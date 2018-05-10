#include <unordered_map>
#include <vector>
#include <stdio.h>
#include "headers/parameters.h"
#include "headers/scanning.h"
#include "headers/boundary_conditions.h"
#include "headers/helper.h"
#include "headers/init.h"



//************************************************************************

void PrecomputeWallBC(int component,
                      int coordinate,
                      struct WallBC &element) {
    int num_lattices = parameters.num_lattices;
    int num_directions = parameters.discretization;

    int i = coords[component];
    int j = coords[component + num_directions];

    int inverse_component = inverse_indices[component];

    int shift = GetIndex(i, j);
    element.target_index = coordinate + inverse_component * num_lattices;
    element.source_index = (coordinate + shift) + component * num_lattices;
}

void PrecomputeMovingWallBC(int component,
                            int coordinate,
                            struct MovingWallBC &element) {
    int num_lattices = parameters.num_lattices;
    int num_directions = parameters.discretization;

    int i = coords[component];
    int j = coords[component + num_directions];

    int inverse_component = inverse_indices[component];

    int shift = GetIndex(i, j);
    
    element.scalar_target_index = coordinate;
    element.target_index = coordinate + inverse_component * num_lattices;
    element.source_index = (coordinate + shift) + component * num_lattices;

    real dot_product_cu = constants.one * (i * boundary_info.wall_velocity_x +
                                           j * boundary_info.wall_velocity_y);

    element.precomputed_data = 2.0 * weights[component] * dot_product_cu;
}



void ScanFlagField(int *flag_field,
                   ptr_boundary_func **boundary_update,
                   int **boundary_coords,
                   int &num_boundaries,
                   BoundaryConditionsHandler &bc_handler) {
    int num_directions = parameters.discretization;
    std::unordered_map<int, std::vector<struct InfoBC>> container;
    
    std::vector<struct WallBC> wall_elements;
    std::vector<struct MovingWallBC> mov_wall_elements;

    for (int j = 0; j < parameters.height; ++j) {
        for (int i = 0; i < parameters.width; ++i) {
            int scalar_self_index = GetIndex(i, j);

            if (flag_field[scalar_self_index] == FLUID) {
                for (int component = 0; component < num_directions; ++component) {
                    int ii = coords[component];
                    int jj = coords[component + num_directions];

                    int scalar_neighbour_index = GetIndex(i + ii, j + jj);
                    int shift = component * parameters.num_lattices;
                    
                    int neighbour_flag = flag_field[scalar_neighbour_index];

                    if (neighbour_flag != FLUID) {
                        struct InfoBC boundary_lattice;
                        boundary_lattice.component = component;

                        switch (neighbour_flag) {
                            case WALL: boundary_lattice.function = ApplyNonSlipBC;
                                struct WallBC wall;
                                PrecomputeWallBC(component,
                                                 scalar_self_index,
                                                 wall);

                                wall_elements.push_back(wall);
                                break;
                            case MOVING_WALL: boundary_lattice.function = ApplyMovingWallBC;
                                struct MovingWallBC moving_wall;
                                PrecomputeMovingWallBC(component,
                                                       scalar_self_index,
                                                       moving_wall);

                                mov_wall_elements.push_back(moving_wall);
                                break;
                            case INFLOW: boundary_lattice.function = ApplyInflowBC;
                                break;
                            case OUTFLOW: boundary_lattice.function = ApplyOutflowBC;
                                break;
                            default:
                                boundary_lattice.function = SkipBoundary;
                        }

                        auto search = container.find(scalar_neighbour_index);

                        if (search != container.end()) {
                            container[scalar_neighbour_index].push_back(boundary_lattice);
                        } else {
                            std::vector<struct InfoBC> vector;
                            container[scalar_neighbour_index] = vector;
                            container[scalar_neighbour_index].push_back(boundary_lattice);
                        }
                    }
                }
            }
        }
    }

    num_boundaries = container.size();
    (*boundary_coords) = (int*)calloc(num_boundaries, sizeof(int));
    (*boundary_update) = (ptr_boundary_func*)calloc(num_boundaries * num_directions,
                                                    sizeof(ptr_boundary_func));

    InitArray<ptr_boundary_func>(*boundary_update,
                                 SkipBoundary,
                                 num_boundaries * num_directions);

    int counter = 0;
    for (auto iterator = container.begin(); iterator != container.end(); ++iterator) {
        (*boundary_coords)[counter] = iterator->first;

        std::vector<struct InfoBC> vector = iterator->second;

        for (auto element = vector.begin(); element != vector.end(); ++element) {
            int shift = element->component * num_boundaries;
            (*boundary_update)[counter + shift] = element->function;
        }
        ++counter;
    }

    bc_handler.SetNonSlipBC(wall_elements);
    bc_handler.SetSlipBC(mov_wall_elements);
}
