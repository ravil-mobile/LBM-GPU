#include <unordered_map>
#include <vector>
#include "headers/parameters.h"
#include "headers/scanning.h"
#include "headers/boundary_conditions.h"
#include "headers/helper.h"
#include <stdio.h>

struct InfoBC {
    int component;
    ptr_boundary_func function;
};

void ScanFlagField(int *flag_field,
                   ptr_boundary_func **boundary_update,
                   int **boundary_coords) {
    int num_directions = parameters.discretization;
    std::unordered_map<int, std::vector<struct InfoBC>> container;

    for (int j = 0; j < parameters.height; ++j) {
        for (int i = 0; i < parameters.width; ++i) {
            int scalar_self_index = GetIndex(i, j);

            if (flag_field[scalar_self_index] == FLUID) {
                for (int component = 0; component < num_directions; ++component) {
                    int ii = coords[component];
                    int jj = coords[component + num_directions];

                    int scalar_neighbour_index = GetIndex(i + ii, j + jj);
                    int neighbour_flag = flag_field[scalar_neighbour_index];
                    

                    if (neighbour_flag != FLUID) {
                        struct InfoBC boundary_lattice;
                        boundary_lattice.component = component;

                        switch (neighbour_flag) {
                            case WALL: boundary_lattice.function = ApplyNonSlipBC;
                                break;
                            case MOVING_WALL: boundary_lattice.function = ApplyMovingWallBC;
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
                            container[scalar_self_index].push_back(boundary_lattice);
                        }
                    }   
                }
            }
        }
    }

    int num_boundaries = container.size();
    *boundary_coords = (int*)calloc(num_boundaries, sizeof(int));
    *boundary_update = (ptr_boundary_func*)calloc(num_boundaries, sizeof(ptr_boundary_func));

    int counter = 0;
    for (auto iterator = container.begin(); iterator != container.end(); ++iterator) {
        (*boundary_coords)[counter++] = iterator->first;
        
        std::vector<struct InfoBC> vector = iterator->second;

        for (auto element = vector.begin(); element > vector.end(); ++element) {
            int shift = element->component * num_boundaries;
            (*boundary_update)[counter + shift] = element->function;
        }
    }
}
