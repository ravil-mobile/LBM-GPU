#include <vector>
#include <stdio.h>
#include "headers/parameters.h"
#include "headers/scanning.h"
#include "headers/boundary_conditions.h"
#include "headers/helper.h"
#include "headers/init.h"

void ScanFlagField(int *flag_field,
                   BoundaryConditionsHandler &bc_handler) {
    int num_directions = parameters.discretization;
   
    std::vector<struct WallBC> wall_elements;
    std::vector<struct MovingWallBC> mov_wall_elements;
    std::vector<struct InflowBC> inflow_elements;
    std::vector<struct OutflowBC> outflow_elements;

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

                        switch (neighbour_flag) {
                            case WALL: 
                                struct WallBC wall;
                                PrecomputeWallBC(component,
                                                 scalar_self_index,
                                                 wall);

                                wall_elements.push_back(wall);
                                break;
                            case MOVING_WALL: 
                                struct MovingWallBC moving_wall;
                                PrecomputeMovingWallBC(component,
                                                       scalar_self_index,
                                                       moving_wall);

                                mov_wall_elements.push_back(moving_wall);
                                break;
                            case INFLOW: 
                                struct InflowBC inlet_lattice;
                                
                                PrecomputeInflowBC(component,
                                                   scalar_self_index,
                                                   inlet_lattice);

                                inflow_elements.push_back(inlet_lattice); 
                                break;
                            case OUTFLOW: 
                                struct OutflowBC outlet_lattice;
                                PrecomputeOutflowBC(component,
                                                    scalar_self_index,
                                                    outlet_lattice);

                                outflow_elements.push_back(outlet_lattice);
                                break;
                            default:
                                break;
                        }

                    }
                }
            }
        }
    }

    bc_handler.SetNonSlipBC(wall_elements);
    bc_handler.SetSlipBC(mov_wall_elements);
    bc_handler.SetInflowBC(inflow_elements);
    bc_handler.SetOutflowBC(outflow_elements);
}
