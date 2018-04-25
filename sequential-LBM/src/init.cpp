#include <math.h>
#include "headers/init.h"
#include "headers/parameters.h"
#include "headers/helper.h"
#include "headers/collision.h"

void InitPopulationField(real* population) {
    int shift = 0;

    for (int component = 0; component < parameters.discretization; ++component) {
        shift = parameters.num_lattices * component;
        for (int i = 0; i < parameters.num_lattices; ++i) {
            population[i + shift] = weights[component];
        }
    }
}

void InitFlagFieldStub(int* flag_field,
                       ptr_update_func *update_density,
                       ptr_update_func *update_velocity,
                       ptr_stream_func *stream_element,
                       char* grid_file) {
    int most_left_index = 0;
    int most_right_index = parameters.width - 1;

    for (int i = 0; i < parameters.height; ++i) {
        // init left wall
        int index = GetIndex(most_left_index, i);
        flag_field[index] = INFLOW;
        update_density[index] = UpdateDensityBC;
        update_velocity[index] = UpdateVelocityBC;
        stream_element[index] = StreamBC; 

        // init right wall
        index = GetIndex(most_right_index, i);
        flag_field[index] = OUTFLOW;
        update_density[index] = UpdateDensityBC;
        update_velocity[index] = UpdateVelocityBC;
        stream_element[index] = StreamBC; 
    }

    int bottom_index = 0;
    int top_index = parameters.height - 1;
    for (int i = 0; i < parameters.width; ++i) {
        // init top (moving) wall
        int index = GetIndex(i, top_index);
        flag_field[index] = WALL;
        update_density[index] = UpdateDensityBC;
        update_velocity[index] = UpdateVelocityBC;
        stream_element[index] = StreamBC; 

        // init bottom wall
        index = GetIndex(i, bottom_index);
        flag_field[index] = WALL;
        update_density[index] = UpdateDensityBC;
        update_velocity[index] = UpdateVelocityBC;
        stream_element[index] = StreamBC; 
    }

    /*
    int obstacle[] = {10, 15, 10, 15};
    for (int j = obstacle[2]; j < obstacle[3]; ++j) {
        for (int i = obstacle[0]; i < obstacle[1]; ++i) {
            int index = GetIndex(i, j);
            flag_field[index] = WALL;
            update_density[index] = UpdateDensityBC;
            update_velocity[index] = UpdateVelocityBC;
            stream_element[index] = StreamBC; 
        }
    }
    */

    int center[] = {315, 315};
    int radius = 30;
    for (int j = (center[1] - radius); j < (center[1] + radius); ++j) {
        for (int i = (center[0] - radius); i < (center[0] + radius); ++i) {
            real delta_x = real(center[0] - i);
            real delta_y = real(center[1] - j);
            real distance = sqrt(delta_x * delta_x + delta_y * delta_y);
            if (real(radius) > distance) {
                int index = GetIndex(i, j);
                index = GetIndex(i, j);
                flag_field[index] = WALL;
                update_density[index] = UpdateDensityBC;
                update_velocity[index] = UpdateVelocityBC;
                stream_element[index] = StreamBC;     
            }
        }
    }
}
