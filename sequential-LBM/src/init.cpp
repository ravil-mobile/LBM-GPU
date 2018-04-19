#include "headers/init.h"
#include "headers/parameters.h"
#include "headers/helper.h"

void init_population_field(real* population) {
    unsigned shift = 0;

    for (unsigned component = 0; component < parameters.discretization; ++component) {
        shift = parameters.num_lattices * component;
        for (unsigned i = 0; i < parameters.num_lattices; ++i) {
            population[i + shift] = weights[component];
        }
    }
}

void InitFlagField(int* flag_field, char* grid_file) {
    unsigned most_left_index = 0;
    unsigned most_right_index = parameters.width - 1;

    unsigned dim = 1;
    int index = 0;

    for (int i = 0; i < parameters.height; ++i) {
        // init left wall
        index = GetIndex(most_left_index, i, dim);
        flag_field[index] = WALL;

        // init right wall
        index = GetIndex(most_right_index, i, dim);
        flag_field[index] = WALL;
    }

    int bottom_index = 0;
    int top_index = parameters.height - 1;
    for (int i = 0; i < parameters.width; ++i) {
        // init top (moving) wall
        index = GetIndex(i, top_index, dim);
        flag_field[i] = MOVING_WALL;

        // init bottom wall
        index = GetIndex(i, bottom_index, dim);
        flag_field[i] = WALL;
    }

    int obstacle[] = {10, 15, 10, 15};
    for (int j = obstacle[2]; j < obstacle[3]; ++j) {
        for (int i = obstacle[0]; i < obstacle[1]; ++i) {
            flag_field[i] = WALL;
        }
    }
}
