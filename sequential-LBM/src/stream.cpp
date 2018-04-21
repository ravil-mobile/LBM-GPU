#include "headers/parameters.h"
#include "headers/helper.h"
#include "headers/stream.h"

void Stream(real *population,
            real *swap_buffer,
            ptr_stream_func *stream_element) {
    for (int j = 0; j < parameters.height; ++j) {
        for (int i = 0; i < parameters.width; ++i) {
            int scalar_index = GetIndex(i, j);
            stream_element[scalar_index](i, j, population, swap_buffer);
        }
    }
}

void StreamFluid(int i, int j, real *population, real *swap_buffer) {
    int self_index = GetIndex(i, j);
    int num_directions = parameters.discretization;
    int num_lattices = parameters.num_lattices;

    for (int component = 0; component < num_directions; ++component) {
        int ii = coords[component];
        int jj = coords[component + num_directions];
        
        int neighbour_index = GetIndex(i + ii, j + jj);
        int shift = component * num_lattices;

        swap_buffer[neighbour_index + shift] = population[self_index + shift];
    }
}

void StreamBC(int i, int j, real *population, real *swap_buffer) {
}
