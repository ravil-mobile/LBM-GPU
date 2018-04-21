#include "headers/parameters.h"
#include "headers/helper.h"
#include "headers/boundary_conditions.h"

void TreatBoundary(ptr_boundary_func *boundary_update,
                   int *boundary_coords,
                   int num_boundaries,
                   real *population,
                   real *velocity,
                   real *density) {
    for (int component = 0; component < parameters.discretization; ++component) {
        for (int index = 0; index < num_boundaries; ++index) {
            int shift = component * num_boundaries;
            boundary_update[index + shift](component,
                                           boundary_coords[index],
                                           population,
                                           velocity,
                                           density);
        }
    }
}

#include <stdio.h>
void SkipBoundary(int component,
                  int coordinate,
                  real *population,
                  real *velocity,
                  real *density) {
}

void ApplyNonSlipBC(int component,
                    int coordinate,
                    real *population,
                    real *velocity,
                    real *density) {
    int num_lattices = parameters.num_lattices;
    int num_directions = parameters.discretization;

    int i = -1 * coords[component];
    int j = -1 * coords[component + num_directions];

    int inverse_component = inverse_indices[component];

    int shift = GetIndex(i, j);
    int self_index = coordinate + component * num_lattices;
    int neighbour_index = (coordinate + shift) + inverse_component * num_lattices;
    
    population[neighbour_index] = population[self_index];
}

void ApplyMovingWallBC(int component,
                       int coordinate,
                       real *population,
                       real *velocity,
                       real *density) {
    int num_lattices = parameters.num_lattices;
    int num_directions = parameters.discretization;

    int i = -1 * coords[component];
    int j = -1 * coords[component + num_directions];

    int inverse_component = inverse_indices[component];
    
    int shift = GetIndex(i, j);
    int self_index = coordinate + component * num_lattices;
    int neighbour_index = (coordinate + shift) + inverse_component * num_lattices;

    int scalar_neighbour_index = coordinate + shift;

    real const_one = constants.one;
    real const_two = constants.two;
    real const_three = constants.three;
    
    real dot_product_cu = const_one * (real(i) * boundary_info.wall_velocity_x + 
                                       real(j) * boundary_info.wall_velocity_y);
  
    real complement = 2.0 * weights[component] * density[scalar_neighbour_index] * dot_product_cu; 
    
    population[neighbour_index] = population[self_index] + complement; 
}

void ApplyInflowBC(int component,
                   int coordinate,
                   real *population,
                   real *velocity,
                   real *density) {

}

void ApplyOutflowBC(int component,
                    int coordinate,
                    real *population,
                    real *velocity,
                    real *density) {

}
