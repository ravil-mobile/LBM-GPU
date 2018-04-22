#include "headers/parameters.h"
#include "headers/collision.h"
#include "headers/helper.h"
#include <stdio.h>

void UpdateDensityField(real *density,
                        real *population,
                        ptr_update_func *update_function) {
    real* null_ptr = 0;
    for (int j = 0; j < parameters.height; ++j) {
        for (int i = 0; i < parameters.width; ++i) {
            int scalar_index = GetIndex(i, j);
            update_function[scalar_index](scalar_index, 
                                          density,
                                          population,
                                          null_ptr);
        }
    }
}
void UpdateDensityFluid(int index, 
                        real *density, 
                        real *population, 
                        real *velocity) {
    real lattice_density = 0.0;
    int num_directions = parameters.discretization;
    int num_lattices = parameters.num_lattices;

    for (int component = 0; component < num_directions; ++component) {
        lattice_density += population[index + component * num_lattices];
    }
    density[index] = lattice_density;
}

void UpdateDensityBC(int index, 
                     real *density, 
                     real *population, 
                     real* velocity) {
}

void UpdateVelocityField(real *velocity,
                         real *population,
                         real *density,
                         ptr_update_func *update_function) {
    for (int j = 0; j < parameters.height; ++j) {
        for(int i = 0; i < parameters.width; ++i) {
            int scalar_index = GetIndex(i, j);
            update_function[scalar_index](scalar_index,
                                          density,
                                          population,
                                          velocity);
        }
    }
}
#include <stdio.h>
void UpdateVelocityFluid(int index,
                         real *density,
                         real *population,
                         real *velocity) {
    real lattice_velocity_x = 0.0;
    real lattice_velocity_y = 0.0;
    int num_directions = parameters.discretization;
    int num_lattices = parameters.num_lattices;
    
    for (int component = 0; component < num_directions; ++component) {
        real distribution = population[index + component * num_lattices];
        lattice_velocity_x += coords[component] * distribution;
        lattice_velocity_y += coords[component + num_directions] * distribution;
    }
    
    real inverse_density = 1.0 / density[index];
    velocity[index] = inverse_density * lattice_velocity_x;
    velocity[index + num_lattices] = inverse_density * lattice_velocity_y;
}

void UpdateVelocityBC(int index,
                      real *density,
                      real *population,
                      real *velocity) {
}

void UpdatePopulationField(real *velocity,
                           real *population,
                           real *density) {
    int num_directions = parameters.discretization;
    int num_lattices = parameters.num_lattices;
    int height = parameters.height;
    int width = parameters.width;
    real relaxation = parameters.relaxation;
    real const_one = constants.one;
    real const_two = constants.two;
    real const_three = constants.three;

    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {
            int scalar_index = GetIndex(i, j);
            
            real local_velocity_x = velocity[scalar_index];
            real local_velocity_y = velocity[scalar_index + num_lattices];
            
            real dot_product_uu = local_velocity_x * local_velocity_x
                                + local_velocity_y * local_velocity_y;

            for (int component = 0; component < num_directions; ++component) {
                real vector_component_x = coords[component];
                real vector_component_y = coords[component + num_directions];
                
                real dot_product_cu = vector_component_x * local_velocity_x
                                    + vector_component_y * local_velocity_y;
                
                real vector_expansion = const_one * dot_product_cu
                                      + const_two * dot_product_cu * dot_product_cu
                                      - const_three * dot_product_uu
                                      + 1.0;
                
                real equilibrium = weights[component] 
                                 * density[scalar_index] 
                                 * vector_expansion;
            
                int shift = component * num_lattices;
                population[scalar_index + shift] -= (relaxation
                                                  * (population[scalar_index + shift] - equilibrium)); 
            }
        }
    }
}
