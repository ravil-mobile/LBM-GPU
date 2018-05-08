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
    int num_lattices = parameters.num_lattices;
    int num_directions = parameters.discretization;

    int i = -1 * coords[component];
    int j = -1 * coords[component + num_directions];

    int inverse_component = inverse_indices[component];

    int shift = GetIndex(i, j);
    // int self_index = coordinate + component * num_lattices;
    int scalar_neighbour_index = coordinate + shift;
    int neighbour_index = scalar_neighbour_index + inverse_component * num_lattices;

    real dot_product_uu = boundary_info.velocity_inflow_x * boundary_info.velocity_inflow_x
                        + boundary_info.velocity_inflow_y * boundary_info.velocity_inflow_y;

    real dot_product_cu = i * boundary_info.velocity_inflow_x
                        + j * boundary_info.velocity_inflow_y;

    real const_one = constants.one;
    real const_two = constants.two;
    real const_three = constants.three;

    real velocity_expansion = const_one * dot_product_cu
                            + const_two * dot_product_cu * dot_product_cu
                            - const_three * dot_product_uu
                            + 1.0;

    real equilibrium = weights[inverse_component]
                     * density[scalar_neighbour_index]
                     * velocity_expansion;

    population[neighbour_index] = equilibrium;
}

void ApplyOutflowBC(int component,
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

    int scalar_neighbour_index = coordinate + shift;
    int neighbour_index = scalar_neighbour_index + inverse_component * num_lattices;

    real neighbour_velocity_x = velocity[scalar_neighbour_index];
    real neighbour_velocity_y = velocity[scalar_neighbour_index + num_lattices];

    real dot_product_uu = neighbour_velocity_x * neighbour_velocity_x
                        + neighbour_velocity_y * neighbour_velocity_y;

    real dot_product_cu = -i * neighbour_velocity_x - j * neighbour_velocity_y;
    real dot_product_cu_inv = i * neighbour_velocity_x + j * neighbour_velocity_y;

    real const_one = constants.one;
    real const_two = constants.two;
    real const_three = constants.three;

    real velocity_expansion = (const_one * dot_product_cu)
                            + (const_two * dot_product_cu * dot_product_cu)
                            - (const_three * dot_product_uu)
                            + 1.0;

    real velocity_expansion_inv = (const_one * dot_product_cu_inv)
                                + (const_two * dot_product_cu_inv * dot_product_cu_inv)
                                - (const_three * dot_product_uu)
                                + 1.0;

    real equilibrium = weights[component]
                     * density[scalar_neighbour_index]
                     * velocity_expansion;

    real equilibrium_inv = weights[inverse_component]
                         * density[scalar_neighbour_index]
                         * velocity_expansion_inv;

    population[neighbour_index] = equilibrium
                                + equilibrium_inv
                                - population[self_index];
}
