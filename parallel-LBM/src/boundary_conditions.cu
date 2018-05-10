#include "headers/parameters.h"
#include "headers/helper.h"
#include "headers/boundary_conditions.h"
#include "headers/kernels.h"

BoundaryConditionsHandler::BoundaryConditionsHandler() {
    dev_boundary_conditions.bc_wall_indices = 0;
    dev_boundary_conditions.num_wall_elements = 0;

    dev_boundary_conditions.bc_moving_wall_indices = 0;
    dev_boundary_conditions.bc_moving_wall_data = 0;
    dev_boundary_conditions.num_moving_wall_elements = 0;
}

BoundaryConditionsHandler::~BoundaryConditionsHandler() {
    if (dev_boundary_conditions.bc_wall_indices != NULL) {
        HANDLE_ERROR(cudaFree(dev_boundary_conditions.bc_wall_indices));
    }
        
    if (dev_boundary_conditions.bc_moving_wall_data != NULL) {
        HANDLE_ERROR(cudaFree(dev_boundary_conditions.bc_moving_wall_data));
    }
        
    if (dev_boundary_conditions.bc_moving_wall_indices != NULL) {
        HANDLE_ERROR(cudaFree(dev_boundary_conditions.bc_moving_wall_indices));
    }
}

void BoundaryConditionsHandler::SetNonSlipBC(std::vector<struct WallBC> elements) {
    
    // concert data from AoS to SoA 
    int size = elements.size();
    int *temp_indices = (int*)calloc(2 * size, sizeof(int));

    for (int i = 0; i < size; ++i) {
        temp_indices[i] = elements[i].source_index;
        temp_indices[i + size] = elements[i].target_index;
    } 
    
    // allocate memory on the DEVICE
    HANDLE_ERROR(cudaMalloc(&dev_boundary_conditions.bc_wall_indices,
                            2 * size * sizeof(int)));
    
    // copy data from HOST to DEVICE
    HANDLE_ERROR(cudaMemcpy(dev_boundary_conditions.bc_wall_indices,
                            temp_indices,
                            2 * size * sizeof(int),
                            cudaMemcpyHostToDevice));
    
    dev_boundary_conditions.num_wall_elements = size;
    
    // free resources on the HOST 
    free(temp_indices);
}

void BoundaryConditionsHandler::SetSlipBC(std::vector<struct MovingWallBC> elements) {
   
    // concert data from AoS to SoA 
    int size = elements.size(); 
    int *temp_indices = (int*)calloc(3 * size, sizeof(int));
    real* temp_data = (real*)calloc(size, sizeof(real));

    for (int i = 0; i < size; ++i) {
        temp_indices[i] = elements[i].source_index;
        temp_indices[i + size] = elements[i].target_index;
        temp_indices[i + 2 * size] = elements[i].scalar_target_index;

        temp_data[i] = elements[i].precomputed_data;
    }
    
    // allocate memory on the DEVICE
    HANDLE_ERROR(cudaMalloc(&dev_boundary_conditions.bc_moving_wall_indices,
                            3 * size * sizeof(int)));

    HANDLE_ERROR(cudaMalloc(&dev_boundary_conditions.bc_moving_wall_data,
                            size * sizeof(real)));

    // copy data from HOST to DEVICE
    HANDLE_ERROR(cudaMemcpy(dev_boundary_conditions.bc_moving_wall_indices,
                            temp_indices,
                            3 * size * sizeof(int),
                            cudaMemcpyHostToDevice));

    HANDLE_ERROR(cudaMemcpy(dev_boundary_conditions.bc_moving_wall_data,
                            temp_data,
                            size * sizeof(real),
                            cudaMemcpyHostToDevice));
    
    dev_boundary_conditions.num_moving_wall_elements = size;

    // free resources on the HOST 
    free(temp_indices);
    free(temp_data);
}

const BoundaryConditions * BoundaryConditionsHandler::GetDeviceData() {
    return &dev_boundary_conditions;
}


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
/*
    int num_lattices = parameters.num_lattices;
    int num_directions = parameters.discretization;

    int i = -1 * coords[component];
    int j = -1 * coords[component + num_directions];

    int inverse_component = inverse_indices[component];

    int shift = GetIndex(i, j);
    int self_index = coordinate + component * num_lattices;
    int neighbour_index = (coordinate + shift) + inverse_component * num_lattices;

    population[neighbour_index] = population[self_index];
*/
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

