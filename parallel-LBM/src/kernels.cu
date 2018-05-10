#include <stdio.h>
#include "headers/kernels.h"
#include "headers/parameters.h"
#include "headers/boundary_conditions.h"

__constant__ struct SimulationParametes parameters_device;
__constant__ struct Constants constants_device;
__constant__ struct BoundaryInfo boundary_info_device;
__constant__ int coords_device[18];
__constant__ real weights_device[9];

void CUDA_CHECK_ERROR() {
    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess) {
        printf("CUDA ERROR: %s\n", cudaGetErrorString(error));
    }
}

void HANDLE_ERROR(cudaError_t error) {
    if (error != cudaSuccess) {
        printf("CUDA STATUS: %s\n", cudaGetErrorString(error));
        exit(1);
    }
}

void CopyConstantsToDevice(struct SimulationParametes parameters,
                           struct Constants constants,
                           struct BoundaryInfo boundary_info,
                           int *coords,
                           real *weights) {
    HANDLE_ERROR(cudaMemcpyToSymbol(parameters_device,
                                    &parameters,
                                    sizeof(struct SimulationParametes)));

    HANDLE_ERROR(cudaMemcpyToSymbol(constants_device,
                                    &constants,
                                    sizeof(struct Constants)));

    HANDLE_ERROR(cudaMemcpyToSymbol(boundary_info_device,
                                    &boundary_info,
                                    sizeof(struct BoundaryInfo)));

    HANDLE_ERROR(cudaMemcpyToSymbol(coords_device,
                                    coords,
                                    18 * sizeof(int)));

    HANDLE_ERROR(cudaMemcpyToSymbol(weights_device,
                                    weights,
                                    9 * sizeof(real)));
} 

__global__ void CheckConstMemoryCopy() {
    int thread_id = threadIdx.x + blockIdx.x * blockDim.x;
    if (thread_id == MASTER) {
        printf(" --- Values allocated in the constant memory--- \n");
        printf(" -> num_lattices %d\n", parameters_device.num_lattices);
        printf(" -> width %d\n", parameters_device.width);

        printf(" -> constant one %f\n", constants_device.one);
        printf(" -> wall velocity %f\n", boundary_info_device.wall_velocity_x);
        printf(" -> weight %f\n\n", weights_device[0]);
    }
}


__device__ int GetIndexDevice(int index_i, int index_j) {
    return index_i + index_j * parameters_device.width; 
}

__global__ void SwapFields(real *a, real *b) {
    int thread_id = threadIdx.x + blockIdx.x * blockDim.x;
    if (thread_id == MASTER) {
        real *temp = a;
        a = b;
        b = temp;
    }
} 

__global__ void InitArrayDevice(real *array, real init_value, int size) {
    int thread_id = threadIdx.x + blockIdx.x * blockDim.x; 
    while (thread_id < size) { 
        array[thread_id] = init_value;
        thread_id += blockDim.x * gridDim.x;
    }
}

__global__ void StreamDevice(real *population,
                             real *swap_buffer,
                             int *flag_field) {
    int num_lattices = parameters_device.num_lattices;
    int num_directions = parameters_device.discretization; 
    int thread_id = threadIdx.x + blockIdx.x * blockDim.x;
    
    while (thread_id < num_lattices) {
        if (flag_field[thread_id] == FLUID) {
            for (int component = 0; component < num_directions; ++component) {
                int ii = coords_device[component];
                int jj = coords_device[num_directions + component];
                int neighbour_index = thread_id + GetIndexDevice(ii, jj);
                int shift = component * num_lattices;

                swap_buffer[shift + neighbour_index] = population[shift + thread_id];
            }
        }
        thread_id += blockDim.x * gridDim.x;
    }
}


__global__ void UpdateDensityFieldDevice(real *density,
                                         real *population,
                                         int *flag_field) {
    int num_lattices = parameters_device.num_lattices;
    int num_directions = parameters_device.discretization; 
    
    int thread_id = threadIdx.x + blockIdx.x * blockDim.x;
    
    while (thread_id < num_lattices) {
        if (flag_field[thread_id] == FLUID) {
            real lattice_density = 0.0;
            for (int component = 0; component < num_directions; ++component) {
                lattice_density += population[component * num_lattices + thread_id];
            }
            density[thread_id] = lattice_density;
        }
        thread_id += blockDim.x * gridDim.x; 
    }
}

__global__ void UpdateVelocityFieldDevice(real *velocity,
                                          real *population,
                                          real *density,
                                          int *flag_field) {
    int num_lattices = parameters_device.num_lattices;
    int num_directions = parameters_device.discretization; 
    
    int thread_id = threadIdx.x + blockIdx.x * blockDim.x;
    
    while (thread_id < num_lattices) {
        if (flag_field[thread_id] == FLUID) {
            real lattice_velocity_x = 0.0;
            real lattice_velocity_y = 0.0;
            
            for (int component = 0; component < num_directions; ++component) {
                real distribution = population[component * num_lattices + thread_id];
                lattice_velocity_x += coords_device[component] * distribution;
                lattice_velocity_y += coords_device[num_directions + component] * distribution;
            }

            real inverse_density = 1.0 / density[thread_id];
            velocity[thread_id] = inverse_density * lattice_velocity_x;
            velocity[num_lattices + thread_id] = inverse_density * lattice_velocity_y;
        }
        thread_id += blockDim.x * gridDim.x; 
    }
}

__global__ void UpdatePopulationFieldDevice(real *velocity,
                                            real *population,
                                            real *density) {
    int num_lattices = parameters_device.num_lattices;
    int num_directions = parameters_device.discretization;  

    real relaxation = parameters_device.relaxation;
    real const_one = constants_device.one;
    real const_two = constants_device.two;
    real const_three = constants_device.three;

    int thread_id = threadIdx.x + blockIdx.x * blockDim.x;
    
    while (thread_id < num_lattices) {
        real local_velocity_x = velocity[thread_id];
        real local_velocity_y = velocity[num_lattices + thread_id];

        real dot_product_uu = local_velocity_x * local_velocity_x
                            + local_velocity_y * local_velocity_y;

        for (int component = 0; component < num_directions; ++component) {
            real vector_component_x = coords_device[component];
            real vector_component_y = coords_device[num_directions + component];

            real dot_product_cu = vector_component_x * local_velocity_x
                                + vector_component_y * local_velocity_y;

            real velocity_expansion = (const_one * dot_product_cu)
                                    + (const_two * dot_product_cu * dot_product_cu)
                                    - (const_three * dot_product_uu)
                                    + 1.0;

            real equilibrium = weights_device[component]
                             * density[thread_id]
                             * velocity_expansion;

            int shift = component * num_lattices;
            population[shift + thread_id] -= (relaxation
                                           * (population[shift + thread_id] - equilibrium));
        }

        thread_id += blockDim.x * gridDim.x;
    }
} 

__global__ void PrintBC(struct BoundaryConditions *boundary_conditions) {
    int thread_id = threadIdx.x + blockIdx.x * blockDim.x;
    if (thread_id == MASTER) {
        printf("GPU: num walls %d\n", boundary_conditions->num_wall_elements);
        printf("GPU: num moving walls %d\n", boundary_conditions->num_moving_wall_elements);

    }
}

__global__ void TreatNonSlipBC(int *indices,
                               real *population,
                               int size) {
    int thread_id = threadIdx.x + blockIdx.x * blockDim.x;
    
    while (thread_id < size) {
        int source = indices[thread_id];
        int target = indices[size + thread_id];
        population[target] = population[source];

        thread_id += blockDim.x * gridDim.x;
    }
}

__global__ void TreatSlipBC(int *indices,
                            real *data,
                            real *density,
                            real *population,
                            int size) {
    int thread_id = threadIdx.x + blockIdx.x * blockDim.x;
    
    while (thread_id < size) {
        int source = indices[thread_id];
        int target = indices[size + thread_id];
        int index = indices[2 * size + thread_id];

        population[target] = population[source] + data[thread_id] * density[index];

        thread_id += blockDim.x * gridDim.x;
    }
}

