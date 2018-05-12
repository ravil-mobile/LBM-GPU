#include "parameters.h"
#include "boundary_conditions.h"
#ifndef SEQUENTIAL_LBM_SRC_HEADERS_KERNELS_H_
#define SEQUENTIAL_LBM_SRC_HEADERS_KERNELS_H_

void CUDA_CHECK_ERROR(); 
void HANDLE_ERROR(cudaError_t);

void CopyConstantsToDevice(struct SimulationParametes parameters,
                           struct Constants constants,
                           struct BoundaryInfo boundary_info,
                           int *coords,
                           real *weights);

__global__ void CheckConstMemoryCopy();

__global__ void InitArrayDevice(real *array, real init_value, int size);

__global__ void StreamDevice(real *population,
                             real *swap_buffer,
                             int *flag_field);

__global__ void UpdateDensityFieldDevice(real *density,
                                         real *population,
                                         int *flag_field);

__global__ void UpdateVelocityFieldDevice(real *velocity,
                                          real *population,
                                          real *density,
                                          int *flag_field);

__global__ void UpdatePopulationFieldDevice(real *velocity,
                                            real *population,
                                            real *density);

__global__ void PrintBC(struct BoundaryConditions* boundary_conditions);

__global__ void TreatNonSlipBC(int *indices,
                               real *population,
                               int size);

__global__ void TreatSlipBC(int *indices,
                            real *data,
                            real *density,
                            real *population,
                            int size); 

__global__ void TreatInflowBC(int *indices,
                              real *data,
                              real *density,
                              real *population,
                              int size);

__global__ void TreatOutflowBC(int *indices,
                               real *velocity,
                               real *density,
                               real *population,
                               int size);

__global__ void ComputeVelocityMagnitude(real *velocity,
                                         real *velocity_magnitude);

#endif  // SEQUENTIAL_LBM_SRC_HEADERS_KERNELS_H_
