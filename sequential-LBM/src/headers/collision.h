#include "parameters.h"

#ifndef SEQUENTIAL_LBM_SRC_HEADERS_COLLISION_H_
#define SEQUENTIAL_LBM_SRC_HEADERS_COLLISION_H_

typedef void (*ptr_update_func)(int, real*);

void UpdateDensityField(real *density,
                        real *population,
                        ptr_update_func *update_function);

void UpdateDensityFluid(int index, real *population);
void UpdateDensityBC(int index, real *population);

void UpdateVelocityField(real *velocity, 
                         real *population,
                         real *density,
                         ptr_update_func *update_function);


void UpdateVelocityFluid(int index, real *population);
void UpdateVelocityBC(int index, real *population);

void UpdatePopulationField(real *velocity,
                           real *population,
                           real *density);

#endif  // SEQUENTIAL_LBM_SRC_HEADERS_COLLISION_H_
