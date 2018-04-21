#include "parameters.h"

#ifndef SEQUENTIAL_LBM_SRC_HEADERS_COLLISION_H_
#define SEQUENTIAL_LBM_SRC_HEADERS_COLLISION_H_

typedef void (*ptr_update_func)(int, real*, real*, real*);

void UpdateDensityField(real *density,
                        real *population,
                        ptr_update_func *update_function);

void UpdateDensityFluid(int index, real *density, real *population, real *velocity = 0);
void UpdateDensityBC(int index, real* density, real *population, real *velocity = 0);

void UpdateVelocityField(real *velocity, 
                         real *population,
                         real *density,
                         ptr_update_func *update_function);


void UpdateVelocityFluid(int index, real *density, real *population, real *velocity);
void UpdateVelocityBC(int index, real *density, real *population, real *velocity);

void UpdatePopulationField(real *velocity,
                           real *population,
                           real *density);

#endif  // SEQUENTIAL_LBM_SRC_HEADERS_COLLISION_H_
