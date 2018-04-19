#include "headers/parameters.h"
#include "headers/collision.h"

typedef void (*ptr_update_func)(int, real*);

void UpdateDensityField(real *density,
                        real *population,
                        ptr_update_func *update_function) {
}

void UpdateDensityFluid(int index, real *population) {
}

void UpdateDensityBC(int index, real *population) {
}

void UpdateVelocityField(real *velocity, 
                         real *population,
                         real *density,
                         ptr_update_func *update_function) {
}


void UpdateVelocityFluid(int index, real *population) {
}

void UpdateVelocityBC(int index, real *population) {
}

void UpdatePopulationField(real *velocity,
                           real *population,
                           real *density) {

}


