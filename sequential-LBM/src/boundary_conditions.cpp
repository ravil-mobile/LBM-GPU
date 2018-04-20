#include "headers/parameters.h"
#include "headers/boundary_conditions.h"

void TreatBoundary(ptr_boundary_func *boundary_func,
                   int *boundary_coords,
                   real *population,
                   real *velocity,
                   real *density) {


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

}

void ApplyMovingWallBC(int component,
                       int coordinate,
                       real *population,
                       real *velocity,
                       real *density) {

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
