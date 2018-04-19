#include "parameters.h"

#ifndef SEQUENTIAL_LBM_SRC_HEADERS_BOUNDARY_CONDITIONS_H_
#define SEQUENTIAL_LBM_SRC_HEADERS_BOUNDARY_CONDITIONS_H_

typedef void (*ptr_boundary_func)(int, int, real*, real*, real*);

void TreatBoundary(ptr_boundary_func *boundary_func,
                   int *boundary_coords,
                   real *population,
                   real *velocity,
                   real *density);

void SkipBoundary(int component,
                  int coordinate,
                  real *population,
                  real *velocity,
                  real *density);

void ApplyNonSlipBC(int component,
                    int coordinate,
                    real *population,
                    real *velocity,
                    real *density);

void ApplyMovingWallBC(int component,
                       int coordinate,
                       real *population,
                       real *velocity,
                       real *density);

void ApplyInflowBC(int component,
                   int coordinate,
                   real *population,
                   real *velocity,
                   real *density);

void ApplyOutflowBC(int component,
                    int coordinate,
                    real *population,
                    real *velocity,
                    real *density);

#endif  // SEQUENTIAL_LBM_SRC_HEADERS_BOUNDARY_CONDITIONS_H_
