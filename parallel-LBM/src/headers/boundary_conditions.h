#include <vector>
#include "parameters.h"

#ifndef SEQUENTIAL_LBM_SRC_HEADERS_BOUNDARY_CONDITIONS_H_
#define SEQUENTIAL_LBM_SRC_HEADERS_BOUNDARY_CONDITIONS_H_

typedef void (*ptr_boundary_func)(int, int, real*, real*, real*);

struct InfoBC {
    int component;
    ptr_boundary_func function;
};

struct WallBC {
    int source_index;
    int target_index;
};

struct MovingWallBC {
    int source_index;
    int target_index;
    int scalar_target_index;
    double precomputed_data;
};

struct BoundaryConditions {
    int *bc_wall_indices;
    int num_wall_elements;

    int *bc_moving_wall_indices;
    real *bc_moving_wall_data;
    int num_moving_wall_elements;
};


class BoundaryConditionsHandler {
public:
    BoundaryConditionsHandler();
    ~BoundaryConditionsHandler();

    void SetNonSlipBC(std::vector<struct WallBC> elements);
    void SetSlipBC(std::vector<struct MovingWallBC> elements);    
    const BoundaryConditions * GetDeviceData();
private:
    BoundaryConditions dev_boundary_conditions;
};




void TreatBoundary(ptr_boundary_func *boundary_update,
                   int *boundary_coords,
                   int num_boundaries,
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
