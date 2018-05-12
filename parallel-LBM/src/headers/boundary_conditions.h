#include <vector>
#include "parameters.h"

#ifndef SEQUENTIAL_LBM_SRC_HEADERS_BOUNDARY_CONDITIONS_H_
#define SEQUENTIAL_LBM_SRC_HEADERS_BOUNDARY_CONDITIONS_H_

typedef void (*ptr_boundary_func)(int, int, real*, real*, real*);

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

struct InflowBC {
    int source_index;
    int target_index;
    int scalar_target_index;
    double precomputed_data;
};

struct OutflowBC {
    int source_index;
    int component;
    int target_index;
    int scalar_target_index;
};

struct BoundaryConditions {
    int *bc_wall_indices;
    int num_wall_elements;

    int *bc_moving_wall_indices;
    real *bc_moving_wall_data;
    int num_moving_wall_elements;

    int *bc_inflow_indices;
    real *bc_inflow_data;
    int num_inflow_elements;

    int *bc_outflow_indices;
    int num_outflow_elements;
};


class BoundaryConditionsHandler {
public:
    BoundaryConditionsHandler();
    ~BoundaryConditionsHandler();

    void SetNonSlipBC(std::vector<struct WallBC> elements);
    void SetSlipBC(std::vector<struct MovingWallBC> elements);    
    void SetInflowBC(std::vector<struct InflowBC> elements); 
    void SetOutflowBC(std::vector<struct OutflowBC> elements); 
    const BoundaryConditions * GetDeviceData();
private:
    BoundaryConditions dev_boundary_conditions;
};

void PrecomputeWallBC(int component,
                      int coordinate,
                      struct WallBC &element);

void PrecomputeMovingWallBC(int component,
                            int coordinate,
                            struct MovingWallBC &element);

void PrecomputeInflowBC(int component,
                        int coordinate,
                        struct InflowBC &element);

void PrecomputeOutflowBC(int component,
                         int coordinate,
                         struct OutflowBC &element);

#endif  // SEQUENTIAL_LBM_SRC_HEADERS_BOUNDARY_CONDITIONS_H_
