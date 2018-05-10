#include "parameters.h"
#include "boundary_conditions.h"

#ifndef SEQUENTIAL_LBM_SRC_HEADERS_SCANNING_H_
#define SEQUENTIAL_LBM_SRC_HEADERS_SCANNING_H_

void ScanFlagField(int *flag_field,
                   ptr_boundary_func **boundary_update,
                   int **boundary_coords,
                   int &num_boundaries,
                   BoundaryConditionsHandler &bc_handler);

#endif  // SEQUENTIAL_LBM_SRC_HEADERS_SCANNING_H_
