#include "parameters.h"
#include "boundary_conditions.h"

#ifndef SEQUENTIAL_LBM_SRC_HEADERS_SCANNING_H_
#define SEQUENTIAL_LBM_SRC_HEADERS_SCANNING_H_

void ScanFlagField(int *flag_field,
                   BoundaryConditionsHandler &bc_handler);

#endif  // SEQUENTIAL_LBM_SRC_HEADERS_SCANNING_H_
