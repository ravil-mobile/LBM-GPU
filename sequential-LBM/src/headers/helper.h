#include <string>
#include "parameters.h"
#include "gnuplot_i.h"

#ifndef SEQUENTIAL_LBM_SRC_HEADERS_HELPER_H_
#define SEQUENTIAL_LBM_SRC_HEADERS_HELPER_H_

int GetIndex(int index_i, int index_j, int dim = 1);

real ComputeVectorMagnitude(real a, real b);

void SetupGnuPlots(gnuplot_ctrl *velocity_frame,
                   gnuplot_ctrl *density_frame);

void DisplayResults(real *velocity, gnuplot_ctrl *velocity_frame,
                    real *density = 0, gnuplot_ctrl *density_frame = 0);

#endif  // SEQUENTIAL_LBM_SRC_HEADERS_HELPER_H_
