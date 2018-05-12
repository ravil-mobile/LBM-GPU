#include <string>
#include "parameters.h"
#include "gnuplot_i.h"

#include <stdio.h>
#ifndef SEQUENTIAL_LBM_SRC_HEADERS_HELPER_H_
#define SEQUENTIAL_LBM_SRC_HEADERS_HELPER_H_

int GetIndex(int index_i, int index_j, const int width);

void SetupGnuPlots(gnuplot_ctrl *velocity_frame,
                   gnuplot_ctrl *density_frame,
                   const struct SimulationParametes &parameters);

void DisplayResults(real *velocity,
                    gnuplot_ctrl *velocity_frame,
                    const struct SimulationParametes &parameters,
                    real *density = 0, 
                    gnuplot_ctrl *density_frame = 0);

template<typename T>
void PrintArray(T *array, int size) {
    for (int i = 0; i < size; ++i) {
        printf("element %d value -> %f\n", i, array[i]);
    }
}

template<typename T>
void VerifyArray(T *array, T ref_value, int size) {
    for (int i = 0; i < size; ++i) {
        if (array[i] != ref_value) {
            printf("VERIFICATION: failure\n");
            return;
        }
    }
    printf("VERIFICATION: success\n");
}

template<typename T>
void CheckArrayConsistency(T *array, T* ref_array, int size) {
    for (int i = 0; i < size; ++i) {
        if (array[i] != ref_array[i]) {
             printf("VERIFICATION: failure\n");
            return;
        }
    }
    printf("CONSISTENCY: success\n");
}
#endif  // SEQUENTIAL_LBM_SRC_HEADERS_HELPER_H_
