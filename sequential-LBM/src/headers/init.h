#include "parameters.h"

#ifndef SEQUENTIAL_LBM_SRC_HEADERS_INIT_H_
#define SEQUENTIAL_LBM_SRC_HEADERS_INIT_H_

template<typename T>
void InitArray(T* array, T init_value, int size) {
    for (int i = 0; i < size; ++i) {
        array[i] = init_value;
    }
}

void InitPopulationField(real* population);

void InitFlagFieldStub(int* flag_field, char* grid_file);

#endif  // SEQUENTIAL_LBM_SRC_HEADERS_INIT_H_

