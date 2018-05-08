#include "parameters.h"
#include "collision.h"
#include "stream.h"

#ifndef SEQUENTIAL_LBM_SRC_HEADERS_INIT_H_
#define SEQUENTIAL_LBM_SRC_HEADERS_INIT_H_

template<typename T>
void InitArray(T* array, T init_value, int size) {
    for (int i = 0; i < size; ++i) {
        array[i] = init_value;
    }
}

void InitPopulationField(real* population);

void InitFlagFieldStub(int* flag_field,
                       ptr_update_func *update_density,
                       ptr_update_func *update_velocity,
                       ptr_stream_func *stream_element,
                       char* grid_file);

#endif  // SEQUENTIAL_LBM_SRC_HEADERS_INIT_H_

