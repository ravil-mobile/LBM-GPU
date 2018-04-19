#include "parameters.h"

#ifndef SEQUENTIAL_LNB_INIT_H_
#define SEQUENTIAL_LNB_INIT_H_

template<typename T>
void init_array(T* array, T init_value, int size) {
    for(int i = 0; i < size; ++i) {
        array[i] = init_value;
    }
}

void init_population_field(real* population);

#endif  // SEQUENTIAL_LNB_INIT_H_

