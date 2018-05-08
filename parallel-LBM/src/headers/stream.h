#include "parameters.h"

#ifndef SEQUENTIAL_LBM_SRC_HEADERS_STREAM_H_
#define SEQUENTIAL_LBM_SRC_HEADERS_STREAM_H_

typedef void (*ptr_stream_func)(int, int, real*, real*);

void Stream(real *population,
            real *swap_buffer,
            ptr_stream_func *stream_element);

void StreamFluid(int i, int j, real *population, real *swap_buffer);
void StreamBC(int i, int j, real *population, real *swap_buffer);
#endif  // SEQUENTIAL_LBM_SRC_HEADERS_STREAM_H_
