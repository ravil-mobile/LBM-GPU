#include "parameters.h"

#ifndef SEQUENTIAL_LBM_SRC_HEADERS_HELPER_H_
#define SEQUENTIAL_LBM_SRC_HEADERS_HELPER_H_

int GetIndex(int index_i, int index_j, unsigned dim) {
    return dim * (index_i + index_j * parameters.width);
}

#endif  // SEQUENTIAL_LBM_SRC_HEADERS_HELPER_H_
