#include <string>
#include "parameters.h"

#ifndef SEQUENTIAL_LBM_SRC_HEADERS_HELPER_H_
#define SEQUENTIAL_LBM_SRC_HEADERS_HELPER_H_

int GetIndex(int index_i, int index_j, int dim = 1);

void WriteVectorFieldToFile(real *data, std::string file_name);

#endif  // SEQUENTIAL_LBM_SRC_HEADERS_HELPER_H_
