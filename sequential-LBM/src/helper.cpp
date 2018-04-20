#include "headers/parameters.h"
#include "headers/helper.h"

int GetIndex(int index_i, int index_j, int dim) {
    return dim * (index_i + index_j * parameters.width);
}
