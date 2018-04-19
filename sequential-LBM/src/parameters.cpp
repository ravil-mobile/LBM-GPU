#include "parameters.h"

struct SimulationParametes parameters;
struct Constants constants;
struct BoundaryInfo boundary_info;

int coords[] = {0, 1, 0, -1,  0, 1, -1, -1,  1,
                0, 0, 1,  0, -1, 1,  1, -1, -1};

real weights[] = {4.0 / 9.0,
                  1.0 / 9.0,
                  1.0 / 9.0,
                  1.0 / 9.0,
                  1.0 / 9.0,
                  1.0 / 36.0,
                  1.0 / 36.0,
                  1.0 / 36.0,
                  1.0 / 36.0};

int inverse_indices[] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
