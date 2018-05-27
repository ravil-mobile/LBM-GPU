#include <string>
#include "parameters.h"

#include <stdio.h>
#ifndef SEQUENTIAL_LBM_SRC_HEADERS_HELPER_H_
#define SEQUENTIAL_LBM_SRC_HEADERS_HELPER_H_

int GetIndex(int index_i, int index_j, const int width);

struct Arguments {
    char *parameters_file;
    char *boundary_cond_file;
    char *mesh_file;
    char *threads_distr_file;
};

int parse_opt(int key, char *arg, struct argp_state *state);

bool isSpace (const std::string& str);

bool predicate (const std::string& line);

void ReadParameterFile(const char *parameter_file,
                       struct SimulationParametes &parameters);

void ReadBoundaryFile(const char *boundary_file,
                      struct BoundaryInfo &boundary_info);

void ReadMeshFile(const char *mesh_file,
                  int *flag_field,
                  struct SimulationParametes &parameters);

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
