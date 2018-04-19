#include <stdio.h>
#include <stdlib.h>
#include "parameters.h"
#include "init.h"
#include "stub.h"

typedef void (*func_ptr)();

void dummy() {};

int main() {
    char parameter_file[] = "parameter.txt";
    char boundary_file[] = "boundary.txt";
    char grid_file[] = "grid.txt";

    ReadInputFilesStub(parameter_file,
                       boundary_file);

    int *flag_field = (int*)calloc(parameters.num_lattices, sizeof(int));
    real *density = (real*)calloc(parameters.num_lattices, sizeof(real));
    real *velocity = (real*)calloc(parameters.dimension * parameters.num_lattices, sizeof(real));
    real *population = (real*)calloc(parameters.discretization * parameters.num_lattices, sizeof(real));
    real *swap_buffer = (real*)calloc(parameters.discretization * parameters.num_lattices, sizeof(real));

    init_array<real>(density, 1.0, parameters.num_lattices);
    init_population_field(population);
    init_population_field(swap_buffer);

    func_ptr *functions = (func_ptr*)calloc(parameters.num_lattices, sizeof(func_ptr));
    init_array<func_ptr>(functions, dummy, parameters.num_lattices);

    printf("num of lattices: %i\n", parameters.num_lattices);
        
    free(flag_field);
    free(density);
    free(velocity);
    free(population);
    free(swap_buffer);

    //free(functions);
    return 0;
}
