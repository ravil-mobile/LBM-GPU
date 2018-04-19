#include <stdio.h>
#include <stdlib.h>
#include "headers/parameters.h"
#include "headers/init.h"
#include "headers/stub.h"
#include "headers/collision.h"
#include "headers/stream.h"
#include "headers/scanning.h"
#include "headers/boundary_conditions.h"

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

    InitArray<real>(density, 1.0, parameters.num_lattices);
    InitPopulationField(population);
    InitPopulationField(swap_buffer);

    ptr_update_func *update_density = (ptr_update_func*)calloc(parameters.num_lattices, 
                                                               sizeof(ptr_update_func));
    
    ptr_update_func *update_velocity = (ptr_update_func*)calloc(parameters.num_lattices,
                                                                sizeof(ptr_update_func));

    ptr_stream_func *stream_element = (ptr_stream_func*)calloc(parameters.num_lattices,
                                                                sizeof(ptr_stream_func));
    
    InitArray<ptr_update_func>(update_density, UpdateDensityFluid, parameters.num_lattices);
    InitArray<ptr_update_func>(update_velocity, UpdateVelocityFluid, parameters.num_lattices);
    InitArray<ptr_stream_func>(stream_element, StreamFluid, parameters.num_lattices);
    
    InitFlagFieldStub(flag_field, grid_file);
    

    ptr_boundary_func *boundary_update = 0;
    int *boundary_coords = 0;

    ScanFlagField(flag_field,
                  &boundary_update,
                  &boundary_coords);



    printf("num of lattices: %i\n", parameters.num_lattices);



    free(flag_field);
    free(density);
    free(velocity);
    free(population);
    free(swap_buffer);

    free(update_density);
    free(update_velocity);
    free(stream_element);
    free(boundary_update);
    free(boundary_coords);
    return 0;
}
