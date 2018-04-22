#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <stdio.h>
#include <time.h>

#include "headers/parameters.h"
#include "headers/init.h"
#include "headers/stub.h"
#include "headers/collision.h"
#include "headers/stream.h"
#include "headers/scanning.h"
#include "headers/boundary_conditions.h"
#include "headers/gnuplot_i.h"
#include "headers/helper.h"

int main() {
    char parameter_file[] = "parameter.txt";
    char boundary_file[] = "boundary.txt";
    char grid_file[] = "grid.txt";
    ReadInputFilesStub(parameter_file,
                       boundary_file);
   
    gnuplot_ctrl *velocity_frame;
    gnuplot_ctrl *density_frame;

    velocity_frame = gnuplot_init();
    density_frame = gnuplot_init();

    SetupGnuPlots(velocity_frame, density_frame);

    int steps_per_report = 10;


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
    

    InitFlagFieldStub(flag_field, 
                      update_density,
                      update_velocity,
                      stream_element,        
                      grid_file);
    
    int num_boundaries = 0;
    ptr_boundary_func *boundary_update = 0;
    int *boundary_coords = 0;

    ScanFlagField(flag_field,
                  &boundary_update,
                  &boundary_coords,
                  num_boundaries);

    clock_t begin = clock();    
    for (int time = 0; time < parameters.num_time_steps; ++time) {
        Stream(population, swap_buffer, stream_element);
        
        TreatBoundary(boundary_update,
                      boundary_coords,
                      num_boundaries,
                      swap_buffer,
                      velocity,
                      density);
        
        std::swap(population, swap_buffer);
        
        UpdateDensityField(density, population, update_density);
        
        UpdateVelocityField(velocity,
                            population,
                            density,
                            update_velocity);
        
        
        UpdatePopulationField(velocity,
                              population,
                              density); 
                
#ifdef DEBUG
        real max_density = *std::max_element(density, 
                                    density + parameters.num_lattices);
        real min_density = *std::min_element(density, 
                                density + parameters.num_lattices);
       

        std::cout << "time step: " << time << "; ";
        std::cout << "max density: " << max_density << "; ";
        std::cout << "min density "  << min_density << std::endl;
#endif

#ifdef GRAPHICS
        if ((time % steps_per_report) == 0) { 
            DisplayResults(velocity, velocity_frame,
                           density, density_frame);

        }
#endif
    }
    clock_t end = clock() - begin;
    double consumed_time = (double)(end - begin) / CLOCKS_PER_SEC;
    double MLUPS = ( parameters.num_lattices * parameters.num_time_steps) 
                 / (consumed_time * 1e6 );

    printf("MLUPS: %4.6f\n", MLUPS);

    gnuplot_close(velocity_frame); 
    gnuplot_close(density_frame);

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
