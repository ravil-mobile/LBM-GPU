#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <iostream>

#include "headers/parameters.h"
#include "headers/init.h"
#include "headers/stub.h"
#include "headers/collision.h"
#include "headers/stream.h"
#include "headers/scanning.h"
#include "headers/boundary_conditions.h"
#include "headers/gnuplot_i.h"
#include "headers/helper.h"
#include "headers/kernels.h"
#include "headers/domain.h"

int main() {

    // choose a proper gpu device with max num threads per block
    cudaDeviceProp property;
    int gpu_count = 0;
    int gpu_device = 0;
    int max_thread_num_per_block = 0;

    HANDLE_ERROR(cudaGetDeviceCount(&gpu_count));
    printf("number of gpu devices detected: %d\n", gpu_count);

    for (int gpu_instance = 0; gpu_instance < gpu_count; ++ gpu_instance) {
        HANDLE_ERROR(cudaGetDeviceProperties(&property, gpu_instance));
#ifdef DEBUG
        printf(" --- General Information for device %d ---\n", gpu_instance);
        printf("name: %s\n", property.name);
        
        printf("warp size: %d\n", property.warpSize);
        printf("max threads per block: %d\n", property.maxThreadsPerBlock);
        printf("size of constant memory: %d\n", property.totalConstMem);
#endif
        if (property.maxThreadsPerBlock > max_thread_num_per_block) {
            gpu_device = gpu_instance;
            max_thread_num_per_block = property.maxThreadsPerBlock;
        }
    }

    // read the enviroment variable "NUM_THREADS_PER_BLOCK"
    // use the maximum value provide by the DEVICE if the variable has not been defined
    char* str_num_threads = getenv ("NUM_THREADS_PER_BLOCK");
    int num_threads = atoi(str_num_threads);
    if (num_threads != 0) {
        max_thread_num_per_block = num_threads;
    }    

    HANDLE_ERROR(cudaSetDevice(gpu_device));
    HANDLE_ERROR(cudaGetDeviceProperties(&property, gpu_device));
    printf("\n --- %s: device has been chosen --- \n", property.name);
    printf(" --- Number threads per block: %d --- \n", max_thread_num_per_block);
        

    // read input data
    char parameter_file[] = "parameter.txt";
    char boundary_file[] = "boundary.txt";
    char grid_file[] = "grid.txt";
    ReadInputFilesStub(parameter_file,
                       boundary_file);
    
    // define cuda grid parameters
    const int NUM_THREADS = max_thread_num_per_block;
    const int NUM_BLOCKS = (parameters.num_lattices + NUM_THREADS) / NUM_THREADS;

#ifdef DEBUG
    printf(" --- num elements: %d --- \n", parameters.num_lattices);
    printf(" --- #threads %d: #blocks: %d --- \n", NUM_THREADS, NUM_BLOCKS);
#endif

    // allocate constant data into the DEVICE
    CopyConstantsToDevice(parameters,
                          constants,
                          boundary_info,
                          coords,
                          weights);

#ifdef DEBUG    
    CheckConstMemoryCopy<<<1,1>>>();
    cudaDeviceSynchronize();
#endif

    gnuplot_ctrl *velocity_frame;
    gnuplot_ctrl *density_frame;

    velocity_frame = gnuplot_init();
    density_frame = gnuplot_init();

    SetupGnuPlots(velocity_frame, density_frame);

    // allocate memory in the HOST
    int *flag_field = (int*)calloc(parameters.num_lattices, sizeof(int));
    real *density = (real*)calloc(parameters.num_lattices, sizeof(real));
    real *velocity = (real*)calloc(parameters.dimension * parameters.num_lattices, sizeof(real));
    real *population = (real*)calloc(parameters.discretization * parameters.num_lattices, sizeof(real));
    real *swap_buffer = (real*)calloc(parameters.discretization * parameters.num_lattices, sizeof(real));

    ptr_update_func *update_density = (ptr_update_func*)calloc(parameters.num_lattices,
                                                               sizeof(ptr_update_func));

    ptr_update_func *update_velocity = (ptr_update_func*)calloc(parameters.num_lattices,
                                                                sizeof(ptr_update_func));

    ptr_stream_func *stream_element = (ptr_stream_func*)calloc(parameters.num_lattices,
                                                                sizeof(ptr_stream_func));
    
    // init memory in the host HOST
    InitArray<real>(density, 1.0, parameters.num_lattices);
    InitPopulationField(population);
    InitPopulationField(swap_buffer);
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
    
    // allocate and init DOMAIN on the DEVICE
    DomainHandler domain_handler;
    domain_handler.InitDomainOnDevice(parameters,
                                      flag_field,
                                      NUM_THREADS,
                                      NUM_BLOCKS);

    Domain *domain = domain_handler.GetDeviceData();

    // allocate and init BOUNDARIES on the DEVICE
    BoundaryConditionsHandler bc_handler;
    
    ScanFlagField(flag_field,
                  &boundary_update,
                  &boundary_coords,
                  num_boundaries,
                  bc_handler);

    const BoundaryConditions* boundaries = bc_handler.GetDeviceData();

    int threads = 0;
    int blocks = 0;

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    // START of algorighm
    cudaEventRecord(start, 0);
    for (int time = 0; time < parameters.num_time_steps; ++time) {
 
        // perform streaming step 
        StreamDevice<<<NUM_BLOCKS, NUM_THREADS>>>(domain->dev_population,
                                                  domain->dev_swap_buffer,
                                                  domain->dev_flag_field);
        CUDA_CHECK_ERROR(); 

        // apply boundary consitions
        if (boundaries->num_wall_elements != 0) {
            threads = min(boundaries->num_wall_elements, NUM_THREADS);
            blocks = (parameters.num_lattices + threads) / threads;
            TreatNonSlipBC<<<blocks, threads>>>(boundaries->bc_wall_indices,
                                                domain->dev_swap_buffer,
                                                boundaries->num_wall_elements); 
            CUDA_CHECK_ERROR();
        }

        if (boundaries->num_moving_wall_elements != 0) {
            threads = min(boundaries->num_moving_wall_elements, NUM_THREADS);
            blocks = (parameters.num_lattices + threads) / threads;
            TreatSlipBC<<<blocks, threads>>>(boundaries->bc_moving_wall_indices,
                                         boundaries->bc_moving_wall_data,
                                         domain->dev_density,
                                         domain->dev_swap_buffer,
                                         boundaries->num_moving_wall_elements);
            CUDA_CHECK_ERROR();
        }

        
        HANDLE_ERROR(cudaDeviceSynchronize());        
        std::swap(domain->dev_population, domain->dev_swap_buffer);
       
        // perform collision step 
        UpdateDensityFieldDevice<<<NUM_BLOCKS, NUM_THREADS>>>(domain->dev_density,
                                                              domain->dev_population,
                                                              domain->dev_flag_field);

        CUDA_CHECK_ERROR(); 

        
        UpdateVelocityFieldDevice<<<NUM_BLOCKS, NUM_THREADS>>>(domain->dev_velocity,
                                                               domain->dev_population,
                                                               domain->dev_density,
                                                               domain->dev_flag_field);
        CUDA_CHECK_ERROR(); 

        
        threads = 960;
        blocks = (parameters.num_lattices + threads) / threads;
        UpdatePopulationFieldDevice<<<blocks, threads>>>(domain->dev_velocity,
                                                         domain->dev_population,
                                                         domain->dev_density);
        CUDA_CHECK_ERROR(); 


#ifdef DEBUG
        HANDLE_ERROR(cudaMemcpy(density,
                                domain->dev_density,
                                parameters.num_lattices * sizeof(real),
                                cudaMemcpyDeviceToHost));

        real max_density = *std::max_element(density,
                                    density + parameters.num_lattices);
        real min_density = *std::min_element(density,
                                density + parameters.num_lattices);


        std::cout << "time step: " << time << "; ";
        std::cout << "max density: " << max_density << "; ";
        std::cout << "min density "  << min_density << std::endl;
#endif

#ifdef GRAPHICS
        HANDLE_ERROR(cudaMemcpy(velocity,
                                domain->dev_velocity,
                                parameters.num_lattices * parameters.dimension * sizeof(real),
                                cudaMemcpyDeviceToHost));

        if ((time % parameters.steps_per_report) == 0) {
            DisplayResults(velocity, velocity_frame);
            // DisplayResults(velocity, velocity_frame,
            //               density, density_frame);
        }
#endif
    }
    cudaEventRecord(stop, 0);
    
    cudaEventSynchronize(start);
    cudaEventSynchronize(stop);
    float elapsed_time = 0;
    cudaEventElapsedTime(&elapsed_time, start, stop); 
    
    
    double MLUPS = (parameters.num_lattices * parameters.num_time_steps)
                 / (elapsed_time * 1e3);

    printf("MLUPS: %4.6f\n", MLUPS);

    // end of algorithm

    // free HOST recourses
    gnuplot_close(velocity_frame);
    gnuplot_close(density_frame);

    
    // free DEVICE resources
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
