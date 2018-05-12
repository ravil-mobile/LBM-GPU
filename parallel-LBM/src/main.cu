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
    int max_num_threads_per_block = 0;
    int max_num_registers_per_block = 0;

    HANDLE_ERROR(cudaGetDeviceCount(&gpu_count));
    printf("number of gpu devices detected: %d\n", gpu_count);

    for (int gpu_instance = 0; gpu_instance < gpu_count; ++ gpu_instance) {
        HANDLE_ERROR(cudaGetDeviceProperties(&property, gpu_instance));
#ifdef DEBUG
        printf(" --- General Information for device %d ---\n", gpu_instance);
        printf("name: %s\n", property.name);
        
        printf("warp size: %d\n", property.warpSize);
        printf("max num. threads per block: %d\n", property.maxThreadsPerBlock);
        printf("max num. registers per block: %d\n", property.regsPerBlock);
        printf("size of constant memory: %d\n", property.totalConstMem);
#endif
        if (property.maxThreadsPerBlock > max_num_threads_per_block) {
            gpu_device = gpu_instance;
            max_num_threads_per_block = property.maxThreadsPerBlock;
            max_num_registers_per_block = property.regsPerBlock;
        }
    }

    // read the enviroment variable "MAX_NUM_THREADS_PER_BLOCK"
    // use the maximum value provide by the DEVICE if the variable has not been defined
    char* str_num_threads = getenv ("MAX_NUM_THREADS_PER_BLOCK");
    if (str_num_threads != NULL) {
        int num_threads = atoi(str_num_threads);
        if (num_threads != 0) {
            max_num_threads_per_block = num_threads;
        }
    } 

    HANDLE_ERROR(cudaSetDevice(gpu_device));
    HANDLE_ERROR(cudaGetDeviceProperties(&property, gpu_device));
    printf("\n --- %s: device has been chosen --- \n", property.name);
    printf(" --- Number threads per block: %d --- \n", max_num_threads_per_block);
    printf(" --- Number registers per block: %d --- \n", property.regsPerBlock);

    struct SimulationParametes parameters;
    struct BoundaryInfo boundary_info;
    struct Constants constants;

    // read input data
    char parameter_file[] = "parameter.txt";
    char boundary_file[] = "boundary.txt";
    char grid_file[] = "grid.txt";
    ReadInputFilesStub(parameters,
                       boundary_info,
                       constants,
                       parameter_file,
                       boundary_file);
    
    // define cuda grid parameters
    const int MAX_NUM_THREADS = max_num_threads_per_block;
    const int MAX_NUM_BLOCKS = (parameters.num_lattices + MAX_NUM_THREADS) / MAX_NUM_THREADS;

    const int MAX_NUM_USED_REGISTERS_PER_WARP = 35;
    const int MIN_NUM_THREADS = max_num_registers_per_block / MAX_NUM_USED_REGISTERS_PER_WARP;
    const int MIN_NUM_BLOCKS = (parameters.num_lattices + MIN_NUM_THREADS) / MIN_NUM_THREADS;

#ifdef DEBUG
    printf(" --- num elements: %d --- \n", parameters.num_lattices);
    printf(" --- max #threads %d: max #blocks: %d --- \n", MAX_NUM_THREADS, MAX_NUM_BLOCKS);
    printf(" --- min #threads %d: min #blocks: %d --- \n", MIN_NUM_THREADS, MIN_NUM_BLOCKS);

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

    SetupGnuPlots(velocity_frame, density_frame, parameters);

    // allocate memory in the HOST
    int *flag_field = (int*)calloc(parameters.num_lattices, sizeof(int));
    real *density = (real*)calloc(parameters.num_lattices, sizeof(real));
    real *velocity_magnitude = (real*)calloc(parameters.num_lattices, sizeof(real));

    InitFlagFieldStub(flag_field,
                      grid_file,
                      parameters);
    
    // allocate and init DOMAIN on the DEVICE
    DomainHandler domain_handler;
    domain_handler.InitDomainOnDevice(parameters,
                                      flag_field,
                                      MAX_NUM_THREADS,
                                      MAX_NUM_BLOCKS);

    const Domain *domain = domain_handler.GetDeviceData();

    // allocate and init BOUNDARIES on the DEVICE
    BoundaryConditionsHandler bc_handler;
    
    ScanFlagField(flag_field,
                  bc_handler,
                  parameters,
                  constants,
                  boundary_info);

    const BoundaryConditions *boundaries = bc_handler.GetDeviceData();

    int threads = 0;
    int blocks = 0;

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    // START of algorighm
    cudaEventRecord(start, 0);
    for (int time = 0; time < parameters.num_time_steps; ++time) {
        
        // perform streaming step 
        StreamDevice<<<MAX_NUM_BLOCKS, MAX_NUM_THREADS>>>(domain->dev_population,
                                                          domain->dev_swap_buffer,
                                                          domain->dev_flag_field);
        //CUDA_CHECK_ERROR(); 
        
        // apply boundary consitions
        if (boundaries->num_wall_elements != 0) {
            threads = min(boundaries->num_wall_elements, MAX_NUM_THREADS);
            blocks = (parameters.num_lattices + threads) / threads;
            TreatNonSlipBC<<<blocks, threads>>>(boundaries->bc_wall_indices,
                                                domain->dev_swap_buffer,
                                                boundaries->num_wall_elements); 
            //CUDA_CHECK_ERROR();
        }

        if (boundaries->num_moving_wall_elements != 0) {
            threads = min(boundaries->num_moving_wall_elements, MAX_NUM_THREADS);
            blocks = (parameters.num_lattices + threads) / threads;
            TreatSlipBC<<<blocks, threads>>>(boundaries->bc_moving_wall_indices,
                                         boundaries->bc_moving_wall_data,
                                         domain->dev_density,
                                         domain->dev_swap_buffer,
                                         boundaries->num_moving_wall_elements);
            //CUDA_CHECK_ERROR();
        }

        if (boundaries->num_inflow_elements != 0) {
            threads = min(boundaries->num_inflow_elements, MAX_NUM_THREADS);
            blocks = (parameters.num_lattices + threads) / threads;
            TreatInflowBC<<<blocks, threads>>>(boundaries->bc_inflow_indices,
                                               boundaries->bc_inflow_data,
                                               domain->dev_density,
                                               domain->dev_swap_buffer,
                                               boundaries->num_inflow_elements);
            //CUDA_CHECK_ERROR();
        }

        if (boundaries->num_outflow_elements != 0) {
            threads = min(boundaries->num_outflow_elements, MAX_NUM_THREADS);
            blocks = (parameters.num_lattices + threads) / threads;
            TreatOutflowBC<<<blocks, threads>>>(boundaries->bc_outflow_indices,
                                                domain->dev_velocity,
                                                domain->dev_density,
                                                domain->dev_swap_buffer,
                                                boundaries->num_outflow_elements);
            //CUDA_CHECK_ERROR();
        }


        HANDLE_ERROR(cudaDeviceSynchronize());
        domain_handler.SwapPopulationFields(); 
       
        // perform collision step 
        UpdateDensityFieldDevice<<<MAX_NUM_BLOCKS, MAX_NUM_THREADS>>>(domain->dev_density,
                                                                      domain->dev_population,
                                                                      domain->dev_flag_field);

        //CUDA_CHECK_ERROR(); 

        
        UpdateVelocityFieldDevice<<<MAX_NUM_BLOCKS, MAX_NUM_THREADS>>>(domain->dev_velocity,
                                                                       domain->dev_population,
                                                                       domain->dev_density,
                                                                       domain->dev_flag_field);
        //CUDA_CHECK_ERROR(); 

        
        UpdatePopulationFieldDevice<<<MIN_NUM_BLOCKS, MIN_NUM_THREADS>>>(domain->dev_velocity,
                                                                         domain->dev_population,
                                                                         domain->dev_density);
        //CUDA_CHECK_ERROR(); 
        
        /*
        threads = 468;
        blocks = (parameters.num_lattices + threads) / threads; 
        UpdatePopulationFieldDevice<<<blocks, threads>>>(domain->dev_velocity,
                                                         domain->dev_population,
                                                         domain->dev_density);
        */

#ifdef DEBUG



        if ((time % parameters.steps_per_report) == 0) {

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
        }
#endif

#ifdef GRAPHICS
        
        if ((time % parameters.steps_per_report) == 0) {
        
            ComputeVelocityMagnitude<<<MAX_NUM_BLOCKS, MAX_NUM_THREADS>>>(domain->dev_velocity,
                                                                          domain->dev_velocity_magnitude);
 

            HANDLE_ERROR(cudaMemcpy(velocity_magnitude,
                                    domain->dev_velocity_magnitude,
                                    parameters.num_lattices * sizeof(real),
                                    cudaMemcpyDeviceToHost));

            DisplayResults(velocity_magnitude, velocity_frame, parameters);
            // DisplayResults(velocity, velocity_frame,
            //               density, density_frame);
        }
#endif
    }
    // END of algorithm
    cudaEventRecord(stop, 0);
    
    cudaEventSynchronize(start);
    cudaEventSynchronize(stop);
    float elapsed_time = 0;
    cudaEventElapsedTime(&elapsed_time, start, stop); 
    
    double MLUPS = (parameters.num_lattices * parameters.num_time_steps)
                 / (elapsed_time * 1e3);

    printf("MLUPS: %4.6f\n", MLUPS);

    // free HOST recourses
    gnuplot_close(velocity_frame);
    gnuplot_close(density_frame);

    
    // free DEVICE resources
    free(flag_field);
    free(density);
    free(velocity_magnitude);
    return 0;
}
