#include "headers/parameters.h"
#include "headers/helper.h"
#include "headers/boundary_conditions.h"
#include "headers/kernels.h"

BoundaryConditionsHandler::BoundaryConditionsHandler() {
    dev_boundary_conditions.bc_wall_indices = 0;
    dev_boundary_conditions.num_wall_elements = 0;

    dev_boundary_conditions.bc_moving_wall_indices = 0;
    dev_boundary_conditions.bc_moving_wall_data = 0;
    dev_boundary_conditions.num_moving_wall_elements = 0;

    dev_boundary_conditions.bc_inflow_indices = 0;
    dev_boundary_conditions.bc_inflow_data = 0;
    dev_boundary_conditions.num_inflow_elements = 0;

    dev_boundary_conditions.bc_outflow_indices = 0;
    dev_boundary_conditions.num_outflow_elements = 0;

}

BoundaryConditionsHandler::~BoundaryConditionsHandler() {
    if (dev_boundary_conditions.bc_wall_indices != NULL) {
        HANDLE_ERROR(cudaFree(dev_boundary_conditions.bc_wall_indices));
    }
        
    if (dev_boundary_conditions.bc_moving_wall_data != NULL) {
        HANDLE_ERROR(cudaFree(dev_boundary_conditions.bc_moving_wall_data));
    }
        
    if (dev_boundary_conditions.bc_moving_wall_indices != NULL) {
        HANDLE_ERROR(cudaFree(dev_boundary_conditions.bc_moving_wall_indices));
    }

    if (dev_boundary_conditions.bc_inflow_data != NULL) {
        HANDLE_ERROR(cudaFree(dev_boundary_conditions.bc_inflow_data));
    }
        
    if (dev_boundary_conditions.bc_inflow_indices != NULL) {
        HANDLE_ERROR(cudaFree(dev_boundary_conditions.bc_inflow_indices));
    }
    
    if (dev_boundary_conditions.bc_outflow_indices != NULL) {
        HANDLE_ERROR(cudaFree(dev_boundary_conditions.bc_outflow_indices));
    }
}

void BoundaryConditionsHandler::SetNonSlipBC(std::vector<struct WallBC> elements) {
    
    // concert data from AoS to SoA 
    int size = elements.size();
    int *temp_indices = (int*)calloc(2 * size, sizeof(int));

    for (int i = 0; i < size; ++i) {
        temp_indices[i] = elements[i].source_index;
        temp_indices[i + size] = elements[i].target_index;
    } 
    
    // allocate memory on the DEVICE
    HANDLE_ERROR(cudaMalloc(&dev_boundary_conditions.bc_wall_indices,
                            2 * size * sizeof(int)));
    
    // copy data from HOST to DEVICE
    HANDLE_ERROR(cudaMemcpy(dev_boundary_conditions.bc_wall_indices,
                            temp_indices,
                            2 * size * sizeof(int),
                            cudaMemcpyHostToDevice));
    
    dev_boundary_conditions.num_wall_elements = size;
    
    // free resources on the HOST 
    free(temp_indices);
}

void BoundaryConditionsHandler::SetSlipBC(std::vector<struct MovingWallBC> elements) {
   
    // concert data from AoS to SoA 
    int size = elements.size(); 
    int *temp_indices = (int*)calloc(3 * size, sizeof(int));
    real* temp_data = (real*)calloc(size, sizeof(real));

    for (int i = 0; i < size; ++i) {
        temp_indices[i] = elements[i].source_index;
        temp_indices[i + size] = elements[i].target_index;
        temp_indices[i + 2 * size] = elements[i].scalar_target_index;

        temp_data[i] = elements[i].precomputed_data;
    }
    
    // allocate memory on the DEVICE
    HANDLE_ERROR(cudaMalloc(&dev_boundary_conditions.bc_moving_wall_indices,
                            3 * size * sizeof(int)));

    HANDLE_ERROR(cudaMalloc(&dev_boundary_conditions.bc_moving_wall_data,
                            size * sizeof(real)));

    // copy data from HOST to DEVICE
    HANDLE_ERROR(cudaMemcpy(dev_boundary_conditions.bc_moving_wall_indices,
                            temp_indices,
                            3 * size * sizeof(int),
                            cudaMemcpyHostToDevice));

    HANDLE_ERROR(cudaMemcpy(dev_boundary_conditions.bc_moving_wall_data,
                            temp_data,
                            size * sizeof(real),
                            cudaMemcpyHostToDevice));
    
    dev_boundary_conditions.num_moving_wall_elements = size;

    // free resources on the HOST 
    free(temp_indices);
    free(temp_data);
}

void BoundaryConditionsHandler::SetInflowBC(std::vector<struct InflowBC> elements) {

    // concert data from AoS to SoA 
    int size = elements.size(); 
    int *temp_indices = (int*)calloc(2 * size, sizeof(int));
    real* temp_data = (real*)calloc(size, sizeof(real));

    for (int i = 0; i < size; ++i) {
        temp_indices[i] = elements[i].target_index;
        temp_indices[i + size] = elements[i].scalar_target_index;

        temp_data[i] = elements[i].precomputed_data;
    }
    
    // allocate memory on the DEVICE
    HANDLE_ERROR(cudaMalloc(&dev_boundary_conditions.bc_inflow_indices,
                            2 * size * sizeof(int)));

    HANDLE_ERROR(cudaMalloc(&dev_boundary_conditions.bc_inflow_data,
                            size * sizeof(real)));

    // copy data from HOST to DEVICE
    HANDLE_ERROR(cudaMemcpy(dev_boundary_conditions.bc_inflow_indices,
                            temp_indices,
                            2 * size * sizeof(int),
                            cudaMemcpyHostToDevice));

    HANDLE_ERROR(cudaMemcpy(dev_boundary_conditions.bc_inflow_data,
                            temp_data,
                            size * sizeof(real),
                            cudaMemcpyHostToDevice));
    
    dev_boundary_conditions.num_inflow_elements = size;

    // free resources on the HOST 
    free(temp_indices);
    free(temp_data);
}


void BoundaryConditionsHandler::SetOutflowBC(std::vector<struct OutflowBC> elements) {

    // concert data from AoS to SoA 
    int size = elements.size(); 
    int *temp_indices = (int*)calloc(4 * size, sizeof(int));

    for (int i = 0; i < size; ++i) {
        temp_indices[i] = elements[i].source_index;
        temp_indices[i + size] = elements[i].scalar_target_index;
        temp_indices[i + 2 * size] = elements[i].target_index;
        temp_indices[i + 3 * size] = elements[i].component;
    }
    
    // allocate memory on the DEVICE
    HANDLE_ERROR(cudaMalloc(&dev_boundary_conditions.bc_outflow_indices,
                            4 * size * sizeof(int)));

    // copy data from HOST to DEVICE
    HANDLE_ERROR(cudaMemcpy(dev_boundary_conditions.bc_outflow_indices,
                            temp_indices,
                            4 * size * sizeof(int),
                            cudaMemcpyHostToDevice));

    dev_boundary_conditions.num_outflow_elements = size;

     // free resources on the HOST 
    free(temp_indices);
}

    
const BoundaryConditions * BoundaryConditionsHandler::GetDeviceData() {
    return &dev_boundary_conditions;
}


void TreatBoundary(ptr_boundary_func *boundary_update,
                   int *boundary_coords,
                   int num_boundaries,
                   real *population,
                   real *velocity,
                   real *density) {
    for (int component = 0; component < parameters.discretization; ++component) {
        for (int index = 0; index < num_boundaries; ++index) {
            int shift = component * num_boundaries;
            boundary_update[index + shift](component,
                                           boundary_coords[index],
                                           population,
                                           velocity,
                                           density);
        }
    }
}


void PrecomputeWallBC(int component,
                      int coordinate,
                      struct WallBC &element) {
    int num_lattices = parameters.num_lattices;
    int num_directions = parameters.discretization;

    int i = coords[component];
    int j = coords[component + num_directions];

    int inverse_component = inverse_indices[component];

    int shift = GetIndex(i, j);
    element.target_index = coordinate + inverse_component * num_lattices;
    element.source_index = (coordinate + shift) + component * num_lattices;
}

void PrecomputeMovingWallBC(int component,
                            int coordinate,
                            struct MovingWallBC &element) {
    int num_lattices = parameters.num_lattices;
    int num_directions = parameters.discretization;

    int i = coords[component];
    int j = coords[component + num_directions];

    int inverse_component = inverse_indices[component];

    int shift = GetIndex(i, j);
    
    element.scalar_target_index = coordinate;
    element.target_index = coordinate + inverse_component * num_lattices;
    element.source_index = (coordinate + shift) + component * num_lattices;

    real dot_product_cu = constants.one * (i * boundary_info.wall_velocity_x +
                                           j * boundary_info.wall_velocity_y);

    element.precomputed_data = 2.0 * weights[component] * dot_product_cu;
}

void PrecomputeInflowBC(int component,
                        int coordinate,
                        struct InflowBC &element) {
    int num_lattices = parameters.num_lattices;
    int num_directions = parameters.discretization;

    int i = coords[component];
    int j = coords[component + num_directions];
    
    int inverse_component = inverse_indices[component];
    int shift = GetIndex(i, j);
    
    element.scalar_target_index = coordinate;
    element.target_index = coordinate + inverse_component * num_lattices;

    real dot_product_uu = boundary_info.velocity_inflow_x * boundary_info.velocity_inflow_x
                        + boundary_info.velocity_inflow_y * boundary_info.velocity_inflow_y;

    real dot_product_cu = -1 * (i * boundary_info.velocity_inflow_x + 
                                j * boundary_info.velocity_inflow_y);

    real const_one = constants.one;
    real const_two = constants.two;
    real const_three = constants.three;

    real velocity_expansion = const_one * dot_product_cu
                            + const_two * dot_product_cu * dot_product_cu
                            - const_three * dot_product_uu
                            + 1.0;

    element.precomputed_data = weights[component] * velocity_expansion;
}

void PrecomputeOutflowBC(int component,
                         int coordinate,
                         struct OutflowBC &element) {

    int num_lattices = parameters.num_lattices;
    int num_directions = parameters.discretization;

    int i = coords[component];
    int j = coords[component + num_directions];
    
    int inverse_component = inverse_indices[component];
    int shift = GetIndex(i, j);
        
    element.source_index = (coordinate + shift) + component * num_lattices;
    element.scalar_target_index = coordinate;
    element.target_index = coordinate + inverse_component * num_lattices;
    element.component = component;
}

