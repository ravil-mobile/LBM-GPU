#include <math.h>
#include "headers/parameters.h"
#include "headers/stub.h"

void ReadInputFilesStub(char *parameter_file,
                        char *boundary_file) {
    InitParametersStub();
    InitBoundaryConditionStub();
}

void InitParametersStub() {
    parameters.simulation_time = 1.0;
    parameters.num_time_steps = 5000;
    parameters.dimension = 2;
    parameters.discretization = 9;
    parameters.delta_x = 0.577 * 1e-3;
    parameters.delta_t = 1e-3;
    parameters.speed_of_sound = parameters.delta_x / parameters.delta_t;
    parameters.viscosity = 85.06e-6;
    parameters.tau = 0.5 * (1.0 + 6.0 * (parameters.viscosity *
                                         parameters.delta_t) /
                                       pow(parameters.delta_x, 2.0));
    parameters.relaxation = 1.0 / parameters.tau;
    parameters.width = 80;
    parameters.height = 25;
    parameters.num_lattices = parameters.width * parameters.height;

    constants.one = 1.0 / pow(parameters.speed_of_sound, 2.0);
    constants.two = 0.5 * pow(constants.one, 2.0);
    constants.three = 0.5 * constants.one;
}

void InitBoundaryConditionStub() {
    boundary_info.wall_velocity_x = 0.05;
    boundary_info.wall_velocity_y = 0.00;

    boundary_info.density_inflow = 1.00;
    boundary_info.velocity_inflow_x = 0.26;
    boundary_info.velocity_inflow_y = 0.00;

    boundary_info.density_outflow = 1.00;
    boundary_info.velocity_outflow_x = 0.26;
    boundary_info.velocity_outflow_y = 0.00;
}
