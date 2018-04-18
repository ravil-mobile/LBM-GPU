#ifndef SEQUENTIAL_LNB_PARAMETERS_H_
#define SEQUENTIAL_LNB_PARAMETERS_H_

typedef double real;

struct SimulationParametes {
    real simulation_time;
    unsigned num_time_steps;
    unsigned dimension;
    unsigned discretization;
    real delta_x;
    real delta_t;
    real speed_of_sound;
    real viscosity;
    real tau;
    real relaxation;
    unsigned width;
    unsigned height;
    unsigned num_lattices;
} parameters;

struct Constants {
    real one;
    real two;
    real three;
} constants;

struct BoundaryInfo {
    real wall_velocity_x;
    real wall_velocity_y;
    real density_inflow;
    real velocity_inflow_x;
    real velocity_inflow_y;
    real density_outflow;
    real velocity_outflow_x;
    real velocity_outflow_y;
} boundary_info;

enum flags {FLUID, WALL, MOVING_WALL, INFLOW, OUTFLOW};
extern int coords[];
extern real weights[];
extern int inverse_indices[]; 

#endif  // SEQUENTIAL_LNB_PARAMETERS_H_
