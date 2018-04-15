import numpy as np


class SimulationParametes:

    simulation_time = 0.0
    num_time_steps = 0
    dimension = 2
    discretization = 9

    delta_x = 0.0
    delta_t = 0.0
    speed_of_sound = 0.0
    viscosity = 0.0
    tau = 0.0
    relaxation = 0.0

    width = 0
    height = 0
    num_lattices = 0

    obstacle = 0


class BoundaryInfo:

    wall_velocity = [0.0, 0.0]

    density_inflow = 0.0
    velocity_inflow = [0.0, 0.0]
    inflow_direction = 0

    density_outflow = 0.0
    velocity_outflow = [0.0, 0.0]
    outflow_direction = 0


class Constants:

    const_1 = 0.0
    const_2 = 0.0
    const_3 = 0.0

    coords = np.array([0, 1, 0, -1,  0, 1, -1, -1, 1,
                       0, 0, 1,  0, -1, 1,  1, -1,-1])

    '''
    coords = np.array([ 0,  0,   # 0
                        1,  0,   # 1
                        0,  1,   # 2
                       -1,  0,   # 3
                        0, -1,   # 4
                        1,  1,   # 5
                       -1,  1,   # 6
                       -1, -1,   # 7
                        1, -1])  # 8
    '''

    weights = np.array([4.0 / 9.0,   # 0
                        1.0 / 9.0,   # 1
                        1.0 / 9.0,   # 2
                        1.0 / 9.0,   # 3
                        1.0 / 9.0,    # 4
                        1.0 / 36.0,   # 5
                        1.0 / 36.0,   # 6
                        1.0 / 36.0,   # 7
                        1.0 / 36.0])  # 8

    inverse_indices = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6])

    flags = {"fluid": 0, "wall": 1, "mov_wall": 2, "inflow": 3, "outflow": 4}
