import numpy as np
import matplotlib.pyplot as plt

from helper import *
from update_density import *
from update_velocity import *
from update_population import *

# NOTATION:
#   | i - index
#   |
#   |
#   |__________ j - index
#

def main():
    # initialization of input parameters, configuration and constants
    param = {}
    param["simulation_time"] = 1.0
    param["time_step"] = 0.001
    param["dimension"] = 2
    param["discretization"] = 9
    param["tau"] = 0.8
    param["relaxation"] = 1.0 / param["tau"]
    param["speed_of_sound"] = 1.0 / np.sqrt(3.0)
    param["viscosity"] = param["speed_of_sound"]**2 * (param["tau"] - 0.5)
    param["width"] = 100
    param["height"] = 100
    param["num_lattices"] = param["width"] * param["height"]

    precomputed_constants = {}
    precomputed_constants["const_1"] = 1.0 / (param["speed_of_sound"] * param["speed_of_sound"])
    precomputed_constants["const_2"] = 0.5 * precomputed_constants["const_1"] * precomputed_constants["const_1"]
    precomputed_constants["const_3"] = 0.5 * precomputed_constants["const_1"]

    flags = {"fluid": 0, "wall": 1, "mov_wall": 2, "inflow": 3, "outflow": 4}

    # INFO: keep coordinate as 1D array instead of 2D
    # to speed up computations of dot products
    coords = np.array([0, 0,    #0
                       1, 0,    #1
                       0, 1,    #2
                       -1, 0,   #3
                       0, -1,   #4
                       1, 1,    #5
                       -1, 1,   #6
                       -1, -1,  #7
                       1, -1])  #8

    weights = np.array([4.0 / 9.0,   #0
                        1.0 / 9.0,   #1
                        1.0 / 9.0,   #2
                        1.0 / 9.0,   #3
                        1.0 / 9.0,   #4
                        1.0 / 36.0,  #5
                        1.0 / 36.0,  #6
                        1.0 / 36.0,  #7
                        1.0 / 36.0]) #8

    # arrays in c-style
    flag_filed = np.zeros(param["num_lattices"], dtype=np.uint8)
    density = np.ones(param["num_lattices"], dtype=np.float32)
    velocity = np.zeros(param["dimension"] * param["num_lattices"], dtype=np.float32)
    population = np.ones(param["num_lattices"] * param["discretization"], dtype=np.float32)
    swap_buffer = np.ones(param["num_lattices"] * param["discretization"], dtype=np.float32)

    # list of functions
    update_density_funcs = [update_density_fluid_cell] * param["num_lattices"]
    update_velocity_funcs = [update_velocity_fluid_cell] * param["num_lattices"]

    init_cavity_flag_field(flag_filed,
                           update_density_funcs,
                           update_velocity_funcs,
                           flags,
                           param)


    # run simulation
    update_density_field(density, population, update_density_funcs, param)
    update_velocity_field(velocity, population, density, coords, update_velocity_funcs, param)
    update_population_field(velocity, population, density, coords, weights, precomputed_constants, param)

    display(field = flag_filed, dim=1, parameters=param)


    print_dict(param)


def collide():
    pass


def stream():
    pass


def init_cavity_flag_field(field, update_funcs, update_velocity_funcs, flags, parameters):

    bottom_index = 0
    most_left_index = 0
    top_index = parameters["height"] - 1
    most_right_index = parameters["width"] - 1
    for i in range(parameters["width"]):
        # init top (moving) wall
        index = get_index(index_i=top_index, index_j=i, param=parameters, dim=1)
        field[index] = flags["mov_wall"]
        update_funcs[index] = update_density_bc_cell
        update_velocity_funcs[index] = update_velocity_bc_cell

        # init bottom wall
        index = get_index(index_i=bottom_index, index_j=i, param=parameters, dim=1)
        field[index] = flags["wall"]
        update_funcs[index] = update_density_bc_cell
        update_velocity_funcs[index] = update_velocity_bc_cell

        # init left wall
        index = get_index(index_i=i, index_j=most_left_index, param=parameters, dim=1)
        field[index] = flags["wall"]
        update_funcs[index] = update_density_bc_cell
        update_velocity_funcs[index] = update_velocity_bc_cell

        # init right wall
        index = get_index(index_i=i, index_j=most_right_index, param=parameters, dim=1)
        field[index] = flags["wall"]
        update_funcs[index] = update_density_bc_cell
        update_velocity_funcs[index] = update_velocity_bc_cell


def display(field, parameters, dim=1, shift=0):
    screen = np.zeros((parameters["height"], parameters["width"]))
    for i in range(parameters["height"]):
        for j in range(parameters["width"]):
            screen[i][j] = field[get_index(i, j, parameters, dim) + shift]

    img = plt.imshow(screen, origin='lower')
    plt.colorbar(img)
    plt.show()


def convert_to_2d():
    pass


def print_dict(dictionary):
    print("parameters of the problem:")
    for key in dictionary.keys():
        print(key, ":", dictionary[key])

main()