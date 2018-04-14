import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import time

from helper import *
from update_density import *
from update_velocity import *
from update_population import *
from streaming import *
from boundary_conditions import *
from scanning_boundaries import *

# NOTATION:
#   | j - index
#   |
#   |
#   |__________ i - index
#

#      indices scheme
#
#        6   2   5
#         \  |  /
#          \ | /
#   3 ______\|/______ 1
#           /|\
#          / | \
#         /  |  \
#        7   4   8


def main():
    # initialization of input parameters, configuration and constants
    param = {}
    param["simulation_time"] = 1.0
    param["num_time_steps"] = 500
    param["dimension"] = 2
    param["discretization"] = 9

    param["tau"] = 1.5
    param["relaxation"] = 1.0 / param["tau"]
    param["speed_of_sound"] = 1.0 / np.sqrt(3.0)
    param["viscosity"] = param["speed_of_sound"]**2 * (param["tau"] - 0.5)

    param["width"] = 100
    param["height"] = 25
    param["num_lattices"] = param["width"] * param["height"]

    param["obstacle"] = [20, 25, 10, 15]

    param["wall_velocity"] = [0.05, 0.0]

    param["density_inflow"] = 1.0
    param["velocity_inflow"] = [0.1, 0.0]
    param["inflow_direction"] = 1

    param["density_outflow"] = 1.0
    param["velocity_outfow"] = [0.1, 0.0]
    param["outflow_direction"] = 1

    bc_conditions = {}
    bc_conditions["lid_velocity"] = 0.1

    precomputed_constants = {}
    precomputed_constants["const_1"] = 1.0 / (param["speed_of_sound"] * param["speed_of_sound"])
    precomputed_constants["const_2"] = 0.5 * precomputed_constants["const_1"] * precomputed_constants["const_1"]
    precomputed_constants["const_3"] = 0.5 * precomputed_constants["const_1"]

    flags = {"fluid": 0, "wall": 1, "mov_wall": 2, "inflow": 3, "outflow": 4}

    # INFO: keep coordinate as 1D array instead of 2D
    # to speed up computations of dot products
    coords = np.array([ 0,  0,   #0
                        1,  0,   #1
                        0,  1,   #2
                       -1,  0,   #3
                        0, -1,   #4
                        1,  1,   #5
                       -1,  1,   #6
                       -1, -1,   #7
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

    inverse_indices = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6])

    # arrays in c-style
    flag_field = np.zeros(param["num_lattices"], dtype=np.uint8)
    density = np.ones(param["num_lattices"], dtype=np.float32)
    velocity = np.zeros(param["dimension"] * param["num_lattices"], dtype=np.float32)
    population = np.ones(param["num_lattices"] * param["discretization"], dtype=np.float32)
    swap_buffer = np.ones(param["num_lattices"] * param["discretization"], dtype=np.float32)

    # list of functions
    update_density_funcs = [update_density_fluid_cell] * param["num_lattices"]
    update_velocity_funcs = [update_velocity_fluid_cell] * param["num_lattices"]
    stream_funcs = [stream_fluid_cell] * param["num_lattices"]

    # initialization of all fields
    init_cavity_flag_field(flag_field,
                           flags,
                           update_density_funcs,
                           update_velocity_funcs,
                           stream_funcs,
                           param)


    init_cavity_population_field(population, weights, param)

    boundaries, boundary_coordinates = scanning(flag_field, flags, coords, param)

    fig = plt.figure(1)
    fig.clf()
    ax = fig.add_subplot(1, 1, 1)

    # run simulation
    start = time.time()
    for step in range(param["num_time_steps"]):

        # perform streaming
        stream(population, swap_buffer, coords, stream_funcs, param)

        treat_boundaries(boundaries,
                         bc_conditions,
                         boundary_coordinates,
                         swap_buffer,
                         velocity,
                         density,
                         coords,
                         inverse_indices,
                         weights,
                         precomputed_constants,
                         param)


        # swap buffers
        temp = population
        population = swap_buffer
        swap_buffer = temp

        # computer collision
        update_density_field(density, population, update_density_funcs, param)
        update_velocity_field(velocity, population, density, coords, update_velocity_funcs, param)
        update_population_field(velocity, population, density, coords, weights, precomputed_constants, param)


        #display_scalar_field(field=density, figure=fig, axis=ax, dim=1, parameters=param, shift=0)
        #display_vector_magnitude_2d(field=velocity, figure=fig, axis=ax, parameters=param)
        display_vector_field_2d(field=velocity, figure=fig, axis=ax, parameters=param)
        draw_obstacle(ax, param)
        plt.pause(0.01)

        print("iteration step: %i; density: max = %f; min = %f" % (step, np.max(density), np.min(density)))


    elapsed_time = time.time()- start
    MLUPS = (param["num_lattices"] * param["num_time_steps"]) / (elapsed_time * 1e6)
    print("MLUPS: ", MLUPS)
    print("elapsed time", elapsed_time)



def treat_boundaries(boundaries,
                     bc_conditions,
                     boundary_coordinates,
                     population,
                     velocity,
                     density,
                     coords,
                     inverse_indices,
                     weights,
                     constants,
                     parameters):

    num_components = parameters["discretization"]

    for counter, boundary_function in enumerate(boundaries):
        component = counter % num_components
        scalar_counter = counter // num_components

        boundary_function(component,
                          bc_conditions,
                          boundary_coordinates[scalar_counter],
                          population,
                          velocity,
                          density,
                          coords,
                          inverse_indices,
                          weights,
                          constants,
                          parameters)


def init_cavity_flag_field(field, flags, update_funcs, update_velocity_funcs, stream_funcs, parameters):

    most_left_index = 0
    most_right_index = parameters["width"] - 1

    for j in range(parameters["height"]):

        # init left wall
        index = get_index(index_i=most_left_index, index_j=j, param=parameters, dim=1)
        field[index] = flags["inflow"]
        update_funcs[index] = update_density_bc_cell
        update_velocity_funcs[index] = update_velocity_bc_cell
        stream_funcs[index] = stream_fluid_bc

        # init right wall
        index = get_index(index_i=most_right_index, index_j=j, param=parameters, dim=1)
        field[index] = flags["outflow"]
        update_funcs[index] = update_density_bc_cell
        update_velocity_funcs[index] = update_velocity_bc_cell
        stream_funcs[index] = stream_fluid_bc

    bottom_index = 0
    top_index = parameters["height"] - 1

    for i in range(parameters["width"]):
        # init top (moving) wall
        index = get_index(index_i=i, index_j=top_index, param=parameters, dim=1)
        field[index] = flags["wall"]
        update_funcs[index] = update_density_bc_cell
        update_velocity_funcs[index] = update_velocity_bc_cell
        stream_funcs[index] = stream_fluid_bc

        # init bottom wall
        index = get_index(index_i=i, index_j=bottom_index, param=parameters, dim=1)
        field[index] = flags["wall"]
        update_funcs[index] = update_density_bc_cell
        update_velocity_funcs[index] = update_velocity_bc_cell
        stream_funcs[index] = stream_fluid_bc


    for j in range(parameters["obstacle"][2], parameters["obstacle"][3]):
        for i in range(parameters["obstacle"][0], parameters["obstacle"][1]):

            index = get_index(index_i=i, index_j=j, param=parameters, dim=1)
            field[index] = flags["wall"]
            update_funcs[index] = update_density_bc_cell
            update_velocity_funcs[index] = update_velocity_bc_cell
            stream_funcs[index] = stream_fluid_bc


def init_cavity_population_field(population, weights, parameters):

    num_components = parameters["discretization"]

    for j in range(parameters["height"]):
        for i in range(parameters["width"]):
            vector_population_index = get_index(i, j, parameters, dim=num_components)

            for component in range(num_components):
                population[vector_population_index + component] = weights[component]


def display_scalar_field(field, figure, axis, parameters, dim=1, shift=0):

    screen = np.zeros((parameters["height"], parameters["width"]))

    for j in range(parameters["height"]):
        for i in range(parameters["width"]):
            screen[j][i] = field[get_index(i, j, parameters, dim) + shift]

    img = axis.imshow(screen, origin='lower')


def display_vector_magnitude_2d(field, figure, axis, parameters):

    screen = np.zeros((parameters["height"], parameters["width"]))

    for j in range(parameters["height"]):
        for i in range(parameters["width"]):
            index = get_index(i, j, parameters, dim=2)
            screen[j][i] = np.sqrt(field[index] * field[index] + field[index + 1] * field[index + 1])

    axis.imshow(screen, origin='lower')


def display_vector_field_2d(field, figure, axis, parameters):
    axis.clear()

    X, Y = np.mgrid[0:parameters["height"], 0:parameters["width"]]

    U = np.zeros((parameters["height"], parameters["width"]))
    V = np.zeros((parameters["height"], parameters["width"]))
    R = np.zeros((parameters["height"], parameters["width"]))

    for j in range(parameters["height"]):
        for i in range(parameters["width"]):
            index = get_index(i, j, parameters, dim=2)
            R[j][i] = np.sqrt(field[index] * field[index] + field[index + 1] * field[index + 1])
            U[j][i] = field[index]
            V[j][i] = field[index + 1]

    axis.quiver(Y, X, U, V, R, alpha=.5)


def draw_obstacle(axis, parameters):
    x_origin = parameters["obstacle"][0]
    y_origin = parameters["obstacle"][2]
    width = parameters["obstacle"][1] - parameters["obstacle"][0]
    height = parameters["obstacle"][3] - parameters["obstacle"][2]
    axis.add_patch(
        patches.Rectangle( (x_origin, y_origin),  # (x,y)
                           width,  # width
                           height,  # height
                           fill=False
        )
    )

def print_dict(dictionary):
    print("parameters of the problem:")
    for key in dictionary.keys():
        print(key, ":", dictionary[key])

main()