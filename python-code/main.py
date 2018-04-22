import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import time

from parameters import SimulationParametes as parameters
from parameters import BoundaryInfo as boundary_info
from parameters import Constants as constants

from update_density import *
from update_velocity import *
from update_population import *
from streaming import *
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
    parameters.simulation_time = 1.0
    parameters.num_time_steps = 100000
    parameters.dimension = 2
    parameters.discretization = 9

    parameters.delta_x = 0.577 * 1e-3
    parameters.delta_t = 1e-3
    parameters.speed_of_sound = parameters.delta_x / parameters.delta_t
    parameters.viscosity = 85.06e-6
    parameters.tau = 0.5 * (1 + (6.0 * parameters.viscosity * parameters.delta_t) / (parameters.delta_x**2))
    parameters.relaxation = 1.0 / parameters.tau

    parameters.width = 30
    parameters.height = 30
    parameters.num_lattices = parameters.width * parameters.height

    #parameters.obstacle = [10, 15, 10, 15]
    parameters.obstacle = [0, 0, 0, 0]

    boundary_info.wall_velocity = [0.05, 0.0]

    boundary_info.density_inflow = 1.0
    boundary_info.velocity_inflow = [0.26, 0.0]
    boundary_info.inflow_direction = 1

    boundary_info.density_outflow = 1.0
    boundary_info.velocity_outflow = [0.26, 0.0]
    boundary_info.outflow_direction = 1

    constants.const_1 = 1.0 / (parameters.speed_of_sound**2)
    constants.const_2 = 0.5 * (constants.const_1**2)
    constants.const_3 = 0.5 * (constants.const_1)

    # arrays in c-style
    flag_field = np.zeros(parameters.num_lattices, dtype=np.uint8)
    density = np.ones(parameters.num_lattices, dtype=np.float32)
    velocity = np.zeros(parameters.dimension * parameters.num_lattices, dtype=np.float32)
    population = np.ones(parameters.discretization * parameters.num_lattices, dtype=np.float32)
    swap_buffer = np.ones(parameters.discretization * parameters.num_lattices, dtype=np.float32)

    # list of functions
    update_density_funcs = [update_density_fluid_cell] * parameters.num_lattices
    update_velocity_funcs = [update_velocity_fluid_cell] * parameters.num_lattices
    stream_funcs = [stream_fluid_cell] * parameters.num_lattices

    # initialization of all fields
    init_cavity_flag_field(flag_field,
                           update_density_funcs,
                           update_velocity_funcs,
                           stream_funcs)

    init_cavity_population_field(population, swap_buffer)

    boundaries, boundary_coordinates = scanning(flag_field)

    fig = plt.figure(1)
    fig.clf()
    ax = fig.add_subplot(1, 1, 1)

    stream_time = 0.0
    treat_boundaries_time = 0.0
    update_density_time = 0.0
    update_velocity_time = 0.0
    update_population_time = 0.0

    # run simulation
    start = time.time()
    for step in range(parameters.num_time_steps):

        # perform streaming
        start_function_call = time.time()
        stream(population, swap_buffer, stream_funcs)
        stream_time += (time.time() - start_function_call)

        start_function_call = time.time()
        treat_boundaries(boundaries,
                         boundary_coordinates,
                         swap_buffer,
                         velocity,
                         density)
        treat_boundaries_time += (time.time() - start_function_call)

        # swap buffers
        temp = population
        population = swap_buffer
        swap_buffer = temp

        # computer collision
        start_function_call = time.time()
        update_density_field(density, population, update_density_funcs)
        update_density_time += (time.time() - start_function_call)

        start_function_call = time.time()
        update_velocity_field(velocity, population, density, update_velocity_funcs)
        update_velocity_time += (time.time() - start_function_call)

        start_function_call = time.time()
        update_population_field(velocity, population, density)
        update_population_time += (time.time() - start_function_call)

        if (step % 250) == 0 :
            # display_scalar_field(field=density, figure=fig, axis=ax, dim=1, shift=0)
            # display_vector_magnitude_2d(field=velocity, figure=fig, axis=ax)
            display_vector_field_2d(field=velocity, figure=fig, axis=ax)
            draw_obstacle(ax)
            plt.pause(0.01)


            print("iteration step: %i; density: max = %f; min = %f" % (step, np.max(density), np.min(density)))

    elapsed_time = time.time() - start
    MLUPS = (parameters.num_lattices * parameters.num_time_steps) / (elapsed_time * 1e6)

    delimiter = "=" * 80
    print(delimiter)
    print("MLUPS: ", MLUPS)
    print("elapsed time", elapsed_time)

    print(delimiter)
    print("stream_time: %f percent" % (100.0 * (stream_time / elapsed_time)))
    print("treat_boundaries_time: %f percent" % (100.0 * (treat_boundaries_time / elapsed_time)))
    print("update_density_time: %f percent" % (100.0 * (update_density_time / elapsed_time)))
    print("update_velocity_time: %f percent" % (100.0 * (update_velocity_time / elapsed_time)))
    print("update_population_time: %f percent" % (100.0 * (update_population_time / elapsed_time)))


def treat_boundaries(boundaries,
                     boundary_coordinates,
                     population,
                     velocity,
                     density):

    num_components = parameters.discretization

    for counter, boundary_function in enumerate(boundaries):
        component = counter % num_components
        scalar_counter = counter // num_components

        boundary_function(component,
                          boundary_coordinates[scalar_counter],
                          population,
                          velocity,
                          density)


def init_cavity_flag_field(field, update_funcs, update_velocity_funcs, stream_funcs):

    most_left_index = 0
    most_right_index = parameters.width - 1

    for j in range(parameters.height):

        # init left wall
        index = get_index(index_i=most_left_index, index_j=j, dim=1)
        field[index] = constants.flags["wall"]
        update_funcs[index] = update_density_bc_cell
        update_velocity_funcs[index] = update_velocity_bc_cell
        stream_funcs[index] = stream_fluid_bc

        # init right wall
        index = get_index(index_i=most_right_index, index_j=j, dim=1)
        field[index] = constants.flags["wall"]
        update_funcs[index] = update_density_bc_cell
        update_velocity_funcs[index] = update_velocity_bc_cell
        stream_funcs[index] = stream_fluid_bc

    bottom_index = 0
    top_index = parameters.height - 1

    for i in range(parameters.width):
        # init top (moving) wall
        index = get_index(index_i=i, index_j=top_index, dim=1)
        field[index] = constants.flags["mov_wall"]
        update_funcs[index] = update_density_bc_cell
        update_velocity_funcs[index] = update_velocity_bc_cell
        stream_funcs[index] = stream_fluid_bc

        # init bottom wall
        index = get_index(index_i=i, index_j=bottom_index, dim=1)
        field[index] = constants.flags["wall"]
        update_funcs[index] = update_density_bc_cell
        update_velocity_funcs[index] = update_velocity_bc_cell
        stream_funcs[index] = stream_fluid_bc

    for j in range(parameters.obstacle[2], parameters.obstacle[3]):
        for i in range(parameters.obstacle[0], parameters.obstacle[1]):

            index = get_index(index_i=i, index_j=j, dim=1)
            field[index] = constants.flags["wall"]
            update_funcs[index] = update_density_bc_cell
            update_velocity_funcs[index] = update_velocity_bc_cell
            stream_funcs[index] = stream_fluid_bc


def init_cavity_population_field(population, swap_buffer):

    num_components = parameters.discretization

    for j in range(parameters.height):
        for i in range(parameters.width):
            vector_population_index = get_index(i, j, dim=num_components)

            for component in range(num_components):
                population[vector_population_index + component] = constants.weights[component]
                swap_buffer[vector_population_index + component] = constants.weights[component]


def display_scalar_field(field, figure, axis, dim=1, shift=0):

    screen = np.zeros((parameters.height, parameters.width))

    for j in range(parameters.height):
        for i in range(parameters.width):
            screen[j][i] = field[get_index(i, j, dim) + shift]

    axis.imshow(screen, origin='lower')


def display_vector_magnitude_2d(field, figure, axis,):

    screen = np.zeros((parameters.height, parameters.width))

    for j in range(parameters.height):
        for i in range(parameters.width):
            index = get_index(i, j, dim=2)
            screen[j][i] = np.sqrt(field[index] * field[index] + field[index + 1] * field[index + 1])

    axis.imshow(screen, origin='lower')


def display_vector_field_2d(field, figure, axis):
    axis.clear()

    X, Y = np.mgrid[0:parameters.height, 0:parameters.width]

    U = np.zeros((parameters.height, parameters.width))
    V = np.zeros((parameters.height, parameters.width))
    R = np.zeros((parameters.height, parameters.width))

    for j in range(parameters.height):
        for i in range(parameters.width):
            index = get_index(i, j, dim=2)
            R[j][i] = np.sqrt(field[index] * field[index] + field[index + 1] * field[index + 1])
            U[j][i] = field[index]
            V[j][i] = field[index + 1]

    axis.quiver(Y, X, U, V, R, alpha=.5)


def draw_obstacle(axis):
    x_origin = parameters.obstacle[0]
    y_origin = parameters.obstacle[2]
    width = parameters.obstacle[1] - parameters.obstacle[0]
    height = parameters.obstacle[3] - parameters.obstacle[2]
    axis.add_patch(
        patches.Rectangle((x_origin, y_origin),
                          width,
                          height,
                          fill=False))


def print_dict(dictionary):
    print("parameters of the problem:")
    for key in dictionary.keys():
        print(key, ":", dictionary[key])


main()
