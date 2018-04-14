from helper import *
from parameters import SimulationParametes as parameters
from parameters import Constants as constants

def stream(population, swap_buffer, stream_funcs):

    height = parameters.height
    width = parameters.width

    # call update functions according the boundary conditions
    for j in range(height):
        for i in range(width):
            scalar_index = get_index(i, j, dim=1)
            stream_funcs[scalar_index](i, j, population, swap_buffer)


def stream_fluid_cell(i, j, population, swap_buffer):

    dim = parameters.dimension
    num_directions = parameters.discretization

    vector_self_index = get_index(i, j, dim=num_directions)

    for component in range(num_directions):
        ii = constants.coords[dim * component]
        jj = constants.coords[dim * component + 1]
        vector_neighbour_index = get_index(i + ii, j + jj, dim=num_directions)
        swap_buffer[vector_neighbour_index + component] = population[vector_self_index + component]


def stream_fluid_bc(i, j, population, swap_buffer):
    pass
