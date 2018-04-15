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

    num_directions = parameters.discretization
    num_lattices = parameters.num_lattices

    self_index = get_index(i, j, dim=1)

    for component in range(num_directions):

        ii = constants.coords[component]
        jj = constants.coords[component + num_directions]

        neighbour_index = get_index(i + ii, j + jj, dim=1)

        shift = component * num_lattices
        swap_buffer[neighbour_index + shift] = population[self_index + shift]


def stream_fluid_bc(i, j, population, swap_buffer):
    pass
