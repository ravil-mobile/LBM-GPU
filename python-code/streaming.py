from helper import *


def stream(population, swap_buffer, coords, stream_funcs, parameters):

    height = parameters["height"]
    width = parameters["width"]

    # call update functions according the boundary conditions
    for j in range(height):
        for i in range(width):
            scalar_index = get_index(i, j, parameters, dim=1)
            stream_funcs[scalar_index](i, j, population, swap_buffer, coords, parameters)


def stream_fluid_cell(i, j, population, swap_buffer, coords, parameters):

    dim = parameters["dimension"]
    num_directions = parameters["discretization"]

    vector_self_index = get_index(i, j, parameters, dim=num_directions)

    for component in range(num_directions):
        ii = coords[dim * component]
        jj = coords[dim * component + 1]
        vector_neighbour_index = get_index(i + ii, j + jj, parameters, dim=num_directions)
        swap_buffer[vector_neighbour_index + component] = population[vector_self_index + component]


def stream_fluid_bc(i, j, population, swap_buffer, coords, parameters):
    pass
