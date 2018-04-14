from helper import *
from parameters import SimulationParametes as parameters


def update_density_field(density, population, update_functions):

    num_directions = parameters.discretization

    height = parameters.height
    width = parameters.width

    # call update functions according the boundary conditions
    for j in range(height):
        for i in range(width):
            scalar_index = get_index(i, j, dim=1)

            vector_population_index = num_directions * scalar_index

            density[scalar_index] = update_functions[scalar_index](vector_population_index,
                                                                   population,
                                                                   num_directions)


def update_density_fluid_cell(index, population, num_directions):
    lattice_density = 0.0
    for i in range(index, index + num_directions):
        lattice_density += population[i]
    return lattice_density


def update_density_bc_cell(index, population, num_directions):
    return 1.0
