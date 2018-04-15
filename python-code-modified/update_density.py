from helper import *
from parameters import SimulationParametes as parameters


def update_density_field(density, population, update_functions):

    height = parameters.height
    width = parameters.width

    # call update functions according the boundary conditions
    for j in range(height):
        for i in range(width):
            scalar_index = get_index(i, j, dim=1)

            density[scalar_index] = update_functions[scalar_index](scalar_index,
                                                                   population)


def update_density_fluid_cell(index, population):

    lattice_density = 0.0
    num_lattices = parameters.num_lattices

    for component in range(parameters.discretization):

        lattice_density += population[index + component * num_lattices]

    return lattice_density


def update_density_bc_cell(index, population):
    return 1.0
