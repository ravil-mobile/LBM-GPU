from helper import *
from parameters import SimulationParametes as parameters
from parameters import Constants as constants


def update_velocity_field(velocity, population, density, update_functions):

    num_lattices = parameters.num_lattices
    height = parameters.height
    width = parameters.width

    # call update functions according the boundary conditions
    for j in range(height):
        for i in range(width):

            # compute index for the scalar field
            scalar_index = get_index(i, j, dim=1)
            inverse_lattice_density = 1.0 / density[scalar_index]

            lattice_velocity = update_functions[scalar_index](scalar_index,
                                                              population)

            velocity[scalar_index] = inverse_lattice_density * lattice_velocity[0]
            velocity[scalar_index + num_lattices] = inverse_lattice_density * lattice_velocity[1]


def update_velocity_fluid_cell(index, population):

    lattice_velocity = [0.0, 0.0]
    num_directions = parameters.discretization
    num_lattices = parameters.num_lattices

    for component in range(num_directions):

        lattice_velocity[0] += (constants.coords[component] *
                                population[index + component * num_lattices])

        lattice_velocity[1] += (constants.coords[component + num_directions] *
                                population[index + component * num_lattices])

    return lattice_velocity


def update_velocity_bc_cell(index, population):
    return [0.0, 0.0]
