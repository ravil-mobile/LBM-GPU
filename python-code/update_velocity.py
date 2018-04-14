from helper import *
from parameters import SimulationParametes as parameters
from parameters import Constants as constants


def update_velocity_field(velocity, population, density, update_functions):

    dim = parameters.dimension
    num_directions = parameters.discretization

    height = parameters.height
    width = parameters.width

    # call update functions according the boundary conditions
    for j in range(height):
        for i in range(width):

            # compute index for the scalar field
            scalar_index = get_index(i, j, dim=1)
            inverse_lattice_density = 1.0 / density[scalar_index]

            vector_population_index = num_directions * scalar_index

            lattice_velocity = update_functions[scalar_index](vector_population_index,
                                                              population)

            vector_velocity_index = dim * scalar_index

            velocity[vector_velocity_index] = inverse_lattice_density * lattice_velocity[0]
            velocity[vector_velocity_index + 1] = inverse_lattice_density * lattice_velocity[1]


def update_velocity_fluid_cell(index, population):

    lattice_velocity = [0.0, 0.0]

    for direction, i in enumerate(range(index, index + parameters.discretization)):

        lattice_velocity[0] += constants.coords[parameters.dimension * direction] * population[i]
        lattice_velocity[1] += constants.coords[parameters.dimension * direction + 1] * population[i]

    return lattice_velocity


def update_velocity_bc_cell(index, population):
    return [0.0, 0.0]
