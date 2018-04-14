from helper import *

def update_velocity_field(velocity, population, density, coords, update_functions, parameters):

    dim = parameters["dimension"]
    num_directions = parameters["discretization"]

    height = parameters["height"]
    width = parameters["width"]

    # call update functions according the boundary conditions
    for j in range(height):
        for i in range(width):

            # compute index for the scalar field
            scalar_index = get_index(i, j, parameters, dim=1)
            inverse_lattice_density = 1.0 / density[scalar_index]

            vector_population_index = num_directions * scalar_index

            lattice_velocity = update_functions[scalar_index](vector_population_index,
                                                              population,
                                                              coords,
                                                              dim,
                                                              num_directions)

            vector_velocity_index = dim * scalar_index

            velocity[vector_velocity_index] = inverse_lattice_density * lattice_velocity[0]
            velocity[vector_velocity_index + 1] = inverse_lattice_density * lattice_velocity[1]


def update_velocity_fluid_cell(index, population, coords, dim, num_directions):

    lattice_velocity = [0.0, 0.0]

    for direction, i in enumerate(range(index, index + num_directions)):

        lattice_velocity[0] += coords[dim * direction] * population[i]
        lattice_velocity[1] += coords[dim * direction + 1] * population[i]

    return lattice_velocity


def update_velocity_bc_cell(index, population, coords, dim, num_directions):
    return [0.0, 0.0]