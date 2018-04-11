from helper import *

def update_velocity_field(velocity, population, density, coords, update_functions, parameters):
    # call update functions according the boundary conditions
    for i in range(parameters["height"]):
        for j in range(parameters["width"]):
            # compute index for the scalar field
            scalar_index = get_index(i, j, parameters, dim=1)
            inverse_lattice_density = 1.0 / density[scalar_index]

            vector_population_index = parameters["discretization"] * scalar_index

            lattice_velocity = update_functions[scalar_index](vector_population_index,
                                                              population,
                                                              coords,
                                                              parameters["dimension"],
                                                              parameters["discretization"])

            vector_velocity_index = parameters["dimension"] * scalar_index

            for component, k in enumerate(range(vector_velocity_index,
                                                vector_velocity_index + parameters["dimension"])):
                velocity[k] = inverse_lattice_density * lattice_velocity[component]


def update_velocity_fluid_cell(index, population, coords, dim, num_directions):

    lattice_velocity = [0.0] * dim

    for direction, i in enumerate(range(index, index + num_directions)):
        for component in range(dim):
            lattice_velocity[component] += population[i] * coords[dim*direction + component]

    return lattice_velocity


def update_velocity_bc_cell(index, population, coords, dim, num_directions):
    return [0.0] * dim