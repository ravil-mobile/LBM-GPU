from helper import *

def update_density_field(density, population, update_functions, parameters):
    # call update functions according the boundary conditions
    for i in range(parameters["height"]):
        for j in range(parameters["width"]):
            scalar_index = get_index(i, j, parameters, dim=1)
            vector_population_index = parameters["discretization"] * scalar_index
            density[scalar_index] = update_functions[scalar_index](vector_population_index,
                                                                   population,
                                                                   parameters["discretization"])


def update_density_fluid_cell(index, population, num_directions):
    lattice_density = 0.0
    for i in range(index, index + num_directions):
        lattice_density += population[i]
    return lattice_density


def update_density_bc_cell(index, population, num_directions):
    return 1.0