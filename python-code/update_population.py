from helper import *
from parameters import SimulationParametes as parameters
from parameters import Constants as constants


def update_population_field(velocity, population, density):

    dim = parameters.dimension
    num_directions = parameters.discretization
    relaxation = parameters.relaxation

    height = parameters.height
    width = parameters.width

    # call update functions according the boundary conditions
    for j in range(height):
        for i in range(width):
            scalar_index = get_index(i, j, dim=1)
            vector_population_index = num_directions * scalar_index
            vector_velocity_index = dim * scalar_index

            # Compute dot product - u.T * u
            dot_product_UU = 0.0
            for k in range(vector_velocity_index, vector_velocity_index + dim):
                dot_product_UU += velocity[k] * velocity[k]

            # Compute update of population for each direction within the lattice
            for component, k in enumerate(range(vector_population_index,
                                                vector_population_index + num_directions)):

                # Compute coordinate
                vector_component_x = constants.coords[dim * component]
                vector_component_y = constants.coords[dim * component + 1]

                # Compute velocity expantion for a particular direction within the lattice
                dot_product_CU = (vector_component_x * velocity[vector_velocity_index] +
                                  vector_component_y * velocity[vector_velocity_index + 1])

                velocity_expansion = (constants.const_1 * dot_product_CU +
                                      constants.const_2 * dot_product_CU * dot_product_CU -
                                      constants.const_3 * dot_product_UU +
                                      1.0)

                # Compute equilibrium distribution
                equilibrium = constants.weights[component] * density[scalar_index] * velocity_expansion

                # Perform update
                population[k] -= relaxation * (population[k] - equilibrium)
