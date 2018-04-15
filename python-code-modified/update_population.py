from helper import *
from parameters import SimulationParametes as parameters
from parameters import Constants as constants


def update_population_field(velocity, population, density):

    num_directions = parameters.discretization
    num_lattices = parameters.num_lattices
    relaxation = parameters.relaxation

    height = parameters.height
    width = parameters.width

    # call update functions according the boundary conditions
    for j in range(height):
        for i in range(width):

            scalar_index = get_index(i, j, dim=1)
            local_velocity_x = velocity[scalar_index]
            local_velocity_y = velocity[scalar_index + num_lattices]

            # Compute dot product - u.T * u
            dot_product_UU = (local_velocity_x * local_velocity_x +
                              local_velocity_y * local_velocity_y)

            # Compute update of population for each direction within the lattice
            for component in range(num_directions):

                # Compute coordinate
                vector_component_x = constants.coords[component]
                vector_component_y = constants.coords[component + num_directions]

                # Compute velocity expantion for a particular direction within the lattice
                dot_product_CU = (vector_component_x * local_velocity_x +
                                  vector_component_y * local_velocity_y)

                velocity_expansion = (constants.const_1 * dot_product_CU +
                                      constants.const_2 * dot_product_CU * dot_product_CU -
                                      constants.const_3 * dot_product_UU +
                                      1.0)

                # Compute equilibrium distribution
                equilibrium = constants.weights[component] * density[scalar_index] * velocity_expansion

                shift = component * num_lattices

                # Perform update
                population[scalar_index + shift] -= relaxation * (population[scalar_index + shift] - equilibrium)
