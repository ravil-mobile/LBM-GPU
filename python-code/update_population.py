from helper import *


def update_population_field(velocity, population, density, coords, weights, constants, parameters):

    dim = parameters["dimension"]

    # call update functions according the boundary conditions
    for i in range(parameters["height"]):
        for j in range(parameters["width"]):
            scalar_index = get_index(i, j, parameters, dim=1)
            vector_population_index = parameters["discretization"] * scalar_index
            vector_velocity_index = dim * scalar_index

            # Compute dot product - u.T * u
            dot_product_UU = 0.0
            for k in range(vector_velocity_index, vector_velocity_index + parameters["dimension"]):
                dot_product_UU += velocity[k] * velocity[k]

            # Compute update of population for each direction within the lattice
            for component, k in enumerate(range(vector_population_index,
                                                vector_population_index + parameters["discretization"])):

                # Compute velocity expantion for a particular direction within the lattice
                dot_product_CU = 0.0
                for vector_component, index in enumerate(range(vector_velocity_index,
                                                               vector_velocity_index + parameters["dimension"])):

                    dot_product_CU += coords[dim * component + vector_component] * velocity[index]


                velocity_expansion = constants["const_1"] * dot_product_CU \
                                   + constants["const_2"] * dot_product_CU * dot_product_CU \
                                   - constants["const_3"] * dot_product_UU \
                                   + 1.0

                # Compute equilibrium distribution
                equilibrium = weights[component] * density[scalar_index] * velocity_expansion

                # Perform update
                population[k] -= parameters["relaxation"] * (population[k] - equilibrium)
