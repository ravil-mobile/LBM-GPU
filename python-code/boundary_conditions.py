from helper import *

def get_bc_function(flag):

    if flag == "wall":
        return apply_non_slip_bc
    elif flag == "mov_wall":
        return apply_mov_wall_bc
    else:
        return skip_boundary


def skip_boundary(component,
                  bc_conditions,
                  coordinate,
                  population,
                  velocity,
                  density,
                  coords,
                  inverse_indices,
                  weights,
                  constants,
                  parameters):
    pass


def apply_non_slip_bc(component,
                      bc_conditions,
                      coordinate,
                      population,
                      velocity,
                      density,
                      coords,
                      inverse_indices,
                      weights,
                      constants,
                      parameters):

    # Compute inverse coordinate
    i = -1 * coords[parameters["dimension"] * component]
    j = -1 * coords[parameters["dimension"] * component + 1]

    # Find inverse index
    inverse_component = inverse_indices[component]

    # compute the shift to the target fluid lattice
    shift = get_index(i, j, parameters, dim=parameters["discretization"])

    self_index = coordinate + component
    neightbour_index = coordinate + shift + inverse_component

    population[neightbour_index] = population[self_index]


def apply_mov_wall_bc(component,
                      bc_conditions,
                      coordinate,
                      population,
                      velocity,
                      density,
                      coords,
                      inverse_indices,
                      weights,
                      constants,
                      parameters):


    # Compute inverse coordinate
    i = -1 * coords[parameters["dimension"] * component]
    j = -1 * coords[parameters["dimension"] * component + 1]

    # Find inverse index
    inverse_component = inverse_indices[component]

    # compute the dot product Cinv.T * Uwall
    dot_product_CU = constants["const_1"] \
                   * (i * parameters["wall_velocity"][0] + j * parameters["wall_velocity"][1])

    # compute the shift to the target fluid lattice
    shift = get_index(i, j, parameters, dim=parameters["discretization"])

    self_index = coordinate + component
    neightbour_index = coordinate + shift + inverse_component
    scalar_sneightbour_index = (coordinate + shift) // parameters["discretization"]

    population[neightbour_index] = population[self_index] \
                                 + 2.0 * weights[component] * density[scalar_sneightbour_index] * dot_product_CU


def apply_inflow_bc(component,
                    bc_conditions,
                    coordinate,
                    population,
                    velocity,
                    density,
                    coords,
                    inverse_indices,
                    weights,
                    constants,
                    parameters):

    # Compute inverse coordinate
    i = -1 * coords[parameters["dimension"] * component]
    j = -1 * coords[parameters["dimension"] * component + 1]

    # Find inverse index
    inverse_component = inverse_indices[component]

    # compute the shift to the target fluid lattice
    shift = get_index(i, j, parameters, dim=parameters["discretization"])

    self_index = coordinate + component
    neightbour_index = coordinate + shift + inverse_component


    dot_product_UU = parameters["velocity_inflow"][0] * parameters["velocity_inflow"][0] \
                   + parameters["velocity_inflow"][1] * parameters["velocity_inflow"][1]

    dot_product_CU = i * parameters["velocity_inflow"][0] + j * parameters["velocity_inflow"][1]


    velocity_expansion = constants["const_1"] * dot_product_CU \
                         + constants["const_2"] * dot_product_CU * dot_product_CU \
                         - constants["const_3"] * dot_product_UU \
                         + 1.0

    # Compute equilibrium distribution
    equilibrium = weights[inverse_component] * parameters["density_inflow"] * velocity_expansion


    population[neightbour_index] = equilibrium


def apply_outflow_bc(component,
                     bc_conditions,
                     coordinate,
                     population,
                     velocity,
                     density,
                     coords,
                     inverse_indices,
                     weights,
                     constants,
                     parameters):

    # Compute coordinate
    i = coords[parameters["dimension"] * component]
    j = coords[parameters["dimension"] * component + 1]

    # Compute inverse coordinate
    i_inv = -1 * coords[parameters["dimension"] * component]
    j_inv = -1 * coords[parameters["dimension"] * component + 1]

    # Find inverse index
    inverse_component = inverse_indices[component]

    # compute the shift to the target fluid lattice
    shift = get_index(i_inv, j_inv, parameters, dim=parameters["discretization"])

    self_index = coordinate + component
    neightbour_index = coordinate + shift + inverse_component

    neightbour_vector_index = ((coordinate  + shift ) // parameters["discretization"]) * parameters["dimension"]
    neightbour_velocity_x = velocity[neightbour_vector_index]
    neightbour_velocity_y = velocity[neightbour_vector_index + 1]

    dot_product_UU = neightbour_velocity_x * neightbour_velocity_x \
                   + neightbour_velocity_y * neightbour_velocity_y

    dot_product_CU = i * neightbour_velocity_x + j * neightbour_velocity_y
    dot_product_CU_inv = i_inv * neightbour_velocity_x + j_inv * neightbour_velocity_y

    velocity_expansion = constants["const_1"] * dot_product_CU \
                       + constants["const_2"] * dot_product_CU * dot_product_CU \
                       - constants["const_3"] * dot_product_UU \
                       + 1.0

    velocity_expansion_inv = constants["const_1"] * dot_product_CU_inv \
                           + constants["const_2"] * dot_product_CU_inv * dot_product_CU_inv \
                           - constants["const_3"] * dot_product_UU \
                           + 1.0


    # Compute equilibrium distribution
    equilibrium = weights[component] * parameters["density_outflow"] * velocity_expansion
    equilibrium_inv = weights[inverse_component] * parameters["density_outflow"] * velocity_expansion_inv

    population[neightbour_index] = equilibrium + equilibrium_inv - population[self_index]
