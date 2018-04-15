from helper import *
from parameters import SimulationParametes as parameters
from parameters import Constants as constants
from parameters import BoundaryInfo as boundary_info


def get_bc_function(flag):

    if flag == "wall":
        return apply_non_slip_bc
    elif flag == "mov_wall":
        return apply_mov_wall_bc
    else:
        return skip_boundary


def skip_boundary(component,
                  coordinate,
                  population,
                  velocity,
                  density):
    pass


def apply_non_slip_bc(component,
                      coordinate,
                      population,
                      velocity,
                      density):

    num_lattices = parameters.num_lattices
    num_directions = parameters.discretization

    # Compute inverse coordinate
    i = -1 * constants.coords[component]
    j = -1 * constants.coords[component + num_directions]

    # Find inverse index
    inverse_component = constants.inverse_indices[component]

    # compute the gap to the target fluid lattice
    shift = get_index(i, j, dim=1)

    self_index = (coordinate) + component * num_lattices
    neightbour_index = (coordinate + shift) + inverse_component * num_lattices

    population[neightbour_index] = population[self_index]


def apply_mov_wall_bc(component,
                      coordinate,
                      population,
                      velocity,
                      density):

    num_directions = parameters.discretization
    num_lattices = parameters.num_lattices

    # Compute inverse coordinate
    i = -1 * constants.coords[component]
    j = -1 * constants.coords[component + num_directions]

    # Find inverse index
    inverse_component = constants.inverse_indices[component]

    # compute the dot product Cinv.T * Uwall
    dot_product_CU = constants.const_1 * (i * boundary_info.wall_velocity[0] +
                                          j * boundary_info.wall_velocity[1])

    # compute the shift to the target fluid lattice
    shift = get_index(i, j, dim=1)

    self_index = (coordinate) + component * num_lattices
    neightbour_index = (coordinate + shift) + inverse_component * num_lattices

    scalar_neighbour_index = coordinate + shift

    population[neightbour_index] = (population[self_index] +
                                    (2.0 * constants.weights[component] *
                                     density[scalar_neighbour_index] *
                                     dot_product_CU))

def apply_inflow_bc(component,
                    coordinate,
                    population,
                    velocity,
                    density):

    num_directions = parameters.discretization
    num_lattices = parameters.num_lattices

    # Compute inverse coordinate
    i = -1 * constants.coords[component]
    j = -1 * constants.coords[component + num_directions]

    # Find inverse index
    inverse_component = constants.inverse_indices[component]

    # compute the shift to the target fluid lattice
    shift = get_index(i, j, dim=1)

    self_index = (coordinate) + component * num_lattices
    neightbour_index = (coordinate + shift) + inverse_component * num_lattices

    dot_product_UU = (boundary_info.velocity_inflow[0] * boundary_info.velocity_inflow[0] +
                      boundary_info.velocity_inflow[1] * boundary_info.velocity_inflow[1])

    dot_product_CU = i * boundary_info.velocity_inflow[0] + j * boundary_info.velocity_inflow[1]


    velocity_expansion = (constants.const_1 * dot_product_CU +
                          constants.const_2 * dot_product_CU * dot_product_CU -
                          constants.const_3 * dot_product_UU +
                          1.0)

    # Compute equilibrium distribution
    equilibrium = constants.weights[inverse_component] * boundary_info.density_inflow * velocity_expansion

    population[neightbour_index] = equilibrium


def apply_outflow_bc(component,
                     coordinate,
                     population,
                     velocity,
                     density):


    num_directions = parameters.discretization
    num_lattices = parameters.num_lattices

    # Compute inverse coordinate
    i = 1 * constants.coords[component]
    j = 1 * constants.coords[component + num_directions]

    i_inv = -i
    j_inv = -j

    # Find inverse index
    inverse_component = constants.inverse_indices[component]

    # compute the shift to the target fluid lattice
    shift = get_index(i_inv, j_inv, dim=1)

    vector_self_index = (coordinate) + component * num_lattices
    vector_neightbour_index = (coordinate + shift) + inverse_component * num_lattices
    scalar_neightbour_index = coordinate + shift

    neightbour_velocity_x = velocity[scalar_neightbour_index]
    neightbour_velocity_y = velocity[scalar_neightbour_index + num_lattices]

    dot_product_UU = (neightbour_velocity_x * neightbour_velocity_x +
                      neightbour_velocity_y * neightbour_velocity_y)

    dot_product_CU = i * neightbour_velocity_x + j * neightbour_velocity_y
    dot_product_CU_inv = i_inv * neightbour_velocity_x + j_inv * neightbour_velocity_y


    velocity_expansion = (constants.const_1 * dot_product_CU +
                          constants.const_2 * dot_product_CU * dot_product_CU -
                          constants.const_3 * dot_product_UU +
                          1.0)

    velocity_expansion_inv = (constants.const_1 * dot_product_CU_inv +
                              constants.const_2 * dot_product_CU_inv * dot_product_CU_inv -
                              constants.const_3 * dot_product_UU +
                              1.0)


    # Compute equilibrium distribution
    equilibrium = (constants.weights[component] *
                   boundary_info.density_outflow *
                   velocity_expansion)

    equilibrium_inv = (constants.weights[inverse_component] *
                       boundary_info.density_outflow *
                       velocity_expansion_inv)

    population[vector_neightbour_index] = (equilibrium +
                                           equilibrium_inv -
                                           population[vector_self_index])
