from helper import *
from boundary_conditions import *


def scanning(flag_field, flags, coords, parameters):

    fluid_flag = flags["fluid"]
    dim = parameters["dimension"]

    num_directions = parameters["discretization"]

    height = parameters["height"]
    width = parameters["width"]

    container = {}

    # go through the domain and examine each lattice
    for j in range(height):
        for i in range(width):

            scalar_self_index = get_index(i, j, parameters, dim=1)

            # if the lattice is fluid, examine its surroundings
            if flag_field[scalar_self_index] == fluid_flag:

                for component in range(num_directions):
                    ii = coords[dim * component]
                    jj = coords[dim * component + 1]

                    scalar_neighbour_index = get_index(i + ii, j + jj, parameters, dim=1)
                    neighbour_flag = flag_field[scalar_neighbour_index]

                    # if the lattice neighbour is not fluid, add it
                    # to the container together with their relative positions
                    if neighbour_flag != fluid_flag:

                        function = skip_boundary

                        # Choose appropriate function according to the boundary conditions
                        if neighbour_flag == flags["wall"]:
                            function = apply_non_slip_bc

                        elif neighbour_flag == flags["mov_wall"]:
                            function = apply_mov_wall_bc

                        elif neighbour_flag == flags["inflow"]:
                            function = apply_inflow_bc

                        elif neighbour_flag == flags["outflow"]:
                            function = apply_outflow_bc


                        boundary_lattice = (component, function)

                        # add boundary lattice to the set
                        if scalar_neighbour_index in container:
                            container[scalar_neighbour_index].append(boundary_lattice)
                        else:
                            container[scalar_neighbour_index] = [boundary_lattice]


    num_bloundaries = len(container.keys())
    boundary_coordinates = [0] * num_bloundaries
    boundaries = [skip_boundary] * num_bloundaries * num_directions

    # Go through the entire container with boundary ellements
    for counter, neighbour_index in enumerate(container.keys()):
        boundary_index = num_directions * counter

        boundary_coordinates[counter] = num_directions * neighbour_index
        for i in range(len(container[neighbour_index])):

            component, function = container[neighbour_index][i]
            boundaries[boundary_index + component] = function

    return boundaries, boundary_coordinates