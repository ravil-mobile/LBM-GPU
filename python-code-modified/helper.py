from parameters import SimulationParametes as parameters


def get_index(index_i, index_j, dim):
    return dim * (index_i + index_j * parameters.width)
