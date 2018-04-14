import numpy as np


def get_index(index_i, index_j, param, dim):
    return dim * (index_i + index_j * param["width"])
