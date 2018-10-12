#!/usr/bin/env python3
# coding: utf-8

import numpy as np
from scipy import sparse

DEFAULT_BINNING = 1000


def tsv2csv(tsv, binning=DEFAULT_BINNING, output=None):

    matrix = np.genfromtxt(tsv, dtype=None, comments=None, delimiter="\t")
    (_, _, pos1, *_, pos2, _, _, _, _) = zip(*matrix)

    positions1 = np.array(pos1) // binning
    positions2 = np.array(pos2) // binning

    minimum = min(np.amin(positions1), np.amin(positions2))
    positions1 -= minimum
    positions2 -= minimum

    n = int(max(np.amax(positions1), np.amax(positions2))) + 1
    try:
        assert len(positions1) == len(positions2)
    except AssertionError:
        print(len(positions1), len(positions2))
        raise
    sparse_matrix = sparse.coo_matrix(
        (np.ones(len(positions1)), (positions1, positions2)), shape=(n, n)
    )

    dense_matrix = np.array(sparse_matrix.todense(), dtype=np.int32)
    if output is not None:
        np.savetxt(output, dense_matrix, fmt="%i")

    return dense_matrix
