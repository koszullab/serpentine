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
    assert len(positions1) == len(
        positions2
    ), "Mismatch between lengths {} and {}".format(
        len(positions1), len(positions2)
    )
    sparse_matrix = sparse.coo_matrix(
        (np.ones(len(positions1)), (positions1, positions2)), shape=(n, n)
    )

    dense_matrix = np.array(sparse_matrix.todense(), dtype=np.int32)
    if output is not None:
        np.savetxt(output, dense_matrix, fmt="%i")

    return dense_matrix


def misha2csv(misha=None, binning=DEFAULT_BINNING, output=None):

    from rpy2.robjects import r

    r_library_expression = """
    library("shaman");
    library("misha")
    """

    if misha is None:
        r_import_expression = """
        gsetroot(shaman_get_test_track_db());
        contact_map <- gextract("hic_obs", gintervals.2d(2, 175e06, 
        178e06, 2, 175e06, 178e06), colnames="score")
        """
    else:
        r.assign("path", misha)
        r_import_expression = """
        contact_map <- gextract("hic_obs", gintervals.2d(2, 0, 
        178e06, 2, 175e06, 178e06), colnames="score")
        """

    r(r_library_expression)
    r(r_import_expression)
    r("write.table(contact_map, 'exported_map.csv')")

    matrix = np.genfromtxt("exported_map.csv", dtype=None, skip_header=True)

    (_, _, start1, end1, _, start2, end2, contacts, _) = zip(*matrix)

    pos1 = (np.array(start1) + np.array(end1)) // 2
    pos2 = (np.array(start2) + np.array(end2)) // 2

    positions1 = np.array(pos1) // binning
    positions2 = np.array(pos2) // binning

    minimum = min(np.amin(positions1), np.amin(positions2))
    positions1 -= minimum
    positions2 -= minimum

    n = int(max(np.amax(positions1), np.amax(positions2))) + 1
    assert len(positions1) == len(
        positions2
    ), "Mismatch between lengths {} and {}".format(
        len(positions1), len(positions2)
    )
    sparse_matrix = sparse.coo_matrix(
        (contacts, (positions1, positions2)), shape=(n, n)
    )

    dense_matrix = np.array(sparse_matrix.todense(), dtype=np.int32)
    if output is not None:
        np.savetxt(output, dense_matrix, fmt="%i")

    return dense_matrix
