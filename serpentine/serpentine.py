#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Serpentine binning

An implementation of the so-called 'serpentine binning' procedure described
in Baudry et al.

Command line::

    Usage:
        serpentine.py [<matrixA>] [<matrixB>] [--threshold=auto] [--verbose]
                      [--min-threshold=auto] [--trend=high] [--triangular]
                      [--limit=3] [--demo] [--demo-size=500]

    Arguments:
        matrixA                         The first input matrix, in plain text
                                        CSV format. Optional in demo mode.
        matrixB                         The second input matrix, in plain text
                                        CSV format. Optional in demo mode or
                                        single binning mode.

    Options:
        -h, --help                      Display this help message.
        --version                       Display the program's current version.
        -t auto, --threshold auto       Threshold value to trigger binning.
                                        [default: auto]
        -m auto, --min-threshold auto   Minimum value to force trigger binning
                                        in either matrix. [default: auto]
        --trend high                    Trend to subtract to the differential
                                        matrix, possible values are "mean":
                                        equal amount of positive and negative
                                        differences, and "high": normalize
                                        at the regions with higher coverage.
                                        [default: high]
        --triangular                    Treat the matrix as triangular,
                                        useful when plotting matrices adjacent
                                        to the diagonal. [default: False]
        --limit 3                       Set the z-axis limit on the
                                        plot of the differential matrix.
                                        [default: 3]
        --demo                          Run a demo on randomly generated
                                        matrices. [default: False]
        --demo-size 500                 Size of the test matrix for the demo.
                                        [default: 500]
        -v, --verbose                   Show verbose output. [default: False]
"""

import sys
import numpy as _np
import pandas as _pd
import docopt as _doc
import itertools as _it
import collections as _col
from matplotlib import pyplot as _plt
from matplotlib import colors as _cols
import warnings as _warns
from random import choice as _choice
import multiprocessing as _mp
from datetime import datetime as _datetime
from serpentine.version import __version__
from typing import Tuple, Optional
import functools

_warns.filterwarnings(action="ignore")

DEFAULT_MIN_THRESHOLD = 10.0
DEFAULT_THRESHOLD = 50.0
DEFAULT_ITERATIONS = 10.0
DEFAULT_SIZE = 300.0
DEFAULT_PRECISION = 0.05


ASCII_SNAKE = """

                            Serpentine binning

        .';.;olc:'
    ;xOOOxkOkkkkkl.
        ...';lxOOOO;            ,cooc,                 ...
                .lk00o        .lkolllllxc.           'okxdxxd:.
                    ;dK0d:,,;cxOl:'......:dd;..   ..cxo,......:ld,
                    .c0XXKKKKd;;.        .;x0K00000o,.         .:xo,.
                    ..oxkxd:,;.           .,clooc'              .;oOOxoooxo
                        ...'.                                       .',;;,.
"""

def alternate_print(logfile):
    return functools.partial(print, file=open(logfile, "a"))

def serpentin_iteration_multi(
    M: _np.ndarray,
    threshold: float = DEFAULT_THRESHOLD,
    minthreshold: float = DEFAULT_MIN_THRESHOLD,
    triangular: bool = False,
    verbose: bool = True,
    offset: int = 0,
    get_bins: bool = False,
) -> _np.ndarray:

    """Perform a single iteration of serpentin binning, multiple matrices version

    Each serpentin binning is generally executed in multiple
    iterations in order to smooth the random variability in the
    bin aggregation. This funciton performs a single iteration.

    Parameters
    ----------
    M : array_like
        The matrices to be compared, as stacked 2D matrices
    threshold : float, optional
        The threshold of rebinning for the highest coverage matrix.
    minthreshold : float, optional
        The threshold for both matrices
    triangular : bool, optional
        Set triangular if you are interested in rebin only half of the
        matrix (for instance in the case of matrices which are
        already triangular, default is false)
    verbose : bool, optional
        Set it false if you are annoyed by the printed output.
    offset : int, optional
        Diagonals to ignore when performing the binning.
    get_bins : bool, optional
        Whether to return the identified bins, in which
        case it will be returned as tuple indexes. One element of the tuple
        per serpentine.
        Each serpentine can be used as an 2D index for each stack in M,
        to slice the values relative to its bins.
        Note: this options consumes a significant amount of memory.
        Default is False.

    Returns
    -------
    D : array_like
        A 4D matrix where the first two dimensions are indexes, containing:
        the rebinned matrices on the indices-diagonal,
        the log-ratio matrices out of diagonal, expressed in base 2.
        Attention, the log-ratio matrices needs to be individually normalized by subtracting
        an appropriate value for the zero (MDbefore or numpy.mean functions are there
        to help you in this task).
    bins : Tuple, optional
        A tuple containing the serpentines.
        Only returned if the supplied 'bins' parameter is True.
    """

    try:
        assert len(M.shape) == 3
    except:
            raise ValueError(
                "M should be stacked 2D arrays"
            )

    dim0, dim1, dim2 = M.shape

    if triangular:
        try:
            assert dim1 == dim2
        except AssertionError:
            raise ValueError(
                "Matrices must be square"
            )

    try:
        assert minthreshold < threshold
    except AssertionError:
        raise ValueError("Minimal threshold should be lower than maximal")

    def pixel_neighs_triangular(i, j, size):

        if i > 0:
            if i - 1 >= j:
                yield (i - 1, j)
        if i < size - 1:
            if i + 1 >= j:
                yield (i + 1, j)
        if j > 0:
            if i >= j - 1:
                yield (i, j - 1)
        if j < size - 1:
            if i >= j + 1:
                yield (i, j + 1)

    def pixel_neighs(i, j, w, h):

        if i > 0:
            yield (i - 1, j)
        if i < w - 1:
            yield (i + 1, j)
        if j > 0:
            yield (i, j - 1)
        if j < h - 1:
            yield (i, j + 1)

    U = _np.copy(M)
    U = U.reshape((dim0, dim1 * dim2))

    if triangular:
        pixels = [
            _np.array([i * dim2 + j], dtype=_np.int32)
            for (i, j) in _it.product(range(dim1), range(dim2))
            if i >= j and abs(i - j) >= offset
        ]

        neighs = [
            set(
                int((a * (a + 1) / 2)) + b
                for (a, b) in pixel_neighs_triangular(i, j, dim2)
            )
            for (i, j) in _it.product(range(dim1), range(dim2))
            if i >= j and abs(i - j) >= offset
        ]
        start = int(dim1 * (dim1 + 1) / 2)
        tot = start

    else:
        pixels = [
            _np.array([i * dim2 + j], dtype=_np.int32)
            for (i, j) in _it.product(range(dim1), range(dim2))
            if abs(i - j) >= offset
        ]
        neighs = [
            set(
                (a * dim2) + b
                for (a, b) in pixel_neighs(i, j, dim1, dim2)
            )
            for (i, j) in _it.product(range(dim1), range(dim2))
            if abs(i - j) >= offset
        ]
        start = dim1 * dim2
        tot = start

    previous_existent = 0
    current_existent = 1

    def print_iteration(i, tot, start, verbose):

        if not verbose:
            return
        percent = 100 * float(tot) / start
        iteration_string = "{}\t Total serpentines: {} ({} %)".format(
            i, tot, percent
        )
        print(iteration_string)

    # merger
    i = 0
    while current_existent != previous_existent:
        print_iteration(i, tot, start, verbose)
        i += 1
        tot = 0
        for serp in _np.random.permutation(range(len(pixels))):
            if pixels[serp] is not None:
                tot += 1
                # choose where to expand
                if pixels[serp].size == 1:  # Optimization for performances
                    summes = U[:,(pixels[serp])[0]]
                else:
                    summes = _np.sum(U[:,pixels[serp]],axis=1)

                thresh = _np.all(summes < threshold)
                # note, i do not understand why minthreshold if dim0 > 2
                minthresh = _np.any(summes < minthreshold)

                if thresh or minthresh:
                    try:
                        min_neigh = _choice(tuple(neighs[serp]))
                    except IndexError:
                        break

                    # Merge pixels
                    pixels[serp] = _np.concatenate(
                        (pixels[serp], pixels[min_neigh]), axis=0
                    )
                    # Merge neighbours (and remove self)
                    neighs[serp].remove(min_neigh)
                    neighs[min_neigh].remove(serp)
                    neighs[serp].update(neighs[min_neigh])

                    # Update neighbours of merged
                    for nneigh in neighs[min_neigh]:
                        neighs[nneigh].remove(min_neigh)
                        neighs[nneigh].add(serp)

                    # Delete merged serpentin
                    pixels[min_neigh] = None
                    neighs[min_neigh] = None

        previous_existent = current_existent
        current_existent = sum((serp is not None for serp in pixels))

    print_iteration(i, tot, start, verbose)
    if verbose:
        print("{}\t Over: {}".format(i, _datetime.now()))

    pix = (p for p in pixels if p is not None)

    U = U.astype(_np.float32)
    for serp in pix:
        U[:,serp] = _np.sum(U[:,serp],axis=1).reshape((dim0,1)) * 1.0 / len(serp)
    U = U.reshape((dim0, dim1, dim2))

    D = _np.zeros((dim0, dim0, dim1, dim2))

    if triangular:
        trili = _np.tril_indices(dim1)
        for i in range(dim0):
            D[i,i] = (
                _np.tril(U[i])
                + _np.transpose(_np.tril(U[i]))
                - _np.diag(_np.diag(_np.tril(U[i])))
            )
            for j in range(dim0):
                if i != j:
                    D[i,j,trili] = U[i,trili] * 1.0 / U[j,trili]
                    D[i,j,trili] = _np.log2(D[i,j,trili])

    else:
        D[_np.eye(dim0, dtype=bool)] = U
        for i in range(dim0):
            for j in range(dim0):
                D[i,j] = _np.log2(U[i] * 1.0 / U[j])

    if get_bins:
        # convert the bins in 2D coordinates
        bins = []
        pix = (p for p in pixels if p is not None)

        for p in pix:
            x = p // dim2
            y = p - dim2 * x
            bins.append((tuple(x), tuple(y)))

        return D, tuple(bins)
    else:
        return D

def serpentin_iteration(
    A: _np.ndarray,
    B: _np.ndarray,
    threshold: float = DEFAULT_THRESHOLD,
    minthreshold: float = DEFAULT_MIN_THRESHOLD,
    triangular: bool = False,
    verbose: bool = True,
    offset: int = 0,
    get_bins: bool = False,
) -> Tuple[_np.ndarray, _np.ndarray, _np.ndarray]:

    """Perform a single iteration of serpentin binning

    Each serpentin binning is generally executed in multiple
    iterations in order to smooth the random variability in the
    bin aggregation. This funciton performs a single iteration.

    Parameters
    ----------
    A, B : array_like
        The matrices to be compared.
    threshold : float, optional
        The threshold of rebinning for the highest coverage matrix.
    minthreshold : float, optional
        The threshold for both matrices
    triangular : bool, optional
        Set triangular if you are interested in rebin only half of the
        matrix (for instance in the case of matrices which are
        already triangular, default is false)
    verbose : bool, optional
        Set it false if you are annoyed by the printed output.
    offset : int, optional
        Diagonals to ignore when performing the binning.
    get_bins : bool, optional
        Whether to return the identified bins, in which
        case it will be returned as tuple indexes. One element of the tuple
        per serpentine.
        Each serpentine can be used as an 2D index for each stack in M,
        to slice the values relative to its bins.
        Note: this options consumes a significant amount of memory.
        Default is False.

    Returns
    -------
    Amod, Bmod : array_like
        The rebinned matrices
    D : array_like
        The log-ratio matrix, expressend in base 2. Attention, the
        matrix need to be normalized by subtractiong an appropriate
        value for the zero (MDbefore or numpy.mean functions are there
        to help you in this task).
    bins : Tuple, optional
        A tuple containing the serpentines.
        Only returned if the supplied 'bins' parameter is True.
    """

    try:
        assert(A.shape == B.shape)
    except AssertionError:
        raise ValueError(
            "Matrices must have identical shape"
        )

    M = _np.stack((A,B))

    if get_bins:
        sM, bins = serpentin_iteration_multi(M,
            threshold, minthreshold, triangular,
            verbose, offset, get_bins
        )
    else:
        sM = serpentin_iteration_multi(M,
            threshold, minthreshold, triangular,
            verbose, offset, get_bins
        )

    Amod = sM[0,0]
    Bmod = sM[1,1]
    D = sM[0,1]
    
    if get_bins:
        return Amod, Bmod, D, bins
    else:
        return Amod, Bmod, D

def _serpentin_iteration_multi_mp(value):
    return serpentin_iteration_multi(*value)


def serpentin_binning_multi(
    M: _np.ndarray,
    threshold: float = DEFAULT_THRESHOLD,
    minthreshold: float = DEFAULT_MIN_THRESHOLD,
    iterations: float = DEFAULT_ITERATIONS,
    precision: float = DEFAULT_PRECISION,
    triangular: bool = False,
    force_symmetric: bool = False,
    verbose: bool = True,
    parallel: int = 4 ,
    offset: int = 0,
    get_bins: bool = False, 
) -> Tuple:
    """Perform the serpentin binning, multi array API

    The function will perform the algorithm to serpentin bin two
    or more matrices, iterations can be done in series or in parallel,
    convinient for multi-processor machines.

    Parameters
    ----------
    M : array_like
        The matrices to be compared, as stacked 2D matrices.
    threshold : float, optional
        The threshold of rebinning for the highest coverage matrix. Default is
        set by the DEFAULT_THRESHOLD parameter, which is 50 if unchanged.
    minthreshold : float, optional
        The threshold for both matrices. Default is set by the
        DEFAULT_MIN_THRESHOLD parameter, which is 10 if unchanged.
    iterations : int, optional
        The number of iterations requested, more iterations will
        consume more time, but also will result in better and smoother
        results. If 0 or negative, iterations will continue until the average
        matrix after one more iteration doesn't differ by more than the
        precision parameter. Default is 10.
    precision : float, optional
        If the iterations parameter is 0 or negative, the iterations will
        continue until the average matrix after one more iteration doesn't
        differ by more than the precision parameter. Default is 0.05.
    triangular : bool, optional
        Set triangular if you are interested in rebinning only half of the
        matrix (for instance in the case of matrices which are
        already triangular, default is False).
    force_symmetric : bool, optional
    verbose : bool, optional
        Whether to print additional output during the computation. Default is
        False.
    parallel : int, optional
        Set it to the number of your processor if you want to attain
        maximum speeds. Default is 4.
    offset : int, optional
        Number of diagonals to ignore when performing the binning.
    get_bins : bool, optional
        Whether to report the identified bins, in which
        case it will be returned as tuple of a tuple of indexes. The
        outer tuple contains a tuple for each iteration, inside that, one element
        per serpentine.
        Each serpentine can be used as an 2D index for each stack in M,
        to slice the values relative to its bins.
        Note: this options consumes a significant amount of memory.
        Default is False.

    Returns
    -------
    sM : array_like
        A 4D matrix where the first two dimensions are indexes, containing:
        the rebinned matrices on the indices-diagonal,
        the log-ratio matrices out of diagonal, expressed in base 2.
        Attention, the log-ratio matrices needs to be individually normalized by subtracting
        an appropriate value for the zero (MDbefore or numpy.mean functions are there
        to help you in this task).
    bins : Tuple, optional
        A tuple containing the tuple of serpentines for each iteration. Only
        returned if the supplied 'bins' parameter is True.
    """

    try:
        assert len(M.shape) == 3
    except:
            raise ValueError(
                "M should be stacked 2D arrays"
            )

    dim0, dim1, dim2 = M.shape

    if triangular:
        try:
            assert dim1 == dim2
        except AssertionError:
            raise ValueError(
                "Matrices must be square"
            )

    try:
        assert minthreshold < threshold
    except AssertionError:
        raise ValueError("Minimal threshold should be lower than maximal")

    iterations = int(iterations)

    sM = _np.zeros((dim0, dim0, dim1, dim2))
    if get_bins:
        bins = []

    serp_size_distribution = _col.Counter()

    if parallel > 1:
        if verbose:
            print(
                "Starting {} binning processes in batches of {}...".format(
                    iterations, parallel
                )
            )
        p = _mp.Pool(parallel)
        iterator = (
            (M, threshold, minthreshold, triangular, verbose, offset, get_bins)
            for x in range(iterations)
        )
        res = p.map(_serpentin_iteration_multi_mp, iterator)

        for r in res:
            if get_bins:
                Mt, binst = r
                bins.append(binst)
            else:
                Mt = r
            sM = sM + Mt

    else:
        if verbose:
            print(
                "{} Starting {} binning processes...".format(
                    _datetime.now(), iterations
                )
            )
        if iterations > 0:
            for _ in range(int(iterations)):
                if get_bins:
                    Mt, binst = serpentin_iteration_multi(
                        M,
                        threshold=threshold,
                        minthreshold=minthreshold,
                        triangular=triangular,
                        verbose=verbose,
                        offset=offset,
                        get_bins=get_bins
                    )
                    bins.append(binst)
                else:
                    Mt = serpentin_iteration_multi(
                        M,
                        threshold=threshold,
                        minthreshold=minthreshold,
                        triangular=triangular,
                        verbose=verbose,
                        offset=offset,
                        get_bins=get_bins
                    )
                sM = sM + Mt
        else:
            iterations = 1
            current_diff = float("inf")
            while current_diff < precision:
                if get_bins:
                    Mt, binst = serpentin_iteration_multi(
                        M,
                        threshold=threshold,
                        minthreshold=minthreshold,
                        triangular=triangular,
                        verbose=verbose,
                        offset=offset,
                        get_bins=get_bins
                    )
                    bins.append(binst)
                else:
                    Mt = serpentin_iteration_multi(
                        M,
                        threshold=threshold,
                        minthreshold=minthreshold,
                        triangular=triangular,
                        verbose=verbose,
                        offset=offset,
                        get_bins=get_bins
                    )
                new_sM = sM + Mt
                sM_diff = _np.abs(
                    (new_sM / (iterations + 1)) - (sM / iterations)
                )
                if (
                    _np.amax(sM_diff) < precision
                ):
                    break

                else:
                    sM = new_sM
                    iterations += 1

    sM = sM * 1.0 / iterations

    if force_symmetric:
        for i in dim0:
            for j in dim1:
                sM = _np.tril(sM[i,j]) + _np.tril(sM[i,j]).T - _np.diag(_np.diag(sM[i,j]))

    if get_bins:
        return sM, tuple(bins)
    else:
        return sM


def serpentin_binning(
    A: _np.ndarray,
    B: _np.ndarray,
    threshold: float = DEFAULT_THRESHOLD,
    minthreshold: float = DEFAULT_MIN_THRESHOLD,
    iterations: float = DEFAULT_ITERATIONS,
    precision: float = DEFAULT_PRECISION,
    triangular: bool = False,
    force_symmetric: bool = False,
    verbose: bool = True,
    parallel: int = 4 ,
    offset: int = 0,
    get_bins: bool = False,
) -> Tuple:
    """Perform the serpentin binning

    The function will perform the algorithm to serpentin bin two
    matrices, iterations can be done in series or in parallel,
    convinient for multi-processor machines.

    Parameters
    ----------
    A, B : array_like
        The matrices to be compared.
    threshold : float, optional
        The threshold of rebinning for the highest coverage matrix. Default is
        set by the DEFAULT_THRESHOLD parameter, which is 50 if unchanged.
    minthreshold : float, optional
        The threshold for both matrices. Default is set by the
        DEFAULT_MIN_THRESHOLD parameter, which is 10 if unchanged.
    iterations : int, optional
        The number of iterations requested, more iterations will
        consume more time, but also will result in better and smoother
        results. If 0 or negative, iterations will continue until the average
        matrix after one more iteration doesn't differ by more than the
        precision parameter. Default is 10.
    precision : float, optional
        If the iterations parameter is 0 or negative, the iterations will
        continue until the average matrix after one more iteration doesn't
        differ by more than the precision parameter. Default is 0.05.
    triangular : bool, optional
        Set triangular if you are interested in rebinning only half of the
        matrix (for instance in the case of matrices which are
        already triangular, default is False).
    force_symmetric : bool, optional
    verbose : bool, optional
        Whether to print additional output during the computation. Default is
        False.
    parallel : int, optional
        Set it to the number of your processor if you want to attain
        maximum speeds. Default is 4.
    offset : int, optional
        Number of diagonals to ignore when performing the binning.
    get_bins : bool, optional
        Whether to report the identified bins, in which
        case it will be returned as tuple of a tuple of indexes. The
        outer tuple contains a tuple for each iteration, inside that, one element
        per serpentine.
        Each serpentine can be used as an index for the inputs A and B,
        to slice the values relative to its bins.
        Note: this options consumes a significant amount of memory.
        Default is False.

    Returns
    -------
    sA, sB : array_like
        The rebinned matrices
    sK : array_like
        The log-ratio matrix, expressend in base 2. Attention, the
        matrix needs to be normalized by subtracting an appropriate
        value for the zero (MDbefore or numpy.mean functions are there
        to help you in this task).
    bins : Tuple, optional
        A tuple containing the tuple of serpentines for each iteration. Only
        returned if the supplied 'bins' parameter is True.
    """

    try:
        assert(A.shape == B.shape)
    except AssertionError:
        raise ValueError(
            "Matrices must have identical shape"
        )

    M = _np.stack((A,B))

    if get_bins:
        sM, bins = serpentin_binning_multi(M,
            threshold=threshold,
            iterations=iterations,
            precision=precision,
            minthreshold=minthreshold,
            triangular=triangular,
            force_symmetric=force_symmetric,
            verbose=verbose,
            parallel=parallel,
            offset=offset,
            get_bins=get_bins
        )

        sA = sM[0,0]
        sB = sM[1,1]
        sK = sM[0,1]
        return sA, sB, sK, bins

    else:
        sM = serpentin_binning_multi(M,
            threshold=threshold,
            iterations=iterations,
            precision=precision,
            minthreshold=minthreshold,
            triangular=triangular,
            force_symmetric=force_symmetric,
            verbose=verbose,
            parallel=parallel,
            offset=offset,
            get_bins=get_bins
        )

        sA = sM[0,0]
        sB = sM[1,1]
        sK = sM[0,1]
        return sA, sB, sK


def _MDplot(ACmean, ACdiff, trans, xlim=None, ylim=None):
    _plt.scatter(ACmean, ACdiff - trans)
    _plt.xlabel("Log2 Mean contact number")
    _plt.ylabel("Log2 ratio")
    if xlim is not None:
        _plt.xlim(xlim[0], xlim[1])
    if ylim is not None:
        _plt.ylim(ylim[0], ylim[1])


def mad(data: _np.ndarray, axis: Optional[int] = None) -> _np.float64:

    """Median absolute deviation

    Calculates the median absolute deviation of data

    Parameters
    ----------
    data : array_like
        The dataset
    axis : int, optional
        The axis over which perform the numpy.median function

    Returns
    -------
    float
        The median absolute deviation
    """

    return _np.median(_np.absolute(data - _np.median(data, axis)), axis)


def outstanding_filter(X: _np.ndarray, sds: Optional[float] = 2.0) -> _np.ndarray:

    """Generate filtering index that removes outstanding values (three standard
    MADs above or below the median).

    Parameters
    ----------
    X : array_like
        The dataset

    Returns
    -------
    array_like
        The boolean filter

    Example
    -------

        >>> import numpy as np
        >>> X = np.arange(25).reshape((5,5))
        >>> X += X.T
        >>> X[2,2] += 10000
        >>> print(X)
        [[    0     6    12    18    24]
         [    6    12    18    24    30]
         [   12    18 10024    30    36]
         [   18    24    30    36    42]
         [   24    30    36    42    48]]
        >>> O = outstanding_filter(X)
        >>> print(O)
        [False False  True False False]
    """

    with _np.errstate(divide="ignore", invalid="ignore"):
        norm = _np.log10(_np.sum(X, axis=0))
        median = _np.median(norm)
        sigma = 1.4826 * mad(norm)

    return (norm < median - sds * sigma) + (norm > median + sds * sigma)


def fltmatr(X: _np.ndarray, flt: _np.ndarray) -> _np.ndarray:

    """Filter a 2D matrix in both dimensions according to an boolean
    vector.


    Parameters
    ----------
    X : array_like
        The input matrix
    flr : array_like
        The boolean filter

    Returns
    -------
    array_like
        The filtered matrix

    Example
    -------

        >>> import numpy as np
        >>> M = np.ones((5, 5))
        >>> M[2:4, 2:4] = 2
        >>> flt = M.sum(axis=1) > 5
        >>> X = fltmatr(M, flt)
        >>> print([X[0,0], X[0,1], X[1,0], X[1,1]])
        [2.0, 2.0, 2.0, 2.0]

    """

    X = _np.copy(X)
    X = X[flt, :]
    X = X[:, flt]
    return X


def _madplot(
    ACmean, ACdiff, s=10, xlim=None, ylim=None, showthr=True, show=True
):
    df = _pd.DataFrame({"m": ACmean, "d": ACdiff})

    with _warns.catch_warnings():
        _warns.simplefilter("ignore")

        df = df[_np.logical_not(_np.isinf(abs(df["m"])))]
        df = df[_np.logical_not(_np.isinf(abs(df["d"])))]

        df = df.sort_values(by="m", ascending=False)
        x = _np.zeros(s)
        y1 = _np.zeros(s)
        y2 = _np.zeros(s)
        q = (max(df["m"]) - min(df["m"])) / (s)

        k = 0
        for i in _np.arange(min(df["m"]), max(df["m"]), q):
            r = df[(df["m"] > i) * (df["m"] < i + q)]
            x[k] = _np.median(r["m"])
            y1[k] = _np.median(r["d"])
            y2[k] = 1.4826 * mad(r["d"])
            if _np.isnan(y2[k]):
                y2[k] = 0
            k = k + 1

    y1lim = _np.mean(y1[-3:])
    y2lim = _np.mean(y2[-3:])
    y2limv = _np.std(y2[-3:])
    if show:
        _MDplot(ACmean, ACdiff, y1lim, xlim=xlim, ylim=ylim)
        _plt.plot(x, y1 - y1lim, color="y")
        _plt.plot(x, y2, color="r")

    if _np.sum(_np.abs(y2 - y2lim) > y2limv * 2) > 0:
        xa = (x[(_np.abs(y2 - y2lim) > y2limv * 2)])[-1]
    else:
        xa = _np.percentile(ACmean[ACmean > 0], 99)

    if _np.isnan(xa) or _np.isinf(xa) or xa < _np.log2(25.0):
        xa = _np.log2(DEFAULT_THRESHOLD)

    if showthr and show:
        _plt.axvline(x=xa)
    return y1lim, 2 ** xa


def MDbefore(XA, XB, s=10, xlim=None, ylim=None, triangular=False, show=True):

    """MD plot before binning

    The MD plot is the main metric provided by the serpentin binning
    package, it visualized the variability in function of the mean
    coverage of a couple of matrices to be compared. This version is
    optimized to be run before the serpentin binning. The return
    values of this funciton will be hints on the values to use to
    normalize the differential matrix and as a threshold for the
    serpentin binning algorithm.

    Parameters
    ----------
    XA, XB : array_like
        The input matrices
    s : int, optional
        How many point use for the trend lines, depends on statistics
        and range
    xlim : float, optional
        Limits for the x-axis
    ylim : float, optional
        Limits for the y-axis
    triangular : bool, optional
        Set triangular if you are interested in rebin only half of the
        matrix (for instance in the case of matrices which are
        already triangular, default is false)
    show : bool, optional
        Set it to false if you are not interested in the graph but
        only in the return values of this function.

    Returns
    -------
    trend : float
        This is the extrapolated trend of the two matrices, it is
        measured by the ratio of the most covered part of the two
        matrices, use it to set the zero to the output of the
        differential analysis if you think this works better than the
        np.mean() function.
    threshold : float
        If the statistics permits the automatical extimation of the
        threshold, this will be a good parameter to use as input to
        the serpentin binning algorithm
    """

    if triangular:
        with _np.errstate(divide="ignore", invalid="ignore"):
            ACmean = _np.log2(
                (_np.tril(XA).flatten() + _np.tril(XB).flatten()) * 1.0 / 2
            )
            ACdiff = _np.log2(_np.tril(XB).flatten() / _np.tril(XA).flatten())
    else:
        with _np.errstate(divide="ignore", invalid="ignore"):
            ACmean = _np.log2(XA.flatten() + XB.flatten()) * 1.0 / 2
            ACdiff = _np.log2(XB.flatten() / XA.flatten())

    return _madplot(ACmean, ACdiff, s, xlim, ylim, show=show)


def MDafter(
    XA, XB, XD, s=10, xlim=None, ylim=None, triangular=False, show=True
):

    """MD plot after binning

    The MD plot is the main metric provided by the serpentin binning
    package, it visualized the variability in function of the mean
    coverage of a couple of matrices to be compared. This version is
    optimized to be run after the serpentin binning. The return values
    of this function should be generally ignored.

    Parameters
    ----------
    XA, XB : array_like
        The input matrices
    XD : array_like
        The differential matrix obtained by serpentin binnning
    s : int, optional
        How many point use for the trend lines, depends on statistics
        and range
    xlim : float, optional
        Limits for the x-axis
    ylim : float, optional
        Limits for the y-axis
    triangular : bool, optional
        Set triangular if you are interested in rebin only half of the
        matrix (for instance in the case of matrices which are
        already triangular, default is false)
    show : bool, optional
        Set it to false if you are not interested in the graph but
        only in the return values of this function.

    Returns
    -------
    trend, threshold : float
        Normally something which should not bother you
    """

    if triangular:
        with _np.errstate(divide="ignore", invalid="ignore"):
            ACmean = _np.log2(
                (_np.tril(XA).flatten() + _np.tril(XB).flatten()) * 1.0 / 2
            )
            ACdiff = _np.tril(XD).flatten()
    else:
        with _np.errstate(divide="ignore", invalid="ignore"):
            ACmean = _np.log2(XA.flatten() + XB.flatten()) * 1.0 / 2
            ACdiff = XD.flatten()

    return _madplot(ACmean, ACdiff, s, xlim, ylim, showthr=False, show=show)


def dshow(
    dif, trend, limit=3, triangular=False, colorbar=True, cmap=None, ax=_plt
):

    """Show differential matrix

    A boilerplate around the imshow matplotlib function to show the
    differential matrix

    Parameters
    ----------
    dif : array_like
        The differential matrix
    trend : float
        The value of the zero, please use either the output of
        MDbefore function or the value of md.mean(dif)
    limit : float, optional
        The colorscale limit of the log-ratio, setting it to 2 or 3
        seems like a sensible choice. Defaults to 3
    triangular : bool, optional
        Set triangular if you are interested in rebin only half of the
        matrix (for instance in the case of matrices which are
        already triangular, default is false)
    colorbar: bool, optional
        Whether to include a colorbar in the plot display. Default is True.
    cmap: str, optional
        Color map of the plotted matrix. Should be ideally diverging, default
        is sesismic.
    ax: optional
        Set axis, defaults to matplotlib library

    Returns
    -------
    The plot
    """

    if cmap is None:
        colors = [
            (120.0 / 350 / 2, 180.0 / 350 / 2, 230.0 / 350 / 2),
            (179.0 / 255, 205.0 / 255, 227.0 / 255),
            (1, 1, 1),
            (251.0 / 255, 180.0 / 255, 174.0 / 255),
            (248.0 / 350 / 2, 120.0 / 350 / 2, 109.0 / 350 / 2),
        ]
        cmap_name = "pastelpentine"
        cm = _cols.LinearSegmentedColormap.from_list(cmap_name, colors, N=21)
    else:
        cm = cmap

    if triangular:
        trili = _np.tril_indices(_np.int(_np.sqrt(dif.size)))
        triui = _np.triu_indices(_np.int(_np.sqrt(dif.size)))
        diagi = _np.diag_indices(_np.int(_np.sqrt(dif.size)))
        plotta = _np.copy(dif)
        plotta[trili] = plotta[trili] - trend
        plottadiag = plotta[diagi]
        plotta[triui] = _np.nan
        plotta[diagi] = plottadiag
    else:
        plotta = _np.copy(dif) - trend

    im = ax.imshow(
        plotta, vmin=-limit, vmax=limit, cmap=cm, interpolation="none"
    )
    if colorbar:
        _plt.colorbar(im)


def mshow(XX, subplot=_plt, colorbar=True, triangular=False):

    """Boilerplate around the imshow function to show a matrix.
    """

    del triangular
    colors = [(1, 1, 1), (1, 0, 0), (0, 0, 0)]
    cmap_name = "radios"
    cm2 = _cols.LinearSegmentedColormap.from_list(cmap_name, colors, N=64)

    with _np.errstate(divide="ignore", invalid="ignore"):
        im = subplot.imshow(_np.log10(XX), cmap=cm2, interpolation="none")
        if colorbar:
            _plt.colorbar(im)

    return im


def fromupdiag(filename):

    """Load a DADE matrix into a numpy array
    """

    result = None
    with open(filename) as f:
        header = f.readline().split()
        header.pop(0)
        total = len(header)
        result = _np.zeros((total, total))
        count = 0
        for line in f:
            data = line.split()
            data.pop(0)
            len(data)
            result[count, count:total] = data
            count += 1

    result = result + _np.transpose(result) - _np.diag(_np.diag(result))
    return result


def _plot(U, V, W, cmap=None, colorbar=True, triangular=False, limit=3):

    fig = _plt.figure()

    ax1 = fig.add_subplot(2, 2, 1)
    im1 = ax1.imshow(
        U,
        interpolation="none",
        clim=(0, 27),
        vmax=_np.percentile(U, 99),
        cmap="Reds",
    )
    if colorbar:
        _plt.colorbar(im1)

    ax2 = fig.add_subplot(2, 2, 2)
    im2 = ax2.imshow(
        V,
        interpolation="none",
        clim=(0, 27),
        vmax=_np.percentile(U, 99),
        cmap="Reds",
    )
    if colorbar:
        _plt.colorbar(im2)

    ax3 = fig.add_subplot(2, 2, 3)
    dshow(
        W,
        0,
        limit=limit,
        triangular=triangular,
        colorbar=colorbar,
        cmap=cmap,
        ax=ax3,
    )


def _demo(
    threshold=DEFAULT_THRESHOLD,
    minthreshold=DEFAULT_MIN_THRESHOLD,
    size=DEFAULT_SIZE,
    triangular=True,
    limit=3,
    trend="mean",
    verbose=True,
):

    """Perform binning on a random matrix
    """

    _np.random.seed(15)

    A = _np.random.random((size, size)) * 10

    _np.random.seed(80)

    B = _np.random.random((size, size)) * 10
    A = A + A.T
    B = B + B.T

    if threshold == "auto" or trend == "high":
        mdthreshold, mdtrend = MDbefore(
            A, B, s=10, triangular=triangular, show=False
        )
    if threshold == "auto":
        threshold = mdthreshold
    if minthreshold == "auto":
        minthreshold = threshold / 5.0
    if trend == "high":
        trend = mdtrend

    U, V, W = serpentin_binning(
        A,
        B,
        threshold=threshold,
        minthreshold=minthreshold,
        parallel=4,
        triangular=triangular,
        verbose=verbose,
    )
    _plot(U, V, W, triangular=triangular, limit=limit)
    _plt.show()


def extract_serpentines(M):
    """Extract serpentine structure

    Isolate serpentines based on shared pixel values
    in a contact map.

    Parameters
    ----------
    M : numpy.ndarray
        The (serpentinized) input contact map
    """

    values = set(_it.chain(*M))
    for value in values:
        serpentine_mask = _np.nonzero(M == value)
        yield list(zip(*serpentine_mask))


def barycenter(serp, weights=None):
    """Compute weighted serpentine barycenter

    Compute the (optionally weighted) barycenter of a
    serpentine, where the weights would be equal to the values
    of the map itself.

    Parameters
    ----------
    serp : iterable
        An iterable of serpentine coordinates
    weights : numpy.ndarray or None, optional
        If None, the barycenter is unweighted. Otherwise, if it
        is a contact map, the barycenter is weighted by the values
        of the map at the serpentine's coordinates. Default is None.

    Returns
    -------
    bary : tuple
        The barycenter coordinates
    """

    bary = _np.zeros(2)
    serp_total = 0
    for coord in serp:
        if weights is not None:
            bary += _np.array(coord) * weights[coord]
            serp_total += weights[coord]
        else:
            bary += _np.array(coord)
            serp_total += 1

    return tuple(bary / serp_total)


def all_barycenters(M_serp, weights=None):
    """Compute all serpentine barycenters

    Extract all serpentines from a serpentinized matrix,
    then compute all serpentine barycenters, optionally
    weigthed by the non-serpentinized matrix.

    Parameters
    ----------
    M_serp : numpy.ndarray
        The serpentinized matrix
    weights : numpy.ndarray or None, optional
        The non-serpentinized (original) matrix acting
        as weights. Default is None.
    """

    serps = extract_serpentines(M_serp)

    for serp in serps:
        yield barycenter(serp, weights=weights)


def _main():

    arguments = _doc.docopt(__doc__, version=__version__)

    inputA = arguments["<matrixA>"]
    inputB = arguments["<matrixB>"]
    threshold = arguments["--threshold"]
    minthreshold = arguments["--min-threshold"]
    size = int(arguments["--demo-size"])
    triangular = arguments["--triangular"]
    limit = int(arguments["--limit"])
    trend = arguments["--trend"]
    is_demo = int(arguments["--demo"])
    verbose = arguments["--verbose"]

    if threshold != "auto":
        threshold = float(threshold)
    if minthreshold != "auto":
        minthreshold = float(minthreshold)
    if trend != "mean" and trend != "high":
        print('Error! The --trend option accepts only "mean" or "high" values')
        return

    if is_demo:
        if verbose:
            print(ASCII_SNAKE)
        _demo(
            threshold=threshold,
            minthreshold=minthreshold,
            size=size,
            triangular=triangular,
            limit=limit,
            verbose=verbose,
            trend=trend,
        )

    elif inputA and inputB:
        if verbose:
            print(ASCII_SNAKE)
        try:
            A = fromupdiag(inputA)
            B = fromupdiag(inputB)
            triangular = True
        except Exception:
            try:
                A = _np.loadtxt(inputA)
                B = _np.loadtxt(inputB)
            except Exception as e:
                print(e)
                print(
                    "Error when processing {} or {}, please check that the "
                    "files exist with reading permissions and that they have "
                    "the correct format.".format(inputA, inputB)
                )
                return

        if threshold == "auto" or trend == "high":
            mdtrend, mdthreshold = MDbefore(
                A, B, s=10, triangular=triangular, show=False
            )
        if threshold == "auto":
            threshold = mdthreshold
        if minthreshold == "auto":
            minthreshold = threshold / 5.0
        if trend == "high":
            trend = mdtrend

        U, V, W = serpentin_binning(
            A,
            B,
            threshold=threshold,
            minthreshold=minthreshold,
            triangular=triangular,
            verbose=verbose,
        )
        if trend == "mean":
            if triangular:
                trili = _np.tril_indices(_np.int(_np.sqrt(W.size)))
                trend = _np.mean(W[trili])
            else:
                trend = _np.mean(W)

        _plot(
            A,
            B,
            _np.log2(B / A) - _np.log2(_np.mean(B) / _np.mean(A)),
            triangular=triangular,
            limit=limit,
        )
        _plot(U, V, W - trend, triangular=triangular, limit=limit)
        _plt.show()

    else:
        print(
            "Error: there must be either two matrix files as inputs (in csv "
            "format) as arguments or the --demo flag must be supplied.\n"
            "See 'serpentine --help' for a list of all options and arguments."
        )


if __name__ == "__main__":

    _main()
