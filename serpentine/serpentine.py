#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Serpentine binning

An implementation of the so-called 'serpentine binning' procedure described
in Scolari et al.

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
from matplotlib import pyplot as _plt
from matplotlib import colors as _cols
import warnings as _warns
from random import choice as _choice
import multiprocessing as _mp
from datetime import datetime as _datetime
from serpentine.version import __version__
from typing import Tuple, Optional

_warns.filterwarnings(action="ignore")

DEFAULT_MIN_THRESHOLD = 10.
DEFAULT_THRESHOLD = 50.
DEFAULT_ITERATIONS = 10.
DEFAULT_SIZE = 300.

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

sys.stdout = sys.stderr


def serpentin_iteration(
    A: _np.ndarray,
    B: _np.ndarray,
    threshold: float = DEFAULT_THRESHOLD,
    minthreshold: float = DEFAULT_MIN_THRESHOLD,
    triangular: bool = False,
    verbose: bool = True,
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

    Returns
    -------
    Amod, Bmod : array_like
        The rebinned matrices
    D : array_like
        The log-ratio matrix, expressend in base 2. Attention, the
        matrix need to be normalized by subtractiong an appropriate
        value for the zero (MDbefore or numpy.mean functions are there
        to help you in this task).
    """

    if triangular:
        try:
            assert A.shape == B.shape
            assert len(A.shape) == 2
            assert min(A.shape) == max(A.shape)
        except AssertionError:
            raise ValueError(
                "Matrices must be square and have identical shape"
            )
    else:
        try:
            assert A.shape == B.shape
            assert len(A.shape) == 2
        except AssertionError:
            raise ValueError("Matrices must have identical shape")

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

    size = A.shape
    U = _np.copy(A)
    V = _np.copy(B)
    U = U.reshape((size[0] * size[1]))
    V = V.reshape((size[0] * size[1]))

    if triangular:
        pixels = [
            _np.array([i * size[0] + j], dtype=_np.int32)
            for (i, j) in _it.product(range(size[0]), range(size[0]))
            if i >= j
        ]

        neighs = [
            set(
                int((a * (a + 1) / 2)) + b
                for (a, b) in pixel_neighs_triangular(i, j, size[0])
            )
            for (i, j) in _it.product(range(size[0]), range(size[0]))
            if i >= j
        ]
        start = int(size[0] * (size[0] + 1) / 2)
        tot = start

    else:
        pixels = [
            _np.array([i * size[1] + j], dtype=_np.int32)
            for (i, j) in _it.product(range(size[0]), range(size[1]))
        ]
        neighs = [
            set(
                (a * size[1]) + b
                for (a, b) in pixel_neighs(i, j, size[0], size[1])
            )
            for (i, j) in _it.product(range(size[0]), range(size[1]))
        ]
        start = size[0] * size[1]
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
                    a = U[(pixels[serp])[0]]
                    b = V[(pixels[serp])[0]]
                else:
                    a = _np.sum(U[pixels[serp]])
                    b = _np.sum(V[pixels[serp]])

                thresh = a < threshold and b < threshold
                minthresh = a < minthreshold or b < minthreshold
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
    V = V.astype(_np.float32)
    for serp in pix:
        U[serp] = _np.sum(U[serp]) * 1. / len(serp)
        V[serp] = _np.sum(V[serp]) * 1. / len(serp)
    U = U.reshape((size[0], size[1]))
    V = V.reshape((size[0], size[1]))

    if triangular:
        Amod = (
            _np.tril(U)
            + _np.transpose(_np.tril(U))
            - _np.diag(_np.diag(_np.tril(U)))
        )
        Bmod = (
            _np.tril(V)
            + _np.transpose(_np.tril(V))
            - _np.diag(_np.diag(_np.tril(V)))
        )
        trili = _np.tril_indices(_np.int(_np.sqrt(Bmod.size)))
        D = _np.zeros_like(Bmod)
        D[trili] = V[trili] * 1. / U[trili]
        D[trili] = _np.log2(D[trili])

    else:
        Amod = U
        Bmod = V
        D = V * 1. / U
        D = _np.log2(D)

    return (Amod, Bmod, D)


def _serpentin_iteration_mp(value):
    return serpentin_iteration(*value)


def serpentin_binning(
    A: _np.ndarray,
    B: _np.ndarray,
    threshold: float = DEFAULT_THRESHOLD,
    minthreshold: float = DEFAULT_MIN_THRESHOLD,
    iterations: float = DEFAULT_ITERATIONS,
    triangular: bool = False,
    verbose: bool = True,
    parallel: int = 4,
) -> Tuple[_np.ndarray, _np.ndarray, _np.ndarray]:

    """Perform the serpentin binning

    The function will perform the algorithm to serpentin bin two
    matrices, iterations can be done in series or in parallel,
    convinient for multi-processor machines.

    Parameters
    ----------
    A, B : array_like
        The matrices to be compared.
    threshold : float, optional
        The threshold of rebinning for the highest coverage matrix.
    minthreshold : float, optional
        The threshold for both matrices
    iterations : int, optional
        The number of iterations requested, more iterations will
        consume more time, but also will result in better and smoother
        results
    triangular : bool, optional
        Set triangular if you are interested in rebinning only half of the
        matrix (for instance in the case of matrices which are
        already triangular, default is False)
    verbose : bool, optional
        Set it false if you are annoyed by the printed output.
    parallel : int, optional
        Set it to the number of your processor if you want to attain
        maximum speeds

    Returns
    -------
    sA, sB : array_like
        The rebinned matrices
    sK : array_like
        The log-ratio matrix, expressend in base 2. Attention, the
        matrix needs to be normalized by subtracting an appropriate
        value for the zero (MDbefore or numpy.mean functions are there
        to help you in this task).
    """

    if triangular:
        try:
            assert A.shape == B.shape
            assert len(A.shape) == 2
            assert min(A.shape) == max(A.shape)
        except AssertionError:
            raise ValueError(
                "Matrices must be square and have identical shape"
            )
    else:
        try:
            assert A.shape == B.shape
            assert len(A.shape) == 2
        except AssertionError:
            raise ValueError("Matrices must have identical shape")
    try:
        assert minthreshold < threshold
    except AssertionError:
        raise ValueError("Minimal threshold should be lower than maximal")

    iterations = int(iterations)

    sK = _np.zeros_like(A)
    sA = _np.zeros_like(A)
    sB = _np.zeros_like(A)

    if parallel > 1:
        if verbose:
            print(
                "Starting {} binning processes in batches of {}...".format(
                    iterations, parallel
                )
            )
        p = _mp.Pool(parallel)
        iterator = (
            (A, B, threshold, minthreshold, triangular, verbose)
            for x in range(iterations)
        )
        res = p.map(_serpentin_iteration_mp, iterator)

        for r in res:
            At, Bt, Kt = r
            sK = sK + Kt
            sA = sA + At
            sB = sB + Bt

    else:
        if verbose:
            print(
                "{} Starting {} binning processes...".format(
                    _datetime.now(), iterations
                )
            )
        for _ in range(iterations):
            At, Bt, Kt = serpentin_iteration(
                A,
                B,
                threshold=threshold,
                minthreshold=minthreshold,
                triangular=triangular,
                verbose=verbose,
            )
            sK = sK + Kt
            sA = sA + At
            sB = sB + Bt

    sK = sK * 1. / iterations
    sA = sA * 1. / iterations
    sB = sB * 1. / iterations

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


def outstanding_filter(X: _np.ndarray) -> _np.ndarray:

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

    return (norm < median - 3 * sigma) + (norm > median + 3 * sigma)


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
        >>> print(M)
        [[ 1.  1.  1.  1.  1.]
         [ 1.  1.  1.  1.  1.]
         [ 1.  1.  2.  2.  1.]
         [ 1.  1.  2.  2.  1.]
         [ 1.  1.  1.  1.  1.]]
        >>> flt = M.sum(axis=1) > 5
        >>> X = fltmatr(M, flt)
        >>> print(X)
        [[ 2.  2.]
         [ 2.  2.]]

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

    if _np.isnan(xa) or _np.isinf(xa) or xa < _np.log2(25.):
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
                (_np.tril(XA).flatten() + _np.tril(XB).flatten()) * 1. / 2
            )
            ACdiff = _np.log2(_np.tril(XB).flatten() / _np.tril(XA).flatten())
    else:
        with _np.errstate(divide="ignore", invalid="ignore"):
            ACmean = _np.log2(XA.flatten() + XB.flatten()) * 1. / 2
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
                (_np.tril(XA).flatten() + _np.tril(XB).flatten()) * 1. / 2
            )
            ACdiff = _np.tril(XD).flatten()
    else:
        with _np.errstate(divide="ignore", invalid="ignore"):
            ACmean = _np.log2(XA.flatten() + XB.flatten()) * 1. / 2
            ACdiff = XD.flatten()

    return _madplot(ACmean, ACdiff, s, xlim, ylim, showthr=False, show=show)


def dshow(dif, trend, limit=3, triangular=False, cmap=None, ax=_plt):

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
            (120. / 350 / 2, 180. / 350 / 2, 230. / 350 / 2),
            (179. / 255, 205. / 255, 227. / 255),
            (1, 1, 1),
            (251. / 255, 180. / 255, 174. / 255),
            (248. / 350 / 2, 120. / 350 / 2, 109. / 350 / 2),
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


def _plot(U, V, W, cmap=None, triangular=False, limit=3):

    fig = _plt.figure()

    ax1 = fig.add_subplot(2, 2, 1)
    im1 = ax1.imshow(
        U,
        interpolation="none",
        clim=(0, 27),
        vmax=_np.percentile(U, 99),
        cmap="Reds",
    )
    _plt.colorbar(im1)

    ax2 = fig.add_subplot(2, 2, 2)
    im2 = ax2.imshow(
        V,
        interpolation="none",
        clim=(0, 27),
        vmax=_np.percentile(U, 99),
        cmap="Reds",
    )
    _plt.colorbar(im2)

    ax3 = fig.add_subplot(2, 2, 3)
    dshow(W, 0, limit=limit, triangular=triangular, cmap=cmap, ax=ax3)


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
        minthreshold = threshold / 5.
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
            minthreshold = threshold / 5.
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
