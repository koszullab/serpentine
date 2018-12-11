#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Serpentine binning for Hi-C contact maps

Provides:
    -An implementation of the serpentine binning procedure described in
    Scolari et al., usable on single or multiple contact maps
    -A set of functions for quickly visualizing the effects of serpentine
    binning in terms of contact distribution (median absolute deviation, etc.)

"""

from .serpentine import (
    dshow,
    fltmatr,
    fromupdiag,
    mad,
    MDafter,
    MDbefore,
    mshow,
    outstanding_filter,
    serpentin_binning,
    serpentin_iteration,
)
from .benchmarks import *

from .version import __version__ as version

__author__ = "Cluster Buster (scovit, a.k.a. Vittore F. Scolari), \
              Lyamovich (baudrly, a.k.a. Lyam Baudry)"
__copyright__ = "Copyright © 2017-2018, Institut Pasteur, Paris, France"
__credits__ = ["Cluster Buster", "Lyamovich"]
__license__ = "Artistic"
__maintainer__ = "Cluster Buster"
__email__ = "vittore.scolari@pasteur.fr"
__status__ = "Pre-Alpha"
__version__ = version
