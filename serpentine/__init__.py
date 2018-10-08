#!/usr/bin/env python3
# coding: utf-8

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
