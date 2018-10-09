#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Serpentine unit testing

Basic tests for the serpentine binning functions.

"""

import numpy as np
import pytest
import serpentine as serp

SIZE_PARAMETERS = ("matrix_size", [5, 10, 20, 50, 100])


@pytest.mark.parametrize(*SIZE_PARAMETERS)
def test_single_iteration_identical(matrix_size):
    """Test if two identical input matrices are binned the same exact way.
    """
    inputA = np.random.random((matrix_size, matrix_size))
    inputB = np.copy(inputA)
    outputA, outputB, diff = serp.serpentin_iteration(inputA, inputB)
    assert np.isclose(outputA, outputB).all()
    assert (diff == 0).all()


@pytest.mark.parametrize(*SIZE_PARAMETERS)
def test_total_smearing(matrix_size):
    """Test if the matrices are totally smeared into identical values if the
    threshold is greater than or equal to their summed values.
    """

    inputA = np.random.random((matrix_size, matrix_size))
    inputB = np.random.random((matrix_size, matrix_size))
    max_value = max(np.sum(inputA), np.sum(inputB))
    outputA, outputB, _ = serp.serpentin_iteration(
        inputA, inputB, threshold=max_value
    )
    assert np.isclose(
        outputA / np.average(outputA), np.ones((matrix_size, matrix_size))
    ).all()
    assert np.isclose(
        outputB / np.average(outputB), np.ones((matrix_size, matrix_size))
    ).all()


@pytest.mark.parametrize(*SIZE_PARAMETERS)
def test_no_binning(matrix_size):
    """Test if the matrices are not binned at all if the threshold value is
    lower than or equal to their minimum value.    
    """
    inputA = np.random.random((matrix_size, matrix_size))
    inputB = np.random.random((matrix_size, matrix_size))
    thresh_value = min(np.min(inputA), np.min(inputB))
    outputA, outputB, _ = serp.serpentin_iteration(
        inputA, inputB, threshold=thresh_value, minthreshold=thresh_value / 2.
    )
    assert np.isclose(outputA, inputA).all()
    assert np.isclose(outputB, inputB).all()


@pytest.mark.parametrize(*SIZE_PARAMETERS)
def test_symmetrical(matrix_size):
    """Test if output symmetry is respected when the inputs are symmetric.
    """
    inputA = np.random.random((matrix_size, matrix_size))
    inputB = np.random.random((matrix_size, matrix_size))
    outputA, outputB, _ = serp.serpentin_iteration(
        inputA, inputB, triangular=True
    )
    assert np.isclose(np.triu(outputA), np.tril(outputA).T).all()
    assert np.isclose(np.triu(outputB), np.tril(outputB).T).all()
