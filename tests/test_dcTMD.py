# -*- coding: utf-8 -*-
"""Tests for the clustering module."""

import numpy as np
from os.path import dirname, join
import pytest
from dcTMD import storing, dcTMD


VELOCITY = 0.001
RESOLUTION = 1
VERBOSE = True
TEMPERATURE = 300
INDICES = np.array([1, 3, 5])
SIGMA = 0.1
MODE = 'reflect'
HERE = dirname(__file__)
TEST_FILE_DIR = join(HERE, 'testdata')


@pytest.fixture
def ref_workestimator(scope="session"):
    workset_name = f'{TEST_FILE_DIR}/workset'
    workset = storing.load(filename=workset_name)
    estimator = dcTMD.WorkEstimator(temperature=TEMPERATURE)
    estimator.fit(workset)
    estimator.smooth_friction(SIGMA, MODE)
    return estimator


@pytest.fixture
def ref_forceestimator(scope="session"):
    forceset_name = f'{TEST_FILE_DIR}/forceset'
    forceset = storing.load(filename=forceset_name)
    estimator = dcTMD.ForceEstimator(temperature=TEMPERATURE)
    estimator.fit(forceset)
    estimator.smooth_friction(SIGMA, MODE)
    return estimator


def assert_estimator_equality(estimator1, estimator2):
    """Compare to WorkEstimator or ForceEstimator instances via asserts."""
    assert type(estimator1) is type(estimator2)
    assert estimator1.temperature == estimator2.temperature
    assert estimator1.verbose == estimator2.verbose
    np.testing.assert_almost_equal(
        estimator1.W_mean_,
        estimator2.W_mean_,
    )
    np.testing.assert_almost_equal(
        estimator1.W_diss_,
        estimator2.W_diss_,
    )
    np.testing.assert_almost_equal(
        estimator1.dG_,
        estimator2.dG_,
    )
    np.testing.assert_almost_equal(
        estimator1.friction_,
        estimator2.friction_,
        decimal=6,
    )
    np.testing.assert_almost_equal(
        estimator1.friction_smooth_,
        estimator2.friction_smooth_,
        decimal=6,
    )


def test_WorkEstimator(ref_workestimator):
    workestimator_name = f'{TEST_FILE_DIR}/workestimator'
    estimator = storing.load(filename=workestimator_name)
    assert_estimator_equality(estimator, ref_workestimator)


def test_ForceEstimator(ref_forceestimator):
    forceestimator_name = f'{TEST_FILE_DIR}/forceestimator'
    estimator = storing.load(filename=forceestimator_name)
    assert_estimator_equality(estimator, ref_forceestimator)
