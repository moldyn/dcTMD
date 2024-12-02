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
MODE = 'nearest'
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
    """Compare two WorkEstimator or ForceEstimator instances via asserts."""
    assert type(estimator1) is type(estimator2)
    assert estimator1.temperature == estimator2.temperature
    assert estimator1.verbose == estimator2.verbose
    np.testing.assert_almost_equal(
        estimator1.W_mean_,
        estimator2.W_mean_,
        err_msg='Mean work test failed.',
    )
    np.testing.assert_almost_equal(
        estimator1.W_diss_,
        estimator2.W_diss_,
        err_msg='Dissipative work test failed.',
    )
    np.testing.assert_almost_equal(
        estimator1.dG_,
        estimator2.dG_,
        err_msg='Free energy test failed.',
    )
    np.testing.assert_almost_equal(
        estimator1.friction_,
        estimator2.friction_,
        decimal=6,
        err_msg='Friction test failed.',
    )
    np.testing.assert_almost_equal(
        estimator1.friction_smooth_,
        estimator2.friction_smooth_,
        decimal=6,
        err_msg='Smoothed friction test failed.',
    )


def test_WorkEstimator(ref_workestimator):
    workestimator_name = f'{TEST_FILE_DIR}/workestimator'
    estimator = storing.load(filename=workestimator_name)
    assert_estimator_equality(estimator, ref_workestimator)


def test_ForceEstimator(ref_forceestimator):
    forceestimator_name = f'{TEST_FILE_DIR}/forceestimator'
    estimator = storing.load(filename=forceestimator_name)
    assert_estimator_equality(estimator, ref_forceestimator)

    # assert that memory kernel does not raise error
    ref_forceestimator.memory_kernel([0])


@pytest.mark.parametrize(
    "set, Estimator",
    [('workset', dcTMD.WorkEstimator),
     ('forceset', dcTMD.ForceEstimator)],
)
def test_estimator_reset(set, Estimator):
    """Test if _reset is called by fit to delete friction_smooth_.

    When refitting a WorkEstimator with a new WorkSet, many attributes are
    calculated directly by the fit function. However, the smoothed friction
    is not directly calculated. As it no longer corresponds to the fitted
    WorkSet, it is deleted. This is what is tested for, here, and analogously
    for the ForceEstimator."
    """
    set_name = f'{TEST_FILE_DIR}/{set}'
    set = storing.load(filename=set_name)
    estimator = Estimator(temperature=TEMPERATURE)
    estimator.fit(set)
    estimator.smooth_friction(SIGMA, MODE)
    # Check if attribute exists
    assert hasattr(estimator, 'friction_smooth_')  # noqa: WPS421
    estimator.fit(set)
    # Check if attribute no longer exists
    assert not hasattr(estimator, 'friction_smooth_')  # noqa: WPS421
