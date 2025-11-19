# -*- coding: utf-8 -*-
"""Tests for the clustering module."""

import numpy as np
from os.path import dirname, join
import pytest
from dcTMD import storing, dcTMD


VELOCITY = 0.001
RESOLUTION = 1
RESOLUTION10 = 10
VERBOSE = True
TEMPERATURE = 300
SIGMA = 0.1
MODE = 'nearest'
HERE = dirname(__file__)
TEST_FILE_DIR = join(HERE, 'testdata')


# test memory kernel method of ForceEstimator
def _valid_indices(forceestimator):
    n = len(forceestimator.force_set.time_)
    # Choose three well-spaced, guaranteed-valid indices.
    # Ensure they are >=1 and <= n-2 to avoid boundary weirdness.
    return np.array([1, max(2, n // 3), n - 2], dtype=int)


@pytest.fixture
def ref_workestimator(scope="session"):
    workset_name = f'{TEST_FILE_DIR}/workset'
    workset = storing.load(filename=workset_name)
    estimator = dcTMD.WorkEstimator(temperature=TEMPERATURE)
    estimator.fit(workset)
    estimator.smooth_friction(SIGMA, MODE)
    return estimator


@pytest.fixture
def ref_workestimator_res10(scope="session"):
    workset_name = f'{TEST_FILE_DIR}/workset_res10'
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


@pytest.fixture
def ref_forceestimator_res10(scope="session"):
    forceset_name = f'{TEST_FILE_DIR}/forceset_res10'
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
    np.testing.assert_allclose(
        estimator1.W_mean_,
        estimator2.W_mean_,
        rtol=1e-06,
        err_msg='Mean work test failed.',
    )
    np.testing.assert_allclose(
        estimator1.W_diss_,
        estimator2.W_diss_,
        rtol=1e-06,
        err_msg='Dissipative work test failed.',
    )
    np.testing.assert_allclose(
        estimator1.dG_,
        estimator2.dG_,
        rtol=1e-06,
        err_msg='Free energy test failed.',
    )
    np.testing.assert_allclose(
        estimator1.friction_,
        estimator2.friction_,
        rtol=1e-06,
        err_msg='Friction test failed.',
    )
    np.testing.assert_allclose(
        estimator1.friction_smooth_,
        estimator2.friction_smooth_,
        rtol=1e-06,
        err_msg='Smoothed friction test failed.',
    )


def test_WorkEstimator(ref_workestimator):
    workestimator_name = f'{TEST_FILE_DIR}/workestimator'
    estimator = storing.load(filename=workestimator_name)
    assert_estimator_equality(estimator, ref_workestimator)


def test_WorkEstimator_res10(ref_workestimator_res10):
    workestimator_name = f'{TEST_FILE_DIR}/workeestimator_res10'
    estimator = storing.load(filename=workestimator_name)
    assert_estimator_equality(estimator, ref_workestimator_res10)


def test_ForceEstimator(ref_forceestimator):
    forceestimator_name = f'{TEST_FILE_DIR}/forceestimator'
    estimator = storing.load(filename=forceestimator_name)
    assert_estimator_equality(estimator, ref_forceestimator)


def test_ForceEstimator_res10(ref_forceestimator_res10):
    forceestimator_name = f'{TEST_FILE_DIR}/forceestimator_res10'
    estimator = storing.load(filename=forceestimator_name)
    assert_estimator_equality(estimator, ref_forceestimator_res10)


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


def test_memory_kernel(ref_forceestimator):
    forceestimator_name = f'{TEST_FILE_DIR}/forceestimator'
    estimator = storing.load(filename=forceestimator_name)

    index = _valid_indices(ref_forceestimator)
    ref_forceestimator.memory_kernel(index=index)

    np.testing.assert_allclose(
        ref_forceestimator.memory_kernel_,
        estimator.memory_kernel_,
        rtol=1e-06,
        err_msg='Memeory kernel using integer as index failed.',
    )


def test_memory_kernel_int_vs_singleton_list_equal(ref_forceestimator):
    fe = ref_forceestimator
    idx = _valid_indices(ref_forceestimator)[1]

    # calculate kenel with int index
    fe.memory_kernel(index=idx)
    kernel_from_int = fe.memory_kernel_
    idx_int = fe.memory_kernel_index_

    # list index
    fe.memory_kernel(index=np.array([idx]))
    kernel_from_list = fe.memory_kernel_
    idx_list = fe.memory_kernel_index_

    # int vs [int] should be identical
    assert idx_int == [idx]
    assert idx_int == idx_list
    np.testing.assert_allclose(kernel_from_int, kernel_from_list, rtol=1e-12)

    # check shape
    assert kernel_from_int.shape == (1, len(fe.force_set.time_))
    assert kernel_from_list.shape == (1, len(fe.force_set.time_))


def test_memory_kernel_ndx_striding(ref_forceestimator):
    fe = ref_forceestimator
    n = len(fe.force_set.time_)
    stride = max(1, n // 10)
    fe.memory_kernel(ndx_striding=stride)

    # Expected indices per implementation:
    expected_indices = np.arange(stride, n - 1, stride, dtype=int)

    # Stored indices and shapes
    assert np.array_equal(fe.memory_kernel_index_, expected_indices)


def test_memory_kernel_deterministic_for_same_input(ref_forceestimator):
    indices = _valid_indices(ref_forceestimator)
    out1 = ref_forceestimator.memory_kernel(index=indices)
    out2 = ref_forceestimator.memory_kernel(index=indices)
    np.testing.assert_allclose(out1, out2, rtol=1e-12)


def test_memory_kernel_errors(ref_forceestimator):
    n = len(ref_forceestimator.force_set.time_)
    idx = max(1, n // 5)
    with pytest.raises(
        ValueError,
        match='Only index or ndx_resolution can be given.'
    ):
        ref_forceestimator.memory_kernel(index=idx, ndx_striding=10)
    with pytest.raises(
        ValueError,
        match='Either index or ndx_resolution must be given.'
    ):
        ref_forceestimator.memory_kernel()
    with pytest.raises(
        ValueError,
        match='Index values must be less than length of data.'
    ):
        ref_forceestimator.memory_kernel(index=np.array([n]))
