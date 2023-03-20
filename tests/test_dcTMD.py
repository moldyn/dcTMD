# -*- coding: utf-8 -*-
"""
Tests for the clustering module.

MIT License
Copyright (c) 2021-2022, Victor Tänzel
All rights reserved.
"""

import numpy as np
import pytest
import os
import dcTMD


VELOCITY = 0.001
RESOLUTION = 1
VERBOSE = True
TEMPERATURE = 300
INDICES = np.array([1, 3, 5])
HERE = os.path.dirname(__file__)


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
    )


@pytest.fixture
def ref_workestimator(scope="session"):
    estimator = dcTMD.dcTMD.WorkEstimator(temperature=TEMPERATURE)

    estimator.fit(ref_workset)
    return estimator


@pytest.fixture
def ref_forceestimator(scope="session"):
    estimator = dcTMD.dcTMD.ForceEstimator(temperature=TEMPERATURE)

    estimator.fit(ref_forceset)
    return estimator



import numpy as np
import pytest
from dcTMD.dcTMD import ForceSet, ForceEstimator


@pytest.fixture
def force_set():
    time = np.array([0.0, 0.1, 0.2, 0.3, 0.4])
    position = np.array([0.0, 0.1, 0.2, 0.3, 0.4])
    velocity = 1.0
    force = np.array([[0.0, 0.1, 0.2, 0.3, 0.4],
                      [1.0, 1.1, 1.2, 1.3, 1.4],
                      [2.0, 2.1, 2.2, 2.3, 2.4]])
    resolution = 2
    return ForceSet(time=time, position=position, velocity=velocity,
                    force=force, resolution=resolution)


def test_force_set_init(force_set):
    assert np.allclose(force_set.time_, [0.0, 0.2, 0.4])
    assert np.allclose(force_set.position_, [0.0, 0.2, 0.4])
    assert force_set.velocity == 1.0
    assert np.allclose(force_set.force_,
                       [[0.0, 0.2, 0.4], [1.0, 1.2, 1.4], [2.0, 2.2, 2.4]])
    assert force_set.resolution == 2


def test_force_set_slice(force_set):
    sliced = force_set[1:3]
    assert np.allclose(sliced.time_, [0.2])
    assert np.allclose(sliced.position_, [0.2])
    assert sliced.velocity == 1.0
    assert np.allclose(sliced.force_, [[1.2], [2.2]])
    assert sliced.resolution == 2


@pytest.fixture
def force_estimator(force_set):
    return ForceEstimator(temperature=300.0)


def test_force_estimator_fit(force_estimator, force_set):
    fitted_estimator = force_estimator.fit(force_set)
    assert fitted_estimator is force_estimator
    assert np.allclose(force_estimator.dG_, [0.0, -0.1, -0.2])
    assert np.allclose(force_estimator.friction_,
                       [0.0, -0.01512851, -0.03025702])


import pytest
import numpy as np
from dcTMD.dcTMD import WorkEstimator
from dcTMD.storing import load
from dcTMD._typing import Float1DArray

@pytest.fixture
def work_estimator():
    work = load('my_work_set')
    work_estimator = WorkEstimator(temperature=290.15)
    work_estimator.fit(work)
    return work_estimator

def test_estimate_free_energy(work_estimator):
    W_mean, W_diss, dG = work_estimator.estimate_free_energy()
    assert isinstance(W_mean, np.ndarray)
    assert isinstance(W_diss, np.ndarray)
    assert isinstance(dG, np.ndarray)
    assert len(W_mean) == len(W_diss) == len(dG) == work_estimator.work_set.n_windows

def test_estimate_free_energy_errors(work_estimator):
    n_resamples = 100
    mode = 'resample'
    seed = 42
    W_mean_resampled, W_diss_resampled, dG_resampled = work_estimator.estimate_free_energy_errors(n_resamples, mode, seed)
    assert isinstance(W_mean_resampled, np.ndarray)
    assert isinstance(W_diss_resampled, np.ndarray)
    assert isinstance(dG_resampled, np.ndarray)
    assert len(W_mean_resampled) == len(W_diss_resampled) == len(dG_resampled) == work_estimator.work_set.n_windows
    assert all(work_estimator.W_mean_resampled_.shape == W_mean_resampled.shape)
    assert all(work_estimator.W_diss_resampled_.shape == W_diss_resampled.shape)
    assert all(work_estimator.dG_resampled_.shape == dG_resampled.shape)

def test_transform(work_estimator):
    dG, friction = work_estimator.transform(None)
    assert isinstance(dG, np.ndarray)
    assert isinstance(friction, np.ndarray)
    assert len(dG) == len(friction) == work_estimator.work_set.n_windows

def test_fit(work_estimator):
    assert work_estimator.W_mean_.shape == (work_estimator.work_set.n_windows,)
    assert work_estimator.W_diss_.shape == (work_estimator.work_set.n_windows,)
    assert work_estimator.dG_.shape == (work_estimator.work_set.n_windows,)
    assert isinstance(work_estimator.friction_, np.ndarray)
    assert work_estimator.s_W_mean_.shape == (work_estimator.work_set.n_windows,)
    assert work_estimator.s_W_diss_.shape == (work_estimator.work_set.n_windows,)
    assert work_estimator.s_dG_.shape == (work_estimator.work_set.n_windows,)
    assert isinstance(work_estimator.W_mean_resampled_, np.ndarray)
    assert isinstance(work_estimator.W_diss_resampled_, np.ndarray)
    assert isinstance(work_estimator.dG_resampled_, np.ndarray)
    assert all(work_estimator.W_mean_resampled_.shape == work_estimator.W_mean_.shape)
    assert all(work_estimator.W_diss_resampled_.shape == work_estimator.W_diss_.shape)
    assert all(work_estimator.dG_resampled_.shape == work_estimator.dG_.shape)
