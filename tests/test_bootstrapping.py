# -*- coding: utf-8 -*-
"""Tests for the utils smoothing module."""

import pytest
import numpy as np
from dcTMD.utils import bootstrapping
from dcTMD.dcTMD import WorkEstimator
from dcTMD.storing import load
from os.path import dirname, join


VELOCITY = 0.001
RESOLUTION = 1
VERBOSE = True
TEMPERATURE = 300
INDICES = np.array([1, 3, 5])
HERE = dirname(__file__)
TEST_FILE_DIR = join(HERE, 'testdata')


@pytest.fixture
def ref_workestimator(scope="session"):
    workset_name = f'{TEST_FILE_DIR}/workset'
    workset = load(filename=workset_name)
    estimator = WorkEstimator(temperature=TEMPERATURE)
    return estimator.fit(workset)


def test_bootstrapping(ref_workestimator):
    def func(work_set):
        return ref_workestimator.estimate_free_energy(work_set)

    # test 'std' mode
    descriptor = {
        'mode': 'std',
        'n_resamples': 100,
        'seed': 42,
    }
    s_quantity, quantity_resampled = bootstrapping(ref_workestimator,
                                                   func,
                                                   descriptor)
    # Check if the output is of the expected shape
    assert s_quantity.shape == (1, 3, len(ref_workestimator.position_))
    assert quantity_resampled.shape == (descriptor['n_resamples'],
                                        3,
                                        len(ref_workestimator.position_))

    # Check if the output is of the expected type
    assert isinstance(s_quantity, np.ndarray)
    assert isinstance(quantity_resampled, np.ndarray)

    # test confidence interval mode
    descriptor = {
        'mode': 0.5,
        'n_resamples': 100,
        'seed': 42,
    }
    s_quantity, quantity_resampled = bootstrapping(ref_workestimator,
                                                   func,
                                                   descriptor)
    # Check if the output is of the expected shape
    assert s_quantity.shape == (1, 2, 3, len(ref_workestimator.position_))
    # Check if the output is of the expected type
    assert isinstance(s_quantity, np.ndarray)
    assert isinstance(quantity_resampled, np.ndarray)
