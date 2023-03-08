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

