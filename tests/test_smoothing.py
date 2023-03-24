# -*- coding: utf-8 -*-
"""Tests for the smoothing module.

MIT License
Copyright (c) 2021-2023, Miriam Jäger
All rights reserved.
"""

import pytest
import numpy as np
from os.path import dirname, join
from dcTMD.dcTMD import WorkEstimator
from dcTMD.storing import load
from dcTMD.utils.smoothing import gaussfilter_friction


VELOCITY = 0.001
RESOLUTION = 1
VERBOSE = True
TEMPERATURE = 300
INDICES = np.array([1, 3, 5])
HERE = dirname(__file__)
TESTFILE_DIR = join(HERE, 'testdata')


@pytest.fixture
def ref_workestimator(scope="session"):
    workset_name = join(TESTFILE_DIR, 'workset')
    workset = load(filename=workset_name)
    estimator = WorkEstimator(temperature=TEMPERATURE)
    return estimator.fit(workset)


@pytest.fixture
def ref_friction(scope="session"):
    filename = join(TESTFILE_DIR, 'workset.smoothed_friction')
    return np.loadtxt(filename)


def test_gaussfilter_friction(ref_workestimator, ref_friction):
    smooth_friction = gaussfilter_friction(
        ref_workestimator.friction_, ref_workestimator.position_, 0.1,
    )
    np.testing.assert_array_almost_equal(smooth_friction, ref_friction)
