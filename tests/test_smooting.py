# -*- coding: utf-8 -*-
"""
Tests for the smoothing module.

MIT License
Copyright (c) 2021-2022, Miriam Jäger
All rights reserved.
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose
import matplotlib.pyplot as plt
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
TEST_FILE_DIR = join(HERE, 'testdata')


@pytest.fixture
def ref_workestimator(scope="session"):
    workset_name = f'{TEST_FILE_DIR}/workset'
    workset = load(filename=workset_name)
    estimator = WorkEstimator(temperature=TEMPERATURE)
    return estimator.fit(workset)


def test_gaussfilter_friction(ref_workestimator):
    smooth_friction_ = gaussfilter_friction(ref_workestimator.friction_,
                                            ref_workestimator.position_,
                                            0.1,
                                            )
    
