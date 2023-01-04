# -*- coding: utf-8 -*-
"""
Tests for the work module.

MIT License
Copyright (c) 2021-2022, Victor Tänzel
All rights reserved.
"""

import dcTMD
import pytest
import os
import glob
import numpy as np

VELOCITY = 0.001
RESOLUTION = 1
VERBOSE = True
INDICES = np.array([1, 3, 5])
HERE = os.path.dirname(__file__)


def filenames():
    path = os.path.join(HERE, 'testdata')
    files = glob.glob(f'{path}/t_middle_*_pullf.xvg')
    return sorted(files)


def assert_sets_equality(set1, set2):
    """Compare to WorkSet or ForceSet instances via asserts."""
    assert type(set1) is type(set2)  # noqa: WPS516
    assert set1.velocity == set2.velocity
    assert set1.resolution == set2.resolution
    assert set1.verbose == set2.verbose
    assert (set1.names_ == set2.names_).all()
    np.testing.assert_almost_equal(
        set1.work_,
        set2.work_,
    )
    np.testing.assert_almost_equal(
        set1.time_,
        set2.time_,
    )
    np.testing.assert_almost_equal(
        set1.position_,
        set2.position_,
    )


@pytest.fixture
def ref_workset():
    workset = dcTMD.storing.WorkSet(
        velocity=VELOCITY,
        resolution=RESOLUTION,
        verbose=VERBOSE,
    )
    workset.fit(filenames())
    return workset


@pytest.fixture
def red_workset():
    workset = dcTMD.storing.WorkSet(
        velocity=VELOCITY,
        resolution=RESOLUTION,
        verbose=VERBOSE,
    )
    workset.fit(filenames())
    red_workset = workset.reduce(INDICES)
    return red_workset


@pytest.mark.parametrize(
    'filename', ['testdata/workset'],
)
def test_load(filename, ref_workset):
    workset = dcTMD.storing.load(
        filename=f'{HERE}/{filename}'
    )
    # compare workset to reference
    assert_sets_equality(workset, ref_workset)


@pytest.mark.parametrize(
    'filename', ['testdata/red_workset'],
)
def test_WorkSet(filename, red_workset):
    workset = dcTMD.storing.load(
        filename=f'{HERE}/{filename}'
    )
    # compare workset to reference
    assert_sets_equality(workset, red_workset)
