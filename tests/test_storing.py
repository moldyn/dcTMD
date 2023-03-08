# -*- coding: utf-8 -*-
"""
Tests for the storing module.

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
TEMPERATURE = 300
INDICES = np.array([1, 3, 5])
HERE = os.path.dirname(__file__)
TEST_FILE_DIR = os.path.join(HERE, 'testdata')


def filenames(scope="session"):
    files = glob.glob(f'{TEST_FILE_DIR}/t_middle_*_pullf.xvg')
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
def ref_workset(scope="session"):
    workset = dcTMD.storing.WorkSet(
        velocity=VELOCITY,
        resolution=RESOLUTION,
        verbose=VERBOSE,
    )
    workset.fit(filenames())
    return workset


@pytest.fixture
def red_workset(scope="session"):
    workset = dcTMD.storing.WorkSet(
        velocity=VELOCITY,
        resolution=RESOLUTION,
        verbose=VERBOSE,
    )
    workset.fit(filenames())
    red_workset = workset.reduce(INDICES)
    return red_workset


@pytest.fixture
def ref_forceset(scope="session"):
    forceset = dcTMD.storing.ForceSet(
        velocity=VELOCITY,
        resolution=RESOLUTION,
        verbose=VERBOSE,
    )
    forceset.fit(filenames())
    return forceset


def test_save_load_forceset(ref_forceset, tmpdir):
    """Test save/load."""
    filename = str(tmpdir.mkdir('sub').join('load_set_test'))
    dcTMD.storing.save(filename, ref_forceset)
    dataset = dcTMD.storing.load(filename=filename)
    # compare set to reference
    assert_sets_equality(dataset, ref_forceset)


def test_save_load_workset(ref_workset, tmpdir):
    """Test save/load."""
    filename = str(tmpdir.mkdir('sub').join('load_set_test'))
    dcTMD.storing.save(filename, ref_workset)
    dataset = dcTMD.storing.load(filename=filename)
    # compare set to reference
    assert_sets_equality(dataset, ref_workset)


@pytest.mark.parametrize(
    'filename', ['workset'],
)
def test_WorkSet(filename, ref_workset):
    workset = dcTMD.storing.load(
        filename=f'{TEST_FILE_DIR}/{filename}'
    )
    # compare workset to reference
    assert_sets_equality(workset, ref_workset)


@pytest.mark.parametrize(
    'filename', ['forceset'],
)
def test_ForceSet(filename, ref_forceset):
    workset = dcTMD.storing.load(
        filename=f'{TEST_FILE_DIR}/{filename}'
    )
    # compare workset to reference
    assert_sets_equality(workset, ref_forceset)
