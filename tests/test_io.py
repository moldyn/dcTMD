# -*- coding: utf-8 -*-
"""
Tests for the io module.

MIT License
Copyright (c) 2021-2022, Victor Tänzel
All rights reserved.
"""

import numpy as np
import os
from os.path import basename, dirname, join
import pytest
from dcTMD import storing, dcTMD
from dcTMD.io import load_pullf, write_output

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
    workset = storing.load(filename=workset_name)
    print(workset)
    estimator = dcTMD.WorkEstimator(temperature=TEMPERATURE)
    estimator.fit(workset)
    return estimator


@pytest.fixture
def ref_forceestimator(scope="session"):
    forceset_name = f'{TEST_FILE_DIR}/forceset'
    forceset = storing.load(filename=forceset_name)
    estimator = dcTMD.ForceEstimator(temperature=TEMPERATURE)
    estimator.fit(forceset)
    return estimator


@pytest.fixture
def filenames_list():
    return [f'{TEST_FILE_DIR}/t_middle_01_pullf.xvg',
            f'{TEST_FILE_DIR}/t_middle_03_pullf.xvg',
            f'{TEST_FILE_DIR}/t_middle_04_pullf.xvg',
            ]


def test_load_pullf(filenames_list):
    # Test case for loading filenames from a file
    print(f'{TEST_FILE_DIR}/pullf_filenames.dat')
    filenames = load_pullf(f'{TEST_FILE_DIR}/pullf_filenames.dat')
    assert len(filenames) == 18
    assert basename(filenames[0]) == basename(filenames_list[0])
    assert basename(filenames[1]) == basename(filenames_list[1])
    assert basename(filenames[2]) == basename(filenames_list[2])
    
    # Test case for globbing filenames
    files = load_pullf(f'{TEST_FILE_DIR}/t_middle_*_pullf.xvg')
    filenames = sorted(files)
    assert len(filenames) == 18
    assert filenames[0] in filenames
    assert filenames[1] in filenames
    assert filenames[2] in filenames
        

def assert_npzfile_equality(npzfile, estimator):
    """Compare saved .npz files and Work/ForceEstimator"""
    np.testing.assert_almost_equal(
        npzfile["x"],
        estimator.position_,
    )
    np.testing.assert_almost_equal(
        npzfile["Wmean"],
        estimator.W_mean_,
    )
    np.testing.assert_almost_equal(
        npzfile["Wdiss"],
        estimator.W_diss_,
    )
    np.testing.assert_almost_equal(
        npzfile["dG"],
        estimator.dG_,
    )
    np.testing.assert_almost_equal(
        npzfile["Gamma"],
        estimator.friction_,
    )


def assert_datfile_equality(datfile, estimator):
    """Compare saved .dat files and Work/ForceEstimator"""
    np.testing.assert_almost_equal(
        datfile[:, 0],
        estimator.position_,
    )
    np.testing.assert_almost_equal(
        datfile[:, 1],
        estimator.W_mean_,
    )
    np.testing.assert_almost_equal(
        datfile[:, 2],
        estimator.W_diss_,
    )
    np.testing.assert_almost_equal(
        datfile[:, 3],
        estimator.dG_,
    )
    np.testing.assert_almost_equal(
        datfile[:, 4],
        estimator.friction_,
    )


def test_write_output_dat_workestimator(
        ref_workestimator,
        tmp_path,
    ):
    out = tmp_path / "test"
    estimator = ref_workestimator
    write_output(str(out), estimator, filetype=['.dat'])
    n_traj = len(estimator.names_)
    test_file = tmp_path / f"test_N{n_traj}_dG.dat"
    assert (test_file).is_file()
    datfile = np.loadtxt(test_file)
    print(datfile.shape)    
    assert_datfile_equality(datfile, estimator)


def test_write_output_dat_forceestimator(
        ref_forceestimator,
        tmp_path,
    ):
    out = tmp_path / "test"
    estimator = ref_forceestimator
    write_output(str(out), estimator, filetype=['.dat'])
    n_traj = len(estimator.names_)
    test_file = tmp_path / f"test_N{n_traj}_dG.dat"
    print(test_file)
    assert (test_file).is_file()


def test_write_output_npz_workestimator(
        ref_workestimator,
        tmp_path,
    ):
    out = tmp_path / "test"
    estimator = ref_workestimator
    write_output(str(out), estimator, filetype=['.npz'])
    n_traj = len(estimator.names_)
    test_file = tmp_path / f"test_N{n_traj}_dG.npz"
    assert (test_file).is_file()
    npzfile = np.load(test_file)
    assert_npzfile_equality(npzfile, estimator)


def test_write_output_npz_forceestimator(
        ref_forceestimator,
        tmp_path,
    ):
    out = tmp_path / "test"
    estimator = ref_forceestimator
    write_output(str(out), estimator, filetype=['.npz'])
    n_traj = len(estimator.names_)
    test_file = tmp_path / f"test_N{n_traj}_dG.npz"
    assert (test_file).is_file()
    npzfile = np.load(test_file)
    assert_npzfile_equality(npzfile, estimator)

