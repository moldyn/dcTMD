# -*- coding: utf-8 -*-
"""Tests for the io module."""

import numpy as np
from os.path import basename, dirname, join
import pytest
from dcTMD import storing, dcTMD
from dcTMD.io import load_pullf, write_output, load_output


VELOCITY = 0.001
RESOLUTION = 1
VERBOSE = True
TEMPERATURE = 300
INDICES = np.array([1, 3, 5])
HERE = dirname(__file__)
TEST_FILE_DIR = join(HERE, 'testdata')


@pytest.fixture
def ref_workestimator(scope='session'):
    workset_name = f'{TEST_FILE_DIR}/workset'
    workset = storing.load(filename=workset_name)
    print(workset)
    estimator = dcTMD.WorkEstimator(temperature=TEMPERATURE)
    estimator.fit(workset)
    return estimator


@pytest.fixture
def ref_forceestimator(scope='session'):
    forceset_name = f'{TEST_FILE_DIR}/forceset'
    forceset = storing.load(filename=forceset_name)
    estimator = dcTMD.ForceEstimator(temperature=TEMPERATURE)
    estimator.fit(forceset)
    return estimator


@pytest.fixture
def filenames_list():
    return [f'{TEST_FILE_DIR}/t_middle_01_pullf.xvg',
            f'{TEST_FILE_DIR}/t_middle_03_pullf.xvg',
            f'{TEST_FILE_DIR}/t_middle_04_pullf.xvg',]


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
    np.testing.assert_allclose(
        npzfile['x'],
        estimator.position_,
        rtol=1e-06,
    )
    np.testing.assert_allclose(
        npzfile['Wmean'],
        estimator.W_mean_,
        rtol=1e-06,
    )
    np.testing.assert_allclose(
        npzfile['Wdiss'],
        estimator.W_diss_,
        rtol=1e-05,
    )
    np.testing.assert_allclose(
        npzfile['dG'],
        estimator.dG_,
        rtol=1e-06,
    )
    np.testing.assert_allclose(
        npzfile['Gamma'],
        estimator.friction_,
        rtol=1e-06,
    )


def assert_datfile_equality(datfile, estimator):
    """Compare saved .dat files and Work/ForceEstimator"""
    np.testing.assert_allclose(
        datfile[:, 0],
        estimator.position_,
        rtol=1e-06,
    )
    np.testing.assert_allclose(
        datfile[:, 1],
        estimator.W_mean_,
        rtol=1e-06,
    )
    np.testing.assert_allclose(
        datfile[:, 2],
        estimator.W_diss_,
        rtol=1e-05,
    )
    np.testing.assert_allclose(
        datfile[:, 3],
        estimator.dG_,
        rtol=1e-06,
    )
    np.testing.assert_allclose(
        datfile[:, 4],
        estimator.friction_,
        rtol=1e-06,
    )


def test_write_output_dat_workestimator(ref_workestimator, tmp_path):
    out = tmp_path / 'test'
    estimator = ref_workestimator
    write_output(str(out), estimator, filetype=['dat'])
    n_traj = len(estimator.names_)
    test_file = tmp_path / f'test_N{n_traj}.dat'
    assert (test_file).is_file()
    datfile = np.loadtxt(test_file)
    assert_datfile_equality(datfile, estimator)


def test_write_output_dat_forceestimator(ref_forceestimator, tmp_path):
    out = tmp_path / 'test'
    estimator = ref_forceestimator
    write_output(str(out), estimator, filetype=['dat'])
    n_traj = len(estimator.names_)
    test_file = tmp_path / f'test_N{n_traj}.dat'
    assert (test_file).is_file()
    datfile = np.loadtxt(test_file)
    assert_datfile_equality(datfile, estimator)


def test_write_output_npz_workestimator(ref_workestimator, tmp_path):
    out = tmp_path / 'test'
    estimator = ref_workestimator
    write_output(str(out), estimator, filetype=['npz'])
    n_traj = len(estimator.names_)
    test_file = tmp_path / f'test_N{n_traj}.npz'
    assert (test_file).is_file()
    npzfile = np.load(test_file)
    assert_npzfile_equality(npzfile, estimator)


def test_write_output_npz_forceestimator(ref_forceestimator, tmp_path):
    out = tmp_path / 'test'
    estimator = ref_forceestimator
    write_output(str(out), estimator, filetype=['npz'])
    n_traj = len(estimator.names_)
    test_file = tmp_path / f'test_N{n_traj}.npz'
    assert (test_file).is_file()
    npzfile = np.load(test_file)
    assert_npzfile_equality(npzfile, estimator)


def test_load_output_dat_workestimator_fromfile(ref_workestimator):
    filepath = f'{TEST_FILE_DIR}/workeestimator_dcTMDresults_N18.dat'
    estimator = ref_workestimator
    res_dict = load_output(filepath)
    assert_npzfile_equality(res_dict, estimator)


def test_load_output_npz_workestimator_fromfile(ref_workestimator):
    filepath = f'{TEST_FILE_DIR}/workeestimator_dcTMDresults_N18.npz'
    estimator = ref_workestimator
    res_dict = load_output(filepath)
    assert_npzfile_equality(res_dict, estimator)


def test_load_output_dat_workestimator(ref_workestimator, tmp_path):
    out = tmp_path / 'test'
    estimator = ref_workestimator
    write_output(str(out), estimator, filetype=['dat'])
    n_traj = len(estimator.names_)
    test_file = tmp_path / f'test_N{n_traj}.dat'
    res_dict = load_output(str(test_file))
    assert_npzfile_equality(res_dict, estimator)


def test_load_output_dat_forceestimator(ref_forceestimator, tmp_path):
    out = tmp_path / 'test'
    estimator = ref_forceestimator
    write_output(str(out), estimator, filetype=['dat'])
    n_traj = len(estimator.names_)
    test_file = tmp_path / f'test_N{n_traj}.dat'
    res_dict = load_output(str(test_file))
    assert_npzfile_equality(res_dict, estimator)


def test_load_output_npz_workestimator(ref_workestimator, tmp_path):
    out = tmp_path / 'test'
    estimator = ref_workestimator
    write_output(str(out), estimator, filetype=['npz'])
    n_traj = len(estimator.names_)
    test_file = tmp_path / f'test_N{n_traj}.npz'
    res_dict = load_output(str(test_file))
    print(res_dict)
    assert_npzfile_equality(res_dict, estimator)


def test_load_output_npz_forceestimator(ref_forceestimator, tmp_path):
    out = tmp_path / 'test'
    estimator = ref_forceestimator
    write_output(str(out), estimator, filetype=['npz'])
    n_traj = len(estimator.names_)
    test_file = tmp_path / f'test_N{n_traj}.npz'
    res_dict = load_output(str(test_file))
    assert_npzfile_equality(res_dict, estimator)
