import os
import numpy as np
import pytest
import shutil
from dcTMD.featureset import FeatureSet
from os.path import dirname, join
from pathlib import Path

HERE = dirname(__file__)
TEST_FILE_DIR = join(HERE, 'testdata')


# test 2D case
@pytest.fixture
def filenames_2d():
    fname = os.path.join(TEST_FILE_DIR, "feature_2D_filenamesname.txt")
    names = np.loadtxt(fname, dtype=str).tolist()
    fnames = []
    for f in names:
        fnames.append(f'{TEST_FILE_DIR}/{f}')

    return fnames


@pytest.fixture
def feature_array_3d():
    fname = os.path.join(TEST_FILE_DIR, "feature_array_3D.npy")
    return np.load(fname)


@pytest.fixture
def traj_feature_set3d(filenames_2d):
    return FeatureSet(filenames=filenames_2d, verbose=False)


def test_fill_array(traj_feature_set3d, feature_array_3d):
    """
    Test that fill_array loads data into an array with
    the expected shape and values.
    """
    result_array = traj_feature_set3d.fill_array()
    assert result_array.shape == feature_array_3d.shape
    np.testing.assert_allclose(result_array, feature_array_3d)


# test 2D case
@pytest.fixture
def filenames_1d():
    fname = os.path.join(TEST_FILE_DIR, "feature_1D_filenamesname.txt")
    names = np.loadtxt(fname, dtype=str).tolist()
    fnames = []
    for f in names:
        fnames.append(f'{TEST_FILE_DIR}/{f}')

    return fnames


@pytest.fixture
def feature_array_2d():
    fname = os.path.join(TEST_FILE_DIR, "feature_array_2D.npy")
    return np.load(fname)


@pytest.fixture
def traj_feature_set2d(filenames_1d):
    return FeatureSet(filenames=filenames_1d, verbose=False)


def test_fill_array2d(traj_feature_set2d, feature_array_2d):
    """
    Test that fill_array loads data into an array with
    the expected shape and values.
    """
    result_array = traj_feature_set2d.fill_array()
    assert result_array.shape == feature_array_2d.shape
    np.testing.assert_allclose(result_array, feature_array_2d)


# test the other functions
def test_get_filenames(filenames_2d):
    """
    Test that _get_filenames works correctly with filenameprefix and wildcard.
    """
    wildcard = os.path.join(TEST_FILE_DIR, "feature_2D_{}.txt")
    tfs = FeatureSet(
        filenameprefix=filenames_2d,
        wildcard=wildcard,
        verbose=False
    )
    # Check that each file found exists.
    for fname in tfs.filenames:
        assert os.path.exists(fname), f"File {fname} does not exist."


def test_invalid_initialization():
    """
    Test that initializing without either filenames or
    filenameprefix (with wildcard) raises ValueError.
    """
    with pytest.raises(ValueError):
        _ = FeatureSet()


def test_shape_mismatch():
    """
    Test that if one file has a shape mismatch, it
    gets skipped and the corresponding slice remains zeros.
    """
    # Create a temporary directory with two files:
    d = Path(TEST_FILE_DIR) / "temp_data"
    d.mkdir(exist_ok=True)
    file1 = os.path.join(d, "file_01.txt")
    file2 = os.path.join(d, "file_02.txt")

    arr1 = np.random.rand(100, 10)
    np.savetxt(file1, arr1)

    arr2 = np.random.rand(50, 10)
    np.savetxt(file2, arr2)

    tfs = FeatureSet(filenames=[str(file1), str(file2)], verbose=False)
    result_array = tfs.fill_array()

    # The expected shape is (2, 100, 10), since the first file sets the shape.
    assert result_array.shape == (2, 100, 10)
    # The data from file1 should match, while the slice for file2 is zero.
    np.testing.assert_allclose(result_array[0], arr1)
    np.testing.assert_allclose(result_array[1], np.zeros((100, 10)))
    # Remove temporary directory
    shutil.rmtree(d)


def test_read_testfile(tmp_path):
    """
    Test that _read_testfile correctly sets self.fileshape
    based on the contents of the first file.
    """
    data = np.arange(12).reshape(3, 4)
    d = Path(TEST_FILE_DIR) / "temp_data"
    d.mkdir(exist_ok=True)

    test_file = d / "test_file.txt"
    np.savetxt(str(test_file), data)

    # Initialize the FeatureSet instance with test_file.
    tfs = FeatureSet(filenames=[str(test_file)], verbose=False)
    tfs._read_testfile()
    assert tfs.fileshape == data.shape
    # Remove temporary directory
    shutil.rmtree(d)


def test_file_not_found():
    """
    Test that attempting to load a non-existent file raises a ValueError.
    """
    non_existent = os.path.join(TEST_FILE_DIR, "non_existent.txt")
    tfs = FeatureSet(filenames=[str(non_existent)], verbose=False)
    with pytest.raises(FileNotFoundError):
        _ = tfs.fill_array()
