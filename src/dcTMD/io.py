# -*- coding: utf-8 -*-
"""Module handling the input/output functions."""
__all__ = ['write_output']

import numpy as np
from beartype import beartype

from dcTMD._typing import (
    ArrayLikeStr,
    Float1DArray,
    Int,
    Str,
)


@beartype
def _read_testfile(
    file_names: ArrayLikeStr,
    verbose: bool = False,
) -> Float1DArray:
    """Read test force file to determine t and its length.

    Parameters
    ----------
    file_names : list of str or np.ndarray of str
        Contains all force file names
    verbose: bool, optional
        Enables verbose mode.

    Returns
    -------
    t : np.ndarray
        contains time
    """
    if verbose:
        print(f'using {file_names[0]} to initialize arrays')
    test_file = np.loadtxt(file_names[0], comments=['@', '#'])
    return test_file[:, 0]


@beartype
def _load_pullf(
    pullf_glob_pattern: Str,
    pullf_files: Str,
) -> ArrayLikeStr:
    """Return force file names via glob or from file."""
    if pullf_glob_pattern is not None:
        import glob
        files = glob.glob(pullf_glob_pattern)
    if pullf_files is not None:
        files = np.loadtxt(pullf_files, dtype=str)  # TODO reshape input?
    return files


@beartype
def write_output(
    out: Str,
    N_traj: Int,
    **kwargs,
) -> None:
    """Take all calculated quantities and save them."""
    print('writing dG...')
    results = [(key, item) for key, item in kwargs.items() if key != 'errors']
    results = np.asarray(results, dtype=object)
    header = results.T[0]
    arrays = np.vstack(results.T[1])
    if kwargs['errors']:
        header = np.append(
            header, [f's_{name}' for name in header if name != 'x'],
        )
        arrays = np.vstack([arrays, kwargs['errors']])
    np.savetxt(
        f'{out}_{N_traj}_dG.dat',
        arrays.T,
        fmt='%20.8f',  # noqa: WPS323
        header='    '.join(header),
    )
    np.savez(
        f'{out}_{N_traj}_dG.npz',
        kwargs,
    )
