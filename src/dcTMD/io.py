# -*- coding: utf-8 -*-
"""
"""
__all__ = []

import numpy as np
from beartype import beartype
from beartype.typing import Any, Optional
from typing import List


@beartype
def _read_testfile(file_names: List[str],
                   verbose: bool = False,
                   ) -> np.ndarray:
    """ 
    Read test force file to determine t and its length.

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
def load_pullf(pullf_glob_pattern: str,
               pullf_files: str,
               ) -> np.ndarray:
    """Return force file names via glob or from file."""
    if pullf_glob_pattern != None:
        import glob
        files = glob.glob(pullf_glob_pattern)
    if pullf_files != None:
        files = np.loadtxt(pullf_files, dtype=str)  # TODO reshape input?
    return files


@beartype
def write_output(o: str, 
                 N: int, 
                 **kwargs,
                 ) -> None:
    print('writing dG...')
    results = np.asarray([(name, a) for name, a in kwargs.items()
                          if name != 'errors'], dtype=object)
    header = results.T[0]
    arrays = np.vstack(results.T[1])
    if kwargs["errors"]:
        header = np.append(
            header, ['s_' + name for name in header if name != 'x'])
        arrays = np.vstack([arrays, kwargs['errors']])
    np.savetxt(o + "_" + str(N) + '_dG.dat',
               arrays.T,
               fmt='%20.8f',
               header='    '.join(header),
               )
    
