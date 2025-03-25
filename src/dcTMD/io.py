# -*- coding: utf-8 -*-
# MIT License
# Copyright (c) 2021-2023, Miriam Jäger
# All rights reserved.
"""Submodule handling the input/output operations."""

__all__ = ['write_output']

import glob
from pathlib import Path
import numpy as np
from beartype import beartype

from dcTMD._typing import (
    Str,
    Str1DArray,
    List,
)


@beartype
def load_pullf(pullf_files: Str) -> (List[Str], Str1DArray):
    """Load filenames and resturns them as list.

    Filenames can be taken from file or glob them from globpattern.

    Parameters
    ----------
    pullf_files :
        file which contains `pullf` filenames or globpattern

    Returns
    -------
    filenames:
        list/array of filenames

    Examples
    --------
    >>> from dcTMD.io import load_pullf
    >>> # load filenames form file
    >>> filenames = load_pullf('pullf_files.txt')
    >>> # load filenames with glob pattern
    >>> filenames = load_pullf('data/*.pullf.xvg')
    """
    try:
        filenames = np.loadtxt(pullf_files, dtype='str')
    except FileNotFoundError as fnf_error:
        print(f'file {fnf_error} using glob.glob({pullf_files})')
        filenames = glob.glob(pullf_files)

    if not len(filenames):
        raise ValueError('No constraint force files found.')

    return filenames


@beartype
def write_output(
    out: Str,
    estimator,
    filetype=('dat', 'npz'),
) -> None:
    """Take all calculated quantities and save them.

    Parameters
    ----------
    out :
        Output name. By default f'{out}_N{n_traj}'
    estimator :
        Either a ForceEstimator or WorkEstimator instance.
    filetype:
        Output filetype, either 'dat', 'npz' or both ('dat', 'npz').

    Examples
    --------
    >>> from dcTMD.storing import load
    >>> from dcTMD.io import write_output
    >>> from dcTMD.dcTMD import WorkEstimator
    >>> # Save the results from WorkEstimator
    >>> # (ForceEstimator works similarly)
    >>> # calculate dcTMD results from workset
    >>> work = load('my_work_set')
    >>> work_estimator = WorkEstimator(temperature=290.15)
    >>> work_estimator.fit(work)
    >>> out = 'my_dcTMD_results'
    >>> # save results as 'npz' file
    >>> write_output(out, work_estimator, filetype='npz')
    >>> # results saves as 'my_dcTMD_results_N100.npz'
    >>> # save results as 'dat' file
    >>> write_output(out, work_estimator, filetype='dat')
    >>> # save results as 'dat' and 'npz' file
    >>> write_output(out, work_estimator, filetype=('dat', 'npz'))
    """
    n_traj = len(estimator.names_)
    out = f'{out}_N{n_traj}'

    results_dict = {
        'x': estimator.position_,
        'Wmean': estimator.W_mean_,
        'Wdiss': estimator.W_diss_,
        'dG': estimator.dG_,
        'Gamma': estimator.friction_,
    }
    if hasattr(estimator, 's_dG_'):
        results_dict.update({
            's_W_mean': estimator.s_W_mean_,
            's_W_diss': estimator.s_W_diss_,
            's_dG': estimator.s_dG_,
        })
    if hasattr(estimator, 'friction_smooth_'):
        results_dict['Gamma_smooth'] = estimator.friction_smooth_
    if hasattr(estimator, 's_friction_'):
        results_dict['s_Gamma'] = estimator.s_friction_

    if 'dat' in filetype:
        header = list(results_dict.keys())
        arrays = np.vstack(list(results_dict.values()))
        print(f'save file {out}.dat')
        np.savetxt(
            f'{out}.dat',
            arrays.T,
            fmt='%20.8f',  # noqa: WPS323
            header='    '.join(header),
        )
    if 'npz' in filetype:
        print(f'save file {out}.npz')
        np.savez(
            f'{out}.npz',
            **results_dict,
        )


@beartype
def load_output(filepath: str) -> dict[str, np.ndarray]:
    """
    Load file produced by `write_output`.

    Parameters
    ----------
    filepath : str
        Path to the .dat or .npz file to be loaded.

    Returns
    -------
    Dict[str, Any]
        A dictionary containing the loaded data. Keys correspond to the
        quantities saved in the .npz file, such as 'x', 'Wmean', 'Wdiss',
        'dG', 'Gamma', and optionally 's_W_mean', 's_W_diss', 's_dG',
        'Gamma_smooth', 's_Gamma' if they were present during saving.

    Examples
    --------
    >>> from dcTMD.storing import load
    >>> from dcTMD.io import write_output, load_output
    >>> from dcTMD.dcTMD import WorkEstimator
    >>> # Save the results from WorkEstimator
    >>> work = load('my_work_set')
    >>> work_estimator = WorkEstimator(temperature=290.15)
    >>> work_estimator.fit(work)
    >>> # save results as '.npz' file
    >>> out = 'my_dcTMD_results'
    >>> write_output(out, work_estimator, filetype='npz')
    >>> # load results
    >>> results = load_output('my_dcTMD_results_N100.npz')
    >>> positions = results['x']
    >>> mean_work = results['Wmean']
    >>> dG = results['dG']
    """
    if not Path(filepath).is_file():
        raise FileNotFoundError(
            f'The file "{filepath}" does not exist.',
        )
    if '.npz' in filepath:
        with np.load(filepath) as npzfile:
            res_dict = {key: npzfile[key] for key in npzfile.files}
        print(f'Loaded data from {filepath}')
        return res_dict
    elif '.dat' in filepath:
        with open(filepath, 'r') as datfile:
            # Read the header line, which starts with '#'
            for line in datfile:
                if line.startswith('#', 0, 5):
                    headers = line.split()
                    headers.remove('#')
                    break
        dctmd_results = np.loadtxt(filepath)
        res_dict = {}
        for idx, key in enumerate(headers):
            res_dict[key] = dctmd_results[:, idx]
        print(f'Loaded data from {filepath}')
        return res_dict
    else:
        print('Could not load file.')
        print(f'{filepath} needs to be .dat or .npz file.')
