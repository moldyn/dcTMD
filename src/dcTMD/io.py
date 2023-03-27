# -*- coding: utf-8 -*-
# MIT License
# Copyright (c) 2021-2023, Victor Tänzel, Miriam Jäger
# All rights reserved.
"""Submodule handling the input/output operations."""

__all__ = ['write_output']

import glob

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
    filetype=('.dat', '.npz'),
) -> None:
    """Take all calculated quantities and save them.

    Parameters
    ----------
    out :
        Output name. By default f'{out}_N{n_traj}_dG'
    estimator :
        Either a ForceEstimator or WorkEstimator instance.
    filetype:
        Output filetype, either '.dat', '.npz' or ['.dat', '.npz'].

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
    >>> # save results as '.npz' file
    >>> write_output(out, work_estimator, filetype='.npz')
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

    if '.dat' in filetype:
        header = list(results_dict.keys())
        arrays = np.vstack(list(results_dict.values()))
        print(f'save file {out}_dG.dat')
        np.savetxt(
            f'{out}_dG.dat',
            arrays.T,
            fmt='%20.8f',  # noqa: WPS323
            header='    '.join(header),
        )
    if '.npz' in filetype:
        print(f'save file {out}_dG.npz')
        np.savez(
            f'{out}_dG.npz',
            **results_dict,
        )
