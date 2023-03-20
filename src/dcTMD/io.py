# -*- coding: utf-8 -*-
"""Module handling the input/output functions."""
__all__ = ['write_output']

import numpy as np
from beartype import beartype

from dcTMD._typing import (
    Str,
    Str1DArray,
    List,
)


@beartype
def load_pullf(
    pullf_files: Str
) -> (List[Str], Str1DArray):
    """Load filenames from file or glob them from globpattern."""
    try:
        filenames = np.loadtxt(pullf_files, dtype='str')
    except FileNotFoundError as fnf_error:
        print(f'file {fnf_error} using glob.glob({pullf_files})')
        import glob
        filenames = glob.glob(pullf_files)

    if len(filenames) == 0:
        print("No constraint force files found.")
        exit()

    return filenames


@beartype
def write_output(
    out: Str,
    estimator,
    filetype=['.dat', '.npz']
) -> None:
    """
    Takes all calculated quantities and saves them

    Parameters
    ----------
    out :
        Output name. By default f'{out}_N{n_traj}_dG'
    estimator :
        Either a ForceEstimator or WorkEstimator instance.
    filetype:
        Output filetype, either '.dat', '.npz' or ['.dat', '.npz'].
    
    Examples :
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
    if hasattr(estimator, "s_dG_"):
        results_dict = {
            "x": estimator.position_,
            "Wmean": estimator.W_mean_,
            "Wdiss": estimator.W_diss_,
            "dG": estimator.dG_,
            "Gamma": estimator.friction_,
            "s_W_mean": estimator.s_W_mean_,
            "s_W_diss": estimator.s_W_diss_,
            "s_dG": estimator.s_dG_,
        }
    else:
        results_dict = {
            "x": estimator.position_,
            "Wmean": estimator.W_mean_,
            "Wdiss": estimator.W_diss_,
            "dG": estimator.dG_,
            "Gamma": estimator.friction_,
        }

    if '.dat' in filetype:
        results = [(key, item) for key, item in results_dict.items()]
        results = np.asarray(results, dtype=object)
        header = results.T[0]
        arrays = np.vstack(results.T[1])
        print(f'save file {out}_N{n_traj}_dG.dat')
        np.savetxt(
            f'{out}_N{n_traj}_dG.dat',
            arrays.T,
            fmt='%20.8f',  # noqa: WPS323
            header='    '.join(header),
        )
    if '.npz' in filetype:
        print(f'dave file {out}_N{n_traj}_dG.npz')
        np.savez(
            f'{out}_N{n_traj}_dG.npz',
            **results_dict,
        )
