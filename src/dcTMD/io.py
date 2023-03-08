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
    # ToDo: add exception for Union typing
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
    dataset,
    estimator,
    filetype=['.dat', '.npz']
) -> None:
    """Take all calculated quantities and save them."""
    n_traj = len(dataset.names_)
    if hasattr(estimator, "s_dG_"):
        results_dict = {
            "x": dataset.position_,
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
            "x": dataset.position_,
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
        np.savetxt(
            f'{out}_N{n_traj}_dG.dat',
            arrays.T,
            fmt='%20.8f',  # noqa: WPS323
            header='    '.join(header),
        )
    if '.npz' in filetype:
        np.savez(
            f'{out}_N{n_traj}_dG.npz',
            **results_dict,
        )
