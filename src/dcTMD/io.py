# -*- coding: utf-8 -*-
"""Module handling the input/output functions."""
__all__ = ['write_output']

import numpy as np
from beartype import beartype

from dcTMD._typing import (
    Int,
    Str,
)


@beartype
def write_output(
    out: Str,
    n_traj: Int,
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
        f'{out}_{n_traj}_dG.dat',
        arrays.T,
        fmt='%20.8f',  # noqa: WPS323
        header='    '.join(header),
    )
    np.savez(
        f'{out}_{n_traj}_dG.npz',
        kwargs,
    )
