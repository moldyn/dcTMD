# -*- coding: utf-8 -*-
"""
Functions used for bootstrapping.

MIT License
Copyright (c) 2022, Victor Tänzel
All rights reserved.
"""

__all__ = ['bootstrapping']

import numpy as np
from beartype import beartype
from beartype.typing import Union, Callable, Optional
from sklearn.base import BaseEstimator, TransformerMixin

from dcTMD._typing import (
    Int,
    Float,
    Float2DArray,
    Float3DArray,
    StrStd,
    NumInRange0to1,
)


@beartype
def _bootstrap_reducer(
    descriptor: dict,
    *args: Union[Float2DArray, Float3DArray],
):
    """
    Perform the last bootstrap step.

    Reduces sequences of a resampled statistics by calculating standard
    deviations or confidence intervals, depending on the 'mode' attribute of
    the given object obj. If mode is a real number in [0, 1), confidence
    intervals will be computed. If it is the string 'std', then the standard
    deviation is calculated.

    Parameters
    ----------
    obj :
        Instance of WorkEstimator which provides its 'mode' attribute.
    *args :
        Resampled statistics, to be reduced.

    Returns
    -------
    tuple :
        Tuple of reduced statistics.
    """
    if descriptor['mode'] == 'std':
        def reducer(x):
            return np.std(x, axis=0)
    else:
        # Implement the confidence intervals here.
        def reducer(x):
            pass
    return tuple(reducer(arg) for arg in args)


@beartype
def bootstrapping(
    obj,
    func: Callable,
    descriptor: dict,
    verbose: bool = False,
):
    """Bootstrapping."""
    import tqdm
    # Set up an array in which the resampled quantities are saved.
    # How many quantities are returned by func? If a tuple is returned,
    # use its length. Else, only one quantitity is returned, thus
    # len_of_return is 1.
    N_traj, length_data = np.shape(obj.work_set.work_)
    probe_return = func(obj.work_set.work_)
    if type(probe_return) == tuple:
        len_of_return = len(probe_return)
    else:
        len_of_return = 1
    quantity_resampled = np.empty((
        descriptor['n_resamples'],
        len_of_return,
        length_data,
    ))
    # Initialize RNG.
    seed = descriptor['seed']
    rng = np.random.default_rng(seed)
    for ind in tqdm.tqdm(
        range(descriptor['n_resamples']),
        desc='Bootstrapping progress',
    ):
        # Draw random work time traces
        random_indices = rng.integers(0, N_traj, N_traj)
        work_set_resampled = obj.work_set.work_[random_indices]

        # Calculate and save the relevant statistic
        quantity = func(work_set_resampled)
        quantity_resampled[ind] = quantity
    # There are now boostrapped quantities in the '_resampled' variables.
    # We are interested in the element-wise distributions and thus
    # calculate (1) the standard distribution of the resampled quantity
    # at all points or (2) confidence intervals.
    if verbose:
        print('Finished resampling, starting reduction.')
    s_quantity = _bootstrap_reducer(
        descriptor,
        quantity_resampled,
    )
    s_quantity = np.array(s_quantity)
    # The distributions of the '_resampled' variables must be inspected
    # and are thus also returned.
    return s_quantity, quantity_resampled
