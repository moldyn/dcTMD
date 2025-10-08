# -*- coding: utf-8 -*-
# MIT License
# Copyright (c) 2022, Victor Tänzel, Miriam Jäger
# All rights reserved.
"""Functions used for bootstrapping."""

__all__ = ['bootstrapping']

import sys
import tqdm
import numpy as np
from beartype import beartype
from beartype.typing import Callable, Optional, Union
from dcTMD._typing import (
    Float2DArray,
    Float3DArray,
)


@beartype
def _bootstrap_reducer(
    descriptor: dict,
    *args: Union[Float2DArray, Float3DArray],
):
    """
    Perform the last bootstrap step.

    Reduces arrays of resampled statistics by calculating standard deviations
    or confidence intervals, depending on the `mode` attribute of the given
    object obj. If mode is a real number in [0, 1), confidence intervals will
    be computed. If it is the string 'std', then the standard deviation is
    calculated.

    Parameters
    ----------
    descriptor :
        Dictionary of the estimator which provides its `mode` attribute.
    args :
        Resampled statistics, to be reduced.

    Returns
    -------
    tuple :
        Tuple of reduced statistics.
    """
    if descriptor['mode'] == 'std':
        # Calculate the standard deviation of the resampled quantity.
        def reducer(resampled_quantity):  # noqa: WPS440
            return np.std(resampled_quantity, axis=0)
    else:
        def reducer(resampled_quantity):
            confidence_level = descriptor['mode']
            # Calculate the lower and upper bounds of the confidence interval
            lower_bound = np.percentile(
                resampled_quantity,
                (1 - confidence_level) / 2 * 100,
                axis=0,
            )
            upper_bound = np.percentile(
                resampled_quantity,
                (1 + confidence_level) / 2 * 100,
                axis=0,
            )
            # Return and the confidence interval as a tuple
            return lower_bound, upper_bound
    reduced_args = tuple(reducer(arg) for arg in args)
    return np.asarray(reduced_args)


@beartype
def _bootstrap_resampling(  # noqa: WPS210
    estimator,
    func: Callable,
    descriptor: dict,
):
    """Perform the resampling step of the bootstrap algorithm."""
    # Set up an array in which the resampled quantities are saved.
    # How many quantities are returned by func? If a tuple is returned,
    # use its length. Else, only one quantitity is returned.
    n_traj, length_data = np.shape(estimator.work_set.work_)
    probe_return = func(estimator.work_set.work_)
    if isinstance(probe_return, tuple):
        len_of_return = len(probe_return)
    else:
        len_of_return = 1

    quantity_resampled = np.empty((
        descriptor['n_resamples'],
        len_of_return,
        length_data,
    ))
    # Initialize RNG.
    rng = np.random.default_rng(descriptor['seed'])
    for idx in tqdm.tqdm(
        range(descriptor['n_resamples']),
        desc='Bootstrapping progress',
    ):
        # Draw random work time traces
        random_indices = rng.integers(0, n_traj, n_traj)
        work_set_resampled = estimator.work_set.work_[random_indices]

        # Calculate and save the relevant statistic
        quantity = func(work_set_resampled)
        quantity_resampled[idx] = quantity
    return quantity_resampled


@beartype
def bootstrapping(
    estimator,
    func: Callable,
    descriptor: dict,
    verbose: Optional[bool] = False,
):
    """
    Perform a bootstrapping error analysis.

    The bootstrapping error analysis is performed using a given function `func`
    by drawing random work trajectories from the `estimator`'s WorkSet
    instance with replacement. The quantity of interest is then calculated for
    the new sample and stored in `quantity_resampled`. This process is repeated
    `n_resamples` times, which is a key of the `descriptor` dictionary. The
    random number generator can be fed a `seed`, the second key of the
    `descriptor` which is optional. Thirdly, a `mode` must be in the
    `descriptor`, which can either be the string 'std' for a standard
    distribution of the resampled quantity or a number in the interval [0, 1)
    which yields confidence intervals instead.

    Parameters
    ----------
    estimator :
        Instance of a WorkEstimator.
    func :
        Function which takes a WorkSet instance as single argument and returns
        the the quantitity for which the bootstrapping error analysis is
        performed.
    descriptor :
        Dictionary of the estimator which provides `mode`, `n_resampled`
        and `seed` as keys.
    verbose :
        Enables verbose mode.

    Returns
    -------
    s_quantity :
        Estimated error for the quantity returned by `func`.
    quantity_resampled :
        Quantities returned by `func` for the resampled work trajectories.
    """

    """
    n_traj, length_data = np.shape(estimator.work_set.work_)
    probe_return = func(estimator.work_set.work_)
    if isinstance(probe_return, tuple):
        len_of_return = len(probe_return)
    else:
        len_of_return = 1

    quantity_resampled = np.empty((
        descriptor['n_resamples'],
        len_of_return,
        length_data,
    ))
    # Initialize RNG.
    rng = np.random.default_rng(descriptor['seed'])
    for idx in tqdm.tqdm(
        range(descriptor['n_resamples']),
        desc='Bootstrapping progress',
    ):
        # Draw random work time traces
        random_indices = rng.integers(0, n_traj, n_traj)
        work_set_resampled = estimator.work_set.work_[random_indices]

        # Calculate and save the relevant statistic
        quantity = func(work_set_resampled)
        quantity_resampled[idx] = quantity
    """
    quantity_resampled = _bootstrap_resampling(
        estimator,
        func,
        descriptor,
    )
    # There are now bootstrapped quantities in the '_resampled' variables.
    # We are interested in the element-wise distributions and thus
    # calculate (1) the standard distribution of the resampled quantity
    # at all points or (2) confidence intervals.
    if verbose:
        sys.stdout.write('Finished resampling, starting reduction.')
    s_quantity = _bootstrap_reducer(
        descriptor,
        quantity_resampled,
    )
    # The distributions of the '_resampled' variables must be inspected
    # and are thus also returned.
    return s_quantity, quantity_resampled
