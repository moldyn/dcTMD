# -*- coding: utf-8 -*-
"""
Classes calculating the dcTMD quantities.

MIT License
Copyright (c) 2022, Victor Tänzel, Miriam Jäger
All rights reserved.
"""

__all__ = ['Boris']

import numpy as np
from beartype import beartype
from beartype.typing import Union, Tuple
from sklearn.base import BaseEstimator, TransformerMixin

from dcTMD._typing import (
    Int,
    Float,
    Float1DArray,
    Float2DArray,
    StrStd,
    NumInRange0to1,
)


@beartype
def _bootstrap_mode_reducer(
    obj,
    *args: Float2DArray,
) -> Tuple[
    Union[Float1DArray, Tuple[Float1DArray, Float1DArray]],
    ...,
]:
    """
    Perform the last bootstrap step.

    Reduces sequences of a resampled statistics by calculating standard
    deviations or confidence intervals, depending on the 'mode' attribute of
    the given object obj. If obj.mode is a real number in [0, 1), confidence
    intervals will be computed. If it is the string 'std', then the standard
    deviation is calculated.

    Parameters
    ----------
    obj :
        Instance of Boris which provides its 'mode' attribute.
    *args :
        Resampled statistics, to be reduced.

    Returns
    -------
    tuple :
        Tuple of reduced statistics.
    """
    if obj.mode_ == 'std':
        def reducer(x):
            return np.std(x, axis=0)
    else:
        import scipy.stats as st

        def reducer(x):
            return st.t.interval(
                obj.mode_,
                len(x)-1,
                loc=np.mean(x),
                scale=st.sem(x),
            )
    return tuple(reducer(arg) for arg in args)


class WorkEstimator(TransformerMixin, BaseEstimator):
    """
    Class for performing dcTMD analysis on a work set.

    Parameters
    ----------
    temperature :
        Temperature at which the simulations were carried out, in K.
    verbose :
        Enables verbose mode.

    Attributes
    ----------
    W_mean_ :
        Mean work, in kJ/mol.
    W_diss_ :
        Dissipative work, in kJ/mol.
    dG_ :
        Free energy estimate, in kJ/mol.
    mode_ :
        Parameter of estimate_free_energy_errors(). Decides how the
        bootstrapping errors are calculated.        
    s_W_mean_ :
        Boostrapping error of the mean work. Calculated via 
        estimate_free_energy_errors().
    s_W_diss_ :
        Boostrapping error of the dissipative work. Calculated via
        estimate_free_energy_errors().
    s_dG_ :
        Boostrapping error of the free energy estimate. Calculated via
        estimate_free_energy_errors().
    W_mean_resampled_ :
        Resampled mean work, needed to inspect its distribution. Calculated
        via estimate_free_energy_errors().
    W_diss_resampled_ :
        Resampled dissipative work, needed to inspect its distribution.
        Calculated via estimate_free_energy_errors().
    dG_resampled_ :
        Resampled free energy estimate, needed to inspect its distribution.
        Calculated via estimate_free_energy_errors().

    Examples :


    """
    @beartype
    def __init__(
        self,
        temperature: Float,
        verbose: bool = False,
    ) -> None:
        """Initialize class."""
        self.temperature = temperature
        self.verbose = verbose

    @beartype
    def fit(
        self,
        work_set,
    ):
        """
        Estimate free energy and friction.

        Parameters
        ----------
        work_set :
            Instance of WorkSet containing constraint forces, for which the
            free energy and friction are estimated.

        Returns
        -------
        self :
            Fitted estimator.
        """
        self.work_set = work_set
        self.estimate_free_energy()
        self.estimate_friction()
        return self

    @beartype
    def estimate_free_energy(self, work_set=None):
        """
        Estimate free energy.

        Parameters
        ----------
        work_set : optional
            Instance of a WorkSet containing constraint forces, for which the 
            free energy and friction are estimated.

        Returns
        -------
        W_mean : 1D np.array
            Mean work, in kJ/mol.
        W_diss : 1D np.array
            Dissipative work, in kJ/mol.
        dG : 1D np.array
            Free energy estimate, in kJ/mol.
        """
        # Besides calculating dcTMD quantitites for the class, this function
        # is also called from the bootstrapping routine. In the latter case,
        # which comes with a passed work_set parameter, attributes should not
        # be overwritten.
        if work_set is None:
            work_set = self.work_set.work_
            is_bootstrapping = False
        else:
            is_bootstrapping = True

        from scipy.constants import R  # noqa: WPS347
        RT = R * self.temperature / 1e3

        W_mean = np.mean(work_set, axis=0)
        W_var = np.var(work_set, axis=0)
        W_diss = 1 / (2 * RT) * W_var
        dG = W_mean - W_diss

        if not is_bootstrapping:
            self.W_mean_ = W_mean
            self.W_diss_ = W_diss
            self.dG_ = dG
        return W_mean, W_diss, dG

    @beartype
    def estimate_free_energy_errors(
        self,
        N_resamples: Int,
        mode: Union[StrStd, NumInRange0to1],
    ):
        """
        Estimate bootstrapping errors for the free energy estimate.

        Boostrapping errors are calculated for the free energy estimate and the
        related quantities mean and dissipative work. Return matches the one of
        estimate_free_energy().

        Parameters
        ----------
        N_resamples :
            Number of drawn resamples for bootstrapping error analysis.
        mode :
            Chooses between reducing the resampled statistic via (1) 'std' the
            element-wise calculation of standard deviations or (2) confidence
            intervals if 'mode' is a float in [0, 1).

        Returns
        -------
        s_W_mean_ :
            Error estimate of the mean work.
        s_W_diss_ :
            Error estimate of the mean work.
        s_dG_ :
            Error estimate of free energy.
        """
        self.N_resamples = N_resamples
        self.mode_ = mode
        self._bootstrap_free_energy()
        return self.s_W_mean_, self.s_W_diss_, self.s_dG_

    @beartype
    def _bootstrap_free_energy(
        self,
    ):
        """Perform bootstrapping for the free energy."""
        import tqdm
        N_traj, length_data = np.shape(self.work_set.work_)

        W_mean_resampled = np.empty((self.N_resamples, length_data))
        W_diss_resampled = np.empty((self.N_resamples, length_data))
        dG_resampled = np.empty((self.N_resamples, length_data))

        for ind in tqdm.tqdm(
            range(self.N_resamples),
            desc='Bootstrapping progress',
        ):
            # Draw random work time traces
            random_indices = np.random.randint(0, N_traj, N_traj)
            work_set_resampled = self.work_set.work_[random_indices]
            # Calculate and save the relevant statistics
            W_mean, W_diss, dG = self.estimate_free_energy(
                work_set=work_set_resampled,
            )
            W_mean_resampled[ind] = W_mean
            W_diss_resampled[ind] = W_diss
            dG_resampled[ind] = dG
        # There are now boostrapped quantities in the '_resampled' variables.
        # We are interested in the element-wise distributions and thus
        # calculate (1) the standard distribution of the resampled quantity
        # at all points or (2) confidence intervals.
        s_W_mean, s_W_diss, s_dG,  = _bootstrap_mode_reducer(
            self,
            W_mean_resampled,
            W_diss_resampled,
            dG_resampled,
        )
        self.s_W_mean_ = s_W_mean
        self.s_W_diss_ = s_W_diss
        self.s_dG_ = s_dG
        # The distributions of the '_resampled' variables must be inspected
        # and are thus saved as attributes.
        self.W_mean_resampled_ = W_mean_resampled
        self.W_diss_resampled_ = W_diss_resampled
        self.dG_resampled_ = dG_resampled

    @beartype
    def estimate_friction(self, W_diss=None):
        """Estimate bootstrapping errors for the friction."""
        # Besides calculating dcTMD quantitites to the class, this function
        # is also called from the bootstrapping routine. In the latter case,
        # which comes with a passed work_set parameter, attributes should not
        # be overwritten.
        if W_diss is None:
            is_bootstrapping = False
            try:
                W_diss = self.W_diss_
            except:
                self.estimate_free_energy()
                W_diss = self.W_diss_
        else:
            is_bootstrapping = True

        delta_x = self.work_set.position_[1] - self.work_set.position_[0]
        if self.verbose:
            print(f'calculating friction, delta_x: {delta_x}nm')
        friction = np.diff(
            W_diss, prepend=W_diss[0]
        ) / (delta_x * self.work_set.velocity)
        if not is_bootstrapping:
            self.friction = friction
        return friction

    @beartype
    def estimate_friction_errors(
        self,
        N_resamples: Int,
        mode: Union[StrStd, NumInRange0to1],
    ):
        """Estimate bootstrapping errors for """
        self.N_resamples = N_resamples
        self.mode_ = mode
        self._bootstrap_friction()
        return self.s_friction

    @beartype
    def _bootstrap_friction(
        self,
    ):
        """Perform bootstrapping for the free energy."""
        import tqdm
        N_traj, length_data = np.shape(self.work_set.work_)
        friction_resampled = np.empty((self.N_resamples, length_data))

        for ind in tqdm.tqdm(
            range(self.N_resamples),
            desc='Bootstrapping progress',
        ):
            # Draw random work time traces
            random_indices = np.random.randint(0, N_traj, N_traj)
            work_set_resampled = self.work_set.work_[random_indices]
            # Calculate and save the relevant statistics
            _, W_diss, _ = self.estimate_free_energy(
                work_set=work_set_resampled,
            )
            friction = self.estimate_friction(W_diss=W_diss)
            friction_resampled[ind] = friction
        # There are now boostrapped quantities in the '_resampled' variables.
        # We are interested in the element-wise distributions and thus
        # calculate (1) the standard distribution of the resampled quantity
        # at all points or (2) confidence intervals.
        s_friction,  = _bootstrap_mode_reducer(
            self,
            friction,
        )
        self.s_friction_ = s_friction
        # The distributions of the '_resampled' variables must be inspected
        # and are thus saved as attributes.
        self.friction_resampled_ = friction_resampled
