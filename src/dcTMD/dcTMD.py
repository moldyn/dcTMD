# -*- coding: utf-8 -*-
"""
Classes calculating the dcTMD quantities.


MIT License
Copyright (c) 2022, Victor Tänzel, Miriam Jäger
All rights reserved.
"""

__all__ = []

import numpy as np
from beartype import beartype
from beartype.typing import Optional, Union
from sklearn.base import BaseEstimator, TransformerMixin

from dcTMD._typing import (
    Int,
    Float,
    ArrayLikeStr,
    Float1DArray,
    Float2DArray,
    StrStd,
    NumInRange0to1,
)

@beartype
def _bootstrap_mode_reducer(obj, *args):
    """
    Perform the last bootstrap step.
    
    Reduces sequences of a resampled statistics by calculating standard
    deviations or confidence intervals, depending on the 'mode' attribute of
    the given object.    
    """
    if obj.mode == 'std':
        def reducer(x):
            return np.std(x, axis=0)
    else:
        import scipy.stats as st

        def reducer(x):
            return st.t.interval(
                obj.mode,
                len(x)-1,
                loc=np.mean(x),
                scale=st.sem(x),
            )
    return tuple(reducer(arg) for arg in args)



def Boris(TransformerMixin, BaseEstimator):
    @beartype
    def __init__(
        self,
        temperature: Float,
        verbose: bool = False,
    ) -> None:
        """Initialize class."""
        self.temperature = temperature

    @beartype
    def fit(
        self,
        work_set: Float2DArray,
    ):
        self.work_set = work_set
        self.estimate_free_energy()
        return self

    @beartype
    def estimate_free_energy(self, work_set=None):
        if work_set is None:
            work_set = self.work_set
        from scipy.constants import R  # noqa: WPS347
        RT = R * self.temperature / 1e3

        self.W_mean_ = np.mean(work_set, axis=0)
        W_var = np.var(work_set, axis=0)
        self.W_diss_ = 1 / (2 * RT) * W_var
        self.dG_ = self.W_mean_ - self.W_diss_
        return self.W_mean_, self.W_diss_, self.dG_

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
        self.mode = mode
        self._bootstrap_free_energy(mode)
        return self.s_W_mean_, self.s_W_diss_, self.s_dG_

    @beartype
    def _bootstrap_free_energy(
        self,
    ):
        """Perform bootstrapping for the free energy."""
        import tqdm
        N_traj, length_data = np.shape(self.work_set)

        W_mean_resampled = np.empty((self.N_resamples, length_data))
        W_diss_resampled = np.empty((self.N_resamples, length_data))
        dG_resampled = np.empty((self.N_resamples, length_data))

        for ind in tqdm.tqdm(
            range(self.N_resamples),
            desc='Bootstrapping progress',
        ):
            random_indices = np.random.randint(0, N_traj, N_traj)
            work_set_resampled = self.work_set[random_indices]
            W_mean, W_diss, dG = self.estimate_free_energy(
                work_set=work_set_resampled,
            )
            W_mean_resampled[ind] = W_mean
            W_diss_resampled[ind] = W_diss
            dG_resampled[ind] = dG
        s_W_mean, s_W_diss, s_dG,  = _bootstrap_mode_reducer(
            W_mean_resampled,
            W_diss_resampled,
            dG_resampled,
        )
        self.s_W_mean_ = s_W_mean
        self.s_W_diss_ = s_W_diss
        self.s_dG_ = s_dG
        self.W_mean_resampled_ = W_mean_resampled
        self.W_diss_resampled_ = W_diss_resampled
        self.dG_resampled_ = dG_resampled
