# -*- coding: utf-8 -*-
"""
Classes calculating the dcTMD quantities.

MIT License
Copyright (c) 2022, Victor Tänzel, Miriam Jäger
All rights reserved.
"""

__all__ = ['WorkEstimator', 'ForceEstimator']

import numpy as np
from beartype import beartype
from beartype.typing import Union, Tuple, Optional
from sklearn.base import BaseEstimator, TransformerMixin

from dcTMD.utils import bootstrapping
from dcTMD._typing import (
    Int,
    Float,
    Str,
    StrStd,
    NumInRange0to1,
    Float1DArray,
    Float2DArray,
)


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
    position_ :
        Positions time trace, product of time trace and velocity, in nm.
    W_mean_ :
        Mean work, in kJ/mol.
    W_diss_ :
        Dissipative work, in kJ/mol.
    dG_ :
        Free energy estimate, in kJ/mol.
    friction_:
        Friction factor in kJ/mol/(nm^2/ps).
    mode_ :
        Parameter of estimate_free_energy_errors(). Decides how the
        bootstrapping errors are calculated.
    s_W_mean_ :
        Bootstrapping error of the mean work. Calculated via
        estimate_free_energy_errors().
    s_W_diss_ :
        Bootstrapping error of the dissipative work. Calculated via
        estimate_free_energy_errors().
    s_dG_ :
        Bootstrapping error of the free energy estimate. Calculated via
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

    Examples
    --------
    >>> from dcTMD.dcTMD import WorkEstimator
    >>> from dcTMD.storing import load
    >>> work = load('my_work_set')
    >>> # Instantiate a WorkEstimator instance and fit it with the WorkSet
    >>> # instance work
    >>> work_estimator = WorkEstimator(temperature=290.15)
    >>> work_estimator.fit(work)
    >>> work_estimator.dG_
    array([..., ])
    """

    @beartype
    def __init__(
        self,
        temperature: (Float, Int),
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
            Instance of WorkSet containing constraint force work time traces,
            for which the free energy and friction are estimated.

        Returns
        -------
        self :
            Fitted estimator.
        """
        self.work_set = work_set
        self.position_ = work_set.position_
        self.names_ = work_set.names_
        self.estimate_free_energy()
        self.estimate_friction()
        return self

    @beartype
    def transform(
        self,
        X,
        y=None,
    ) -> Tuple[Float1DArray, Float1DArray]:
        """Return free energy and friction estimates."""
        return self.dG_, self.friction_

    @beartype
    def estimate_free_energy(
        self,
        work_set=None,
    ) -> Tuple[Float1DArray, Float1DArray, Float1DArray]:
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
        dG_ : 1D np.array
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
        n_resamples: Int,
        mode: Union[StrStd, NumInRange0to1],
        seed: Optional[Int] = None,
    ) -> Tuple[Union[Float1DArray, Float2DArray],
               Union[Float1DArray, Float2DArray],
               Union[Float1DArray, Float2DArray],
               ]:
        """
        Estimate bootstrapping errors for the free energy estimate.

        Bootstrapping errors are calculated for the free energy estimate and
        the related quantities mean and dissipative work. Return matches the
        one of [dcTMD.dcTMD.WorkEstimator.estimate_free_energy][].

        Parameters
        ----------
        n_resamples :
            Number of drawn resamples for bootstrapping error analysis.
        mode :
            Chooses between reducing the resampled statistic via (1) 'std' the
            element-wise calculation of standard deviations or (2) confidence
            intervals if `mode` is a float in [0, 1).
        seed :
            Seed for the random number generator.

        Returns
        -------
        s_W_mean_ :
            Error estimate of the mean work.
        s_W_diss_ :
            Error estimate of the mean work.
        s_dG_ :
            Error estimate of free energy.

        Examples
        --------
        >>> from dcTMD.dcTMD import WorkEstimator
        >>> work_estimator.estimate_free_energy_errors(1000, mode='std')
        Bootstrapping progress: 100%|██████████| 1000/1000 [00:00<00:00, 12797.15it/s]
        >>> work_estimator.s_dG_
        array([..., ])
        """
        self.free_energy_error_ = {
            'mode': mode,
            'n_resamples': n_resamples,
            'seed': seed,
        }
        self._bootstrap_free_energy()
        return self.s_W_mean_, self.s_W_diss_, self.s_dG_

    @beartype
    def _bootstrap_free_energy(
        self,
    ) -> None:
        """Use utils/bootstrapper.py for error estimation."""
        # Prepare and run bootstrapper
        def func(my_work_set):
            return self.estimate_free_energy(my_work_set)
        s_quantity, quantity_resampled = bootstrapping.bootstrapping(
            self,
            func=func,
            descriptor=self.free_energy_error_,
        )
        # Save error estimates and bootstrapped quantities
        self.W_mean_resampled_ = quantity_resampled[:, 0]
        self.W_diss_resampled_ = quantity_resampled[:, 1]
        self.dG_resampled_ = quantity_resampled[:, 2]
        
        if self.free_energy_error_['mode'] == 'std':
            self.s_W_mean_ = s_quantity[0, 0]
            self.s_W_diss_ = s_quantity[0, 1]
            self.s_dG_ = s_quantity[0, 2]
        else:
            self.s_W_mean_ = s_quantity[0, :, 0]
            self.s_W_diss_ = s_quantity[0, :, 1]
            self.s_dG_ = s_quantity[0, :, 2]
            

    @beartype
    def estimate_friction(
        self,
        W_diss: Optional[Float1DArray] = None,
    ) -> Float1DArray:
        """
        Estimate bootstrapping errors for the friction.

        Besides calculating dcTMD quantitites to the class, this function
        is also called from the bootstrapping routine. In the latter case,
        which comes with a passed `W_diss` parameter, attributes should not
        be overwritten.

        Parameters
        ----------
        W_diss :
            Dissipative work, in kJ/mol. Is passed if this function is used
            during bootstrapping.

        Returns
        -------
        friction :
            Friction factor in kJ/mol/(nm^2/ps).
        """
        if W_diss is None:
            is_bootstrapping = False
            try:
                W_diss = self.W_diss_
            except AttributeError:
                self.estimate_free_energy()
                W_diss = self.W_diss_
        else:
            is_bootstrapping = True

        delta_x = self.work_set.position_[1] - self.work_set.position_[0]
        if self.verbose:
            print(f'calculating friction, delta_x: {delta_x}nm')
        friction = np.diff(
            W_diss, prepend=W_diss[0],
        ) / (delta_x * self.work_set.velocity)
        if not is_bootstrapping:
            self.friction_ = friction
        return friction

    @beartype
    def estimate_friction_errors(
        self,
        n_resamples: Int,
        mode: Union[StrStd, NumInRange0to1],
        seed: Optional[Int] = None,
    ) -> Float1DArray:
        """
        Estimate bootstrapping errors for the free energy estimate.

        Bootstrapping errors are calculated for the free energy estimate and
        the related quantities mean and dissipative work. Return matches the
        one of estimate_free_energy().

        Parameters
        ----------
        n_resamples :
            Number of drawn resamples for bootstrapping error analysis.
        mode :
            Chooses between reducing the resampled statistic via
                1.  'std' the element-wise calculation of standard deviations,
                2.  confidence intervals if `mode` is a float in [0, 1).
        seed :
            Seed for the random number generator.

        Returns
        -------
        s_friction_ :
            Bootstrap error of the friction factor in kJ/mol/(nm^2/ps).

        Examples
        --------
        >>> from dcTMD.dcTMD import WorkEstimator
        >>> work_estimator.estimate_friction_errors(1000, mode='std')
        Bootstrapping progress: 100%|██████████| 1000/1000 [00:00<00:00, 10245.63it/s]
        >>> work_estimator.s_friction_
        array([..., ])
        """
        self.friction_error_ = {
            'n_resamples': n_resamples,
            'mode': mode,
            'seed': seed,
        }
        self._bootstrap_friction()
        return self.s_friction_

    @beartype
    def _bootstrap_friction(
        self,
    ) -> None:
        """Use utils/bootstrapper.py for error estimation."""
        # Prepare and run bootstrapper
        def func(my_work_set):
            return self.estimate_friction(
                self.estimate_free_energy(my_work_set)[1],
            )
        s_quantity, quantity_resampled = bootstrapping.bootstrapping(
            self,
            func=func,
            descriptor=self.friction_error_,
        )
        # Save error estimates and bootstrapped quantities
        if self.free_energy_error_['mode'] == 'std':
            self.s_friction_ = s_quantity[0, 0]
        else:
            self.s_friction_ = s_quantity[0, :, 0]
        self.friction_resampled_ = quantity_resampled[:, 0]

    @beartype
    def smooth_friction(
        self,
        sigma: Float,
        mode: Str = "reflect",
    ) -> Float1DArray:
        """Smooth friction with gaussain kernel.

        Parameters
        ----------
        sigma:
            standard deviation of gaussian kernel in nm
        mode:
            options: ‘reflect’, ‘constant’, ‘nearest’, ‘mirror’, ‘wrap’
            The mode parameter determines how the input array is
            extended beyond its boundaries. Default is ‘reflect’.
            Behavior for each option see scipy.ndimage.gaussian_filter1d

        Returns
        -------
        friction_smooth_ : 1d np.array
            Smoothed friction.
        """
        from dcTMD.utils.smoothing import gaussfilter_friction
        self.friction_smooth_ = gaussfilter_friction(self.friction_,
                                                     self.position_,
                                                     sigma=0.1,
                                                     )
        return self.friction_smooth_


class ForceEstimator(TransformerMixin, BaseEstimator):
    """
    Class for performing dcTMD analysis on a force set.

    Parameters
    ----------
    temperature :
        Temperature at which the simulations were carried out, in K.
    verbose :
        Enables verbose mode.

    Attributes
    ----------
    position_ :
        Positions time trace, product of time trace and velocity, in nm.
    W_mean_ :
        Mean work, in kJ/mol.
    W_diss_ :
        Dissipative work, in kJ/mol.
    dG_ :
        Free energy estimate, in kJ/mol.
    friction_ :
        Friction factor in kJ/mol/(nm^2/ps).

    Examples
    --------
    >>> from dcTMD.dcTMD import ForceEstimator
    >>> from dcTMD.storing import load
    >>> force = load('my_force_set')
    >>> # Instantiate a ForceEstimator instance and fit it with the
    >>> # ForceSet instance
    >>> force_estimator = ForceEstimator(temperature=290.15)
    >>> force_estimator.fit(force)
    >>> force_estimator.dG_
    array([..., ])
    """

    @beartype
    def __init__(
        self,
        temperature: (Float, Int),
        verbose: bool = False,
    ) -> None:
        """Initialize class."""
        self.temperature = temperature
        self.verbose = verbose

    @beartype
    def fit(
        self,
        force_set,
    ):
        """
        Estimate free energy and friction.

        Parameters
        ----------
        force_set :
            Instance of ForceSet containing constraint forces, for which the
            free energy and friction are estimated.

        Returns
        -------
        self :
            Fitted estimator.
        """
        self.force_set = force_set
        self.names_ = force_set.names_
        self.estimate_free_energy_friction()
        return self

    @beartype
    def transform(
        self,
        X,
        y=None,
    ) -> Tuple[Float1DArray, Float1DArray]:
        """Return free energy and friction estimates."""
        return self.dG_, self.friction_

    @beartype
    def estimate_free_energy_friction(
        self,
    ) -> Tuple[Float1DArray, Float1DArray, Float1DArray, Float1DArray]:
        """
        Estimate free energy and friction.
        From force auto corrleation

        Returns
        -------
        W_mean : 1D np.array
            Mean work, in kJ/mol.
        W_diss : 1D np.array
            Dissipative work, in kJ/mol.
        dG_ : 1D np.array
            Free energy estimate, in kJ/mol.
        friction_ : 1D np.array
            Friction factor in kJ/mol/(nm^2/ps).
        """
        from scipy.constants import R  # noqa: WPS347
        from scipy.integrate import cumulative_trapezoid
        RT = R * self.temperature / 1e3

        # average and variance over all trajectories in each time step
        # shape: (length_data)
        force_mean = np.mean(self.force_set.force_, axis=0)

        # calculate $\delta f_c(t) = f_c(t) - \left< f_c (t) \right>_N$
        self.delta_force_array = self.force_set.force_ - force_mean

        # integrate over time
        int_delta_force = cumulative_trapezoid(
            self.delta_force_array,
            self.force_set.time_,
            axis=-1,
            initial=0,
        )
        intcorr = np.multiply(
            self.delta_force_array,
            int_delta_force,
        )
        friction_ = np.mean(intcorr, axis=0) / RT

        print('Calculating dissipative work...')
        W_mean_ = cumulative_trapezoid(
            force_mean,
            self.force_set.position_,
            initial=0,
        )
        W_diss = cumulative_trapezoid(
            friction_,
            self.force_set.position_,
            initial=0,
        ) * self.force_set.velocity

        # Reduce resolution
        self.position_ = self.force_set.position_[::self.force_set.resolution]  # noqa: E501
        self.W_mean_ = W_mean_[::self.force_set.resolution]
        self.W_diss_ = W_diss[::self.force_set.resolution]
        self.dG_ = self.W_mean_ - self.W_diss_
        self.friction_ = friction_[::self.force_set.resolution]

        return self.W_mean_, self.W_diss_, self.dG_, self.friction_

    @beartype
    def smooth_friction(
        self,
        sigma: Float,
        mode: Str = "reflect",
    ) -> Float1DArray:
        """Smooth friction with gaussain kernel.

        Parameters
        ----------
        sigma:
            standard deviation of gaussian kernel in nm
        mode:
            options: ‘reflect’, ‘constant’, ‘nearest’, ‘mirror’, ‘wrap’
            The mode parameter determines how the input array is
            extended beyond its boundaries. Default is ‘reflect’.
            Behavior for each option see scipy.ndimage.gaussian_filter1d

        Returns
        -------
        friction_smooth_ : 1d np.array
            Smoothed friction.
        """
        from dcTMD.utils.smoothing import gaussfilter_friction
        self.friction_smooth_ = gaussfilter_friction(self.friction_,
                                                     self.position_,
                                                     sigma=0.1,
                                                     )
        return self.friction_smooth_

    def memory_kernel(
        self,
        x_indices: Float1DArray,
    ) -> Float1DArray:
        """
        Calculate memory kernel at positions X "forward" in time
        from fluctuation-dissipation.
        see e.g. R. Zwanzig, “Nonequilibrium statistical mechanics”,
        Oxford University Press (2001).

        Parameters
        ----------
        x_indices : np.ndarray
            Indices at which memory kernel is calculated.

        Returns
        -------
        corr_set : np.ndarray
            shape: (len(X), length_data)
        NaN are set to zero
        """
        _, length_data = self.delta_force_.shape
        corr_set = np.zeros((len(x_indices), length_data))

        for ind, tt in enumerate(range(length_data)):
            entries = self.delta_force_[:, tt:-2] * \
                self.delta_force_[:, tt + 1:-1]
            corr_set[ind, tt:-2] = np.mean(
                entries,
                axis=0,
            )
        return corr_set
