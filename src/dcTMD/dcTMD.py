# -*- coding: utf-8 -*-
# MIT License
# Copyright (c) 2022, Victor Tänzel, Miriam Jäger
# All rights reserved.
"""Classes `WorkEstimator`, `ForceEstimator` calculating the dcTMD quantities.

This submodule contains two classes, WorkEstimator and ForceEstimator, which
are used for the dcTMD analysis of constraint force time traces. Both class can
be used to calculate the mean work, dissipative work, free energy and friction
estimate of a set of constraint force time traces.
"""

__all__ = ['WorkEstimator', 'ForceEstimator']

import numpy as np
from abc import ABC
from beartype import beartype
from beartype.typing import Union, Tuple, Optional
from sklearn.base import BaseEstimator, TransformerMixin

from dcTMD import utils
from dcTMD._typing import (
    Int,
    Float,
    Str,
    StrStd,
    NumInRange0to1,
    Float1DArray,
    Float2DArray,
    Index1DArray,
)


class _SmoothBasisEstimator(ABC):
    """Class with the smoothing method for `WorkEstimator`, `ForceEstimator`.
    """
    @beartype
    def smooth_friction(
        self,
        sigma: Float,
        mode: Str = 'reflect',
    ) -> Float1DArray:
        """Smooth friction with gaussian kernel.

        Parameters
        ----------
        sigma:
            standard deviation of gaussian kernel in nm
        mode:
            options: `reflect`, `constant`, `nearest`,`mirror`, `wrap`
            The mode parameter determines how the input array is
            extended beyond its boundaries. Default is `reflect`.
            Behavior for each option see scipy.ndimage.gaussian_filter1d.

        Returns
        -------
        friction_smooth_ : 1d np.array
            Smoothed friction.
        """
        self.friction_smooth_ = utils.gaussfilter_friction(
            self.friction_,
            self.position_,
            sigma=sigma,
            mode=mode,
        )
        return self.friction_smooth_

    @beartype
    def _reset(self) -> None:
        """Reset friction_smooth_ attribute."""
        if hasattr(self, 'friction_smooth_'):  # noqa: WPS421
            del self.friction_smooth_  # noqa: WPS420


class WorkEstimator(
    _SmoothBasisEstimator, TransformerMixin, BaseEstimator,
):  # noqa: WPS230, WPS214
    """Class for performing dcTMD analysis on a work set.

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
        Parameter of [WorkEstimator.estimate_free_energy_errors][dcTMD.dcTMD.\
WorkEstimator.estimate_free_energy_errors]. Decides how the bootstrapping
        errors are calculated.
    s_W_mean_ :
        Bootstrapping error of the mean work. Calculated via [WorkEstimator.\
estimate_free_energy_errors][dcTMD.dcTMD.WorkEstimator.\
estimate_free_energy_errors].
    s_W_diss_ :
        Bootstrapping error of the dissipative work. Calculated via
        [WorkEstimator.estimate_free_energy_errors][dcTMD.dcTMD.WorkEstimator.\
estimate_free_energy_errors].
    s_dG_ :
        Bootstrapping error of the free energy estimate. Calculated via
        [WorkEstimator.estimate_free_energy_errors][dcTMD.dcTMD.WorkEstimator.\
estimate_free_energy_errors].
    W_mean_resampled_ :
        Resampled mean work, needed to inspect its distribution. Calculated
        via [WorkEstimator.estimate_free_energy_errors][dcTMD.dcTMD.\
WorkEstimator.estimate_free_energy_errors].
    W_diss_resampled_ :
        Resampled dissipative work, needed to inspect its distribution.
        Calculated via [WorkEstimator.estimate_free_energy_errors][dcTMD.\
dcTMD.WorkEstimator.estimate_free_energy_errors].
    dG_resampled_ :
        Resampled free energy estimate, needed to inspect its distribution.
        Calculated via [WorkEstimator.estimate_free_energy_errors][dcTMD.\
dcTMD.WorkEstimator.estimate_free_energy_errors].

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
        temperature: Union[Float, Int],
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
        self._reset()
        self.work_set = work_set
        self.position_ = work_set.position_
        self.names_ = work_set.names_
        self.estimate_free_energy()
        self.estimate_friction()
        return self

    @beartype
    def transform(
        self, X, y=None,
    ) -> Tuple[Float1DArray, Float1DArray]:  # noqa: WPS111
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
    def estimate_free_energy_errors(  # noqa: WPS320
        self,
        n_resamples: Int,
        mode: Union[StrStd, NumInRange0to1],
        seed: Optional[Int] = None,
    ) -> Tuple[Union[Float1DArray, Float2DArray],
               Union[Float1DArray, Float2DArray],  # noqa: WPS318
               Union[Float1DArray, Float2DArray],
               ]:
        """
        Estimate bootstrapping errors for the free energy estimate.

        Bootstrapping errors are calculated for the free energy estimate and
        the related quantities mean and dissipative work. Return matches the
        one of [WorkEstimator.estimate_free_energy][dcTMD.dcTMD.WorkEstimator.\
estimate_free_energy].

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
        Bootstrapping progress: 100%|██████████| 1000/1000 [00:00<00:00, 12797.15it/s]  # noqa
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
        s_quantity, quantity_resampled = utils.bootstrapping(
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

        Besides calculating dcTMD quantities to the class, this function
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
        one of [WorkEstimator.estimate_free_energy][dcTMD.dcTMD.WorkEstimator.\
estimate_free_energy].

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
        Bootstrapping progress: 100%|██████████| 1000/1000 [00:00<00:00, 10245.63it/s]  # noqa
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
        s_quantity, quantity_resampled = utils.bootstrapping(
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


class ForceEstimator(
    _SmoothBasisEstimator, TransformerMixin, BaseEstimator,
):  # noqa: WPS230
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
        temperature: Union[Float, Int],
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
        self._reset()
        self.force_set = force_set
        self.names_ = force_set.names_
        self.estimate_free_energy_friction()
        return self

    @beartype
    def transform(
        self, X, y=None,
    ) -> Tuple[Float1DArray, Float1DArray]:  # noqa: WPS111
        """Return free energy and friction estimates."""
        return self.dG_, self.friction_

    @beartype
    def estimate_free_energy_friction(
        self,
    ) -> Tuple[Float1DArray, Float1DArray, Float1DArray, Float1DArray]:
        """
        Estimate free energy and friction from force auto correlation.

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

        # average over all trajectories in each time step
        force_mean = np.mean(self.force_set.force_, axis=0)

        # calculate $\delta f(t) = f(t) - \left< f(t) \right>_N$
        self.delta_force_array = self.force_set.force_ - force_mean

        # integrate f(t) over time
        int_delta_force = cumulative_trapezoid(
            self.delta_force_array,
            self.force_set.time_,
            axis=-1,
            initial=0,
        )
        # multiply $\delta f(t) \int_0^t \delta f(t') dt'$ for each N
        intcorr = np.multiply(
            self.delta_force_array,
            int_delta_force,
        )
        friction = np.mean(intcorr, axis=0) / RT

        W_mean = cumulative_trapezoid(
            force_mean,
            self.force_set.position_,
            initial=0,
        )
        W_diss = cumulative_trapezoid(
            friction,
            self.force_set.position_,
            initial=0,
        ) * self.force_set.velocity

        # Reduce resolution
        self.position_ = self.force_set.position_[::self.force_set.resolution]
        self.W_mean_ = W_mean[::self.force_set.resolution]
        self.W_diss_ = W_diss[::self.force_set.resolution]
        self.dG_ = self.W_mean_ - self.W_diss_
        self.friction_ = friction[::self.force_set.resolution]

        return self.W_mean_, self.W_diss_, self.dG_, self.friction_

    @beartype
    @staticmethod
    def kernel_at_ndx(
        delta_force_array: Float2DArray,
        ndx: Int
    ) -> Float1DArray:
        """
        Calculate the kernel at a specific index.

        Args:
            ndx (Int): Index at which the kernel is calculated.

        Returns:
            Float1DArray: The mean force correlation
            < df(t(x)) df(t) >_N at the given index.
        """
        delta_force_point = delta_force_array[:, ndx]
        force_correlation_at_ndx = (
            delta_force_array.T * delta_force_point
        ).T
        return np.mean(force_correlation_at_ndx, axis=0)

    @beartype
    def memory_kernel(
        self,
        index: Union[Int, Index1DArray, None] = None,
        ndx_striding: Union[Int, None] = None,
    ) -> Float2DArray:
        """
        Calculate memory kernel at index.

        Either give a index (as Int or an array of indices) or
        ndx_striding as an argument. The latter creates the index array
        from the force data with striding.

        Parameters
        ----------
        x_indices :
            Indices at which the memory kernel is calculated.
            If None, indices will
            be generated based on `ndx_striding`. Default is None.
        ndx_striding:
            Resolution for creating index array.
            If provided, indices will be
            generated at intervals of `ndx_resolution`. Default is None.

        Returns
        -------
        numpy.ndarray
            - A 2D NumPy array containing the memory kernel values.

        Examples
        --------
        >>> from dcTMD.dcTMD import ForceEstimator
        >>> from dcTMD.storing import load
        >>> force = load('my_force_set')
        >>> # Instantiate a ForceEstimator instance and fit it with the
        >>> # ForceSet instance
        >>> force_estimator = ForceEstimator(temperature=290.15)

        # Example usage with specific indices:

        >>> kernel = force_estimator.memory_kernel(index=[10, 20, 30])
        >>> print(kernel)

        # Example usage with resolution:

        >>> force_estimator.memory_kernel(ndx_resolution=1000)
        >>> print(force_estimator.memory_kernel_index_)
        >>> print(force_estimator.memory_kernel_)
        """
        # Argument validation (mutual exclusion / requirement)
        if index is not None and ndx_striding is not None:
            raise ValueError(
                'Only index or ndx_resolution can be given.'
            )
        if index is None and ndx_striding is None:
            raise ValueError(
                'Either index or ndx_resolution must be given.'
            )
        # read in index or create index array
        if index is not None:
            if isinstance(index, (int, np.integer)):
                print('create index with single int')
                index = np.array([index])
            if np.any(index >= len(self.force_set.time_)):
                raise ValueError(
                    'Index values must be less than length of data.'
                )
        elif ndx_striding is not None:
            print('create index with ndx_resolution')
            index = np.arange(
                ndx_striding,
                len(self.force_set.time_) - 1,
                ndx_striding,
                dtype=int
            )
        # calculate memory kernel at given indices
        correlation_set = np.zeros((len(index), len(self.force_set.time_)))
        for i, ndx in enumerate(index):
            correlation_set[i] = self.kernel_at_ndx(
                self.delta_force_array,
                ndx
            )
        self.memory_kernel_ = correlation_set
        self.memory_kernel_index_ = index
        return correlation_set
