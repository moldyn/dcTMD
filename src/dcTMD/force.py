# -*- coding: utf-8 -*-
"""Calculate dcTMD quantities via the autocovariance function of the force."""
__all__ = []


import numpy as np
from beartype import beartype
from beartype.typing import Tuple
from dcTMD._typing import (
    Int,
    Float,
    ArrayLikeStr,
    Float1DArray,
    Float2DArray,
)


@beartype
def _fill_force_array(
    time: Float1DArray,
    file_names: ArrayLikeStr,
    verbose: bool = False,
) -> Tuple[Float2DArray, ArrayLikeStr]:
    """
    Fill the force array by reading in force files.

    Parameters
    ----------
    time : np.ndarray
        Time trace in ps.
    file_names: array_like
        Contains force file names.

    Returns
    -------
    force_array: np.ndarray
        Array containing constraint forces with shape=(len(file_names),
        len(t)).
    force_array_names : list
        List with force file names.
    """
    # allocate memory
    force_array = np.zeros((len(file_names), len(time)))
    force_array_names = []

    # read in data and fill force_array
    for ind, current_file_name in enumerate(file_names):
        if verbose:
            print(f'reading file {current_file_name}')
        input_data = np.loadtxt(current_file_name, comments=('@', '#'))
        # test if inpufile is corrupted
        if input_data[:, 0].shape != time.shape:
            print(f'skip file {current_file_name}\n')
            print(f'shape is {input_data.shape}')
            continue
        force_array[ind, :] = input_data[:, 1]
        force_array_names.append(current_file_name)

    # removing rows with only zero
    force_array = force_array[~np.all(force_array == 0, axis=1)]

    return force_array, force_array_names


def pullf_to_force_array(
    file_names: ArrayLikeStr,
    verbose: bool = False,
    res: Int = 1,
) -> Tuple[Float2DArray, Float1DArray, ArrayLikeStr]:
    """
    Write data of GROMACS pullf.xvg in np.ndarray.

    Parameters
    ----------
    file_names : array_like
        Contains force file names.
    verbose : bool, optional.
        Enables verbose mode.
    res : int, optional.
        Unused and kept for compatibiltiy.

    Returns
    -------
    force_array : np.ndarray
        Contains constraint forces for all simulations.
    x : 1D np.ndarray
        Contains positions related to the time trace, in nm.
    time_step: float
        Timestep of the pullf.xvg file
    force_array_names : 1D np.ndarray
        Contains the file names corresponding to the force_array.
    """
    # read testfile and determine t
    from dcTMD.io import _read_testfile
    time = _read_testfile(file_names, verbose)

    # fill arrays with data
    force_array, force_array_names = _fill_force_array(
        time, file_names, verbose,
    )

    return force_array, time, force_array_names


def calc_dG(
    force_array: Float1DArray,
    temp: Float,
    time: Float1DArray,
    vel: Float,
    res: Int = 1,
) -> Tuple[Float1DArray, Float1DArray, Float1DArray]:
    """
    Estimate free energy via the dissipation correction.

    Three steps:
        1) calculate force ACF
        2) calculate W_diss from integrate force ACF
        3) calculate dG from  dG = <W> - W_diss

    Parameters
    ----------
    force_array : np.array
        Contains time traces of the constraint forces for each simulation.
    temp : float
        Simulation temperature in K.
    time : 1D np.array
        Times related to the entries in the force time trace in ps.
    vel : float
        Pulling velocity in nm/ps.
    res : int, optional
        Striding for returned quantities.

    Returns
    -------
    W_mean : 1D np.array
        Average work in kJ/mol.
    W_diss : 1D np.array
        Dissipative work in kJ/mol
    W_mean-W_diss : 1D np.array
        Free energy estimate in kJ/mol.
    """
    from scipy.integrate import cumulative_trapezoid
    R = 0.0083144598  # noqa: WPS111
    RT = R * temp  # [R]=[kJ/(mol K)]
    pos = time * vel
    """
    * force average: calculate < f_c (t) >_N.
    **Important:** this is an ensemble average over the trajectory ensemble N,
    not the time average over t
    """
    # average and variance over all trajectories in each time step
    force_mean = np.mean(force_array, axis=0)  # shape: (length_data)
    W_mean = cumulative_trapezoid(force_mean, pos, initial=0)
    # calculate $\delta f_c(t) = f_c(t) - \left< f_c (t) \right>_N$ for all t
    delta_force_array = force_array - force_mean
    """
    * optimized algorithm for numerical evaluation:
    * integrate: $\int_0^t dt' \delta f_c(t')$ for all $t'$
    * multiply by $\delta f_c(t)$ to yield $\int_0^t dt'\delta f_c(t) \delta
    * f_c(t')$ for $t$ with all $t' \leq t$ each
    * then calculate the ensemble average $\left< \int_0^t dt' \delta f_c(t)
    * \delta f_c(t') \right>$
    """
    int_delta_force_array = cumulative_trapezoid(
        delta_force_array,
        time,
        axis=-1,
        initial=0,
    )
    intcorr = np.multiply(delta_force_array, int_delta_force_array)
    """
    # same as:
    for n in range(N):
      for i in range(length_data):
          intcorr[n,i] = delta_force_array[n,i]*int_delta_force_array[n,i]
    """
    gamma = np.mean(intcorr, axis=0) / RT
    # * $W_{diss}$ from integration:
    print('calculating dissipative work...')
    W_diss = cumulative_trapezoid(gamma, pos, initial=0) * vel
    dG = W_mean[::res] - W_diss[::res]
    return W_mean[::res], W_diss[::res], dG


def calc_dG_and_friction(
    force_array: Float2DArray,
    temp: Float,
    time: Float1DArray,
    vel: Float,
    sigma: Float,
    res: Int = 1,
) -> Tuple[Float1DArray, ...]:
    """
    Calculate friction Gamma and free energy dG by dissipation correction.

    Steps:
        1) calculate force ACF
        2) calculate Gamma
        3) calculate W_diss from integrate force ACF
        4) calculate dG from  dG = <W> - W_diss

    force_array : np.ndarray
        Contains time traces of the constraint forces for each simulation.
    temp : float
        Simulation temperature in K.
    time : 1D np.ndarray
        Times related to the entries in the force time trace in ps.
    vel : float
        Pulling velocity in nm/ps.
    sigma : float, optional
        Width in nm of a Gaussian filter to smooth the friction factor.
    res : int
        Striding for returned quantities.

    Returns
    -------
    W_mean : 1D np.ndarray
        Average work in kJ/mol.
    W_diss : 1D np.ndarray
        Dissipative work in kJ/mol.
    W_mean-W_diss : 1D np.ndarray
        Free energy estimate in kJ/mol.
    gamma: 1D np.ndarray
        Fricion in kJ/mol/(nm^2/ps).
    gamma_smooth
        Smoothed fricion in kJ/mol/(nm^2/ps).
    """
    from dcTMD.utils import gaussfilter_friction
    from scipy.integrate import cumulative_trapezoid

    pos = time * vel
    x_length = pos[1] - pos[0]
    R = 0.0083144598  # noqa: WPS111
    RT = R * temp  # [R]=[kJ/(mol K)]
    """
    * force average: calculate < f_c (t) >_N.
    **Important:** this is an ensemble average over the trajectory ensemble N,
    not the time average over t
    """
    # average and variance over all trajectories in each time step
    force_mean = np.mean(force_array, axis=0)  # shape: (length_data)
    W_mean = cumulative_trapezoid(force_mean, pos, initial=0)
    # calculate $\delta f_c(t) = f_c(t) - \left< f_c (t) \right>_N$ for all t
    delta_force_array = force_array - force_mean  # shape: (N, length_data)

    # ~~~ evaluation
    """
    * optimized algorithm for numerical evaluation:
    * integrate: $\int_0^t dt' \delta f_c(t')$ for all $t'$
    * multiply by $\delta f_c(t)$ to yield $\int_0^t dt'\delta f_c(t)
    * \delta f_c(t')$ for $t$ with all $t' \leq t$ each then calculate the
    * ensemble average $\left< \int_0^t dt' \delta f_c(t) \delta f_c(t')
    * \right>$
    """
    int_delta_force_array = cumulative_trapezoid(
        delta_force_array,
        time,
        axis=-1,
        initial=0,
    )
    intcorr = np.multiply(delta_force_array, int_delta_force_array)
    """
    # similar to :
    for n in range(N):
      for i in range(length_data):
          intcorr[n,i] = delta_force_array[n,i]*int_delta_force_array[n,i]
    """
    gamma = np.mean(intcorr, axis=0) / RT
    gamma_smooth = gaussfilter_friction(gamma, sigma, x_length)
    """
    # * autocorrelation function evaluation:
    # * calculate $\left< \delta f_c(t) \delta f_c(t') \right>$ for the last
    # * $t$
    corr_set = np.zeros(np.shape(force_array))
    print("calculating and processing ACF...\n")
    for n in range(N):
        corr_set[n, :] = delta_force_array[n, :]*delta_force_array[n, -1]
    autocorr_set = np.mean(corr_set, axis=0)
    """
    # * $W_{diss}$ from integration:
    print('Calculating dissipative work...')
    W_diss = cumulative_trapezoid(gamma, pos, initial=0) * vel

    # Reduce resolution
    W_mean = W_mean[::res]
    W_diss = W_diss[::res]
    dG = W_mean - W_diss
    gamma = gamma[::res]
    gamma_smooth = gamma_smooth[::res]
    return W_mean, W_diss, dG, gamma, gamma_smooth


def memory_kernel(
    delta_force_array: Float2DArray,
    x_indices: Float1DArray,
) -> Float1DArray:
    """
    Calculate memory kernel at positions X "forward" in time.

    Calculate memory kernel at positions X "forward" in time
    from fluctuation-dissipation.
    see e.g. R. Zwanzig, “Nonequilibrium statistical mechanics”,
    Oxford University Press (2001).
    corr_set[n,i] = delta_force_array[n,i]*delta_force_array[n,int
    ((length_data*args.x)-1)]


    Parameters
    ----------
    delta_force_array : np.ndarray
        Calculted via delta_force_array = force_array - force_mean
        with force_mean = np.mean(force_array, axis=0), an ensemble everage
        in each time step.
    x_indices : np.ndarray
        Indices at which memory kernel is calculated.

    Returns
    -------
    corr_set : np.ndarray
        shape: (len(X), length_data)
        NaN are set to zero
    """
    _, length_data = delta_force_array.shape
    corr_set = np.zeros((len(x_indices), length_data))

    for ind, tt in enumerate(range(length_data)):
        entries = delta_force_array[:, tt:-2] * \
            delta_force_array[:, tt + 1:-1]  # noqa: WPS221
        corr_set[ind, tt:-2] = np.mean(
            entries,
            axis=0,
        )
    return corr_set
