# -*- coding: utf-8 -*-
"""
"""
__all__ = []


import numpy as np
from beartype import beartype
from beartype.typing import Any, Optional
from typing import List



@beartype
def _fill_force_array(t: np.ndarray,
                      file_names: List[str],
                      verbose=False,
                      ) -> Any:
    """ 
    Fill the force array by reading in force files.

    Parameters
    ----------
    t : np.ndarray
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
    force_array = np.zeros((len(file_names), len(t)))
    force_array_names = []

    # read in data and fill force_array
    for i, current_file_name in enumerate(file_names):
        if verbose:
            print("reading file {}".format(current_file_name))
        input_data = np.loadtxt(current_file_name, comments=("@", "#"))
        # test if inpufile is corrupted
        if input_data[:, 0].shape != t.shape:
            print("skip file {}\n".format(current_file_name))
            print("shape is {}".format(input_data.shape))
            continue
        force_array[i, :] = input_data[:, 1]
        force_array_names.append(current_file_name)

    # removing rows with only zero
    force_array = force_array[~np.all(force_array == 0, axis=1)]

    return force_array, force_array_names


def pullf_to_force_array(file_names: List[str],
                         verbose=False,
                         res=1,
                         ) -> Any:
    """ 
    Write data of GROMACS pullf.xvg in np.ndarray.

    Parameters
    ----------
    file_names : array_like
        Contains force file names.
    verbose : bool, optional.
        Enables verbose mode.
    res : int, optional.
        Striding for returned quantities.

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
    t = _read_testfile(file_names, verbose)

    if verbose:
        print('length of pullf file is {}'.format(len(t)))
        print('output length is {}'.format(len(t[::res])))

    # fill arrays with data
    force_array, force_array_names = _fill_force_array(t, file_names, verbose)

    return force_array, t, force_array_names


def calc_dG(force_array: np.ndarray,
            T: float,
            t: np.ndarray,
            vel: float,
            res: int,
            ) -> Any:
    """ 
    Estimates free energy via the dissipation correction in three steps:
        1) calculate force ACF
        2) calculate W_diss from integrate force ACF
        3) calculate dG from  dG = <W> - W_diss

    Parameters
    ----------
    force_array : np.array
        Contains time traces of the constraint forces for each simulation.
    T : float
        Simulation temperature in K.
    t : 1D np.array
        Times related to the entries in the force time trace in ps.
    vel : float
        Pulling velocity in nm/ps.
    res : int
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
    RT = 0.0083144598*T  # [R]=[kJ/(mol K)]
    x = t * vel
    """
    * force average: calculate < f_c (t) >_N.
    **Important:** this is an ensemble average over the trajectory ensemble N,
    not the time average over t
    """
    # average and variance over all trajectories in each time step
    force_mean = np.mean(force_array, axis=0)  # shape: (length_data)
    W_mean = cumulative_trapezoid(force_mean, x, initial=0)
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
    int_delta_force_array = cumulative_trapezoid(delta_force_array, t,
                                                 axis=-1, initial=0)
    intcorr = np.multiply(delta_force_array, int_delta_force_array)
    """
    # same as:
    for n in range(N):
      for i in range(length_data):
          intcorr[n,i] = delta_force_array[n,i]*int_delta_force_array[n,i]
    """
    gamma = np.mean(intcorr, axis=0) / RT
    # * $W_{diss}$ from integration:
    print("calculating dissipative work...\n")
    W_diss = cumulative_trapezoid(gamma, x, initial=0) * vel
    return W_mean[::res], W_diss[::res], W_mean[::res] - W_diss[::res]


def calc_dG_and_friction(force_array: np.ndarray,
                         T: float,
                         t: np.ndarray,
                         vel: float,
                         sigma: float,
                         res: int = 1,
                         ) -> Any:
    """ calculate Gamma and dG dissipation correction
        1) calculate force ACF
        2) calculate Gamma
        3) calculate W_diss from integrate force ACF
        4) calculate dG from  dG = <W> - W_diss

    force_array : np.ndarray
        Contains time traces of the constraint forces for each simulation.
    T : float
        Simulation temperature in K.
    t : 1D np.ndarray
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
    from .utils import gaussfilter_friction
    from scipy.integrate import cumulative_trapezoid

    x = t * vel
    x_length = x[1]-x[0]

    RT = 0.0083144598*T  # [R]=[kJ/(mol K)]
    """
    * force average: calculate < f_c (t) >_N.
    **Important:** this is an ensemble average over the trajectory ensemble N,
    not the time average over t
    """
    # average and variance over all trajectories in each time step
    force_mean = np.mean(force_array, axis=0)  # shape: (length_data)
    W_mean = cumulative_trapezoid(force_mean, x, initial=0)
    # calculate $\delta f_c(t) = f_c(t) - \left< f_c (t) \right>_N$ for all t
    delta_force_array = force_array - force_mean  # shape: (N, length_data)

    # ~~~ evaluation
    """
    * optimized algorithm for numerical evaluation:
    * integrate: $\int_0^t dt' \delta f_c(t')$ for all $t'$
    * multiply by $\delta f_c(t)$ to yield $\int_0^t dt'\delta f_c(t) \delta f_c(t')$ for $t$ with all $t' \leq t$ each
    * then calculate the ensemble average $\left< \int_0^t dt' \delta f_c(t) \delta f_c(t') \right>$
    """
    int_delta_force_array = cumulative_trapezoid(delta_force_array, t,
                                                 axis=-1, initial=0)
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
    # * calculate $\left< \delta f_c(t) \delta f_c(t') \right>$ for the last $t$
    corr_set = np.zeros(np.shape(force_array))
    print("calculating and processing ACF...\n")
    for n in range(N):
        corr_set[n, :] = delta_force_array[n, :]*delta_force_array[n, -1]
    autocorr_set = np.mean(corr_set, axis=0)
    """
    # * $W_{diss}$ from integration:
    print("calculating dissipative work...\n")
    W_diss = cumulative_trapezoid(gamma, x, initial=0) * vel

    return W_mean[::res], W_diss[::res], W_mean[::res]-W_diss[::res], \
        gamma[::res], gamma_smooth[::res]


def memory_kernel(delta_force_array: np.ndarray,
                  X: np.ndarray,
                  ) -> np.ndarray:
    """
    Calculate memory kernel at positions X "forward" in time
    from fluctuation-dissipation 
    see e.g. R. Zwanzig, “Nonequilibrium statistical mechanics”, 
    Oxford University Press (2001).

    Parameters
    ----------
    delta_force_array : np.ndarray
        Calculted via delta_force_array = force_array - force_mean
        with force_mean = np.mean(force_array, axis=0), an ensemble everage 
        in each time step.
    X : np.ndarray
        Indices at which memory kernel is calculated.

    Returns
    -------
    corr_set : np.ndarray
        shape: (len(X), length_data)
        NaN are set to zero

        corr_set[n,i] = delta_force_array[n,i]*delta_force_array[n,int
        ((length_data*args.x)-1)]
    """
    N, length_data = delta_force_array.shape
    corr_set = np.zeros((len(X), length_data))

    for i, t in enumerate(range(length_data)):
        corr_set[i, t:-2] = np.mean(delta_force_array[:, t:-2]
                                    * delta_force_array[:, t+1:-1],
                                    axis=0)
    return corr_set
