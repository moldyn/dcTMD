# -*- coding: utf-8 -*-

__all__ = ['pullf_to_work_array', 'calc_dG', 'calc_friction']

import math
import numpy as np
import pandas as pd
from beartype import beartype
from beartype.typing import Any, Dict, Optional


@beartype
def _read_testfile(file_names: list[str],
                   verbose: bool = False,
                   ) -> np.ndarray:
    """ 
    Read test force file to determine t and its length.

    Parameters
    ----------
    file_names : list or 1D np.array 
        Contains all force file names
    verbose: bool, optional
        Enables verbose mode.

    Returns
    -------
    t : np.array
        contains time
    """
    if verbose:
        print(f'using {file_names[0]} to initialize arrays')
    test_file = np.loadtxt(file_names[0], comments=['@', '#'])
    return test_file[:, 0]


@beartype
def _fill_work_array(t: np.ndarray,
                     file_names: list,
                     vel: float,
                     verbose: bool = False,
                     res: int = 1,
                     ) -> Any:
    """
    Helper that integrates the force files.

    Parameters
    ----------
    t[1] : float
        Timestep of the pullf.xvg file in ps,
        this information is needed later for friction smoothing.
    file_names : list or 1D np.array 
        Constains file names.
    vel : float
        pulling velocity in nm/ps
    verbose: bool, optional
        Enables verbose mode.
    res : int, optional
        Striding for returned quantities.
    Returns
    ------- 
    work_array : np.arrray 
        contains cumulative work for each trajectory
        shape (number_of_files, length_data[::res]) 
    work_array_names : 1D np.array 
        Contains the file names corresponding to the work_array.     
    """
    from scipy.integrate import cumulative_trapezoid
    # allocate memory
    work_array = np.zeros((len(file_names), len(t[::res])))
    work_array_names = []
    x = t * vel
    # read in data and fill work_array
    for i, current_file_name in enumerate(file_names):
        if verbose:
            print(f'reading file {current_file_name}')
        current_file_data = np.loadtxt(current_file_name, comments=['@', '#'])
        # test if file is corrupted, else add it
        if current_file_data[:, 0].shape != t.shape:
            print(f'skip file {current_file_name}')
            print(f'shape is {current_file_data.shape}')
        else:  # force in units kJ/(mol*nm), work in units kJ/mol
            work_array[i, :] = cumulative_trapezoid(current_file_data[:, 1],
                                                    x,
                                                    initial=0,
                                                    )[::res]
            work_array_names.append(current_file_name)
    # removing rows with only zero
    work_array = work_array[~np.all(work_array == 0, axis=1)]
    return work_array, work_array_names


@beartype
def pullf_to_work_array(file_names: list,
                        vel: float,
                        verbose: bool = False,
                        res: int = 1
                        ) -> Any:
    """ 
    Writes data of GROMACS pullf.xvg in np.array full_force_set

    Parameters
    ----------
    file_names : list or 1D np.array 
        containing file names
    vel : float
        pulling velocity in nm/ps
    verbose : bool, optional.
        Enables verbose mode.
    res : int, optional.
        striding for returned quantities

    Returns
    ------- 
    work_array : np.arrray 
        contains cumulative work for each trajectory
        shape (number_of_files, length_data[::res]) 
    x : 1D np.array
        in [nm]; len(t) = length_data[::res]
    t[1] : float
        timestep of the pullf.xvg file [ps]
        this information is needed later for gamma smoothing 
    work_array_names : 1D np.array 
        Contains the file names corresponding to the work_array. 
    """
    # read testfile and determine t
    t = _read_testfile(file_names, verbose)
    if verbose:
        print(f'length of pullf file is {len(t)}')
        print(f'output length is {len(t[::res])}')
    # fill arrays with data
    work_array, work_array_names = _fill_work_array(t,
                                                    file_names,
                                                    vel,
                                                    verbose,
                                                    res,
                                                    )
    return work_array, t[::res], work_array_names


@beartype
def calc_dG(work_set: np.ndarray,
            T: float,
            errors: bool = False,
            N_resamples: Optional[int] = None,
            verbose: bool = False,
            ) -> Any:
    """ 
    Calculate dG from dG = <W> - 1/(2kbT) <dW^2> 
    with optional bootstrap errors.

    For a bootstrapping error estimation, set errors=True and pass a 
    N_resamples. The latter defines how many times trajectories are randomly
    drawn with replacement. Hence, they may be drawn multiple times or not at 
    all. Note: n_samples=None means that the new sample is as large as the 
    original sample. The standard deviation of the bootstrap sample results 
    serves as an estimated error.
    (1) Manually check convergence!
    (2) Supply sufficient data (many trajectories)!

    Parameters
    ----------
    work_set: np.array
        Numpy array of shape (trajectories, time), containing the 
        constraint work.
    T: float
        Simulation temperature in K.
    errors: bool, optional
        Decides if bootstrap errors are calculated and returned.
    N_resamples: int, optional
        Number of drawn resamples for bootstrapping error analysis.
    verbose: bool, optional 
        Enables verbose mode.

    Returns
    -------
    W_mean: 1D np.array
        average work in kJ/mol
    W_diss: 1D np.array
        dissipative work in kJ/mol
    W_mean-W_diss: 1D np.array
        dissipation corrected work or \Delta G in kJ/mol
    np.std(s_W_mean, axis=0): 1D np.array, optional
        bootstrap error of average work in kJ/mol
    np.std(s_W_diss, axis=0): 1D np.array, optional
        bootstrap error of dissipative work in kJ/mol
    np.std(s_dG, axis=0): 1D np.array, optional
        bootstrap error of \Delta G in kJ/mol
    """
    try:
        from scipy.constants import R
        RT = R*T/1e3  # [R]=[kJ/(mol K)]
    except:
        RT = 0.0083144598*T

    N, length_data = np.shape(work_set)

    W_mean = np.mean(work_set, axis=0)  # shape: (length_data); kJ/mol
    W_var = np.var(work_set, axis=0)    # shape: (length_data); (kJ/mol)^2
    W_diss = 1/(2*RT)*W_var
    if not (errors):
        return W_mean, W_diss, W_mean-W_diss
    else:
        s_W_mean, s_W_diss, s_dG = _calc_dG_errors(work_set,
                                                   T,
                                                   N_resamples,
                                                   verbose
                                                   )
        return W_mean, W_diss, W_mean-W_diss, s_W_mean, s_W_diss, s_dG


@beartype
def _calc_dG_errors(work_set: np.ndarray,
                    T: float,
                    N_resamples: int,
                    verbose: bool = False,
                    ) -> Any:
    """
    Helper that contains the bootstrapping loop.

    Parameters
    ----------
    work_set : np.array
        Numpy array of shape (trajectories, time), containing the 
        constraint work.
    T : float
        Simulation temperature in K.
    N_resamples : int
        Number of drawn resamples for bootstrapping error analysis.
    verbose : bool, optional
        Enables verbose mode.

    Returns
    -------
    np.std(s_W_mean, axis=0) : 1D np.array, optional
        bootstrap error of average work in kJ/mol
    np.std(s_W_diss, axis=0) : 1D np.array, optional
        bootstrap error of dissipative work in kJ/mol
    np.std(s_dG, axis=0) : 1D np.array, optional
        bootstrap error of \Delta G in kJ/mol
    """
    N, length_data = np.shape(work_set)

    s_W_mean = np.empty((N_resamples, length_data))
    s_W_diss = np.empty((N_resamples, length_data))
    s_dG = np.empty((N_resamples, length_data))

    for i in range(N_resamples):
        if verbose and np.mod(i, 200) == 0:
            print(i, end=' ', flush=True)
        random_indices = np.random.randint(0, N, N)
        re_work_set = work_set[random_indices]
        s_W_mean[i], s_W_diss[i], s_dG[i] = calc_dG(re_work_set,
                                                    T,
                                                    errors=False,
                                                    )
    return np.std(s_W_mean, axis=0), np.std(s_W_diss, axis=0), \
        np.std(s_dG, axis=0)


@beartype
def calc_friction(W_diss: np.ndarray,
                  vel: float,
                  time_step: float,
                  errors: bool = False,
                  N_resamples: int = None,
                  workset: np.array = None,
                  verbose: bool = False,
                  ):
    """
    Calculate friction from the dissipative work via a derivative.
    TODO: check if time step implementation is coorect!!!
            delta_x instead of time_step in  function!

    Parameters
    ----------
    W_diss : 1D np.array
        Dissipative work from calc_dG()
    vel : float
        Pulling velocity in nm/ps.
    time_step: float
        Time step of pullf.xvg files in ps.
    errors: bool, optional
        Decides if bootstrap errors are calculated and returned.
    N_resamples: int, optional
        Number of drawn resamples for bootstrapping error analysis.
    work_set: np.array, optional
        Numpy array of shape (trajectories, time), containing the 
        constraint work.
    verbose: bool, optional
        Enables verbose mode.

    Returns
    -------
    Gamma : 1D np.array
        Friction factor in kJ/mol/(nm^2/ps).
    np.std(s_gamma, axis=0): 1D np.array, optional
        Bootstrap error of the friction factor Gamma in kJ/mol/(nm^2/ps).
    """
    # time_step in ps; x in nm, d/dx --> delta x
    delta_x = time_step * vel
    if verbose:
        print(f'calculating friction\n\
            timestep = {time_step}ps\n\
            delta_x = {delta_x}nm')
    if errors == False:
        return np.diff(W_diss, prepend=W_diss[0]) / (delta_x * vel)
    else:
        _calc_friction_errors()


@beartype
def _calc_friction_errors(W_diss: np.ndarray,
                          vel: float,
                          time_step: float,
                          N_resamples: int,
                          verbose: bool = False,
                          ):
    """
    Calculate friction from the dissipative work via a derivative.
    TODO: check if time step implementation is coorect!!!
            delta_x instead of time_step in  function!

    Parameters
    ----------
    work_set: np.array, optional
        Numpy array of shape (trajectories, time), containing the 
        constraint work.
    vel : float
        Pulling velocity in nm/ps.
    time_step : float
        Time step of pullf.xvg files in ps.
    N_resamles : int
        Number of drawn resamples for bootstrapping error analysis.
    verbose : bool, optional
        Enables verbose mode.

    Returns
    -------
    Gamma : 1D np.array
        Friction factor in kJ/mol/(nm^2/ps).
    np.std(s_gamma, axis=0): 1D np.array, optional
        Bootstrap error of the friction factor Gamma in kJ/mol/(nm^2/ps).
    """
    N, length_data = np.shape(work_set)

    s_W_diss = np.empty_like((N_resamples, length_data))
    s_gamma = np.empty_like(s_W_diss)
    s_gamma_smooth = np.empty_like(s_W_diss)

    for i in range(N_resamples):
        if np.mod(i, 100) == 0:
            print(i, end=' ', flush=True)
        random_indices = np.random.randint(0, N, N)
        re_work_set = work_set[random_indices]
        W_mean_re, W_diss_re, dG_re, gamma_re, gamma_smooth_re = \
            calc_dG_and_friction(re_work_set, T,
                                 vel, time_step,
                                 sigma, errors=False)
        s_W_mean[i] = W_mean_re
        s_W_diss[i] = W_diss_re
        s_dG[i] = dG_re
        s_gamma[i] = gamma_re
        s_gamma_smooth[i] = gamma_smooth_re

    return W_mean, W_diss, W_mean-W_diss, gamma, gamma_smooth,\
        np.std(s_W_mean, axis=0), np.std(s_W_diss, axis=0), \
        np.std(s_dG, axis=0), np.std(s_gamma, axis=0), \
        np.std(s_gamma_smooth, axis=0)


def calc_dG_and_friction(work_set: np.ndarray,
                         T: float,
                         vel: float,
                         time_step: float,
                         sigma: float,
                         errors: bool = False,
                         N_resamples: int = None
                         ) -> Any:
    """
    Calculate \Delta G and friction from work set.
    dG = <W> - 1/2kbT <dW^2> 
    gamma = d/dx 1/2kbT <dW^2> = d/dx W_diss(x)  
    with optional bootstrap errors

    For a bootstrapping error estimation, set errors=True and pass a 
    N_resamples. The latter defines how many times trajectories are randomly
    drawn with replacement. Hence, they may be drawn multiple times or not at 
    all. Note: n_samples=None means that the new sample is as large as the 
    original sample. The standard deviation of the bootstrap sample results 
    serves as an estimated error.
    (1) Manually check convergence!
    (2) Supply sufficient data (many trajectories)!

    Parameters
    ----------
    W_diss: 1D np.array
        dissipative work
    T: float
        simulation temperature in K
    vel: float
        pulling velocity in nm/ps
    time_step: float
        time step of pullf.xvg files in ps
    errors: bool, optional
        decides if bootstrap errors are calculated and returned
    N_resamples: int, optional
        number of drawn resamples for bootstrapping error analysis

    Returns
    -------
    W_mean: 1D np.array
        average work in kJ/mol
    W_diss: 1D np.array
        dissipative work in kJ/mol
    W_mean-W_diss: 1D np.array
        dissipation corrected work or \Delta G in kJ/mol
    Gamma: 1D np.array
        in kJ/mol/(nm^2/ps)
    np.std(s_W_mean, axis=0): 1D np.array, optional
        bootstrap error of average work in kJ/mol
    np.std(s_W_diss, axis=0): 1D np.array, optional
        bootstrap error of dissipative work in kJ/mol
    np.std(s_dG, axis=0): 1D np.array, optional
        bootstrap error of \Delta G in kJ/mol
    np.std(s_gamma, axis=0): 1D np.array, optional
        bootstrap error of gamma in kJ/mol/(nm^2/ps)
    """
    from scipy.constants import R
    # RT = R*T/1e3  # TODO test this
    RT = 0.0083144598*T  # [R] = kJ/(mol K)
    N, length_data = np.shape(work_set)

    # calculate dG
    W_mean = np.mean(work_set, axis=0)  # shape: (length_data); kJ/mol
    W_var = np.var(work_set, axis=0)    # shape: (length_data); (kJ/mol)^2
    W_diss = 1/(2*RT)*W_var

    # calculate friction
    # time_step in ps; x in nm, d/dx --> delta x
    delta_x = time_step * vel
    gamma = np.diff(W_diss, prepend=W_diss[0]) / (delta_x * vel)
    gamma_smooth = gausfilter_friction(gamma, sigma, time_step*vel)
    if not (errors):
        return W_mean, W_diss, W_mean-W_diss, gamma, gamma_smooth
    else:
        from sklearn.utils import resample
        s_W_mean = np.empty((N_resamples, length_data))
        s_W_diss = np.empty_like(s_W_mean)
        s_dG = np.empty_like(s_W_mean)
        s_gamma = np.empty_like(s_W_mean)
        s_gamma_smooth = np.empty_like(s_W_mean)

        n_samples = math.ceil(0.9 * N)
        print(f'start resampling; n_sample={n_samples}')
        for i in range(N_resamples):
            if np.mod(i, 100) == 0:
                print(i, end=' ', flush=True)

            re_work_set = resample(work_set, n_samples=n_samples)
            # random_indices = np.random.randint(0, N, N)
            # re_work_set = work_set[random_indices]
            W_mean_re, W_diss_re, dG_re, gamma_re, gamma_smooth_re = \
                calc_dG_and_friction(re_work_set, T,
                                     vel, time_step,
                                     sigma, errors=False)
            s_W_mean[i] = W_mean_re
            s_W_diss[i] = W_diss_re
            s_dG[i] = dG_re
            s_gamma[i] = gamma_re
            s_gamma_smooth[i] = gamma_smooth_re

        return W_mean, W_diss, W_mean-W_diss, gamma, gamma_smooth,\
            np.std(s_W_mean, axis=0), np.std(s_W_diss, axis=0), \
            np.std(s_dG, axis=0), np.std(s_gamma, axis=0), \
            np.std(s_gamma_smooth, axis=0)


def gausfilter_friction(frict, sigma, x_length):
    """
    smoothes friction data from NEQGamma
    """
    from math import ceil
    from scipy.ndimage.filters import gaussian_filter
    blur = ceil(sigma/x_length)
    return gaussian_filter(frict, sigma=blur, mode='nearest')
