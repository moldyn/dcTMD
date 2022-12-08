# -*- coding: utf-8 -*-

__all__ = ['pullf_to_work', 'calc_dG', 'calc_friction']

import numpy as np
from beartype import beartype
from beartype.typing import Any, Optional
from typing import List


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
    file_names : array_like
        Constains file names.
    vel : float
        Pulling velocity in nm/ps.
    verbose : bool, optional
        Enables verbose mode.
    res : int, optional
        Striding for returned quantities.

    Returns
    ------- 
    work_array : np.ndarrray 
        contains cumulative work for each trajectory
        shape (number_of_files, length_data[::res]) 
    work_array_names : 1D np.ndarray 
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
    file_names : list or 1D np.ndarray 
        Contains file names.
    vel : float
        Pulling velocity in nm/ps.
    verbose : bool, optional.
        Enables verbose mode.
    res : int, optional.
        Striding for returned quantities.

    Returns
    ------- 
    work_array : np.ndarrray 
        Contains cumulative work for each trajectory,
        shape (number_of_files, length_data[::res]).
    x : 1D np.ndarray
        Contains positions related to the time trace, in nm.
    t[1] : float
        Timestep of the pullf.xvg file ps, this information is needed later 
        for friction smoothing.
    work_array_names : 1D np.ndarray 
        Contains the file names corresponding to the work_array. 
    """
    # read testfile and determine t
    from dcTMD.io import _read_testfile
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
    work_set: np.ndarray
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
        Average work in kJ/mol.
    W_diss: 1D np.array
        Dissipative work in kJ/mol.
    W_mean-W_diss: 1D np.array
        Free energy estimate in kJ/mol.
    np.std(s_W_mean, axis=0): 1D np.array, optional
        Bootstrap error of average work in kJ/mol.
    np.std(s_W_diss, axis=0): 1D np.array, optional
        Bootstrap error of dissipative work in kJ/mol.
    np.std(s_dG, axis=0): 1D np.array, optional
        Bootstrap error of free energy in kJ/mol.
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
    work_set : np.ndarray
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
    np.std(s_W_mean, axis=0) : 1D np.ndarray
        Bootstrap error of average work in kJ/mol.
    np.std(s_W_diss, axis=0) : 1D np.ndarray
        Bootstrap error of dissipative work in kJ/mol.
    np.std(s_dG, axis=0) : 1D np.ndarray
        Bootstrap error of free energy in kJ/mol.
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
                  t: np.ndarray,
                  sigma: Optional[float] = None,
                  errors: bool = False,
                  T: Optional[float] = None,
                  N_resamples: Optional[int] = None,
                  work_set: Optional[np.array] = None,
                  verbose: bool = False,
                  ) -> Any:
    """
    Calculate friction from the dissipative work via a derivative. Smoothing 
    with a Gaussian kernel is possible via the sigma parameter.

    For a bootstrapping error estimation, set errors=True and pass N_resamples.
    The latter defines how many random sets of trajectories are drawn with 
    replacement. The new samples have the size of the original work set. Hence,
    in a given set, some trajectories may be drawn multiple times and others
    not at all. The standard deviation of the bootstrap sample results serves 
    as an estimated error.
    (1) Please ensure convergence.
    (2) Supply sufficient data, i.e. many trajectories.

    TODO: check if time step implementation is coorect!!!
            delta_x instead of time_step in  function!

    Parameters
    ----------
    W_diss : 1D np.ndarray
        Dissipative work from calc_dG()
    vel : float
        Pulling velocity in nm/ps.
    t : np.ndarray
        Time trace corresponding to force file entries.
    sigma : float, optional
        Width in nm of a Gaussian filter to smooth the friction factor.
    errors : bool, optional
        Decides if bootstrap errors are calculated and returned.
    T : float, optional
        Simulation temperature in K, needed for bootstrapping error analysis.
    N_resamples : int, optional
        Number of drawn resamples for bootstrapping error analysis.
    work_set : np.array, optional
        Numpy array of shape (trajectories, time), containing the 
        constraint work.
    verbose : bool, optional
        Enables verbose mode.

    Returns
    -------
    friction : 1D np.ndarray
        Friction factor in kJ/mol/(nm^2/ps).
    smooth_friction : 1D np.ndarray, optional
        Smoothed friction factor in kJ/mol/(nm^2/ps).
    s_friction : 1D np.ndarray, optional
        Bootstrap error of the friction factor in kJ/mol/(nm^2/ps).
    s_smooth_friction : 1D np.ndarray, optional
        Bootstrap error of the friction factor in kJ/mol/(nm^2/ps).
    """
    # time_step in ps; x in nm, d/dx --> delta x
    x = t * vel
    delta_x = x[1] - x[0]
    if verbose:
        print(f'calculating friction, delta_x: {delta_x}nm')
    friction = np.diff(W_diss, prepend=W_diss[0]) / (delta_x * vel)

    if sigma != None:
        from .utils import gaussfilter_friction
        smooth_friction = gaussfilter_friction(friction, x, sigma)
    if errors == True:
        if sigma == None:
            s_friction = _calc_friction_errors(work_set,
                                               T,
                                               vel,
                                               t,
                                               N_resamples,
                                               sigma,
                                               verbose,
                                               )
            return friction, s_friction
        else:
            s_friction, s_smooth_friction = _calc_friction_errors()
            return friction, s_friction, smooth_friction, s_smooth_friction


@beartype
def _calc_friction_errors(work_set: np.ndarray,
                          T: float,
                          vel: float,
                          t: np.ndarray,
                          N_resamples: int,
                          sigma: Optional[float] = None,
                          verbose: bool = False,
                          ) -> Any:
    """
    Calculate friction from the dissipative work via a derivative.
    TODO: check if time step implementation is coorect!!!
            delta_x instead of time_step in  function!

    Parameters
    ----------
    work_set : np.ndarray
        Numpy array of shape (trajectories, time), containing the 
        constraint work.
    T : float
        Simulation temperature in K.
    vel : float
        Pulling velocity in nm/ps.
    t : np.ndarray
        Time trace corresponding to force file entries.
    N_resamles : int
        Number of drawn resamples for bootstrapping error analysis.
    sigma : float, optional
        Standard deviation of gaussian kernel in nm, used for error of smoothed
        friction.
    verbose : bool, optional
        Enables verbose mode.

    Returns
    -------
    np.std(s_friction, axis=0) : 1D np.ndarray
        Bootstrap error of the friction in kJ/mol/(nm^2/ps).
    np.std(s_smooth_friction, axis=0) : 1D np.ndarray, optional
        Bootstrap error of the smoothed friction in kJ/mol/(nm^2/ps).
    """
    N, length_data = np.shape(work_set)

    s_friction = np.empty((N_resamples, length_data))
    if sigma != None:
        s_smooth_friction = np.empty_like(s_friction)

    for i in range(N_resamples):
        if verbose and np.mod(i, 200) == 0:
            print(i, end=' ', flush=True)
        random_indices = np.random.randint(0, N, N)
        re_work_set = work_set[random_indices]
        _, W_diss, _ = calc_dG(re_work_set,
                               T,
                               errors=False,
                               verbose=False,
                               )
        friction_s = calc_friction(W_diss,
                                   vel,
                                   t,
                                   sigma,
                                   errors=False,
                                   verbose=False,
                                   )
        if sigma == None:
            s_friction[i] = friction_s
        else:
            s_friction[i], s_smooth_friction[i] = friction_s
    if sigma == None:
        return np.std(s_friction, axis=0)
    else:
        return np.std(s_friction, axis=0), np.std(s_smooth_friction, axis=0)


@beartype
def calc_dG_and_friction(work_set: np.ndarray,
                         T: float,
                         vel: float,
                         t: np.ndarray,
                         sigma: Optional[float] = None,
                         verbose: bool = False,
                         ) -> Any:
    """
    Wrapper that calculates free energy and friction. May smooth the friction
    if sigma is provided.

    Parameters
    ----------
    work_set : np.ndarray, optional
        Numpy array of shape (trajectories, time), containing the 
        constraint work.
    T : float
        Simulation temperature in K.
    vel : float
        Pulling velocity in nm/ps.
    t : float
        Time trace corresponding to force file entries.
    sigma : float, optional
        Width in nm of a Gaussian filter to smooth the friction factor.
    verbose : bool, optional
        Enables verbose mode.

    Returns
    -------
    dG : np.ndarray   
        Free energy estimate in kJ/mol.
    friction : 1D np.ndarray
        Friction factor in kJ/mol/(nm^2/ps).
    smooth_friction : 1D np.ndarray, optional
        Smoothed friction factor in kJ/mol/(nm^2/ps).    
    """
    _, W_diss, dG = calc_dG(work_set,
                            T,
                            verbose,
                            )
    friction_s = calc_friction(W_diss,
                               vel,
                               t,
                               sigma,
                               verbose=verbose,
                               )
    if sigma == None:
        return dG, friction_s
    else:
        return dG, friction_s[0], friction_s[1]
