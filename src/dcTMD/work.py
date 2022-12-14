# -*- coding: utf-8 -*-
"""Calculate dcTMD quantities via the cumulants of the work."""

__all__ = ['pullf_to_work_array', 'calc_dG', 'calc_friction']

import numpy as np
from beartype import beartype
from beartype.typing import Union, Optional, Tuple

from dcTMD._typing import (
    Int,
    Float,
    ArrayLikeStr,
    Float1DArray,
    Float2DArray,
)


@beartype
def _fill_work_array(
    time: Float1DArray,
    file_names: ArrayLikeStr,
    vel: Float,
    verbose: bool = False,
    res: Int = 1,
) -> Tuple[Float2DArray, ArrayLikeStr]:
    """
    Help integrate the force files.

    Force in units kJ/(mol*nm), work in units kJ/mol.

    Parameters
    ----------
    time : np.ndarray
        Time trace corresponding to force file entries.
    file_names : array_like
        Contains file names.
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
    import tqdm
    from scipy.integrate import cumulative_trapezoid
    work_array = np.zeros((len(file_names), len(time[::res])))
    work_array_names = []
    pos = time * vel
    # read in data and fill work_array
    for ind, current_file_name in tqdm.tqdm(
        enumerate(file_names),
        total=len(file_names),
        desc='Loading & integrating force files',
    ):
        if verbose:
            print(f'reading file {current_file_name}')
        current_file_data = np.loadtxt(current_file_name, comments=['@', '#'])
        # test if file is corrupted, else add it
        if current_file_data[:, 0].shape == time.shape:
            work_array[ind, :] = cumulative_trapezoid(
                current_file_data[:, 1],
                pos,
                initial=0,
            )[::res]
            work_array_names.append(current_file_name)
        else:
            print(f'skip file {current_file_name}')
            print(f'shape is {current_file_data.shape}')
    # removing rows with only zero
    work_array = work_array[~np.all(work_array == 0, axis=1)]
    return work_array, work_array_names


@beartype
def pullf_to_work_array(
    file_names: ArrayLikeStr,
    vel: Float,
    verbose: bool = False,
    res: Int = 1,
) -> Tuple[Float2DArray, Float1DArray, ArrayLikeStr]:
    """
    Integrate forces from GROMACS pullf.xvg files.

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
    time : np.ndarray
        Time trace corresponding to force file entries.
    work_array_names : 1D np.ndarray
        Contains the file names corresponding to the work_array.
    """
    # read testfile and determine t
    from dcTMD.io import _read_testfile
    time = _read_testfile(file_names, verbose)
    if verbose:
        time_length = len(time)
        time_length_reduced = len(time[::res])
        print(f'length of pullf file is {time_length}')
        print(f'reduced length is {time_length_reduced}')
    # fill arrays with data
    work_array, work_array_names = _fill_work_array(
        time,
        file_names,
        vel,
        verbose,
        res,
    )
    return work_array, time[::res], work_array_names


@beartype
def calc_dG(
    work_set: Float2DArray,
    temp: Float,
    N_resamples: Optional[Int] = None,
    verbose: bool = False,
) -> Tuple[Float1DArray, ...]:
    """
    Estimate free energy with cumulant expansion & optional bootstrapping.

    For a bootstrapping error estimation, set pass N_resamples. The latter
    defines how many random sets of trajectories are drawn with replacement.
    The new samples have the size of the original work set. Hence, in a given
    set, some trajectories may be drawn multiple times and others not at all.
    The standard deviation of the bootstrap sample results serves as an
    estimated error.
    (1) Please ensure convergence.
    (2) Supply sufficient data, i.e. many trajectories.
    (3) Inspect the distributions.

    Parameters
    ----------
    work_set: np.ndarray
        Numpy array of shape (trajectories, time), containing the
        constraint work.
    temp: float
        Simulation temperature in K.
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
        from scipy.constants import R  # noqa: WPS347
    except ModuleNotFoundError:
        RT = 8.3144598
    RT = R * temp / 1e3

    W_mean = np.mean(work_set, axis=0)  # shape: (length_data); kJ/mol
    W_var = np.var(work_set, axis=0)    # shape: (length_data); (kJ/mol)^2
    W_diss = 1 / (2 * RT) * W_var
    dG = W_mean - W_diss
    if N_resamples is not None:
        s_W_mean, s_W_diss, s_dG = _calc_dG_errors(
            work_set,
            temp,
            N_resamples,
            verbose,
        )
        return W_mean, W_diss, dG, s_W_mean, s_W_diss, s_dG  # noqa: WPS227
    return W_mean, W_diss, dG


@beartype
def _calc_dG_errors(
    work_set: Float2DArray,
    temp: Float,
    N_resamples: Int,
    verbose: bool = False,
) -> Tuple[Float1DArray, ...]:
    """
    Perform the bootstrapping error estimation.

    Parameters
    ----------
    work_set : np.ndarray
        Numpy array of shape (trajectories, time), containing the
        constraint work.
    temp : float
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
    import tqdm
    N_traj, length_data = np.shape(work_set)

    s_W_mean = np.empty((N_resamples, length_data))
    s_W_diss = np.empty((N_resamples, length_data))
    s_dG = np.empty((N_resamples, length_data))

    for ind in tqdm.tqdm(
        range(N_resamples),
        desc='Bootstrapping progress',
    ):
        random_indices = np.random.randint(0, N_traj, N_traj)
        re_work_set = work_set[random_indices]
        W_mean, W_diss, dG = calc_dG(
            re_work_set,
            temp,
            N_resamples=None,
        )
        s_W_mean[ind] = W_mean
        s_W_diss[ind] = W_diss
        s_dG[ind] = dG
    return np.std(s_W_mean, axis=0), \
        np.std(s_W_diss, axis=0), \
        np.std(s_dG, axis=0)


@beartype
def calc_friction(
    W_diss: Float1DArray,
    vel: Float,
    time: Float1DArray,
    sigma: Optional[Float] = None,
    N_resamples: Optional[Int] = None,
    temp: Optional[Float] = None,
    work_set: Optional[Float2DArray] = None,
    verbose: bool = False,
) -> Union[Tuple[Float1DArray, ...], Float1DArray]:
    """
    Calculate (& smooth) friction from the dissipative work via its derivative.

    For a bootstrapping error estimation, set pass N_resamples. The latter
    defines how many random sets of trajectories are drawn with replacement.
    The new samples have the size of the original work set. Hence, in a given
    set, some trajectories may be drawn multiple times and others not at all.
    The standard deviation of the bootstrap sample results serves as an
    estimated error.
    (1) Please ensure convergence.
    (2) Supply sufficient data, i.e. many trajectories.
    (3) Inspect the distributions.

    TODO: check if time step implementation is coorect!!!
            delta_x instead of time_step in  function!

    Parameters
    ----------
    W_diss : 1D np.ndarray
        Dissipative work from calc_dG()
    vel : float
        Pulling velocity in nm/ps.
    time : np.ndarray
        Time trace corresponding to force file entries.
    sigma : float, optional
        Width in nm of a Gaussian filter to smooth the friction factor.
    N_resamples : int, optional
        Number of drawn resamples for bootstrapping error analysis.
    temp : float, optional
        Simulation temperature in K, needed for bootstrapping error analysis.
    work_set : np.array, optional
        Numpy array of shape (trajectories, time), containing the
        constraint work, needed for bootstrapping error analysis.
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
    returns = ()
    # time_step in ps; pos=x in nm, d/dx --> delta x
    pos = time * vel
    delta_x = pos[1] - pos[0]
    if verbose:
        print(f'calculating friction, delta_x: {delta_x}nm')
    friction = np.diff(W_diss, prepend=W_diss[0]) / (delta_x * vel)
    returns = (*returns, friction)

    if sigma is not None:
        from dcTMD.utils import gaussfilter_friction
        returns = (*returns, gaussfilter_friction(friction, pos, sigma))
    if N_resamples is not None:
        returns = (*returns, *_calc_friction_errors(
            work_set,
            temp,
            vel,
            time,
            N_resamples,
            sigma,
            verbose,
        ),
        )
    return returns if len(returns) > 1 else returns[0]


@beartype
def _calc_friction_errors(
    work_set: Float2DArray,
    temp: Float,
    vel: Float,
    time: Float1DArray,
    N_resamples: Int,
    sigma: Optional[Float] = None,
    verbose: bool = False,
) -> Tuple[Float1DArray, ...]:
    """
    Calculate friction from the dissipative work via a derivative.

    TODO: check if time step implementation is coorect!!!
    delta_x instead of time_step in  function!

    Parameters
    ----------
    work_set : np.ndarray
        Numpy array of shape (trajectories, time), containing the
        constraint work.
    temp : float
        Simulation temperature in K.
    vel : float
        Pulling velocity in nm/ps.
    time : np.ndarray
        Time trace corresponding to force file entries.
    N_resamples : int
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
    import tqdm
    N_traj, length_data = np.shape(work_set)

    s_friction = np.empty((N_resamples, length_data))
    if sigma is not None:
        s_smooth_friction = np.empty_like(s_friction)

    for ind in tqdm.tqdm(
        range(N_resamples),
        desc='Bootstrapping progress',
    ):
        random_indices = np.random.randint(0, N_traj, N_traj)
        re_work_set = work_set[random_indices]
        _, W_diss, _ = calc_dG(
            re_work_set,
            temp,
            N_resamples=None,
            verbose=False,
        )
        frictions = calc_friction(
            W_diss,
            vel,
            time,
            sigma,
            N_resamples=None,
            verbose=False,
        )
        if sigma is None:
            s_friction[ind] = frictions
        else:
            s_friction[ind], s_smooth_friction[ind] = frictions  # noqa: WPS414
    if sigma is None:
        return (np.std(s_friction, axis=0), )
    return np.std(s_friction, axis=0), np.std(s_smooth_friction, axis=0)


@beartype
def calc_dG_and_friction(
    work_set: Float2DArray,
    temp: Float,
    vel: Float,
    time: Float1DArray,
    sigma: Optional[Float] = None,
    verbose: bool = False,
) -> Tuple[Float1DArray, ...]:
    """
    Calculate (& smooth) free energy and friction.

    May smooth the friction if sigma is provided.

    Parameters
    ----------
    work_set : np.ndarray, optional
        Numpy array of shape (trajectories, time), containing the
        constraint work.
    temp : float
        Simulation temperature in K.
    vel : float
        Pulling velocity in nm/ps.
    time : float
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
    _, W_diss, dG = calc_dG(
        work_set,
        temp,
        verbose,
    )
    friction_s = calc_friction(
        W_diss,
        vel,
        time,
        sigma,
        verbose=verbose,
    )
    if sigma is None:
        return dG, friction_s
    return dG, friction_s[0], friction_s[1]
