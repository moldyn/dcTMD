# -*- coding: utf-8 -*-
"""
changelog:
    - rename Gamma to gamma
    - work_array_names added
    - error estimation for dG and gamma added
"""

__all__ = []


from typing import Optional
import numpy as np
import pandas as pd
import math
from scipy.integrate import cumtrapz
from typing import Optional
import gc


def pullf_to_work_array(file_names: list, vel: float, skip: int =17, 
                        verbose=False, res=1):
    """ Writes data of GROMACS pullf.xvg in np.array full_force_set
    
    input:
        
        file_names: 
            list or 1D np.array containing file_names
            len(file_names) = number_of_files
        
        skip: int
            number of lines to skip in pullf.xvg file
        
        vel: float
            pulling velocity in nm/ps
        
        res: int
            striding

    out: 
        
        work_array: 
            np.arrray which contains force data 
            shape (number_of_files, length_data[::res]) 

        x: 1D np.array
            in [nm]; len(t) = length_data[::res]

        t[1]: float
            timestep of the pullf.xvg file [ps]
            this information is needed later for gamma smoothing 
        
        work_array_names: 
            1D np.array which contains the file_names coorrespoonding
            to full_force_set. This is useful when a pathseparation needs
            to be performed.

    """
    # prepare array and testfile
    number_of_files = len(file_names)
    print("{} files found".format(number_of_files))
    if verbose:
        print("reading file {}".format(file_names[0]))
    test_file = pd.read_csv(file_names[0],
                            sep='\s+', 
                            header=None, 
                            skiprows=skip, 
                            dtype=float,
                            usecols=[0]
                            ).to_numpy()

    t = test_file
    x = t * vel
    
    length_data = len(t)
    red_length_data = len(t[::res])

    if verbose:
        print('length of pullf file is {}'.format(length_data))
        print('output length is {}'.format(red_length_data))
    #~~~~~~~~~~~~

    work_array = np.zeros((number_of_files, red_length_data))
    work_array_names = []
    
    # read in data and fill work_array
    for i, current_file_name in enumerate(file_names):
        if verbose:
            print("reading file {}".format(current_file_name))
        input_data = pd.read_csv(current_file_name,
                                 sep='\s+',
                                 header=None,
                                 skiprows=skip,
                                 dtype=float,
                                 #usecols=[1],
                                 ).to_numpy().T
        # test if file is corrupted
        if input_data.shape != test_file.shape:
            print("skip file {}".format(current_file_name))
            print("shape is {}".format(input_data.shape))
            continue
        # force in units kJ/(mol*nm), work in units kJ/mol
        work_array[i, :] = cumtrapz(input_data[1], x, initial=0)[::res]
        work_array_names.append(current_file_name)

    # removing rows with only zero
    work_array = work_array[~np.all(work_array == 0, axis=1)]

    return work_array, t[::res], work_array_names


def calc_dG(work_set: np.ndarray, T: float, errors=False, 
                N_resamples:Optional[int]=None):
    """ calculate dG from dG = <W> - 1/(2kbT) <dW^2> 
        with optional bootstrap errors

    For a bootstrapping error estimation, set errors=True and pass a 
    N_resamples. The latter defines how many times trajectories are randomly
    drawn with replacement. Hence, they may be drawn multiple times or not at 
    all. This is implemented via sklearn.utils.resample but may be replaced
    easily. Note: n_samples=None means that the new sample is as large as the 
    original sample. The standard deviation of the bootstrap sample results 
    serves as an estimated error.
    (1) Manually check convergence!
    (2) Supply sufficient data (many trajectories)!

    input:
        work_set: np.array
        
        T: float
            simulation temperature in [K]
        
        optional:
        errors: bool
        
        N_resamples: int

    out:
        W_mean: 1D np.array
            average work in [kJ/mol]

        W_diss: 1D np.array
            dissipative work in [kJ/mol]

        W_mean-W_diss: 1D np.array
            dissipation corrected work or \Delta G in [kJ/mol]

        optional:
        np.std(s_W_mean, axis=0): 1D np.array
            bootstrap error of average work in [kJ/mol]

        np.std(s_W_diss, axis=0): 1D np.array
            bootstrap error of dissipative work in [kJ/mol]
        
        np.std(s_dG, axis=0): 1D np.array
            bootstrap error of \Delta G in [kJ/mol]
    
    """
    from scipy.constants import R
    RT = R*T/1e3            #[R]=[kJ/(mol K)]
                    
    N, length_data = np.shape(work_set)

    W_mean = np.mean(work_set, axis=0)  # shape: (length_data); kJ/mol
    W_var = np.var(work_set, axis=0)    # shape: (length_data); (kJ/mol)^2
    W_diss =  1/(2*RT)*W_var
    if not(errors):
        return W_mean, W_diss, W_mean-W_diss
    else:
        from sklearn.utils import resample
        s_W_mean = np.empty((N_resamples, length_data))
        s_W_diss = np.empty((N_resamples, length_data))
        s_dG = np.empty((N_resamples, length_data))
        
        for i in range(N_resamples):
            if np.mod(i, 100)==0:
                print(i, end=' ', flush=True)
            re_work_set = resample(work_set, n_samples=None)
            # random_indices = np.random.randint(0, N, N)
            # re_work_set = work_set[random_indices]

            W_mean_re, W_diss_re, dG_re = calc_dG(re_work_set, T, errors=False)
            s_W_mean[i] = W_mean_re
            s_W_diss[i] = W_diss_re
            s_dG[i] = dG_re

        return W_mean, W_diss, W_mean-W_diss, np.std(s_W_mean, axis=0), \
            np.std(s_W_diss, axis=0), np.std(s_dG, axis=0)


def calc_friction(W_diss, vel, time_step): #calc_gamma_from_W_diss
    """
    TODO: check if time step implementation is coorect!!!

    Calculate friction.
    gamma = d/dx W_diss(x)

    input:
        W_diss: 1D np.array
            dissipative work from calc_dG()

        vel: float
            pulling velocity in [nm/ps]

        time_step: float
            time step of pullf.xvg files in [ps]

    out: 
        Gamma: 1D np.array
            in [kJ/mol/(nm^2/ps)]

    """
    #time_step in ps; x in nm, d/dx --> delta x
    delta_x = time_step * vel
    print('calculating friction\n\
            timestep = {}ps \ndelta_x = {}nm'.format(time_step, delta_x))
    gamma = np.diff(W_diss, prepend=W_diss[0]) / (delta_x * vel)

    return gamma


def calc_dG_and_friction(work_set: np.ndarray, T: float, vel: float,
                        time_step: float, sigma: float, errors=False, N_resamples=None):
    """
    Calculate \Delta G and friction from work set.
    dG = <W> - 1/2kbT <dW^2> 
    gamma = d/dx 1/2kbT <dW^2> = d/dx W_diss(x)  
    with optional bootstrap errors

    For a bootstrapping error estimation, set errors=True and pass a 
    N_resamples. The latter defines how many times trajectories are randomly
    drawn with replacement. Hence, they may be drawn multiple times or not at 
    all. This is implemented via sklearn.utils.resample but may be replaced
    easily. Note: n_samples=None means that the new sample is as large as the 
    original sample. The standard deviation of the bootstrap sample results 
    serves as an estimated error.
    (1) Manually check convergence!
    (2) Supply sufficient data (many trajectories)!

    input:
        W_diss: 1D np.array
            dissipative work from calc_dG()

        T: float
            simulation temperature in [K]

        vel: float
            pulling velocity in [nm/ps]

        time_step: float
            time step of pullf.xvg files in [ps] 

        optional:
        errors: bool
        
        N_resamples: int

    out: 
        W_mean: 1D np.array
            average work in [kJ/mol]

        W_diss: 1D np.array
            dissipative work in [kJ/mol]

        W_mean-W_diss: 1D np.array
            dissipation corrected work or \Delta G in [kJ/mol]
        
        Gamma: 1D np.array
            in [kJ/mol/(nm^2/ps)]

        optional:
        np.std(s_W_mean, axis=0): 1D np.array
            bootstrap error of average work in [kJ/mol]

        np.std(s_W_diss, axis=0): 1D np.array
            bootstrap error of dissipative work in [kJ/mol]
        
        np.std(s_dG, axis=0): 1D np.array
            bootstrap error of \Delta G in [kJ/mol]

        np.std(s_gamma, axis=0): 1D np.array
            bootstrap error of gamma in [kJ/mol/(nm^2/ps)]
    """
    from scipy.constants import R
    #RT = R*T/1e3  # TODO test this
    RT = 0.0083144598*T                     #[R]=[kJ/(mol K)]
    N, length_data = np.shape(work_set)
    
    # calculate dG
    W_mean = np.mean(work_set, axis=0)  # shape: (length_data); kJ/mol
    W_var = np.var(work_set, axis=0)    # shape: (length_data); (kJ/mol)^2
    W_diss =  1/(2*RT)*W_var

    # calculate friction
    # time_step in ps; x in nm, d/dx --> delta x
    delta_x = time_step * vel
    #print('calculating friction\n\
    #        timestep = {}ps \ndelta_x = {}nm'.format(time_step, delta_x))
    gamma = np.diff(W_diss, prepend=W_diss[0]) / (delta_x * vel)
    gamma_smooth = gausfilter_friction(gamma, sigma, time_step*vel)
    if not(errors):
        return W_mean, W_diss, W_mean-W_diss, gamma, gamma_smooth
    else:
        from sklearn.utils import resample
        s_W_mean = np.empty((N_resamples, length_data))
        s_W_diss = np.empty_like(s_W_mean)
        s_dG = np.empty_like(s_W_mean)
        s_gamma = np.empty_like(s_W_mean)
        s_gamma_smooth = np.empty_like(s_W_mean)
        
        n_samples = math.ceil(0.9 * N)
        print("start resampling; n_sample={}".format(n_samples))
        for i in range(N_resamples):
            if np.mod(i, 100)==0:
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
