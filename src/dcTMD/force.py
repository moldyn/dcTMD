# -*- coding: utf-8 -*-
"""
"""
__all__ = []


import scipy.integrate
from scipy.integrate import cumtrapz
import numpy as np
import pandas as pd



def pullf_to_force_array(file_names: list, vel: float, skip: int =17,
                             verbose=False, res=1):
    """ Writes data of GROMACS pullf.xvg in np.array force_array
    
    input:
        
        file_names: 
            list or 1D np.array containing file_names
            len(file_names) = number_of_files
        
        skip: 
            number of lines to skip in pullf.xvg file

    out: 
        
        force_array: 
            np.arrray which contains force data 
            shape: (number_of_files, length_data)

        x: 1D np.array
            in [nm]; len(t) = length_data[::res]

        time_step: 
            timestep of the pullf.xvg file
        
        force_array_names: 
            1D np.array which contains the file_names coorrespoonding
            to force_array. This is useful when a pathseparation needs
            to be performed.

    """
    # prepare np.array
    number_of_files = len(file_names)
    print("{} files found".format(number_of_files))
    if verbose:
        print("reading file {}".format(file_names[0]))
    test_file = pd.read_csv(file_names[0],
                            sep='\s+', 
                            header=None, 
                            skiprows=skip, 
                            dtype=float
                            ).to_numpy().T

    t = test_file[0]
    length_data = len(t)

    if verbose:
        print('length of pullf file is {}'.format(length_data))
    
    force_array = np.zeros((number_of_files, length_data))
    force_array_names = np.empty(number_of_files, dtype=str)

    # read in data and fill force_array
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
        # test if inpufile is corrupted
        if input_data.shape != test_file.shape:
            print("skip file {}\n".format(current_file_name))
            print("shape is {}".format(input_data.shape))
            continue
        force_array[i, :] = input_data[1]
        force_array_names[i] = current_file_name

    # removing rows with only zero
    force_array = force_array[~np.all(force_array == 0, axis=1)]
    # TODO: fix this
    #force_array_names = force_array_names[~np.all(force_array_names == '')]

    return force_array, t, force_array_names


def calc_dG(force_array: np.ndarray, T: float, t, vel: float, res:int):
    """ calculate dG dissipation correction
        1) calculate force ACF
        2) calculate W_diss from integrate force ACF
        3) calculate dG from  dG = <W> - W_diss
       
    input:
        force_array: np.array

        T: float
            simulation temperature in [K]

        t: 1D np.array

        vel: float
            pulling velocity in [nm/ps]

        res: int
            striding
   
    out:
        W_mean: 1D np.array
            average work in [kJ/mol]

        W_diss: 1D np.array
            dissipative work in [kJ/mol]

        W_mean-W_diss: 1D np.array
            dissipation corrected work or \Delta G in [kJ/mol]
    """

    from scipy.constants import R
    #RT = R*T/1e3  # TODO test this
    RT = 0.0083144598*T                     #[R]=[kJ/(mol K)]
    x = t * vel
    """
    * force average: calculate < f_c (t) >_N.
    **Important:** this is an ensemble average over the trajectory ensemble N,
    not the time average over t
    """
    # average and variance over all trajectories in each time step
    force_mean = np.mean(force_array, axis=0)  # shape: (length_data)
    W_mean = cumtrapz(force_mean, x, initial=0)
    # calculate $\delta f_c(t) = f_c(t) - \left< f_c (t) \right>_N$ for all t
    delta_force_array = force_array - force_mean

    """
    * optimized algorithm for numerical evaluation:
    * integrate: $\int_0^t dt' \delta f_c(t')$ for all $t'$
    * multiply by $\delta f_c(t)$ to yield $\int_0^t dt'\delta f_c(t) \delta f_c(t')$ for $t$ with all $t' \leq t$ each
    * then calculate the ensemble average $\left< \int_0^t dt' \delta f_c(t) \delta f_c(t') \right>$
    """
    int_delta_force_array = cumtrapz(delta_force_array, t,
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
    W_diss = cumtrapz(gamma, x, initial=0) * vel
    
    return W_mean[::res], W_diss[::res], W_mean[::res] - W_diss[::res]


def calc_dG_and_friction(force_array: np.ndarray, T: float, t, vel: float,
                        sigma: float, res: int=1):
    """ calculate Gamma and dG dissipation correction
        1) calculate force ACF
        2) calculate Gamma
        3) calculate W_diss from integrate force ACF
        4) calculate dG from  dG = <W> - W_diss
       
    input:
        force_array: np.array

        T: float
            simulation temperature in [K]

        t: 1D np.array

        vel: float
            pulling velocity in [nm/ps]
   
    out:
        W_mean: 1D np.array
            average work in [kJ/mol]

        W_diss: 1D np.array
            dissipative work in [kJ/mol]

        W_mean-W_diss: 1D np.array
            dissipation corrected work or \Delta G in [kJ/mol]

        Gamma: 1D np.array
            in [kJ/mol/(nm^2/ps)]
    """

    x = t * vel
    x_length = x[1]-x[0]
    N, length_data = np.shape(force_array)
    from scipy.constants import R
    #RT = R*T/1e3  # TODO test this
    RT = 0.0083144598*T                     #[R]=[kJ/(mol K)]
    """
    * force average: calculate < f_c (t) >_N.
    **Important:** this is an ensemble average over the trajectory ensemble N,
    not the time average over t
    """
    # average and variance over all trajectories in each time step
    force_mean = np.mean(force_array, axis=0)  # shape: (length_data)
    W_mean = scipy.integrate.cumtrapz(force_mean, x, initial=0)
    # calculate $\delta f_c(t) = f_c(t) - \left< f_c (t) \right>_N$ for all t
    delta_force_array = force_array - force_mean

    # ~~~ evaluation
    """
    * optimized algorithm for numerical evaluation:
    * integrate: $\int_0^t dt' \delta f_c(t')$ for all $t'$
    * multiply by $\delta f_c(t)$ to yield $\int_0^t dt'\delta f_c(t) \delta f_c(t')$ for $t$ with all $t' \leq t$ each
    * then calculate the ensemble average $\left< \int_0^t dt' \delta f_c(t) \delta f_c(t') \right>$
    """
    int_delta_force_array = scipy.integrate.cumtrapz(delta_force_array, t,
                                                    axis=-1, initial=0)
    intcorr = np.multiply(delta_force_array, int_delta_force_array)
    """
    # similar to :
    for n in range(N):
      for i in range(length_data):
          intcorr[n,i] = delta_force_array[n,i]*int_delta_force_array[n,i]
    """
    gamma = np.mean(intcorr, axis=0) / RT
    gamma_smooth = gausfilter_friction(gamma, sigma, x_length)
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
    W_diss = scipy.integrate.cumtrapz(gamma, x, initial=0) * vel
    
    return W_mean[::res], W_diss[::res], W_mean[::res]-W_diss[::res], \
            gamma[::res], gamma_smooth[::res]


def memory_kernel():
    """
    TODO
    calculate memory kernel at position x
    corr_set[n,i] = delta_force_array[n,i]*delta_force_array[n,int((length_data*args.x)-1)]
    """

def gausfilter_friction(frict, sigma, x_length):
    """
    smoothes friction data from NEQGamma
    """
    import math
    blur = math.ceil(sigma/x_length)
    print(blur)
    blurred = np.zeros(frict.shape)
    #print(input_file_data[:, 0:10])
    from scipy.ndimage.filters import gaussian_filter
    blurred = gaussian_filter(frict, sigma=blur, mode='reflect')
    #blurred = np.where(blurred[1]<0, 0, blurred[1])

    return blurred