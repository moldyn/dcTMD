# -*- coding: utf-8 -*-
"""
"""
__all__ = []

import numpy as np
import matplotlib.pyplot as plt

def gausfilter_friction(x_length, frict, sigma):
    """
    smoothes friction data from NEQGamma
    """
    from math import ceil
    blur = ceil(sigma/x_length)
    blurred = np.zeros(frict.shape)
    #print(input_file_data[:, 0:10])
    from scipy.ndimage.filters import gaussian_filter
    blurred = gaussian_filter(frict, sigma=blur, mode='reflect')
    #blurred = np.where(blurred[1]<0, 0, blurred[1])

    return blurred


def smooth_friction(x, frict, sigma=False):
    from scipy.ndimage.filters import gaussian_filter
    """
    smoothes friction data from NEQGamma
    input:
        x: 1d np array

        frict: 1d np array

        sigma: float
            window size for smoothing in [nm] e.g. 0.1
        if sigma=False iterates over ascending sigma values until 
        frict > 0 or until sigma_max=0.1x 
        Sets valus below 0 if sigma_max is reached.'

    out:
        frict_smooth
    """    
    if sigma:
        frict_smooth = gausfilter_friction(x, frict, sigma)            
    else:
        j = [0.001, 0.005, 0.01, 0.05, 0.1]
        sigma = j * x[-1]
        frict_smooth = gausfilter_friction(x, frict, sigma[0])
        for s in sigma:
            frict_smooth = gausfilter_friction(x, frict, s)
            if frict_smooth > 0:
                sigma = s
                break
                
        if frict_smooth < 0:
            print('sigma={} is no sufficient. Setting values<0 to 0'.format(s))
            frict_smooth = np.where(frict_smooth<0, 0, frict_smooth)
            sigma = s
            
    return frict_smooth, sigma


def plot_dG(ax1, x, dG, dG_color='black'):
    # TODO: test this an debug..
    ax1.plot(x, dG, lw=1, color=dG_color, label=r'PMF$_{tot}$')
    ax1.set_ylabel(r'$\Delta$G [kJ/mol]')
    ax1.set_xlabel('pulling distance [nm]')
    ax1.grid(False)
    ax1.legend()
    return


def plot_Gamma_twinax(ax1, x, frict, frict_color='blue'):
    # TODO: test this an debug..
    ax2 = ax1.twinx()
    ax2.plot(x, frict, c=frict_color, lw=1, label=r'gauss filtered $\Gamma$')
    ax2.spines['right'].set_color(frict_color)
    ax2.set_ylabel(r'$\Gamma$ [kg mol$^{-1}$ ns$^{-1}]$',
                   color=frict_color)
    ax2.tick_params(axis='y', color=frict_color, labelcolor=frict_color)
    ax2.grid(False)
    return ax1, ax2




