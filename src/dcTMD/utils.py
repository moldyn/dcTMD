# -*- coding: utf-8 -*-
"""
"""
__all__ = []

import numpy as np


def gaussfilter_friction(friction: np.ndarray,
                         x: np.ndarray,
                         sigma: float,
                         ) -> np.ndarray:
    """
    Smoothes friction with a gaussian kernel and reflective borders.

    Parameters
    ----------
    friction : 1d np.array
        Array that contains the friction.
    x : 1d np.array 
        Positions corresponding to entries in friction array in nm.
    sigma : float
        Standard deviation of gaussian kernel in nm.

    Returns
    -------
    blurred : 1d np.array
        Smoothed friction.
    """
    from scipy.ndimage.filters import gaussian_filter
    delta_x = x[1] - x[0]
    blur = np.ceil(sigma / delta_x).astype(int)
    blurred_friction = gaussian_filter(friction, sigma=blur, mode='reflect')
    return blurred_friction


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
