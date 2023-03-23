# -*- coding: utf-8 -*-
"""Smooth the friction."""
__all__ = ['gaussfilter_friction']

import numpy as np

from dcTMD._typing import (
    Float,
    Float1DArray,
    Str
)


def gaussfilter_friction(
    friction: Float1DArray,
    pos: Float1DArray,
    sigma: Float,
    mode: Str = "reflect",
) -> Float1DArray:
    """
    Smoothes friction with a gaussian kernel and 'nearest' borders.

    Parameters
    ----------
    friction : 1d np.array
        Array that contains the friction.
    pos : 1d np.array
        Positions corresponding to entries in friction array in nm.
    sigma : float
        Standard deviation of gaussian kernel in nm.
    mode:
        options: ‘reflect’, ‘constant’, ‘nearest’, ‘mirror’, ‘wrap’
        The mode parameter determines how the input array is
        extended beyond its boundaries. Default is ‘reflect’.
        Behavior for each option see scipy.ndimage.gaussian_filter1d

    Returns
    -------
    blurred_friction : 1d np.array
        Smoothed friction.
    """
    from scipy.ndimage import gaussian_filter1d
    delta_x = pos[1] - pos[0]
    blur = np.ceil(sigma / delta_x).astype(int)
    return gaussian_filter1d(friction, sigma=blur, mode=mode)
