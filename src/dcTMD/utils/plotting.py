# -*- coding: utf-8 -*-
# MIT License
# Copyright (c) 2021-2023, Miriam Jäger
# All rights reserved.
"""Simple plot functions for dcTMD results."""

import sys
from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import norm


# Constants
LABEL_POSITION_NM = r'$x$ [nm]'
LABEL_FRICTION = r'$\Gamma$ [kJ nm$^2$/(mol ps)]'
UNIT_ENERGY = '[kJ/mol]'
LABEL_dG = r'$\Delta G$'
ALPHA_VALUE = 0.3  # transparency
WORKLINES_WIDTH = 0.5  # line width for work lines


def plot_dcTMD_results(
    estimator,
    friction=None,
    x=None,
    figsize=(4, 4),
):
    """Plot dcTMD results overview in two subplots.
    This function generates a figure with two subplots.
    The top subplot contains free energy, dissipative work, and mean work.
    The bottom subplot contains friction vs. position.

    Args:
        estimator:
            dcTMD work or force estimator object,
            containing the dcTMD results.
        friction (array-like, optional):
            Friction values to plot. If not provided,
            the function uses `estimator.friction_`. If `estimator` has a
            `friction_smooth_` attribute, it will use that instead.
        x (array-like, optional):
            The x-axis positions for the plots. If not
            provided, `estimator.position_` is used.
        figsize (tuple, optional):
            Size of the figure. Default is (4, 4).

    Returns:
        tuple: containing fig and axs

    Example:
        >>> plot_dcTMD_results(
        ...     estimator=my_estimator,
        ...     figsize=(6, 6),
        ... )
        >>> plt.show()
    """
    fig, axs = plt.subplots(ncols=1,
                            nrows=2,
                            sharex=True,
                            figsize=figsize,
                            )
    if x is None:
        x = estimator.position_
    plot_dG_Wdiss(estimator, axs[0], x=x)
    if friction is None:
        friction = estimator.friction_
        if hasattr(estimator, 'friction_smooth_'):
            friction = estimator.friction_smooth_
    plot_Gamma(x, friction, axs[1])
    axs[0].legend(
        loc='lower left',
        mode='expand',
        bbox_to_anchor=(0, 1.05, 1, 0.2),
        bbox_transform=axs[0].transAxes,
        frameon=False,
        ncol=3,
    )
    axs[0].set_xlabel("")
    plt.tight_layout()
    return fig, axs


def plot_dG_Wdiss(workestimator, ax, x=None):
    """Plot free energy, dissipative work and mean work vs position.

    Args:
        workestimator:
            dcTMD workestimator object,
            containing the dcTMD results.
        ax (matplotlib.axes.Axes):
            The axes object where the plot will be drawn.
        x (array-like, optional):
            The x-axis positions for the plot. If not provided,
            `workestimator.position_` is used.

    Returns:
        None

    Notes:
        - The x-axis represents the position along the coordinate (in nm).
        - The y-axis represents the energy values (in kJ/mol).

    Example:
        >>> fig, ax = plt.subplots()
        >>> plot_dG_Wdiss(workestimator=my_estimator, ax=ax)
        >>> plt.show()
    """
    if x is None:
        x = workestimator.position_
    ax.plot(x, workestimator.dG_, label=LABEL_dG)
    ax.plot(x, workestimator.W_mean_, label=r'W$_{\mathrm{mean}}$')
    ax.plot(x, workestimator.W_diss_, label=r'W$_{\mathrm{diss}}$')
    ax.set(xlabel=LABEL_POSITION_NM,
           ylabel=UNIT_ENERGY,
           xlim=[min(x), max(x)],
           )


def plot_Gamma(x, friction, ax, label=None):
    """Plot friction factor (Γ) vs position.

    Args:
        x (array-like):
            Positions along the coordinate (in nm).
        friction (array-like):
            Friction factor values (in kJ nm²/(mol ps)).
        ax (matplotlib.axes.Axes):
            The axes object where the plot will be drawn.
        label (str, optional):
            Label for the friction curve. Default is None.

    Example:
        >>> fig, ax = plt.subplots()
        >>> plot_Gamma(positions, friction_values, ax)
        >>> plt.legend()
        >>> plt.show()
    """
    ax.plot(x, friction, label=label)
    ax.set(xlabel=LABEL_POSITION_NM,
           ylabel=LABEL_FRICTION,
           xlim=[min(x), max(x)],
           )


def plot_dG(x, dG, ax, label=None, color=None):
    """Plot free energy vs position.

    Args:
        x (array-like):
            Positions along the coordinate (in nm).
        dG (array-like):
            Free energy values (in kJ/mol).
        ax (matplotlib.axes.Axes):
            The axes object where the plot will be drawn.
        label (str, optional):
            Label for the free energy curve. Default is None.
        color (str, optional):
            Color of the free energy curve. Default is None.

    Returns:
        matplotlib.lines.Line2D:
            The line object representing the plotted curve.

    Example:
        >>> fig, ax = plt.subplots()
        >>> plot_dG(positions, free_energy_values, ax)
        >>> plt.legend()
        >>> plt.show()
    """
    line = ax.plot(x, dG, label=label, color=color)[0]
    ax.set(xlabel=LABEL_POSITION_NM,
           ylabel=r'$\Delta G$ [kJ/mol]',
           xlim=[min(x), max(x)],
           )
    return line


def plot_dG_werrors(
    workestimator,
    ax,
    labeldG=None,
    color=None,
):
    """Plot free energy with errors against position.
    This function generates a plot of the free energy change (ΔG)
    as a function of position (x), including error bands if available.

    Args:
        workestimator:
            dcTMD workestimator object.
            It is expected to have the following attributes:
            `position_`, `dG_`, `s_dG_`
        ax (matplotlib.axes.Axes):
            The axes object where the plot will be drawn.
        labeldG (str, optional):
            Label for the free energy curve. Default is None.
        color (str, optional):
            Color of the free energy curve. Default is None.

    Notes:
        - If `s_dG_` is not available in `workestimator`,
            the function will print an errors message.
        - The error band is plotted as a shaded region
            around the free energy curve.

    Example:
        >>> fig, ax = plt.subplots()
        >>> plot_dG_werrors(workestimator=my_estimator, ax=ax)
        >>> plt.show()
    """
    if hasattr(workestimator, 's_dG_'):
        x = workestimator.position_
        dG = workestimator.dG_
        sdG = workestimator.s_dG_
        line = plot_dG(x, dG, ax, label=labeldG, color=color)
        color = line.get_color()
        if len(sdG) == 2:
            ax.fill_between(
                x,
                sdG[0],
                sdG[1],
                facecolor=color,
                alpha=ALPHA_VALUE,
            )
        else:
            ax.fill_between(
                x,
                dG - sdG,
                dG + sdG,
                facecolor=color,
                alpha=ALPHA_VALUE,
            )
    else:
        sys.stdout.write(
            f'no errors are determined for {workestimator}'
        )
        sys.stdout.write(
            'use estimate_free_energy_errors() to determine errors'
        )
        return


def plot_worklines(workset, ax, x=None, color='#777', res=1):
    """Line plots of the individual work trajectories
    in the workset.

    Args:
        workset:
            workset object
        ax (matplotlib.axes.Axes):
            The axes object where the plot will be drawn.
        x (array-like, optional):
            The x-axis positions for the plot. If not provided,
            `workset.position_` is used.
        color (str, optional):
            Color of the work lines. Default is '#777'.
        res (int, optional):
            Resolution for downsampling the data.
            Only every `res`-th point
            will be plotted. Default is 1 (no downsampling).

    Notes:
        - The x-axis represents the position (in nm).
        - The y-axis represents the work values (in kJ/mol).
        - Each trajectory is plotted as a semi-transparent line.

    Example:
        >>> fig, ax = plt.subplots()
        >>> plot_worklines(workset=my_workset, ax=ax)
        >>> plt.show()
    """
    if x is None:
        x = workset.position_
    for w in workset.work_:
        w_reduced = w[::res]
        x_reduced = x[::res]
        ax.plot(
            x_reduced,
            w_reduced,
            color=color,
            alpha=ALPHA_VALUE,
            lw=WORKLINES_WIDTH,
        )
    ax.set(xlabel=LABEL_POSITION_NM,
           ylabel=r'work $W$ [kJ/mol]',
           xlim=[min(x), max(x)],
           )


def plot_histo_normaldist(histodata, ax, color='None', label=None):
    """Plots a histogram of the input data and overlays a
    probability density function (PDF) of a normal distribution
    fitted to the data.

    Args:
        histodata (array-like):
            The data to be plotted.
            It will be flattened if not already 1D.
        ax (matplotlib.axes.Axes):
            The axes object where the plot will be drawn.
        color (str, optional):
            Color of the histogram and the fitted
            normal distribution curve.
            Default is 'None'.
        label (str, optional):
            Label for the histogram. Default is None.

    Notes:
        - Bin width is estimated with the Freedman-Diaconis rule (bins='fd').
        - The normal distribution is fitted using `scipy.stats.norm.fit`.
        - The x-axis represents the probability density (P),
        - The y-axis represents the work values (W) in kJ/mol.

    Example:
        >>> fig, ax = plt.subplots()
        >>> plot_histo_normaldist(my_data, ax=ax)
        >>> plt.legend()
        >>> plt.show()
    """
    histo = histodata.ravel()
    ax.hist(
        histo,
        bins='fd',
        density=True,
        histtype='stepfilled',
        align='mid',
        alpha=0.5,
        orientation='horizontal',
        color=color,
        label=label,
        ec=color,
        zorder=3,
    )
    mu, std = norm.fit(histo)
    # Plot the PDF.
    y = np.linspace(np.min(histo) - 10,
                    np.max(histo) + 10,
                    100,
                    )
    p = norm.pdf(y, mu, std)
    ax.plot(p, y, color=color, zorder=1)
    ax.set(xlabel=r'$P$',
           ylabel=r'$W$ [kJ/mol]',
           )
