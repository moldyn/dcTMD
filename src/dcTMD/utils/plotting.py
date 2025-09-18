# -*- coding: utf-8 -*-
# MIT License
# Copyright (c) 2021-2023, Miriam Jäger
# All rights reserved.
"""Simple plot functions for dcTMD results."""


import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import probplot, norm


"""
def fig_sizeA4width():
    # Convert cm to inches
    # A4 width
    fig_width_cm = 21
    fig_height_cm = 14.8 / 2
    inches_per_cm = 1 / 2.54
    # width in inches
    fig_width = fig_width_cm * inches_per_cm
    # height in inches
    fig_height = fig_height_cm * inches_per_cm
    fig_size = (fig_width, fig_height)
    return fig_size


def fig_sizehalfA4width():
    # Convert cm to inches
    fig_width_cm = 21 / 2
    fig_height_cm = 14.8 / 2
    inches_per_cm = 1 / 2.54
    # width in inches
    fig_width = fig_width_cm * inches_per_cm
    # height in inches
    fig_height = fig_height_cm * inches_per_cm
    fig_size = (fig_width, fig_height)
    return fig_size
"""


def plot_dcTMD_results(
    estimator,
    friction=None,
    x=None,
    figsize=(4, 4),
):
    """Plot dcTMD results overview in two subplots.

    Top subplot constains free energy, dissipative work and mean work.
    Bottom subplot contains friction vs. position.
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
    """Plot free energy, dissipative work and mean work vs position."""
    if x is None:
        x = workestimator.position_
    ax.plot(x, workestimator.dG_, label=r'$\Delta G$')
    ax.plot(x, workestimator.W_mean_, label=r'W$_{\mathrm{mean}}$')
    ax.plot(x, workestimator.W_diss_, label=r'W$_{\mathrm{diss}}$')
    ax.set(xlabel=r'$x$ [nm]',
           ylabel=r'[kJ/mol]',
           xlim=[min(x), max(x)],
           )


def plot_Gamma(x, friction, ax, label=None):
    """Plot friction factor vs position."""
    ax.plot(x, friction, label=label)
    ax.set(xlabel=r'$x$ [nm]',
           ylabel=r'$\Gamma$ [kJ nm$^2$/(mol ps)]',
           xlim=[min(x), max(x)],
           )


def plot_dG(x, dG, ax, label=None):
    """Plot free energy vs position."""
    line, = ax.plot(x, dG, label=label)
    ax.set(xlabel=r'$x$ [nm]',
           ylabel=r'$\Delta G$ [kJ/mol]',
           xlim=[min(x), max(x)],
           )
    return line


def plot_dG_werrors(workestimator, ax, labeldG=None):
    """Plot free energy with errors against position."""
    if hasattr(workestimator, 's_dG_'):
        x = workestimator.position_
        dG = workestimator.dG_
        sdG = workestimator.s_dG_
        line = plot_dG(x, dG, ax, label=labeldG)
        color = line.get_color()
        if len(sdG) == 2:
            ax.fill_between(
                x,
                sdG[0],
                sdG[1],
                facecolor=color,
                alpha=0.3,
            )
        else:
            ax.fill_between(
                x,
                dG - sdG,
                dG + sdG,
                facecolor=color,
                alpha=0.3,
            )
    else:
        print(f'no errors are determined for {workestimator}')
        print('use estimate_free_energy_errors() to determine errors')
        return


def plot_worklines(workset, ax, x=None, color='#777', res=1):
    """Line plots of work of the individual trajectories
    in the workset."""
    if x is None:
        x = workset.position_
    for w in workset.work_:
        ax.plot(x[::res], w[::res], color=color, alpha=.3, lw=.5)

    ax.set(xlabel=r'$x$ [nm]',
           ylabel=r'work $W$ [kJ/mol]',
           xlim=[min(x), max(x)],
           )


def plot_histo_normaldist(data, ax, color='blue', label=None):
    """Plots a histogram of the data and
    a normal distribution fitted to the data."""
    data = data.flatten()
    ax.hist(data,
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
    mu, std = norm.fit(data)
    # Plot the PDF.
    y = np.linspace(np.min(data) - 10,
                    np.max(data) + 10,
                    100
                    )
    p = norm.pdf(y, mu, std)
    ax.plot(p, y, color=color, zorder=1)
    ax.set(xlabel=r'$P$',
           ylabel=r'$W$ [kJ/mol]',
           )


def plot_worknormalitychecks(
    workset,
    index,
    x=None,
    colors=None,
    figsize=(6, 2),
    axs=None,
    res=10,
):
    """Plots the work values of trajectories individually.

    Also adds histograms and normality plots for the indices given in `index`.
    """
    if axs is None:
        print('No axs given. Create figure.')
        fig, axs = plt.subplots(
            ncols=3,
            nrows=1,
            figsize=figsize,
        )
    if x is None:
        x = workset.position_
    plot_worklines(workset, axs[0], x=x, res=res)

    if not colors:
        cmap = plt.get_cmap('Dark2')
        colors = cmap.colors

    for j, idx in enumerate(index):
        work = workset.work_[:, idx].flatten()
        axs[1].set_title(r'Histogram at $x$')
        plot_histo_normaldist(work, axs[1], colors[j])
        axs[0].axvline(x[idx],
                       color=colors[j],
                       zorder=3,
                       label=rf'$x={x[idx]:.2f}$',
                       )

        probplot(work, plot=axs[2], fit=True)
        axs[2].get_lines()[j * 2].set_color(colors[j])
        axs[2].set_title('Normality plot')

    axs[0].legend()
    return axs
