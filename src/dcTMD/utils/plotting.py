# -*- coding: utf-8 -*-
"""
Simple plot functions for dcTMD results.

MIT License
Copyright (c) 2021-2022, Miriam Jäger
All rights reserved.
"""


import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import probplot, norm


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


def plot_dcTMD_results(x, workestimator, friction):
    """Plot dcTMD results in two subplots."""
    fig, axs = plt.subplots(ncols=1,
                            nrows=2,
                            sharex=True,
                            figsize=fig_sizehalfA4width(),
                            )
    plot_dG_Wdiss(x, workestimator, axs[0])
    plot_Gamma(x, friction, axs[1])
    axs[0].legend(loc='lower left', mode='expand',
                  bbox_to_anchor=(0, 0.9, 1, 0.2),
                  frameon=False,
                  ncol=3,)
    axs[0].set_xlabel("")
    plt.tight_layout()
    return fig, axs


def plot_dG_Wdiss(x, workestimator, ax):
    """Plot free energy, dissipative work and mean work vs x."""
    ax.plot(x, workestimator.dG_, label=r'$\Delta G$')
    ax.plot(x, workestimator.W_mean_, label=r'W$_{\mathrm{mean}}$')
    ax.plot(x, workestimator.W_diss_, label=r'W$_{\mathrm{diss}}$')
    ax.set(xlabel=r'position $x$ [nm]',
           ylabel=r'[kJ/mol]',
           xlim=[min(x), max(x)],
           )


def plot_Gamma(x, friction, ax, label=None):
    """Plot friction vs x."""
    ax.plot(x, friction, label=rf"{label}")
    ax.set(xlabel=r'position $x$ [nm]',
           ylabel=r'$\Gamma$ [kJ/mol/(nm$^2$/ps)]',
           xlim=[min(x), max(x)],
           )


def plot_dG(x, dG, ax, label=None):
    """Plot dG vs x."""
    ax.plot(x, dG, label=rf"{label}")
    ax.set(xlabel=r'position $x$ [nm]',
           ylabel=r'$\Delta G$ [kJ/mol]',
           xlim=[min(x), max(x)],
           )


def plot_worklines(x, workset, ax):
    """Plot work values of trajectories individually."""
    for w in workset:
        ax.plot(x, w, color='#777', alpha=.5, lw=.5)

    ax.set(xlabel=r'position $x$ [nm]',
           ylabel=r'work $W$ [kJ/mol]',
           xlim=[min(x), max(x)],
           )


def plot_histo_normaldist(data, ax, color='blue', label=None):
    data = data.flatten()
    ax.hist(data,
            bins='fd',
            density=True,
            histtype='stepfilled',
            align='mid',
            alpha=0.5,
            orientation='horizontal',
            color=color,
            label=rf'{label}',
            ec=color,
            zorder=3,
            )
    mu, std = norm.fit(data)
    # Plot the PDF.
    y = np.linspace(np.min(data)-10,
                    np.max(data)+10,
                    100
                    )
    p = norm.pdf(y, mu, std)
    ax.plot(p, y, color=color, zorder=1)
    ax.set(xlabel=r'$P$',
           ylabel=r'$W$ [kJ/mol]',
           )


def plot_worknormalitychecks(x, workset, index, colors=None):
    """Plots the work values of trajectories individually.
    And adds histograms and normality plots for the indices given in index."""
    fig, axs = plt.subplots(ncols=3,
                            nrows=1,
                            figsize=fig_sizeA4width()
                            )
    plot_worklines(x, workset, axs[0])

    if not colors:
        from matplotlib.cm import get_cmap
        cmap = get_cmap('Dark2')
        colors = cmap.colors

    for j, idx in enumerate(index):
        data = workset[:, idx].flatten()
        axs[1].set_title('Histrogram at x')
        plot_histo_normaldist(data, axs[1], colors[j])
        axs[0].axvline(x[idx],
                       color=colors[j],
                       zorder=3,
                       label=rf'$x={x[idx]}$nm',
                       )

        probplot(data, plot=axs[2], fit=True)
        axs[2].get_lines()[j*2].set_color(colors[j])

    axs[0].legend()
    plt.tight_layout()

    return fig, axs
