# -*- coding: utf-8 -*-
# MIT License
# Copyright (c) 2021-2023, Victor Tänzel, Miriam Jaeger, Steffen Wolf.
# All rights reserved.
"""CLI of dcTMD."""

import click
from dcTMD.storing import WorkSet, ForceSet
from dcTMD.dcTMD import WorkEstimator, ForceEstimator
from dcTMD.io import load_pullf, write_output
from dcTMD.storing import save
from dcTMD.utils import plotting
import matplotlib.pyplot as plt

MODES = ('work', 'force')


@click.command(no_args_is_help=True)
@click.option(
    '-m',
    '--mode',
    type=click.Choice(MODES, case_sensitive=True),
    default=MODES[0],
    show_default=True,
    required=True,
    help='Use either work or force autocovariance function to calculate '
    + 'dcTMD quantities.',
)
@click.option(
    '-f',
    '--file',
    'pullf_files',
    required=True,
    help='Input: File containing list of all constraint force file names '
    + 'or glob pattern e.g."*.xvg" to generate a list of all constraint '
    + 'force files using glob.glob()',
)
@click.option(
    '-o',
    '--outname',
    type=click.Path(),
    default='./',
    help='Output: Path/prefix of output names.',
)
@click.option(
    '-T',
    '--temperature',
    type=float,
    required=True,
    help='Simulation temperature in K.',
)
@click.option(
    '-vel',
    '--velocity',
    type=float,
    required=True,
    help='Pulling velocity in nm/ps.',
)
@click.option(
    '--res',
    type=int,
    default=1,
    show_default=True,
    required=False,
    help='Striding to reduce size of returned free energy and friction.',
)
@click.option(
    '-s',
    '--sigma',
    type=float,
    default=None,
    required=False,
    help='Standard deviation of gaussian filter in nm.',
)
@click.option(
    '-v',
    '--verbose',
    is_flag=True,
    default=False,
    show_default=True,
    help='Enable verbose mode.',
)
@click.option(
    '-p',
    '--plot',
    is_flag=True,
    default=False,
    show_default=True,
    help='Plots free energy and smoothed friction.',
)
@click.option(
    '-sd',
    '--save_dataset',
    is_flag=True,
    default=False,
    show_default=True,
    help='Save the Work/ForceSet instance to file.',
)
def main(  # noqa: WPS211, WPS216
    mode,
    pullf_files,
    outname,
    temperature,
    velocity,
    res,
    sigma,
    verbose,
    plot,
    save_dataset,
) -> None:
    """Calculate free energy and friction for given constraint force files.

    \b
    -------------------------
    |         dcTMD         |
    -------------------------

    Analysis tools for dissipation-corrected targeted molecular dynamics, which
    is an enhanced sampling method to enforce rare events in biomolecular
    systems. When publishing results gained with this python package, please
    cite the following publications:
    (1) Tänzel, Victor and Jäger, Miriam and Wolf, Steffen in preparation.
    (2) Wolf, Steffen, and Gerhard Stock. 'Targeted molecular dynamics
    calculations of free energy profiles using a nonequilibrium friction
    correction.' Journal of chemical theory and computation 14.12 (2018): 6175-
    6182.
    """
    # Click aftercare
    if verbose:
        click.echo(
            'Input:\n '
            f'mode: {mode}\n '
            f'file: {pullf_files}\n '
            f'outname: {outname}\n '
            f'temperature: {temperature}\n '
            f'velocity: {velocity}\n '
            f'resolution: {res}\n '
            f'sigma: {sigma}\n '
            f'verbose: {verbose}\n '
            f'plot: {plot}\n '
            f'save dataset: {save_dataset}\n'
        )

    # Set up mode
    if mode == 'work':
        dataset = WorkSet(
            velocity=velocity,
            resolution=res,
            verbose=verbose,
        )
        estimator = WorkEstimator(temperature)
    elif mode == 'force':
        dataset = ForceSet(
            velocity=velocity,
            resolution=res,
            verbose=verbose,
        )
        estimator = ForceEstimator(temperature)

    # Loading constraint force files
    filenames = load_pullf(pullf_files)

    # Generate work/force set
    dataset.fit(filenames)
    if save_dataset:
        out = f'{outname}_{mode}set'
        save(out, dataset)
    # Calculate Wmean, Wdiss, dG and friction factor
    estimator.fit(dataset)

    # Smooth friction
    if sigma:
        estimator.smooth_friction(sigma=0.1, mode='reflect')

    # save data as .npz and .dat file
    outname = f'{outname}_{mode}'
    write_output(outname, estimator)

    if plot:
        plotting.plot_dcTMD_results(
            estimator,
        )
        plt.savefig(f'{outname}.png')


if __name__ == '__main__':
    main()  # pragma: no cover
