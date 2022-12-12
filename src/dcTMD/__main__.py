#!/usr/bin/env python3
"""
CLI of dcTMD.

MIT License
Copyright (c) 2021, Victor Tänzel, Miriam Jaeger, Steffen Wolf.
All rights reserved.
"""

import click
import dcTMD

import numpy as np

from dcTMD._clickextensions import Mutex

MODES = ('work', 'force')


@click.command(no_args_is_help=True)
@click.option(
    '-dG',
    '--calc_dG',
    'calc_dG',
    is_flag=True,
    type=bool,
    default=False,
    help='Set flag to calculate dG',
)
@click.option(
    '-friction',
    '--calc_friction',
    'calc_friction',
    is_flag=True,
    type=bool,
    default=False,
    help='Set flag to calculate friction',
)
@click.option(
    '-m',
    '--mode',
    type=click.Choice(MODES, case_sensitive=True),
    default=MODES[0],
    show_default=True,
    required=True,
    help='Use either work or force autocovariance function to calculate' +
    'dcTMD quantities.',
)
@click.option(
    '-f',
    '--file',
    'pullf_files',
    help='Input: File containing list of all constraint force file names',
    cls=Mutex,
    not_required_if=['pullf_glob_pattern'],
)
@click.option(
    '-g',
    '--glob',
    'pullf_glob_pattern',
    help='Input: Glob pattern generating a list of all constraint force file \
names',
    cls=Mutex,
    not_required_if=['pullf_files'],
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
    'T',
    type=float,
    required=True,
    help='Simulation temperature in K.',
)
@click.option(
    '-vel',
    '--velocity',
    'vel',
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
    help='Striding to reduce size of returned free energy and friction',
)
@click.option(
    '-s',
    '--sigma',
    type=float,
    required=False,
    help='Standard deviation of gaussian filter in nm.',
)
@click.option(
    '-N',
    '--N_resamples',
    type=float,
    required=False,
    help='Number of resamples used in optional bootstrapping.',
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
    'p',
    '--plot',
    is_flag=True,
    default=False,
    show_default=True,
    help='Plots free energy and smoothed friction.',
)
def main(
    calc_dG,
    calc_friction,
    mode,
    pullf_files,
    pullf_glob_pattern,
    outname,
    temp,
    vel,
    res,
    sigma,
    N_resamples,
    verbose,
    plot,
) -> None:
    r"""
    \b
    -------------------------
    |         dcTMD         |
    -------------------------
    Calculate free energy and friction for given constraint force files.

    Analysis tools for dissipation-corrected targeted molecular dynamics, which
    is an enhanced sampling method to enforce rare events in biomolecular
    systems. When publishing results gained with this python package, please
    cite the following publications:
    (1) Tänzel, Victor and Jäger, Miriam and Wolf, Steffen in preparation.
    (2) Wolf, Steffen, and Gerhard Stock. "Targeted molecular dynamics
    calculations of free energy profiles using a nonequilibrium friction
    correction." Journal of chemical theory and computation 14.12 (2018): 6175-
    6182.
    """
    # Click aftercare
    if verbose:
        click.echo(f'Input:\n calc_dG: {calc_dG};\n calc_friction: \
{calc_friction};\n mode: {mode}\n file: {pullf_files}\n glob: \
{pullf_glob_pattern}\n outname: {outname}\n temperature: {temp}\n \
velocity: {vel}\n resolution: {res}\n sigma: {sigma}\n verbose: \
{verbose}\n, plot: {plot}\n, N_resamples: {N_resamples}\n')

    if not calc_dG and not calc_friction:
        print("No '-dG' or '-friction' flag found, please specify what to do.")
        exit()

    if not (pullf_files or pullf_glob_pattern):
        print("Please provide constraint force files via '-f' or '-g'.")
        exit()

    # Set up mode
    if mode == 'work':
        calc_dG = dcTMD.work.calc_dG
        calc_dG_and_friction = dcTMD.work.calc_dG_and_friction
        pullf_to_dataset = dcTMD.work.pullf_to_work_array
        outname += '_from_work'
    if mode == 'force':
        calc_dG = dcTMD.force.calc_dG
        calc_dG_and_friction = dcTMD.force.calc_dG_and_friction
        pullf_to_dataset = dcTMD.force.pullf_to_force_array
        outname += '_from_forceacf'

    # Loading constraint force files
    files = dcTMD.io._load_pullf(pullf_glob_pattern, pullf_files)

    # Generate work/force set
    dataset, time, filenames = pullf_to_dataset(files, vel, verbose, res)

    pos = time * vel
    N_traj, _ = np.shape(dataset)
    timestep = time[1] - time[0]

    if mode == 'work':
        if calc_friction:
            print('calculate free energy and friction')
            Wmean, Wdiss, dG, Gamma, Gamma_smooth, *args = \
                calc_dG_and_friction(
                    dataset,
                    temp,
                    vel,
                    timestep,
                    sigma,
                    N_resamples,
                )
            outname += f'_sig{sigma}nm'
        else:
            print('calculate free energy')
            Wmean, Wdiss, dG, *args = calc_dG(
                dataset,
                temp,
                N_resamples,
            )

    if mode == 'force':
        if calc_friction:
            print('calculate free energy and friction')
            Wmean, Wdiss, dG, Gamma, Gamma_smooth, *args = \
                calc_dG_and_friction(dataset, temp, time, vel, sigma, res)

            outname += f'_sig{sigma}nm'
            print(f'saving output {outname}')
        else:
            print('calculate free energy')
            Wmean, Wdiss, dG, *args = calc_dG(dataset, temp, time, vel, res)
            # striding
        pos = pos[::res]

    if calc_friction:
        dcTMD.io.write_output(
            outname,
            N_traj,
            pos=pos,
            Wmean=Wmean,
            Wdiss=Wdiss,
            dG=dG,
            Gamma=Gamma,
            Gamma_smooth=Gamma_smooth,
            errors=args,
        )
    else:
        dcTMD.io.write_output(
            outname,
            N_traj,
            pos=pos,
            Wmean=Wmean,
            Wdiss=Wdiss,
            dG=dG,
            errors=args,
        )


if __name__ == '__main__':
    main()
