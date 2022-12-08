#!/usr/bin/env python3
"""
CLI of dcTMD.

MIT License
Copyright (c) 2021, Victor Taenzel, Miriam Jaeger
All rights reserved.

"""

import click
import dcTMD

import numpy as np
import matplotlib.pyplot as plt
from .utils import plot_dG, plot_Gamma_twinax
from ._clickextensions import Mutex

MODES = ['work', 'force']


@click.command(
    no_args_is_help=True
)
@click.option('-dG',
              '--calc_dG',
              'calc_dG',
              is_flag=True,
              type=bool,
              default=False,
              help='Set flag to calculate dG',
              )
@click.option('-friction',
              '--calc_friction',
              'calc_friction',
              is_flag=True,
              type=bool,
              default=False,
              help='Set flag to calculate friction',
              )
@click.option('-m',
              '--mode',
              type=click.Choice(MODES, case_sensitive=True),
              default=MODES[0],
              show_default=True,
              required=True,
              help='Use either work or force autocovariance function to \
calculate dG/friction.',
              )
@click.option('-f',
              '--file',
              'pullf_files',
              help='Input: File containing list of all constraint force \
file names',
              cls=Mutex,
              not_required_if=['pullf_glob_pattern']
              )
@click.option('-g',
              '--glob',
              'pullf_glob_pattern',
              help='Input: Glob pattern generating a list of all constraint \
force file names',
              cls=Mutex,
              not_required_if=['pullf_files']
              )
@click.option('-o',
              '--outname',
              type=click.Path(),
              default='./',
              help='Output: Path/prefix of output names.',
              )
@click.option('-T',
              '--temperature',
              'T',
              type=float,
              required=True,
              help='Simulation temperature in K.'
              )
@click.option('-vel',
              '--velocity',
              'vel',
              type=float,
              required=True,
              help='Pulling velocity in nm/ps.',
              )
@click.option('--res',
              type=int,
              default=1,
              show_default=True,
              required=False,
              help='Striding to reduce size of returned free energy and \
friction',
              )
@click.option('-s',
              '--sigma',
              type=float,
              required=False,
              help='Standard deviation of gaussian filter in nm.',
              )
@click.option('-v',
              '--verbose',
              is_flag=True,
              default=False,
              show_default=True,
              help='Enable verbose mode.'
              )
@click.option('p',
              '--plot',
              is_flag=True,
              default=False,
              show_default=True,
              help='Plots free energy and smoothed friction.'
              )
def main(calc_dG,
         calc_friction,
         mode,
         pullf_files,
         pullf_glob_pattern,
         outname,
         T,
         vel,
         res,
         sigma,
         verbose,
         plot,
         ) -> None:
    """
    \b
    -------------------------
    |         dcTMD         |
    -------------------------
    TODO: lizenz and paper
    Calculates dG and friction for given constraint force files.\\
        use pip -e to install package
    """
    # Click aftercare
    if verbose:
        click.echo(f'Input:\n calc_dG: {calc_dG};\n calc_friction: \
{calc_friction};\n mode: {mode}\n file: {pullf_files}\n glob: \
{pullf_glob_pattern}\n outname: {outname}\n temperature: {T}\n \
velocity: {vel}\n resolution: {res}\n sigma: {sigma}\n verbose: \
{verbose}\n, plot: {plot}\n')
    if not calc_dG and not calc_friction:
        print("No '-dG' or '-friction' flag found, please specify what to do.")
        exit()

    if not (pullf_files or pullf_glob_pattern):
        print("Please provide constraint force files via '-f' or '-g'.")
        exit()

    # Set up mode
    if mode == 'work':
        pass
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
    files = dcTMD.io.load_pullf(pullf_glob_pattern, pullf_files)

    # Generate work/force set
    dataset, t, filenames = pullf_to_dataset(files, vel, verbose, res)

    x = t * vel
    N, length_data = np.shape(dataset)
    timestep = t[1] - t[0]
    errors = False

    if mode == 'work':
        if calc_friction:
            print('calculate free energy and friction')
            Wmean, Wdiss, dG, Gamma, Gamma_smooth, *args = \
                calc_dG_and_friction(dataset, T, vel, timestep,
                                     sigma, errors, 10000)
            outname += '_sig{}nm'.format(sigma)
        else:
            print('calculate free energy')
            Wmean, Wdiss, dG, *args = calc_dG(dataset, T, errors, 10000)

    if mode == 'force':
        if calc_friction:
            print('calculate free energy and friction')
            Wmean, Wdiss, dG, Gamma, Gamma_smooth, *args = \
                calc_dG_and_friction(dataset, T, t, vel, sigma, res)

            outname += '_sig{}nm'.format(sigma)
            print("saving output {}".format(outname))
        else:
            print('calculate free energy')
            Wmean, Wdiss, dG, *args = calc_dG(dataset, T, t, vel, res)
            # striding
        x = x[::res]

    if not calc_friction:
        dcTMD.io.write_output(outname, N, x=x, Wmean=Wmean,
                              Wdiss=Wdiss, dG=dG, errors=args)
    else:
        dcTMD.io.write_output(outname, N, x=x, Wmean=Wmean,
                              Wdiss=Wdiss, dG=dG, Gamma=Gamma,
                              Gamma_smooth=Gamma_smooth, errors=args)

    if plot:
        fig, ax1 = plt.subplots(figsize=(5, 2.2))
        plot_dG(ax1, x, dG)
        if calc_friction:
            plot_Gamma_twinax(ax1, x, Gamma_smooth)

        plt.tight_layout()
        plt.savefig(outname + '.pdf')


if __name__ == '__main__':
    main()
