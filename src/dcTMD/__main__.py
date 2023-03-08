#!/usr/bin/env python3
"""
CLI of dcTMD.

MIT License
Copyright (c) 2021, Victor Tänzel, Miriam Jaeger, Steffen Wolf.
All rights reserved.
"""

import click
import dcTMD

MODES = ('work', 'force')


@click.command(no_args_is_help=True)
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
    required=True,
    help='Input: File containing list of all constraint force file names' +
    'or glob pattern e.g."*.xvg" to generate a list of all constraint ' +
    'force files using glob.glob()',
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
    '--resamples',
    type=int,
    required=False,
    help='Number of resamples used in optional bootstrapping. ' +
    'This is only available in mode work',
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
    help='Save the Work/ForceSet class to file.',
)
def main(
    mode,
    pullf_files,
    outname,
    temperature,
    velocity,
    res,
    sigma,
    resamples,
    verbose,
    plot,
    save_dataset,
) -> None:
    r"""
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
    # ToDo: discuss if we want to include error estimation
    if verbose:
        click.echo(f'Input:\n mode: {mode}\n file: {pullf_files}\n \
        outname: {outname}\n temperature: {temperature}\n \
        velocity: {velocity}\n resolution: {res}\n sigma: {sigma}\n verbose: \
        {verbose}\n, plot: {plot}\n, N_resamples: {resamples}\n\
        save dataset: {save_dataset}\n')

    # Set up mode
    if mode == 'work':
        dataset = dcTMD.storing.WorkSet(velocity=velocity,
                                        resolution=res,
                                        verbose=verbose,
                                        )
        estimator = dcTMD.dcTMD.WorkEstimator(temperature)
    if mode == 'force':
        dataset = dcTMD.storing.ForceSet(velocity=velocity,
                                         resolution=res,
                                         verbose=verbose,
                                         )
        estimator = dcTMD.dcTMD.ForceEstimator(temperature)

    # Loading constraint force files
    filenames = dcTMD.io.load_pullf(pullf_files)

    # Generate work/force set
    dataset.fit(filenames)
    if save_dataset:
        out = outname + f"_{len(dataset.names_)}_{mode}set"
        dcTMD.storing.save(out, dataset)
    # Calculate Wmean, Wdiss, dG and friction factor
    estimator.fit(dataset)

    """
    # calculate errors
    if resamples:
        estimator.estimate_free_energy_errors(N_resamples)
    """

    # save data as .npz and .dat file
    outname += f'_{mode}'
    dcTMD.io.write_output(outname, dataset, estimator)


if __name__ == '__main__':
    main()
