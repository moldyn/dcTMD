## Introduction
This python package aids with the analysis of targeted molecular dynamics (TMD) trajectories according to   dissipation-corrected TMD method. TMD simulations enforce rare unbinding events of ligands from proteins via a constraint pulling bias. With dcTMD, free energy profiles and friction factors are estimated along the unbinding coordinate. For a methodological overview, see our [article](https://pubs.acs.org/doi/full/10.1021/acs.jctc.8b00835)


## Disclaimer
> **Warning**
> This package is still in beta stage. Please open an issue if you encounter
> any bug/error.

> S. Wolf, and G. Stock,
> *Targeted molecular dynamics calculations of free energy profiles using a nonequilibrium friction correction.*,
> Journal of chemical theory and computation (2018)

This package will be published soon:

> V. Tänzel, and M. Jäger, and S. Wolf,
> *Dissipation Corrected Targeted Molecular Dynamics*,
> in preparation (2022)

We kindly ask you to cite these articles in case you use this software package for published works.


## Installation
The package will be available on pipy and conda. Until then, install it via:
```bash
python3 -m pip install git+ssh://git@github.com/moldyn/dcTMD.git
```



## Classes and Modules
The module is structured into the following submodules:

    io: This submodule provides functions for handling input and output operations in dcTMD package. Specifically, the two functions: load_pullf and write_output.

    storing: This modules creates force or work sets from .pullf files which are needed for further analysis.
    This submodule provides two classes, WorkSet and ForceSet, that store constraint force data as work or force time traces, respectively.

    dcTMD: 

    utils: bootstrapping, plotting and smoothing

### Methods:

    __init__(self, temperature: (Float, Int), verbose: bool = False): Constructor for the class, which takes the temperature of the simulation and a verbose flag as input parameters.
    fit(self, work_set): Method to estimate the free energy and friction for a given work set. It takes a work set as input and returns a fitted estimator.
    transform(self, X, y=None): Method to return the free energy and friction estimates as a tuple.
    estimate_free_energy(self, work_set=None): Method to estimate the free energy for a given work set. It takes an optional work set as input and returns the mean work, dissipative work, and free energy estimate as a tuple.
    estimate_free_energy_errors(self, n_resamples: Int, mode: Union[StrStd, NumInRange0to1], seed=None): Method to estimate the errors in the free energy estimates using bootstrapping. It takes the number of resamples, the mode of calculation, and the random seed as input parameters and returns the resampled mean work, resampled dissipative work, and resampled free energy estimate as a tuple.

Attributes:

    temperature: The temperature of the simulation.
    verbose: The verbose flag, which determines whether to print verbose output.
    W_mean_: The mean work, in kJ/mol.
    W_diss_: The dissipative work, in kJ/mol.
    dG_: The free energy estimate, in kJ/mol.
    friction_: The friction factor in kJ/mol/(nm^2/ps).
    mode_: The parameter used to estimate free energy errors.
    s_W_mean_: The bootstrapping error of the mean work.
    s_W_diss_: The bootstrapping error of the dissipative work.
    s_dG_: The bootstrapping error of the free energy estimate.
    W_mean_resampled_: The resampled mean work, needed to inspect its distribution.
    W_diss_resampled_: The resampled dissipative work, needed to inspect its distribution.
    dG_resampled_: The resampled free energy estimate, needed to inspect its distribution.

This code defines a class ForceEstimator for performing dcTMD (discrete coordinate thermodynamics) analysis on a force set. The class has several methods and attributes, including fit, transform, and estimate_free_energy_friction.

The fit method takes in a force_set object and estimates the free energy and friction. The transform method returns the estimated free energy and friction. The estimate_free_energy_friction method estimates the free energy and friction using force auto-correlation.

The memory_kernel method calculates the memory kernel at positions x_indices "forward" in time from fluctuation-dissipation. It returns a 2D array corr_set with shape (len(x_indices), length_data).

Overall, this class seems to be part of a larger codebase for analyzing molecular dynamics simulations.