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

## Sections:
- [**Theoretical Background:**](tutorials/theory.md): Here, you will learn the basic theory behind dcTMD. Including Jarzinskys equality, the derivation of the free energy and friction estimate as well as the main assumptions made on the way.

- [**Create pulling trajectories with Gromacs:**](tutorials/Gromacs.md): Here, you will learn how you can set up constraint targeted MD simulations using the pull code implemented in Gromacs. 

- **dcTMD Analysis:** In section [dcTMD via Work](tutorials/work.ipynb) and [dcTMD via Force](tutorials/force.ipynb) you will learn how to analyze the constraint pulling trajectories with dcTMD as described in [Theory](tutorials/theory.md).

- [**Command Line Interface:**](tutorials/CLI.ipynb) In this section, we will provide a short guide to the command line interface of `dcTMD`, which provides some common analysis and visualization functionality.