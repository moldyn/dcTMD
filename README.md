<div align="center">
  <p>
    <a href="https://github.com/wemake-services/wemake-python-styleguide" alt="wemake-python-styleguide">
        <img src="https://img.shields.io/badge/style-wemake-000000.svg" /></a>
    <a href="https://moldyn.github.io/dcTMD" alt="Docs">
        <img src="https://img.shields.io/badge/mkdocs-Documentation-brightgreen" /></a>
    <a href="https://github.com/moldyn/dcTMD/blob/main/LICENSE" alt="License">
        <img src="https://img.shields.io/github/license/moldyn/dcTMD" /></a>
  </p>

  <p>
    <a href="#features">Features</a> •
    <a href="#installation">Installation</a> •
    <a href="#usage">Usage</a>
  </p>
</div>


# dcTMD

> **Warning**
> This package is still in beta stage. Please open an issue if you encounter
> any bug/error.

This package aids in the analysis of dissipation-corrected targeted molecular dynamics (dcTMD) simulations. The simulaiton method enforces rare unbinding events of ligands from proteins via a constraint pulling bias. Subsequently, free energy profiles and friction factors are estimated along the unbinding coordinate. For a methological overview, see our [article](https://pubs.acs.org/doi/full/10.1021/acs.jctc.8b00835)

> S. Wolf, and G. Stock,
> *Targeted molecular dynamics calculations of free energy profiles using a nonequilibrium friction correction.*,
> Journal of chemical theory and computation (2018)

This package will be published soon:

> V. Tänzel, and M. Jäger, and S. Wolf,
> *Dissipation Corrected Targeted Molecular Dynamics
> in preparation

We kindly ask you to cite these articles in case you use this software package for published works.

## Features
- Intuitive usage via module and via CI
- Sklearn-style API for fast integration into your Python workflow
- Supports Python 3.6-3.10
- Multitude of [publications](https://www.moldyn.uni-freiburg.de/publications.html) with dcTMD

## Implemented Key Functionalities
- Estimation of free energy profiles and friction factors along the unbinding coordinate of ligands as decribed by [Wolf and Stock 2018](https://pubs.acs.org/doi/full/10.1021/acs.jctc.8b00835).
- Analysis of separate unbinding pathways as decribed by [Wolf et al. 2018](https://arxiv.org/abs/2212.07154).

## Installation
The package will be available on pipy and conda. Until then, install it via:
```bash
python3 -m pip install git+ssh://git@github.com/moldyn/dcTMD.git
```

## Usage
Check out the documentation for an overview over all modules as well as the tutorials.

## Roadmap

Tier 1:
- Finish ForceSet, ForceEstimator
- Rewrite IO of results.
- Save & load Estimator instances.
- Rewrite CLI.
- Smoothing of friction within Estimators.

- Tests
- Linting

- Documentation
- Tutorial
- Pathway separation tutorial
- Github page

Tier 2:
- RAM estimator schreiben
- Set up package
- Gaussian error estimation
- WorkSet plots: 1d, 2d distributions
- Estimator plots: free energy, friction & both
- Confidence intervals
- Check if smoothing width works as expected
- Discuss gaussian kernel borders
- Normality plot
- Exponential estimator

Tier 3:
- Design logo.
- Publish to pipy & conda.
- Publish publication note.
