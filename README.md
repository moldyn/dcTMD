<div align="center">
    <img  style="width: 400px;" src="https://raw.githubusercontent.com/moldyn/dcTMD/main/docs/dctmd_logo.svg" />
  <p>
    <a href="https://github.com/wemake-services/wemake-python-styleguide" alt="wemake-python-styleguide">
        <img src="https://img.shields.io/badge/style-wemake-000000.svg" /></a>
    <a href="https://beartype.rtfd.io" alt="bear-ified">
        <img src="https://raw.githubusercontent.com/beartype/beartype-assets/main/badge/bear-ified.svg" /></a>
    <a href="https://moldyn.github.io/dcTMD" alt="Docs">
        <img src="https://img.shields.io/badge/mkdocs-Documentation-brightgreen" /></a>
    <a href="https://pypi.org/project/dcTMD" alt="PyPI">
        <img alt="PyPI" src="https://img.shields.io/pypi/v/dcTMD" /></a>
    <a href="https://anaconda.org/conda-forge/dcTMD" alt="conda version">
        <img src="https://img.shields.io/conda/vn/conda-forge/dcTMD" /></a>
    <a href="https://pepy.tech/project/dcTMD" alt="Downloads">
        <img src="https://pepy.tech/badge/dcTMD" /></a>
    <a href="https://github.com/moldyn/dcTMD/blob/main/LICENSE" alt="License">
        <img src="https://img.shields.io/github/license/moldyn/dcTMD" /></a>
    <a href="https://github.com/moldyn/dcTMD/actions/workflows/codeql.yml" alt="CodeQL">
        <img src="https://github.com/moldyn/dcTMD/actions/workflows/codeql.yml/badge.svg?branch=main" /></a>
    <a href="https://github.com/moldyn/dcTMD/actions/workflows/pytest.yml" alt="GitHub Workflow Status">
        <img src="https://img.shields.io/github/actions/workflow/status/moldyn/dcTMD/pytest.yml?branch=main"></a>
    <a href="https://codecov.io/gh/moldyn/dcTMD" >
        <img src="https://codecov.io/gh/moldyn/dcTMD/branch/main/graph/badge.svg?token=XMLP2VUU33"/></a>
  </p>

  <p>
    <a href="#features">Features</a> •
    <a href="#installation">Installation</a> •
    <a href="https://moldyn.github.io/dcTMD/getting_started/">Tutorials</a> •
    <a href="https://moldyn.github.io/dcTMD/">Docs</a>
  </p>
</div>


# dcTMD

This package aids in the analysis of dissipation-corrected targeted molecular dynamics (dcTMD) simulations. The method enforces rare unbinding events of ligands from proteins via a constraint pulling bias. Subsequently, free energy profiles and friction factors are estimated along the unbinding coordinate. For a methodological overview, see our [article](https://pubs.acs.org/doi/full/10.1021/acs.jctc.8b00835).

> S. Wolf, and G. Stock,  
> *Targeted molecular dynamics calculations of free energy profiles using a nonequilibrium friction correction.*,
> **J. Chem. Theory Comput.** 2018 14 (12), 6175-6182,  
> doi: [10.1021/acs.jctc.8b00835](https://pubs.acs.org/doi/10.1021/acs.jctc.8b00835)

This package will be published soon:

> V. Tänzel, M. Jäger, D. Nagel, and S. Wolf,  
> *Dissipation Corrected Targeted Molecular Dynamics*,  
> in preparation 2023

We kindly ask you to cite these articles in case you use this software package for published works.

## Features
- Intuitive usage via module and CI
- Sklearn-style API for fast integration into your Python workflow
- Supports Python 3.8-3.10
- Multitude of [publications](https://www.moldyn.uni-freiburg.de/publications.html) with dcTMD

## Implemented Key Functionalities
- Estimation of free energy profiles and friction factors along the unbinding coordinate of ligands as described by [Wolf and Stock 2018](https://pubs.acs.org/doi/full/10.1021/acs.jctc.8b00835).
- Analysis of separate unbinding pathways as described by [Wolf et al. 2022](https://arxiv.org/abs/2212.07154).

## Installation
The package will be available on PiPY and conda. Until then, install it via:
```bash
python3 -m pip install git+ssh://git@github.com/moldyn/dcTMD.git
```

## Usage
Check out the documentation for an overview over all modules as well as the tutorials.

## Roadmap

- [ ] New Features:
    - [ ] Gaussian error estimation
    - [ ] 2d distribution WorkSet plots
    - [x] Estimator plots: free energy, friction & both
    - [x] Normality plot
    - [x] Confidence intervals
    - [ ] Exponential estimator class
- [ ] Discuss gaussian kernel borders

