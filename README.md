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

This package aids in the analysis of dissipation-corrected targeted molecular dynamics (dcTMD) simulations. The method enforces rare unbinding events of ligands from proteins via a constraint pulling bias. Subsequently, free energy profiles and friction factors are estimated along the unbinding coordinate. For a methodological overview, see our [article](https://pubs.acs.org/doi/full/10.1021/acs.jctc.8b00835)

> S. Wolf, and G. Stock,
> *Targeted molecular dynamics calculations of free energy profiles using a nonequilibrium friction correction.*,
> Journal of chemical theory and computation (2018)

This package will be published soon:

> V. Tänzel, and M. Jäger, and S. Wolf,
> *Dissipation Corrected Targeted Molecular Dynamics*,
> in preparation (2022)

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
The package will be available on pipy and conda. Until then, install it via:
```bash
python3 -m pip install git+ssh://git@github.com/moldyn/dcTMD.git
```

## Usage
Check out the documentation for an overview over all modules as well as the tutorials.

## Roadmap

Next steps:
- [ ] Finish ForceSet, ForceEstimator.
- [ ] Rewrite IO of results (discuss this).
- [ ] Rewrite CLI and provide only simplest case: Estimation of free energy and friction from list of files plus saving.
- [ ] Smoothing of friction within Estimators via self.smooth_friction() methods. Better(?) even: Write abstract parent class which provides this smoothing method and inherit it.
- [ ] Tests for all functionality. So far, only the WorkSet has two tests. Add tests for the ForceSet, for the Estimators WorkEstimator and ForceEstimator as well as all their methods. Also provide tests for the io and for the utils, so smoothing and bootstrapping.
- [ ] Linting using [wemake](https://github.com/wemake-services/wemake-python-styleguide). Provide short explanation for ignores.
- [ ] Documentation: (1) Go through the MkDoc pages and check the automatically generated documentation based on the docstrings. Often, there are minor mistakes in content and formatting. (2) Create index pages, which are so far empty. (3) Link the tutorials in the navigation bar by adding notebook support.
- [ ] Tutorials: (1) For dcTMD analysis with this package, (2) for pathway separation following the work for trypsin.
- [ ] Set repo to public and make sure the documentation page works, as well as the tests. Enforce high test coverage.

Secondary importance:
- [ ] New Features: 
    - [ ] RAM estimator
    - [ ] Gaussian error estimation
    - [ ] WorkSet plots: 1d, 2d distributions
    - [ ] Estimator plots: free energy, friction & both
    - [ ] Normality plot
    - [ ] Confidence intervals
    - [ ] Exponential estimator class
- [ ] Check if smoothing width works as expected
- [ ] Discuss gaussian kernel borders

Finishing steps:
- [ ] Design logo.
- [ ] Publish to pypi & conda.
- [ ] Publish publication note.
