---
title: "dcTMD: a python package for performing dissipation-corrected targeted molecular dynamics"
authors:
  - name: Miriam Jäger
    affiliation: 1
  - name: Victor Tänzel
    affiliation: 1
  - name: Daniel Nagel
    affiliation: 2
  - name: Steffen Wolf
    affiliation: 1
affiliations:
  - index: 1
    name: University of Freiburg, Institute of Physics, Freiburg, Germany
  - index: 2
    name: University of Heidelberg, Institute of Theoretical Physics, Heidelberg, Germany
date: "2025-11-27"
bibliography: paper.bib
---


# Summary
`dcTMD` is a Python package designed to extract free-energy and nonequilibrium friction estimates from targeted molecular dynamics (TMD) simulations [@Schlitter1994_TMD]. The method implemented here is called *dissipation-corrected targeted molecular dynamics* (dcTMD) by @Wolf2018_dcTMD.
Given a set of non-equilibrium simulations where, for example, a ligand is pulled from a binding site via a velocity constraint, this tool performs automated post-processing of the bias-force time traces to estimate the underlying free-energy landscape and the friction (dissipation) along the unbinding coordinate.  

The method is based on a second-order cumulant expansion of *Jarzynski’s equality* [@Jarzynski1997_NonEq], which connects nonequilibrium work distributions to equilibrium free-energy differences. Combined with a Markovian Langevin Equation, dcTMD further allows the extraction of *position- and velocity-dependent friction coefficients* from the same nonequilibrium data. This approach has been succesfully applied in multiple studies [@Wolf2020_MultisecondDissociation @Jaeger2022_IonChannelConductance @Bray_2022_JCIM @Post2022_JCTC @Cai_2023_anisotropic @Taenzel2024_CommunityPaths @Jaeger2025_SimilarityMeasures @Milster_2025_NEQ].
The resulting free-energy and friction profiles can subsequently be used to estimate *binding and unbinding rate constants* following @Wolf2020_MultisecondDissociation.  

The software is intended for molecular dynamics practitioners interested in ligand–protein unbinding, mechanistic interpretation of binding kinetics, and quantitative modeling of non-equilibrium effects in soft condensed matter and biomolecular systems.  

# Statement of need
Ligand unbinding from proteins is of fundamental interest in computational biophysics [@Copeland_2016_drugtarget @Schuetz_2017_kinetics]. In many cases, the unbinding event is rare and requires enhanced-sampling or biased-simulation strategies to observe within computationally feasible timescales. 
dcTMD-based workflows have been shown to yield accurate free energy and non-equilibrium friction coefficients from velocity-constrained pulling simulations. The dcTMD package builds on this work by offering a unified, documented, and extensible implementation (including pathway separation analysis, see @Wolf_2023_path) that is currently not available, thereby lowering the barrier for applying dcTMD to new biomolecular systems and for reproducing published dcTMD studies.
​The *dcTMD* tool offers:  
- automatic parsing of Gromacs pulling simulation outputs,  
- estimation of work distributions from trajectory ensembles,  
- estimation of free-energy profiles along the biasing coordinate,  
- estimation of non-equilibrium friction coefficients along the same coordinate,  
- force autocorrelation analysis.

By providing a dedicated Python framework with an `scikit-learn`-style API, `dcTMD` enables users to integrate dissipation-corrected analysis into existing workflows, ensuring reproducibility and broad accessibility.
The software has already been successfully applied in several studies (e.g. [@Taenzel2024_CommunityPaths], [@Jaeger2025_SimilarityMeasures]), and is expected to promote the wider adoption of the dissipation-corrected targeted MD approach in computational chemistry and biophysics.  

# Implementation and architecture
The code is written in Python (versions 3.9–3.14) and is available under the MIT license.  
**Repository:** [https://github.com/moldyn/dcTMD](https://github.com/moldyn/dcTMD)

Key architectural features include:  
- A modular API following `fit`/`transform` conventions familiar from *scikit-learn*, easing integration into analysis pipelines.
- Input support for *GROMACS* pulling trajectories.
- Core functionality for computing free energy and non-equilibrium friction profiles along the biasing coordinate.  
- Support for analysis of multiple unbinding pathways.  
- Force correlation analysis for non-equilibrium friction analysis.  
- Continuous integration and testing via GitHub Actions; documentation hosted at [https://moldyn.github.io/dcTMD](https://moldyn.github.io/dcTMD).  

## Use case
A typical workflow begins with the user performing at least 100 independent velocity-constraint pulling simulations. `dcTMD` provides two analysis routes, both following the same workflow pattern:

1. Work-based analysis
   using a `WorkSet` and `WorkEstimator`

2. Force-correlation analysis
   using a `ForceSet` and `ForceEstimator`

Both methods yield free-energy and friction profiles but differ in how these properties are estimated.

#### 1. Load trajectories into a `WorkSet` or `ForceSet`
The user loads all pulling trajectories into an appropriate container:

* `WorkSet` for the work-based route, which is computationally cheaper, as the resolution of the trajectories can be reduced after integration.
* `ForceSet` for the force-correlation route.

Each trajectory contains the constraint force $f(t)$ along a pulling coordinate $x(t) = x_0 + vt$ with a constraint velocity $v$. The work along each trajectory is computed as $W(x) = \int_{x_0}^{x} \mathrm{d}x'\, f(x')$.

#### 2. Perform dcTMD analysis via an estimator
* **Work-based estimator (`WorkEstimator`)**

    The free-energy profile is estimated as $\Delta G(x) = \langle W(x) \rangle - \frac{\beta}{2}\langle \delta W(x)^2 \rangle$, with $\delta W = W - \langle W \rangle$, $\beta = (k_\mathrm{B}T)^{-1}$, and $\langle . \rangle$ denoting a mean over the ensemble of trajectories. The dissipated work is $W_\mathrm{diss}(x) = \frac{\beta}{2}\langle \delta W(x)^2 \rangle$. The non-equilibrium position-dependent friction is obtained from its derivative as $\Gamma(x) = \frac{1}{v} \, \frac{\mathrm{d}}{\mathrm{d}x} W_\mathrm{diss}(x)$.

* **Force-correlation-based estimator (`ForceEstimator`)**
  In this approach the $\Delta G$ and $\Gamma$ are computed directy from the force data as $\Delta G(x) = \int_{x_0}^{x} \mathrm{d}x' \,\langle f(x') \rangle \,-\, v \int_{x_0}^{x} \mathrm{d}x' \,\Gamma(x')$ and $\Gamma(x) = \beta \int_0^{t(x)} \mathrm{d}\tau \,  \langle \delta f(t(x)) \, \delta f(\tau) \rangle$. Furthermore, the two-time force autocorrelation function $C_t(\tau) = \langle \delta f(t(x)) \, \delta f(\tau) \rangle$, which corresponds to the ''Memory Kernel'' of the dcTMD methodology @Post2022, can be plotted for selected values of $t(x)$ to gain insight into timescales within the ''bath'' degrees of freedom constituting dissipation channels and thus friction sources.

#### 3. Visualize and interpret results
`dcTMD` provides plotting and export tools for:
* free-energy profiles $\Delta G(x)$
* friction profiles $\Gamma(x)$
* work distributions

In addition, the tool supports trajectory separation, enabling analysis of different dissociation routes (see @Wolf_2023_path).

#### 4. Example
![Figures crated using data taken from  @Wolf2020_MultisecondDissociation of trypsin-benzamidine unbinding.  a)-c) work distribution analysis. d) Decomposition of mean work $W_{\rm mean} = \langle W(x) \rangle$ into free energy $\Delta G(x)$ and dissipation work $W_{\rm diss}(x)$. e) non-equilibrium friction coefficient $\Gamma (x)$ along the pulling coordinate $x$.](figures/image2.png){width=\linewidth}

Figure 1 displays a common analysis of a set of unbinding trajectories from TMD simulations of the trypsin-benzamidine complex [@Wolf2020_MultisecondDissociation]. The analysis of the work distribution displays good agreement with a normal distribution at two different exaluated positions of the pulling coordinate $x$. The mean work $W_{\rm mean} = \langle W(x) \rangle$, which shows no features on its own, yields a free energy profile $\Delta G(x)$, which displays a clearly defined transition state at $x \approx 0.45$ nm as well as a bound state in form of a free energy minimum at $x \approx 0.0$ nm and an unbound continuum for $x > 0.6$ nm. The maximum in friction $\Gamma$ around $x = 0.5$ nm is indicative of changes in the hydration of both ligand and binding site.

# Impact
By providing an open and reproducible implementation of the dcTMD methodology, the software lowers the barrier for researchers to apply dissipation-corrected targeted molecular dynamics to ligand unbinding problems as well as other condensed soft matter systems. This enables broader exploration of nonequilibrium binding kinetics, supports mechanistic interpretation of frictional contributions, and provides access to position-dependent friction from targeted MD trajectories.  We anticipate the tool will be adopted in academic molecular simulation groups and in pharmaceutical research exploring unbinding free energies and kinetics.

# Acknowledgements
The implementation of dcTMD builds on the Python scientific stack, relying on **NumPy** for numerical operations, **Matplotlib** for visualization, and **Click** for the command-line interface.  We thank Gerhard Stock, Matthias Post and Georg Diez for valuable discussions, and Fabian Rohrbach and Leo Küchler for testing the software. This work has been supported by the Deutsche Forschungsgemeinschaft (DFG) via the Research Unit FOR 5099 *“Reducing complexity of nonequilibrium systems”* (project No. 431945604).  The authors acknowledge support by the High Performance and Cloud Computing Group at the Zentrum für Datenverarbeitung of the University of Tübingen and the Rechenzentrum of the University of Freiburg, as well as the state of Baden-Württemberg through bwHPC and the DFG through Grants Nos. INST 37/935-1 FUGG and INST 39/963-1 FUGG.

# References
