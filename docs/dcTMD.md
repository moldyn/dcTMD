Analysis tools for dissipation-corrected targeted molecular dynamics, which is an enhanced sampling method to enforce rare events in biomolecular systems.

The module is structured into the following submodules:

- [**io:**][dcTMD.io] This submodule contains all methods related to reading data from text files and writing data to text files, including helpful header comments.

- [**storing:**][dcTMD.storing] This submodule offers techniques for the analysis of state trajectories&mdash;commonly known as Molecular Dynamics (MD)&mdash;without relying on Markov state models. It encompasses functions for determining timescales, recognizing significant events, correcting dynamical anomalies, and evaluating various state discretization methods.  These functions provide a comprehensive solution for analyzing time-series data and understanding the underlying dynamics of complex systems.

- [**dcTMD:**][dcTMD.dcTMD] This submodule contains two classes, WorkEstimator and ForceEstimator, which are used for the dcTMD analysis of constraint force time traces. Both class can be used to calculate the mean work, dissipative work, free energy and friction estimate of a set of constraint force time traces.

- [**utils:**][dcTMD.utils] This submodule provides utility functions such as smoothing, plotting and error estimation via bootstrapping. The functions in this submodule can be used in conjunction with other parts of the software.