Analysis tools for dissipation-corrected targeted molecular dynamics, which is an enhanced sampling method to enforce rare events in biomolecular systems.

The module is structured into the following submodules:

- [**io:**][dcTMD.io] This submodule contains all methods related to reading data from text files and writing data to text files, including helpful header comments.

- [**storing:**][dcTMD.storing] This submodule creates force or work sets from .pullf files (force time traces genertaed by gromacs) which are needed for further analysis. It provides two classes, WorkSet and ForceSet, that store constraint force data as work or force time traces, respectively.

- [**dcTMD:**][dcTMD.dcTMD] This submodule contains two classes, WorkEstimator and ForceEstimator, which are used for the dcTMD analysis of constraint force time traces. Both class can be used to calculate the mean work, dissipative work, free energy and friction estimate of a set of constraint force time traces.

- [**utils:**][dcTMD.utils] This submodule provides utility functions such as smoothing, plotting and error estimation via bootstrapping. The functions in this submodule can be used in conjunction with other parts of the software.
