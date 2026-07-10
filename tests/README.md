# Test and tutorial data

## Origin

The `testdata/t_middle_*_pullf.xvg` files contain force time traces from
targeted molecular dynamics (dcTMD) simulations of the unbinding of 
benzamidine from trypsin. They are a reduced subset of the
simulation ensemble published in

> S. Wolf, B. Lickert, S. Bray, and G. Stock,
> *Multisecond ligand dissociation dynamics from atomistic simulations*,
> Nat. Commun. **11**, 2918 (2020), doi: 10.1038/s41467-020-16655-1

## Simulation setup

For details have a look at the supplementary information. 

| Property | Value |
| --- | --- |
| System | trypsin, PDB ID 3PTB, in complex with benzamidine |
| Box | dodecahedral, 7.5 x 7.5 x 5.3 nm, 8971 water molecules |
| Ions | 16 Na+ and 25 Cl-, charge neutral at 0.1 M salt |
| MD engine | Gromacs v2018 |
| Force field | Amber99SB* for protein and ions, GAFF for benzamidine |
| Ligand charges | RESP, from HF/6-31G* calculations |
| Water model | TIP3P |
| Temperature | 290.15 K |
| Pressure | 1 bar, Parrinello-Rahman barostat during pulling |
| Bias | moving distance constraint, Gromacs PULL code, SHAKE |
| Pulling coordinate | distance between the center of mass of the benzamidine heavy atoms and that of the Calpha atoms of the central beta sheet of trypsin |
| Time step for integration | 1 fs |
| Number of steps that elapse between writing forces to the output trajectory | 1 |
| Pulling velocity | 1 m/s, equivalently 0.001 nm/ps |
| Pulling distance | 2 nm, reached after 2000 ps |

## File format

Each `.xvg` file holds two columns in Gromacs internal units. The first column
is the simulation time in ps, the second the constraint force in kJ/(mol nm).
Header lines starting with `#` or `@` are ignored by `dcTMD`. Every file
contains 20001 frames, one every 0.1 ps, covering the full 2000 ps of the
pulling simulation.

## Scope and limitations

The `testdata/` directory contains 18 force trajectories, whereas the published 
work used 400 statistically independent trajectories. Free energy and friction
estimates from dcTMD require an ensemble of this size to converge. The bundled
files are intended to exercise the code and to let users run the analysis
tutorials on a laptop without access to HPC resources. They are not sufficient
to reproduce the free energy and friction profiles reported in the reference
above, and results obtained from them should not be interpreted physically.