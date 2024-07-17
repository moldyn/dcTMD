### Input files

Download the [tutorial\_files.tar.gz](https://github.com/moldyn/dcTMD/blob/484da088ad8e5e91886d4d057b7c452e7b7d9aab/docs/tutorials/tutorial_files.tar.gz) and unpack via

```console
tar -xzvf ./tutorial_files.tar
```
You will find a folder with the following files:


```console
ls

3ptb_AMBER99SB_ben_pushEQUIBRUN.mdp 
3ptb_AMBER99SB_ben_pushRUN_v0.001.mdp
3ptb_AMBER99SB_ben.gro
3ptb_AMBER99SB_ben.top
3ptb_AMBER99SB_ben.ndx
3ptb_AMBER99SB_ben_Protein_chain_A.itp
3ptb_AMBER99SB_ben_Ion_chain_B.itp
3ptb_ben_H2_GMX_RESP.itp
posre_Protein_chain_A.itp
posre_Ion_chain_B.itp
posre_ben.itp

```

* Files for run input commands:
    + The *pushEQUIBRUN.mdp* file is an initial equilibration file for generating start simulation files with different initial velocity distributions.

    + The *pushRUN_v0.001.mdp* file is the respective command input for the non-equilibrium pulling. 
    
    + **3ptb** refers to the protein data base code of a trypsin structure, **AMBER99SB** to the employed force field, and **ben** to the benzamidine ligand.


* Structure file: 
    + *3ptb_AMBER99SB_ben.gro* in which the trypsin-benzamidine complex is equilibrated in TIP3P water with a physiological NaCl concentration and a single Ca2+ ion.

* Topologies and position restraint files:
    + *3ptb_AMBER99SB_ben.top*
    + *3ptb_AMBER99SB_ben_Protein_chain_A.itp*
    + *3ptb_AMBER99SB_ben_Ion_chain_B.itp*
    + *3ptb_ben_H2_GMX_RESP.itp*
    + *posre_Protein_chain_A.itp*
    + *posre_Ion_chain_B.itp*
    + *posre_ben.itp*

* Index file:
    + *3ptb_AMBER99SB_ben.ndx*
    + **Important:** the index file needs to include an anchor group from whose center of mass the ligand is pulled away (in this case: the group *[sheet]* containing C-alpha atoms from the central beta-sheet) and the ligand itself (or better, the heavy atoms of the ligand, here group *[ BEN_heavy ]*). If you want to create a respective anchor index for your own simulation problem, choose an anchor group that is tighly connected to the remainder of the protein (such as C-alpha atoms in alpha-helices and beta-sheets). The vector connecting the centers of mass of anchor and ligand needs to roughly point into the direction of a putative unbinding path.


### Carrying out pulling MD simulations

For the generation of the input structure in your own project, we advise you to carry out an initial NPT equilibration simulation of at least 10 ns length. Here, we have done this already for you and generated an equilibrated structure.

You will need a number (optimally between 100-200, but here we restrict ourselves to 10) of equilibrated trajectories with different initial velocity distributions. For this, generate an initial equilibration folder and the simulation start TPR files using:

```console
mkdir equib
cd equib/

for i in {000..009}
do
gmx grompp -f ../3ptb_AMBER99SB_ben_pushEQUIBRUN.mdp -c ../3ptb_AMBER99SB_ben.gro -r ../3ptb_AMBER99SB_ben.gro -p ../3ptb_AMBER99SB_ben.top -n ../3ptb_AMBER99SB_ben.ndx -o 3ptb_AMBER99SB_ben_pushEQUIBRUN_"$i".tpr -maxwarn 1 
done
```

and run the individual simulations via e.g.
```console
gmx mdrun -v -deffnm 3ptb_AMBER99SB_ben_pushEQUIBRUN_001
```

As these initial runs only require simulations of 0.1 ns length, they should be ready within a reasonably short time, i.e., some minutes.

When all equilibration simulations have been carried out, prepare a separate directory for the pulling simulations and the individual pulling input TPR files via:
```console
cd ..
mkdir v0.001
cd v0.001/

for i in {000..009}
do
gmx grompp -f ../3ptb_AMBER99SB_ben_pushRUN_v0.001.mdp -c ../equib/3ptb_AMBER99SB_ben_pushEQUIBRUN_"$i".gro -p ../3ptb_AMBER99SB_ben.top -n ../3ptb_AMBER99SB_ben.ndx -o 3ptb_AMBER99SB_ben_pushRUN_0.001_"$i".tpr
done
```
Note that the notation ***\_0.001\_*** stands for a velocity in Gromacs units of 0.001~nm/ps, i.e., 1~m/s. To our current experience, this is a sweet-spot velocity with the best trade-off between slow pulling and minimal computational effort. Run the simulations via e.g.:
```console
gmx mdrun -v -deffnm 3ptb_AMBER99SB_ben_pushRUN_0.001_000
```
These simulations will each require 1-2 hours on a modern workstation, so you better run them in parallel on a HPC cluster of your choice.

For all further analysis, you require the **3ptb\_AMBER99SB\_ben\_pushRUN\_0.001\_*\_pullf.xvg** files (with * denoting the respective run number).
