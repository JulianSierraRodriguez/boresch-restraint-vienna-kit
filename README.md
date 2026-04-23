# About this module

This module was developed during Julian Sierra's research stay on March-April 2026 with Prof. Stefan Boresch. This module is focused on testing and comparing two anchor search algorithm for Boresch restraints using OpenMM, MDAnalysis and OpenFE.

# How to install this module

To be able to use this module the best thing is to create a conda environment with the OpenFE 1.9.1 and python 3.13, which was the one existing when this module was being developed.

```bash 
conda create -c conda-forge -n openfe openfe=1.9.1 python=3.13
conda activate openfe
```

Once we are in the conda environment, we install this module using pip inside the folder:

```bash 
git clone https://github.com/JulianSierraRodriguez/boresch-restraint-vienna-kit.git
cd boresch-restraint-vienna-kit/
pip install -e .
```

# Possible modules

There several possible functions divided into different modules:

## simulations_openMM

## restraints_openfe

## restraints_william

## plot restraints_mda

## drawing_boresch_restraints


# Examples

In the same folder where you perform the pip install command, there is a directory named "examples", inside here we have several with the needed inputs to run them in the directory.

## simple_example

This test is a simple example test, where we start with a CIF and SDF from the Protein Data Bank and perform a short simulation, then it searches for anchors and draws them using RDKit.

## comparison_search_algorithms

Inside the directory "comparison_search_algorithms", we have a script that performs ten times the MD procedure, and then uses two search algorithms on the same trajectories to compare the results. The algorithms are:

-  OpenFE 1.9.1 boresch restraint search from the ABFE workflow.
-  Boresch restraint search by Wu et al. (https://doi.org/10.1021/acs.jctc.5c00861)

This script performs the search of both algorithms and then plots the restraints on the residues involved using RDKit.

To run the script get inside the comparison_search_algorithms and run:

```bash 
python comparison_W_openfe.py
```

This will make 10 new directories where it will run the simulations and save the figures, and at the end will produce several spreadsheets with the results, to compare both algorithms.

## long_simulation

```python 
your_code = do_some_stuff
```