# About this module

This module was developed during Julian Sierra's research stay on (March-April 2026) with Prof. Stefan Boresch. It focuses on implementing and comparing two anchor-search algorithm for Boresch restraints using OpenMM, MDAnalysis and OpenFE.

The goal of this toolkit is to support reproducible testing and benchmarking of restraint placement strategies in moleculer dynamic simulations. 

# How to install this module

To use this module, it is recomended to create a dedicated Conda environment using OpenFE 1.9.1 and Python 3.13, which were the versions used during development:

```bash 
conda create -c conda-forge -n openfe openfe=1.9.1 python=3.13
conda activate openfe
```

Once the environment is activated, install package from the source:

```bash 
git clone https://github.com/JulianSierraRodriguez/boresch-restraint-vienna-kit.git
cd boresch-restraint-vienna-kit/
pip install -e .
```

# Available modules

The package is organized into several functional modules:

## simulations_openMM

This module containts all utilities required for system preparation and molecular dynamics simulations using OpenMM.

It provides three main workflows:

* ```preparing_from_PDBank``` : prepares a structure from the Protein Data Bank.
* ```preparation_simulations``` : generates OpenMM system and runs preparation steps including minimization, NVT heating and NPT equilibration.
* ```production_simulations```, performs production molecular dynamics using checkpointing for restartability.

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