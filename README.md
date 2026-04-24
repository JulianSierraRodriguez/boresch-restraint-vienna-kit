# About this package

This module was developed during Julian Sierra's research stay on (March-April 2026) with Prof. Stefan Boresch. It focuses on implementing and comparing two anchor-search algorithm for Boresch restraints using OpenMM, MDAnalysis and OpenFE.

The goal of this toolkit is to support reproducible testing and benchmarking of restraint placement strategies in moleculer dynamic simulations. 

# Installation

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

This module implements a function that uses OpenFE's anchor-search algorithm for Boresch restraints in their ABFE workflow.

After execution, it returns the six anchors atoms for the Boresch restraints ```[h0,h1,h2,g0,g1,g2]```, along with the equilibrium parameters required to define the restraints.

## restraints_william

This module containts utilities required to run the anchor-search algorith from Wu et al. (https://doi.org/10.1021/acs.jctc.5c00861) in the workflow ```restraint_search_william```.

After execution, it returns the six anchors atoms for the Boresch restraints ```[h0,h1,h2,g0,g1,g2]```, along with the equilibrium parameters required to define the restraints.

## plot restraints_mda

This modules provides utilities to compute and visualize Boresch restraint Parameters along a trajectory.

Using MDAnalysis, it evaluates restraint geometries over time and generates plots using Matplotlib.

## drawing_boresch_restraints

This modules provides tools to visualize selected anchor atoms in a 2D representation of ligand-Protein complex using RDKit.

It helps map the selected anchors onto chemical structures for interpretability.

# Examples

The ```examples/``` directory (located in the repository root after installation) contains ready-to-run workflows with the required input files.

## simple_example

This is a minimal test case that starts from CIF and SDF files obtained from the Protein Data Bank. It performs a short simulation, runs anchor selection, and visualizes the resulting anchors using RDKit.

To run the simple test:
```bash 
python simple_example.py
```

## comparison_search_algorithms

This directory contains, we have a script that performs ten times the MD procedure, and then uses the two anchor-search algorithms on the same trajectories to compare the results.

This script performs the search of both algorithms and then plots the restraints on the residues involved using RDKit.

To run the comparison (Warning, it takes time):

```bash 
python comparison_W_openfe.py
```

This generates 10 simulation runs in separate directories, saves visualizations, and produces summary spreadsheets for comparison of both methods.

## long_simulation

This directory contains scripts for longer and more detailed simulations:
In the 'long_simulation' directory we have several scripts:

* ```long_simulation.py``` : runs extended MD simulations using checkpointing vi ```simulation_openMM```.
* ```analysis_long_sim.py``` : extract the first 1ns of the simulation to find anchors with both algorithms, draw the found anchors with RDKit and plot the Boresch Parameters throughout the trajectory using:
  * ```restraints_william```
  * ```restraints_openfe```
  * ```drawing_boresch_restraints```
  * ```plot_restraints_mda```
* ```prepare_sim_for_VMD.py``` : combines trajectories into a single file, reduces size, removes the waters and centers the protein for visualization in VMD.

To run the long simulation (Warning, it takes time):

```bash 
python long_simulation.py
```

To run the analysis for the long simulation (Warning, it takes time and you need the outputs from long_simulation.py):

```bash 
python analysis_long_sim.py
```

To run the preparation of VMD for the long simulation (Warning, you need the outputs from long_simulation.py):

```bash 
python prepare_sim_for_VMD.py
```