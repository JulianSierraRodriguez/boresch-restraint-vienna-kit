# About this module



# How to install this module

To be able to use this module the best thing is to create a conda environment with the OpenFE 1.9.1 and python 3.13, which was the one existing when this module was being developed.

```bash 
conda create -c conda-forge -n openfe openfe=1.9.1 python=3.13
conda activate openfe
```
Once we are in the conda environment, we install this module using pip in the folder with pyproject.toml:

```bash 
git clone https://github.com/JulianSierraRodriguez/boresch-restrain-vienna-kit.git
cd boresch-restrain-vienna-kit/boresch-restraint-vienna-kit/
pip install -e .
```

# Possible modules

## simulation_julian

## restraints_openfe

## restraints_william

## plot restraints_mda

## drawing_boresch_restraints


# Examples

## comparison_search_algorithms

## long_simulation

## simple_example

```python 
your_code = do_some_stuff
```