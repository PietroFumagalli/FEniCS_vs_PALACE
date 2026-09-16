# Electromagnetic FEM Benchmark: FEniCS vs Palace

This repository contains a robust, reproducible computational benchmarking suite comparing two Finite Element Method (FEM) solvers: **FEniCS** and **Palace**. 

The benchmark validates the time-harmonic Helmholtz equation by simulating a transverse electromagnetic wave propagating through a 3D coaxial cable. The numerical results (computed using first-order Nédélec edge elements) are validated against the exact closed-form analytical solution.

## Repository Structure
* `run.py`: The main orchestrator script that runs the meshing and solver pipelines.
* `coaxial_cable.py`: Generates the parameterized 3D tetrahedral mesh using Gmsh.
* `testcase_fenics.py`: The FEniCS solver implementation, including the weak form setup and $L^2$ error computation.
* `environment.yml`: Conda environment file containing the exact dependency pins required for reproducibility.

## Installation & Setup

**1. Prerequisites**
You must have [Miniconda](https://docs.anaconda.com/free/miniconda/) or Anaconda installed on your system. 

**2. Installation**
Download the project to your local machine and navigate into the directory, then create the isolated environment using the provided configuration file. 
```bash
conda env create -f environment.yml
```
Then activate the environment:
```bash
conda activate thesis_env
```
And then you are able to run the experiment by launching
```bash
python run.py
```
