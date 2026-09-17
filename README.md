# Electromagnetic FEM Benchmark: FEniCS vs Palace

This repository contains a reproducible computational benchmarking of two Finite Element Method solvers: **FEniCS** and **Palace**. 

The benchmark validates the time-harmonic Helmholtz equation by simulating a transverse electromagnetic wave propagating through a 3D coaxial cable. The numerical results (computed using first-order Nédélec edge elements) are validated against the exact closed-form analytical solution.

## Repository Structure
* `run.py`: The main orchestrator script that runs the meshing and solver pipelines.
* `coaxial_cable.py`: Generates the parameterized 3D tetrahedral mesh using Gmsh.
* `testcase_fenics.py`: The FEniCS solver implementation, including the weak form setup and $L^2$ error computation.
* `environment.yml`: Conda environment file containing the exact dependency required for reproducibility.
* `palace_config.json`: The JSON configuration file defining the boundary conditions, materials, and solver settings for Palace.

## Installation & Setup

**1. Prerequisites.**
You must have [Miniconda](https://docs.anaconda.com/free/miniconda/) or Anaconda installed on your system. 

**2. FEniCS & Python Dependencies.**
Download the project to your local machine and navigate into the directory, then create the isolated environment using the provided configuration file. 
```bash
conda env create -f environment.yml
```
Then activate the environment:
```bash
conda activate thesis_env
```

**2. Palace Installation.**
Palace is a compiled, high-performance C++ solver and must be installed separately outside of the Conda environment.
You can download a pre-compiled binary or build it from source by following the instructions on the [official Palace GitHub repository](https://github.com/awslabs/palace). Once installed, open `run.py` and update the `palace_bin` variable to point to the exact path of your local Palace executable:
```python
palace_bin = "/path/to/your/palace/build/bin/palace"
```

**3. Run the experiment.**
Now you are able to run the experiment by launching
```bash
python run.py
```
