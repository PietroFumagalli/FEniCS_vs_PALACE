import os
import sys
import subprocess
import time
import os

print("COAXIAL CABLE FEM BENCHMARK (Fenics vs Palace)")

# 1. MESH GENERATION
print("(1/3) Generating Mesh ...")
t0_mesh = time.perf_counter()
subprocess.run([sys.executable, "coaxial_cable.py"], check=True)
t1_mesh = time.perf_counter()
print(f"Mesh generated in {t1_mesh - t0_mesh:.4f} seconds.\n")

# 2. FENICS SOLVER
print("(2/3) Running Fenics...")
t0_fenics = time.perf_counter()
result_fenics = subprocess.run([sys.executable, "testcase_fenics.py"], check=True)
t1_fenics = time.perf_counter()
fenics_time = t1_fenics - t0_fenics
print(f"Finished in {fenics_time:.4f} seconds.")

# 3. PALACE SOLVER
print("\n(3/3) Palace Solver: Pending.")