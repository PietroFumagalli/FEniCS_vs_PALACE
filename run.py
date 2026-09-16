import os
import sys
import subprocess
import time
import os
import result

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
print(f"Finished in {t1_fenics - t0_fenics:.4f} seconds.")

# 3. PALACE SOLVER
print("\n(3/3) Running Palace...")
palace_bin = os.path.expanduser("~/Desktop/temp/palace_source/build/bin/palace")
t0_palace = time.perf_counter()
subprocess.run([palace_bin, "palace_config.json"], check=True)
t1_palace = time.perf_counter()
print(f"Finished in {t1_palace - t0_palace:.4f} seconds.")