from mpi4py import MPI
import ufl
from dolfinx import mesh, fem, io
from dolfinx.fem.petsc import LinearProblem
import numpy as np
from dolfinx.io import gmsh as d_gmsh
from dolfinx.io import XDMFFile
import resource
import gmsh
from petsc4py import PETSc

# 1. load the mesh
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0) # suppress terminal output
mesh_data = d_gmsh.read_from_msh("cable.msh", MPI.COMM_WORLD, gdim=3)
domain = mesh_data.mesh
facet_tags = mesh_data.facet_tags

# 2. define the Nédélec edge element space of order 1 (forced to complex)
dtype_complex = np.complex128
V = fem.functionspace(domain, ("N1curl", 1))

# parameters of the problem
frequency = 30e6                # in this way the wavelength is 10 m
omega = 2.0 * np.pi * frequency
mu_r = 1.0
epsilon_r = 1.0
c0 = 3e8                        # speed of light in vacuum
k0 = omega / c0
k_squared = k0**2 * epsilon_r

# trial and test functions
E = ufl.TrialFunction(V)
v = ufl.TestFunction(V)

# integration measure for the boundary facets
ds = ufl.Measure("ds", domain=domain, subdomain_data=facet_tags)

# 4. weak form definition (Helmholtz Equation)
# (\nabla \times E, \nabla \times v) - k^2 (\epsilon_r E, v) = (F, v) + boundary terms (absorbing boundary conditions)
a = ufl.inner(ufl.curl(E), ufl.curl(v)) * ufl.dx - k_squared * ufl.inner(E, v) * ufl.dx
im_k0 = fem.Constant(domain, 1j * k0)
a += im_k0 * ufl.inner(E, v) * ds(3)
f = fem.Constant(domain, np.array([0.0, 0.0, 0.0], dtype=dtype_complex))
L = ufl.inner(f, v) * ufl.dx

print("Finite element space and weak form configured successfully!")

# 5. inlet boundary condition

# define the analytical solution for the trnasverse electromagnetic mode in a coaxial cable,
#  this is the electrostatic field between two concentric cylinders, which solves the 2D Laplace equation
def tem_mode(x):
    r_in = 1.0
    r_out = 3.0
    V0 = 1.0
    amplitude = V0 / np.log(r_out / r_in)
    
    r_sq = x[0]**2 + x[1]**2
    E_val = np.zeros((3, x.shape[1]), dtype=dtype_complex)
    E_val[0] = amplitude * (x[0] / r_sq)
    E_val[1] = amplitude * (x[1] / r_sq)
    return E_val

# interpolate the analytical solution on the inlet boundary as a Nédélec function, which will be used to enforce the Dirichlet boundary condition on the inlet
u_in = fem.Function(V)
u_in.interpolate(tem_mode)

# apply Dirichlet boundary condition on the inlet boundary
fdim = domain.topology.dim - 1
inlet_facets = facet_tags.find(2)
inlet_dofs = fem.locate_dofs_topological(V, fdim, inlet_facets)
bc_inlet = fem.dirichletbc(u_in, inlet_dofs)

# 6. Dirichlet boundary condition on the tangential component of the electric field on the boundary
inner_wall_facets = facet_tags.find(4)
outer_wall_facets = facet_tags.find(5)

wall_facets = np.concatenate([inner_wall_facets, outer_wall_facets])
wall_dofs = fem.locate_dofs_topological(V, fdim, wall_facets)

u_pec = fem.Function(V)
u_pec.x.array[:] = 0.0
bc_pec = fem.dirichletbc(u_pec, wall_dofs)

# 7. solve the problem 
problem = LinearProblem(a, L, bcs=[bc_pec, bc_inlet], petsc_options_prefix="helmholtz_solver")

# force the existing solver to use GMRES
problem.solver.setType(PETSc.KSP.Type.GMRES)
problem.solver.getPC().setType(PETSc.PC.Type.ILU)
problem.solver.setTolerances(rtol=1e-8, max_it=1000)

# Print the residual at each iteration like Palace
def gmres_monitor(ksp, its, rnorm):
    print(f"  FEniCS GMRES Iteration {its} | Residual: {rnorm:.4e}")

problem.solver.setMonitor(gmres_monitor)

E_h = problem.solve()

print("Weak form solved successfully!")

num_dofs = V.dofmap.index_map.size_global * V.dofmap.index_map_bs
mem_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)

print(f"fenics_DOFS: {num_dofs}")
print(f"fenics_RAM: {mem_mb:.2f}")

'''
# Paraview export
# interpolate the Nédélec solution onto a node-based Lagrange space
W = fem.functionspace(domain, ("Lagrange", 1, (3,)))
E_W_complex = fem.Function(W)
E_W_complex.interpolate(E_h)

# Create a Real-valued function specifically for XDMF export
E_vis = fem.Function(W, dtype=np.float64)
E_vis.name = "Electric_Field_Real"

# Copy only the real part of the interpolated data for visualization
E_vis.x.array[:] = np.real(E_W_complex.x.array)

with XDMFFile(domain.comm, "electric_field.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    xdmf.write_function(E_vis)
'''

'''
# Computation of the L2 error with respect to the analytical solution
# define the spatial coordinates
x = ufl.SpatialCoordinate(domain)
r_sq = x[0]**2 + x[1]**2

# construct the exact 3D phasor: e^{-j * k0 * z} * (x/r^2, y/r^2, 0)
V0 = 1.0
r_in = 1.0
r_out = 3.0
amplitude = V0 / np.log(r_out / r_in)
phase = ufl.cos(k0 * x[2]) - 1j * ufl.sin(k0 * x[2])
E_exact = ufl.as_vector([
    amplitude * phase * (x[0] / r_sq),
    amplitude * phase * (x[1] / r_sq),
    fem.Constant(domain, 0.0 + 0.0j)
])

# formulate the L2 error integral
error_form = fem.form(ufl.inner(E_h - E_exact, E_h - E_exact) * ufl.dx)
exact_norm_form = fem.form(ufl.inner(E_exact, E_exact) * ufl.dx)

# 4. Assemble the scalar integrals across the mesh
local_error_sq = fem.assemble_scalar(error_form)
local_exact_sq = fem.assemble_scalar(exact_norm_form)

global_error_sq = domain.comm.allreduce(local_error_sq, op=MPI.SUM)
global_exact_sq = domain.comm.allreduce(local_exact_sq, op=MPI.SUM)

abs_L2_error = np.sqrt(np.real(global_error_sq))
rel_L2_error = abs_L2_error / np.sqrt(np.real(global_exact_sq))

print(f"Absolute L2 Error: {abs_L2_error:.4e}")
print(f"Relative L2 Error: {rel_L2_error * 100:.2f} %")
'''