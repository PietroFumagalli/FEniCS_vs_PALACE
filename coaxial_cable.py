import gmsh

gmsh.initialize()
gmsh.option.setNumber("Mesh.Algorithm3D", 1)
gmsh.option.setNumber("General.Terminal", 0) # suppress terminal output
gmsh.model.add("coaxial_cable")

# Cap the maximum element size for wave accuracy (lambda/10)
gmsh.option.setNumber("Mesh.MeshSizeMax",0.2)

# Allow Gmsh to adaptively use smaller elements down to this size for curved surfaces
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.1)

L = 10.0      # cable length
r_in = 1.0    # inner radius
r_out = 3.0   # outer radius

ext = gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, L, r_out)
int = gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, L, r_in)
gmsh.model.occ.cut([(3, ext)], [(3, int)])

gmsh.model.occ.synchronize()

# Define physical group for the dielectric volume
volumes = gmsh.model.getEntities(3)
gmsh.model.addPhysicalGroup(3, [volumes[0][1]], tag=1)
gmsh.model.setPhysicalName(3, 1, "Dielectric")

# Sort the surfaces to assign the correct boundary tags
surfaces = gmsh.model.getEntities(2)
inlet_surfs = []
outlet_surfs = []
inner_surfs = []
outer_surfs = []

for dim, tag in surfaces:
    com = gmsh.model.occ.getCenterOfMass(dim, tag)
    if abs(com[2] - 0.0) < 1e-3:       
        inlet_surfs.append(tag)        # z = 0
    elif abs(com[2] - L) < 1e-3:       
        outlet_surfs.append(tag)       # z = L
    else:
        # If it's not the ends, it's a cylinder wall
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
        if abs(xmax - r_out) < 1e-3:
            outer_surfs.append(tag)    # r = 3.0
        else:
            inner_surfs.append(tag)    # r = 1.0

# 3. Apply the 2D Physical Groups
gmsh.model.addPhysicalGroup(2, inlet_surfs, tag=2)
gmsh.model.setPhysicalName(2, 2, "Inlet")

gmsh.model.addPhysicalGroup(2, outlet_surfs, tag=3)
gmsh.model.setPhysicalName(2, 3, "Outlet")

gmsh.model.addPhysicalGroup(2, inner_surfs, tag=4)
gmsh.model.setPhysicalName(2, 4, "Inner_PEC")

gmsh.model.addPhysicalGroup(2, outer_surfs, tag=5)
gmsh.model.setPhysicalName(2, 5, "Outer_PEC")

# Generate and save the mesh
gmsh.model.mesh.generate(3)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write("cable.msh")

# gmsh.fltk.run()
gmsh.finalize()
