import meshio
msh = meshio.read("data/mesh1.msh")

def create_mesh(mesh, cell_type, prune_z=False):
    
    cells = mesh.get_cells_type(cell_type)
    cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
    out_mesh = meshio.Mesh(
        mesh.points[:,:2], cells={cell_type: cells}, 
        cell_data={"name_to_read":[cell_data]}
    )
    return out_mesh

mesh = create_mesh(msh, "triangle", prune_z=True)
meshio.write("data/mesh1.xdmf", mesh)