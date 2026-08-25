"""Validate the vortex STLs: manifold, watertight, dimensions."""
import pymeshlab
import os

folder = os.path.dirname(os.path.abspath(__file__))
files = ["vortex_base.stl", "vortex_top_collar.stl", "vortex_tube_printable.stl",
         "wind_base.stl", "wind_venturi_head.stl", "wind_venturi_hat.stl"]

for f in files:
    path = os.path.join(folder, f)
    ms = pymeshlab.MeshSet()
    ms.load_new_mesh(path)
    m = ms.current_mesh()
    topo = ms.get_topological_measures()
    geo = ms.get_geometric_measures()
    bb = m.bounding_box()
    dim = bb.max() - bb.min()
    print(f"=== {f} ===")
    print(f"  vertices={m.vertex_number()}  faces={m.face_number()}")
    print(f"  bbox: {dim[0]:.2f} x {dim[1]:.2f} x {dim[2]:.2f} mm")
    print(f"  boundary_edges={topo.get('boundary_edges')}")
    print(f"  non_two_manifold_edges={topo.get('non_two_manifold_edges')}")
    print(f"  non_two_manifold_vertices={topo.get('non_two_manifold_vertices')}")
    print(f"  connected_components={topo.get('connected_components_number')}")
    print(f"  genus={topo.get('genus')}")
    vol = geo.get('mesh_volume')
    if vol is not None:
        print(f"  volume={vol/1000:.1f} cm^3")
    print()
