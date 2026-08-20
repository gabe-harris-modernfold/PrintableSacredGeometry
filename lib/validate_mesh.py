#!/usr/bin/env python3
"""Validate an STL for 3D printing using PyMeshLab.

Usage:  python validate_mesh.py [path-to.stl]   (default: vortex_horn_torus.stl)

Checks: watertight (closed, 0 boundary edges), 2-manifold, single component,
no self-intersecting faces, and outward-facing normals (positive volume).

NOTE (Windows / Microsoft Store Python): pymeshlab ships its Qt/MeshLab DLLs in
the package root but doesn't always register that dir on the DLL search path
before loading its plugins, so a bare `import pymeshlab` can fail with
"Cannot load library ... io_e57.dll". We add the package dir explicitly first.
"""
import sys, os, importlib.util

spec = importlib.util.find_spec("pymeshlab")
if spec and spec.submodule_search_locations:
    os.add_dll_directory(spec.submodule_search_locations[0])
import pymeshlab as ml

path = sys.argv[1] if len(sys.argv) > 1 else "vortex_horn_torus.stl"
ms = ml.MeshSet()
ms.load_new_mesh(path)
m = ms.current_mesh()

topo = ms.apply_filter("get_topological_measures")
geo  = ms.apply_filter("get_geometric_measures")
ms.apply_filter("compute_selection_by_self_intersections_per_face")
self_int = ms.current_mesh().selected_face_number()

vol      = geo.get("mesh_volume", 0.0)
manifold = bool(topo.get("is_mesh_two_manifold", False))
boundary = topo.get("boundary_edges", -1)
nonman_e = topo.get("non_two_manifold_edges", -1)
nonman_v = topo.get("non_two_manifold_vertices", -1)
comps    = topo.get("connected_components_number", -1)

print(f"file:        {path}")
print(f"vertices:    {m.vertex_number()}    faces: {m.face_number()}")
print(f"two-manifold:{manifold}")
print(f"boundary edges:        {boundary}")
print(f"non-manifold edges:    {nonman_e}")
print(f"non-manifold vertices: {nonman_v}")
print(f"connected components:  {comps}")
print(f"self-intersecting faces: {self_int}")
print(f"volume (mm^3): {vol:.1f}    surface area (mm^2): {geo.get('surface_area',0):.1f}")

checks = {
    "watertight (0 boundary edges)": boundary == 0,
    "2-manifold":                    manifold,
    "single component":              comps == 1,
    "no self-intersections":         self_int == 0,
    "outward normals (vol > 0)":     vol > 0,
}
print("\n" + "-" * 40)
for name, ok in checks.items():
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
print("-" * 40)
print("RESULT:", "PRINT-READY" if all(checks.values()) else "NEEDS ATTENTION")
sys.exit(0 if all(checks.values()) else 1)
