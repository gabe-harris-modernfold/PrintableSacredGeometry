"""
Horn torus: the special torus where the tube radius equals the revolution
radius (R == r), so the inner hole closes to a single point at the center.

Parametrization (u around the axis, v around the tube):
    x = (R + r cos v) cos u
    y = (R + r cos v) sin u
    z =       r sin v
With R == r, the inner equator (v = pi) collapses to the origin -> the
characteristic pinch point of a horn torus.
"""

import numpy as np
import trimesh

# ---- parameters ----------------------------------------------------------
R          = 30.0   # revolution radius == tube radius for a horn torus
NU         = 200    # segments around the main axis
NV         = 120    # segments around the tube
OUT        = "horn_torus.stl"
# --------------------------------------------------------------------------


def main():
    r = R
    u = np.linspace(0, 2 * np.pi, NU, endpoint=False)
    v = np.linspace(0, 2 * np.pi, NV, endpoint=False)
    uu, vv = np.meshgrid(u, v, indexing="ij")

    ring = R + r * np.cos(vv)
    x = ring * np.cos(uu)
    y = ring * np.sin(uu)
    z = r * np.sin(vv)
    verts = np.column_stack([x.ravel(), y.ravel(), z.ravel()])

    # faces over the (NU x NV) grid, wrapping in both directions
    faces = []
    for i in range(NU):
        for j in range(NV):
            a = i * NV + j
            b = ((i + 1) % NU) * NV + j
            c = ((i + 1) % NU) * NV + (j + 1) % NV
            d = i * NV + (j + 1) % NV
            faces.append([a, b, c])
            faces.append([a, c, d])
    faces = np.array(faces)

    mesh = trimesh.Trimesh(vertices=verts, faces=faces, process=False)
    mesh.merge_vertices()                       # weld the collapsed inner ring
    mesh.update_faces(mesh.nondegenerate_faces())  # drop zero-area pinch faces
    mesh.remove_unreferenced_vertices()
    mesh.fix_normals()

    print(f"horn torus R=r={R}: verts={len(mesh.vertices)}, faces={len(mesh.faces)}")
    print(f"watertight={mesh.is_watertight}, euler={mesh.euler_number}")
    mesh.export(OUT)
    sz = np.round(mesh.bounds[1] - mesh.bounds[0], 1)
    print(f"exported {OUT}: size={sz} mm")


if __name__ == "__main__":
    main()
