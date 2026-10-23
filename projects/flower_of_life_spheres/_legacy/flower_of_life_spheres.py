"""
Flower of Life built from spheres, exported to STL.

The Flower of Life is circles on a triangular (hexagonal) lattice:
a center, a ring of 6, then a ring of 12 -> 19 circles for 2 rings.
Neighboring centers are spaced one radius apart in the classic figure,
which makes the circles overlap and form the petal/vesica pattern.

Here each circle becomes a SPHERE on that same lattice. Two modes:

  SPACING_FACTOR = 1.0  -> centers one radius apart, spheres overlap
                          (true Flower of Life look; union -> solid).
  SPACING_FACTOR = 2.0  -> centers two radii apart, spheres just touch
                          (close-packed "bag of marbles" look).
"""

import numpy as np
import trimesh

# ---- parameters ----------------------------------------------------------
MODE            = "cubic"  # "hex" flat disk, "fcc" 3D ball, "cubic" NxNxN cube
SPHERE_RADIUS   = 10.0   # mm
RINGS           = 2      # hex mode: 1->7, 2->19, 3->37 spheres
SPACING_FACTOR  = 1.0    # 1.0 = each sphere's surface passes through neighbor centers (Flower of Life)
SUBDIVISIONS    = 3      # icosphere subdivisions; higher = smoother + bigger file
UNION           = True   # boolean-union into one watertight solid
# cubic-mode only:
CUBE_N          = 4      # spheres per edge (N x N x N)
# fcc-mode only:
FCC_RADIUS      = 3.5    # cluster radius, in lattice nearest-neighbor units
ENVELOPE        = "cube"    # "sphere" = rounded ball, "cube" = rounded cube
OUTPUT          = "flower_of_life.stl"
# --------------------------------------------------------------------------


def hex_lattice_points(rings, spacing):
    """Center points of a hex-packed disk: center + `rings` surrounding rings."""
    pts = []
    for q in range(-rings, rings + 1):
        for r in range(-rings, rings + 1):
            s = -q - r
            if max(abs(q), abs(r), abs(s)) <= rings:   # hex (axial) distance
                x = spacing * (q + r / 2.0)
                y = spacing * (r * np.sqrt(3) / 2.0)
                pts.append((x, y, 0.0))
    return np.array(pts)


def fcc_lattice_points(extent, spacing, envelope):
    """FCC close-packing: integer sites with (i+j+k) even, clipped to envelope.

    Nearest-neighbor distance equals `spacing`. `extent` is the cluster
    half-size in nearest-neighbor units.
    """
    a = spacing * np.sqrt(2.0)          # conventional cube edge for this nn dist
    n = int(np.ceil(extent * np.sqrt(2.0))) + 1
    cutoff = extent * spacing
    pts = []
    for i in range(-n, n + 1):
        for j in range(-n, n + 1):
            for k in range(-n, n + 1):
                if (i + j + k) % 2 != 0:
                    continue
                p = 0.5 * a * np.array([i, j, k])
                if envelope == "sphere":
                    if np.linalg.norm(p) <= cutoff + 1e-9:
                        pts.append(p)
                else:  # cube
                    if np.max(np.abs(p)) <= cutoff + 1e-9:
                        pts.append(p)
    return np.array(pts)


def cubic_lattice_points(n, spacing):
    """Simple-cubic N x N x N grid, centered at the origin."""
    offs = (np.arange(n) - (n - 1) / 2.0) * spacing
    xx, yy, zz = np.meshgrid(offs, offs, offs, indexing="ij")
    return np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])


def build():
    spacing = SPHERE_RADIUS * SPACING_FACTOR
    if MODE == "hex":
        centers = hex_lattice_points(RINGS, spacing)
    elif MODE == "cubic":
        centers = cubic_lattice_points(CUBE_N, spacing)
    else:
        centers = fcc_lattice_points(FCC_RADIUS, spacing, ENVELOPE)
    print(f"{MODE}: {len(centers)} spheres, r={SPHERE_RADIUS}, spacing={spacing:.2f}")

    spheres = []
    for c in centers:
        s = trimesh.creation.icosphere(subdivisions=SUBDIVISIONS, radius=SPHERE_RADIUS)
        s.apply_translation(c)
        spheres.append(s)

    if UNION:
        mesh = trimesh.boolean.union(spheres)
        print(f"union -> watertight={mesh.is_watertight}, "
              f"verts={len(mesh.vertices)}, faces={len(mesh.faces)}")
    else:
        mesh = trimesh.util.concatenate(spheres)
        print(f"concatenated {len(spheres)} meshes (not unioned)")

    mesh.export(OUTPUT)
    b = mesh.bounds
    print(f"exported {OUTPUT}  size={b[1]-b[0]} mm")


if __name__ == "__main__":
    build()
