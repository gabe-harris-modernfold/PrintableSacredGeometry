"""
Sponge: the CONNECTED negative space of a near-tangent FCC sphere packing.

Same idea as negative_space.py, but the spheres are spaced just under tangent
(FACTOR between sqrt(3)~1.732 and 2.0). Above the 1.732 sealing threshold the
triangular windows between spheres stay OPEN, so all the tetrahedral and
octahedral pockets link into a single connected porous solid -> one printable
piece, like a sponge / open-cell foam.
"""

import numpy as np
import trimesh

# ---- parameters ----------------------------------------------------------
R          = 30.0   # sphere radius (mm)
FACTOR     = 1.90   # center spacing / R. 2.0 = tangent; 1.732..2.0 = open windows (sponge)
EXTENT     = 2      # FCC cube half-extent, in nearest-neighbor units
SUB        = 3      # icosphere subdivisions
OUTPUT     = "sponge.stl"
# --------------------------------------------------------------------------


def fcc_centers(extent, spacing):
    a = spacing * np.sqrt(2.0)
    n = int(np.ceil(extent * np.sqrt(2.0))) + 1
    cutoff = extent * spacing
    pts = []
    for i in range(-n, n + 1):
        for j in range(-n, n + 1):
            for k in range(-n, n + 1):
                if (i + j + k) % 2:
                    continue
                p = 0.5 * a * np.array([i, j, k])
                if np.max(np.abs(p)) <= cutoff + 1e-9:
                    pts.append(p)
    return np.array(pts)


def main():
    spacing = R * FACTOR
    centers = fcc_centers(EXTENT, spacing)

    spheres = []
    for c in centers:
        s = trimesh.creation.icosphere(subdivisions=SUB, radius=R)
        s.apply_translation(c)
        spheres.append(s)
    solid = trimesh.boolean.union(spheres)
    print(f"{len(centers)} spheres, spacing={spacing:.2f}, overlap={2*R-spacing:.2f} mm")

    lo, hi = centers.min(0), centers.max(0)
    box = trimesh.creation.box(extents=(hi - lo))
    box.apply_translation((lo + hi) / 2.0)

    sponge = trimesh.boolean.difference([box, solid])

    comps = sponge.split(only_watertight=False)
    print(f"connected components: {len(comps)} (1 = single sponge), "
          f"watertight={sponge.is_watertight}")
    sponge.export(OUTPUT)
    sz = np.round(sponge.bounds[1] - sponge.bounds[0], 1)
    print(f"exported {OUTPUT}: size={sz} mm, verts={len(sponge.vertices)}, "
          f"faces={len(sponge.faces)}")


if __name__ == "__main__":
    main()
