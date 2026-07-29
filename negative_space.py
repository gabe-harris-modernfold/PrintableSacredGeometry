"""
Negative space (interstitial voids) of an FCC close-packing of spheres.

In FCC packing the gaps between spheres are the classic interstitial holes:
- TETRAHEDRAL voids: enclosed by 4 mutually-touching spheres (small, sharp).
- OCTAHEDRAL voids: enclosed by 6 spheres (larger).

We build the sphere union, subtract it from a box clipped to the lattice
interior, then split the result into connected components. Spheres are given
a slight overlap (FACTOR < 2.0) so each interior void seals into its own
closed concave solid instead of leaking through the tangent pinch-points.
"""

import numpy as np
import trimesh

# ---- parameters ----------------------------------------------------------
R          = 30.0   # sphere radius (mm)
FACTOR     = 1.70   # center spacing / R. voids seal into separate solids when < sqrt(3)~1.732
EXTENT     = 2      # FCC cube half-extent, in nearest-neighbor units
SUB        = 3      # icosphere subdivisions
OUTPUT     = "negative_space.stl"
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

    neg = trimesh.boolean.difference([box, solid])

    comps = neg.split(only_watertight=False)
    eps = 1e-3
    interior = [c for c in comps
                if np.all(c.bounds[0] > lo + eps) and np.all(c.bounds[1] < hi - eps)]
    print(f"{len(comps)} total components, {len(interior)} fully-interior voids")

    vols = sorted(c.volume for c in interior)
    if vols:
        print("interior void volumes (mm^3), sorted:",
              ", ".join(f"{v:.1f}" for v in vols))

    out = trimesh.util.concatenate(interior) if interior else neg
    out.export(OUTPUT)
    print(f"exported {OUTPUT}: verts={len(out.vertices)}, faces={len(out.faces)}")


if __name__ == "__main__":
    main()
