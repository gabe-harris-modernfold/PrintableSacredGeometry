"""
Take the 78 disconnected void shapes from negative_space.stl and tie them
together into one printable piece: a base plate with a vertical post rising
into each shape's centroid. Shapes stacked in the same column share a post
(skewer), so everything ends up as a single connected solid.
"""

import numpy as np
import trimesh

# ---- parameters ----------------------------------------------------------
SRC          = "negative_space.stl"
OUT          = "negative_space_posts.stl"
POST_R       = 1.5    # post radius (mm)
BASE_THICK   = 5.0    # base plate thickness (mm)
BASE_MARGIN  = 6.0    # base overhang beyond the shape footprint (mm)
MIN_POST     = 1.0    # skip posts shorter than this (shape already on the base)
# --------------------------------------------------------------------------


def main():
    m = trimesh.load(SRC)
    comps = m.split(only_watertight=False)
    lo, hi = m.bounds
    z0 = lo[2]                       # top surface of the base = lowest shape point

    parts = list(comps)

    # base plate
    ext = np.array([hi[0] - lo[0] + 2 * BASE_MARGIN,
                    hi[1] - lo[1] + 2 * BASE_MARGIN,
                    BASE_THICK])
    base = trimesh.creation.box(extents=ext)
    base.apply_translation([(lo[0] + hi[0]) / 2,
                            (lo[1] + hi[1]) / 2,
                            z0 - BASE_THICK / 2])
    parts.append(base)

    # one post per shape, base up to the shape's centroid
    n_posts = 0
    for c in comps:
        cx, cy, cz = c.centroid
        if cz - z0 < MIN_POST:
            continue
        post = trimesh.creation.cylinder(
            radius=POST_R, segment=[[cx, cy, z0 - 0.1], [cx, cy, cz]])
        parts.append(post)
        n_posts += 1

    result = trimesh.boolean.union(parts)
    ncomp = len(result.split(only_watertight=False))
    print(f"{len(comps)} shapes + {n_posts} posts + base -> "
          f"{ncomp} connected component(s), watertight={result.is_watertight}")

    result.export(OUT)
    sz = np.round(result.bounds[1] - result.bounds[0], 1)
    print(f"exported {OUT}: size={sz} mm, verts={len(result.vertices)}, "
          f"faces={len(result.faces)}")


if __name__ == "__main__":
    main()
