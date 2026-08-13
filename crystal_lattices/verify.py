#!/usr/bin/env python3
"""Independent check of every STL `build.py` wrote, using trimesh rather than the code
that generated them.

What "valid" means here
-----------------------
These parts are unions of overlapping watertight solids -- a hub plus its sleeves, a cell
plus its embossed glyph prisms -- because this environment has no mesh-boolean backend to
fuse them with. So `trimesh.load(path).is_watertight` is False for a labelled part *by
construction*, and checking it on the whole file would be checking the wrong thing.

The right check is per connected body: split the mesh and require every body to be
watertight, consistently wound and of positive volume. A single leaking body is a real
defect; several disjoint sealed bodies is the intended output, and every slicer unions
them on import.

Usage: python verify.py [--verbose]
"""

import os
import re
import sys

import numpy as np
import trimesh

HERE = os.path.dirname(os.path.abspath(__file__))
STL = os.path.join(HERE, "stl")
BED = 320.0
PETG = 1.27          # g/cm3

#: Fraction of the solid volume actually extruded, at 3 walls + 10% gyroid. Two figures,
#: because the sets differ by an order of magnitude in wall-to-volume ratio: a 60 mm solid
#: cell is nearly all void, whereas a Ø5 strut is 2.4 mm of wall across a 5 mm diameter and
#: barely hollow at all. Estimates -- the slicer is the authority.
HOLLOW = 0.26        # bulky pieces: solid cells, Wigner-Seitz cells
SLENDER = 0.75       # strut-and-hub pieces: kit, fused frames, blocks
SLENDER_SETS = ("kit", "assembled", "supercell")

#: Kit filenames carry their quantity, e.g. node_A_3way_x8.stl -- so the file count is not
#: the part count, and summing file volumes understates the kit by a factor of about five.
QTY = re.compile(r"_x(\d+)\.stl$")


def bodies(m):
    """Connected components of the mesh, split by shared vertex index.

    Deliberately not `trimesh.split`: that routes through `submesh`, which calls
    `fill_holes` on any component that is not already watertight, which needs networkx --
    absent here. So a genuinely broken body made the checker crash with ModuleNotFoundError
    instead of reporting the defect. scipy's connected_components has no such dependency."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    f = np.asarray(m.faces)
    n = len(m.vertices)
    rows = np.concatenate([f[:, 0], f[:, 1], f[:, 2]])
    cols = np.concatenate([f[:, 1], f[:, 2], f[:, 0]])
    g = coo_matrix((np.ones(len(rows), np.int8), (rows, cols)), shape=(n, n))
    _, label = connected_components(g, directed=False)
    fl = label[f[:, 0]]
    return [trimesh.Trimesh(m.vertices, f[fl == k], process=False)
            for k in np.unique(fl)]


def check(path):
    m = trimesh.load(path, process=True)
    parts = bodies(m)
    issues = []
    for i, b in enumerate(parts):
        if not b.is_watertight:
            issues.append(f"body {i} not watertight")
        if not b.is_winding_consistent:
            issues.append(f"body {i} winding inconsistent")
        if b.volume <= 0.0:
            issues.append(f"body {i} volume {b.volume:.3f} <= 0")
    ext = m.extents
    if ext.max() > BED:
        issues.append(f"{ext.max():.1f} mm exceeds the {BED:.0f} mm bed")
    return {"bodies": len(parts), "tris": len(m.faces), "extents": ext,
            "volume": float(sum(b.volume for b in parts)), "issues": issues}


def main(verbose=False):
    if not os.path.isdir(STL):
        print("no stl/ directory -- run build.py first")
        return 1

    paths = []
    for root, _, files in os.walk(STL):
        paths += [os.path.join(root, f) for f in sorted(files) if f.endswith(".stl")]
    paths.sort()

    bad, tris = [], 0
    groups = {}
    for p in paths:
        r = check(p)
        m = QTY.search(os.path.basename(p))
        r["qty"] = int(m.group(1)) if m else 1
        tris += r["tris"]
        rel = os.path.relpath(p, STL).replace("\\", "/")
        groups.setdefault(rel.split("/")[0], []).append((rel, r))
        if r["issues"]:
            bad.append((rel, r["issues"]))
        if verbose:
            print(f"{rel:62} x{r['qty']:<3} {r['bodies']:3} bodies {r['tris']:7} tris "
                  f"{r['extents'][0]:6.1f} x{r['extents'][1]:6.1f} x{r['extents'][2]:6.1f}")

    print(f"{'set':16} {'files':>6} {'parts':>6} {'bodies':>7} {'triangles':>10} "
          f"{'cm3':>9} {'PETG':>8}")
    total_g = 0.0
    for g, items in sorted(groups.items()):
        parts = sum(r["qty"] for _, r in items)
        v = sum(r["volume"] * r["qty"] for _, r in items) / 1000.0
        grams = v * PETG * (SLENDER if g in SLENDER_SETS else HOLLOW)
        total_g += grams
        print(f"{g:16} {len(items):6} {parts:6} "
              f"{sum(r['bodies'] for _, r in items):7} "
              f"{sum(r['tris'] for _, r in items):10} {v:9.1f} {grams:7.0f} g")

    print(f"\n{len(paths)} STLs, {tris} triangles, "
          f"{sum(r['qty'] for items in groups.values() for _, r in items)} printed parts")
    print(f"PETG total about {total_g:.0f} g at 3 walls + 10% infill "
          f"(x{SLENDER:g} for strut sets, x{HOLLOW:g} for bulky ones)")

    if bad:
        print(f"\nFAILED ({len(bad)}):")
        for rel, issues in bad:
            print(f"  {rel}: {'; '.join(issues)}")
        return 1
    print("\nevery connected body is watertight, consistently wound and positive-volume")
    return 0


if __name__ == "__main__":
    sys.exit(main("--verbose" in sys.argv))
