#!/usr/bin/env python3
"""Web assets for the site page web/from-one-circle.html: the matched set, gathered.

Reads the STLs the five set projects build (run their scripts first) and writes:
  site/assets/one-circle/<key>.nsm    preview meshes for the in-page 3D viewer, in the
                                      grid pose (three-fold axis up, as each piece sits on
                                      its plate); .nsm layout as in negative-space/web_assets.py
  site/assets/one-circle/<key>.webp   studio renders, the no-WebGL fallback
  site/assets/one-circle/geometry.js  window.OC: plate figures, solid corners and edges,
                                      sphere centres and every file's size, for the SVGs
  site/stl/one-circle/*.stl           the printable files, unchanged
  site/stl/one-circle-set.zip         all of them plus PRINTING.txt

Run from this directory: python web_assets.py"""

import json
import os
import shutil
import struct
import sys
import tempfile
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "lib"))

import numpy as np
import trimesh
from PIL import Image

import render_mesh as rm
import sacred_grid as sg

PROJECTS = os.path.join(HERE, "..")
SITE = os.path.join(HERE, "..", "..", "site")
ASSETS = os.path.join(SITE, "assets", "one-circle")
STL_DIR = os.path.join(SITE, "stl", "one-circle")
ZIP = os.path.join(SITE, "stl", "one-circle-set.zip")
PREVIEW_TRIS = 60_000

# key -> (project file, viewer source); viewer source "grid" regenerates the grid pose
PIECES = {f"plate-{k:02d}": (f"sacred_geometry_plates/{k:02d}_{n}.stl", "stl")
          for k, n in enumerate(sg.PLATES, 1)}
PIECES["pin"] = ("sacred_geometry_plates/plate_pin.stl", None)
PIECES.update({s: (f"platonic_solid_frames/{s}_frame.stl", "grid") for s in sg.SOLIDS})
PIECES.update({f"stellated-{s}": (f"stellated_platonic_solids/stellated_{s}_frame.stl", "grid") for s in sg.SOLIDS})
PIECES["five-tetrahedra"] = ("five_intersecting_tetrahedra/five_intersecting_tetrahedra.stl", "bodies")
PIECES["five-tetrahedra-rigid"] = ("five_intersecting_tetrahedra/five_intersecting_tetrahedra_rigid.stl", "stl")
PIECES["egg-of-life"] = ("flower_of_life_spheres/egg_of_life.stl", "grid")
PIECES["flower-of-life-spheres"] = ("flower_of_life_spheres/flower_of_life_spheres.stl", "grid")

PRINTING = """FROM ONE CIRCLE -- the matched sacred-geometry set
https://promethean-theoria.netlify.app/from-one-circle.html

Every piece is sized from one circle of radius R = 15 mm. PETG, 0.4 mm nozzle.

PLATES (01-08)  3 mm thick, 170 mm across. Print flat, no supports.
  Stack in order on three pins: print plate_pin.stl three times (it lies on its
  flat side), or use 3 mm dowels. Holes are 3.4 mm. The three holes are spaced
  unevenly (0, 120, 180 degrees), so a plate only seats one way round.

EDGE FRAMES  (Platonic solids and their stellations)  5 mm struts, 7 mm nodes.
  Each file is already seated on its widest face with a 0.6 mm flat footprint.
  Supports on (tree supports work best), brim recommended. 15% infill.

FIVE TETRAHEDRA  five_intersecting_tetrahedra.stl prints in place as five free,
  interlocked frames (0.8 mm clearance): supports only from the bed, or the
  frames fuse. five_intersecting_tetrahedra_rigid.stl is the same compound welded
  into one solid, for a sturdier print.

SPHERES  egg_of_life.stl, flower_of_life_spheres.stl sit on a 1.5 mm footprint.
  Supports under the lower spheres. Set either on the Seed or Flower plate,
  corner down, and its outline covers the plate's circles.
"""


def write_nsm(path, meshes):
    lo = np.min([m.bounds[0] for m in meshes], axis=0)
    hi = np.max([m.bounds[1] for m in meshes], axis=0)
    span = np.where(hi > lo, hi - lo, 1.0)
    out = [b"NSM1", struct.pack("<I", len(meshes)), struct.pack("<6f", *lo, *hi)]
    for m in meshes:
        q = np.round((m.vertices - lo) / span * 65535).astype("<u2")
        pos = q.tobytes()
        pos += b"\0" * (-len(pos) % 4)
        out += [struct.pack("<II", len(m.vertices), len(m.faces)), pos, m.faces.astype("<u4").tobytes()]
    data = b"".join(out)
    with open(path, "wb") as f:
        f.write(data)
    return len(data)


def grid_mesh(key):
    if key.startswith("stellated-"):
        return sg.to_trimesh(sg.grid_view("stellated_" + key[len("stellated-"):]))
    return sg.to_trimesh(sg.grid_view(key.replace("-", "_")))


def render(mesh, out, az=35, el=24, size=640):
    m = sg.decimated(mesh, 24_000)
    tmp = os.path.join(tempfile.gettempdir(), "one_circle_render.png")
    rm.render([(m.vertices, m.faces, "frame")], az, el, w=size, h=size, out=tmp)
    Image.open(tmp).save(out, "WEBP", quality=84, method=6)


def r2(a):
    return np.round(np.asarray(a, float), 2).tolist()


def main():
    os.makedirs(ASSETS, exist_ok=True)
    os.makedirs(STL_DIR, exist_ok=True)
    files = {}
    for key, (rel, view) in PIECES.items():
        src = os.path.join(PROJECTS, rel)
        stl = trimesh.load(src)
        name = os.path.basename(rel)
        shutil.copyfile(src, os.path.join(STL_DIR, name))
        files[key] = {"file": name, "kb": round(os.path.getsize(src) / 1024),
                      "mm": r2(stl.extents), "g": round(stl.volume / 1000 * 1.27, 1)}
        if view is None:
            continue
        if view == "grid":
            meshes = [grid_mesh(key)]
        elif view == "bodies":
            meshes = stl.split(only_watertight=False)
        else:
            meshes = [stl]
        meshes = [sg.decimated(m, PREVIEW_TRIS) for m in meshes]
        c = np.vstack([m.vertices for m in meshes])
        mid = (c.min(axis=0) + c.max(axis=0)) / 2
        for m in meshes:
            m.vertices = m.vertices - mid
        n = write_nsm(os.path.join(ASSETS, f"{key}.nsm"), meshes)
        files[key]["tris"] = int(sum(len(m.faces) for m in meshes))
        if not key.startswith("plate-"):
            render(trimesh.util.concatenate(meshes), os.path.join(ASSETS, f"{key}.webp"))
        print(f"{key:24s} {name:42s} {files[key]['kb']:6d} KB stl  {n // 1024:5d} KB nsm  {files[key]['mm']}")

    geom = {
        "R": sg.R, "rays": list(sg.RAYS), "stroke": sg.STROKE, "frameIn": sg.FRAME_IN, "frameOut": sg.FRAME_OUT,
        "pins": [r2(sg.ray(a, sg.PIN_AT)) for a in sg.PIN_ANGLES], "pinHole": sg.PIN_HOLE,
        "plates": {}, "solids": {}, "stellations": {}, "spheres": {}, "files": files,
    }
    for name in sg.PLATES:
        f = sg.figure(name)
        geom["plates"][name] = {"circles": r2(f["circles"]), "lines": [[r2(p), r2(q), w] for p, q, w in f["lines"]],
                                "discs": r2(f["discs"]) if f["discs"] else [],
                                "spokes": [[r2(p), r2(q)] for p, q in sg.spokes(name)]}
    for name in sg.SOLIDS + ("five_tetrahedra",):
        V, E = sg.solid(name)
        geom["solids"][name] = {"V": r2(V), "E": [list(map(int, e)) for e in E]}
    for name in sg.SOLIDS:
        V, E, true = sg.stellation(name)
        geom["stellations"][name] = {"V": r2(V), "E": [list(map(int, e)) for e in E], "true": true,
                                     "label": sg.STELLATED[name]}
    for name, n in (("egg_of_life", 2), ("flower_of_life_spheres", 3)):
        geom["spheres"][name] = r2((sg.Q @ sg.sphere_centres(n).T).T)
    with open(os.path.join(ASSETS, "geometry.js"), "w", encoding="utf-8") as f:
        f.write("window.OC=" + json.dumps(geom, separators=(",", ":")) + ";\n")
    print("wrote geometry.js", os.path.getsize(os.path.join(ASSETS, "geometry.js")) // 1024, "KB")

    with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED) as z:
        for key, info in files.items():
            z.write(os.path.join(STL_DIR, info["file"]), f"one-circle-set/{info['file']}")
        z.writestr("one-circle-set/PRINTING.txt", PRINTING)
    print(f"wrote {os.path.relpath(ZIP, HERE)}  {os.path.getsize(ZIP) / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
