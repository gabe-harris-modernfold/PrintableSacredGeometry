#!/usr/bin/env python3
"""Negative Space, icosahedral: the desk ribbon (negative_space_desk.py) wrapped
around an invisible icosahedron instead of an invisible sphere.

Two opposite vertices sit on the axis, where the sphere's poles were, so the rod
runs vertex to vertex and the existing rod + eye (negative_space_desk_rod.stl) fits
unchanged. The same 4.12-turn path is projected radially from the centre onto
the icosahedron: the ribbon's inner face lies flat on the faces and creases at
every edge it crosses. Its width is re-solved along the faces, so the band is
still as wide (across the surface) as the temperature says -- a plain radial
projection would stretch it 0.8-1.3x depending on where on a face it falls.
Thickness is normal to the faces: the outer face lies on the same icosahedron
grown by THICK.

Writes negative_space_ico.stl (ribbon + both vertex caps, print upright with tree
supports) and negative_space_ico.png."""

import numpy as np
import pyvista as pv
from scipy.spatial import ConvexHull

import manifold3d as mf
import negative_space_desk as nd
import preview_render as pr
import sag_check

pv.OFF_SCREEN = True

R_C = nd.R                      # circumradius: vertices where the sphere's poles were
N_ACROSS = 32                   # grid cells across the band, ~1 mm, so creases stay crisp
CAP_PROUD = 0.1                 # mm the vertex caps stand off each ribbon face


def icosahedron(rc):
    """Vertices with one on +z, and the hull's outward face normals and inradius."""
    k = np.arange(5)
    ring = lambda z, phase: np.column_stack([
        2 / np.sqrt(5) * np.cos(phase + 2 * np.pi * k / 5),
        2 / np.sqrt(5) * np.sin(phase + 2 * np.pi * k / 5), np.full(5, z)])
    v = rc * np.vstack([[0, 0, 1], ring(1 / np.sqrt(5), 0.0),
                        ring(-1 / np.sqrt(5), np.pi / 5), [0, 0, -1]])
    hull = ConvexHull(v)
    return v, hull, hull.equations[:, :3], -hull.equations[0, 3]


VERTS, HULL, NORMALS, R_IN = icosahedron(R_C)


def project(d, r):
    """Unit directions d (..., 3) out to the icosahedron of inradius r."""
    return d * (r / (d @ NORMALS.T).max(axis=-1))[..., None]


def direction(lat, lon):
    return pr.sphere_pt(lat, lon, 1.0)


def band(rib):
    """The ribbon on the icosahedron as a closed shell. Each row's across-band
    offsets are solved so that distance measured over the faces from the
    centreline is +-w/2, then pushed to the inner and outer icosahedra."""
    tv = nd.band_t(rib)
    lat, lon = pr.lat_lon(tv)
    w = rib.width(tv)
    m = 129                                     # trial offsets per row
    d_max = 0.5 * w / (nd.R * 0.6)              # generous: the surface scale is >= 0.63
    trial = d_max[:, None] * np.linspace(-1.0, 1.0, m)[None, :]
    pts = project(direction(lat[:, None] + trial, np.broadcast_to(lon[:, None], trial.shape)),
                  R_IN)
    step = np.linalg.norm(np.diff(pts, axis=1), axis=2)
    dist = np.concatenate([np.zeros((len(tv), 1)), np.cumsum(step, axis=1)], axis=1)
    dist -= dist[:, [m // 2]]                   # signed distance over the surface
    target = np.linspace(-0.5, 0.5, N_ACROSS + 1)[None, :] * w[:, None]
    off = np.array([np.interp(target[i], dist[i], trial[i]) for i in range(len(tv))])
    d = direction(lat[:, None] + off, np.broadcast_to(lon[:, None], off.shape))
    return nd.shell(project(d, R_IN), project(d, R_IN + nd.THICK))


def vertex_caps():
    """Caps 0.1 mm proud of the ribbon on both faces: exactly coplanar faces would
    leave slivers when the ribbon's buried end is unioned in."""
    grow = (R_IN + nd.THICK + CAP_PROUD) / R_IN
    shrink = (R_IN - CAP_PROUD) / R_IN
    surface = mf.Manifold.hull_points(VERTS * grow) - mf.Manifold.hull_points(VERTS * shrink)
    return nd.caps(surface, reach=R_C * grow + 1.0)


def edges():
    """The invisible icosahedron's 30 edges, for the explanatory panel."""
    pairs = set()
    for tri in HULL.simplices:
        for a, b in ((0, 1), (1, 2), (2, 0)):
            pairs.add(tuple(sorted((tri[a], tri[b]))))
    lines = np.hstack([[2, a, b] for a, b in sorted(pairs)])
    return pv.PolyData(VERTS, lines=lines)


def render(meshes, out):
    lo = np.min([m.bounds[0] for m in meshes], axis=0)
    box = pv.Box(bounds=(-nd.BED / 2, nd.BED / 2, -nd.BED / 2, nd.BED / 2, lo[2],
                         lo[2] + nd.BED))
    views = [("Three-quarter view, inside the 320 mm build volume", (1.0, 0.55, 0.45), "box"),
             ("The icosahedron it implies (not printed)", (1.0, 0.25, 0.2), "ghost"),
             ("Close-up: top vertex cap, ribbon creasing over the edges", None, None)]
    pl = pv.Plotter(off_screen=True, shape=(1, 3), window_size=(2700, 1100), border=False)
    for col, (title, direction_, extra) in enumerate(views):
        pl.subplot(0, col)
        pl.set_background("#1c1e23")
        for m in meshes:
            pl.add_mesh(pv.wrap(m), color="#e8dfcf", specular=0.15, ambient=0.18)
        if extra == "box":
            pl.add_mesh(box, style="wireframe", color="#6f7782", line_width=1.5)
        if extra == "ghost":
            pl.add_mesh(edges(), color="#7fb3d5", line_width=2.0, opacity=0.8)
        pl.camera.up = (0, 0, 1)
        pl.camera.view_angle = 30.0
        if direction_ is None:
            p = np.array([0.0, 0.0, R_C - 12.0])
            d = np.array([1.0, 0.5, 0.8])
            pl.camera.focal_point = tuple(p)
            pl.camera.position = tuple(p + 190.0 * d / np.linalg.norm(d))
            pl.reset_camera_clipping_range()
        else:
            d = np.asarray(direction_) / np.linalg.norm(direction_)
            pl.camera.focal_point = (0.0, 0.0, 10.0)
            pl.camera.position = tuple(np.array([0.0, 0.0, 10.0]) + 1500.0 * d)
            pl.reset_camera()
            pl.camera.zoom(1.1)
        pl.add_text(title, position="upper_left", font_size=13, color="#d6d0c4")
    pl.screenshot(out)


def main():
    rib = nd.make_ribbon()
    cap_parts, bore = vertex_caps()
    ribbon = nd.to_trimesh(mf.Manifold.batch_boolean([band(rib)] + cap_parts, mf.OpType.Add)
                           - bore)
    ribbon.export("negative_space_ico.stl")
    rod = nd.to_trimesh(nd.rod_part())

    print(f"icosahedron: circumradius {R_C:.0f} mm, inradius {R_IN:.1f} mm, edge "
          f"{np.linalg.norm(VERTS[0] - VERTS[1]):.1f} mm")
    nd.report("negative_space_ico.stl", ribbon)
    down = ribbon.face_normals[:, 2] < -np.cos(np.radians(45.0))
    off_bed = ribbon.triangles_center[:, 2] > ribbon.bounds[0, 2] + 0.5
    print(f"  {100 * ribbon.area_faces[down & off_bed].sum() / ribbon.area:.0f} % of the surface "
          f"overhangs past 45 deg (needs support)")
    lo = np.minimum(ribbon.bounds[0], rod.bounds[0])
    hi = np.maximum(ribbon.bounds[1], rod.bounds[1])
    print(f"assembled with negative_space_desk_rod.stl: {(hi - lo).round(1)} mm")

    place = lambda la, lo_, off: project(direction(la, lo_), R_IN + off)
    face_normal = lambda p: NORMALS[np.argmax(p @ NORMALS.T, axis=-1)]
    for name, kw in [("sphere", {}), ("icosahedron", dict(place=place, normal=face_normal))]:
        s = sag_check.sag(rib, nd.THICK, 2 * nd.ROD_R, **kw)
        print(f"  sag on the {name}: ribbon drops {s['sag']:.1f} mm, rod bows "
              f"{s['rod_bow']:.1f} mm, ribbon {s['mass']:.0f} g incl. rod")
    render([ribbon, rod], "negative_space_ico.png")
    print("wrote negative_space_ico.png")


if __name__ == "__main__":
    main()
