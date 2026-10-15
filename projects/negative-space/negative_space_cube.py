#!/usr/bin/env python3
"""Negative Space, tilted cube: the desk ribbon (negative_space_desk.py) wrapped
around an invisible cube balanced on one corner, its body diagonal vertical.

Two opposite corners sit on the axis where the sphere's poles were, so the desk
rod and eye (negative_space_desk_rod.stl) fit unchanged. As on the pyramid
(negative_space_pyramid.py) the ribbon is mapped by height: each latitude of the
sphere is a level of the cube, linear from the top corner to the bottom one, so
every turn is a level triangle near the corners and a level hexagon round the
middle, and never rises again. The cube's six other corners split the diagonal
into thirds, and so the year: Jan - Apr wraps the three upper faces only, May - Aug
all six, Sep - Dec the three lower ones. Every face leans 35.3 deg off vertical,
so no surface of the band is a ceiling.

Thickness is normal to the faces; how a crease is filled depends on the edge. The
three edges from each axis corner lie in the cube's mirror planes, so there each
face carries its own level lines on to the bisector plane and its mirror image
meets it: a true mitre. The six zig-zag middle edges are not mirrors -- the band
crosses them at 45 deg and runs straight on over them unfolded, and for that no
sharp mitre exists (the two faces' band edges would reach the outer corner line
12 mm apart). There the fold is filled square over the stretch the band covers,
as a thick strip bent over the edge would be, and so are middle corners the band
covers. The caps are the two axis corners of the shell itself, cut off along a
level line 20 mm in, so the level band runs straight out of them; the bottom one
is cut flat 4 mm up, for a first layer.

  negative_space_cube.stl   ribbon + both axis corners, printed corner-down as it
                            hangs (brim, a few tree supports for the tall
                            narrow start); takes negative_space_desk_rod.stl
Also writes negative_space_cube.png."""

import math

import manifold3d as mf
import numpy as np
import pyvista as pv
import trimesh
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union

import negative_space_desk as nd
import negative_space_pyramid as npyr
import preview_render as pr
import sag_check

pv.OFF_SCREEN = True

Z_TOP = nd.R                    # corners on the axis at +-R: the desk rod fits
SIDE = 2.0 * Z_TOP / math.sqrt(3.0)
R_IN = SIDE / 2.0               # every face R_IN from the centre
GROW = (R_IN + nd.THICK) / R_IN
SLOPE = math.atan(math.sqrt(2.0))       # every face rises 54.7 deg from horizontal
SCALE = Z_TOP / (math.pi / 2) / math.sin(SLOPE)   # mm of fall line per radian of lat
FOOT = nd.TUBE_BOTTOM - 1.0     # bottom corner cut flat here; rod end 1 mm up inside
H_CAP = 20.0                    # caps: each axis corner of the shell, cut level this far in
H_TRIM = 12.0                   # band ends this far from the axis corners, inside the caps
MIN_HOLE = 100.0                # mm2: smaller holes in a face's footprint are slivers, filled
REACH = 15.0                    # face prisms drawn this far past their edges, then trimmed
WINDOW = math.radians(85.0)     # a face's own level lines are used within this azimuth


def cube():
    """Corners and faces of the cube with its body diagonal along z. Faces 0-2 lean
    outward-up (azimuths 0, 120, 240), 3-5 outward-down (180, 300, 60)."""
    ez = np.array([1.0, 1.0, 1.0]) / math.sqrt(3.0)
    ex = np.array([2.0, -1.0, -1.0]) / math.sqrt(6.0)
    rot = np.vstack([ex, np.cross(ez, ex), ez])
    axes = np.vstack([np.eye(3), -np.eye(3)])
    normals = axes @ rot.T
    signs = np.array(np.meshgrid([-1, 1], [-1, 1], [-1, 1])).T.reshape(-1, 3)
    verts = R_IN * signs @ rot.T
    faces = []
    for a, n in zip(axes, normals):
        ring = [i for i, s in enumerate(signs) if s @ a > 0]
        o, u, v = frame(n)
        ang = [math.atan2((verts[i] - o) @ v, (verts[i] - o) @ u) for i in ring]
        faces.append([ring[k] for k in np.argsort(ang)])
    edges = []
    for f in range(6):
        for g in range(f + 1, 6):
            if abs(normals[f] @ normals[g]) < 1e-9:
                i, j = sorted(set(faces[f]) & set(faces[g]))
                edges.append((i, j, f, g))
    return verts, normals, faces, edges


def frame(n):
    a = np.array([0.0, 0.0, 1.0]) if abs(n[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    u = np.cross(a, n)
    u /= np.linalg.norm(u)
    return R_IN * n, u, np.cross(n, u)


VERTS, NORMALS, FACES, EDGES = cube()
FRAMES = [frame(n) for n in NORMALS]
AZ = [math.atan2(n[1], n[0]) for n in NORMALS]


def mirror_edge(f, g):
    """Both faces lean the same way: an edge from an axis corner, in a mirror plane."""
    return NORMALS[f][2] * NORMALS[g][2] > 0


def height(lat):
    return Z_TOP * np.asarray(lat) / (np.pi / 2)


def surf(lat, lon):
    """(lat, lon) to the cube: the level at height(lat), met at azimuth lon."""
    lat, lon = np.broadcast_arrays(lat, lon)
    z = height(lat)
    d = np.stack([np.cos(lon), np.sin(lon)], axis=-1)
    h = d @ NORMALS[:, :2].T
    lim = np.where(h > 1e-9, (R_IN - z[..., None] * NORMALS[:, 2]) / np.maximum(h, 1e-9),
                   np.inf)
    rho = np.maximum(lim.min(axis=-1), 0.0)
    return np.stack([rho * np.cos(lon), rho * np.sin(lon), z], axis=-1)


def face_plane(f, lat, lon):
    """Face f's own plane carried on past its edges: each latitude's level line on
    it, met at azimuth lon. On the face itself this is surf()."""
    z = height(lat)
    n = NORMALS[f]
    rho = (R_IN - z * n[2]) / (math.hypot(n[0], n[1]) * np.cos(lon - AZ[f]))
    return np.stack([rho * np.cos(lon), rho * np.sin(lon), z], axis=-1)


def band_t(rib):
    """Path parameters ~1 mm apart, ending H_TRIM from each axis corner."""
    lat, lon = pr.lat_lon(rib.tv)
    p = surf(lat, lon)
    on = np.nonzero(Z_TOP - np.abs(p[:, 2]) >= H_TRIM)[0]
    s = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(p, axis=0), axis=1))]
    return np.interp(np.linspace(s[on[0]], s[on[-1]], int(s[on[-1]] - s[on[0]]) + 1), s, rib.tv)


def walls(rib, tv):
    """Latitudes of the band's two edges, solved so it is w wide square to its own
    direction, measured over the faces (as negative_space_pyramid.walls)."""
    lat, lon = pr.lat_lon(tv)
    w = rib.width(tv)
    trial = (0.65 * w / SCALE)[:, None] * np.linspace(-1.0, 1.0, npyr.M_TRIAL)[None, :]
    la = np.clip(lat[:, None] + trial, -np.pi / 2, np.pi / 2)
    lo = np.broadcast_to(lon[:, None], la.shape)
    p = surf(la, lo)
    dt = 1e-6
    tan = surf(np.clip(la - np.pi * dt, -np.pi / 2, np.pi / 2),
               lo + 2 * np.pi * pr.N_TURNS * dt) - p
    tan /= np.maximum(np.linalg.norm(tan, axis=-1, keepdims=True), 1e-12)
    t_mid = tan[:, 1:] + tan[:, :-1]
    t_mid /= np.maximum(np.linalg.norm(t_mid, axis=-1, keepdims=True), 1e-12)
    across = np.linalg.norm(np.cross(np.diff(p, axis=1), t_mid), axis=-1)
    dist = np.concatenate([np.zeros((len(tv), 1)), np.cumsum(across, axis=1)], axis=1)
    dist -= dist[:, [npyr.M_TRIAL // 2]]
    off = np.array([np.interp([-0.5 * wi, 0.5 * wi], d, tr)
                    for wi, d, tr in zip(w, dist, trial)])
    return lat[:, None] + off, lon


def to_2d(points, f):
    o, u, v = FRAMES[f]
    q = np.atleast_2d(points) - o
    return np.column_stack([q @ u, q @ v])


def to_3d(xy, f):
    o, u, v = FRAMES[f]
    return o + xy[:, :1] * u + xy[:, 1:] * v


def own_band(f, lat_w, lon):
    """Face f's band along its own level lines, in its 2-D frame."""
    ok = np.cos(lon - AZ[f]) > math.cos(WINDOW)
    polys = []
    for i0, i1 in npyr.runs(ok):
        if i1 - i0 < 2:
            continue
        ring = np.vstack([face_plane(f, lat_w[i0:i1, 0], lon[i0:i1]),
                          face_plane(f, lat_w[i0:i1, 1], lon[i0:i1])[::-1]])
        poly = Polygon(to_2d(ring, f))
        polys.append(poly if poly.is_valid else poly.buffer(0))
    return unary_union(polys)


def edge_zone(f, i, j):
    """The strip of face f's plane beyond its edge i-j, square to it, REACH deep and
    OVERLAP longer than the edge at each end."""
    a, b = to_2d(VERTS[[i, j]], f)
    along = (b - a) / np.linalg.norm(b - a)
    out = np.array([-along[1], along[0]])
    if out @ (a - np.asarray(Polygon(to_2d(VERTS[FACES[f]], f)).centroid.coords)[0]) < 0:
        out = -out
    a, b = a - npyr.OVERLAP * along, b + npyr.OVERLAP * along
    return Polygon([a - npyr.OVERLAP * out, b - npyr.OVERLAP * out,
                    b + REACH * out, a + REACH * out])


def beyond_level(f, z0, up):
    """The part of face f's plane above (up) or below level z0, in its 2-D frame."""
    o, u, v = FRAMES[f]
    grad = np.array([u[2], v[2]])                # height rises along this in the frame
    at = (z0 - o[2]) * grad / (grad @ grad)      # a point on the level line
    along = np.array([-grad[1], grad[0]]) / np.linalg.norm(grad)
    side = (1.0 if up else -1.0) * grad / np.linalg.norm(grad)
    big = 1e3
    return Polygon([at - big * along, at + big * along,
                    at + big * along + big * side, at - big * along + big * side])


def footprint(f, own):
    """Where face f's prism carries plastic, in its 2-D frame: the band, and the caps
    (the face beyond a level line H_CAP from each axis corner), on the face itself
    and, past each mirror edge, along its own level lines on to the bisector. Holes
    a band edge leaves where it grazes a cap's rim are filled."""
    reach = Polygon(to_2d(VERTS[FACES[f]], f)).buffer(npyr.OVERLAP, join_style="mitre")
    for i, j, fa, fb in EDGES:
        if f in (fa, fb) and mirror_edge(fa, fb):
            reach = reach.union(edge_zone(f, i, j))
    caps = beyond_level(f, Z_TOP - H_CAP, True).union(beyond_level(f, H_CAP - Z_TOP, False))
    region = own[f].union(caps).intersection(reach)
    return unary_union([Polygon(g.exterior, [r for r in g.interiors
                                              if Polygon(r).area > MIN_HOLE])
                        for g in getattr(region, "geoms", [region])
                        if g.geom_type == "Polygon"])


def covered(own_f, i, j, f):
    """Stretches of edge i-j the band covers (from face f's band), each OVERLAP longer."""
    from shapely.geometry import LineString
    seg = to_2d(VERTS[[i, j]], f)
    d = seg[1] - seg[0]
    cover = own_f.intersection(LineString(seg))
    for g in getattr(cover, "geoms", [cover]):
        if g.geom_type == "LineString" and g.length > 1e-6:
            s = np.sort((np.asarray(g.coords)[[0, -1]] - seg[0]) @ d / (d @ d))
            s += np.array([-1.0, 1.0]) * npyr.OVERLAP / math.sqrt(d @ d)
            yield VERTS[i] + s[0] * (VERTS[j] - VERTS[i]), VERTS[i] + s[1] * (VERTS[j] - VERTS[i])


def prism(region, f):
    n, (o, u, v) = NORMALS[f], FRAMES[f]
    contours = []
    for g in getattr(region, "geoms", [region]):
        if g.geom_type != "Polygon" or g.area < 1e-6:
            continue
        contours.append(np.asarray(g.exterior.coords)[:-1])
        contours += [np.asarray(r.coords)[:-1] for r in g.interiors]
    solid = mf.Manifold.extrude(mf.CrossSection(contours, mf.FillRule.EvenOdd),
                                nd.THICK + 1.2).translate((0.0, 0.0, -0.2))
    solid = solid.transform(np.column_stack([u, v, n, o]))
    for g in range(6):                           # back to just past each mirror bisector
        if abs(n @ NORMALS[g]) < 1e-9 and mirror_edge(f, g):
            solid = solid.trim_by_plane(npyr.unit(n - NORMALS[g]), -0.03)
    return solid


def ribbon_solid(lat_w, lon):
    """Ribbon + both corner caps: face prisms (caps included, see footprint),
    middle-edge folds and middle-corner fills, clipped to the shell between the
    cube and its THICK-grown copy, bored for the rod and footed. A round cap, as on
    the sphere, would let the level band -- twice as far out at an edge as
    mid-face -- spike out of it at each edge. Also returns the band's area on the
    faces."""
    shell = mf.Manifold.hull_points(VERTS * GROW) - mf.Manifold.hull_points(VERTS)
    own = [own_band(f, lat_w, lon) for f in range(6)]
    area = sum(own[f].intersection(Polygon(to_2d(VERTS[FACES[f]], f))).area for f in range(6))
    parts = [prism(footprint(f, own), f) for f in range(6)]
    for i, j, fa, fb in EDGES:
        if not mirror_edge(fa, fb):
            parts += [npyr.fold(p0, p1, NORMALS[fa], NORMALS[fb])
                      for p0, p1 in covered(own[fa], i, j, fa)]
    for v in range(len(VERTS)):
        fs = [f for f in range(6) if v in FACES[f]]
        if abs(VERTS[v][2]) < Z_TOP - 1.0 and own[fs[0]].contains(Point(to_2d(VERTS[v], fs[0])[0])):
            parts.append(npyr.corner(VERTS[v], [NORMALS[f] for f in fs]))
    z0, z1 = -Z_TOP * GROW - 5.0, Z_TOP * GROW + 5.0
    bore = mf.Manifold.cylinder(z1 - z0, nd.ROD_R + nd.ROD_CLEAR,
                                circular_segments=nd.SEG).translate((0.0, 0.0, z0))
    solid = (shell ^ mf.Manifold.batch_boolean(parts, mf.OpType.Add)) - bore
    return solid.trim_by_plane((0.0, 0.0, 1.0), FOOT), area


def edges():
    lines = np.hstack([[2, i, j] for i, j, _, _ in EDGES])
    return pv.PolyData(VERTS, lines=lines)


def render(meshes, out):
    lo = np.min([m.bounds[0] for m in meshes], axis=0)
    hi = np.max([m.bounds[1] for m in meshes], axis=0)
    box = pv.Box(bounds=(-nd.BED / 2, nd.BED / 2, -nd.BED / 2, nd.BED / 2, lo[2],
                         lo[2] + nd.BED))
    views = [("Three-quarter view, inside the 320 mm build volume", (1.0, 0.55, 0.45), "box"),
             ("The cube it implies (not printed)", (1.0, 0.25, 0.15), "ghost"),
             ("Close-up: the top corner cap, rod and eye", (1.0, 0.4, 0.5), None)]
    pl = pv.Plotter(off_screen=True, shape=(1, 3), window_size=(2700, 1100), border=False)
    c = 0.5 * (lo + hi)
    for col, (title, direction, extra) in enumerate(views):
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
        if extra is None:                        # close-up on the top corner
            p = np.array([0.0, 0.0, Z_TOP - 15.0])
            pl.camera.focal_point = tuple(p)
            pl.camera.position = tuple(p + 170.0 * npyr.unit(np.asarray(direction, float)))
            pl.reset_camera_clipping_range()
        else:
            d = npyr.unit(np.asarray(direction, float))
            pl.camera.focal_point = tuple(c)
            pl.camera.position = tuple(c + 1500.0 * d)
            pl.reset_camera()
            pl.camera.zoom(1.1)
        pl.add_text(title, position="upper_left", font_size=13, color="#d6d0c4")
    pl.screenshot(out)


def main():
    rib = nd.make_ribbon()
    tv = band_t(rib)
    lat_w, lon = walls(rib, tv)
    solid, area = ribbon_solid(lat_w, lon)
    print(f"solid: genus {solid.genus()} (2 = the rod bore through both caps)")
    npyr.clean(nd.to_trimesh(solid)).export("negative_space_cube.stl")
    ribbon = trimesh.load("negative_space_cube.stl")
    rod = nd.to_trimesh(nd.rod_part())

    edge = np.linalg.norm(VERTS[FACES[0][0]] - VERTS[FACES[0][1]])
    print(f"cube: side {edge:.1f} mm, body diagonal {2 * Z_TOP:.0f} mm, faces "
          f"{math.degrees(SLOPE):.1f} deg from horizontal")
    lat, lon_c = pr.lat_lon(tv)
    s = np.linalg.norm(np.diff(surf(lat, lon_c), axis=0), axis=1)
    w = rib.width(tv)
    print(f"  band area on the faces {area / 100:.0f} cm2 vs width x length "
          f"{(0.5 * (w[1:] + w[:-1]) * s).sum() / 100:.0f} cm2 (caps excluded)")
    nd.report("negative_space_cube.stl", ribbon)
    print(f"  euler {ribbon.euler_number} (genus {(2 - ribbon.euler_number) // 2}); foot "
          f"{FOOT:.0f} mm, rod end {nd.TUBE_BOTTOM:.0f} mm")
    off_bed = ribbon.triangles_center[:, 2] > ribbon.bounds[0, 2] + 0.5
    for deg in (45.0, 60.0):
        down = ribbon.face_normals[:, 2] < -math.sin(math.radians(deg))
        print(f"  {100 * ribbon.area_faces[down & off_bed].sum() / ribbon.area:.1f} % of the "
              f"surface overhangs past {deg:.0f} deg from vertical")
    print(f"  rod + eye vs ribbon: {(solid ^ nd.rod_part()).volume():.2f} mm3 overlap")
    lo = np.minimum(ribbon.bounds[0], rod.bounds[0])
    hi = np.maximum(ribbon.bounds[1], rod.bounds[1])
    print(f"assembled with negative_space_desk_rod.stl: {(hi - lo).round(1)} mm")

    def place(la, lo_, off):
        p = surf(la, lo_)
        return p + off * NORMALS[np.argmax(p @ NORMALS.T, axis=-1)]

    face_normal = lambda p: NORMALS[np.argmax(p @ NORMALS.T, axis=-1)]
    date = lambda t: pr.datetime(pr.YEAR, 1, 1) + pr.timedelta(days=int(t * pr.DAYS))
    for name, kw in [("sphere", {}), ("tilted cube", dict(place=place, normal=face_normal))]:
        r = sag_check.sag(rib, nd.THICK, 2 * nd.ROD_R, **kw)
        print(f"  sag hanging, {name}: ribbon drops {r['sag']:.1f} mm at "
              f"{date(r['t_sag']):%b %d}, rod bows {r['rod_bow']:.1f} mm, peak stress "
              f"{r['sigma']:.2f} MPa, {r['mass']:.0f} g incl. rod")
    render([ribbon, rod], "negative_space_cube.png")
    print("wrote negative_space_cube.png")


if __name__ == "__main__":
    main()
