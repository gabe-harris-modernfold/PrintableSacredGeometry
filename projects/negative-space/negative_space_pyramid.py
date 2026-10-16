#!/usr/bin/env python3
"""Negative Space, pyramidal: the desk ribbon (negative_space_desk.py) wrapped
around an invisible square pyramid with the Great Pyramid's proportions -- the
Kepler triangle, slant height / half-base = phi, faces rising at 51.8 deg.

The ribbon is mapped onto it by height rather than by angle. Each latitude of the
sphere becomes a level of the pyramid, linear along the fall lines: on the four
faces a latitude is a horizontal square, on the base a smaller square around the
centre. So every turn down the faces is a level square that never rises again,
the turns stay evenly spaced, and the year divides between faces and base as the
fall line does -- in the golden ratio, 1/phi of it (Jan 1 .. Aug 14) down the
faces, the rest spiralling in across the base, which lies flat on the bed.
(A radial projection from the centre, as negative_space_ico.py uses, would sag
every turn into a swag at each face centre and leave islands to support.)

Thickness is normal to the faces: each face's share of the band is extruded
straight out of that face, the creases are filled (see ribbon_solid), and the
whole is clipped to the pyramid shell, cut square just outside the base edges.
Rays from the centre meet these big faces at up to ~70 deg, so radial thickness
would shear the band's edges by up to ~17 mm.

  negative_space_pyramid.stl      ribbon + apex and base caps, printed base down
  negative_space_pyramid_rod.stl  rod + eye, printed flat; push it down through
                                  the apex hole into the base cap and glue
Also writes negative_space_pyramid.png."""

import math

import manifold3d as mf
import numpy as np
import pyvista as pv
import trimesh
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import unary_union

import negative_space_desk as nd
import preview_render as pr
import sag_check

pv.OFF_SCREEN = True

PHI = (1.0 + math.sqrt(5.0)) / 2.0
BASE = 290.0                    # inner base side; the grown pyramid's base is ~315
A = BASE / 2.0
H = A * math.sqrt(PHI)          # Kepler triangle: height^2 = phi * half-base^2
R_IN = A * H / (A + math.hypot(A, H))   # insphere radius; the origin is its centre
Z_BASE, Z_APEX = -R_IN, H - R_IN
LAT_B = math.pi / 2 - math.pi / PHI     # latitude of the base edge (1/phi of the year)
SCALE = PHI * PHI * A / math.pi         # mm of fall line per radian of latitude
GROW = (R_IN + nd.THICK) / R_IN         # the outer pyramid, faces THICK out
PLINTH = nd.THICK               # outer shell cut square this far outside the base edges
TRIM = 6.0                      # ribbon ends this far inside the cap rim
M_TRIAL = 129                   # trial offsets per across-band row
EPS_FACE = 0.1                  # projected onto a face only where dir . n > this
LEAN = 0.02                     # base folds reach a hair past each face normal
OVERLAP = 0.05                  # face prisms reach this far past their face's edges
BASE_FACE = 4
TUBE_BOTTOM = Z_BASE - nd.THICK + 1.0   # rod end 1 mm up inside the base cap
TUBE_TOP = Z_APEX * GROW + 6.0          # rod top: where the eye starts


def pyramid():
    """Apex, then base corners at 45, 135, 225, 315 deg; four sides facing +x, +y,
    -x, -y, then the base; every face plane at R_IN from the origin."""
    th = math.atan2(H, A)                       # face slope = its normal's tilt
    az = np.radians(45.0 + 90.0 * np.arange(4))
    verts = np.vstack([[0.0, 0.0, Z_APEX],
                       np.column_stack([A * math.sqrt(2) * np.cos(az),
                                        A * math.sqrt(2) * np.sin(az), np.full(4, Z_BASE)])])
    k = np.radians(90.0 * np.arange(4))
    normals = np.vstack([np.column_stack([math.sin(th) * np.cos(k), math.sin(th) * np.sin(k),
                                          np.full(4, math.cos(th))]), [0.0, 0.0, -1.0]])
    faces = [[0, 1 + (j - 1) % 4, 1 + j] for j in range(4)] + [[1, 2, 3, 4]]
    edges = ([(0, 1 + j, j, (j + 1) % 4) for j in range(4)] +          # slant
             [(1 + (j - 1) % 4, 1 + j, j, 4) for j in range(4)])       # base
    corners = [(1 + j, (j, (j + 1) % 4, 4)) for j in range(4)]
    for f, n in zip(faces, normals):
        assert np.allclose(verts[f] @ n, R_IN)
    return verts, normals, faces, edges, corners


VERTS, NORMALS, FACES, EDGES, CORNERS = pyramid()


def frame(n):
    """A face's plane: its tangent point on the insphere and two in-plane axes."""
    a = np.array([0.0, 0.0, 1.0]) if abs(n[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    u = np.cross(a, n)
    u /= np.linalg.norm(u)
    return R_IN * n, u, np.cross(n, u)


FRAMES = [frame(n) for n in NORMALS]


def unit(v):
    return v / np.linalg.norm(v)


def surf(lat, lon):
    """(lat, lon) to the pyramid: down the fall line from the apex to the base edge
    while lat > LAT_B, then in across the base to its centre, both linear in lat."""
    lat, lon = np.broadcast_arrays(lat, lon)
    c, s = np.cos(lon), np.sin(lon)
    edge = A / np.maximum(np.abs(c), np.abs(s))         # base square at this azimuth
    u = (np.pi / 2 - lat) / (np.pi / 2 - LAT_B)
    rho = np.where(u <= 1.0, u, (np.pi / 2 + lat) / (np.pi / 2 + LAT_B)) * edge
    z = np.where(u <= 1.0, Z_APEX + u * (Z_BASE - Z_APEX), Z_BASE)
    return np.stack([rho * c, rho * s, z], axis=-1)


def band_t(rib):
    """Path parameters ~1 mm apart, trimmed TRIM inside each cap's rim."""
    lat, lon = pr.lat_lon(rib.tv)
    p = surf(lat, lon)
    on = np.nonzero(np.hypot(p[:, 0], p[:, 1]) >= nd.CAP_R - TRIM)[0]
    s = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(p, axis=0), axis=1))]
    s0, s1 = s[on[0]], s[on[-1]]
    return np.interp(np.linspace(s0, s1, int(s1 - s0) + 1), s, rib.tv)


def walls(rib, tv):
    """The band's two edges on the pyramid. Each row runs down its meridian; the
    offsets are solved so the band is w wide square to its own direction, measured
    over the faces (rows meet the band up to ~30 deg off square near the corners).
    Returns both edges' points and their latitudes (rows, 2)."""
    lat, lon = pr.lat_lon(tv)
    w = rib.width(tv)
    trial = (0.65 * w / SCALE)[:, None] * np.linspace(-1.0, 1.0, M_TRIAL)[None, :]
    la = np.clip(lat[:, None] + trial, -np.pi / 2, np.pi / 2)
    lo = np.broadcast_to(lon[:, None], la.shape)
    p = surf(la, lo)
    dt = 1e-6                                    # the band's direction at each trial point
    tan = surf(np.clip(la - np.pi * dt, -np.pi / 2, np.pi / 2),
               lo + 2 * np.pi * pr.N_TURNS * dt) - p
    tan /= np.maximum(np.linalg.norm(tan, axis=-1, keepdims=True), 1e-12)
    t_mid = tan[:, 1:] + tan[:, :-1]
    t_mid /= np.maximum(np.linalg.norm(t_mid, axis=-1, keepdims=True), 1e-12)
    across = np.linalg.norm(np.cross(np.diff(p, axis=1), t_mid), axis=-1)
    dist = np.concatenate([np.zeros((len(tv), 1)), np.cumsum(across, axis=1)], axis=1)
    dist -= dist[:, [M_TRIAL // 2]]
    off = np.array([np.interp([-0.5 * wi, 0.5 * wi], d, tr)
                    for wi, d, tr in zip(w, dist, trial)])
    lat_w = lat[:, None] + off
    return surf(lat_w[:, 0], lon), surf(lat_w[:, 1], lon), lat_w


def runs(mask):
    d = np.diff(np.r_[0, mask.astype(int), 0])
    return zip(np.nonzero(d == 1)[0], np.nonzero(d == -1)[0])


def strip(w0, w1, f):
    """The band projected from the centre onto face f's plane, in its 2-D frame,
    wherever both edges face it. Central projection is one-to-one there, so the
    strip stays a simple polygon; only the part inside the face is used."""
    n, (o, u, v) = NORMALS[f], FRAMES[f]
    ok = np.ones(len(w0), bool)
    for w in (w0, w1):
        ok &= w @ n > EPS_FACE * np.linalg.norm(w, axis=1)
    polys = []
    for i0, i1 in runs(ok):
        if i1 - i0 < 2:
            continue
        ring = np.vstack([w0[i0:i1], w1[i0:i1][::-1]])
        q = ring * (R_IN / (ring @ n))[:, None] - o
        poly = Polygon(np.column_stack([q @ u, q @ v]))
        polys.append(poly if poly.is_valid else poly.buffer(0))
    return unary_union(polys)


def to_2d(points, f):
    o, u, v = FRAMES[f]
    q = np.atleast_2d(points) - o
    return np.column_stack([q @ u, q @ v])


def prism(region, f):
    """Region (2-D, in face f's frame) pushed straight out of the face, from just
    inside it to past the outer face; the shell clip trims both ends."""
    n, (o, u, v) = NORMALS[f], FRAMES[f]
    contours = []
    for g in getattr(region, "geoms", [region]):
        if g.geom_type != "Polygon" or g.area < 1e-6:
            continue
        contours.append(np.asarray(g.exterior.coords)[:-1])
        contours += [np.asarray(r.coords)[:-1] for r in g.interiors]
    if not contours:
        return None
    solid = mf.Manifold.extrude(mf.CrossSection(contours, mf.FillRule.EvenOdd),
                                nd.THICK + 1.2).translate((0.0, 0.0, -0.2))
    return solid.transform(np.column_stack([u, v, n, o]))


def side_plane(k, lat, lon):
    """Side face k's own plane carried on past its edges: each latitude's level line
    on it, met at azimuth lon. On the face itself this is surf()."""
    th = math.atan2(H, A)
    z = Z_APEX + (np.pi / 2 - lat) / (np.pi / 2 - LAT_B) * (Z_BASE - Z_APEX)
    rho = (R_IN - math.cos(th) * z) / (math.sin(th) * np.cos(lon - math.radians(90.0 * k)))
    return np.stack([rho * np.cos(lon), rho * np.sin(lon), z], axis=-1)


def side_strip(k, lat_w, lon):
    """Side face k's band in its 2-D frame, carried past both slant edges along the
    face's own level lines (within 70 deg of the face's azimuth)."""
    ok = np.cos(lon - math.radians(90.0 * k)) > math.cos(math.radians(70.0))
    polys = []
    for i0, i1 in runs(ok):
        if i1 - i0 < 2:
            continue
        ring = np.vstack([side_plane(k, lat_w[i0:i1, 0], lon[i0:i1]),
                          side_plane(k, lat_w[i0:i1, 1], lon[i0:i1])[::-1]])
        poly = Polygon(to_2d(ring, k))
        polys.append(poly if poly.is_valid else poly.buffer(0))
    return unary_union(polys)


def side_region(k):
    """Where side face k's prism may reach in its plane: well past the slant edges
    (the bisector trims it back), only OVERLAP past the base edge."""
    tri = to_2d(VERTS[FACES[k]], k)
    apex, c0, c1 = tri
    out = c0 - apex + c1 - apex                 # away from the apex, square to the base
    out -= (out @ unit(c1 - c0)) * unit(c1 - c0)
    out = unit(out)
    keep = Polygon(np.vstack([c0 + OVERLAP * out - 1e3 * unit(c1 - c0),
                              c1 + OVERLAP * out + 1e3 * unit(c1 - c0),
                              apex - 1e3 * out + 1e3 * unit(c1 - c0),
                              apex - 1e3 * out - 1e3 * unit(c1 - c0)]))
    return Polygon(tri).buffer(15.0, join_style="mitre").intersection(keep)


def fold(p0, p1, na, nb):
    """A base fold where the band covers edge p0..p1: everything whose nearest point
    on the pyramid lies on that stretch of the edge."""
    m = unit(na + nb)
    a, b = unit(na + LEAN * (na - nb)), unit(nb + LEAN * (nb - na))
    t = nd.THICK
    pts = [p + d for p in (p0, p1) for d in (-0.2 * m, 6 * t * a, 6 * t * b, 3 * t * m)]
    return mf.Manifold.hull_points(np.array(pts))


def covered(strip_2d, p0, p1, f):
    """Stretches of edge p0..p1 the band covers, from its strip in face f's frame,
    each reaching OVERLAP past the band's edges."""
    seg = to_2d(np.array([p0, p1]), f)
    d = seg[1] - seg[0]
    cover = strip_2d.intersection(LineString(seg))
    for g in getattr(cover, "geoms", [cover]):
        if g.geom_type == "LineString" and g.length > 1e-6:
            s = np.sort((np.asarray(g.coords)[[0, -1]] - seg[0]) @ d / (d @ d))
            s += np.array([-1.0, 1.0]) * OVERLAP / math.sqrt(d @ d)
            yield p0 + s[0] * (p1 - p0), p0 + s[1] * (p1 - p0)


def corner(p, normals):
    """Everything whose nearest point on the pyramid is corner p: all of it when
    the band runs over the corner."""
    m = unit(sum(normals))
    out = np.linalg.solve(np.array(normals), np.ones(3))   # to the grown corner
    t = nd.THICK
    pts = [p - 0.2 * m, p + 1.5 * t * out]
    pts += [p + 6 * t * unit(n + LEAN * (n - m)) for n in normals]
    return mf.Manifold.hull_points(np.array(pts))


def ribbon_solid(w0, w1, lat_w, lon):
    """Ribbon + caps, clipped to the shell between the pyramid and its THICK-grown
    copy (cut square round the base) and bored for the rod. Each face's band is
    pushed straight out of the face. At a slant edge each side face carries its
    own level lines on to the bisector plane, where its mirror-image neighbour
    meets it: a true mitre for a level band. (Mitring with the neighbour unfolded
    flat instead would bend every band edge 63 deg at the crease -- the pyramid's
    angle deficit -- and stick a spur out at each.) The band crosses the base
    edges almost lengthwise (~5 deg), where a mitre would hang flaps under it
    wherever it runs close by, so a base fold is filled only where the band covers
    the edge, as a thick strip folded square over it would be; corners the band
    covers are filled whole. Also returns the band's area on the faces."""
    # The grown pyramid's base edges are 52 deg knife edges 12 mm out on the bed and
    # its base corners 18 mm spikes; cut square at PLINTH, so every fold keeps
    # THICK of plastic round it and stands on a short vertical wall.
    side = 2.0 * (A + PLINTH)
    plinth = mf.Manifold.cube((side, side, 4.0 * H), center=True)
    shell = (mf.Manifold.hull_points(VERTS * GROW) ^ plinth) - mf.Manifold.hull_points(VERTS)
    strips = [strip(w0, w1, f) for f in range(len(FACES))]
    footprint = 0.0
    parts = []
    for f, face in enumerate(FACES):
        outline = Polygon(to_2d(VERTS[face], f))
        footprint += strips[f].intersection(outline).area
        if f == BASE_FACE:
            # Past the edge by OVERLAP into the fold: a prism side through the
            # pyramid's own edge would meet four surfaces on one line (slivers).
            parts.append(prism(strips[f].intersection(
                outline.buffer(OVERLAP, join_style="mitre")), f))
            continue
        solid = prism(side_strip(f, lat_w, lon).intersection(side_region(f)), f)
        for g in ((f + 1) % 4, (f - 1) % 4):    # back to just past each slant bisector
            solid = solid.trim_by_plane(unit(NORMALS[f] - NORMALS[g]), -0.03)
        parts.append(solid)
    for i, j, fa, fb in EDGES:
        if fb == BASE_FACE:
            parts += [fold(p0, p1, NORMALS[fa], NORMALS[fb])
                      for p0, p1 in covered(strips[fa], VERTS[i], VERTS[j], fa)]
    for v, fs in CORNERS:
        if strips[fs[-1]].contains(Point(to_2d(VERTS[v], fs[-1])[0])):
            parts.append(corner(VERTS[v], [NORMALS[f] for f in fs]))
    z0, z1 = Z_BASE - nd.THICK - 5.0, Z_APEX * GROW + 5.0
    parts.append(mf.Manifold.cylinder(z1 - z0, nd.CAP_R, circular_segments=nd.SEG)
                 .translate((0.0, 0.0, z0)))
    bore = mf.Manifold.cylinder(z1 - z0, nd.ROD_R + nd.ROD_CLEAR,
                                circular_segments=nd.SEG).translate((0.0, 0.0, z0))
    return (shell ^ mf.Manifold.batch_boolean(parts, mf.OpType.Add)) - bore, footprint


def clean(mesh):
    """Collapse the sub-micron slivers the booleans leave where two faces' copies
    of one band edge nearly coincide, so the STL reloads watertight."""
    mesh = mesh.copy()
    mesh.merge_vertices(digits_vertex=4)
    mesh.update_faces(mesh.nondegenerate_faces())
    mesh.update_faces(mesh.unique_faces())
    return mesh


def edges():
    """The invisible pyramid's eight edges, for the explanatory panel."""
    lines = np.hstack([[2, i, j] for i, j, _, _ in EDGES])
    return pv.PolyData(VERTS, lines=lines)


def render(meshes, out):
    lo = np.min([m.bounds[0] for m in meshes], axis=0)
    hi = np.max([m.bounds[1] for m in meshes], axis=0)
    box = pv.Box(bounds=(-nd.BED / 2, nd.BED / 2, -nd.BED / 2, nd.BED / 2, lo[2],
                         lo[2] + nd.BED))
    views = [("Three-quarter view, inside the 320 mm build volume", (1.0, 0.55, 0.45), "box"),
             ("The pyramid it implies (not printed)", (1.0, 0.3, 0.12), "ghost"),
             ("From below: the base spiral and the folds over the base edges",
              (1.0, 0.45, -0.8), None)]
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
        d = np.asarray(direction) / np.linalg.norm(direction)
        pl.camera.focal_point = tuple(c)
        pl.camera.position = tuple(c + 1500.0 * d)
        pl.camera.up = (0, 0, 1)
        pl.camera.view_angle = 30.0
        pl.reset_camera()
        pl.camera.zoom(1.1)
        pl.add_text(title, position="upper_left", font_size=13, color="#d6d0c4")
    pl.screenshot(out)


def main():
    rib = nd.make_ribbon()
    tv = band_t(rib)
    w0, w1, lat_w = walls(rib, tv)
    solid, footprint = ribbon_solid(w0, w1, lat_w, pr.lat_lon(tv)[1])
    print(f"solid: genus {solid.genus()} (2 = the rod bore through both caps)")
    clean(nd.to_trimesh(solid)).export("negative_space_pyramid.stl")
    rod = nd.rod_part(TUBE_BOTTOM, TUBE_TOP)
    rod_flat = nd.lay_flat(rod)
    rod_flat.export("negative_space_pyramid_rod.stl")
    rod_in_place = nd.to_trimesh(rod)

    day_b = pr.datetime(pr.YEAR, 1, 1) + pr.timedelta(days=int(pr.DAYS / PHI))
    edge = math.hypot(A * math.sqrt(2), H)
    print(f"pyramid: base {BASE:.0f} mm, height {H:.1f} mm, faces at "
          f"{math.degrees(math.atan2(H, A)):.2f} deg, slant edges {edge:.1f} mm, insphere "
          f"{R_IN:.1f} mm; base edge reached {day_b:%b %d}")
    lat, lon = pr.lat_lon(tv)
    s = np.linalg.norm(np.diff(surf(lat, lon), axis=0), axis=1)
    w = rib.width(tv)
    print(f"  band footprint {footprint / 100:.0f} cm2 vs width x length "
          f"{(0.5 * (w[1:] + w[:-1]) * s).sum() / 100:.0f} cm2 (caps excluded)")
    ribbon = trimesh.load("negative_space_pyramid.stl")       # check what was written
    nd.report("negative_space_pyramid.stl", ribbon)
    off_bed = ribbon.triangles_center[:, 2] > ribbon.bounds[0, 2] + 0.5
    for deg in (45.0, 55.0):                     # overhang angle from vertical
        down = ribbon.face_normals[:, 2] < -math.sin(math.radians(deg))
        print(f"  {100 * ribbon.area_faces[down & off_bed].sum() / ribbon.area:.1f} % of the "
              f"surface overhangs past {deg:.0f} deg from vertical")
    nd.report("negative_space_pyramid_rod.stl", rod_flat)
    lo = np.minimum(ribbon.bounds[0], rod_in_place.bounds[0])
    hi = np.maximum(ribbon.bounds[1], rod_in_place.bounds[1])
    print(f"assembled: {(hi - lo).round(1)} mm")

    def place(la, lo_, off):
        p = surf(la, lo_)
        return p + off * NORMALS[np.argmax(p @ NORMALS.T, axis=-1)]

    face_normal = lambda p: NORMALS[np.argmax(p @ NORMALS.T, axis=-1)]
    cases = [("sphere, hanging", {}),
             ("pyramid, hanging", dict(place=place, normal=face_normal, top=TUBE_TOP)),
             ("pyramid, standing on its base", dict(place=place, normal=face_normal,
                                                     top=TUBE_TOP, stand=Z_BASE - 1.0))]
    date = lambda t: pr.datetime(pr.YEAR, 1, 1) + pr.timedelta(days=int(t * pr.DAYS))
    for name, kw in cases:
        r = sag_check.sag(rib, nd.THICK, 2 * nd.ROD_R, **kw)
        print(f"  sag, {name}: ribbon drops {r['sag']:.1f} mm at {date(r['t_sag']):%b %d}, "
              f"rod bows {r['rod_bow']:.1f} mm, peak stress {r['sigma']:.2f} MPa, "
              f"{r['mass']:.0f} g incl. rod")
    render([ribbon, rod_in_place], "negative_space_pyramid.png")
    print("wrote negative_space_pyramid.png")


if __name__ == "__main__":
    main()
