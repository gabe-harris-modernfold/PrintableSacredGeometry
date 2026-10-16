#!/usr/bin/env python3
"""Desk-size Negative Space: the "bold" design scaled to the 320 x 320 x 320 mm bed,
all PETG, two prints. Nothing metal, no spokes: the ribbon is held only where it
flares into a small spherical cap at each pole, and a thin rod runs between the
caps up to the hanging eye.

  negative_space_desk.stl      ribbon + both caps, printed upright with tree supports
  negative_space_desk_rod.stl  rod + eye, printed lying flat so its layers run along
                               it; push it down through both cap holes and glue

The rod cannot print in place: with no spokes it would stand ~290 mm free between
the caps, far too slender to survive the nozzle, and bending would load its layers
in peel.

The ribbon is uncut -- no moon holes, day notches or numerals -- and only its width
still carries data (the temperature). Path, width and any cuts come from
preview_render.Ribbon, with the sphere radius and the printing minimums overridden
for this scale; switching APERTURE / NOTCHES / NUMBERS back on restores the cuts as
boolean subtractions (manifold3d). sag_check.py models the same ribbon
(make_ribbon) under its own weight.

Also writes negative_space_desk.png (the assembly)."""

import math

import manifold3d as mf
import numpy as np
import pyvista as pv
import trimesh

import preview_render as pr

pv.OFF_SCREEN = True

BED = 320.0
R = 138.0                       # sphere radius; total height lands at ~312 mm
TURNS = pr.DAYS / (3 * pr.SYNODIC)   # one turn per three lunations, as "bold"
W = (14.0, 32.0)                # band width at the Jan / Jul normal
APERTURE = 0.0                  # moon hole / band width; 0 = no moon holes (was 0.72)
THICK = 6.0                     # set by sag_check.py: self-weight sag vs plastic
NUMBERS, NOTCHES = None, False
PRINT_MINIMUMS = dict(LIGAMENT=1.5, NOTCH_W=1.2, NOTCH_D=2.5, MONTH_NOTCH_W=3.0,
                      END_CLEAR=0.8, MIN_FEATURE=1.2, NUM_GAP=1.5, MIN_MOON=4.0)
ROD_R = 2.5                     # 5 mm rod: sag_check.py --rod 8 5 4 3 for the trade-off
ROD_CLEAR = 0.2                 # radial clearance in the cap holes, taken up by glue
ROD_FLAT = 0.3                  # trimmed off the rod's underside so it sits flat to print
CAP_R = 15.0                    # pole caps: spherical patches this far out from the axis
ROOT = (CAP_R, 30.0, 10.0)      # band flares to 10 mm over the last 30 mm into each cap
TUBE_BOTTOM, TUBE_TOP = -R - 5.0, R + 12.0    # rod ends (TUBE_TOP: where the eye starts)
EYE_R, EYE_HOLE_R, EYE_THICK = 12.0, 5.5, 2 * ROD_R   # eye flat in the rod's print plane
SEG = 48                        # circular segments


def to_manifold(verts, tris):
    return mf.Manifold(mf.Mesh(vert_properties=np.asarray(verts, np.float32),
                               tri_verts=np.asarray(tris, np.uint32)))


def band_t(rib):
    """Path parameters along the ribbon, trimmed 5 mm inside each pole cap's rim so
    its ends are buried in the cap (and never fold over the pole), ~1 mm apart."""
    lat_all, _ = pr.lat_lon(rib.tv)
    on = np.nonzero(R * np.cos(lat_all) >= CAP_R - 5.0)[0]
    t0, t1 = rib.tv[on[0]], rib.tv[on[-1]]
    n_a = int((np.interp(t1, rib.tv, rib.s_v) - np.interp(t0, rib.tv, rib.s_v)) / 1.0)
    return np.linspace(t0, t1, n_a + 1)


def band(rib):
    """The uncut ribbon on the sphere as a closed shell."""
    tv = band_t(rib)
    vv = np.linspace(-0.5, 0.5, 7)
    lat, lon = pr.lat_lon(tv)
    lat_g = lat[:, None] + vv[None, :] * rib.width(tv)[:, None] / R
    inner = pr.sphere_pt(lat_g, np.broadcast_to(lon[:, None], lat_g.shape))
    return shell(inner, inner * (R + THICK) / R)


def shell(inner, outer):
    """Closed solid between two (along, across, 3) grids of matching points: the
    inner and outer faces, both long edges and both end caps."""
    n_a, n_c = inner.shape[0] - 1, inner.shape[1] - 1
    verts = np.vstack([inner.reshape(-1, 3), outer.reshape(-1, 3)])
    n = (n_a + 1) * (n_c + 1)
    idx = np.arange(n).reshape(n_a + 1, n_c + 1)
    I, O = idx, idx + n
    quads = []
    a, b, c, d = I[:-1, :-1], I[1:, :-1], I[1:, 1:], I[:-1, 1:]
    quads += [np.stack([a, d, c, b], -1).reshape(-1, 4)]                 # inner, faces in
    quads += [np.stack([a + n, b + n, c + n, d + n], -1).reshape(-1, 4)]  # outer
    quads += [np.stack([I[:-1, 0], I[1:, 0], O[1:, 0], O[:-1, 0]], -1)]   # south edge
    quads += [np.stack([I[:-1, -1], O[:-1, -1], O[1:, -1], I[1:, -1]], -1)]  # north edge
    quads += [np.stack([I[0, :-1], O[0, :-1], O[0, 1:], I[0, 1:]], -1)]   # start cap
    quads += [np.stack([I[-1, :-1], I[-1, 1:], O[-1, 1:], O[-1, :-1]], -1)]  # end cap
    q = np.vstack(quads)
    tris = np.vstack([q[:, [0, 1, 2]], q[:, [0, 2, 3]]])
    return to_manifold(verts, tris)


def frame(t):
    """Origin on the sphere at path parameter t, unit path tangent T, the
    perpendicular Y toward the January pole, and the outward normal n."""
    lat, lon = pr.lat_lon(np.array([t, t + 1e-6]))
    p = pr.sphere_pt(lat, lon)
    n = p[0] / R
    tan = p[1] - p[0]
    tan -= tan.dot(n) * n
    tan /= np.linalg.norm(tan)
    return p[0], tan, np.cross(n, tan), n


def cutter(poly, t):
    """A 2-D outline in the band's local (along, across) mm, pushed radially
    through the full ribbon thickness at path parameter t."""
    poly = np.asarray(poly, float)
    x, y = poly[:, 0], poly[:, 1]
    if 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y) < 0:
        poly = poly[::-1]
    solid = mf.Manifold.extrude(mf.CrossSection([poly]), THICK + 4.0)
    p, tan, y_ax, n = frame(t)
    base = p - 2.0 * n
    return solid.transform(np.column_stack([tan, y_ax, n, base]))


def rect(x0, x1, y0, y1):
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]


def moon_outline(r, phase, full, k=24):
    th = np.linspace(-0.5 * math.pi, 0.5 * math.pi, k + 1)
    if full:
        th = np.linspace(0.0, 2.0 * math.pi, 2 * k, endpoint=False)
        return np.column_stack([r * np.cos(th), r * np.sin(th)])
    side = 1.0 if phase < math.pi else -1.0      # waxing lit on +x, as seen from outside
    y = r * np.sin(th)
    outer = np.column_stack([side * np.sqrt(r * r - y * y), y])
    term = np.column_stack([side * math.cos(phase) * np.sqrt(r * r - y * y), y])[::-1][1:-1]
    return np.vstack([outer, term])


def t_at(rib, s):
    return float(np.interp(s, rib.s_v, rib.tv))


def cutters(rib):
    cuts = []
    for d in range(pr.DAYS):                     # moons
        r = rib.r_moon[d]
        if r < 1.0 or rib.moon_none[d]:
            continue
        t = t_at(rib, 0.5 * (rib.s_day[d] + rib.s_day[d + 1]))
        out = moon_outline(r, rib.phase[d], rib.moon_full[d])
        cuts.append(cutter(out + [0.0, rib.y_moon[d]], t))
    for b in range(1, pr.DAYS if rib.notches else 1):   # day / month notches, upper edge
        t = b / pr.DAYS
        w = float(rib.width(np.array(t)))
        if w < 2.0:
            continue
        nw = pr.MONTH_NOTCH_W if b in rib.month_starts else pr.NOTCH_W
        depth = min(pr.NOTCH_D, pr.NOTCH_MAX * w)
        cuts.append(cutter(rect(-nw / 2, nw / 2, w / 2 - depth, w / 2 + 3.0), t))
    if not NUMBERS:
        return cuts
    h, s, g = NUMBERS["h"], NUMBERS["stroke"], NUMBERS["bridge"]
    dw, sp = pr.DIGIT_ASPECT * h, pr.DIGIT_SPACE * h
    seg_box = {"a": (s + g, dw - s - g, h - s, h), "b": (dw - s, dw, h / 2 + g / 2, h),
               "c": (dw - s, dw, 0, h / 2 - g / 2), "d": (s + g, dw - s - g, 0, s),
               "e": (0, s, 0, h / 2 - g / 2), "f": (0, s, h / 2 + g / 2, h),
               "g": (s + g, dw - s - g, h / 2 - s / 2, h / 2 + s / 2)}
    for left, nw, dom in rib.labels:             # stencil day numbers on the lower lane
        t = t_at(rib, left + 0.5 * nw)
        y0 = -0.5 * float(rib.width(np.array(t))) + pr.LIGAMENT
        digits = [int(dom) // 10, int(dom) % 10] if dom >= 10 else [int(dom)]
        for i, dig in enumerate(digits):
            u0 = -0.5 * nw + i * (dw + sp)
            for name, on in zip("abcdefg", pr.SEG_TABLE[dig]):
                if on:
                    x0, x1, v0, v1 = seg_box[name]
                    cuts.append(cutter(rect(u0 + x0, u0 + x1, y0 + v0, y0 + v1), t))
    return cuts


def caps(surface=None, reach=R + THICK + 1.0):
    """Pole caps -- the ribbon's surface layer (by default the sphere shell R .. R +
    THICK) inside radius CAP_R of the axis -- and the bore through both for the rod."""
    if surface is None:
        surface = mf.Manifold.sphere(R + THICK, 512) - mf.Manifold.sphere(R, 512)
    out = [surface ^ mf.Manifold.cylinder(reach, CAP_R, circular_segments=SEG)
           .translate((0.0, 0.0, z0)) for z0 in (0.0, -reach)]
    bore = mf.Manifold.cylinder(2 * reach, ROD_R + ROD_CLEAR, circular_segments=SEG)
    return out, bore.translate((0.0, 0.0, -reach))


def rod_part(bottom=None, top=None):
    """Rod and eye as assembled (rod up the axis, eye in the XZ plane), with a flat
    trimmed along the -Y side for printing it lying down. The rod runs from bottom
    to top (default TUBE_BOTTOM .. TUBE_TOP)."""
    bottom = TUBE_BOTTOM if bottom is None else bottom
    top = TUBE_TOP if top is None else top
    rod = mf.Manifold.cylinder(top - bottom, ROD_R, circular_segments=SEG)
    zc = top + EYE_R - 3.0
    eye = mf.Manifold.extrude(mf.CrossSection.circle(EYE_R, 96)
                              - mf.CrossSection.circle(EYE_HOLE_R, 64), EYE_THICK)
    part = (rod.translate((0.0, 0.0, bottom))
            + eye.rotate((90.0, 0.0, 0.0)).translate((0.0, 0.5 * EYE_THICK, zc)))
    keep = mf.Manifold.cube((2000.0, 2000.0, 2000.0)).translate(
        (-1000.0, -ROD_R + ROD_FLAT, -1000.0))
    return part ^ keep


def lay_flat(rod):
    """rod_part() turned onto its trimmed flat and centred on the bed."""
    flat = to_trimesh(rod.rotate((90.0, 0.0, 0.0)))
    flat.apply_translation([-flat.centroid[0], -flat.centroid[1], -flat.bounds[0, 2]])
    return flat


def render(meshes, out):
    """The assembly: three views, the first inside the build volume."""
    lo = np.min([m.bounds[0] for m in meshes], axis=0)
    hi = np.max([m.bounds[1] for m in meshes], axis=0)
    box = pv.Box(bounds=(-BED / 2, BED / 2, -BED / 2, BED / 2, lo[2], lo[2] + BED))
    views = [("Three-quarter view, inside the 320 mm build volume", (1.0, 0.55, 0.45), True),
             ("Side", (1.0, 0.0, 0.0), False),
             ("Close-up: ribbon into the top pole cap, rod through it", None, False)]
    pl = pv.Plotter(off_screen=True, shape=(1, 3), window_size=(2700, 1100), border=False)
    for col, (title, direction, show_box) in enumerate(views):
        pl.subplot(0, col)
        pl.set_background("#1c1e23")
        for m in meshes:
            pl.add_mesh(pv.wrap(m), color="#e8dfcf", specular=0.15, ambient=0.18)
        if show_box:
            pl.add_mesh(box, style="wireframe", color="#6f7782", line_width=1.5)
        if direction is None:                    # close-up on the top pole cap
            p = np.array([0.0, 0.0, R - 5.0])
            d = np.array([1.0, 0.6, 0.7])
            pl.camera.focal_point = tuple(p)
            pl.camera.position = tuple(p + 170.0 * d / np.linalg.norm(d))
            pl.camera.up = (0, 0, 1)
            pl.camera.view_angle = 30.0
            pl.reset_camera_clipping_range()
        else:
            d = np.asarray(direction) / np.linalg.norm(direction)
            c = 0.5 * (lo + hi)
            pl.camera.focal_point = tuple(c)
            pl.camera.position = tuple(c + 1500.0 * d)
            pl.camera.up = (0, 0, 1)
            pl.camera.view_angle = 30.0
            pl.reset_camera()
            pl.camera.zoom(1.1)
        pl.add_text(title, position="upper_left", font_size=13, color="#d6d0c4")
    pl.screenshot(out)


def make_ribbon():
    """preview_render's Ribbon at this scale: radius, turns and printing minimums."""
    pr.R, pr.N_TURNS = R, TURNS
    for k, v in PRINT_MINIMUMS.items():
        setattr(pr, k, v)
    return pr.Ribbon(*W, APERTURE, NUMBERS, root=ROOT, notches=NOTCHES)


def to_trimesh(m):
    mg = m.to_mesh()
    return trimesh.Trimesh(mg.vert_properties[:, :3], mg.tri_verts, process=True)


def report(name, mesh):
    bodies = len(mesh.split(only_watertight=False))
    print(f"{name}: extents {mesh.extents.round(1)} mm (bed {BED:.0f}), watertight "
          f"{mesh.is_watertight}, {bodies} bod{'y' if bodies == 1 else 'ies'}, "
          f"{mesh.volume * pr.PETG_G_MM3:.0f} g PETG")


def main():
    rib = make_ribbon()
    cuts = cutters(rib)
    ribbon = band(rib)
    if cuts:
        ribbon -= mf.Manifold.batch_boolean(cuts, mf.OpType.Add)
    cap_parts, bore = caps()
    ribbon = to_trimesh(mf.Manifold.batch_boolean([ribbon] + cap_parts, mf.OpType.Add) - bore)
    ribbon.export("negative_space_desk.stl")
    rod = rod_part()
    rod_in_place = to_trimesh(rod)
    rod_flat = lay_flat(rod)
    rod_flat.export("negative_space_desk_rod.stl")

    report("negative_space_desk.stl", ribbon)
    down = ribbon.face_normals[:, 2] < -math.cos(math.radians(45.0))
    off_bed = ribbon.triangles_center[:, 2] > ribbon.bounds[0, 2] + 0.5
    print(f"  ribbon {THICK:.1f} mm, {len(cuts)} cuts; "
          f"{100 * ribbon.area_faces[down & off_bed].sum() / ribbon.area:.0f} % of the surface "
          f"overhangs past 45 deg (needs support)")
    report("negative_space_desk_rod.stl", rod_flat)
    print(f"  {2 * ROD_R:.0f} mm rod, prints flat; cap holes {2 * (ROD_R + ROD_CLEAR):.1f} mm")
    lo = np.minimum(ribbon.bounds[0], rod_in_place.bounds[0])
    hi = np.maximum(ribbon.bounds[1], rod_in_place.bounds[1])
    print(f"assembled: {(hi - lo).round(1)} mm")
    if cuts:
        (net, _), (run, run_mm, _) = rib.integrity()
        print(f"  integrity: thinnest cross-section {100 * net:.0f} % of the band, longest "
              f"cut {100 * run:.0f} % ({run_mm:.1f} mm); limits {100 * pr.MIN_NET:.0f} / "
              f"{100 * pr.MAX_CUT:.0f} %")
    render([ribbon, rod_in_place], "negative_space_desk.png")
    print("wrote negative_space_desk.png")


if __name__ == "__main__":
    main()
