#!/usr/bin/env python3
"""Printable involute gears by the generating method (shapely).

External gears are cut by rolling a basic rack cutter (20 deg, addendum 1.25 m
on the cutter so it also cuts the root), which gives undercut and profile shift
correctly without closed-form involute bookkeeping.  The internal ring is cut by
rolling the mating planet, enlarged by the backlash and tip clearance, round
inside it as a shaper cutter, which makes the ring conjugate to the planet as
printed.

Run directly for a mesh self-check:  python gears.py
"""
import math
import numpy as np
from shapely.geometry import Polygon, Point
from shapely.ops import unary_union
from shapely import affinity

PA = math.radians(20)

def _rack(m, n_teeth, thin):
    """Rack cutter in its own frame: pitch line y = 0, cutter teeth pointing to -y."""
    p = math.pi * m
    w0 = p / 4 + thin / 2                      # cutter tooth half-width at the pitch line
    teeth = []
    for k in range(-n_teeth, n_teeth + 1):
        u = k * p + p / 2                       # cutter teeth sit in the gear's spaces
        lo, hi = -1.25 * m, 2.0 * m
        hw = lambda y: w0 + y * math.tan(PA)    # cutter tooth narrows toward its tip (y < 0)
        teeth.append(Polygon([(u - hw(lo), lo), (u + hw(lo), lo), (u + hw(hi), hi), (u - hw(hi), hi)]))
    body = Polygon([(-(n_teeth + 1) * p, 1.5 * m), ((n_teeth + 1) * p, 1.5 * m),
                    ((n_teeth + 1) * p, 4 * m), (-(n_teeth + 1) * p, 4 * m)])
    return unary_union(teeth + [body])

def external(Z, m, x=0.0, backlash=0.15, steps=320):
    """External spur gear outline, one tooth centred on +y."""
    rp = m * Z / 2
    ra = rp + m * (1 + x)
    rack = _rack(m, 4, backlash)
    span = 4.0 * math.pi / Z                                        # the neighbouring spaces get fully cut
    cuts = []
    for phi in np.linspace(-span, span, steps):
        r = affinity.translate(rack, -rp * phi, rp + x * m)       # rack rolls along the pitch circle
        cuts.append(affinity.rotate(r, -phi, origin=(0, 0), use_radians=True))
    blank = Point(0, 0).buffer(ra, 256)
    cut = blank.difference(unary_union(cuts))
    half = math.pi / Z
    wedge = Polygon([(0, 0)] + [(2 * ra * math.cos(math.pi / 2 + a), 2 * ra * math.sin(math.pi / 2 + a))
                                for a in np.linspace(-half, half, 16)])
    tooth = cut.intersection(wedge)
    rf = rp - m * (1.25 - x)
    gear = unary_union([affinity.rotate(tooth, 2 * math.pi * k / Z, origin=(0, 0), use_radians=True)
                        for k in range(Z)] + [Point(0, 0).buffer(rf + 0.05, 128)])
    return _main_body(gear)

def _main_body(g):
    """Largest piece of a generated outline, after closing hairline seams; any
    leftover fragment above 0.05 mm^2 is an error, not something to drop silently."""
    g = g.buffer(0.02).buffer(-0.02)
    if g.geom_type == "Polygon":
        return g.simplify(0.002)
    parts = sorted(g.geoms, key=lambda p: -p.area)
    assert all(p.area < 0.05 for p in parts[1:]), f"gear fragments: {[round(p.area, 3) for p in parts[1:6]]}"
    return parts[0].simplify(0.002)

def internal(Zr, Zp, m, planet, r_outer, x_p=0.0, clearance=0.2, steps=520):
    """Internal ring cut by the planet outline rolling inside it.  `planet` is the
    planet's shapely outline with a tooth centred on +y; x_p its profile shift."""
    a = m * (Zr - Zp) / 2
    rp_p = m * Zp / 2
    cutter = planet.buffer(clearance)                               # backlash + tip clearance
    span = 0.85                                                      # every planet tooth passes the cut sector
    cuts = []
    for psi in np.linspace(-span, span, steps):
        # absolute planet spin is -psi (Zr - Zp) / Zp; the orbit rotation below adds psi
        phi = -psi * Zr / Zp
        c = affinity.rotate(cutter, phi, origin=(0, 0), use_radians=True)
        c = affinity.translate(c, 0, a)
        cuts.append(affinity.rotate(c, psi, origin=(0, 0), use_radians=True))
    r_tip = a + (rp_p - m * (1.25 + x_p)) + 0.4                      # ring tips clear the planet root
    blank = Point(0, 0).buffer(r_outer, 256).difference(Point(0, 0).buffer(r_tip, 256))
    cut = blank.difference(unary_union(cuts))
    half = math.pi / Zr
    wedge = Polygon([(0, 0)] + [(2 * r_outer * math.cos(math.pi / 2 + q), 2 * r_outer * math.sin(math.pi / 2 + q))
                                for q in np.linspace(-half, half, 16)])
    seg = cut.intersection(wedge)
    ring = unary_union([affinity.rotate(seg, 2 * math.pi * k / Zr, origin=(0, 0), use_radians=True)
                        for k in range(Zr)])
    return _main_body(ring)

def mesh_check(sun, planet, ring, Zs, Zp, Zr, m, n_planets=3, samples=24):
    """Worst overlap area (mm^2) between the sun, planets and ring over one tooth
    pitch of carrier motion, with the ring fixed.  Zero means no interference."""
    a = m * (Zs + Zp) / 2
    worst = 0.0
    for k in range(samples):
        th_c = 2 * math.pi / (Zs + Zr) * k / samples                       # one mesh pitch of carrier
        th_s = th_c * (1 + Zr / Zs)                                        # sun angle (ring fixed)
        th_p = -th_c * Zr / Zp                                            # planet spin (absolute)
        s = affinity.rotate(sun, th_s, origin=(0, 0), use_radians=True)
        for j in range(n_planets):
            q = th_c + 2 * math.pi * j / n_planets
            # planet frame: tooth on +y at zero spin; meshing phase with the sun at q = 0
            p = affinity.rotate(planet, th_p + math.pi / Zp * (1 if Zp % 2 == 0 else 0),
                                origin=(0, 0), use_radians=True)
            p = affinity.translate(p, 0, a)
            p = affinity.rotate(p, q, origin=(0, 0), use_radians=True)
            worst = max(worst, s.intersection(p).area, ring.intersection(p).area)
    return worst

if __name__ == "__main__":
    Zs, Zp, Zr, m = 15, 21, 57, 1.5
    sun = external(Zs, m, x=+0.15)
    planet = external(Zp, m, x=-0.15)
    ring = internal(Zr, Zp, m, planet, r_outer=Zr * m / 2 + 7, x_p=-0.15)
    print("areas", round(sun.area), round(planet.area), round(ring.area))
    print("worst overlap mm^2:", round(mesh_check(sun, planet, ring, Zs, Zp, Zr, m), 3))
