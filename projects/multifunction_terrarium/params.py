"""Every dimension in this object is traceable to a physical constant.

That is the form language (IDEA.md sec.6): nothing here is styled. If you change a
number, change the physics comment with it or the object stops meaning anything.
"""

import math

# ---------------------------------------------------------------- fluid constants
SIGMA   = 0.072          # N/m    water surface tension, 20 C
RHO     = 1000.0         # kg/m3
GRAV    = 9.81           # m/s2
MU      = 1.0e-3         # Pa.s
N_PETG  = 1.57           # refractive index

#: the water's own wavelength -- sets BOTH the floor on lattice cell size and the
#: crossover between "texture as effective medium" and "texture as discrete features"
CAPILLARY_LEN = math.sqrt(SIGMA / (RHO * GRAV)) * 1e3        # 2.709 mm

#: max height a sessile puddle can stand before gravity flattens it. Cell walls
#: shorter than this get bulged over and the maze stops existing.
PUDDLE_MAX_H  = 2 * CAPILLARY_LEN * math.sin(math.radians(78) / 2)   # 3.42 mm

#: Rayleigh-Taylor spacing of pendant drops on a wetted ceiling. Drip sites placed
#: closer than lam_c are ours; wider apart and the physics chooses them instead.
RT_LAM_C   = 2 * math.pi * CAPILLARY_LEN                      # 17.02 mm

# ---------------------------------------------------------------- printer / material
NOZZLE      = 0.4
LAYER       = 0.2
BED         = 320.0      # 320 x 320 x 320, PETG
PANE_T      = 0.6        # flat-printed port panes: 3 layers, optically flat both faces
WALL_THIN   = 0.9        # 2 perimeters -- lattice ribs, screen wall
WALL_WET    = 1.6        # anything below a waterline
OVERHANG_OK = 45.0       # deg from vertical, support-free

# ---------------------------------------------------------------- drop production
LIP_W       = NOZZLE     # drip lip / notch root width. Tate: V is linear in this.
def drop_volume(d_mm):
    """Tate's law with the Harkins-Brown correction. mm -> uL."""
    return 0.6 * math.pi * (d_mm * 1e-3) * SIGMA / (RHO * GRAV) * 1e9

DROP_V   = drop_volume(LIP_W)                                 # 5.54 uL
DROP_D   = (6 * DROP_V * 1e-9 / math.pi) ** (1 / 3) * 1e3     # 2.19 mm

#: above this a site stops dripping and becomes a jet -- total loss of drop count.
#: We = 4 at the orifice.
_v_jet     = math.sqrt(4 * SIGMA / (RHO * LIP_W * 1e-3))
JET_LIMIT  = _v_jet * math.pi * (LIP_W * 1e-3 / 2) ** 2 / (DROP_V * 1e-9)   # 19.3 /s

#: a drop must clear the lip before touching down or it bridges and no detachment
#: event happens at all. ~4x drop diameter.
STEP_MIN   = 4 * DROP_D                                       # 8.8 mm
STEP       = 10.0        # terrace riser. drops-in-air ~ 1/sqrt(STEP), so: as short
                         # as STEP_MIN allows.
NOTCH_P    = 6.0         # notch pitch along the lip. > 2*DROP_D (no coalescence),
                         # and small enough to keep per-site rate under JET_LIMIT.

# ---------------------------------------------------------------- lattice
CELL_FREE  = 3.0         # >= CAPILLARY_LEN or a single drop spans the cell and the
                         # maze disappears. This is a physics floor, not a print one.
RIB_W      = WALL_THIN
CELL_PITCH = CELL_FREE + RIB_W                                # 3.9 mm
RIB_H      = 3.5         # >= PUDDLE_MAX_H, else cells bulge over and merge

VOGEL_P    = 0.50        # r = c*k^p. 0.5 = equal area (uniform flow balance);
                         # aperiodicity comes from the 5/6/7 cell defects and the
                         # 46% scatter in edge length, not from grading the areas.
GOLDEN_ANG = math.pi * (3 - math.sqrt(5))

# ---------------------------------------------------------------- optics / viewing
FRESNEL_MAX = 45.0       # deg incidence. <=45 the window reads clear (<12% loss);
                         # by 70 deg it is a mirror.
EYE_STAND   = 1550.0
EYE_SIT     = 1200.0
VIEW_DIST   = 700.0
TABLE_H     = 750.0

# ---------------------------------------------------------------- global form
HEX_R       = 155.0      # circumradius. across-corners 310 <= 320 bed and footprint.
HEX_AF      = HEX_R * math.sqrt(3)                            # 268.5 across flats
N_MODULE    = 3
LEVELS      = 14         # cascade cone in the mid volume
MOD_H       = LEVELS * STEP                                   # 140 mm
TOTAL_H     = N_MODULE * MOD_H                                # 630 mm

SCREEN_HEX   = True      # hexagonal, phase-locked to the shell: a uniform 30 mm
                         # theatre gap all the way round, and 15% more lip than a
                         # circle that has to clear the same across-flats
SHELL_IN     = HEX_R - 3.0                                    # 152 circumradius
THEATRE      = 30.0      # the gap the drops fall through, seen through the ports
SCREEN_R_BOT = 100.0     # cascade cone: narrow at top, widening as it descends so
SCREEN_R_TOP =  60.0     # every lip sits over the next tread and the drop re-forms
TREAD        = (SCREEN_R_BOT - SCREEN_R_TOP) / LEVELS         # 2.86 mm
SCREEN_W     = WALL_THIN

#: travertine lobes. Scalloping the lip buys arc length -- and therefore drip sites --
#: without moving the cone out toward the shell.
LOBES     = 36
LOBE_A    = 8.0


# ---------------------------------------------------------------- flow
Q_TRICKLE = 6.0          # L/h, solar pump through a restrictor
Q_DUMP    = 48.0         # L/h equivalent during a bell-siphon discharge
SIPHON_RATIO = (1, 2, 3) # flood volumes; periods beat on a common 6T

def per_site_rate(q_Lh, n_sites):
    """drops/s at one notch. Must stay under JET_LIMIT at Q_DUMP."""
    return (q_Lh / 3600 * 1e-3) / (DROP_V * 1e-9) / n_sites

def band_tilt(z):
    """Facet tilt that splits the difference between a standing and a seated eye.
    Worst-case incidence over the whole object comes out at 14 deg."""
    ds = math.degrees(math.atan2(EYE_STAND - z, VIEW_DIST))
    di = math.degrees(math.atan2(EYE_SIT   - z, VIEW_DIST))
    return (ds + di) / 2, max(abs(ds - (ds + di) / 2), abs(di - (ds + di) / 2))

if __name__ == "__main__":
    print(f"capillary length   {CAPILLARY_LEN:6.3f} mm   <- cell-size floor")
    print(f"puddle max height  {PUDDLE_MAX_H:6.2f} mm   <- rib height floor")
    print(f"R-T drip spacing   {RT_LAM_C:6.2f} mm")
    print(f"drop               {DROP_V:6.2f} uL, {DROP_D:.2f} mm dia from a {LIP_W} mm lip")
    print(f"jet limit          {JET_LIMIT:6.1f} drops/s per site")
    print(f"terrace step       {STEP:6.1f} mm  (floor {STEP_MIN:.1f})")
    print(f"cell pitch         {CELL_PITCH:6.2f} mm   rib {RIB_W} x {RIB_H} mm")
    print(f"object             {HEX_R*2:.0f} across corners x {TOTAL_H:.0f} tall, "
          f"{N_MODULE} modules of {MOD_H:.0f}")
    print(f"tread              {TREAD:6.2f} mm  (>= capillary {CAPILLARY_LEN:.2f}? "
          f"{'yes' if TREAD >= CAPILLARY_LEN else 'NO'})")


# ---------------------------------------------------------------- vessel form
# Ovoid crown over a slightly tapered 12-sided body. Shared by every module so
# the shell, the cascade and the drip spirals all agree on one surface.
import numpy as _np
N_FACE   = 12
BODY_H   = 360.0
CROWN_H  = 180.0
R_BODY   = 152.0
VESSEL_H = BODY_H + CROWN_H          # 540
SHELL_W  = 2.4                       # structural frame
SPLIT_Z  = (0.0, 180.0, 360.0, VESSEL_H)   # printable modules, each <= 320 tall

Z_RES    = 92.0                      # waterline / top of the reservoir
Z_BED    = 205.0                     # top of the living bed tray

def vessel_r(z):
    z = _np.asarray(z, float)
    body = R_BODY - 5.0 * _np.clip(z, 0, BODY_H) / BODY_H
    t = _np.clip((z - BODY_H) / CROWN_H, 0, 1)
    crown = (R_BODY - 5.0) * _np.sqrt(_np.clip(1 - t * t, 0, 1))
    return _np.where(z <= BODY_H, body, crown)

def hexf(theta, n=N_FACE):
    """Radius multiplier turning a circle into a regular n-gon."""
    a = _np.pi / n
    return _np.cos(a) / _np.cos((_np.asarray(theta) % (2 * a)) - a)


def crown_inner_z(r, wall=SHELL_W, n_face=N_FACE):
    """Height of the shell's INNER face at radius r, on the crown.

    Solved off the same buffered profile vessel._sweep builds, not an analytic
    normal offset. Near the apex the profile turns sharply and the mitred buffer
    departs from a smooth offset by up to 0.76 mm -- more than enough to bury a
    0.4 mm seat clearance and leave the condenser tap embedded in the wall.

    The 12-gon's inner FACES sit at cos(pi/n) of the section radius and a face is
    the closest the wall ever comes to the axis, so that is what a seated part has
    to clear.
    """
    from shapely.geometry import LineString
    zz = _np.linspace(BODY_H, VESSEL_H, 60)
    pts = list(zip(vessel_r(zz), zz))
    dense = []
    for a, b in zip(pts[:-1], pts[1:]):
        a, b = _np.asarray(a, float), _np.asarray(b, float)
        n = max(1, int(_np.ceil(_np.linalg.norm(b - a) / 2.5)))
        dense += [tuple(a + (b - a) * i / n) for i in range(n)]
    dense.append(tuple(pts[-1]))
    poly = LineString(dense).buffer(wall / 2, cap_style=2, join_style=2, mitre_limit=6)
    if poly.geom_type == "MultiPolygon":
        poly = max(poly.geoms, key=lambda g: g.area)
    r_sec = r / math.cos(math.pi / n_face)
    hit = poly.intersection(LineString([(r_sec, BODY_H - 50), (r_sec, VESSEL_H + 50)]))
    if hit.is_empty:
        raise ValueError(f"radius {r} is off the crown")
    return min(c[1] for g in getattr(hit, "geoms", [hit]) for c in g.coords)
