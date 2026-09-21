"""Anterior silhouette of a standing adult, with measured per-segment skin temperature
and measured per-segment heat-transfer coefficients attached to each region.

Every temperature and coefficient in SEGMENTS is measured data, not a modelling choice.
Sources are named per field in SOURCES below; nothing here is invented to make the
simulation behave. Where a segment the solver needs was not measured directly (the
chest and back are absent from the 15-site supine study), the value is taken from the
second mapping study and flagged with `derived=True`.

The rasteriser turns the silhouette into three arrays the LBM solver consumes:

    solid      (ny, nx) bool   -- True inside the body
    seg_id     (ny, nx) int    -- index into SEGMENTS, -1 in the fluid
    theta_wall (ny, nx) float  -- normalised wall temperature at solid nodes

Normalisation is theta = (T_skin - T_ambient) / DELTA_T_REF, where DELTA_T_REF is the
hottest segment's excess over ambient, so theta runs 1.0 at the anterior neck down to
about 0.53 at the sole of the foot.
"""

import numpy as np
from matplotlib.path import Path

# ---------------------------------------------------------------------------
# Measurement conditions
# ---------------------------------------------------------------------------

T_AMBIENT = 24.4        # degC -- ambient of the 30-subject SST study
BODY_HEIGHT = 1.70      # m    -- standing height the silhouette is scaled to

SOURCES = {
    "T_skin": (
        "Regional skin surface temperature, 30 subjects, 24.4 degC ambient, 46% RH, "
        "15 exposed sites, supine (PMC7992731). Ten statistically distinct thermal "
        "zones spanning 29.8-35.0 degC."
    ),
    "T_skin_derived": (
        "Chest/back/pelvis are not among the 15 exposed sites of PMC7992731. Values "
        "follow the Chudecka thermal-mapping result that the trunk runs warmest "
        "overall, placed between the measured cheek (33.2) and upper arm (32.6)."
    ),
    "h_r,h_c": (
        "16-segment articulated thermal manikin in a boundary-layer wind tunnel, "
        "standing posture, still air (v < 0.1 m/s). Radiative and natural-convective "
        "coefficients isolated per body part."
    ),
    "transition": (
        "Natural convection over the standing body stays laminar to about 0.8 m of "
        "body height for 306 K skin against 298 K air, then transitions to turbulent."
    ),
}

# ---------------------------------------------------------------------------
# Segment geometry + measured properties
# ---------------------------------------------------------------------------
# Geometry is in metres, origin at the floor on the body's midline, +y up.
# 'capsule' = (x0, y0, x1, y1, radius); 'disc' = (xc, yc, radius); 'poly' = vertices.
# 'mirror' duplicates the shape at -x (limbs).
#
# T_skin  degC, measured skin surface temperature
# h_r     W/m^2K, measured radiative coefficient
# h_c     W/m^2K, measured natural-convective coefficient, standing, still air

SEGMENTS = [
    dict(name="head", zone="face 34.7-35.0", kind="disc",
         geom=(0.0, 1.600, 0.095), mirror=False,
         T_skin=34.8, h_r=4.1, h_c=3.6),

    dict(name="neck", zone="anterior neck, warmest site", kind="poly",
         geom=[(-0.055, 1.455), (0.055, 1.455), (0.055, 1.545), (-0.055, 1.545)],
         mirror=False,
         T_skin=35.0, h_r=4.1, h_c=3.6),

    dict(name="chest", zone="trunk, warmest overall (mapping study)", kind="poly",
         geom=[(-0.225, 1.430), (0.225, 1.430), (0.160, 1.100), (-0.160, 1.100)],
         mirror=False, derived=True,
         T_skin=33.6, h_r=4.5, h_c=3.0),

    dict(name="abdomen", zone="lower trunk", kind="poly",
         geom=[(-0.160, 1.100), (0.160, 1.100), (0.145, 1.010), (-0.145, 1.010)],
         mirror=False, derived=True,
         T_skin=32.8, h_r=4.4, h_c=2.9),

    dict(name="pelvis", zone="pelvic region", kind="poly",
         geom=[(-0.145, 1.010), (0.145, 1.010), (0.170, 0.930), (-0.170, 0.930)],
         mirror=False, derived=True,
         T_skin=32.2, h_r=4.2, h_c=3.4),

    dict(name="upper_arm", zone="upper arm 32.6", kind="capsule",
         geom=(0.200, 1.395, 0.235, 1.105, 0.048), mirror=True,
         T_skin=32.6, h_r=5.2, h_c=2.9),

    dict(name="forearm", zone="forearm 32.6", kind="capsule",
         geom=(0.235, 1.105, 0.250, 0.860, 0.040), mirror=True,
         T_skin=32.6, h_r=4.9, h_c=3.7),

    dict(name="hand", zone="palm 33.8, second-warmest zone", kind="capsule",
         geom=(0.250, 0.860, 0.252, 0.740, 0.042), mirror=True,
         T_skin=33.8, h_r=4.1, h_c=4.1),

    dict(name="thigh", zone="thigh 31.9, coolest of the proximal sites", kind="capsule",
         geom=(0.085, 0.930, 0.095, 0.505, 0.075), mirror=True,
         T_skin=31.9, h_r=4.3, h_c=4.1),

    dict(name="lower_leg", zone="shin 32.6-32.8, warmer than the thigh above it",
         kind="capsule",
         geom=(0.095, 0.505, 0.085, 0.105, 0.052), mirror=True,
         T_skin=32.7, h_r=5.3, h_c=4.1),

    dict(name="foot", zone="sole 29.8-30.1, coolest zone", kind="poly",
         geom=[(0.035, 0.000), (0.165, 0.000), (0.165, 0.065), (0.035, 0.085)],
         mirror=True,
         T_skin=30.0, h_r=3.9, h_c=5.1),
]

DELTA_T_REF = max(s["T_skin"] for s in SEGMENTS) - T_AMBIENT   # 10.6 K, the neck


# ---------------------------------------------------------------------------
# Air properties and the dimensionless groups they imply
# ---------------------------------------------------------------------------

def air_properties(t_film_c=None):
    """Air properties at the film temperature between skin and ambient."""
    if t_film_c is None:
        t_mean_skin = float(np.mean([s["T_skin"] for s in SEGMENTS]))
        t_film_c = 0.5 * (t_mean_skin + T_AMBIENT)
    t_k = t_film_c + 273.15
    return dict(
        T_film_C=t_film_c,
        T_film_K=t_k,
        beta=1.0 / t_k,             # 1/K, ideal gas
        nu=1.57e-5,                 # m^2/s at ~300 K
        alpha=2.21e-5,              # m^2/s
        Pr=0.71,
        k=0.0263,                   # W/mK
        g=9.81,
    )


def rayleigh(length, delta_t, props=None):
    """Rayleigh number over `length` for a temperature excess `delta_t`."""
    p = props or air_properties()
    return p["g"] * p["beta"] * delta_t * length ** 3 / (p["nu"] * p["alpha"])


def free_fall_velocity(length=BODY_HEIGHT, delta_t=DELTA_T_REF, props=None):
    """Buoyant velocity scale sqrt(g beta dT L) -- the sim's velocity unit."""
    p = props or air_properties()
    return float(np.sqrt(p["g"] * p["beta"] * delta_t * length))


def transition_height(delta_t=8.0, ra_crit=1.0e9, props=None):
    """Height at which the local Rayleigh number reaches `ra_crit`."""
    p = props or air_properties()
    return float((ra_crit * p["nu"] * p["alpha"] / (p["g"] * p["beta"] * delta_t)) ** (1 / 3))


# ---------------------------------------------------------------------------
# Rasteriser
# ---------------------------------------------------------------------------

def _capsule_distance(px, py, x0, y0, x1, y1):
    """Distance from each (px, py) to the segment (x0,y0)-(x1,y1)."""
    dx, dy = x1 - x0, y1 - y0
    denom = dx * dx + dy * dy
    if denom == 0.0:
        return np.hypot(px - x0, py - y0)
    t = ((px - x0) * dx + (py - y0) * dy) / denom
    t = np.clip(t, 0.0, 1.0)
    return np.hypot(px - (x0 + t * dx), py - (y0 + t * dy))


def _shape_mask(seg, X, Y, flip):
    sx = -1.0 if flip else 1.0
    kind, geom = seg["kind"], seg["geom"]
    if kind == "disc":
        xc, yc, r = geom
        return np.hypot(X - sx * xc, Y - yc) <= r
    if kind == "capsule":
        x0, y0, x1, y1, r = geom
        return _capsule_distance(X, Y, sx * x0, y0, sx * x1, y1) <= r
    if kind == "poly":
        verts = [(sx * vx, vy) for vx, vy in geom]
        pts = np.column_stack([X.ravel(), Y.ravel()])
        return Path(verts).contains_points(pts).reshape(X.shape)
    raise ValueError(f"unknown segment kind {kind!r}")


def rasterise(nx, ny, dx, body_x_center, body_y_floor):
    """Rasterise the silhouette onto an (ny, nx) lattice.

    `dx` is metres per lattice cell; the body's midline sits at lattice column
    `body_x_center` and its feet on lattice row `body_y_floor`.

    Returns solid, seg_id, theta_wall.
    """
    jj, ii = np.meshgrid(np.arange(ny), np.arange(nx), indexing="ij")
    X = (ii - body_x_center) * dx
    Y = (jj - body_y_floor) * dx

    solid = np.zeros((ny, nx), dtype=bool)
    seg_id = np.full((ny, nx), -1, dtype=np.int16)

    for idx, seg in enumerate(SEGMENTS):
        m = _shape_mask(seg, X, Y, flip=False)
        if seg.get("mirror"):
            m |= _shape_mask(seg, X, Y, flip=True)
        # First segment to claim a cell keeps it, so limbs do not overwrite the torso.
        fresh = m & ~solid
        seg_id[fresh] = idx
        solid |= m

    theta_wall = np.zeros((ny, nx), dtype=np.float64)
    for idx, seg in enumerate(SEGMENTS):
        theta_wall[seg_id == idx] = (seg["T_skin"] - T_AMBIENT) / DELTA_T_REF

    return solid, seg_id, theta_wall


def segment_table():
    """The measured table, with the dimensionless wall temperature each row implies."""
    rows = []
    for idx, s in enumerate(SEGMENTS):
        rows.append(dict(
            index=idx,
            name=s["name"],
            T_skin=s["T_skin"],
            dT=s["T_skin"] - T_AMBIENT,
            theta=(s["T_skin"] - T_AMBIENT) / DELTA_T_REF,
            h_r=s["h_r"],
            h_c=s["h_c"],
            radiative_share=s["h_r"] / (s["h_r"] + s["h_c"]),
            derived=bool(s.get("derived", False)),
            zone=s["zone"],
        ))
    return rows


if __name__ == "__main__":
    p = air_properties()
    print(f"film temperature      {p['T_film_C']:.1f} degC")
    print(f"ambient               {T_AMBIENT} degC")
    print(f"reference excess dT   {DELTA_T_REF:.1f} K  (anterior neck)")
    print(f"Ra over body height   {rayleigh(BODY_HEIGHT, DELTA_T_REF, p):.3e}")
    print(f"free-fall velocity    {free_fall_velocity():.3f} m/s")
    print(f"buoyant time scale    {BODY_HEIGHT / free_fall_velocity():.2f} s")
    print()
    print("local Rayleigh number along the body, dT = 8 K (manikin conditions):")
    for h in (0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.7):
        print(f"   x = {h:.1f} m   Ra_x = {rayleigh(h, 8.0, p):.3e}")
    print(f"Ra_x = 1e9 reached at x = {transition_height(8.0):.2f} m "
          f"(reported transition: 0.8 m)")
    print(f"Ra_x = 1e8 reached at x = {transition_height(8.0, 1e8):.2f} m")
    print()
    print(f"{'segment':<11}{'T':>6}{'dT':>6}{'theta':>7}{'h_r':>6}{'h_c':>6}{'rad%':>7}")
    for r in segment_table():
        flag = " *" if r["derived"] else ""
        print(f"{r['name']:<11}{r['T_skin']:>6.1f}{r['dT']:>6.1f}{r['theta']:>7.3f}"
              f"{r['h_r']:>6.1f}{r['h_c']:>6.1f}{100 * r['radiative_share']:>6.0f}%{flag}")
    print("  * chest/abdomen/pelvis not among the 15 measured sites -- see SOURCES")
