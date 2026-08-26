"""Ovoid-crowned vessel with a descending drip spiral.

Fed two ways, as asked: condensate running down the crown falls into the open
gutters, and a 10 mm intake tube feeds a four-level splitter that divides into
sixteen 2.5 mm bores.

Why the spiral is an OPEN C-section and not a sealed pipe: at 2.5 mm the bore sits
below the capillary length (Bond 0.85), so a hole in the underside of a closed pipe
needs ~576 Pa to break its meniscus while the bore only supplies ~25 Pa of head. It
would hold every drop and discharge only at the far end. An open gutter has a free
surface, no Laplace barrier, and sheds from a notched lip by Tate's law. It also
catches condensation, which a sealed pipe could not.
"""

import math
import numpy as np
import trimesh
import params as P

BORE      = 2.5          # ID, as specified
GWALL     = 0.9
R_IN, R_OUT = BORE / 2, BORE / 2 + GWALL
LIP_A     = 40.0         # deg above horizontal: the outer spill lip
INNER_A   = 140.0        # inner wall runs higher so water spills outward only
NOTCH_A   = 4.0          # lip drops to this at a notch -> the drip site
NOTCH_P   = 12.0         # mm along the gutter
N_SPIRAL  = 16           # 4 splitter levels; >=9 needed to stay under the jet limit
TURNS     = 1.35
LEVEL_FRAC= 0.72         # of each step that is level gutter; the rest is the drop
Z_TOP, Z_BOT = 500.0, 215.0

VESSEL_H = P.VESSEL_H
vessel_r = P.vessel_r

# ---------------------------------------------------------------- crown collector
# The cone is a CONDENSER TAP, not a finial. Its flat seat touches the dome's inner
# face at r=SEAT_R: that wetted contact line is the collection mechanism. Film
# creeping down the ceiling is drawn into the contact instead of nucleating pendant
# drops at RT_LAM_C (17.0 mm) spacing and raining wherever it likes. Inside SEAT_R
# the seat hangs ~0.9 mm under the apex -- far below CAPILLARY_LEN, so that gap
# fills and feeds the same contact line rather than dripping.
SEAT_R    = 12.0         # contact circle. Larger than the intake OD so the seat is
                         # an annulus around the stem, not a point on the apex.
CONE_R    = 24.0         # base rim: condensate runs down the flank and sheds here,
                         # onto the branches below, which carry it on as surface flow
CONE_H    = 6.0          # short on purpose -- every mm here is a mm the fan-out
                         # loses, and the fan-out is angle-starved. Apex-up, so it
                         # is self-supporting at any half-angle.
STEM_UP   = 26.0         # intake above the seat, through the apex bore
LV_FRAC   = (0.0, 0.10, 0.25, 0.55, 1.0)
                         #: radial schedule. Spread LATE: the last level has 16 short
                         #: branches, the first has 2 long ones, so pushing the travel
                         #: outward buys 3 deg of overhang over an even schedule.
#: first bifurcation. Sits as high as its own node sphere allows under the seat --
#: any higher and the sphere breaks through the cone's seat and out of the dome.
#: Every mm it gains is a mm of fan-out, and the fan-out is what is angle-starved.
def _z_root():
    return crown_in_z(SEAT_R) - (BORE * 4 / 2 + GWALL)


def crown_in_z(r):
    """Height of the crown's inner face at radius r -- see params.crown_inner_z."""
    return P.crown_inner_z(r)


def _teardrop(rb, n=26):
    """Bore section: circle with a 45 deg roof peak.

    The branches run 69 deg off vertical -- 21 deg off HORIZONTAL -- so a round bore
    has a nearly flat roof with nothing under it and sags shut on the first layer
    that has to bridge it. The peak lets every layer self-support."""
    a = np.linspace(math.radians(135), math.radians(405), n)
    return np.vstack([np.c_[rb * np.cos(a), rb * np.sin(a)], [0.0, rb * math.sqrt(2)]])


def r_outer(rb):
    """Outer radius that leaves GWALL over the teardrop PEAK, not over the circle.

    The peak stands at rb*sqrt(2). Sizing the tube to rb+GWALL -- the obvious thing,
    and what the area-conserving formula gives -- puts the peak outside the wall on
    the two trunk levels and slots the tube open along its whole length."""
    return rb * math.sqrt(2) + GWALL


def _bore(p, q, rb, over=(0.6, 0.6)):
    """Bore cutter along p->q. Teardrop when the run is shallow enough to need a
    self-supporting roof; plain round when it is near-vertical and has none."""
    from shapely.geometry import Polygon
    p, q = np.asarray(p, float), np.asarray(q, float)
    d = q - p; L = float(np.linalg.norm(d)); w = d / L
    if abs(w[2]) > math.cos(math.radians(P.OVERHANG_OK)):
        return trimesh.creation.cylinder(radius=rb, sections=28,
                                         segment=[p - w * over[0], q + w * over[1]])
    up = np.array([0.0, 0.0, 1.0]) - w * w[2]
    up /= np.linalg.norm(up)
    m = trimesh.creation.extrude_polygon(Polygon(_teardrop(rb)), L + over[0] + over[1])
    T = np.eye(4); T[:3, 0] = np.cross(up, w); T[:3, 1] = up; T[:3, 2] = w
    T[:3, 3] = p - w * over[0]
    return m.apply_transform(T)


def helix(k, n=1400):
    """One spiral, riding just proud of the vessel, descending as a STAIRCASE.

    A smooth 20-degree helix would race the water to its first notch and the whole
    drip zone would collapse to the top of each spiral. Stepping it into level runs
    with STEP drops between makes each run fill and spill along its whole notched
    lip -- the same accumulate/threshold law as the internal terraces, wrapped
    round the outside."""
    t = np.linspace(0, 1, n)
    nstep = max(1, int(round((Z_TOP - Z_BOT) / P.STEP)))
    f = t * nstep
    k0 = np.floor(f); u = np.clip((f - k0 - LEVEL_FRAC) / (1 - LEVEL_FRAC), 0, 1)
    z = Z_TOP - (Z_TOP - Z_BOT) * (k0 + u * u * (3 - 2 * u)) / nstep
    th = k * 2 * math.pi / N_SPIRAL + 2 * math.pi * TURNS * t
    # INSIDE the shell: the vessel stays sealed, and the gutters catch the
    # condensate running down the crown as well as the trickled-in feed
    # inscribed radius, not circumradius: a 12-sided shell's inner FACE sits at
    # circumradius * cos(pi/12) = 0.966x. Using the circumradius put the gutters
    # 3.4 mm through the body wall and 15.7 mm through the crown.
    r = vessel_r(z) * math.cos(math.pi / P.N_FACE) - P.SHELL_W - 2 * R_OUT - 1.5
    return np.c_[r * np.cos(th), r * np.sin(th), z], nstep

def _frames(pts):
    T = np.gradient(pts, axis=0)
    T /= np.linalg.norm(T, axis=1, keepdims=True)
    rad = pts.copy(); rad[:, 2] = 0
    rad /= np.linalg.norm(rad, axis=1, keepdims=True)
    out = rad - (np.sum(rad * T, axis=1, keepdims=True)) * T
    out /= np.linalg.norm(out, axis=1, keepdims=True)
    up = np.cross(T, out)
    up[up[:, 2] < 0] *= -1
    return out, up

def _section(lip_deg, n_arc=26):
    """Closed C: outer arc from the lip the long way round, inner arc back.

    Forced counter-clockwise. The sweep's side quads and its end-cap fans both
    assume CCW; the raw construction comes out CW, which left the cap/side
    junction inconsistently wound (104 bad edge pairs = 2 x 52 section points)."""
    a0, a1 = math.radians(lip_deg), math.radians(INNER_A) - 2 * math.pi
    A = np.linspace(a0, a1, n_arc)
    sec = np.r_[np.c_[R_OUT * np.cos(A), R_OUT * np.sin(A)],
                np.c_[R_IN * np.cos(A[::-1]), R_IN * np.sin(A[::-1])]]
    area2 = np.sum(sec[:, 0] * np.roll(sec[:, 1], -1)
                   - np.roll(sec[:, 0], -1) * sec[:, 1])
    return sec if area2 > 0 else sec[::-1].copy()

def gutter(pts, notch_p=NOTCH_P):
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))]
    out, up = _frames(pts)
    u = (s % notch_p) / notch_p
    lip = np.where(np.abs(u - 0.5) < 0.06, NOTCH_A, LIP_A)      # sharp dip = drip site
    secs = [_section(l) for l in lip]
    ns = len(secs[0])
    V = np.empty((len(pts), ns, 3))
    for i, (c, o, w, sec) in enumerate(zip(pts, out, up, secs)):
        V[i] = c + sec[:, 0, None] * o + sec[:, 1, None] * w
    nt = len(pts)
    i = np.arange(nt - 1)[:, None]; j = np.arange(ns)[None, :]
    a = i * ns + j; b = (i + 1) * ns + j
    c_ = (i + 1) * ns + (j + 1) % ns; d = i * ns + (j + 1) % ns
    F = [np.stack([a, b, c_], -1).reshape(-1, 3), np.stack([a, c_, d], -1).reshape(-1, 3)]
    for base, flip in ((0, True), ((nt - 1) * ns, False)):        # end caps
        fan = [[base, base + m + 1, base + m] if flip else [base, base + m, base + m + 1]
               for m in range(1, ns - 1)]
        F.append(np.array(fan))
    m = trimesh.Trimesh(vertices=V.reshape(-1, 3), faces=np.concatenate(F), process=False)
    m.merge_vertices()
    trimesh.repair.fix_winding(m)      # safety net: propagate consistency
    trimesh.repair.fix_normals(m)      # and point them outward
    if m.volume < 0: m.invert()
    return m, s[-1], int(s[-1] / notch_p)

def splitter():
    """Condenser tap over an area-conserving tube tree: one 10 mm intake to sixteen
    2.5 mm outlets, every level 78.5 mm2.

    Two supplies, one tree. The PUMP feed arrives inside the bores, entering through
    the stem that crosses the apex. The CONDENSATE arrives outside: the cone's seat
    touches the dome, film wicks onto it, runs down the flank, sheds at the base rim
    and continues as surface flow on the branch exteriors to the same sixteen mouths.

    Angles: the crown gives 34 mm of drop to fan 0 -> 81 mm of radius, so a straight
    line from cone to outlet is already 67.2 deg from vertical. Nothing beats that in
    this envelope. LV_FRAC + centred sectors reach 69.1 against the 67.2 floor; it is
    well past OVERHANG_OK and the part needs support, which is why the bores are
    teardropped rather than round."""
    r_end = (float(vessel_r(Z_TOP)) * math.cos(math.pi / P.N_FACE)
             - P.SHELL_W - 2 * R_OUT - 1.5)      # where a spiral begins
    lv_r = [r_end * f for f in LV_FRAC]

    z_root = _z_root()

    def node(k, i):
        if k == 0:
            return np.array([0.0, 0.0, z_root])
        # centred sectors: children straddle the parent's angle instead of one child
        # inheriting it and the other flying a quarter-turn away. The -pi/N_SPIRAL
        # phase puts level 4 exactly on the spiral start angles.
        a = (i + 0.5) * 2 * math.pi / 2 ** k - math.pi / N_SPIRAL
        return np.array([lv_r[k] * math.cos(a), lv_r[k] * math.sin(a), lv_z[k]])

    # z per level proportional to that level's longest chord -> every branch at the
    # same angle, so no single level is the one that fails.
    lv_z = [z_root] * 5
    chord = []
    for k in range(4):
        c = 0.0
        for i in range(max(1, 2 ** k)):
            p = node(k, i)
            for ch in (2 * i, 2 * i + 1):
                c = max(c, float(np.hypot(*(p - node(k + 1, ch))[:2])))
        chord.append(c)
    drop = z_root - Z_TOP
    z = z_root
    for k, c in enumerate(chord):
        z -= drop * c / sum(chord); lv_z[k + 1] = z

    solid, void = [], []
    for k in range(4):
        r_wall = BORE * (math.sqrt(2) ** (3 - k)) / 2          # bore radius, level k
        r_node = BORE * (math.sqrt(2) ** (4 - k)) / 2          # parent bore radius
        for i in range(max(1, 2 ** k)):
            p = node(k, i)
            for ch in (2 * i, 2 * i + 1):
                q = node(k + 1, ch)
                solid.append(trimesh.creation.cylinder(
                    radius=r_outer(r_wall), segment=[p, q], sections=20))
                # last level overshoots its node so the sixteen outlets open
                solid.append(trimesh.creation.icosphere(subdivisions=2,
                             radius=r_outer(r_wall)).apply_translation(q))
                # the outlet cutter must exit CLEAR of the terminal sphere. Stop it
                # short and it slices a loose cap off each of the sixteen tips
                # (r_outer(1.25)=2.67 detaches at any overshoot >= 2.36).
                void.append(_bore(p, q, r_wall, over=(0.6, 4.0 if k == 3 else 0.6)))
            # r_outer allows for the PARENT teardrop's peak arriving at the node.
            # The root has no parent teardrop -- the stem bore above it is round,
            # near-vertical and needs no roof -- so it only needs a plain wall.
            solid.append(trimesh.creation.icosphere(
                subdivisions=2,
                radius=(r_node + GWALL) if k == 0 else r_outer(r_node)
            ).apply_translation(p))
            void.append(trimesh.creation.icosphere(subdivisions=3, radius=r_node)
                        .apply_translation(p))

    z_seat = crown_in_z(SEAT_R)
    r_in   = BORE * 4 / 2                                   # 10 mm intake bore
    # frustum: SEAT_R on top against the dome, CONE_R at the shedding rim
    cone = trimesh.creation.cylinder(radius=CONE_R, height=CONE_H, sections=48)
    v = cone.vertices.copy()
    top = v[:, 2] > 0
    v[top, :2] *= SEAT_R / CONE_R
    cone = trimesh.Trimesh(v, cone.faces, process=False)
    cone.apply_translation([0, 0, z_seat - CONE_H / 2])
    solid.append(cone)

    stem_top = z_seat + STEM_UP
    solid.append(trimesh.creation.cylinder(
        radius=r_in + GWALL, sections=28,
        segment=[[0, 0, z_root], [0, 0, stem_top]]))
    void.append(_bore([0, 0, z_root], [0, 0, stem_top], r_in, over=(0.0, 2.0)))

    m = trimesh.boolean.union(solid, engine="manifold")
    m = trimesh.boolean.difference([m] + void, engine="manifold")
    # manifold leaves ~190 zero-area faces. They cost nothing in memory and the mesh
    # reports watertight, but STL is float32: on round-trip the merge fails around
    # them and the part loads with 564 broken edges. Drop them before it ever
    # reaches a slicer.
    m.update_faces(m.nondegenerate_faces())
    m.merge_vertices()
    trimesh.repair.fix_normals(m)
    return m


def build(verbose=True, part="all"):
    """part: 'gutters' (285 tall), 'splitter' (90 tall), or 'all' for preview.
    Split because the two together are 421 mm and the bed is 320."""
    gutters, tot_len, tot_notch = [], 0.0, 0
    for k in range(N_SPIRAL):
        path, nstep = helix(k)
        g, L, n = gutter(path); gutters.append(g); tot_len += L; tot_notch += n
    pieces = {'gutters': gutters, 'splitter': [splitter()],
              'all': gutters + [splitter()]}[part]
    m = trimesh.util.concatenate(pieces)
    if verbose:
        q_t = P.Q_TRICKLE / 3600 * 1e-3 / (P.DROP_V * 1e-9)
        q_d = P.Q_DUMP / 3600 * 1e-3 / (P.DROP_V * 1e-9)
        _, nstep = helix(0)
        print(f"  spirals          {N_SPIRAL} x {TURNS} turns, bore {BORE} mm")
        print(f"  stepped helix    {nstep} level runs of "
              f"{tot_len/N_SPIRAL/nstep*LEVEL_FRAC:.0f} mm, {P.STEP:.0f} mm drops")
        print(f"  gutter length    {tot_len/1000:.2f} m total ({tot_len/N_SPIRAL:.0f} mm each)")
        print(f"  drip notches     {tot_notch} at {NOTCH_P:.0f} mm pitch")
        print(f"  per-notch rate   {q_t/tot_notch:.2f}/s trickle   {q_d/tot_notch:.2f}/s dump "
              f"(jet limit {P.JET_LIMIT:.1f})")
        print(f"  per-spiral flow  {P.Q_DUMP/3600*1e3/N_SPIRAL:.2f} mL/s at dump "
              f"(pipe jet limit 1.67) "
              f"{'OK' if P.Q_DUMP/3600*1e3/N_SPIRAL < 1.67 else '*** STREAMS ***'}")
        print(f"  intake           {BORE*4:.0f} mm ID -> 2 x {BORE*2*math.sqrt(2)/2*2:.1f}"
              f" -> 4 x {BORE*2:.1f} -> 8 x {BORE*math.sqrt(2):.1f} -> 16 x {BORE}")
        print(f"  mesh             {len(m.faces):,} faces  extents {np.round(m.extents,1)}")
    return m

if __name__ == "__main__":
    print("drip spiral + splitter:")
    m = build()
    m.export("drip_spiral.stl")
    print("  wrote drip_spiral.stl")
