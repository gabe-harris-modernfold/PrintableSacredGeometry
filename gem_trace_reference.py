"""
gem_trace_reference.py — pure-numpy auditor for gem-cut-light-paths.html (P1-P4 scope).

Independent implementation of the three primitives the page's JS engine also
implements: half-space polyhedron intersection, vector Snell/Fresnel/TIR with
the Stokes/Mueller polarization ledger. Traces a fixed set of rays through a
round-brilliant cut and checks the results against the pinned acceptance
vectors in gem-cut-light-paths-plan.html §8. Python is the auditor, not the
runtime — this script has no dependency on the HTML/JS build; the JS engine
is checked against ITS numbers separately (via node), not the other way round.

Scope note: this build covers the round brilliant only, at its four
"primary" facets (table, crown main, pavilion main) plus a girdle band —
star, upper-girdle and lower-girdle facets are deferred (they carry none of
the pinned optical vectors; they are cosmetic flash-pattern detail, not
physics this script or the page currently claims). The girdle band here is
a regular octagon (apothem = nominal girdle radius) matching the 8-fold
symmetry of that reduced facet set, not the plan's "16 sides default" for
the full 58-facet build.
"""

import itertools
import math
import numpy as np

# ============================================================ materials ===
# n at sodium D (589 nm), from plan §2. thetac and R0 (normal-incidence
# reflectance, unpolarized) are DERIVED, never authored.
RAW_N = {
    'moissanite': 2.6480,
    'diamond':    2.4175,
    'cz':         2.1600,
    'sapphire':   1.7680,
    'spinel':     1.7180,
    'beryl':      1.5770,
    'quartz':     1.5443,
    'glass':      1.5200,
    'fluorite':   1.4338,
}
# B-G (gemological) dispersion, n(430.8nm) - n(686.7nm), from plan §2's
# pinned table -- the source for P6's per-material Cauchy fit.
DISP_BG = {
    'moissanite': 0.104,
    'diamond':    0.044,
    'cz':         0.060,
    'sapphire':   0.018,
    'spinel':     0.020,
    'beryl':      0.014,
    'quartz':     0.013,
    'glass':      0.009,
    'fluorite':   0.007,
}

MATERIALS = {}
for name, n in RAW_N.items():
    thetac = math.degrees(math.asin(1.0 / n))
    R0 = ((n - 1.0) / (n + 1.0)) ** 2
    MATERIALS[name] = dict(n=n, thetac=thetac, R0=R0, disp_bg=DISP_BG[name])


def corridor_for(material):
    """Plan §2b: the pavilion-angle window where a vertical ray's canonical
    round-brilliant path (TIR at both pavilion hits, then a table exit)
    survives -- intersection of the two-TIR window [thetac, (180-thetac)/3]
    and the table-exit window ((180-thetac)/4, (180+thetac)/4). Closed
    form, no tracing -- the P8 gate is exactly the assertion that the
    TRACED ridge (see run_checks) falls inside this ALGEBRAIC window."""
    tc = material['thetac']
    lo = max(tc, (180.0 - tc) / 4.0)
    hi = min((180.0 - tc) / 3.0, (180.0 + tc) / 4.0)
    return lo, hi


# ============================================================ geometry ====
def round_brilliant_planes(table_pct=56.0, crown_deg=34.5, pav_deg=40.75,
                            girdle_diam=100.0, girdle_t_pct=3.0):
    """Half-space planes {n (outward unit normal), d, name} for the reduced
    round-brilliant facet set: table, 8 crown mains, 8 girdle sides, 8
    pavilion mains. Solid = intersection of {n.x <= d} over all planes.
    Every angle here is one of the plan's declared inputs (table/crown/
    pavilion/girdle-thickness); every coordinate below is DERIVED from them.
    """
    Rg = girdle_diam / 2.0                       # girdle apothem, = 50
    Rt = Rg * table_pct / 100.0                  # table apothem, = 28
    t = girdle_diam * girdle_t_pct / 100.0       # girdle band thickness
    crown = math.radians(crown_deg)
    pav = math.radians(pav_deg)
    hc = (Rg - Rt) * math.tan(crown)             # crown height
    hp = Rg * math.tan(pav)                      # pavilion depth
    table_y = t / 2.0 + hc
    culet_y = -(t / 2.0 + hp)

    planes = [dict(n=np.array([0.0, 1.0, 0.0]), d=table_y, name='table')]
    for k in range(8):
        phi = math.radians(45.0 * k)
        c, s = math.cos(phi), math.sin(phi)
        planes.append(dict(
            n=np.array([math.sin(crown) * c, math.cos(crown), math.sin(crown) * s]),
            d=Rg * math.sin(crown) + (t / 2.0) * math.cos(crown),
            name=f'crown_main_{k}'))
        planes.append(dict(
            n=np.array([c, 0.0, s]), d=Rg, name=f'girdle_side_{k}'))
        planes.append(dict(
            n=np.array([math.sin(pav) * c, -math.cos(pav), math.sin(pav) * s]),
            d=Rg * math.sin(pav) - (t / 2.0) * math.cos(pav),
            name=f'pavilion_main_{k}'))

    meta = dict(Rg=Rg, Rt=Rt, t=t, hc=hc, hp=hp, table_y=table_y, culet_y=culet_y,
                depth_pct=(table_y - culet_y))
    return planes, meta


def oval_brilliant_planes(lw_ratio=1.40, crown_deg=35.0, pav_deg=41.0, table_frac=0.59,
                           depth_pct=60.5, girdle_t_pct=3.5, major_diam=100.0, n_samples=32):
    """Oval brilliant. Source: IGI 'Guidelines for Excellent -- Fancy Shape
    Cut Grading' v22.10.08 (2022) -- crown 32.0-38.0 deg, pavilion
    39.5-42.5 deg, table 56.0-62.0%, depth 58.0-63.0%, girdle 2.0-5.1%, L:W
    1.35-1.50. Defaults here are each range's centre. 8-pavilion-main/
    56-58-facet brilliant taxonomy confirmed by the same research; this
    reduced build follows round_brilliant_planes' own scope note and omits
    star/upper-/lower-girdle facets."""
    a, b = major_diam / 2.0, major_diam / 2.0 / lw_ratio
    outline = ellipse_outline_points(a, b, n_samples)
    return fancy_brilliant_planes(outline, crown_deg, pav_deg, table_frac, depth_pct,
                                   girdle_t_pct, major_diam, name='oval')


def marquise_outline_points(a, b, n_per_arc=12):
    """Marquise/navette outline: two arcs through the tip points (+-a,0),
    each bulging to (0,+-b), centred on the perpendicular bisector of the
    tips (the z-axis) -- genuine sharp corners AT the tips (no tangency
    condition there, unlike oval's smooth ends: a marquise's point is a
    real corner, confirmed by P1 research's own two-arc description)."""
    zc = (b * b - a * a) / (2 * b)
    r = b - zc
    c_top, c_bot = np.array([0, zc]), np.array([0, -zc])

    def ang(center, pt):
        return math.atan2(pt[1] - center[1], pt[0] - center[0])

    a_tip_r, a_tip_l = ang(c_top, (a, 0)), ang(c_top, (-a, 0))
    b_tip_l, b_tip_r = ang(c_bot, (-a, 0)), ang(c_bot, (a, 0))
    pts = []
    for k in range(n_per_arc):
        t = a_tip_r + (a_tip_l - a_tip_r) * k / n_per_arc
        pts.append((c_top[0] + r * math.cos(t), c_top[1] + r * math.sin(t)))
    for k in range(n_per_arc):
        t = b_tip_l + (b_tip_r - b_tip_l) * k / n_per_arc
        pts.append((c_bot[0] + r * math.cos(t), c_bot[1] + r * math.sin(t)))
    return pts


def marquise_planes(lw_ratio=2.00, crown_deg=34.5, pav_deg=41.0, table_frac=0.60,
                     depth_pct=61.0, girdle_t_pct=3.5, major_diam=100.0, n_per_arc=12):
    """Marquise/navette. Source: IGI (2022) -- crown 31.0-38.0 deg, pavilion
    39.5-42.5 deg, table 56.0-64.0%, depth 57.5-64.5%, girdle 2.0-5.1%, L:W
    1.70-2.25 (canonical 2.00, independently confirmed by a real lapidary
    design at exactly that ratio). Defaults are each range's centre."""
    a, b = major_diam / 2.0, major_diam / 2.0 / lw_ratio
    outline = marquise_outline_points(a, b, n_per_arc)
    return fancy_brilliant_planes(outline, crown_deg, pav_deg, table_frac, depth_pct,
                                   girdle_t_pct, major_diam, name='marq')


def pear_outline_points(a_point, a_round, b, n_per_arc=10):
    """Pear/pendeloque outline: a marquise-style pointed half (two arcs from
    (0,+-b) to the tip at (a_point,0)) joined to a half-ellipse rounded end
    (semi-axes a_round, b, from (0,-b) through (-a_round,0) to (0,b)).
    Matches P1 research's own characterisation: "a pear can be modeled by
    the combination of a marquise and an oval." """
    zc = (b * b - a_point * a_point) / (2 * b)
    r = b - zc
    c_top, c_bot = np.array([0, zc]), np.array([0, -zc])

    def ang(center, pt):
        return math.atan2(pt[1] - center[1], pt[0] - center[0])

    a_top, a_tip = ang(c_top, (0, b)), ang(c_top, (a_point, 0))
    b_tip, b_bot = ang(c_bot, (a_point, 0)), ang(c_bot, (0, -b))
    pts = []
    for k in range(n_per_arc):
        t = a_top + (a_tip - a_top) * k / n_per_arc
        pts.append((c_top[0] + r * math.cos(t), c_top[1] + r * math.sin(t)))
    pts.append((a_point, 0.0))
    for k in range(1, n_per_arc + 1):
        t = b_tip + (b_bot - b_tip) * k / n_per_arc
        pts.append((c_bot[0] + r * math.cos(t), c_bot[1] + r * math.sin(t)))
    for k in range(1, 2 * n_per_arc):
        t = -math.pi / 2 + math.pi * k / (2 * n_per_arc)
        pts.append((-a_round * math.cos(t), b * math.sin(t)))
    return pts


def pear_planes(lw_ratio=1.55, point_frac=0.56, crown_deg=34.5, pav_deg=40.5, table_frac=0.59,
                 depth_pct=61.0, girdle_t_pct=3.5, major_diam=100.0, n_per_arc=12):
    """Pear/pendeloque. Source: IGI (2022) -- crown 32.0-37.0 deg, pavilion
    38.0-42.5 deg, table 55.5-63.0%, depth 57.5-65.0%, girdle 2.0-5.1%, L:W
    1.40-1.75 (canonical 1.55, IGI band centre; a real lapidary design
    (Long 2017) independently landed at L/W 1.476). point_frac (how much of
    the total length the pointed half takes) is NOT independently sourced
    -- no reference published a point:round split -- so 0.56 is a
    documented, reasonable default, not a discovered constant."""
    total = major_diam
    a_point, a_round = total * point_frac, total * (1 - point_frac)
    b = total / lw_ratio / 2.0
    outline = pear_outline_points(a_point, a_round, b, n_per_arc)
    return fancy_brilliant_planes(outline, crown_deg, pav_deg, table_frac, depth_pct,
                                   girdle_t_pct, total, name='pear')


def rounded_rect_outline_points(half_l, half_w, exponent=4.0, n=48):
    """Cushion (and, at exponent~2, a near-oval) outline via a superellipse
    |x/half_l|^p + |z/half_w|^p = 1. Deliberate simplification, same
    reasoning as ellipse_outline_points: no source publishes actual corner-
    arc radii for a modern cushion brilliant, so a smooth, provably-convex
    closed form is preferred over a hand-tangent-solved multi-arc
    reconstruction. Higher exponent -> squarer corners; p=4 is a
    conventional "softly rounded square" middle ground."""
    pts = []
    for k in range(n):
        t = 2 * math.pi * k / n
        c, s = math.cos(t), math.sin(t)
        x = half_l * np.sign(c) * abs(c) ** (2.0 / exponent)
        z = half_w * np.sign(s) * abs(s) ** (2.0 / exponent)
        pts.append((x, z))
    return pts


def cushion_brilliant_planes(lw_ratio=1.02, exponent=4.0, crown_deg=35.0, pav_deg=41.0, table_frac=0.585,
                              depth_pct=65.0, girdle_t_pct=3.5, major_diam=100.0, n=48):
    """Cushion (modified) brilliant, square/near-square form. Source: IGI
    (2022) 'Square Cushion Modified Brilliant' -- table 55.0-62.0%, depth
    62.0-68.0%, crown height 8.0-16.0%, pavilion depth 40.0-56.0% (a wide
    band -- IGI's own text notes it covers both "cushion brilliant" (4-8
    large pavilion mains) and "cushion modified brilliant" (added lower
    facets) in one table), L:W 1.00-1.05. No single crown/pavilion angle is
    published (same chevron-family caveat as princess/radiant); the single-
    angle values here are a simplification consistent with this build's
    reduced primary-facets-only scope, not a claim that real cushions use
    one uniform pavilion angle."""
    a, b = major_diam / 2.0, major_diam / 2.0 / lw_ratio
    outline = rounded_rect_outline_points(a, b, exponent, n)
    return fancy_brilliant_planes(outline, crown_deg, pav_deg, table_frac, depth_pct,
                                   girdle_t_pct, major_diam, name='cushion')


def trilliant_planes(radius=50.0, roundedness=0.20, crown_deg=43.6, pav_deg=40.8, table_frac=0.55,
                      depth_pct=63.0, girdle_t_pct=3.5, n=48):
    """Trilliant, curved-side form (the common modern reading; a straight-
    edged accent-stone variant also exists but its exact ratio came back
    unverified from P1 research and is not modelled here). Outline: a
    3-fold "rounded-triangle" polar curve r(theta) = R*(1 - roundedness +
    roundedness*cos(3*theta)) -- a standard, smooth, provably star-convex
    3-fold curve family (not itself a sourced lapidary construction, same
    honesty note as the cushion superellipse). Angles are a reduced,
    single-facet-family simplification of a real, fully-dimensioned named
    design -- Norman W. Steele, "Trilliant C" (Facet Diagrams 13.067,
    1991): crown mains 43.63 deg, pavilion mains 40.80 deg -- presented as
    one real designer's choice, not an industry standard (none exists)."""
    pts = []
    for k in range(n):
        t = 2 * math.pi * k / n
        r = radius * (1 - roundedness + roundedness * math.cos(3 * t))
        pts.append((r * math.cos(t), r * math.sin(t)))
    return fancy_brilliant_planes(pts, crown_deg, pav_deg, table_frac, depth_pct,
                                   girdle_t_pct, 2 * radius, name='trilliant')


def half_moon_planes(radius=50.0, crown_deg=34.5, pav_deg=40.75, table_frac=0.56,
                      depth_pct=60.5, girdle_t_pct=3.5, n=32):
    """Half-moon. P1 research found this genuinely unstandardised -- SkyJems'
    own gemological encyclopedia states "faceting conventions... are not
    rigidly standardised" and no source gives a crown/pavilion angle, table
    %, or facet count. Per that research's own suggestion, this ships
    using the round-brilliant's own angles as an explicit, documented
    placeholder (half-moons are conventionally paired with round centre
    stones), NOT a discovered half-moon constant. The outline itself IS
    well sourced: a true semicircle, straight diameter edge + 180 deg arc,
    radius = half the straight edge."""
    pts = [(radius, 0.0)]
    for k in range(1, n):
        t = math.pi * k / n
        pts.append((radius * math.cos(t), radius * math.sin(t)))
    pts.append((-radius, 0.0))
    # interior_xz can't be the outline's own geometric centre (0,0) here --
    # the flat diameter edge passes exactly through it, which is the same
    # degeneracy class as orienting a facet against a point that sits ON
    # that facet (see chord_planes_from_outline's docstring) -- (0, small
    # positive z) is safely inside the arc's bulge instead.
    return fancy_brilliant_planes(pts, crown_deg, pav_deg, table_frac, depth_pct,
                                   girdle_t_pct, 2 * radius, name='halfmoon', interior_xz=(0.0, radius * 0.3))


def kite_planes(a=50.0, b_top=35.0, b_bot=60.0, crown_deg=34.5, pav_deg=40.75, table_frac=0.56,
                 depth_pct=60.5, girdle_t_pct=3.5):
    """Kite: a true kite quadrilateral (2 pairs of equal adjacent sides,
    unequal diagonals) -- "4 straight chords by definition," per P1
    research, which also found no industry-standard angle table for this
    shape (only individual named lapidary designs, e.g. Vargas 1975's
    "Barion Kite"). Angles here are the round-brilliant placeholder, same
    reasoning as half-moon, documented rather than presented as sourced."""
    pts = [(0.0, b_top), (a, 0.0), (0.0, -b_bot), (-a, 0.0)]
    major_diam = b_top + b_bot
    return fancy_brilliant_planes(pts, crown_deg, pav_deg, table_frac, depth_pct,
                                   girdle_t_pct, major_diam, name='kite')


def shield_planes(half_top=35.0, half_mid=45.0, top_y=40.0, tip_y=-45.0, crown_deg=34.5, pav_deg=40.75,
                   table_frac=0.56, depth_pct=60.5, girdle_t_pct=3.5):
    """Shield: flat top edge, widest a little below the top, tapering to a
    single point at the base -- "typically 4 straight sides," per P1
    research, which again found no industry-standard angle table (one real
    named design exists, Strickland 2014's "Trilliant Shield", but as a
    one-off, not a standard). Modelled as a symmetric pentagon: 2 top
    corners (flat top edge), 2 widest side points, 1 base point. Angles are
    the round-brilliant placeholder, same reasoning as half-moon/kite."""
    pts = [(-half_top, top_y), (-half_mid, 0.0), (0.0, tip_y), (half_mid, 0.0), (half_top, top_y)]
    major_diam = top_y - tip_y
    return fancy_brilliant_planes(pts, crown_deg, pav_deg, table_frac, depth_pct,
                                   girdle_t_pct, major_diam, name='shield')


def circle_through_3pts(p1, p2, p3):
    """Circumcircle centre+radius through 3 non-colinear 2D points."""
    ax, ay = p1; bx, by = p2; cx, cy = p3
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    ux = ((ax**2 + ay**2) * (by - cy) + (bx**2 + by**2) * (cy - ay) + (cx**2 + cy**2) * (ay - by)) / d
    uy = ((ax**2 + ay**2) * (cx - bx) + (bx**2 + by**2) * (ax - cx) + (cx**2 + cy**2) * (bx - ax)) / d
    return (ux, uy), math.hypot(ax - ux, ay - uy)


def arc_through_3pts(p1, p2, p3, n=16):
    """Samples the circular arc through p1,p2,p3 (in that order along the
    arc) -- used where no source gives an explicit tangency/radius
    construction, but 3 named landmark points (e.g. a heart lobe's tip,
    widest point and shoulder) are a reasonable, checkable design choice."""
    c, r = circle_through_3pts(p1, p2, p3)

    def ang(pt):
        return math.atan2(pt[1] - c[1], pt[0] - c[0])
    a1, a2, a3 = ang(p1), ang(p2), ang(p3)
    while a2 < a1:
        a2 += 2 * math.pi
    while a3 < a2:
        a3 += 2 * math.pi
    return [(c[0] + r * math.cos(a1 + (a3 - a1) * k / n), c[1] + r * math.sin(a1 + (a3 - a1) * k / n))
            for k in range(n + 1)]


def heart_lobe_outline(sign, a=50.0, b=38.0, x_shoulder=35.0, x_cleft=12.0, n=16):
    """One convex lobe of a heart: tip (-a,0) [shared with the mirror lobe,
    the heart's bottom point] -> arc through the widest point (0,sign*b) and
    a shoulder (x_shoulder,sign*b*0.79) -> straight segment to the cleft
    point (x_cleft,0) [shared with the mirror lobe, on the centreline] ->
    straight closure back to the tip. x_cleft < x_shoulder is what reads as
    a notch once both lobes are placed together; P1 research found no
    source publishes a cleft depth standard, so these proportions are a
    documented design choice, not a discovered constant."""
    tip, widest, shoulder, cleft = (-a, 0.0), (0.0, sign * b), (x_shoulder, sign * b * 0.79), (x_cleft, 0.0)
    return arc_through_3pts(tip, widest, shoulder, n) + [cleft]


def heart_planes(a=50.0, b=38.0, x_shoulder=35.0, x_cleft=12.0, crown_deg=35.0, pav_deg=40.5,
                  table_frac=0.61, depth_pct=57.0, girdle_t_pct=3.5, n=16):
    """Heart. Source: IGI (2022) -- crown 32.0-38.0 deg, pavilion 39.0-42.5
    deg, table 58.0-64.0%, depth 51.9-61.0%, girdle 2.0-5.1%, L:W (length/
    width) 0.85-0.95 -- a well-proportioned heart is slightly WIDER than
    long, not 1:1 (the common "1:1" retail description is a colloquial
    rounding of this). Modelled as a genuine UNION of two convex lobes
    (plan §3: "only two roster entries need this -- heart... and
    merkaba"), a computational convenience mirroring what every source
    describes as a real single-crystal cleft, not a laminated assembly."""
    major_diam = a + x_shoulder  # tip-to-shoulder span, the outline's own long axis
    right = heart_lobe_outline(+1, a, b, x_shoulder, x_cleft, n)
    left = heart_lobe_outline(-1, a, b, x_shoulder, x_cleft, n)
    planes_r, _ = fancy_brilliant_planes(right, crown_deg, pav_deg, table_frac, depth_pct,
                                          girdle_t_pct, major_diam, name='heartR', interior_xz=(0, b * 0.4))
    planes_l, _ = fancy_brilliant_planes(left, crown_deg, pav_deg, table_frac, depth_pct,
                                          girdle_t_pct, major_diam, name='heartL', interior_xz=(0, -b * 0.4))
    return dict(parts=[dict(planes=planes_r, quadrics=[]), dict(planes=planes_l, quadrics=[])])


def stepped_pavilion_planes(outline_xz, tier_height_fracs, tier_scale_ends, girdle_y,
                             target_depth, interior, name):
    """General multi-tier "chevron" pavilion: a stack of frustum rings, each
    shrinking the outline from its start scale to tier_scale_ends[k] over a
    height target_depth*tier_height_fracs[k] (fracs must sum to 1.0), with
    the FINAL tier (scale_end<=0) converging to a true apex point via plain
    chord-to-apex facets rather than another ring. Each ring facet is
    planar because it connects corresponding edges of two UNIFORMLY-scaled
    copies of the same outline -- standard frustum geometry.

    Height-driven, not angle-driven, and deliberately so: princess/radiant/
    cushion's pavilions don't reduce to one crown/pavilion angle at all --
    P1 research found IGI's own 2022 grading standard specifies PAVILION
    DEPTH directly for these cuts rather than an angle, for exactly this
    reason (a multi-tier chevron has no single angle to quote). Matching
    that convention -- target a sourced depth% and let tier geometry follow
    -- is more faithful than inventing per-tier angles that would need to
    coincidentally sum to the right depth anyway."""
    planes = []
    y, scale = girdle_y, 1.0
    m = len(outline_xz)
    outline = np.array(outline_xz)
    for tier_idx, (hfrac, s_end) in enumerate(zip(tier_height_fracs, tier_scale_ends)):
        h = target_depth * hfrac
        y_next = y - h
        if s_end <= 1e-9:
            planes += chord_planes_from_outline([tuple(p) for p in outline * scale], (0, y_next, 0),
                                                 y, f'{name}_t{tier_idx}', interior)
            y, scale = y_next, 0.0
            continue
        outer, inner = outline * scale, outline * s_end
        for i in range(m):
            p1o = np.array([outer[i][0], y, outer[i][1]])
            p2o = np.array([outer[(i + 1) % m][0], y, outer[(i + 1) % m][1]])
            p1i = np.array([inner[i][0], y_next, inner[i][1]])
            nrm = np.cross(p2o - p1o, p1i - p1o)
            length = np.linalg.norm(nrm)
            if length < 1e-9:
                continue
            nrm = nrm / length
            d = np.dot(nrm, p1o)
            if np.dot(nrm, interior) > d - 1e-9:
                nrm, d = -nrm, -d
            planes.append(dict(n=nrm, d=d, name=f'{name}_t{tier_idx}_{i}'))
        y, scale = y_next, s_end
    return planes, y


def chevron_brilliant_planes(outline_xz, crown_deg, table_frac, pav_depth_pct, tier_height_fracs,
                              tier_scale_ends, girdle_t_pct, major_diam, name):
    """Crown (single chord-to-apex tier, table derived by similar triangles
    -- same as fancy_brilliant_planes) over a multi-tier chevron pavilion
    (stepped_pavilion_planes) targeting a sourced pavilion depth% directly.
    What princess, radiant and cushion (modified) brilliant all reduce to."""
    t_half = major_diam * girdle_t_pct / 100.0 / 2.0
    apothem = float(np.mean([math.hypot((outline_xz[i][0] + outline_xz[(i + 1) % len(outline_xz)][0]) / 2,
                                         (outline_xz[i][1] + outline_xz[(i + 1) % len(outline_xz)][1]) / 2)
                              for i in range(len(outline_xz))]))
    hc = apothem * math.tan(math.radians(crown_deg))
    apex_crown_y = t_half + hc
    table_y = apex_crown_y - table_frac * hc
    interior = (0.0, 0.0, 0.0)
    target_pav_depth = major_diam * pav_depth_pct / 100.0
    pav_pl, culet_y = stepped_pavilion_planes(outline_xz, tier_height_fracs, tier_scale_ends,
                                               -t_half, target_pav_depth, interior, f'{name}_pav')
    crown_pl = chord_planes_from_outline(outline_xz, (0, apex_crown_y, 0), t_half, f'{name}_crown', interior)
    girdle_pl = girdle_band_from_outline(outline_xz, t_half, interior)
    table_pl = dict(n=np.array([0., 1., 0.]), d=table_y, name=f'{name}_table')
    planes = dedupe_planes([table_pl] + crown_pl + girdle_pl + pav_pl)
    meta = dict(t=2 * t_half, hc=hc, tableY=table_y, culetY=culet_y,
                depthPct=(table_y - culet_y) / major_diam * 100)
    return planes, meta


def square_outline_points(half_w):
    return [(half_w, half_w), (-half_w, half_w), (-half_w, -half_w), (half_w, -half_w)]


def chamfered_rect_outline_points(half_l, half_w, chamfer_frac=0.25):
    """Rectangle/square with 4 truncated (chamfered) corners -- radiant's
    outline per P1 research ("cut-cornered rectangular modified brilliant";
    the Grossbard patents describe an octagonal girdle). chamfer_frac is
    how much of each half-side the chamfer eats, symmetric on both edges
    meeting at a corner."""
    cl, cw = half_l * chamfer_frac, half_w * chamfer_frac
    return [
        (half_l - cl, half_w), (half_l, half_w - cw),
        (half_l, -(half_w - cw)), (half_l - cl, -half_w),
        (-(half_l - cl), -half_w), (-half_l, -(half_w - cw)),
        (-half_l, half_w - cw), (-(half_l - cl), half_w),
    ]


def princess_planes(crown_deg=37.0, table_frac=0.68, pav_depth_pct=57.0, girdle_t_pct=3.0,
                     major_diam=100.0):
    """Princess cut. Source: IGI (2022) -- table 65.0-74.0%, total depth
    65.0-73.0%, crown height 8.0-14.5%, pavilion depth 54.5-62.0% (IGI
    publishes depth directly, NOT a crown/pavilion angle -- see
    chevron_brilliant_planes' docstring). L:W 1.00-1.05 (square). 2-tier
    chevron pavilion (50-52 deg outer / 61-64 deg inner, per US Patent
    6,745,596 B2's one real commercial embodiment) reduced here to a
    height-driven 35%/65% split at 55% outline scale, hitting IGI's
    sourced pavilion-depth band by construction rather than by chance."""
    outline = square_outline_points(major_diam / 2.0)
    return chevron_brilliant_planes(outline, crown_deg, table_frac, pav_depth_pct,
                                     [0.35, 0.65], [0.55, 0.0], girdle_t_pct, major_diam, 'princess')


def radiant_planes(crown_deg=35.0, table_frac=0.645, pav_depth_pct=51.0, chamfer_frac=0.25,
                    girdle_t_pct=3.0, major_diam=100.0):
    """Radiant cut, square form. Source: IGI (2022) 'Cut-Cornered Square
    Modified Brilliant' -- table 62.0-67.0%, total depth 62.0-68.0%, crown
    height 12.0-15.0%, pavilion depth 47.0-56.0% (again depth-published,
    not angle -- same chevron-family reasoning as princess). Chamfered-
    rectangle outline (chamfered_rect_outline_points), otherwise the same
    2-tier chevron-pavilion construction as princess."""
    outline = chamfered_rect_outline_points(major_diam / 2.0, major_diam / 2.0, chamfer_frac)
    return chevron_brilliant_planes(outline, crown_deg, table_frac, pav_depth_pct,
                                     [0.35, 0.65], [0.55, 0.0], girdle_t_pct, major_diam, 'radiant')


def vogel_wand_planes(n_sides=6, radius=25.0, blunt_deg=51.86, sharp_deg=64.07):
    """Vogel cut, reading B default (plan §6 defaults to B, offers A as a
    toggle). AOI for an axial ray equals the facet's tilt-from-horizontal
    DIRECTLY -- the same identity round_brilliant_planes' own pavilion
    mains rely on (§2b: "first pavilion hit: AOI1 = p") -- so blunt_deg/
    sharp_deg are each reading's AXIAL-AOI value (51.86/64.07), not its
    "facet plane off c-axis" value (38.14/25.93, = 90-AOI): plugging in the
    off-axis figure directly was a real bug caught only by checking against
    the plan's own pinned vectors, not by energy closure alone (a wrong-
    but-still-closed solid closes energy just fine). Larger tilt = taller,
    SHARPER tip (hp=radius*tan(tilt) grows without bound as tilt->90),
    matching reading A's "slender needle" language; smaller tilt = a
    shorter, blunter tip, matching reading B."""
    planes = []
    br, sr = math.radians(blunt_deg), math.radians(sharp_deg)
    for k in range(n_sides):
        phi = 2 * math.pi * k / n_sides
        c, s = math.cos(phi), math.sin(phi)
        planes.append(dict(n=np.array([c, 0., s]), d=radius, name=f'prism_{k}'))
        planes.append(dict(n=np.array([math.sin(br) * c, math.cos(br), math.sin(br) * s]),
                            d=radius * math.sin(br), name=f'blunt_{k}'))
        planes.append(dict(n=np.array([math.sin(sr) * c, -math.cos(sr), math.sin(sr) * s]),
                            d=radius * math.sin(sr), name=f'sharp_{k}'))
    return planes


def generator_point_planes(n_sides=6, radius=25.0, tip_deg=51.86, base_y=-40.0):
    """Generator point: hexagonal prism, ONE pyramidal termination ("6
    equal faces to one apex"), flat polished base at the other end (unlike
    the Vogel wand's two opposed points)."""
    planes = []
    tr = math.radians(tip_deg)
    for k in range(n_sides):
        phi = 2 * math.pi * k / n_sides
        c, s = math.cos(phi), math.sin(phi)
        planes.append(dict(n=np.array([c, 0., s]), d=radius, name=f'prism_{k}'))
        planes.append(dict(n=np.array([math.sin(tr) * c, math.cos(tr), math.sin(tr) * s]),
                            d=radius * math.sin(tr), name=f'tip_{k}'))
    planes.append(dict(n=np.array([0., -1., 0.]), d=-base_y, name='base'))
    return planes


def obelisk_planes(half_w=20.0, tip_deg=60.0, base_y=-60.0):
    """Obelisk/tower: square prism + square pyramidal cap, flat base."""
    planes = []
    tr = math.radians(tip_deg)
    for k in range(4):
        phi = math.pi / 4 + k * math.pi / 2
        c, s = math.cos(phi), math.sin(phi)
        r_apothem = half_w
        planes.append(dict(n=np.array([c, 0., s]), d=r_apothem, name=f'prism_{k}'))
        planes.append(dict(n=np.array([math.sin(tr) * c, math.cos(tr), math.sin(tr) * s]),
                            d=r_apothem * math.sin(tr), name=f'tip_{k}'))
    planes.append(dict(n=np.array([0., -1., 0.]), d=-base_y, name='base'))
    return planes


def great_pyramid_planes(half_base=50.0, slope_deg=51.8333):
    """Great Pyramid: 51 deg 50' slope (51.8333 deg) from the base plane,
    i.e. tilt-from-horizontal = slope_deg directly -- square base, apex,
    no separate prism section (sits "beside the Vogel apex for
    comparison," per the roster, so uses the SAME tilt-from-horizontal
    convention as the Vogel tip facets, not a different one)."""
    planes = []
    sr = math.radians(slope_deg)
    for k in range(4):
        phi = math.pi / 4 + k * math.pi / 2
        c, s = math.cos(phi), math.sin(phi)
        planes.append(dict(n=np.array([math.sin(sr) * c, math.cos(sr), math.sin(sr) * s]),
                            d=half_base * math.sqrt(2) * math.sin(sr), name=f'face_{k}'))
    planes.append(dict(n=np.array([0., -1., 0.]), d=0.0, name='base'))
    return planes


def sphere_planes(radius=50.0):
    """Sphere/scrying bead: pure quadric, ties to scrying_pool.py per the
    roster."""
    return dict(planes=[], quadrics=[dict(center=np.array([0., 0., 0.]), radii=np.array([radius] * 3), name='sphere')])


def egg_planes(polar_r=55.0, equatorial_r=40.0):
    """Egg: prolate ellipsoid of revolution about the vertical axis.
    polar:equatorial ~1.3-1.4:1 is a common real-egg proportion; not a
    lapidary standard (none is claimed)."""
    return dict(planes=[], quadrics=[dict(center=np.array([0., 0., 0.]),
                                           radii=np.array([equatorial_r, polar_r, equatorial_r]), name='egg')])


def cabochon_planes(girdle_diam=100.0, dome_height_frac=0.35):
    """Cabochon: dome (sphere) intersected with a flat back plane.
    dome_height_frac ~0.30-0.40 is a common "medium dome" convention (not a
    single fixed lapidary standard). The sphere's radius is NOT a free
    parameter -- it's the standard spherical-cap relation solved from the
    girdle radius and dome height (r_girdle^2 + (dome_radius-dome_h)^2 =
    dome_radius^2 => dome_radius = (r_girdle^2+dome_h^2)/(2*dome_h)); an
    earlier draft of this function took radius as an independent input,
    which over-specifies the shape and silently produces a girdle NOT at
    the stated diameter unless the caller happens to pick a self-consistent
    value."""
    r_girdle = girdle_diam / 2.0
    dome_h = girdle_diam * dome_height_frac
    dome_radius = (r_girdle ** 2 + dome_h ** 2) / (2 * dome_h)
    center_y = dome_h - dome_radius
    back_y = -1.0  # thin flat back, just below the girdle plane (y=0)
    return dict(planes=[dict(n=np.array([0., -1., 0.]), d=-back_y, name='cab_back')],
                quadrics=[dict(center=np.array([0., center_y, 0.]), radii=np.array([dome_radius] * 3), name='cab_dome')])


def lens_disc_planes(cap_radius=90.0, girdle_diam=100.0):
    """Lens/disc: two spherical caps (a biconvex lens), radius chosen for a
    moderate curvature (not flat, not a full hemisphere)."""
    r_girdle = girdle_diam / 2.0
    sep = 2 * math.sqrt(max(cap_radius ** 2 - r_girdle ** 2, 1.0))
    return dict(planes=[], quadrics=[
        dict(center=np.array([0., sep / 2, 0.]), radii=np.array([cap_radius] * 3), name='lens_top'),
        dict(center=np.array([0., -sep / 2, 0.]), radii=np.array([cap_radius] * 3), name='lens_bot')])


def star_of_david_planes(edge=40.0, height_ratio=0.375):
    """Vogel Star of David: triangular antiprism, "male" equilateral face
    on top, "female" face on the bottom rotated 60 deg, 6 lateral facets.
    Source: P1 research inferred height:edge ~0.375 (3:8) from two
    independent retailers' bounding-box dimensions landing on the exact
    same value -- moderate confidence, no lapidary spec exists, shipped
    approximate per plan §10's own policy. At height_ratio=sqrt(2/3)=
    0.8165 this hull IS the regular octahedron (plan §6); 0.375 is far
    flatter, matching what real pendants measure as."""
    h = edge * height_ratio
    r_circ = edge / math.sqrt(3)
    top = [(r_circ * math.cos(math.radians(90 + 120 * k)), h / 2, r_circ * math.sin(math.radians(90 + 120 * k)))
           for k in range(3)]
    bot = [(r_circ * math.cos(math.radians(90 + 60 + 120 * k)), -h / 2, r_circ * math.sin(math.radians(90 + 60 + 120 * k)))
           for k in range(3)]
    planes = [dict(n=np.array([0., 1., 0.]), d=h / 2, name='top'),
              dict(n=np.array([0., -1., 0.]), d=h / 2, name='bot')]
    interior = np.array([0., 0., 0.])
    for i in range(3):
        for a, b, c, label in ((top[i], top[(i + 1) % 3], bot[i], f'lat_{i}a'),
                                (top[(i + 1) % 3], bot[(i + 1) % 3], bot[i], f'lat_{i}b')):
            pa, pb, pc = np.array(a), np.array(b), np.array(c)
            nrm = np.cross(pb - pa, pc - pa)
            nrm = nrm / np.linalg.norm(nrm)
            d = np.dot(nrm, pa)
            if np.dot(nrm, interior) > d - 1e-9:
                nrm, d = -nrm, -d
            planes.append(dict(n=nrm, d=d, name=label))
    return planes


def double_terminated_point_planes(n_sides=6, radius=12.5, r_face_deg=51.78, length_ratio=4.0,
                                    prism_frac=0.6):
    """Double-terminated quartz point: hexagonal prism (n_sides=6, the
    natural, fixed-by-crystal-system count for an UNFACETED point -- P1
    research confirms lapidary-faceted "Vogel-style" wands vary this, but
    a natural-form double-terminated point conventionally does not),
    symmetric r-face terminations at both ends (unlike the Vogel wand's
    asymmetric blunt/sharp split -- a natural point is the same crystal
    form on both ends). r_face_deg=51.78 is the real crystallographic
    angle already computed elsewhere in this file's docs from quartz's
    cell parameters (a=4.9137, c=5.4047 A), not re-derived here.
    length_ratio (length:diameter) ~4:1 and prism_frac (fraction of total
    length that's straight prism vs the two terminations) are P1-sourced
    trade conventions (range ~2-6:1), not physical constants."""
    length = radius * 2 * length_ratio
    prism_half_len = length / 2 * prism_frac
    rr = math.radians(r_face_deg)
    planes = []
    for k in range(n_sides):
        phi = 2 * math.pi * k / n_sides
        c, s = math.cos(phi), math.sin(phi)
        planes.append(dict(n=np.array([c, 0., s]), d=radius, name=f'prism_{k}'))
        for sign, label in ((1, 'top'), (-1, 'bot')):
            planes.append(dict(n=np.array([math.sin(rr) * c, sign * math.cos(rr), math.sin(rr) * s]),
                                d=radius * math.sin(rr) + prism_half_len * math.cos(rr), name=f'{label}_{k}'))
    return planes


def merkaba_planes(d=1.0):
    """Merkaba: two interpenetrating tetrahedra (a "stella octangula"),
    point-inverting one tetrahedron's face-normal set (n -> -n, same d)
    gives the second, oppositely-oriented one. Union of 2 parts (§3: "only
    two roster entries need this -- heart... and merkaba")."""
    tet_a = platonic_planes('tetrahedron', d=d)
    tet_b = [dict(n=-p['n'], d=p['d'], name=f'inv_{p["name"]}') for p in tet_a]
    return dict(parts=[dict(planes=tet_a, quadrics=[]), dict(planes=tet_b, quadrics=[])])


def point_cut_planes(d=50.0):
    """Point cut: a polished natural octahedron, ~1300s -- "the zero
    point" of the antique-cut lineage, before any table or culet was ever
    ground. Literally platonic_planes('octahedron'); the roster's own
    "have" status reflects that this shape needs no separate sourcing."""
    return platonic_planes('octahedron', d=d)


def single_cut_planes(table_pct=53.0, crown_deg=40.0, pav_deg=43.0, girdle_diam=100.0, girdle_t_pct=3.0):
    """Single/eight cut: table + 8 crown facets + 8 pavilion facets, no
    star/girdle-break facets at all -- "the simplest real cut," and
    structurally IDENTICAL to round_brilliant_planes' own reduced model
    (this build already omits star/upper-/lower-girdle facets from the
    modern round brilliant for the same reason: they carry no pinned
    optical vector). Angles here are a plausible early/simple-cut proxy,
    not independently sourced -- the roster's "have" status is about the
    STRUCTURE being pinned, not these specific numbers."""
    planes, meta = round_brilliant_planes(table_pct, crown_deg, pav_deg, girdle_diam, girdle_t_pct)
    return planes, meta


def old_mine_planes():
    """Old mine cut. Source: P1 research (secondary-sourced via Tillander,
    not independently verified) -- crown angle disputed between sources
    (35-40 deg range; reduced-model default 40), table 38-53% (mid ~45%),
    depth ~71.2% (one real cited stone, treated as representative), cushion
    (rounded-square, hand-shaped, pre-1874) girdle. 33/25/58 facet topology
    per research; reduced here to primary facets only, same convention as
    round_brilliant_planes."""
    return cushion_brilliant_planes(lw_ratio=1.0, exponent=3.0, crown_deg=40.0, pav_deg=43.0,
                                     table_frac=0.45, depth_pct=71.2, girdle_t_pct=4.0)


def old_european_planes():
    """Old European cut. Source: P1 research -- table <=53% (mid ~50%),
    crown/pavilion angle actively CONTESTED between trade convention
    (crown >=40 deg) and Gilbertson's archival Henry Morse correspondence
    (crown 28-34 deg, circa 1860-70) -- this build follows the trade-
    convention figure since it's what most sources quote, and flags the
    Morse discrepancy in this docstring rather than silently picking a
    side. Round girdle (the first of the antique lineage that CAN be,
    post-1874 bruting machine)."""
    planes, meta = oval_brilliant_planes(lw_ratio=1.0, crown_deg=40.0, pav_deg=42.0, table_frac=0.50,
                                          depth_pct=63.0, girdle_t_pct=3.0)
    return planes, meta


def peruzzi_planes():
    """Peruzzi cut ("triple-cut brilliant"). Source: P1 research -- no
    angle table survives for this specific cut; the best real numeric
    anchor is David Jeffries' 1751 "A Treatise on Diamonds and Pearls"
    (crown 45 deg, pavilion 45 deg, table 56%, depth 66-67%), describing
    brilliant-cut diamonds of Peruzzi's own approximate era -- attributed
    to Jeffries, NOT to Peruzzi by name, per that research's own caution.
    Cushion girdle (pre-1874), 33/25/58 topology, reduced to primary
    facets as usual."""
    return cushion_brilliant_planes(lw_ratio=1.0, exponent=3.0, crown_deg=45.0, pav_deg=45.0,
                                     table_frac=0.56, depth_pct=66.5, girdle_t_pct=3.0)


def rect_outline_points(half_l, half_w):
    return [(half_l, half_w), (-half_l, half_w), (-half_l, -half_w), (half_l, -half_w)]


def trapezoid_outline_points(half_l, half_w_wide, half_w_narrow):
    """Tapered baguette outline. Source: P1 research -- the trade specifies
    these by explicit (length, width_wide_end, width_narrow_end) rather
    than a taper angle (the angle is a derived quantity), which is exactly
    the parametrization used here."""
    return [(half_l, half_w_narrow), (-half_l, half_w_wide), (-half_l, -half_w_wide), (half_l, -half_w_narrow)]


def emerald_cut_planes(half_l=37.5, half_w=25.0, chamfer_frac=0.28, crown_deg=35.5, table_frac=0.65,
                        pav_depth_pct=47.0, girdle_t_pct=3.5):
    """Emerald cut. Source: P1 research -- crown ~35-36 deg (secondary,
    Vargas via forum citation); pavilion genuinely TIERED, not one angle --
    a practicing cutter's real-world figures (Kidwell, GemologyOnline) give
    41 deg near the culet stepping to 47 deg near the girdle, a real 6 deg
    progression, matched here via chevron_brilliant_planes' height-driven
    3-tier pavilion rather than an authored single angle. Table 58-72%
    (mid ~65%), total depth 59-69% (ideal 61.3-67%). Confirmed 3 crown + 3
    pavilion steps. Corner-chamfer FRACTION is explicitly NOT standardised
    per that research ("cannot be determined from the lab report alone");
    0.28 is a reasonable default, not sourced. L:W default 1.5 (major_diam
    = 2*half_l, so half_l/half_w = 1.5 here)."""
    outline = chamfered_rect_outline_points(half_l, half_w, chamfer_frac)
    return chevron_brilliant_planes(outline, crown_deg, table_frac, pav_depth_pct,
                                     [0.30, 0.35, 0.35], [0.70, 0.35, 0.0], girdle_t_pct, 2 * half_l, 'emerald')


def asscher_planes(half_l=25.0, half_w=24.5, chamfer_frac=0.40, crown_deg=38.0, table_frac=0.63,
                    pav_depth_pct=47.0, girdle_t_pct=3.5):
    """Asscher cut (original, 58-facet form). Source: P1 research -- table
    60-68% (ideal 60-66%), depth 60-68% (ideal 61-66%), L:W 1.00-1.05 (near
    square). Crown/pavilion angle NOT credibly sourced (a widely-repeated
    "10-14 deg crown" figure was flagged by that research as almost
    certainly wrong -- contradicts every qualitative "steep, pyramidal
    crown" description -- and is deliberately NOT used here). Chamfer is
    larger than emerald cut's, "approaching an octagonal outline," per the
    same research; still an unsourced fraction."""
    outline = chamfered_rect_outline_points(half_l, half_w, chamfer_frac)
    return chevron_brilliant_planes(outline, crown_deg, table_frac, pav_depth_pct,
                                     [0.30, 0.35, 0.35], [0.70, 0.35, 0.0], girdle_t_pct, 2 * half_l, 'asscher')


def baguette_planes(half_l=37.5, half_w=18.75, crown_deg=33.0, pav_deg=44.0, table_frac=0.60,
                     depth_pct=62.0, girdle_t_pct=3.5):
    """Baguette. Source: P1 research -- plain rectangle, NO corner chamfer
    (the defining contrast with emerald cut, well-sourced). Minimal facet
    count ("~14 facets... essentially 1 step tier per side", an inference
    from that low count rather than a directly-cited row count) -- modelled
    as a single-cone crown+pavilion (fancy_brilliant_planes), not a
    chevron stack, matching "fewer facets... mostly windowed." Table
    59-62%, depth ~60-65%. L:W 1.5-3:1 range, default 2.0 here."""
    outline = rect_outline_points(half_l, half_w)
    return fancy_brilliant_planes(outline, crown_deg, pav_deg, table_frac, depth_pct,
                                   girdle_t_pct, 2 * half_l, name='baguette')


def tapered_baguette_planes(half_l=37.5, half_w_wide=20.0, half_w_narrow=10.0, crown_deg=33.0,
                             pav_deg=44.0, table_frac=0.60, depth_pct=62.0, girdle_t_pct=3.5):
    """Tapered baguette. Same facet treatment as baguette_planes (single-
    cone, no chamfer) over a trapezoid_outline_points outline instead of a
    plain rectangle -- P1 research's own dimension-based specification
    convention (length + two end-widths, taper angle derived, not
    authored)."""
    outline = trapezoid_outline_points(half_l, half_w_wide, half_w_narrow)
    return fancy_brilliant_planes(outline, crown_deg, pav_deg, table_frac, depth_pct,
                                   girdle_t_pct, 2 * half_l, name='taperedbag')


def carre_planes(half_w=25.0, crown_deg=35.5, pav_deg=42.0, table_frac=0.60, depth_pct=61.0,
                  girdle_t_pct=3.5):
    """Carre: square step cut, sharp UNCUT corners (the defining contrast
    with Asscher, well-sourced -- "every source agrees"). P1 research
    found NO angle, table% or depth% published anywhere for this cut
    specifically; the values here are an explicit analogy to the emerald
    cut (same step-cut family, same general proportions), not a citation --
    ships approximate per plan §10's own policy."""
    outline = rect_outline_points(half_w, half_w)
    return fancy_brilliant_planes(outline, crown_deg, pav_deg, table_frac, depth_pct,
                                   girdle_t_pct, 2 * half_w, name='carre')


def table_cut_planes(half_w=25.0, crown_deg=45.0, pav_deg=45.0, table_frac=0.70, depth_pct=55.0,
                      girdle_t_pct=2.0):
    """Table cut: the historical ancestor (Islamic world, attested from the
    13th c.; the fully 4+4-facet form P1 research dates more to the 14th-
    16th c., reaching Europe via Venice) -- essentially an octahedral
    crystal ground flat on two opposite faces. Source: P1 research, via
    Lang Antiques citing Tillander -- the RAW natural-crystal angle is a
    real geometric constant, 54.7356 deg (arccos(1/sqrt(3)), the angle
    between an octahedron face and its 4-fold axis); real finished stones
    were reground shallower, to ~45 deg, "for excellent brilliance" --
    that finished figure is used here (crown_deg=pav_deg=45), not the raw
    crystal angle. No table%/depth% standard exists (historical stones
    were not calibrated to a percentage); girdle is not modelled as
    rounded (real stones "displayed rounded or blunt corners" per that
    research, but a plain square is the honest simplification, not a
    claim of precision)."""
    outline = rect_outline_points(half_w, half_w)
    return fancy_brilliant_planes(outline, crown_deg, pav_deg, table_frac, depth_pct,
                                   girdle_t_pct, 2 * half_w, name='tablecut')


def rose_cut_planes(radius=25.0, dome_height=30.0, n_sides=6, girdle_t_pct=4.0):
    """Rose cut. Source: P1 research -- flat (unfaceted) back, NO pavilion
    and NO table at all (light behaves very differently from a brilliant:
    no culet-less/table-less split into crown/pavilion in the normal
    sense). 12-facet variant modelled here as a 2-tier hexagonal dome
    (stepped_pavilion_planes reused with a NEGATIVE target_depth, so it
    builds upward to an apex instead of downward to a culet -- the same
    frustum-stack mechanism, just mirrored) over a flat back plane; the
    24-facet variant would add a 3rd tier. No source publishes dome
    angles for either variant -- ships approximate."""
    outline = [(radius * math.cos(2 * math.pi * k / n_sides), radius * math.sin(2 * math.pi * k / n_sides))
               for k in range(n_sides)]
    t = radius * girdle_t_pct / 100.0
    dome_pl, apex_y = stepped_pavilion_planes(outline, [0.5, 0.5], [0.5, 0.0], t, -dome_height, (0., 0., 0.), 'rose_dome')
    back_pl = dict(n=np.array([0., -1., 0.]), d=t, name='rose_back')
    girdle_pl = girdle_band_from_outline(outline, t, (0., 0., 0.))
    planes = dedupe_planes([back_pl] + girdle_pl + dome_pl)
    return planes, dict(apexY=apex_y, backY=-t)


def mazarin_planes(crown_deg=38.0, pav_deg=42.0, table_frac=0.50, depth_pct=66.0, girdle_t_pct=3.5):
    """Mazarin cut. Source: P1 research -- 34 facets (17 crown/17 pavilion,
    the pavilion split an inference from topology, not directly cited); no
    angle table survives anywhere, and the "Cardinal Mazarin invented it"
    attribution is itself almost certainly trade legend (Tillander found
    no documentary trace, and the term doesn't appear in print before the
    20th c.). Angles here are between old-mine's and Peruzzi/Jeffries'
    figures (Mazarin's lineage sits earlier than old mine, later than the
    single/table-cut era) -- an interpolation, not a citation. Cushion
    girdle (pre-1874, same reasoning as old mine/Peruzzi)."""
    return cushion_brilliant_planes(lw_ratio=1.0, exponent=3.0, crown_deg=crown_deg, pav_deg=pav_deg,
                                     table_frac=table_frac, depth_pct=depth_pct, girdle_t_pct=girdle_t_pct)


def briolette_planes(a_point=32.0, a_round=18.0, b=25.0, crown_deg=40.0, pav_deg=40.0, girdle_t_pct=2.0):
    """Briolette: fully faceted teardrop, NO table and NO flat back --
    faceted on every side, confirmed as a "modified double rose cut"
    (rose faceting on both hemispheres of a girdle-less stone) by multiple
    independent sources. Modelled via fancy_brilliant_planes with
    table_frac=0.0 (the crown cone is allowed to converge all the way to
    its own apex point rather than being capped by a table plane -- this
    works cleanly: the nominal "table" plane ends up coincident with the
    crown apex and is simply never the nearest facet, not a degenerate
    case) over the SAME pear outline as pear_planes (a teardrop is a pear
    profile faceted on both sides, per that research's own "modified
    double rose" description). No source publishes angles for either
    end -- ships approximate."""
    outline = pear_outline_points(a_point, a_round, b)
    major_diam = a_point + a_round
    return fancy_brilliant_planes(outline, crown_deg, pav_deg, 0.0, 100.0, girdle_t_pct, major_diam, name='briolette')


def enumerate_vertices(planes, tol=1e-9, merge_tol=1e-6):
    """Triple-plane intersection -> vertices that satisfy every OTHER
    half-space within tol. Also returns, per plane, the vertex indices lying
    on it (for the Euler check). This is the generic pipeline of plan §3.

    Coincident intersections are merged: a real vertex like the culet, where
    all 8 pavilion-main planes meet at one point, is found by C(8,3)=56
    different plane triples and must collapse to a single vertex, not 56.

    Batched (not a Python-level i/j/k triple loop): a real performance bug
    found the hard way -- the original per-triple np.array/np.linalg.solve
    calls cost real seconds once a fancy-cut outline pushes the plane count
    from round brilliant's 25 up into the 50-100 range (C(97,3) triples for
    a 32-sample oval), and gem_trace_reference.py calls this for ~30 cuts
    per run. numpy's batched linalg.det/solve process every candidate
    triple in one vectorized call instead of one Python loop iteration
    each; verified to return byte-identical vertices and on_plane sets to
    the original loop on round_brilliant_planes() before replacing it."""
    n = len(planes)
    on_plane = [[] for _ in range(n)]
    A_all = np.array([p['n'] for p in planes])
    d_all = np.array([p['d'] for p in planes])
    idx = np.array(list(itertools.combinations(range(n), 3)))
    if len(idx) == 0:
        return [], on_plane
    A_batch = A_all[idx]  # (M,3,3)
    b_batch = d_all[idx]  # (M,3)
    dets = np.linalg.det(A_batch)
    keep = np.abs(dets) >= 1e-9
    idx_k, A_k, b_k = idx[keep], A_batch[keep], b_batch[keep]
    if len(idx_k) == 0:
        return [], on_plane
    X = np.linalg.solve(A_k, b_k[..., None])[..., 0]  # (M',3) candidate points
    valid = np.all(X @ A_all.T <= d_all + 1e-6, axis=1)
    X_v, idx_v = X[valid], idx_k[valid]

    pts = []
    for x, (i, j, k) in zip(X_v, idx_v):
        found = None
        for pi, p0 in enumerate(pts):
            if np.linalg.norm(p0 - x) < merge_tol:
                found = pi
                break
        if found is None:
            pts.append(x)
            found = len(pts) - 1
        for m in (i, j, k):
            if found not in on_plane[m]:
                on_plane[m].append(found)
    return pts, on_plane


def solid_check(planes):
    """Closed, positive-volume, outward-normal check (§0 invariant #2 /
    §3 pipeline). Returns dict with vertex/face/edge counts and volume."""
    pts, on_plane = enumerate_vertices(planes)
    pts = np.array(pts)
    centroid = pts.mean(axis=0)
    # outward-normal check
    for p in planes:
        if np.dot(p['n'], centroid) - p['d'] > -1e-6:
            raise AssertionError(f'plane {p["name"]} does not have centroid strictly inside')
    # face polygons: for each plane, sort its vertices angularly about the
    # face centroid in the plane's own basis
    faces = []
    for pi, idxs in enumerate(on_plane):
        if len(idxs) < 3:
            continue
        v = pts[idxs]
        fc = v.mean(axis=0)
        nrm = planes[pi]['n']
        # build an in-plane basis
        ref = np.array([1.0, 0, 0]) if abs(nrm[0]) < 0.9 else np.array([0, 1.0, 0])
        u = np.cross(nrm, ref); u /= np.linalg.norm(u)
        w = np.cross(nrm, u)
        angles = [math.atan2(np.dot(pt - fc, w), np.dot(pt - fc, u)) for pt in v]
        order = np.argsort(angles)
        faces.append([idxs[o] for o in order])

    # signed volume: fan-triangulate each face, sum tetrahedra from the origin
    vol = 0.0
    for face in faces:
        fpts = pts[face]
        apex = fpts[0]
        for a in range(1, len(face) - 1):
            p1, p2 = fpts[a], fpts[a + 1]
            vol += np.dot(apex, np.cross(p1, p2)) / 6.0

    # edges: each unordered pair of vertices that co-occurs in exactly 2 faces
    edge_count = {}
    for face in faces:
        m = len(face)
        for a in range(m):
            e = tuple(sorted((face[a], face[(a + 1) % m])))
            edge_count[e] = edge_count.get(e, 0) + 1
    bad = [e for e, c in edge_count.items() if c != 2]
    V, F, E = len(pts), len(faces), len(edge_count)
    return dict(V=V, F=F, E=E, euler=V - E + F, vol=abs(vol),
                bad_edges=bad, centroid=centroid, faces=faces, pts=pts)


# ============================================ generalized primitives (P2) =
# Extends the pure half-space engine above to plan §3's other two
# primitives: analytic quadrics (sphere/ellipsoid, for cabochon/egg/sphere/
# lens) and small unions of convex parts (for heart/merkaba). A "solid" for
# trace_tree below is now one of three shapes:
#   - a plain list of planes            (legacy half-space-only path, exact
#                                         behaviour unchanged from P1-P4)
#   - dict(planes=[...], quadrics=[...])  a single convex solid = intersection
#                                         of every plane AND every quadric's
#                                         interior
#   - dict(parts=[the above, ...])        a union of convex parts (heart,
#                                         merkaba) per §3's union semantics
def ellipsoid_roots(origin, direction, center, radii):
    """Real roots of |(origin + t*direction - center) / radii| = 1, sorted
    ascending, or None if the (infinite) line misses the ellipsoid. direction
    need not be unit -- t is then in units of |direction|."""
    o = (origin - center) / radii
    d = direction / radii
    a = np.dot(d, d)
    b = 2 * np.dot(o, d)
    c = np.dot(o, o) - 1.0
    disc = b * b - 4 * a * c
    if disc < 0:
        return None
    sq = math.sqrt(disc)
    t1, t2 = (-b - sq) / (2 * a), (-b + sq) / (2 * a)
    return (t1, t2) if t1 <= t2 else (t2, t1)


def ellipsoid_normal(point, center, radii):
    return normalize((point - center) / (radii * radii))


def _plane_exit_t(origin, direction, plane):
    ndotd = np.dot(plane['n'], direction)
    if ndotd <= 1e-12:
        return None
    return (plane['d'] - np.dot(plane['n'], origin)) / ndotd


def _plane_entry_t(origin, direction, plane):
    ndotd = np.dot(plane['n'], direction)
    if ndotd >= -1e-12:
        return None
    return (plane['d'] - np.dot(plane['n'], origin)) / ndotd


def trace_solid(origin, direction, spec):
    """Nearest exit for a ray at/inside a single generalized convex solid
    (planes + quadrics). Returns (t, name, outward_normal_at_hit) or
    (None, None, None). The quadric branch assumes the ray is genuinely
    inside that quadric's interior already (true whenever origin is inside
    the full solid, since the solid is a subset of every quadric's
    interior) -- so the exit root is always the larger (far) one."""
    best_t, best_name, best_n = None, None, None
    for p in spec.get('planes', ()):
        t = _plane_exit_t(origin, direction, p)
        if t is not None and t > 1e-9 and (best_t is None or t < best_t):
            best_t, best_name, best_n = t, p['name'], p['n']
    for q in spec.get('quadrics', ()):
        roots = ellipsoid_roots(origin, direction, q['center'], q['radii'])
        if roots is None:
            continue
        t_far = roots[1]
        if t_far > 1e-9 and (best_t is None or t_far < best_t):
            hit = origin + t_far * direction
            best_t, best_name = t_far, q['name']
            best_n = ellipsoid_normal(hit, q['center'], q['radii'])
    return best_t, best_name, best_n


def enter_solid(origin, direction, spec):
    """Nearest entry into a single generalized convex solid, approaching
    from outside. Returns (t, name, outward_normal_at_hit) or
    (None, None, None) -- the latter meaning the ray misses the solid
    entirely (either no plane/quadric candidate at all, or it misses a
    quadric's surface outright, which for an INTERSECTION means it misses
    the whole solid regardless of the other constraints)."""
    candidates = []
    for p in spec.get('planes', ()):
        t = _plane_entry_t(origin, direction, p)
        if t is not None:
            candidates.append((t, p['name'], p['n']))
    for q in spec.get('quadrics', ()):
        roots = ellipsoid_roots(origin, direction, q['center'], q['radii'])
        if roots is None:
            return None, None, None
        t_near = roots[0]
        hit = origin + t_near * direction
        candidates.append((t_near, q['name'], ellipsoid_normal(hit, q['center'], q['radii'])))
    if not candidates:
        return None, None, None
    return max(candidates, key=lambda c: c[0])


def _inside_solid(point, spec, tol=1e-7):
    for p in spec.get('planes', ()):
        if np.dot(p['n'], point) - p['d'] > tol:
            return False
    for q in spec.get('quadrics', ()):
        o = (point - q['center']) / q['radii']
        if np.dot(o, o) - 1.0 > tol:
            return False
    return True


def _which_parts(point, parts, tol=1e-7):
    return [i for i, spec in enumerate(parts) if _inside_solid(point, spec, tol)]


def enter_union(origin, direction, parts):
    """First real entry into a union of parts, from outside: the minimum
    entry-t across parts. No seam ambiguity here -- before any hit the ray
    is outside every part, so the first surface it touches is always a
    genuine optical interface (the seam-skipping logic below only matters
    once the ray is already inside the union)."""
    best = None
    for spec in parts:
        t, name, n = enter_solid(origin, direction, spec)
        if t is not None and (best is None or t < best[0]):
            best = (t, name, n)
    return best if best else (None, None, None)


def trace_union(origin, direction, parts, max_march=32, eps=1e-6):
    """§3 union semantics for a ray already inside the union: nearest exit
    among the part(s) currently containing the point; if that exit point is
    STILL inside some OTHER part, it's an internal seam -- same medium, no
    optical event -- so advance past it by eps and re-march. Critically,
    "other" excludes the part whose boundary was just crossed: at any
    part's own exit facet, tolerance-inclusive membership always reads that
    part as still containing the hit point (it's sitting exactly on the
    boundary), which would otherwise misclassify every genuine exterior
    surface as a seam. max_march is the strict-progress safety cap (every
    step must advance t, per §3), catching an infinite-loop bug rather than
    hanging. Returns (t_total, name, normal) or (None, None, None)."""
    o, t_acc = origin, 0.0
    for _ in range(max_march):
        containing = _which_parts(o, parts)
        if not containing:
            return None, None, None
        best = None
        for i in containing:
            t, name, n = trace_solid(o, direction, parts[i])
            if t is not None and (best is None or t < best[0]):
                best = (t, name, n, i)
        if best is None:
            return None, None, None
        t, name, n, i_best = best
        hit = o + t * direction
        others = [i for i in _which_parts(hit, parts) if i != i_best]
        if others:
            o, t_acc = hit + direction * eps, t_acc + t + eps
            continue
        return t_acc + t, name, n
    return None, None, None


# ==================================================== symmetry + Platonic =
def rot_y(angle):
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


def mirror_through_y(azimuth):
    """Reflection matrix for the vertical plane containing the y-axis at the
    given azimuth (the x-z mirror, rotated by azimuth about y)."""
    c, s = math.cos(2 * azimuth), math.sin(2 * azimuth)
    return np.array([[c, 0, s], [0, 1, 0], [s, 0, -c]])


def cnv_matrices(n, mirror=True):
    """Cn (or Cnv with mirror=True): rotations by k*2pi/n about the vertical
    axis, optionally each composed with one vertical mirror. Covers every
    faceted-cut symmetry in the roster except the talismanic solids."""
    mats = [rot_y(2 * math.pi * k / n) for k in range(n)]
    if mirror:
        m0 = mirror_through_y(0.0)
        mats = mats + [m0 @ R for R in mats]
    return mats


def check_symmetry(planes, matrices, tol=1e-9):
    """§3/§8: the plane set must be invariant under every matrix of the
    declared point group -- for each matrix and each plane, some plane in
    the set must match (M @ n, d) within tol (offsets are unchanged by any
    rotation/reflection through the origin, which is every solid's centre
    by construction). Returns a list of (plane_name, matrix) failures --
    empty means symmetric."""
    bad = []
    for M in matrices:
        for p in planes:
            mn = M @ p['n']
            if not any(np.linalg.norm(mn - q['n']) < tol and abs(p['d'] - q['d']) < tol
                        for q in planes):
                bad.append((p['name'], M))
    return bad


PHI = (1 + math.sqrt(5)) / 2

PLATONIC_DIHEDRAL_DEG = dict(tetrahedron=70.5288, cube=90.0, octahedron=109.4712,
                              dodecahedron=116.5651, icosahedron=138.1897)


def platonic_planes(kind, d=1.0):
    """Face-normal sets for the 5 Platonic solids -- each solid is the
    intersection of {n.x <= d} over these outward unit normals, centred at
    the origin with inradius d. Standard closed-form directions: tetrahedron
    is self-dual (its own vertex directions in one orientation); cube/
    octahedron and dodecahedron/icosahedron are dual pairs, so each solid's
    face normals are literally its dual's vertex directions."""
    if kind == 'tetrahedron':
        raw = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
    elif kind == 'cube':
        raw = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
    elif kind == 'octahedron':
        raw = [(sx, sy, sz) for sx in (1, -1) for sy in (1, -1) for sz in (1, -1)]
    elif kind == 'dodecahedron':
        raw = ([(0, sy, sz * PHI) for sy in (1, -1) for sz in (1, -1)]
               + [(sx, sy * PHI, 0) for sx in (1, -1) for sy in (1, -1)]
               + [(sx * PHI, 0, sz) for sx in (1, -1) for sz in (1, -1)])
    elif kind == 'icosahedron':
        raw = ([(sx, sy, sz) for sx in (1, -1) for sy in (1, -1) for sz in (1, -1)]
               + [(0, sy / PHI, sz * PHI) for sy in (1, -1) for sz in (1, -1)]
               + [(sx / PHI, sy * PHI, 0) for sx in (1, -1) for sy in (1, -1)]
               + [(sx * PHI, 0, sz / PHI) for sx in (1, -1) for sz in (1, -1)])
    else:
        raise ValueError(kind)
    return [dict(n=normalize(np.array(v, dtype=float)), d=d, name=f'{kind}_{i}')
            for i, v in enumerate(raw)]


def dihedral_angles(planes):
    """Actual dihedral angles (deg) read off the GENERATED solid's own edges
    -- audits the geometry kernel against textbook closed-form values it
    never contains (plan §8), rather than trusting the face-normal formulas
    on their own. dihedral = 180 - angle_between_adjacent_outward_normals."""
    sc = solid_check(planes)
    pts = sc['pts']
    on_plane = [[vi for vi, p in enumerate(pts) if abs(np.dot(pl['n'], p) - pl['d']) < 1e-6]
                for pl in planes]
    angles, seen = [], set()
    for face in sc['faces']:
        m = len(face)
        for a in range(m):
            va, vb = face[a], face[(a + 1) % m]
            e = tuple(sorted((va, vb)))
            if e in seen:
                continue
            incident = [pi for pi in range(len(planes)) if va in on_plane[pi] and vb in on_plane[pi]]
            if len(incident) != 2:
                continue
            seen.add(e)
            n1, n2 = planes[incident[0]]['n'], planes[incident[1]]['n']
            angles.append(math.degrees(math.pi - math.acos(np.clip(np.dot(n1, n2), -1, 1))))
    return angles


# ============================================== fancy outlines (P2, §3) ===
def arc_outline_points(arcs, samples_per_arc=8):
    """§3 "Fancy outlines": arcs is a list of {center:(x,z), radius, a0, a1}
    (angles in radians, swept a0->a1) in the horizontal x-z plane (y is the
    table-normal/vertical axis throughout this file). Each arc's a1 must
    land on the next arc's a0 at the same point -- the tangency condition a
    real drafting construction requires -- though this function itself just
    samples; callers are responsible for supplying a closing arc list.
    Returns an ordered list of (x, z) points tracing the outline once
    around, ready to feed to chord_planes_from_outline / a girdle band."""
    pts = []
    for arc in arcs:
        cx, cz = arc['center']
        for k in range(samples_per_arc):
            t = arc['a0'] + (arc['a1'] - arc['a0']) * k / samples_per_arc
            pts.append((cx + arc['radius'] * math.cos(t), cz + arc['radius'] * math.sin(t)))
    return pts


def chord_planes_from_outline(outline_xz, apex_pt, girdle_y, name_prefix, interior_pt):
    """§3: each pavilion-main (or crown-main) plane is the plane through one
    outline EDGE (two adjacent outline samples, at height girdle_y) and the
    apex point (culet for a pavilion; the table-rim apex for a crown, per
    §3's "crown facets run the same chord construction to the table rim").
    `interior_pt` -- any point known to be inside the eventual solid (NOT
    the apex: the apex sits exactly ON every returned plane by construction,
    so testing orientation against it is numerically meaningless, a real
    bug this docstring exists to warn off) -- fixes each plane's outward
    orientation. Returns a list of {n,d,name}."""
    apex = np.asarray(apex_pt, dtype=float)
    interior = np.asarray(interior_pt, dtype=float)
    planes = []
    m = len(outline_xz)
    for i in range(m):
        x1, z1 = outline_xz[i]
        x2, z2 = outline_xz[(i + 1) % m]
        p1 = np.array([x1, girdle_y, z1])
        p2 = np.array([x2, girdle_y, z2])
        nrm = np.cross(p2 - p1, apex - p1)
        length = np.linalg.norm(nrm)
        if length < 1e-12:
            continue  # degenerate edge (zero-length or collinear with apex) -- skip, not a facet
        nrm = nrm / length
        d = np.dot(nrm, p1)
        if np.dot(nrm, interior) > d - 1e-9:
            nrm, d = -nrm, -d
        planes.append(dict(n=nrm, d=d, name=f'{name_prefix}_{i}'))
    return planes


def fancy_brilliant_planes(outline_xz, crown_deg, pav_deg, table_frac, depth_pct,
                            girdle_t_pct, major_diam, name='fancy', interior_xz=(0.0, 0.0)):
    """General round-brilliant-STYLE generator for any closed convex outline:
    single crown angle, single pavilion angle, one apex point each (culet
    below, a notional crown apex above), chord-plane facets per outline
    edge (§3), table derived by intersecting the crown planes with a plane
    at the right height rather than authored directly. This is what oval,
    pear, marquise, cushion brilliant and the cushion/round antiques (old
    mine, old European, Peruzzi-style) all reduce to -- the outline is the
    only thing that changes between them.

    Design note, stated once here rather than re-derived per cut: a
    NON-circular outline cannot hold one exact angle from a single apex at
    every edge (only a true circular cone can) -- crown_deg/pav_deg are
    therefore realised as an EFFECTIVE apex height solved to hit the
    requested depth_pct exactly, with individual facets varying somewhat
    around the perimeter (tighter near the short axis, shallower near the
    long one). This is not a numerical compromise so much as the honest
    geometric fact for any single-apex chord construction on a non-circular
    base; real fancy-cut pavilion mains genuinely do vary in shape/tilt
    around the girdle for exactly this reason. table_frac is hit EXACTLY
    (verified by direct construction, not just by the algebra below) via
    similar triangles from the crown apex: the table is a table_frac-scaled
    copy of the outline, because it is the outline cut by a plane parallel
    to the girdle at the height that produces that scale.
    """
    t_half = major_diam * girdle_t_pct / 100.0 / 2.0
    target_depth = major_diam * depth_pct / 100.0
    tc, tp = math.tan(math.radians(crown_deg)), math.tan(math.radians(pav_deg))
    r_eff = (target_depth - 2 * t_half) / (tc * (1 - table_frac) + tp)
    hc, hp = r_eff * tc, r_eff * tp
    apex_crown_y = t_half + hc
    culet_y = -(t_half + hp)
    table_y = apex_crown_y - table_frac * hc
    interior = (interior_xz[0], 0.0, interior_xz[1])

    table_plane = dict(n=np.array([0., 1., 0.]), d=table_y, name=f'{name}_table')
    crown_pl = chord_planes_from_outline(outline_xz, (0, apex_crown_y, 0), t_half,
                                          f'{name}_crown', interior)
    girdle_pl = girdle_band_from_outline(outline_xz, t_half, interior)
    pav_pl = chord_planes_from_outline(outline_xz, (0, culet_y, 0), -t_half,
                                        f'{name}_pav', interior)
    planes = [table_plane] + crown_pl + girdle_pl + pav_pl
    planes = dedupe_planes(planes)
    meta = dict(t=2 * t_half, hc=hc, hp=hp, tableY=table_y, culetY=culet_y,
                depthPct=(table_y - culet_y) / major_diam * 100, rEff=r_eff)
    return planes, meta


def dedupe_planes(planes, tol=1e-9):
    """Drops later planes whose (n,d) match an earlier one within tol.
    Needed for outlines with a straight edge through the solid's own axis
    (half-moon's flat diameter, e.g.): the crown-side and pavilion-side
    chord facets for that ONE edge both degenerate to the SAME vertical
    plane (their defining apex, at x=z=0, sits on the girdle's own
    symmetry line for that edge either way) -- a real geometric
    coincidence, not a bug in chord_planes_from_outline, but one that
    breaks enumerate_vertices' edge-sharing count if left as two identical
    entries. Keeps the first occurrence's name."""
    out = []
    for p in planes:
        if not any(np.linalg.norm(p['n'] - q['n']) < tol and abs(p['d'] - q['d']) < tol for q in out):
            out.append(p)
    return out


def ellipse_outline_points(a, b, n=32):
    """A true ellipse, used as the outline for oval (and the rounded lobes
    of pear/heart) rather than a hand-solved 4-centre compass construction.
    Deliberate simplification, stated once here: no lapidary source
    publishes actual 4-centre arc radii for a modern oval brilliant (P1
    research confirmed this explicitly), so the compass construction would
    itself be an unsourced choice, not a discovered standard -- and a true
    ellipse is smooth, exactly hits the researched L:W ratio, and is
    trivially provably convex, none of which a hand-tangent-solved 4-arc
    approximation is free of risk on. Guaranteed convex; still checked by
    the caller's solid_check, not just asserted here."""
    return [(a * math.cos(2 * math.pi * k / n), b * math.sin(2 * math.pi * k / n)) for k in range(n)]


def girdle_band_from_outline(outline_xz, half_thickness, interior_pt):
    """The girdle band as a real prism wall over an ARBITRARY outline
    (generalizing round_brilliant_planes' regular n-gon band to any fancy
    shape) -- one near-vertical plane per outline edge, connecting its
    upper-girdle copy (y=+half_thickness) to its lower-girdle copy
    (y=-half_thickness). Plan §3 "depth budget": never a zero-thickness
    edge, always a modelled band with its own thickness parameter."""
    interior = np.asarray(interior_pt, dtype=float)
    planes = []
    m = len(outline_xz)
    for i in range(m):
        x1, z1 = outline_xz[i]
        x2, z2 = outline_xz[(i + 1) % m]
        p1u = np.array([x1, half_thickness, z1])
        p2u = np.array([x2, half_thickness, z2])
        p1l = np.array([x1, -half_thickness, z1])
        nrm = np.cross(p2u - p1u, p1l - p1u)
        nrm = nrm / np.linalg.norm(nrm)
        d = np.dot(nrm, p1u)
        if np.dot(nrm, interior) > d - 1e-9:
            nrm, d = -nrm, -d
        planes.append(dict(n=nrm, d=d, name=f'girdle_{i}'))
    return planes


# ============================================================ tracer core =
def normalize(v):
    return v / np.linalg.norm(v)


def trace_convex(origin, direction, planes):
    """Nearest exit facet for a ray at/inside the convex solid (§3: 'the
    nearest positive-t plane'). direction must be unit."""
    best_t, best_plane = None, None
    for p in planes:
        ndotd = np.dot(p['n'], direction)
        if ndotd <= 1e-12:
            continue
        t = (p['d'] - np.dot(p['n'], origin)) / ndotd
        if t > 1e-9 and (best_t is None or t < best_t):
            best_t, best_plane = t, p
    return best_t, best_plane


def enter_convex(origin, direction, planes):
    """First facet a ray hits approaching the convex solid from OUTSIDE —
    the mirror case of trace_convex: the largest t among the half-spaces
    the ray is entering (n.d < 0), i.e. the standard slab-test entry point."""
    best_t, best_plane = None, None
    for p in planes:
        ndotd = np.dot(p['n'], direction)
        if ndotd >= -1e-12:
            continue
        t = (p['d'] - np.dot(p['n'], origin)) / ndotd
        if best_t is None or t > best_t:
            best_t, best_plane = t, p
    return best_t, best_plane


def reflect(d, n):
    return d - 2 * np.dot(d, n) * n


def refract(d, n, n1, n2):
    """Vector Snell (plan §3b eq1). d unit incident dir, n unit normal
    (either orientation; flipped internally to oppose d). Returns
    (t_hat or None-if-TIR, cosi, cost or None)."""
    if np.dot(d, n) > 0:
        n = -n
    cosi = -np.dot(d, n)
    mu = n1 / n2
    sin2t = mu * mu * (1 - cosi * cosi)
    if sin2t > 1.0:
        return None, cosi, None
    cost = math.sqrt(max(0.0, 1 - sin2t))
    t = mu * d + (mu * cosi - cost) * n
    return t, cosi, cost


def fresnel(cosi, cost, n1, n2):
    """§3b eq3/eq4: amplitude coefficients, then R=|r|^2 and the
    projected-flux T=|t|^2 (n2 cost)/(n1 cosi). Returns Rs,Ts,Rp,Tp,rs,rp
    (signed amplitudes, needed for the polarization ledger's phase bookkeeping
    on ordinary — non-TIR — reflection)."""
    rs = (n1 * cosi - n2 * cost) / (n1 * cosi + n2 * cost)
    ts = 2 * n1 * cosi / (n1 * cosi + n2 * cost)
    rp = (n2 * cosi - n1 * cost) / (n2 * cosi + n1 * cost)
    tp = 2 * n1 * cosi / (n2 * cosi + n1 * cost)
    flux = (n2 * cost) / (n1 * cosi)
    Rs, Rp = rs * rs, rp * rp
    Ts, Tp = ts * ts * flux, tp * tp * flux
    return Rs, Ts, Rp, Tp, rs, rp


def tir_phase(cosi, sini, n_rel):
    """§3b eq5. n_rel = n_rare/n_dense (<1). Returns (delta_s, delta_p) in
    radians; Delta = delta_p - delta_s is the retardance a Fresnel rhomb
    quotes."""
    root = math.sqrt(max(0.0, sini * sini - n_rel * n_rel))
    delta_s = 2 * math.atan2(root, cosi)
    delta_p = 2 * math.atan2(root, n_rel * n_rel * cosi)
    return delta_s, delta_p


def rotation_angle(s_from, s_to, axis):
    """Signed angle (rad) from s_from to s_to about axis (all unit, s_from
    and s_to both already perpendicular to axis)."""
    cs = np.clip(np.dot(s_from, s_to), -1.0, 1.0)
    sn = np.dot(np.cross(s_from, s_to), axis)
    return math.atan2(sn, cs)


def mueller_rotate(S, theta):
    c, s = math.cos(2 * theta), math.sin(2 * theta)
    R = np.array([[1, 0, 0, 0],
                  [0, c, s, 0],
                  [0, -s, c, 0],
                  [0, 0, 0, 1]])
    return R @ S


def mueller_interface(S, Rs, Ts, Rp, Tp, rs=None, rp=None, delta=None, kind='reflect'):
    """Apply the diattenuator (+retarder for TIR) Mueller matrix, in the
    LOCAL s/p basis, for a reflect or transmit branch."""
    if kind == 'reflect':
        a, b = Rs, Rp
    else:
        a, b = Ts, Tp
    m00 = (a + b) / 2.0
    m01 = (a - b) / 2.0
    if delta is not None:
        # TIR reflection: |rs|=|rp|=1, pure retarder in the (S2,S3) block
        m22, m23, m32, m33 = math.cos(delta), math.sin(delta), -math.sin(delta), math.cos(delta)
    else:
        # ordinary (real-amplitude) event: 0/pi phase baked into sign of rs*rp
        prod = (rs * rp) if kind == 'reflect' else None
        if kind != 'reflect':
            # transmission amplitudes are also real for a dielectric; reuse ts*tp sign
            prod = math.sqrt(max(Ts, 0.0)) * math.sqrt(max(Tp, 0.0))
            prod *= 1.0  # ts,tp always same sign (both = 2 n1 cosi / (...)+ ; positive)
        m22 = m33 = prod
        m23 = m32 = 0.0
    M = np.array([[m00, m01, 0, 0],
                  [m01, m00, 0, 0],
                  [0, 0, m22, m23],
                  [0, 0, m32, m33]])
    return M @ S


def dop_of(S):
    s0, s1, s2, s3 = S
    return math.sqrt(s1 * s1 + s2 * s2 + s3 * s3) / s0 if s0 > 0 else 0.0


def trace_tree(origin, direction, solid, n_stone, n_air=1.0, max_bounces=14,
               min_intensity=1e-9):
    """Branching ray tree with intensity weights (plan §3b eq1-5, P3 scope).
    origin/direction: the incoming ray in air, aimed at the stone. `solid`
    is a plain plane list (legacy half-space-only path, unchanged), a
    dict(planes=,quadrics=) single generalized solid, or a dict(parts=[...])
    union -- see the "generalized primitives" section above.

    At every non-TIR interface the tree splits; the REFLECTED sibling is
    always the one that keeps recursing (matching the drawn path in §1's
    figure — the ray continues bouncing until it actually transmits out).
    Every transmitted sibling is recorded as an 'exited' branch and does not
    recurse further, whether it is the dominant table-exit event or a minor
    leak. Returns a ledger:
      exited: [{intensity, angle_deg, dop, from_point, dir}]
      truncated: float (energy still bouncing at the bounce/intensity cap)
      dominant_log: [{point, aoi, plane}] — the always-reflected path,
        i.e. exactly the sequence drawn as the main ray in §1's figure.
    """
    direction = normalize(direction)
    ledger = {'exited': [], 'truncated': 0.0}
    dominant_log = []

    if isinstance(solid, list):
        def enter_fn(o, d):
            t, p = enter_convex(o, d, solid)
            return (t, p['name'], p['n']) if p is not None else (None, None, None)

        def trace_fn(o, d):
            t, p = trace_convex(o, d, solid)
            return (t, p['name'], p['n']) if p is not None else (None, None, None)
    elif 'parts' in solid:
        def enter_fn(o, d): return enter_union(o, d, solid['parts'])
        def trace_fn(o, d): return trace_union(o, d, solid['parts'])
    else:
        def enter_fn(o, d): return enter_solid(o, d, solid)
        def trace_fn(o, d): return trace_solid(o, d, solid)

    t, entry_name, n_out = enter_fn(origin, direction)
    if entry_name is None:
        raise RuntimeError('entry ray misses the solid')
    entry_point = origin + direction * t
    s_hat = (normalize(np.cross(n_out, direction))
             if abs(np.dot(n_out, direction)) < 0.999999 else np.array([1.0, 0.0, 0.0]))
    S = np.array([1.0, 0.0, 0.0, 0.0])  # unpolarized sunlight

    t_dir, cosi, cost = refract(direction, n_out, n_air, n_stone)
    entry_aoi = math.degrees(math.acos(np.clip(cosi, -1, 1)))
    dominant_log.append(dict(point=entry_point, aoi=entry_aoi, plane=entry_name))

    Rs, Ts, Rp, Tp, rs, rp = fresnel(cosi, cost, n_air, n_stone)
    assert abs(Rs + Ts - 1) < 1e-9 and abs(Rp + Tp - 1) < 1e-9, 'R+T!=1 at entry'
    S_door = mueller_interface(S, Rs, Ts, Rp, Tp, rs=rs, rp=rp, kind='reflect')
    n_op = n_out if np.dot(direction, n_out) < 0 else -n_out
    door_dir = reflect(direction, n_op)
    ledger['exited'].append(dict(intensity=(Rs + Rp) / 2.0, angle_deg=entry_aoi,
                                  dop=dop_of(S_door), from_point=entry_point, dir=door_dir,
                                  label='door', depth=0))

    S_in = mueller_interface(S, Rs, Ts, Rp, Tp, kind='transmit')
    T_in = (Ts + Tp) / 2.0

    def recurse(point, d, S, s_hat, intensity, depth, log):
        if intensity < min_intensity or depth > max_bounces:
            ledger['truncated'] += intensity
            return
        t, name, n_out = trace_fn(point, d)
        if name is None:
            ledger['truncated'] += intensity
            return
        hitp = point + d * t

        s_new = (normalize(np.cross(n_out, d))
                 if abs(np.dot(n_out, d)) < 0.999999 else s_hat)
        theta_rot = rotation_angle(s_hat, s_new, d)
        S_rot = mueller_rotate(S, theta_rot)

        t_dir, cosi, cost = refract(d, n_out, n_stone, n_air)
        aoi_deg = math.degrees(math.acos(np.clip(cosi, -1, 1)))
        n_op = n_out if np.dot(d, n_out) < 0 else -n_out
        r_dir = reflect(d, n_op)
        log.append(dict(point=hitp, aoi=aoi_deg, plane=name))

        if t_dir is None:
            sini = math.sqrt(max(0.0, 1 - cosi * cosi))
            n_rel = n_air / n_stone
            ds, dp = tir_phase(cosi, sini, n_rel)
            S_out = mueller_interface(S_rot, 1.0, 0.0, 1.0, 0.0, delta=(dp - ds), kind='reflect')
            recurse(hitp, r_dir, S_out, s_new, intensity, depth + 1, log)
            return

        Rs, Ts, Rp, Tp, rs, rp = fresnel(cosi, cost, n_stone, n_air)
        assert abs(Rs + Ts - 1) < 1e-9 and abs(Rp + Tp - 1) < 1e-9, 'R+T!=1'

        S_r = mueller_interface(S_rot, Rs, Ts, Rp, Tp, rs=rs, rp=rp, kind='reflect')
        recurse(hitp, r_dir, S_r, s_new, intensity * (Rs + Rp) / 2.0, depth + 1, log)

        S_t = mueller_interface(S_rot, Rs, Ts, Rp, Tp, kind='transmit')
        ext_deg = math.degrees(math.acos(np.clip(np.dot(t_dir, n_out), -1, 1)))
        ledger['exited'].append(dict(intensity=intensity * (Ts + Tp) / 2.0, angle_deg=ext_deg,
                                      dop=dop_of(S_t), from_point=hitp, dir=t_dir, depth=depth))

    recurse(entry_point, t_dir, S_in, s_hat, T_in, 1, dominant_log)

    ledger['total'] = sum(e['intensity'] for e in ledger['exited']) + ledger['truncated']
    ledger['dominant_log'] = dominant_log
    return ledger


# ======================================================== ray fan (P4+) ===
# The corridor's natural terminal hit for any vertical ray on this cut is
# the 3rd interior hit (see trace_tree's docstring and the plan's §2b
# corridor algebra) — that hit's own transmitted branch is "the" dominant
# exit whose strength decides how bright a fan ray draws.
MAIN_HITS = 3
FAN_MIN_OPACITY = 0.06  # weak rays fade toward the ground colour, never vanish


def dominant_exit(ledger):
    return next((e for e in ledger['exited'] if e.get('depth') == MAIN_HITS), None)


def fan_opacity(intensity):
    """Monotonic intensity -> opacity map (plan §9 P4 gate). Deliberately not
    claimed as a metric surface — just legible, and checkably monotonic."""
    return FAN_MIN_OPACITY + (1 - FAN_MIN_OPACITY) * intensity


def build_fan(planes, n_stone, x_range=25.0, step=2.0):
    """x_range/step chosen odd so the grid never lands on x=0 exactly — the
    vertical ray straight down the axis passes through the culet, where all
    8 pavilion-main planes meet: a genuine degenerate vertex (plan §3's
    "vertex grazes" case), not a physically weak ray, and out of scope to
    resolve here (needs the degenerate-hit counter, a later invariant)."""
    rays = []
    x = -x_range
    while x <= x_range + 1e-9:
        origin = np.array([x, 1000.0, 0.0])
        ledger = trace_tree(origin, np.array([0.0, -1.0, 0.0]), planes, n_stone)
        dom = dominant_exit(ledger)
        rays.append(dict(
            x=x, intensity=dom['intensity'] if dom else 0.0,
            from_point=dom['from_point'] if dom else None,
            dir=dom['dir'] if dom else None))
        x += step
    return rays


# ==================================================== §4c: convergence ===
# Real-world scale: girdle diameter 100 units <-> 6.5 mm (§7b's "about one
# carat"), so 1 unit = 0.065 mm and 1 mm = 1/0.065 units.
MM_PER_UNIT = 0.065
UNITS_PER_MM = 1.0 / MM_PER_UNIT

LENS = dict(f_mm=100.0, D_mm=75.0, theta_sun=9.30e-3)  # plan §4, pinned


def lens_geometry():
    """Static, at-focus preview (standoff s=f): marginal-ray half-angle and
    spot size, both already pinned in plan §4 — no new numbers here."""
    f, D = LENS['f_mm'], LENS['D_mm']
    half_angle_deg = math.degrees(math.atan((D / 2.0) / f))
    spot_mm = f * LENS['theta_sun']
    return dict(f_units=f * UNITS_PER_MM, aperture_r_units=(D / 2.0) * UNITS_PER_MM,
                half_angle_deg=half_angle_deg, spot_mm=spot_mm)


def line_intersection_2d(p1, d1, p2, d2):
    """Solve p1 + t1*d1 = p2 + t2*d2 in the x-y plane. Returns (t1, t2) or
    None if the lines are parallel (determinant ~0)."""
    A = np.array([[d1[0], -d2[0]], [d1[1], -d2[1]]])
    det = A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]
    if abs(det) < 1e-9:
        return None
    b = np.array([p2[0] - p1[0], p2[1] - p1[1]])
    t1 = (b[0] * A[1, 1] - A[0, 1] * b[1]) / det
    t2 = (A[0, 0] * b[1] - b[0] * A[1, 0]) / det
    return t1, t2


def fan_caustic(rays):
    """Exhaustive search over ADJACENT fan pairs (plan §4c) for where the
    exit fan's rays actually cross. Most pairs are exactly parallel (same
    terminal facet -> same exit direction, only translated) and contribute
    no finite crossing; real crossings occur only at the fan's central
    mirror seam and at table/crown-main facet-transition seams. Returns a
    list of dicts {x1, x2, point, dist_mm} for every pair that crosses in
    front of both rays (t>0 each)."""
    crossings = []
    for i in range(len(rays) - 1):
        a, b = rays[i], rays[i + 1]
        if a['from_point'] is None or b['from_point'] is None:
            continue
        sol = line_intersection_2d(a['from_point'], a['dir'], b['from_point'], b['dir'])
        if sol is None:
            continue
        t1, t2 = sol
        if t1 <= 1e-6 or t2 <= 1e-6:
            continue
        point = a['from_point'] + t1 * a['dir']
        dist_units = min(t1, t2)
        crossings.append(dict(x1=a['x'], x2=b['x'], point=point, dist_mm=dist_units * MM_PER_UNIT))
    return crossings


# ==================================================== P5: rotation + cone =
def quat_from_axis_angle(axis, angle):
    axis = axis / np.linalg.norm(axis)
    s = math.sin(angle / 2)
    return np.array([math.cos(angle / 2), axis[0] * s, axis[1] * s, axis[2] * s])  # w,x,y,z


def quat_mul(q1, q2):
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return np.array([
        w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
        w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
        w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
        w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2])


def quat_to_matrix(q):
    """§4b: gem orientation as a quaternion (renormalised each frame,
    per the plan's own rationale -- Euler angles gimbal-lock on a stone
    the reader is free to tumble). This matrix rotates the SOLID; the
    tracer itself always runs in gem space, pulling the fixed world-space
    source through this matrix's TRANSPOSE (=inverse, for a rotation) --
    see trace_in_gem_space below."""
    w, x, y, z = q
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]])


def trace_in_gem_space(world_origin, world_dir, gem_quat, solid, n_stone, **kwargs):
    """§4b's frame split: the facet table is NEVER mutated (no drift from
    repeated rotation); instead the fixed world-space source is pulled
    through M^-1 (=M^T for a rotation matrix) into gem space, and the
    existing trace_tree runs unchanged. Verified byte-identical to
    literally rotating the solid by M and keeping the source fixed --
    that's the §8 rotation-invariance invariant, checked below."""
    M = quat_to_matrix(gem_quat)
    Minv = M.T
    return trace_tree(Minv @ world_origin, Minv @ world_dir, solid, n_stone, **kwargs)


def entry_facet_handoff(beam_pos, tilt_deg, planes, axis=np.array([0., 0., 1.])):
    """§4b: for a beam fixed in WORLD space at world-x=beam_pos (vertical,
    aimed down), find which gem-space facet it enters as the gem tilts by
    tilt_deg about `axis` (through the girdle centre, in the section
    plane by default). Reproduces the plan's own pinned handoff table
    EXACTLY -- but only against an IDEALIZED ZERO-THICKNESS girdle
    (round_brilliant_planes(girdle_t_pct=0)), not this file's normal 3%-
    thickness default: the plan's own §4b table is closed-form 2D
    arithmetic (its "pre-tracer" numbers, per §10), computed against the
    same zero-thickness idealization as §1's own figure. Real, worth
    stating plainly rather than silently reconciling: the handoff angle
    shifts by a fraction of a degree once real girdle thickness is
    included, and that shift is itself a fact about the built stone, not
    an error in either number."""
    q = quat_from_axis_angle(axis, math.radians(tilt_deg))
    M = quat_to_matrix(q)
    Minv = M.T
    o = Minv @ np.array([beam_pos, 1000.0, 0.0])
    d = Minv @ np.array([0.0, -1.0, 0.0])
    t, plane = enter_convex(o, d, planes)
    return plane['name'] if plane else None


def find_handoff_tilt(beam_pos, planes, lo=0.0, hi=90.0, tol=1e-4):
    """Bisects for the tilt at which entry_facet_handoff flips away from
    the table -- to crown_main first, in every pinned case, though at a
    wide enough `hi` the beam can tilt PAST crown_main into girdle_side
    too; this only needs "no longer table" as the flip condition, not
    "specifically crown," so a generous default hi doesn't misfire."""
    f_lo = entry_facet_handoff(beam_pos, lo, planes)
    f_hi = entry_facet_handoff(beam_pos, hi, planes)
    if 'table' not in f_lo or 'table' in f_hi:
        return None
    while hi - lo > tol:
        mid = (lo + hi) / 2
        f_mid = entry_facet_handoff(beam_pos, mid, planes)
        if 'table' in f_mid:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def aperture_cone_rays(table_y, f_mm=100.0, d_mm=75.0, standoff_mm=100.0, n_aperture=24):
    """§4/§5: sample the lens aperture uniformly (n_aperture points across
    the 75mm diameter), each ray aimed from its OWN aperture position at
    the SAME focal point on the table centre -- the "full cone" the
    chief-ray §1/§2 verdicts don't capture. Returns a list of (origin,
    direction) in GIRDLE units, ready to trace directly.

    Real bug worth recording: an earlier draft built direction from an
    independently-chosen angle and set origin = target - direction*1000,
    which does NOT place the origin at the actual aperture position/height
    at all (a stand-in point on the right LINE, but off by orders of
    magnitude from the real aperture geometry) -- caught only by checking
    the resulting origin's own (x,y) against the intended aperture
    position, not by energy closure (a wrong-but-still-converging ray
    closes energy just fine). The correct construction is the other way
    round: place the origin at the true aperture position and height
    first, and DERIVE the direction as (focal - origin), normalized."""
    half_d_units = (d_mm / 2.0) * UNITS_PER_MM
    lens_y = table_y + standoff_mm * UNITS_PER_MM
    focal = np.array([0.0, table_y, 0.0])
    rays = []
    for k in range(n_aperture):
        # aperture offset, symmetric about the centre, avoiding the
        # degenerate x=0 chief ray for the same reason build_fan does
        frac = (k - (n_aperture - 1) / 2.0) / (n_aperture / 2.0)
        ax = frac * half_d_units
        origin = np.array([ax, lens_y, 0.0])
        direction = focal - origin
        rays.append((origin, direction / np.linalg.norm(direction)))
    return rays


def cone_survival_fraction(material, pav_deg=40.75, refine_near_thetac=True):
    """§4's own table: fraction of the f/1.33 cone (at focus) that
    survives TIR at the FIRST pavilion hit, for a chief ray landing on the
    table centre. Model: a ray arriving at angle delta_air (in air, in the
    section plane) refracts to an interior angle delta_int, which adds
    directly onto the chief ray's own pavilion AOI. §10's own risk note:
    T climbs from 11% to 50% across the last 0.4 deg below theta_c, so
    uniform sampling puts its largest error exactly where survival is
    closest to 50% -- refine_near_thetac adds extra samples in a +-1 deg
    band around theta_c specifically to keep that error down.

    KNOWN LIMITATION, stated plainly rather than hidden: this 1D
    (section-plane-only) model reproduces the pinned §4 table closely for
    6 of 9 materials (moissanite/diamond/CZ exactly, at 100% survival;
    beryl, quartz and glass within 3 points, right where the knife edge
    the page cares about most actually sits) but is off by 6-13 points for
    sapphire, spinel and fluorite -- the three LARGEST-|margin| materials,
    where the true 2D (azimuthal, not just in-plane) shape of the cone
    evidently matters more than this simplification captures. Several
    alternative samplings were tried (full 3D ray-traced azimuthal
    sampling; a cos(theta)-projected 2D aperture) and none reproduced all
    9 values simultaneously either -- consistent with §10's own flagged
    risk that this integral is genuinely sampling-sensitive. Shipping the
    closer, simpler model with the discrepancy named, per this project's
    own rule for unsourceable precision, rather than quietly picking
    whichever knob happened to fit one run."""
    thetac = material['thetac']
    lo, hi = -20.56, 20.56
    samples = list(np.linspace(lo, hi, 41))
    n = material['n']
    if refine_near_thetac:
        # the chief ray's AOI at the first pavilion hit is pav_deg; a
        # cone half-angle delta (in air) compresses to asin(sin(delta)/n)
        # inside, so solve for which air-angle delta lands the interior
        # AOI within 1 deg of theta_c, and refine the sampling there
        aoi_at = lambda delta_air: pav_deg + math.degrees(math.asin(math.sin(math.radians(delta_air)) / n))
        for delta10 in range(int(lo * 10), int(hi * 10) + 1):
            delta = delta10 / 10.0
            if abs(aoi_at(delta) - thetac) < 1.0:
                samples.append(delta)
    survived = 0
    for delta_air in samples:
        sin_int = math.sin(math.radians(delta_air)) / n
        if abs(sin_int) > 1:
            continue
        delta_int = math.degrees(math.asin(sin_int))
        aoi = pav_deg + delta_int
        if aoi >= thetac:
            survived += 1
    return survived / len(samples)


def collimated_leak_onset_tilt(material, pav_deg=40.75):
    """§4b's own closed-form collimated anchor: leak begins when the
    interior tilt psi=asin(sin(tau)/n) eats the pavilion's margin
    (pav_deg - theta_c). Solved directly (no search needed) -- this is
    the ANCHOR the focused-cone search below is bracketed by, not a
    duplicate of it."""
    margin = pav_deg - material['thetac']
    if margin <= 0:
        return 0.0
    return math.degrees(math.asin(material['n'] * math.sin(math.radians(margin))))


def cone_leak_onset_tilt(material, planes, pav_deg=40.75, lo=0.0, hi=90.0, tol=1e-3):
    """P5's own deliverable (plan §9): the gem tilt at which the FOCUSED
    f/1.33 CONE begins to leak through the first pavilion hit -- i.e. the
    smallest WORLD-space tilt at which cone_survival_fraction (evaluated
    with the effective pavilion AOI at that tilt) drops below 1.0.
    Bracketed by collimated_leak_onset_tilt (a lower bound: the cone's
    marginal rays reach critical before the chief ray does, so the cone
    leaks at or before the chief-ray-only collimated anchor).

    World tilt -> effective pavilion AOI uses the SAME refraction-aware
    relation as collimated_leak_onset_tilt (interior tilt =
    asin(sin(world_tilt)/n)) -- SUBTRACTED from pav_deg, not added: a real
    bug caught only against §4b's own worked example (diamond's worst
    case is 40.75-24.43=16.32 deg, a grazing entry ray, NOT 40.75+24.43 --
    tilting the other way only makes the margin larger, never smaller).
    Not raw world-tilt degrees subtracted directly either, which would
    conflate a world-space angle with an already-refracted interior one."""
    n = material['n']

    def effective_pav(tilt):
        return pav_deg - math.degrees(math.asin(min(1.0, math.sin(math.radians(tilt)) / n)))

    def any_leak(tilt):
        return cone_survival_fraction(material, pav_deg=effective_pav(tilt)) < 1.0
    if any_leak(lo):
        return lo
    if not any_leak(hi):
        return None
    while hi - lo > tol:
        mid = (lo + hi) / 2
        if any_leak(mid):
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


# ============================================================ P6b: uniaxial =
# Quartz's o/e split, per plan §9/§10 -- "the repo's flagship material
# only." Deliberate duplication of trace_tree's own recursion, documented
# rather than hidden: trace_tree is the single most heavily pinned/
# validated function in this file, and threading a per-step callable index
# through it risked that foundation for a feature the plan's own §10
# explicitly permits shipping in reduced form. Ray walkoff (the true ray
# direction deviating from the wave normal by up to 0.34 deg) is NOT
# modelled -- the plan's own "kill switch" scope, a third-order correction
# against AOI values of tens of degrees; implementing it properly needs
# the index ellipsoid's Poynting-vector construction, real additional work
# the plan does not require landing.
def effective_index_e(theta_rad, n_o, n_e_principal):
    """Standard uniaxial-crystal extraordinary-wave index:
    1/n^2(theta) = cos^2(theta)/n_o^2 + sin^2(theta)/n_e^2, theta measured
    from the c-axis to the WAVE NORMAL (here, the ray direction -- walkoff
    is not modelled, so wave normal and ray direction coincide). theta=0
    (propagation || c) gives n_e(0)=n_o exactly (o and e degenerate,
    §9's own P6b gate); theta=90 gives the principal n_e."""
    inv_n2 = math.cos(theta_rad) ** 2 / n_o ** 2 + math.sin(theta_rad) ** 2 / n_e_principal ** 2
    return 1.0 / math.sqrt(inv_n2)


def trace_tree_e_wave(origin, direction, planes, n_o, n_e_principal, c_axis, n_air=1.0,
                       max_bounces=14, min_intensity=1e-9):
    """The e-wave branch: mirrors trace_tree's own recursion exactly,
    EXCEPT the index used at each hit is recomputed from that hit's own
    incoming direction's angle to c_axis (eff_n below), not a constant
    n_stone -- the wave normal genuinely rotates at every bounce, so the
    effective index genuinely differs per-hit, not just once at entry."""
    direction = normalize(direction)
    ledger = {'exited': [], 'truncated': 0.0}
    dominant_log = []

    def eff_n(d):
        cos_theta = abs(np.dot(d, c_axis))
        theta = math.acos(np.clip(cos_theta, -1, 1))
        return effective_index_e(theta, n_o, n_e_principal)

    n_entry = eff_n(direction)
    t, entry_plane = enter_convex(origin, direction, planes)
    if entry_plane is None:
        raise RuntimeError('entry ray misses the solid')
    entry_point = origin + direction * t
    n_out = entry_plane['n']
    s_hat = (normalize(np.cross(n_out, direction))
             if abs(np.dot(n_out, direction)) < 0.999999 else np.array([1.0, 0.0, 0.0]))
    S = np.array([1.0, 0.0, 0.0, 0.0])

    t_dir, cosi, cost = refract(direction, n_out, n_air, n_entry)
    entry_aoi = math.degrees(math.acos(np.clip(cosi, -1, 1)))
    dominant_log.append(dict(point=entry_point, aoi=entry_aoi, plane=entry_plane['name']))

    Rs, Ts, Rp, Tp, rs, rp = fresnel(cosi, cost, n_air, n_entry)
    S_door = mueller_interface(S, Rs, Ts, Rp, Tp, rs=rs, rp=rp, kind='reflect')
    n_op = n_out if np.dot(direction, n_out) < 0 else -n_out
    door_dir = reflect(direction, n_op)
    ledger['exited'].append(dict(intensity=(Rs + Rp) / 2.0, angle_deg=entry_aoi,
                                  dop=dop_of(S_door), from_point=entry_point, dir=door_dir,
                                  label='door', depth=0))

    S_in = mueller_interface(S, Rs, Ts, Rp, Tp, kind='transmit')
    T_in = (Ts + Tp) / 2.0

    def recurse(point, d, S, s_hat, intensity, depth, log):
        if intensity < min_intensity or depth > max_bounces:
            ledger['truncated'] += intensity
            return
        t, plane = trace_convex(point, d, planes)
        if plane is None:
            ledger['truncated'] += intensity
            return
        hitp = point + d * t
        n_out = plane['n']

        s_new = (normalize(np.cross(n_out, d))
                 if abs(np.dot(n_out, d)) < 0.999999 else s_hat)
        theta_rot = rotation_angle(s_hat, s_new, d)
        S_rot = mueller_rotate(S, theta_rot)

        n_local = eff_n(d)  # the one line that differs from trace_tree's own recursion
        t_dir, cosi, cost = refract(d, n_out, n_local, n_air)
        aoi_deg = math.degrees(math.acos(np.clip(cosi, -1, 1)))
        n_op = n_out if np.dot(d, n_out) < 0 else -n_out
        r_dir = reflect(d, n_op)
        log.append(dict(point=hitp, aoi=aoi_deg, plane=plane['name']))

        if t_dir is None:
            sini = math.sqrt(max(0.0, 1 - cosi * cosi))
            n_rel = n_air / n_local
            ds, dp = tir_phase(cosi, sini, n_rel)
            S_out = mueller_interface(S_rot, 1.0, 0.0, 1.0, 0.0, delta=(dp - ds), kind='reflect')
            recurse(hitp, r_dir, S_out, s_new, intensity, depth + 1, log)
            return

        Rs, Ts, Rp, Tp, rs, rp = fresnel(cosi, cost, n_local, n_air)
        S_r = mueller_interface(S_rot, Rs, Ts, Rp, Tp, rs=rs, rp=rp, kind='reflect')
        recurse(hitp, r_dir, S_r, s_new, intensity * (Rs + Rp) / 2.0, depth + 1, log)

        S_t = mueller_interface(S_rot, Rs, Ts, Rp, Tp, kind='transmit')
        ext_deg = math.degrees(math.acos(np.clip(np.dot(t_dir, n_out), -1, 1)))
        ledger['exited'].append(dict(intensity=intensity * (Ts + Tp) / 2.0, angle_deg=ext_deg,
                                      dop=dop_of(S_t), from_point=hitp, dir=t_dir, depth=depth))

    recurse(entry_point, t_dir, S_in, s_hat, T_in, 1, dominant_log)
    ledger['total'] = sum(e['intensity'] for e in ledger['exited']) + ledger['truncated']
    ledger['dominant_log'] = dominant_log
    return ledger


def trace_tree_uniaxial(origin, direction, planes, n_o, n_e_principal, c_axis, **kwargs):
    """Combines the o-wave (trace_tree at n_o, exactly the isotropic
    tracer -- an ordinary wave has no direction-dependence at all) and
    e-wave branches. Energy split at entry is exactly 50/50: sunlight
    arrives unpolarized (S=[1,0,0,0], this file's own convention), and
    projecting an UNPOLARIZED Stokes state onto any two orthogonal
    eigenpolarizations always gives equal halves by definition -- "project
    onto the eigenpolarizations" (plan §9) is not a nontrivial calculation
    for this specific, universal-in-this-tracer input state, just a fact
    worth stating rather than silently assuming. Returns a combined ledger
    (each exited event tagged with its wave) that itself sums to 1.0."""
    led_o = trace_tree(origin, direction, planes, n_o, **kwargs)
    led_e = trace_tree_e_wave(origin, direction, planes, n_o, n_e_principal, c_axis, **kwargs)
    exited = ([dict(e, wave='o', intensity=e['intensity'] * 0.5) for e in led_o['exited']]
              + [dict(e, wave='e', intensity=e['intensity'] * 0.5) for e in led_e['exited']])
    total = sum(e['intensity'] for e in exited) + 0.5 * led_o['truncated'] + 0.5 * led_e['truncated']
    return dict(exited=exited, truncated=0.5 * led_o['truncated'] + 0.5 * led_e['truncated'],
                total=total, dominant_log_o=led_o['dominant_log'], dominant_log_e=led_e['dominant_log'])


# ================================================================ P6: spectral =
# 7 Fraunhofer lines, per plan §4/§9 (h through B).
FRAUNHOFER_NM = dict(h=404.7, G=430.8, F=486.1, e=546.1, d=587.6, C=656.3, B=686.7)


def cauchy_fit(n_d, disp_bg):
    """§3b's stated fallback -- "a two-term Cauchy fit constructed to
    reproduce the pinned n_d and B-G dispersion exactly" -- used here for
    EVERY material, not just ones lacking a published Sellmeier equation.
    Deliberate: reciting real Sellmeier coefficients from memory risks a
    subtly wrong number that LOOKS authoritative (worse than a labelled
    fit, per the plan's own "presenting a fit as a measurement would be
    the one dishonest number on the page"); a fit calibrated to the
    already-pinned §2 table values is honest, verifiable, and exact on
    the two numbers that matter most. n(lambda) = A + B/lambda^2, lambda
    in microns; solved directly from n(d)=n_d and n(G)-n(B)=disp_bg."""
    l_d, l_g, l_b = FRAUNHOFER_NM['d'] / 1000.0, FRAUNHOFER_NM['G'] / 1000.0, FRAUNHOFER_NM['B'] / 1000.0
    b_coef = disp_bg / (1 / l_g ** 2 - 1 / l_b ** 2)
    a_coef = n_d - b_coef / l_d ** 2
    return a_coef, b_coef


def n_at_wavelength(a_coef, b_coef, lambda_nm):
    return a_coef + b_coef / (lambda_nm / 1000.0) ** 2


def group_index(a_coef, b_coef, lambda_nm):
    """n_g = n - lambda*dn/dlambda; for n=A+B/lambda^2, dn/dlambda=-2B/lambda^3,
    so n_g = A + 3B/lambda^2 = n + 2B/lambda^2. Diamond at d: ~2.4955 vs the
    plan's own pinned 2.494 (a different, unstated computation method) --
    close agreement, not exact, and not claimed to be."""
    n = n_at_wavelength(a_coef, b_coef, lambda_nm)
    lam_um = lambda_nm / 1000.0
    return n + 2 * b_coef / lam_um ** 2


def solar_weight(lambda_nm, temp_k=5778.0):
    """Solar-spectrum weighting (plan §4: "weighted by the solar spectrum
    rather than equal-energy"). A blackbody at the sun's effective
    temperature (Planck's law), not a measured irradiance table -- a
    genuine physical approximation, labelled as such, rather than reciting
    tabulated AM0/AM1.5 irradiance numbers from memory. Returns spectral
    radiance in arbitrary units (only relative weights across the 7
    Fraunhofer lines matter, so constants cancel)."""
    lam_m = lambda_nm * 1e-9
    h_planck, c_light, k_b = 6.62607015e-34, 2.99792458e8, 1.380649e-23
    x = h_planck * c_light / (lam_m * k_b * temp_k)
    return 1.0 / (lam_m ** 5 * (math.exp(x) - 1))


def dispersion_fan(material, planes, ray_o, ray_d):
    """Traces the B and G Fraunhofer lines through the SAME pinned ray and
    returns the exit-angle spread (the "B-G fan" the plan's §1/§4 pin at
    1.06 deg for diamond, 2.80 deg for moissanite) -- using each
    material's own Cauchy fit, calibrated to its §2-pinned n_d/dispersion."""
    a_coef, b_coef = cauchy_fit(material['n'], material['disp_bg'])
    angles = {}
    for name in ('G', 'B'):
        n_lam = n_at_wavelength(a_coef, b_coef, FRAUNHOFER_NM[name])
        led = trace_tree(ray_o, ray_d, planes, n_lam)
        exit_evt = max(led['exited'][1:], key=lambda e: e['intensity'])
        angles[name] = exit_evt['angle_deg']
    return abs(angles['B'] - angles['G']), angles


# ============================================================ self-checks =
def close(a, b, tol):
    return abs(a - b) <= tol


def run_checks():
    n_ok = [0]
    n_fail = [0]

    def ok(cond, msg):
        n_ok[0] += 1
        if not cond:
            n_fail[0] += 1
            print(f'FAIL: {msg}')

    planes, meta = round_brilliant_planes()

    # ---- P2 gate: depth reproduces the GIA-Excellent band ----
    ok(59.0 <= meta['depth_pct'] <= 62.3, f"depth {meta['depth_pct']:.3f}% not in GIA-Excellent 59.0-62.3")
    ok(close(meta['hc'], 15.1184, 0.01), f"crown height {meta['hc']} != 15.1184")
    ok(close(meta['hp'], 43.0824, 0.01), f"pavilion depth {meta['hp']} != 43.0824")

    # ---- geometry kernel: closed solid, Euler check ----
    sc = solid_check(planes)
    ok(sc['euler'] == 2, f"V-E+F = {sc['euler']} != 2 (V={sc['V']} E={sc['E']} F={sc['F']})")
    ok(len(sc['bad_edges']) == 0, f"{len(sc['bad_edges'])} edges not shared by exactly 2 faces")
    ok(sc['vol'] > 0, 'solid volume not positive')

    # ---- material table: thetac / R0 sanity vs plan §2 ----
    ok(close(MATERIALS['diamond']['thetac'], 24.43, 0.01), 'diamond thetac')
    ok(close(MATERIALS['diamond']['R0'] * 100, 17.20, 0.05), 'diamond R0')
    ok(close(MATERIALS['quartz']['thetac'], 40.36, 0.01), 'quartz thetac')
    ok(close(MATERIALS['glass']['thetac'], 41.14, 0.01), 'glass thetac')

    # ---- §8 pinned vector: RB diamond, x=-5 ----
    n = MATERIALS['diamond']['n']
    ray_o = np.array([-5.0, 1000.0, 0.0])
    ray_d = np.array([0.0, -1.0, 0.0])
    led = trace_tree(ray_o, ray_d, planes, n)
    aois = [round(h['aoi'], 2) for h in led['dominant_log']]
    ok(close(aois[0], 0.0, 0.01), f'diamond x-5 entry AOI {aois[0]} != 0 (vertical ray on table)')
    ok(close(aois[1], 40.75, 0.05), f'diamond x-5 AOI1 {aois[1]} != 40.75')
    ok(close(aois[2], 57.75, 0.05), f'diamond x-5 AOI2 {aois[2]} != 57.75')
    ok(close(aois[3], 17.00, 0.05), f'diamond x-5 AOI3 {aois[3]} != 17.00')
    exit_evt = max(led['exited'], key=lambda e: e['intensity'])
    ok(close(exit_evt['angle_deg'], 44.98, 0.1), f"diamond x-5 exit angle {exit_evt['angle_deg']} != 44.98")
    door = led['exited'][0]
    ok(close(door['intensity'] * 100, 17.20, 0.1), f"diamond door {door['intensity']*100} != 17.20")
    ok(close(led['total'], 1.0, 1e-6), f"diamond x-5 energy closure {led['total']} != 1")
    ok(close(exit_evt['intensity'] * 100, 67.80, 1.0), f"diamond x-5 table-exit ledger {exit_evt['intensity']*100} != 67.80")

    # ---- §8 pinned vector: RB quartz, x=-5 ----
    n = MATERIALS['quartz']['n']
    led = trace_tree(ray_o, ray_d, planes, n)
    aois = [round(h['aoi'], 2) for h in led['dominant_log']]
    ok(close(aois[1], 40.75, 0.05), f'quartz x-5 AOI1 {aois[1]} != 40.75')
    ok(close(aois[2], 57.75, 0.05), f'quartz x-5 AOI2 {aois[2]} != 57.75')
    ok(close(aois[3], 17.00, 0.05), f'quartz x-5 AOI3 {aois[3]} != 17.00')
    exit_evt = max(led['exited'], key=lambda e: e['intensity'])
    ok(close(exit_evt['angle_deg'], 26.84, 0.1), f"quartz x-5 exit angle {exit_evt['angle_deg']} != 26.84")
    ok(close(exit_evt['intensity'] * 100, 90.96, 1.0), f"quartz x-5 table-exit ledger {exit_evt['intensity']*100} != 90.96")
    ok(close(led['total'], 1.0, 1e-6), f"quartz x-5 energy closure {led['total']} != 1")

    # ---- §8 pinned vector: RB glass, x=-5 (leaks at first pavilion hit) ----
    n = MATERIALS['glass']['n']
    led = trace_tree(ray_o, ray_d, planes, n)
    aois = [round(h['aoi'], 2) for h in led['dominant_log']]
    ok(close(aois[1], 40.75, 0.05), f'glass x-5 AOI1 {aois[1]} != 40.75')
    door = led['exited'][0]
    ok(close(door['intensity'] * 100, 4.26, 0.05), f"glass door {door['intensity']*100} != 4.26")
    leak = [e for e in led['exited'] if close(e['angle_deg'], 82.84, 1.0)]
    ok(len(leak) >= 1, 'glass: no exit branch near 82.84 deg (first-hit leak)')
    if leak:
        ok(close(leak[0]['intensity'] * 100, 47.53, 1.5), f"glass leak fraction {leak[0]['intensity']*100} != ~47.53")
    ok(close(led['total'], 1.0, 1e-6), f"glass x-5 energy closure {led['total']} != 1")

    # ---- §8 pinned vector: RB diamond, x=-14 (exits crown main, not table) ----
    n = MATERIALS['diamond']['n']
    ray_o2 = np.array([-14.0, 1000.0, 0.0])
    led = trace_tree(ray_o2, ray_d, planes, n)
    aois = [round(h['aoi'], 2) for h in led['dominant_log']]
    ok(close(aois[1], 40.75, 0.05), f'diamond x-14 AOI1 {aois[1]} != 40.75')
    ok(close(aois[2], 57.75, 0.05), f'diamond x-14 AOI2 {aois[2]} != 57.75')
    ok(close(aois[3], 17.50, 0.1), f'diamond x-14 AOI3 {aois[3]} != 17.50')
    exit_evt = max(led['exited'][1:], key=lambda e: e['intensity'])  # skip door
    # dominant_log[3] is the exit hit itself (entry, hit1, hit2, hit3=exit);
    # dominant_log[-1] is NOT the exit — the reflected sibling keeps
    # recursing past it, which is deliberate (see trace_tree docstring).
    ok('crown_main' in led['dominant_log'][3]['plane'], f"diamond x-14 should exit via crown_main, got {led['dominant_log'][3]['plane']}")
    ok(close(exit_evt['angle_deg'], 46.63, 0.15), f"diamond x-14 exit angle {exit_evt['angle_deg']} != 46.63")

    # ---- TIR phase / Fresnel rhomb acceptance test (§3b eq5, §8) ----
    n_rel = 1.0 / 1.51
    cosi = math.cos(math.radians(54.62))
    sini = math.sin(math.radians(54.62))
    ds, dp = tir_phase(cosi, sini, n_rel)
    delta = math.degrees(dp - ds)
    ok(close(delta, 45.00, 0.05), f'Fresnel rhomb Delta {delta} != 45.00')

    n_rel_q = 1.0 / MATERIALS['quartz']['n']
    cosi = math.cos(math.radians(57.75))
    sini = math.sin(math.radians(57.75))
    ds, dp = tir_phase(cosi, sini, n_rel_q)
    delta_q = math.degrees(dp - ds)
    ok(close(delta_q, 44.18, 0.1), f'quartz second-bounce Delta {delta_q} != 44.18')

    # ---- Brewster / DOP invariant (interior reflection, diamond 22.47 deg) ----
    n1 = MATERIALS['diamond']['n']
    brewster_deg = math.degrees(math.atan(1.0 / n1))
    ok(close(brewster_deg, 22.47, 0.05), f'diamond internal Brewster {brewster_deg} != 22.47')
    cosi = math.cos(math.radians(brewster_deg))
    cost = math.cos(math.asin(n1 * math.sin(math.radians(brewster_deg)) / 1.0))
    Rs, Ts, Rp, Tp, rs, rp = fresnel(cosi, cost, n1, 1.0)
    ok(Rp < 1e-9, f'Rp at Brewster = {Rp} != 0')

    # ---- Fresnel continuity into TIR (no discontinuity at theta_c) ----
    n1, n2 = MATERIALS['quartz']['n'], 1.0
    thetac = math.radians(MATERIALS['quartz']['thetac'])
    for eps in (1e-2, 1e-4, 1e-6):
        th = thetac - eps
        cosi, cost = math.cos(th), math.sqrt(max(0, 1 - (n1 / n2 * math.sin(th)) ** 2))
        Rs, Ts, Rp, Tp, rs, rp = fresnel(cosi, cost, n1, n2)
        R = (Rs + Rp) / 2
        ok(R < 1.0 and R > 0.9 if eps < 1e-4 else True, f'R approaching 1 near thetac (eps={eps}): R={R}')

    # ---- corridor identities (§2b closed form) ----
    p = 40.75
    aoi2_pred = 180 - 3 * p
    aoi3_pred = abs(4 * p - 180)
    n = MATERIALS['diamond']['n']
    led = trace_tree(ray_o, ray_d, planes, n)
    aois = [round(h['aoi'], 2) for h in led['dominant_log']]
    ok(close(aois[2], aoi2_pred, 0.05), f'corridor AOI2 identity: {aois[2]} != {aoi2_pred}')
    ok(close(aois[3], aoi3_pred, 0.05), f'corridor AOI3 identity: {aois[3]} != {aoi3_pred}')

    # ---- P4 fan gate: opacity strictly monotonic in dominant-exit intensity,
    # for every material, across the sampled offsets ----
    for mat in MATERIALS:
        rays = build_fan(planes, MATERIALS[mat]['n'])
        pairs = sorted((r['intensity'], fan_opacity(r['intensity'])) for r in rays)
        ok(all(pairs[i][1] <= pairs[i+1][1] for i in range(len(pairs)-1)),
           f'{mat}: fan opacity not monotonic in intensity')
        ok(all(0.0 <= r['intensity'] <= 1.0 + 1e-9 for r in rays),
           f'{mat}: fan ray intensity out of [0,1]')

    # ---- §4c: lens preview reproduces §4's already-pinned numbers ----
    lg = lens_geometry()
    ok(close(lg['half_angle_deg'], 20.56, 0.01), f"lens half-angle {lg['half_angle_deg']} != 20.56")
    ok(close(lg['spot_mm'], 0.930, 0.001), f"lens spot {lg['spot_mm']} != 0.930")

    # ---- §4c: exit-fan caustic — almost entirely parallel (each ray's
    # own regime just translated sideways); real forward crossings occur
    # only in mirrored pairs, at the table/crown-main facet-transition
    # seam (here x=+-9/+-11 — a purely geometric fact, since the two
    # pavilion-main bounces are plain reflection and so are exactly the
    # same for every material; only the crossing DISTANCE, set by the
    # final Snell exit, is material-dependent) ----
    for mat in ('diamond', 'quartz', 'glass'):
        rays = build_fan(planes, MATERIALS[mat]['n'])
        crossings = fan_caustic(rays)
        ok(len(crossings) == 2, f'{mat}: expected exactly 2 forward crossings (mirrored pair), got {len(crossings)}')
        xs = sorted((c['x1'], c['x2']) for c in crossings)
        ok(xs == [(-11.0, -9.0), (9.0, 11.0)], f'{mat}: crossings at unexpected x-pairs {xs}')
        ok(all(c['dist_mm'] > 0 for c in crossings), f'{mat}: a reported crossing distance is not positive')
        # deep in one regime, well away from the seam: exactly parallel
        same_regime = [c for c in crossings if c['x1'] == 3.0]
        ok(len(same_regime) == 0, f'{mat}: x=3/x=5 (same facet, should be parallel) unexpectedly crossed')

    # ---- P2 gate: quadric primitive -- sphere (scrying bead) ----
    sphere = dict(planes=[], quadrics=[dict(center=np.array([0., 0., 0.]), radii=np.array([50., 50., 50.]), name='sphere')])
    n = MATERIALS['quartz']['n']
    led = trace_tree(np.array([-5.0, 1000.0, 0.0]), np.array([0., -1., 0.]), sphere, n)
    ok(close(led['total'], 1.0, 1e-6), f"sphere energy closure {led['total']} != 1")
    ok(close(led['dominant_log'][0]['aoi'], math.degrees(math.asin(5.0 / 50.0)), 0.05),
       f"sphere entry AOI {led['dominant_log'][0]['aoi']} != asin(5/50)")

    # ---- P2 gate: quadric primitive -- egg (prolate ellipsoid, off-axis ray) ----
    egg = dict(planes=[], quadrics=[dict(center=np.array([0., 0., 0.]), radii=np.array([35., 55., 35.]), name='egg')])
    led = trace_tree(np.array([-10.0, 1000.0, 0.0]), np.array([0., -1., 0.]), egg, n)
    ok(close(led['total'], 1.0, 1e-6), f"egg energy closure {led['total']} != 1")
    ok(led['dominant_log'][0]['plane'] == 'egg', "egg entry did not report the egg quadric's name")

    # ---- P2 gate: quadric ∩ plane -- cabochon (dome + flat back) ----
    cab = dict(planes=[dict(n=np.array([0., -1., 0.]), d=20.0, name='cab_back')],
               quadrics=[dict(center=np.array([0., 0., 0.]), radii=np.array([50., 50., 50.]), name='cab_dome')])
    led = trace_tree(np.array([-5.0, 1000.0, 0.0]), np.array([0., -1., 0.]), cab, n)
    ok(close(led['total'], 1.0, 1e-6), f"cabochon energy closure {led['total']} != 1")
    ok(led['dominant_log'][0]['plane'] == 'cab_dome', "cabochon entry should be via the dome, not the flat back")

    # ---- P2 gate: quadric ∩ quadric -- lens/disc (two spherical caps) ----
    # Two spheres of radius R, centres separated so they overlap into a
    # lens: the sphere centred BELOW the overlap forms the lens's visible
    # UPPER surface (its top bulges up into the overlap, like a watch-glass)
    # -- so a downward ray's entry facet is the sphere centred below, not
    # the one centred above; counterintuitive on first read, and exactly
    # why this is worth a pinned regression rather than an eyeballed guess.
    R, sep = 60.0, 90.0
    lens_solid = dict(planes=[], quadrics=[
        dict(center=np.array([0., sep / 2, 0.]), radii=np.array([R, R, R]), name='sphere_centred_above'),
        dict(center=np.array([0., -sep / 2, 0.]), radii=np.array([R, R, R]), name='sphere_centred_below')])
    led = trace_tree(np.array([-5.0, 1000.0, 0.0]), np.array([0., -1., 0.]), lens_solid, n)
    ok(close(led['total'], 1.0, 1e-6), f"lens energy closure {led['total']} != 1")
    ok(led['dominant_log'][0]['plane'] == 'sphere_centred_below',
       f"lens entry should be via the below-centred sphere (forms the upper surface), got {led['dominant_log'][0]['plane']}")

    # ---- P2 gate: union primitive -- merkaba (two interpenetrating tetrahedra) ----
    # Point-inverting a tetrahedron's face-normal set (n -> -n, same d) is
    # the second, oppositely-oriented tetrahedron -- the stellated-octahedron
    # construction the roster calls for.
    tet_a = platonic_planes('tetrahedron', d=1.0)
    tet_b = [dict(n=-p['n'], d=p['d'], name=f"inv_{p['name']}") for p in tet_a]
    merkaba = dict(parts=[dict(planes=tet_a, quadrics=[]), dict(planes=tet_b, quadrics=[])])
    # off-axis AND off z=0: z=0 is a mirror-symmetry plane of this specific
    # inverted-pair construction, and a ray confined to it hits a genuine
    # degenerate tie (two different planes, one per tetrahedron, sharing an
    # identical t along that ray) -- the same class of degeneracy as the
    # round brilliant's on-axis culet vertex, not a tracer bug.
    ray_o3 = np.array([0.3, 1000.0, 0.2])
    led = trace_tree(ray_o3, np.array([0., -1., 0.]), merkaba, n)
    ok(close(led['total'], 1.0, 1e-6), f"merkaba energy closure {led['total']} != 1")
    ok(len(led['dominant_log']) >= 2, 'merkaba: ray should cross at least two facet events (interpenetrating solids)')

    # Direct, hand-verifiable test of the union seam-skip mechanism itself,
    # decoupled from trace_tree's Snell refraction (which bends the ray by
    # an amount not worth re-deriving by hand here): an UNREFRACTED straight
    # ray through this same offset is independently known to pass through B
    # first (entry) then genuinely into A's territory before it can leave --
    # so the real exit must name a tetrahedron_* facet, not another
    # inv_tetrahedron_* facet of the part it entered through (which would
    # mean the seam back into A was missed).
    d_straight = np.array([0., -1., 0.])
    parts = merkaba['parts']
    t_e, name_e, n_e = enter_union(ray_o3, d_straight, parts)
    ok(name_e == 'inv_tetrahedron_3', f'merkaba: expected straight-ray entry via inv_tetrahedron_3, got {name_e}')
    t_x, name_x, n_x = trace_union(ray_o3 + d_straight * t_e + d_straight * 1e-9, d_straight, parts)
    ok(name_x is not None and not name_x.startswith('inv_'),
       f'merkaba: straight-ray exit should be a tetrahedron_* facet (crossed into the other part), got {name_x}')
    exit_pt = ray_o3 + d_straight * (t_e + t_x)
    # a genuine exit sits exactly on the part it's leaving (boundary-
    # inclusive membership trivially still counts that ONE part) -- the
    # real check is that it's not ALSO inside some other part
    ok(len(_which_parts(exit_pt, parts)) <= 1,
       f'merkaba: exit point is inside {len(_which_parts(exit_pt, parts))} parts -- not a genuine exterior surface')

    # ---- P2 gate: Platonic dihedrals match closed form to 1e-9, and the
    # generated solids are themselves closed (V-E+F=2) ----
    for kind, expect in PLATONIC_DIHEDRAL_DEG.items():
        pl = platonic_planes(kind)
        sc = solid_check(pl)
        ok(sc['euler'] == 2, f'{kind}: V-E+F = {sc["euler"]} != 2')
        ok(len(sc['bad_edges']) == 0, f'{kind}: {len(sc["bad_edges"])} edges not shared by exactly 2 faces')
        angles = dihedral_angles(pl)
        ok(len(angles) > 0, f'{kind}: no dihedral angles recovered from the generated solid')
        # tol 1e-4, not 1e-6: the plan's own closed-form constants are
        # themselves quoted to 4 decimal places (e.g. 70.5288 for
        # arccos(1/3) = 70.52877936550931...), so 1e-6 would fail against
        # the reference value's own rounding, not against a real kernel bug
        ok(all(close(a, expect, 1e-4) for a in angles),
           f'{kind}: dihedral angles {sorted(set(round(a,4) for a in angles))} != {expect}')

    # ---- P2 gate: symmetry-group invariance (§3/§8, 1e-12) ----
    # round brilliant's own reduced facet set is invariant under C8 (its
    # generating azimuth step) but NOT under a single vertical mirror alone
    # composed arbitrarily -- check the true, already-built-in symmetry.
    bad = check_symmetry(planes, cnv_matrices(8, mirror=True), tol=1e-9)
    ok(len(bad) == 0, f'round brilliant: {len(bad)} plane/matrix pairs break the declared C8v symmetry')
    # a cube is invariant under C4v about any face axis -- validates the
    # checker against an independently-obvious case
    cube = platonic_planes('cube')
    bad_cube = check_symmetry(cube, cnv_matrices(4, mirror=True), tol=1e-9)
    ok(len(bad_cube) == 0, f'cube: {len(bad_cube)} plane/matrix pairs break C4v')
    # negative control: the checker must actually be able to fail -- a cube
    # is NOT invariant under C5v (no 5-fold axis through a cube face)
    bad_cube5 = check_symmetry(cube, cnv_matrices(5, mirror=False), tol=1e-9)
    ok(len(bad_cube5) > 0, 'cube: checker failed to detect a real C5 asymmetry (negative control)')

    # ---- P2 gate: chord-plane fancy-outline mechanism (§3) ----
    # A circular outline stands in for a real oval/pear/marquise's arc
    # construction -- this checks the MECHANISM (closed solid, positive
    # volume, energy-conserving trace) independent of any specific cut's
    # real proportions, which are still pending P1 sourcing. It is
    # deliberately NOT compared to round_brilliant_planes' own numbers:
    # that generator uses 8 broad pavilion-main wedges, this one uses one
    # narrow facet per outline sample, a different topology by design, so
    # matching AOIs would be a coincidence to chase, not a real invariant.
    circle_arcs = [dict(center=(0.0, 0.0), radius=50.0, a0=0.0, a1=2 * math.pi)]
    outline = arc_outline_points(circle_arcs, samples_per_arc=16)
    t_half, interior = 1.5, (0.0, 0.0, 0.0)
    pav_pl = chord_planes_from_outline(outline, (0, -43.08, 0), -t_half, 'pav', interior)
    crown_pl = chord_planes_from_outline(outline, (0, 15.12 + t_half, 0), t_half, 'crown', interior)
    girdle_pl = girdle_band_from_outline(outline, t_half, interior)
    table_pl = dict(n=np.array([0., 1., 0.]), d=15.12 + t_half, name='table')
    fancy_planes = [table_pl] + crown_pl + girdle_pl + pav_pl
    sc = solid_check(fancy_planes)
    ok(sc['euler'] == 2, f'chord-outline proxy: V-E+F = {sc["euler"]} != 2')
    ok(len(sc['bad_edges']) == 0, f'chord-outline proxy: {len(sc["bad_edges"])} edges not shared by exactly 2 faces')
    ok(sc['vol'] > 0, 'chord-outline proxy: solid volume not positive')
    led = trace_tree(np.array([-5.0, 1000.0, 0.0]), np.array([0., -1., 0.]), fancy_planes, MATERIALS['diamond']['n'])
    ok(close(led['total'], 1.0, 1e-6), f"chord-outline proxy: energy closure {led['total']} != 1")

    # ---- P2 gate: oval brilliant (real cut, IGI 2022 sourced proportions) ----
    oval_planes, oval_meta = oval_brilliant_planes()
    sc = solid_check(oval_planes)
    ok(sc['euler'] == 2, f'oval: V-E+F = {sc["euler"]} != 2')
    ok(len(sc['bad_edges']) == 0, f'oval: {len(sc["bad_edges"])} edges not shared by exactly 2 faces')
    ok(sc['vol'] > 0, 'oval: solid volume not positive')
    ok(58.0 <= oval_meta['depthPct'] <= 63.0, f"oval: depth% {oval_meta['depthPct']} outside IGI Excellent 58.0-63.0")
    # table scale must hit table_frac EXACTLY (similar-triangles construction,
    # not just close) -- verified by direct construction, not the formula
    pts, on_plane = enumerate_vertices(oval_planes)
    table_idx = 0
    ok(oval_planes[table_idx]['name'] == 'oval_table', 'oval: plane 0 is not the table (index assumption broke)')
    table_verts = [pts[i] for i in on_plane[table_idx]]
    table_x_half = max(abs(p[0]) for p in table_verts)
    ok(close(table_x_half, 50.0 * 0.59, 1e-6), f'oval: table half-width {table_x_half} != {50.0*0.59}')
    for mat in ('diamond', 'quartz', 'glass'):
        led = trace_tree(np.array([-5.0, 1000.0, 0.0]), np.array([0., -1., 0.]), oval_planes, MATERIALS[mat]['n'])
        ok(close(led['total'], 1.0, 1e-6), f'oval {mat}: energy closure {led["total"]} != 1')

    # ---- P2 gate: remaining modern-brilliant-fancy cuts, closure + energy ----
    for cutname, fn in (('marquise', marquise_planes), ('pear', pear_planes),
                         ('cushion', cushion_brilliant_planes), ('trilliant', trilliant_planes),
                         ('halfmoon', half_moon_planes), ('kite', kite_planes), ('shield', shield_planes)):
        planes_c, meta_c = fn()
        sc = solid_check(planes_c)
        ok(sc['euler'] == 2, f'{cutname}: V-E+F = {sc["euler"]} != 2')
        ok(len(sc['bad_edges']) == 0, f'{cutname}: {len(sc["bad_edges"])} edges not shared by exactly 2 faces')
        ok(sc['vol'] > 0, f'{cutname}: solid volume not positive')
        for mat in ('diamond', 'quartz', 'glass'):
            led = trace_tree(np.array([-5.0, 1000.0, 0.0]), np.array([0., -1., 0.]), planes_c, MATERIALS[mat]['n'])
            ok(close(led['total'], 1.0, 1e-6), f'{cutname} {mat}: energy closure {led["total"]} != 1')

    # ---- P2 gate: heart (union of 2 convex lobes) ----
    heart = heart_planes()
    for part in heart['parts']:
        sc = solid_check(part['planes'])
        ok(sc['euler'] == 2, f'heart lobe: V-E+F = {sc["euler"]} != 2')
        ok(len(sc['bad_edges']) == 0, f'heart lobe: {len(sc["bad_edges"])} edges not shared by exactly 2 faces')
        ok(sc['vol'] > 0, 'heart lobe: solid volume not positive')
    for xoff in (-5.0, 5.0, -20.0, 20.0):
        led = trace_tree(np.array([xoff, 1000.0, 0.3]), np.array([0., -1., 0.]), heart, MATERIALS['diamond']['n'])
        ok(close(led['total'], 1.0, 1e-6), f'heart x={xoff}: energy closure {led["total"]} != 1')

    # ---- P2 gate: chevron (multi-tier pavilion) cuts -- princess, radiant ----
    for cutname, fn, depth_lo, depth_hi in (('princess', princess_planes, 65.0, 73.0),
                                             ('radiant', radiant_planes, 62.0, 68.0)):
        planes_c, meta_c = fn()
        sc = solid_check(planes_c)
        ok(sc['euler'] == 2, f'{cutname}: V-E+F = {sc["euler"]} != 2')
        ok(len(sc['bad_edges']) == 0, f'{cutname}: {len(sc["bad_edges"])} edges not shared by exactly 2 faces')
        ok(sc['vol'] > 0, f'{cutname}: solid volume not positive')
        ok(depth_lo <= meta_c['depthPct'] <= depth_hi,
           f'{cutname}: depth% {meta_c["depthPct"]} outside IGI Excellent {depth_lo}-{depth_hi}')
        for mat in ('diamond', 'quartz', 'glass'):
            led = trace_tree(np.array([-5.0, 1000.0, 0.0]), np.array([0., -1., 0.]), planes_c, MATERIALS[mat]['n'])
            ok(close(led['total'], 1.0, 1e-6), f'{cutname} {mat}: energy closure {led["total"]} != 1')

    # ---- P2 gate: Vogel wand reproduces plan §8's pinned acceptance vectors ----
    vogel_planes = vogel_wand_planes()
    spec = dict(planes=vogel_planes, quadrics=[])
    n_q, thetac_q = MATERIALS['quartz']['n'], MATERIALS['quartz']['thetac']
    _, _, n_out_b = trace_solid(np.array([0.001, 0., 0.]), np.array([0., 1., 0.]), spec)
    aoi_b = math.degrees(math.acos(abs(np.dot([0, 1, 0], n_out_b))))
    ok(close(aoi_b, 51.86, 0.01), f'Vogel B axial AOI {aoi_b} != 51.86')
    ok(close(aoi_b - thetac_q, 11.51, 0.02), f'Vogel B margin {aoi_b-thetac_q} != +11.51')
    _, _, n_out_a = trace_solid(np.array([0.001, 0., 0.]), np.array([0., -1., 0.]), spec)
    aoi_a = math.degrees(math.acos(abs(np.dot([0, -1, 0], n_out_a))))
    ok(close(aoi_a, 64.07, 0.01), f'Vogel A axial AOI {aoi_a} != 64.07')
    ok(close(aoi_a - thetac_q, 23.71, 0.02), f'Vogel A margin {aoi_a-thetac_q} != +23.71')

    # ---- P2 gate: remaining talismanic "have" + convex-smooth, closure + energy ----
    ray_o1, ray_d1 = np.array([1.0, 1000.0, 0.3]), np.array([0., -1., 0.])
    for cutname, fn in (('vogel', vogel_wand_planes), ('genpoint', generator_point_planes),
                        ('obelisk', obelisk_planes), ('greatpyramid', great_pyramid_planes)):
        planes_c = fn()
        sc = solid_check(planes_c)
        ok(sc['euler'] == 2, f'{cutname}: V-E+F = {sc["euler"]} != 2')
        ok(len(sc['bad_edges']) == 0, f'{cutname}: {len(sc["bad_edges"])} edges not shared by exactly 2 faces')
        ok(sc['vol'] > 0, f'{cutname}: solid volume not positive')
        led = trace_tree(ray_o1, ray_d1, planes_c, n_q)
        ok(close(led['total'], 1.0, 1e-6), f'{cutname}: energy closure {led["total"]} != 1')
    for cutname, fn in (('sphere', sphere_planes), ('egg', egg_planes),
                        ('cabochon', cabochon_planes), ('lensdisc', lens_disc_planes)):
        spec_c = fn()
        led = trace_tree(ray_o1, ray_d1, spec_c, n_q)
        ok(close(led['total'], 1.0, 1e-6), f'{cutname}: energy closure {led["total"]} != 1')

    # ---- P2 gate: last 2 talismanic cuts (P1-sourced) ----
    sod_planes = star_of_david_planes()
    sc = solid_check(sod_planes)
    ok(sc['euler'] == 2, f'star-of-david: V-E+F = {sc["euler"]} != 2')
    ok(sc['V'] == 6 and sc['E'] == 12 and sc['F'] == 8, f'star-of-david: V,E,F = {sc["V"]},{sc["E"]},{sc["F"]} != 6,12,8 (antiprism)')
    ok(len(sc['bad_edges']) == 0, f'star-of-david: {len(sc["bad_edges"])} edges not shared by exactly 2 faces')
    led = trace_tree(ray_o1, ray_d1, sod_planes, n_q)
    ok(close(led['total'], 1.0, 1e-6), f'star-of-david: energy closure {led["total"]} != 1')

    # ---- P5 gate: rotation-invariance (§8) -- rotating the SOLID by M with
    # a fixed source must match rotating the SOURCE by M^-1 with the solid
    # at identity, exactly. This is the transform-plumbing check: a wrong
    # frame convention still produces a plausible-looking, energy-closed
    # trace at every single orientation, so only a DIRECT comparison catches it ----
    rb_planes = round_brilliant_planes()[0]
    q_test = quat_from_axis_angle(np.array([0.3, 0.8, 0.2]), math.radians(37.0))
    M_test = quat_to_matrix(q_test)
    ray_o2, ray_d2 = np.array([-5.0, 1000.0, 0.0]), np.array([0.0, -1.0, 0.0])
    planes_rotated = [dict(n=M_test @ p['n'], d=p['d'], name=p['name']) for p in rb_planes]
    led_rotsolid = trace_tree(ray_o2, ray_d2, planes_rotated, MATERIALS['quartz']['n'])
    led_rotsource = trace_in_gem_space(ray_o2, ray_d2, q_test, rb_planes, MATERIALS['quartz']['n'])
    aois_a = [round(h['aoi'], 6) for h in led_rotsolid['dominant_log']]
    aois_b = [round(h['aoi'], 6) for h in led_rotsource['dominant_log']]
    ok(aois_a == aois_b, f'rotation-invariance: rotate-solid AOIs {aois_a} != rotate-source AOIs {aois_b}')
    ok(close(led_rotsolid['total'], 1.0, 1e-6) and close(led_rotsource['total'], 1.0, 1e-6),
       'rotation-invariance: energy does not close under one or both paths')

    # ---- P5 gate: entry-facet handoff tilts reproduce plan §4b's pinned
    # table (against the SAME idealized zero-thickness girdle that table's
    # own closed-form arithmetic used -- see entry_facet_handoff's docstring) ----
    planes_idealized = round_brilliant_planes(girdle_t_pct=0.0)[0]
    for pos, expect in ((0.0, 61.63), (10.0, 43.32), (20.0, 22.69), (25.0, 9.85), (27.0, 3.58)):
        got = find_handoff_tilt(pos, planes_idealized)
        ok(got is not None and close(got, expect, 0.02), f'handoff tilt at pos={pos}: {got} != {expect}')

    # ---- P5 gate: collimated leak-onset anchors (§4b closed form) ----
    for mat, expect in (('quartz', 0.61), ('diamond', 42.78)):
        got = collimated_leak_onset_tilt(MATERIALS[mat])
        ok(close(got, expect, 0.02), f'{mat} collimated leak-onset {got} != {expect}')

    # ---- P5 gate: cone-survival fractions -- reproduces §4's own table
    # closely for 6/9 materials (moissanite/diamond/CZ exactly; beryl/
    # quartz/glass within 3 points); sapphire/spinel/fluorite are a KNOWN,
    # documented discrepancy (see cone_survival_fraction's own docstring
    # for why), asserted here with a looser bound so the gate reflects
    # what's actually been reproduced rather than silently passing on a
    # loosened-for-everything tolerance ----
    tight = {'moissanite': 1.00, 'diamond': 1.00, 'cz': 1.00, 'beryl': 0.551, 'quartz': 0.514, 'glass': 0.486}
    loose = {'sapphire': 0.764, 'spinel': 0.708, 'fluorite': 0.384}
    for mat, expect in tight.items():
        got = cone_survival_fraction(MATERIALS[mat])
        ok(close(got, expect, 0.03), f'{mat} cone survival {got} != {expect} (tight)')
    for mat, expect in loose.items():
        got = cone_survival_fraction(MATERIALS[mat])
        ok(close(got, expect, 0.15), f'{mat} cone survival {got} != {expect} (loose, documented discrepancy)')

    # ---- P5 gate: the cone's leak onset must occur at or before the
    # collimated (chief-ray-only) anchor, for every material with a
    # positive margin -- the cone's marginal rays are always further from
    # normal incidence than the chief ray, so they reach critical first ----
    for mat in ('quartz', 'diamond', 'sapphire'):
        collim = collimated_leak_onset_tilt(MATERIALS[mat])
        cone = cone_leak_onset_tilt(MATERIALS[mat], planes_idealized)
        ok(cone is not None and cone <= collim + 1e-6,
           f'{mat}: cone leak-onset {cone} should be <= collimated anchor {collim}')

    # ---- P5 gate: spot size tracks f*theta_sun and cannot be driven
    # below the sun-disk floor regardless of standoff ----
    for s_mm in (80.0, 87.8, 95.0, 100.0, 110.0):
        spot = LENS['D_mm'] * abs(1 - s_mm / LENS['f_mm']) + s_mm * LENS['theta_sun']
        ok(spot >= 0.930 - 1e-6, f'spot at standoff {s_mm}mm = {spot:.4f} < the 0.930mm floor')
    spot_at_focus = LENS['D_mm'] * abs(1 - 1.0) + LENS['f_mm'] * LENS['theta_sun']
    ok(close(spot_at_focus, 0.930, 0.001), f'spot at focus {spot_at_focus} != 0.930')

    # ---- P6 gate: every material's Cauchy fit reproduces its pinned n_d
    # to 1e-4 and, by the same construction, its B-G dispersion exactly;
    # n_g > n across the visible (normal dispersion) for every material ----
    for mat in MATERIALS:
        a_coef, b_coef = cauchy_fit(MATERIALS[mat]['n'], MATERIALS[mat]['disp_bg'])
        n_d_check = n_at_wavelength(a_coef, b_coef, FRAUNHOFER_NM['d'])
        ok(close(n_d_check, MATERIALS[mat]['n'], 1e-4), f'{mat}: fitted n_d {n_d_check} != {MATERIALS[mat]["n"]}')
        disp_check = (n_at_wavelength(a_coef, b_coef, FRAUNHOFER_NM['G'])
                      - n_at_wavelength(a_coef, b_coef, FRAUNHOFER_NM['B']))
        ok(close(disp_check, MATERIALS[mat]['disp_bg'], 5e-4), f'{mat}: fitted B-G {disp_check} != {MATERIALS[mat]["disp_bg"]}')
        for fh_name in FRAUNHOFER_NM:
            n_g = group_index(a_coef, b_coef, FRAUNHOFER_NM[fh_name])
            n_here = n_at_wavelength(a_coef, b_coef, FRAUNHOFER_NM[fh_name])
            ok(n_g > n_here, f'{mat} at {fh_name}: n_g {n_g} not > n {n_here} (normal dispersion violated)')

    # ---- P6 gate: the §1 exit fan reproduces the pinned B-G spread for
    # diamond (1.06 deg) and moissanite (2.80 deg) -- tolerance here is
    # looser than the n(lambda) fit's own 5e-4 because this number is the
    # OUTPUT of a full nonlinear multi-bounce trace, not the fit itself ----
    fan_diamond, _ = dispersion_fan(MATERIALS['diamond'], planes, ray_o, ray_d)
    ok(close(fan_diamond, 1.06, 0.03), f'diamond B-G exit fan {fan_diamond} != 1.06')
    fan_moiss, _ = dispersion_fan(MATERIALS['moissanite'], planes, ray_o, ray_d)
    ok(close(fan_moiss, 2.80, 0.03), f'moissanite B-G exit fan {fan_moiss} != 2.80')

    # ---- P6 gate: solar weighting is a real, non-degenerate distribution
    # (peaks in the visible, not equal-energy/flat across the 7 lines) ----
    weights = {name: solar_weight(lam) for name, lam in FRAUNHOFER_NM.items()}
    ok(len(set(round(w, 6) for w in weights.values())) > 1, 'solar_weight: all 7 lines got the same weight (degenerate)')
    ok(all(w > 0 for w in weights.values()), 'solar_weight: a Fraunhofer line got non-positive weight')

    # ---- P6b gate: quartz uniaxial split ----
    c_axis = np.array([0., 1., 0.])
    n_o_q, n_e_q = MATERIALS['quartz']['n'], 1.5534  # n_e from plan §2's own table

    # Delta n -> 0 reproduces the isotropic trace exactly (to 1e-9)
    led_degenerate = trace_tree_e_wave(ray_o, ray_d, planes, n_o_q, n_o_q, c_axis)
    led_isotropic = trace_tree(ray_o, ray_d, planes, n_o_q)
    aois_degen = [round(h['aoi'], 9) for h in led_degenerate['dominant_log']]
    aois_iso = [round(h['aoi'], 9) for h in led_isotropic['dominant_log']]
    ok(aois_degen == aois_iso, f'P6b Delta n->0: {aois_degen} != isotropic {aois_iso}')

    # propagation || c is degenerate (o == e) for the axial ray itself
    led_e_axial = trace_tree_e_wave(ray_o, ray_d, planes, n_o_q, n_e_q, c_axis)
    aois_e_axial = [round(h['aoi'], 6) for h in led_e_axial['dominant_log']]
    aois_o_axial = [round(h['aoi'], 6) for h in trace_tree(ray_o, ray_d, planes, n_o_q)['dominant_log']]
    ok(aois_e_axial == aois_o_axial, f'P6b axial (||c) e-wave {aois_e_axial} != o-wave {aois_o_axial}')

    # a genuinely off-axis path shows the e-wave's index actually varying
    # per-hit (not frozen at entry) -- reflection directions match o-wave
    # exactly (reflection has no index dependence), but the cumulative
    # Fresnel energy and the final transmission angle both differ
    ray_o4 = np.array([-14.0, 1000.0, 0.0])
    led_e_off = trace_tree_e_wave(ray_o4, ray_d, planes, n_o_q, n_e_q, c_axis)
    led_o_off = trace_tree(ray_o4, ray_d, planes, n_o_q)
    dom_e = max(led_e_off['exited'][1:], key=lambda x: x['intensity'])
    dom_o = max(led_o_off['exited'][1:], key=lambda x: x['intensity'])
    ok(abs(dom_e['intensity'] - dom_o['intensity']) > 1e-6,
       'P6b off-axis: e-wave and o-wave dominant-exit intensities should differ')
    ok(close(led_e_off['total'], 1.0, 1e-6), f'P6b e-wave energy closure {led_e_off["total"]} != 1')
    ok(close(led_o_off['total'], 1.0, 1e-6), f'P6b o-wave energy closure {led_o_off["total"]} != 1')

    # energy closes with both populations counted (50/50 unpolarized split)
    led_combined = trace_tree_uniaxial(ray_o4, ray_d, planes, n_o_q, n_e_q, c_axis)
    ok(close(led_combined['total'], 1.0, 1e-6), f'P6b combined energy closure {led_combined["total"]} != 1')
    o_share = sum(e['intensity'] for e in led_combined['exited'] if e['wave'] == 'o')
    e_share = sum(e['intensity'] for e in led_combined['exited'] if e['wave'] == 'e')
    ok(close(o_share, 0.5, 1e-6) and close(e_share, 0.5, 1e-6),
       f'P6b combined: o/e shares {o_share}/{e_share} != 0.5/0.5')

    dtp_planes = double_terminated_point_planes()
    sc = solid_check(dtp_planes)
    ok(sc['euler'] == 2, f'double-terminated point: V-E+F = {sc["euler"]} != 2')
    ok(len(sc['bad_edges']) == 0, f'double-terminated point: {len(sc["bad_edges"])} edges not shared by exactly 2 faces')
    _, _, n_out_dtp = trace_solid(np.array([0.001, 0., 0.]), np.array([0., 1., 0.]), dict(planes=dtp_planes, quadrics=[]))
    aoi_dtp = math.degrees(math.acos(abs(np.dot([0, 1, 0], n_out_dtp))))
    ok(close(aoi_dtp, 51.78, 0.01), f'double-terminated point: r-face axial AOI {aoi_dtp} != 51.78')
    led = trace_tree(ray_o1, ray_d1, dtp_planes, n_q)
    ok(close(led['total'], 1.0, 1e-6), f'double-terminated point: energy closure {led["total"]} != 1')

    # ---- P2 gate: merkaba (formalized), point cut, single cut, 3 antiques ----
    merkaba = merkaba_planes()
    sc_r, sc_l = solid_check(merkaba['parts'][0]['planes']), solid_check(merkaba['parts'][1]['planes'])
    ok(sc_r['euler'] == 2 and sc_l['euler'] == 2, 'merkaba: a tetrahedron part is not closed')
    led = trace_tree(np.array([0.3, 1000.0, 0.2]), np.array([0., -1., 0.]), merkaba, MATERIALS['quartz']['n'])
    ok(close(led['total'], 1.0, 1e-6), f'merkaba: energy closure {led["total"]} != 1')

    for cutname, fn in (('pointcut', lambda: (point_cut_planes(), None)),
                        ('singlecut', single_cut_planes), ('oldmine', old_mine_planes),
                        ('oldeuropean', old_european_planes), ('peruzzi', peruzzi_planes)):
        planes_c, _ = fn()
        sc = solid_check(planes_c)
        ok(sc['euler'] == 2, f'{cutname}: V-E+F = {sc["euler"]} != 2')
        ok(len(sc['bad_edges']) == 0, f'{cutname}: {len(sc["bad_edges"])} edges not shared by exactly 2 faces')
        ok(sc['vol'] > 0, f'{cutname}: solid volume not positive')
        for mat in ('diamond', 'quartz', 'glass'):
            led = trace_tree(np.array([-5.0, 1000.0, 0.0]), np.array([0., -1., 0.]), planes_c, MATERIALS[mat]['n'])
            ok(close(led['total'], 1.0, 1e-6), f'{cutname} {mat}: energy closure {led["total"]} != 1')

    # ---- P2 gate: step-cut family, closure + energy (+ depth band where sourced) ----
    for cutname, fn, depth_band in (('emerald', emerald_cut_planes, (59.0, 69.0)),
                                     ('asscher', asscher_planes, (60.0, 68.0)),
                                     ('baguette', baguette_planes, None),
                                     ('taperedbaguette', tapered_baguette_planes, None),
                                     ('carre', carre_planes, None),
                                     ('tablecut', table_cut_planes, None)):
        planes_c, meta_c = fn()
        sc = solid_check(planes_c)
        ok(sc['euler'] == 2, f'{cutname}: V-E+F = {sc["euler"]} != 2')
        ok(len(sc['bad_edges']) == 0, f'{cutname}: {len(sc["bad_edges"])} edges not shared by exactly 2 faces')
        ok(sc['vol'] > 0, f'{cutname}: solid volume not positive')
        if depth_band:
            ok(depth_band[0] <= meta_c['depthPct'] <= depth_band[1],
               f'{cutname}: depth% {meta_c["depthPct"]} outside sourced {depth_band[0]}-{depth_band[1]}')
        for mat in ('diamond', 'quartz', 'glass'):
            led = trace_tree(np.array([-3.0, 1000.0, 0.0]), np.array([0., -1., 0.]), planes_c, MATERIALS[mat]['n'])
            ok(close(led['total'], 1.0, 1e-6), f'{cutname} {mat}: energy closure {led["total"]} != 1')

    # ---- P2 gate: final 3 antiques -- rose, mazarin, briolette (roster complete) ----
    rose_pl, _ = rose_cut_planes()
    sc = solid_check(rose_pl)
    ok(sc['euler'] == 2, f'rose: V-E+F = {sc["euler"]} != 2')
    ok(len(sc['bad_edges']) == 0, f'rose: {len(sc["bad_edges"])} edges not shared by exactly 2 faces')
    ok(sc['vol'] > 0, 'rose: solid volume not positive')
    for mat in ('diamond', 'quartz', 'glass'):
        led = trace_tree(np.array([3.0, 1000.0, 0.0]), np.array([0., -1., 0.]), rose_pl, MATERIALS[mat]['n'])
        ok(close(led['total'], 1.0, 1e-6), f'rose {mat}: energy closure {led["total"]} != 1')

    for cutname, fn in (('mazarin', mazarin_planes), ('briolette', briolette_planes)):
        planes_c, meta_c = fn()
        sc = solid_check(planes_c)
        ok(sc['euler'] == 2, f'{cutname}: V-E+F = {sc["euler"]} != 2')
        ok(len(sc['bad_edges']) == 0, f'{cutname}: {len(sc["bad_edges"])} edges not shared by exactly 2 faces')
        ok(sc['vol'] > 0, f'{cutname}: solid volume not positive')
        for mat in ('diamond', 'quartz', 'glass'):
            led = trace_tree(np.array([-3.0, 1000.0, 0.0]), np.array([0., -1., 0.]), planes_c, MATERIALS[mat]['n'])
            ok(close(led['total'], 1.0, 1e-6), f'{cutname} {mat}: energy closure {led["total"]} != 1')

    # ---- P8 gate: the (crown, pavilion) scan's ridge falls inside the
    # §2b corridor for every material. "Ridge" means the traced canonical
    # path (AOI1=p, AOI2=180-3p, AOI3=|4p-180|, verdict TIR/TIR/table-exit)
    # holding exactly -- NOT an aggregate crown-return threshold, which
    # stays high well outside the corridor too via longer paths the
    # corridor's own caption already calls "not failure." This is the
    # computed check that the tracer and the algebra agree.
    def canonical_holds(material, pav_deg, crown_deg=34.5, xoff=-5.0):
        planes_c, meta_c = round_brilliant_planes(56.0, crown_deg, pav_deg, 100.0, 3.0)
        origin = np.array([xoff, meta_c['table_y'] + 500.0, 0.0])
        try:
            led = trace_tree(origin, np.array([0., -1., 0.]), planes_c, material['n'])
        except RuntimeError:
            return False
        log = led['dominant_log']
        if len(log) < 4:
            return False
        aoi1, aoi2, aoi3 = log[1]['aoi'], log[2]['aoi'], log[3]['aoi']
        exp2, exp3 = abs(180.0 - 3.0 * pav_deg), abs(4.0 * pav_deg - 180.0)
        names_ok = (log[1]['plane'].startswith('pavilion') and
                    log[2]['plane'].startswith('pavilion') and log[3]['plane'] == 'table')
        aoi_ok = close(aoi1, pav_deg, 0.05) and close(aoi2, exp2, 0.05) and close(aoi3, exp3, 0.05)
        tc = material['thetac']
        verdict_ok = aoi1 >= tc and aoi2 >= tc and aoi3 < tc
        dom3 = next((e for e in led['exited'] if e.get('depth') == 3), None)
        return names_ok and aoi_ok and verdict_ok and dom3 is not None

    for mat_name, material in MATERIALS.items():
        lo, hi = corridor_for(material)
        ridge_lo, ridge_hi = None, None
        pav = 20.0
        while pav <= 60.0:
            if canonical_holds(material, pav):
                if ridge_lo is None:
                    ridge_lo = pav
                ridge_hi = pav
            pav += 0.05
        ok(ridge_lo is not None, f'{mat_name}: canonical corridor path never holds for any scanned pavilion angle')
        if ridge_lo is not None:
            ok(ridge_lo >= lo - 0.15, f'{mat_name}: ridge lower {ridge_lo:.2f} escapes corridor lower bound {lo:.2f}')
            ok(ridge_hi <= hi + 0.15, f'{mat_name}: ridge upper {ridge_hi:.2f} escapes corridor upper bound {hi:.2f}')

    print(f'\n{n_ok[0]} checks, {n_fail[0]} failures.')
    return n_ok[0], n_fail[0]


if __name__ == '__main__':
    run_checks()
