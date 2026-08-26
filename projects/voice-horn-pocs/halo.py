#!/usr/bin/env python3
"""POC-C HALO -- printable STLs for the cycle-11 chord.

BUILD NUMBERS, from pareto_c11.csv (the 6.16-cent member the doc calls for):
    V   = 626.06 / 207.25 / 188.79 / 112.45 cm3
    a_n =   5.51 /  10.94 /  14.88 /  14.42 mm   (necks ascend with pitch:
    L_n = 34.84 mm, slot 4 mm                     constant t60 wants it, since
    chamber = 4.021 L (the model's dia 160 x 200) Q ~ a_n sqrt(f))

Two build decisions the network does not see, and one it does:

* The chamber's cavity VOLUME is the acoustic quantity, not its shape, so the
  lid is replaced by a 45-degree cone: the chamber then prints in one piece,
  upright, with no lid, no support and no bridge over a 168 mm span. Radius and
  cylinder height are solved so cylinder + cone - ribs = 4.021 L exactly.
* The rib array (8 -> 12 mm, RULE S1 -- two orders below lambda/8 for the chord,
  fully visible to breath noise above 3.6 kHz) eats 0.32 L of that cavity, which
  is why the built radius is 84 mm rather than the model's 80.
* The wind cowls become flush GRILLES across each neck mouth and the inlet. A
  mesh dome could not be printed inside a sealed chamber; a grille of 1 mm bars
  on 7 mm centres is 86% open -- transparent to a 220-550 Hz chord, opaque to
  the DC flow that would otherwise make an ocarina of a tuned neck.

Every sphere is split on the vertical plane through its neck so both halves
print dome-up off a flat face, with the half-neck as a groove in that face.

Run from this directory:  python halo.py [--no-skin]
"""
import argparse
import math

import numpy as np
import trimesh

import horn_lib as H

OUT = "stl/"

V_TGT = np.array([626.0571, 207.2529, 188.7921, 112.4453]) * 1000.0   # mm^3
A_NECK = np.array([5.5051, 10.9427, 14.8751, 14.4236])
L_N, SLOT_G = 34.8379, 4.0
FT = np.array([220.0, 330.0, 440.0, 550.0])
V_CH = math.pi * 80.0 ** 2 * 200.0                    # the model's chamber, mm^3
D_CH_MODEL = 160.0

R_CH = 84.0                     # built chamber radius (ribs take the difference)
N_RIB, W_RIB, H_RIB = 48, 4.0, (8.0, 12.0)
WALL = H.WALL
Z_INLET, Z_NECK, Z_SLOT = 35.0, 85.0, 145.0
SLOT_POSTS, POST_W = 4, 10.0
INLET_R, INLET_LEN, CUP_R = 20.0, 70.0, 22.5
TH_PITCH, TH_DEPTH, CLR = 4.0, 1.2, 0.50   # printed-thread fit, PETG
BOSS_LEN, SOCK_WALL = 22.0, 3.0
GRILLE_BAR, GRILLE_GAP = 1.0, 6.0
RHO, C_S, MU = 1.2, 343.0, 1.81e-5


def rib_h(th):
    """Rib height, swept 8 -> 12 mm once round the chamber (a printed rainbow)."""
    return H_RIB[0] + (H_RIB[1] - H_RIB[0]) * (np.asarray(th) % (2 * math.pi)) / (2 * math.pi)


def chamber_geometry():
    """Solve the cylinder height that makes cavity volume = the model's 4.021 L."""
    th_c = np.arange(N_RIB) * 2 * math.pi / N_RIB
    a_rib = float(np.sum(W_RIB * rib_h(th_c)))        # rib cross-section, mm^2
    cone = math.pi * R_CH ** 3 / 3.0                  # 45-degree cap, apex on axis
    h_cyl = (V_CH - cone) / (math.pi * R_CH ** 2 - a_rib)
    return h_cyl, a_rib, cone


def rib_theta():
    """Non-uniform theta: tight samples on each rib edge, a few in each gap."""
    span = 2 * math.pi / N_RIB
    half = 0.5 * W_RIB / R_CH
    out = []
    for k in range(N_RIB):
        c = k * span
        out += [c - half - 1e-4, c - half, c + half, c + half + 1e-4]
        out += list(c + half + 1e-4 + (span - 2 * half - 2e-4) * np.linspace(0.1, 0.9, 4))
    return np.asarray(out) % (2 * math.pi)


def chamber(skin=True):
    h_cyl, a_rib, cone = chamber_geometry()
    th = rib_theta()
    span = 2 * math.pi / N_RIB
    half = 0.5 * W_RIB / R_CH
    k = np.floor((th + half) / span).astype(int) % N_RIB
    dth = (th - k * span + math.pi) % (2 * math.pi) - math.pi
    on_rib = np.abs(dth) <= half + 1e-9
    rib = np.where(on_rib, rib_h(k * span), 0.0)

    z_top = h_cyl + R_CH - 8.0                 # truncate the spire and cap it:
    z = np.unique(np.concatenate([              # a cone run to a point leaves a
        np.linspace(0, h_cyl, 40), [h_cyl - 0.5],   # pinhole in a sealed chamber
        np.linspace(h_cyl, z_top, 40)]))
    inner = np.where(z <= h_cyl, R_CH,
                     np.maximum(R_CH - (z - h_cyl), 0.2))          # 45-deg cone
    ri = inner[:, None] - np.where(z[:, None] <= h_cyl, rib[None, :], 0.0)
    ro = inner[:, None] + WALL * np.where(z[:, None] <= h_cyl, 1.0, math.sqrt(2.0))
    body = H.shell(z, ri, ro, theta=th)
    floor = H.cyl((0, 0, -WALL), (0, 0, 0.5), R_CH + WALL, 192)
    finial = H.cyl((0, 0, z_top - 1.0), (0, 0, z_top + 5.0), 8.0 + WALL, 64)
    body = H.union(body, floor, finial)

    body = H.difference(body, slot_cutter())
    for i in range(4):
        body = H.union(body, boss(i, 2 * math.pi * i / 4))
    body = H.union(body, inlet())
    for i in range(4):
        body = H.difference(body, bore_cut(2 * math.pi * i / 4, A_NECK[i],
                                           Z_NECK, R_CH + L_N + 5))
    body = H.difference(body, bore_cut(math.pi / 4, INLET_R, Z_INLET,
                                       R_CH + INLET_LEN + 30, cup=True))
    for i in range(4):
        body = H.union(body, grille(2 * math.pi * i / 4, A_NECK[i], Z_NECK))
    body = H.union(body, grille(math.pi / 4, INLET_R, Z_INLET))
    if skin:
        body = H.union(body, *chamber_skin(h_cyl))
    return body


def slot_cutter():
    """Annular radiating slot, sized so the OPEN area matches the model's
    pi * 160 * 4 mm2 even after four posts hold the crown on."""
    r_mid = R_CH + 0.5 * WALL
    open_len = 2 * math.pi * r_mid - SLOT_POSTS * POST_W
    g = math.pi * D_CH_MODEL * SLOT_G / open_len
    ring = H.cyl((0, 0, Z_SLOT - 0.5 * g), (0, 0, Z_SLOT + 0.5 * g), R_CH + WALL + 6, 256)
    ring = H.difference(ring, H.cyl((0, 0, Z_SLOT - g), (0, 0, Z_SLOT + g),
                                    R_CH - 14.0, 256))
    posts = []
    for k in range(SLOT_POSTS):
        a = 2 * math.pi * k / SLOT_POSTS + math.pi / SLOT_POSTS
        b = trimesh.creation.box((POST_W, 2 * (WALL + 8), 2 * g))
        T = trimesh.transformations.rotation_matrix(a, (0, 0, 1))
        T[:3, 3] = [(R_CH + WALL) * math.cos(a), (R_CH + WALL) * math.sin(a), Z_SLOT]
        b.apply_transform(T)
        posts.append(b)
    return H.difference(ring, *posts)


def _radial(ang, z, mesh):
    T = trimesh.transformations.rotation_matrix(math.pi / 2, (0, 1, 0))
    mesh.apply_transform(T)                                  # +z -> +x
    T2 = trimesh.transformations.rotation_matrix(ang, (0, 0, 1))
    T2[:3, 3] = [0, 0, z]
    mesh.apply_transform(T2)
    return mesh


def boss(i, ang):
    """Neck boss: external thread for the sphere, teardrop roof so a horizontal
    dia-40 stub does not need support."""
    a = A_NECK[i]
    z = np.linspace(R_CH - 2.0, R_CH + BOSS_LEN, 40)
    th = H.grid_theta(96)
    ri = np.full((len(z), len(th)), a)
    ro = ri + WALL + H.thread(z, th, R_CH + 4.0, R_CH + BOSS_LEN - 2.0,
                              TH_PITCH, TH_DEPTH)
    m = H.shell(z, ri, ro)
    roof = trimesh.creation.extrude_polygon(
        _teardrop(a + WALL + TH_DEPTH), BOSS_LEN + 2.0)
    roof.apply_transform(trimesh.transformations.translation_matrix(
        (0, 0, R_CH - 2.0)))
    return _radial(ang, Z_NECK, H.union(m, roof))


def _teardrop(r):
    from shapely.geometry import Polygon
    return Polygon([(-r, 0), (r, 0), (0, r)])


def inlet():
    z = np.linspace(R_CH - 2.0, R_CH + INLET_LEN, 30)
    th = H.grid_theta(96)
    ri = np.full((len(z), len(th)), INLET_R)
    lip = np.clip((z - (R_CH + INLET_LEN - 12.0)) / 12.0, 0, 1)
    ri = ri + (CUP_R - INLET_R) * lip[:, None]
    ro = ri + WALL
    m = H.shell(z, ri, ro)
    roof = trimesh.creation.extrude_polygon(_teardrop(INLET_R + WALL), INLET_LEN)
    roof.apply_transform(trimesh.transformations.translation_matrix(
        (0, 0, R_CH - 2.0)))
    return _radial(math.pi / 4, Z_INLET, H.union(m, roof))


def bore_cut(ang, r, z, reach, cup=False):
    m = H.cyl((0, 0, R_CH - 30.0), (0, 0, reach), r, 96)
    if cup:
        m = H.union(m, H.cyl((0, 0, reach - 12.0), (0, 0, reach + 2.0), CUP_R, 96))
    return _radial(ang, z, m)


def grille(ang, r, z):
    """Flush bars across a mouth: 86% open to the chord, opaque to DC flow."""
    bars = []
    for x in np.arange(-r + GRILLE_GAP * 0.5, r - 0.5, GRILLE_GAP):
        y = math.sqrt(max(r * r - x * x, 0.0))
        if y < 1.0:
            continue
        b = trimesh.creation.box((GRILLE_BAR, 2 * y, 3.0))
        b.apply_transform(trimesh.transformations.translation_matrix(
            (x, 0, R_CH + 1.5)))
        bars.append(b)
    return _radial(ang, z, H.union(*bars))


def chamber_skin(h_cyl):
    """Exterior gradient: polished equator band where the crown is carried,
    sharp claw density rising toward the cone (RULE S3) -- the crown of the
    crystal, curling upward through the shard field."""
    zp = np.concatenate([np.linspace(0, h_cyl, 60),
                         np.linspace(h_cyl, h_cyl + R_CH - WALL, 60)])
    rp = np.where(zp <= h_cyl, R_CH + WALL,
                  np.maximum(R_CH - (zp - h_cyl), 0.2) + WALL * math.sqrt(2))
    keep = ((zp < Z_INLET - 30) | (zp > Z_INLET + 30))
    P, N, frac, area = H.phyllotaxis(zp[keep], rp[keep], 380,
                                     weight=lambda zz: np.where(zz < h_cyl * 0.75,
                                                                0.25, 1.0),
                                     z_lo=Z_SLOT + 12.0)
    h = 1.2 + 7.0 * frac
    rb = np.clip(h / 7.0, 0.7, 0.42 * math.sqrt(area))
    return H.claws(P, N + np.array([0.0, 0.0, 0.8]), h, rb, r_tip=0.18,
                   curl=np.array([0.0, 0.0, 1.0]), curl_amt=0.6)


# ------------------------------------------------------------------ spheres

def sphere_radius(v_target, subdiv=5):
    """Scale so the FACETED cavity holds the modelled volume, not the ideal one."""
    unit = trimesh.creation.icosphere(subdivisions=subdiv, radius=1.0).volume
    return (v_target / unit) ** (1.0 / 3.0)


def sphere_parts(i, skin=True):
    a, r_in = A_NECK[i], sphere_radius(V_TGT[i])
    ctr = R_CH + L_N + r_in                       # sphere centre, from the axis
    outer = trimesh.creation.icosphere(subdivisions=5, radius=r_in + WALL)
    inner = trimesh.creation.icosphere(subdivisions=5, radius=r_in)
    outer.apply_translation((ctr, 0, 0))
    inner.apply_translation((ctr, 0, 0))

    z = np.linspace(R_CH + 4.0, R_CH + L_N + 6.0, 40)
    th = H.grid_theta(96)
    base = np.full((len(z), len(th)), a + WALL + CLR)
    rid = H.thread(z, th, R_CH + 4.0, R_CH + BOSS_LEN - 2.0, TH_PITCH, TH_DEPTH)
    inside = (z <= R_CH + BOSS_LEN)[:, None]
    ri = np.where(inside, base + rid, a)
    ro = np.full_like(ri, a + WALL + CLR + TH_DEPTH + SOCK_WALL)
    stub = H.shell(z, ri, ro)
    stub = _radial(0.0, 0.0, stub)

    body = H.difference(H.union(outer, stub), inner,
                        _radial(0.0, 0.0, H.cyl((0, 0, R_CH), (0, 0, ctr), a, 64)))
    body = H.difference(body, plug_bore(ctr, r_in, i))
    body = H.union(body, plug_boss(ctr, r_in, i))
    if skin:
        body = H.union(body, *sphere_skin(ctr, r_in, i))
        body = H.difference(body, *pinch_saddles(ctr, r_in))

    tabs, holes = clamp_tabs(ctr, r_in)
    body = H.difference(H.union(body, tabs), holes)
    keep = trimesh.creation.box((600, 300, 600))
    keep.apply_translation((ctr, 150, 0))
    keep2 = trimesh.creation.box((600, 300, 600))
    keep2.apply_translation((ctr, -150, 0))
    return H.intersection(body, keep), H.intersection(body, keep2)


def plug_dims(i, r_in):
    d_p = 0.5 * (2 * r_in)
    travel = 0.2 * V_TGT[i] / (math.pi * (0.5 * d_p) ** 2)
    return 0.5 * d_p, travel


def plug_bore(ctr, r_in, i):
    """Bore for the tuning plug, on +y: perpendicular to the split plane, so it
    lands whole in one half instead of being sawn down the middle."""
    rp, travel = plug_dims(i, r_in)
    return H.cyl((ctr, r_in - 3.0, 0),
                 (ctr, r_in + WALL + 0.5 * travel + 14.0, 0), rp, 96)


def plug_boss(ctr, r_in, i):
    rp, travel = plug_dims(i, r_in)
    z0, z1 = r_in - 2.0, r_in + WALL + 0.5 * travel + 12.0
    z = np.linspace(z0, z1, 40)
    th = H.grid_theta(96)
    ri = np.full((len(z), len(th)), rp) + H.thread(z, th, r_in + WALL, z1 - 2.0,
                                                   TH_PITCH, TH_DEPTH)
    ro = np.full_like(ri, rp + TH_DEPTH + SOCK_WALL + 1.0)
    m = H.shell(z, ri, ro)                                   # built along +z
    m.apply_transform(trimesh.transformations.rotation_matrix(-math.pi / 2,
                                                              (1, 0, 0)))
    m.apply_transform(trimesh.transformations.translation_matrix((ctr, 0, 0)))
    return m


def plug(i, skin=False):
    """Tuning plug: +-10% of cavity volume, i.e. the chord centres anywhere from
    about G3 to B3 by turning four of these (f ~ V^-1/2, Law 14.6)."""
    r_in = sphere_radius(V_TGT[i])
    rp, travel = plug_dims(i, r_in)
    L = travel + 24.0
    z = np.linspace(0, L, 60)
    th = H.grid_theta(96)
    ro = np.full((len(z), len(th)), rp - CLR) + H.thread(z, th, 6.0, L - 8.0,
                                                         TH_PITCH, TH_DEPTH)
    ri = np.full_like(ro, rp - CLR - 3.0)
    ri = np.where((z >= L - 6.0)[:, None], 0.2, ri)          # closed face
    cup = H.shell(z, ri, ro)
    head = H.cyl((0, 0, -8.0), (0, 0, 0.5), rp + 7.0, 6)     # hex grip
    m = H.union(cup, head)
    if not skin:
        return m
    dishes = []                                  # fingertip scallop per hex flat
    r_flat = (rp + 7.0) * math.cos(math.pi / 6.0)
    for k in range(6):
        a = (2 * k + 1) * math.pi / 6.0          # trimesh hex verts sit at k*60,
        dish = trimesh.creation.icosphere(subdivisions=3, radius=9.0)   # flats
        dish.apply_translation(((r_flat - 1.2 + 9.0) * math.cos(a),     # between
                                (r_flat - 1.2 + 9.0) * math.sin(a), -3.75))
        dishes.append(dish)
    return H.difference(m, *dishes)


def clamp_tabs(ctr, r_in):
    """Four M3 clamping ears straddling the split plane -- the halves are bolted
    and sealed rather than glued, since a Helmholtz Q lives on a tight cavity."""
    tabs, holes = [], []
    for k in range(4):
        a = math.pi / 4 + k * math.pi / 2
        c = np.array([ctr + (r_in + WALL + 5.0) * math.cos(a), 0.0,
                      (r_in + WALL + 5.0) * math.sin(a)])
        b = trimesh.creation.box((16.0, 10.0, 16.0))
        T = trimesh.transformations.rotation_matrix(a, (0, 1, 0))
        T[:3, 3] = c
        b.apply_transform(T)
        tabs.append(b)
        holes.append(H.cyl(c + np.array([0, -8, 0]), c + np.array([0, 8, 0]), 1.7, 20))
    return H.union(*tabs), trimesh.util.concatenate(holes)


def sphere_skin(ctr, r_in, i):
    """Sharp claw density rising toward the outer pole -- painful to seize by a
    sphere (RULE S3, literally now), comfortable only at the centre and at the
    two pinch saddles the tuning grip actually uses."""
    n = 130 + 40 * i
    idx = np.arange(n)
    zc = 1 - 2 * (idx + 0.5) / n
    rr = np.sqrt(np.maximum(1 - zc * zc, 0))
    ph = idx * H.GOLDEN
    d = np.stack([rr * np.cos(ph), rr * np.sin(ph), zc], -1)
    rp, _ = plug_dims(i, r_in)
    r_boss = rp + TH_DEPTH + SOCK_WALL + 1.0
    off_axis = r_in * np.hypot(d[:, 0], d[:, 2])
    keep = ((d[:, 0] > -0.25)                                # clear of the neck
            & ~((d[:, 1] > 0) & (off_axis < r_boss + 4.0))   # of the plug boss
            & (np.abs(d[:, 2]) < 0.68))                      # of the saddles
    d = d[keep]
    P = np.array([ctr, 0, 0]) + d * (r_in + WALL - 1.0)
    h = 1.5 + 5.0 * np.clip((d[:, 0] + 0.2) / 1.2, 0, 1)
    rb = np.clip(h / 7.0, 0.7, 2.0)
    pole = np.array([1.0, 0.0, 0.0])                         # curl toward the
    tang = pole[None, :] - d * d[:, :1]                      # outer pole
    tang /= np.maximum(np.linalg.norm(tang, axis=1, keepdims=True), 1e-9)
    return H.claws(P, d, h, rb, r_tip=0.18, curl=tang, curl_amt=0.6, sink=2.0)


def pinch_saddles(ctr, r_in):
    """Two smooth 2 mm dished pits at the +-z poles -- the thumb-and-finger
    pinch that turns a sphere against its thread. They straddle the split seam
    (each half carries half a dish), sit clear of the clamp ears at 45 degrees,
    and are the one polite place on an otherwise armored orb."""
    out = []
    for s in (+1.0, -1.0):
        dish = trimesh.creation.icosphere(subdivisions=4, radius=60.0)
        dish.apply_translation((ctr, 0.0, s * (r_in + WALL - 2.0 + 60.0)))
        out.append(dish)
    return out


# -------------------------------------------------------------- verification

def network(V, a_i, L_n, g_s, frq=np.linspace(150, 700, 1101)):
    """The cycle-11 lumped network, re-run on AS-BUILT numbers.

    Chamber compliance, radiating slot, four shunted Helmholtz branches, glottal
    Norton source -- identical to evolve_halo_c11.py, only the inputs change."""
    om = 2 * math.pi * frq
    Vm, an = np.asarray(V) * 1e-9, np.asarray(a_i) * 1e-3
    A_s = math.pi * 0.160 * g_s * 1e-3
    L_s = 0.005 + 1.7 * math.sqrt(A_s / math.pi) * 0.5
    Y = 1j * om * (V_CH * 1e-9) / (RHO * C_S ** 2)
    Zs = 1j * om * RHO * L_s / A_s + RHO * C_S * (om / C_S) ** 2 / (2 * math.pi)
    Y = Y + 1 / Zs
    Zres = []
    for Vi, a in zip(Vm, an):
        A_n = math.pi * a ** 2
        Le = L_n * 1e-3 + 1.7 * a
        Zi = (2 * Le * np.sqrt(2 * MU * RHO * om) / (A_n * a)
              + 1j * om * RHO * Le / A_n + 1 / (1j * om * Vi / (RHO * C_S ** 2)))
        Zres.append(Zi)
        Y = Y + 1 / Zi
    P = 1.0 / (Y + 1 / 4.0e6)
    out = []
    for i, ft in enumerate(FT):
        U = np.abs(P / Zres[i])
        win = (frq > ft * 0.7) & (frq < ft * 1.35)
        j = int(np.argmax(np.where(win, U, 0)))
        half = U[j] / math.sqrt(2)
        lo, hi = j, j
        while lo > 0 and U[lo] > half:
            lo -= 1
        while hi < len(frq) - 1 and U[hi] > half:
            hi += 1
        bw = max(frq[hi] - frq[lo], 0.5)
        out.append((frq[j], 1200 * math.log2(frq[j] / ft), 2.2 * (frq[j] / bw) / frq[j]))
    return out


def verify():
    h_cyl, a_rib, cone = chamber_geometry()
    V_built = math.pi * R_CH ** 2 * h_cyl + cone - a_rib * h_cyl
    r_mid = R_CH + 0.5 * WALL
    g_open = math.pi * D_CH_MODEL * SLOT_G / (2 * math.pi * r_mid - SLOT_POSTS * POST_W)
    print(f"HALO chamber: r {R_CH:.1f} + cylinder {h_cyl:.1f} + 45-deg cone "
          f"{R_CH:.0f} -> cavity {V_built/1e6:.3f} L (model {V_CH/1e6:.3f} L, "
          f"ribs take {a_rib*h_cyl/1e6:.3f} L)")
    print(f"  slot gap built {g_open:.2f} mm over {SLOT_POSTS} posts = the model's "
          f"{math.pi*D_CH_MODEL*SLOT_G:.0f} mm2 of open area")
    V = [float(trimesh.creation.icosphere(subdivisions=5,
                                          radius=sphere_radius(v)).volume)
         for v in V_TGT]
    rows = network(V, A_NECK, L_N, SLOT_G)
    tot = sum(abs(c) for _, c, _ in rows)
    print("  as-built chord:  " + "  ".join(
        f"{f:.1f}Hz({c:+.1f}c, t60 {t:.2f}s)" for f, c, t in rows))
    print(f"  total tuning error {tot:.1f} cents "
          f"(published cycle-11 figure: 6.2)")
    for i, v in enumerate(V):
        r = sphere_radius(V_TGT[i])
        rp, travel = plug_dims(i, r)
        print(f"  sphere {i+1}: dia {2*r:.1f} mm, cavity {v/1000:.1f} cm3, "
              f"neck a {A_NECK[i]:.2f}, plug dia {2*rp:.0f} x {travel:.0f} mm travel")


def build(skin=True):
    verify()
    rows = [H.report("halo_chamber", chamber(skin=skin), OUT + "halo_chamber.stl")]
    for i in range(4):
        a, b = sphere_parts(i, skin=skin)
        rows.append(H.report(f"halo_sphere{i+1}_a", a, OUT + f"halo_sphere{i+1}_a.stl"))
        rows.append(H.report(f"halo_sphere{i+1}_b", b, OUT + f"halo_sphere{i+1}_b.stl"))
        rows.append(H.report(f"halo_plug{i+1}", plug(i, skin=skin),
                             OUT + f"halo_plug{i+1}.stl"))
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--skin", action="store_true",
                    help="add the ornament pass; default is plain smooth bodies")
    a = ap.parse_args()
    rows = build(skin=a.skin)
    print()
    H.print_table(rows)
