#!/usr/bin/env python3
"""Quartz Flash Windscrew -- the lifter-drum tumbler.

Run from this directory:  python tumbler.py
Prints the lifter tables and writes tumbler.png.  windscrew.py takes the drum
dimensions from here; physics.py takes the drum torque and impact numbers.

A round drum with interchangeable straight lifter ribs and replaceable quartz
liner trays, turned at a governed fraction of critical speed.  Models
(assumptions in README.md):
  1. Critical speed of the rock centres, omega_c = sqrt(g / r_c), with
     r_c = liner surface radius - rock radius.
  2. Release: quasi-static balance on a rock held in the rib/liner corner
     (gravity, centrifugal, Coulomb friction on the rib face).  A rib lower than
     the rock radius acts through its edge, i.e. as a face leaning back by
     gamma = atan((a - h) / sqrt(a^2 - (a - h)^2)).
  3. Slide along the face to the tip, in the drum frame, Coriolis included.
  4. Free flight in the inertial frame until the rock meets the charge bed, the
     liner.  Both are quartz, so both give a flash.  There is no through-axle:
     the drum hangs on two 608s in a hub at its inboard end, because an axle
     intercepts the flights, which cross close to the centre at 60-80 % of critical.
  5. Charge bed: circular segment holding the bulk volume of the stones, tilted
     toward the rising side by a dynamic repose angle.
  6. Drum torque = lifts per rev x energy per lift / 2 pi.  It is nearly flat
     from 50 to 95 % of critical and collapses at critical, so a torque-limited
     clutch cannot hold the speed: the build uses a centrifugal governor.
"""
import math
import numpy as np

G = 9.81

# drum (windscrew.py builds from these)
R_IN = 0.060                # shell inner radius, m
DRUM_L = 0.088              # inner length, m
LINER_T = 0.007             # liner tray (2 mm) + epoxied quartz chips (~5 mm)
R_W = R_IN - LINER_T        # working surface: the quartz liner
N_RAILS = 12                # dovetail rails: ribs on alternate rails (6) or all (12); trays between
N_RIBS = 6
AXLE_R = 0.0                # no through-axle: the shaft ends in the inboard hub boss (see README)
RIBS = {                    # name: (height above the liner m, face lean-back from radial deg)
    "h08_f90": (0.008, 0.0),
    "h14_f90": (0.014, 0.0),
    "h14_f60": (0.014, 30.0),
    "h14_f45": (0.014, 45.0),
}
DEFAULT_RIB = "h14_f60"

# charge
ROCK_A = 0.010              # rock radius (20 mm stones)
M_ROCK = 0.012              # kg
V_ROCK = M_ROCK / 2650.0
N_STONES = 10
PACKING = 0.60              # solid fraction of a loose bed of stones
PSI_BED = math.radians(35)  # bed tilt toward the rising side
CARRY_MAX = 3               # stones one 88 mm rib can hold
MU_RIB = 0.35               # quartz on PETG
PHI_TARGET = 0.75           # governed fraction of critical speed

QUARTZ = ("liner", "bed")

def omega_c(a=ROCK_A):
    return math.sqrt(G / (R_W - a))

def rpm(w):
    return w * 30 / math.pi

def bed(n_stones=N_STONES):
    """(fill fraction of the cross-section, half-angle alpha, chord distance d, axis unit vector)."""
    f = n_stones * V_ROCK / PACKING / (math.pi * R_W**2 * DRUM_L)
    lo, hi = 0.0, math.pi
    for _ in range(60):
        mid = (lo + hi) / 2
        if (mid - math.sin(mid) * math.cos(mid)) / math.pi < f:
            lo = mid
        else:
            hi = mid
    alpha = (lo + hi) / 2
    e = np.array([math.sin(PSI_BED), -math.cos(PSI_BED)])
    return f, alpha, R_W * math.cos(alpha), e

def effective_lean(h, beta):
    if h >= ROCK_A:
        return beta
    gamma = math.atan((ROCK_A - h) / math.sqrt(ROCK_A**2 - (ROCK_A - h)**2))
    return max(beta, gamma)

def lift(rib, phi, n_stones=N_STONES, keep_path=False):
    """One rock carried by one rib at speed phi * omega_c (drum turning CCW).
    Returns a dict, or None when the rock is never released (centrifuging)."""
    h, beta_deg = RIBS[rib]
    beta = effective_lean(h, math.radians(beta_deg))
    r_c = R_W - ROCK_A
    w = phi * omega_c()
    f_vec = np.array([-math.sin(beta), math.cos(beta)])     # drum frame, rib at the bottom
    n_vec = np.array([math.cos(beta), math.sin(beta)])
    P0 = np.array([0.0, -r_c])
    _, alpha, d_bed, e_bed = bed(n_stones)

    def applied(th, P, sdot):
        grav = np.array([-G * math.sin(th), -G * math.cos(th)])
        return grav + w * w * P + 2 * w * sdot * n_vec       # + Coriolis

    # 1. release scan, starting where the rib leaves the bed
    th_r = None
    for th in np.arange(max(PSI_BED + alpha, 0.0), math.pi, math.radians(0.25)):
        A = applied(th, P0, 0.0)
        a_f, N = A @ f_vec, -(A @ n_vec)
        if N <= 0 or a_f > MU_RIB * N:
            th_r = th
            break
    if th_r is None:
        return None
    # 2. slide to the tip
    s, sdot, th, dt = 0.0, 0.0, th_r, 1e-4
    s_tip = max(h - ROCK_A, 0.0)
    while s < s_tip:
        P = P0 + s * f_vec
        A = applied(th, P, sdot)
        a_f, N = A @ f_vec, -(A @ n_vec)
        if N <= 0:
            break
        sdot = max(sdot + (a_f - MU_RIB * N) * dt, 0.0)
        s += sdot * dt
        th += w * dt
    P = P0 + s * f_vec
    R = np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])
    p0 = R @ P
    v0 = w * np.array([-p0[1], p0[0]]) + R @ (sdot * f_vec)
    # 3. flight (vectorised), first contact wins
    t = np.arange(2e-4, 1.5, 2e-4)
    pts = p0 + np.outer(t, v0) + 0.5 * np.outer(t * t, [0.0, -G])
    rad = np.linalg.norm(pts, axis=1)
    hits = {"axle": rad <= AXLE_R + ROCK_A if AXLE_R > 0 else np.zeros(len(t), bool),
            "bed": (pts @ e_bed >= d_bed - ROCK_A) & (t > 0.01),
            "liner": (rad >= r_c) & (t > 0.01)}
    first = {k: (np.argmax(v) if v.any() else len(t) - 1) for k, v in hits.items()}
    hit = min(first, key=first.get)
    i = first[hit]
    p, tt = pts[i], t[i]
    v = v0 + np.array([0.0, -G]) * tt
    surf = w * np.array([-p[1], p[0]])                       # what it lands on turns with the drum
    nrm = e_bed if hit == "bed" else p / np.linalg.norm(p)
    v_rel = v - surf
    E = M_ROCK * (G * (p0[1] + r_c) + 0.5 * (v0 @ v0 - (w * r_c) ** 2))
    return dict(rib=rib, phi=phi, th_release=math.degrees(th_r), hit=hit,
                th_land=math.degrees(math.atan2(p[0], -p[1])),
                v_impact=float(np.linalg.norm(v_rel)), v_normal=float(abs(v_rel @ nrm)),
                drop=float(p0[1] - p[1]), t_flight=float(tt), E_lift=float(E),
                path=np.vstack([p0, pts[:i + 1]]) if keep_path else None)

def carry(n_stones=N_STONES):
    return min(CARRY_MAX, max(1.0, n_stones / 4))

def drum_torque(rib, phi, n_stones=N_STONES):
    """Mean torque the drum needs (N m); zero once rocks are no longer released."""
    L = lift(rib, phi, n_stones)
    if L is None:
        return 0.0
    return N_RIBS * carry(n_stones) * max(L["E_lift"], 0.0) / (2 * math.pi)

def strikes_per_s(phi, n_stones=N_STONES):
    return N_RIBS * carry(n_stones) * phi * omega_c() / (2 * math.pi)

def quartz_strike(rib, phi, n_stones=N_STONES):
    """Normal impact speed when the flight ends on quartz, else 0."""
    L = lift(rib, phi, n_stones)
    return L["v_normal"] if (L and L["hit"] in QUARTZ) else 0.0

if __name__ == "__main__":
    wc = omega_c()
    f, alpha, d, _ = bed()
    print(f"drum R {R_IN*1e3:.0f} mm (liner surface {R_W*1e3:.0f}) x {DRUM_L*1e3:.0f} mm; "
          f"rock centres at r = {(R_W-ROCK_A)*1e3:.0f} mm")
    print(f"critical speed {rpm(wc):.0f} rpm; target {PHI_TARGET:.0%} = {rpm(PHI_TARGET*wc):.0f} rpm")
    print(f"{N_STONES} stones of {2*ROCK_A*1e3:.0f} mm: bulk fill {f:.1%}; bed spans "
          f"{math.degrees(PSI_BED-alpha):.0f} to {math.degrees(PSI_BED+alpha):.0f} deg from the bottom")
    print()
    print(" rib       phi  release  lands   at deg  drop mm  v_normal m/s  torque mN m")
    for rib in RIBS:
        for phi in (0.55, 0.65, 0.75, 0.85, 0.95):
            L = lift(rib, phi)
            if L is None:
                print(f" {rib:8s} {phi:4.2f}   carried round (centrifuging)")
                continue
            print(f" {rib:8s} {phi:4.2f} {L['th_release']:7.0f}  {L['hit']:6s} {L['th_land']:7.0f} "
                  f"{L['drop']*1e3:8.0f} {L['v_normal']:12.2f} {drum_torque(rib, phi)*1e3:12.1f}")
        print()
    print("best quartz strike per rib (normal speed on liner or bed):")
    phis = np.round(np.arange(0.50, 0.96, 0.05), 2)
    for rib in RIBS:
        vs = [quartz_strike(rib, p) for p in phis]
        k = int(np.argmax(vs))
        band = [p for p, v in zip(phis, vs) if v >= 0.8 * vs[k]]
        print(f"  {rib:8s} peak {vs[k]:.2f} m/s at {phis[k]:.0%} of critical; >=80 % of peak over "
              f"{min(band):.0%}-{max(band):.0%}")
    print(f"\nfill sweep, {DEFAULT_RIB} at {PHI_TARGET:.0%} (severity vs count):")
    for n in (4, 6, 8, 10, 14, 18, 24):
        L = lift(DEFAULT_RIB, PHI_TARGET, n)
        print(f"  {n:2d} stones  fill {bed(n)[0]:5.1%}  lands on {L['hit']:5s}  v_normal {L['v_normal']:.2f} m/s  "
              f"strikes/s {strikes_per_s(PHI_TARGET, n):5.1f}")
    T = drum_torque(DEFAULT_RIB, PHI_TARGET)
    print(f"\ndrum torque for {DEFAULT_RIB} at {PHI_TARGET:.0%}: {T*1e3:.1f} mN m")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cols = {"h08_f90": "#7a9cc6", "h14_f90": "#1d3f66", "h14_f60": "#c0552b", "h14_f45": "#e0a030"}
    fig, ax = plt.subplots(2, 2, figsize=(11, 10))
    a0 = ax[0, 0]
    t = np.linspace(0, 2 * math.pi, 300)
    a0.plot(R_IN * 1e3 * np.cos(t), R_IN * 1e3 * np.sin(t), color="0.3", lw=1)
    a0.fill_between(R_W * 1e3 * np.cos(t), 0, 0, color="none")
    a0.plot(R_W * 1e3 * np.cos(t), R_W * 1e3 * np.sin(t), color="#b9a77e", lw=4, alpha=0.6,
            label="quartz liner")
    f, alpha, d, e = bed()
    ang = math.atan2(e[1], e[0])
    seg = np.array([[R_W * math.cos(ang + u), R_W * math.sin(ang + u)] for u in np.linspace(-alpha, alpha, 60)])
    a0.fill(seg[:, 0] * 1e3, seg[:, 1] * 1e3, color="#e8dcc8", label=f"charge bed ({N_STONES} stones)")
    for rib, c in cols.items():
        L = lift(rib, PHI_TARGET, keep_path=True)
        if L is not None:
            a0.plot(L["path"][:, 0] * 1e3, L["path"][:, 1] * 1e3, color=c, label=f"{rib}: lands on {L['hit']}")
    a0.add_patch(plt.Circle((0, 0), AXLE_R * 1e3, color="0.5"))
    a0.set_aspect("equal"); a0.set_xlim(-65, 65); a0.set_ylim(-65, 65)
    a0.set_title(f"Flights at {PHI_TARGET:.0%} of critical (drum turns CCW)")
    a0.legend(fontsize=7, loc="upper left")
    ph = np.linspace(0.45, 1.0, 45)
    for rib, c in cols.items():
        Ls = [lift(rib, p) for p in ph]
        ax[0, 1].plot([p for p, L in zip(ph, Ls) if L], [L["th_release"] for L in Ls if L], color=c, label=rib)
        ax[1, 0].plot(ph, [L["v_normal"] if (L and L["hit"] in QUARTZ) else np.nan for L in Ls], color=c, label=rib)
        ax[1, 1].plot(ph, [drum_torque(rib, p) * 1e3 for p in ph], color=c, label=rib)
    for a in (ax[0, 1], ax[1, 0], ax[1, 1]):
        a.axvspan(0.70, 0.80, color="#f3d9c8"); a.set_xlabel("fraction of critical speed"); a.legend(fontsize=8)
    ax[0, 1].set_ylabel("release angle from bottom, deg"); ax[0, 1].set_title("Release angle (sloped face lets go earlier)")
    ax[1, 0].set_ylabel("normal impact speed on quartz, m/s")
    ax[1, 0].set_title("Strike severity (steps = landing moves from the bed to the liner)")
    ax[1, 1].set_ylabel("drum torque, mN m")
    ax[1, 1].set_title("Drum torque: flat, then collapse at critical")
    fig.tight_layout()
    fig.savefig("tumbler.png", dpi=110)
    print("wrote tumbler.png")
