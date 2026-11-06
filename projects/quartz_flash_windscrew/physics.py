#!/usr/bin/env python3
"""Quartz Flash Windscrew -- the numbers behind the build.

Run from this directory:  python physics.py
Prints the design tables and writes physics.png.  The drum itself (lifter
release, flights, landing, drum torque) is modelled in tumbler.py.

Chain of models (assumptions in README.md):
  1. Rotor     Liam-type conical Archimedes screw, linear torque-coefficient curve
               Ct(lam) = Ct0 (1 - lam/lam0), conservative small-rotor values.
  2. Gearbox   4.8 : 1 planetary (sun in, ring fixed, carrier = drum), efficiency
               eta, plus a constant drag torque at the sun.
  3. Governor  inverse centrifugal clutch.  Three shoes on the sun's cage are held
               onto a hub on the shaft by an O-ring (inward force F each).  Torque
               passed to the sun is 3 mu (F - m w^2 r_s) r_hub, which falls as the sun
               speeds up, so the sun settles where that equals the drum's load:
               w_sun^2 = (F - T_sun / (3 mu r_hub)) / (m r_s).  Adding nuts to the
               shoes (m) lowers the governed speed; wind no longer sets it.
  4. Drum load from tumbler.py (lifts per rev x energy per lift).
  5. Impact    Hertz sphere-on-sphere impact: peak force, contact radius, peak
               pressure, and the rim tension that starts a cone crack.
  6. Sliding   flash temperature of a sliding contact, Tian & Kennedy (1994)
               circular-source formula with Blok heat partition.
  7. Photons   Planck emission from the hot real-contact area, scotopic-weighted, into
               a dark-adapted pupil, against the Hecht-Shlaer-Pirenne threshold -- swept
               over emissivity, the most uncertain number here.
"""
import math
import numpy as np
import tumbler as TB

RHO_AIR = 1.20             # kg/m^3 (ridge, ~200 m ASL, 15 C)
G = 9.81
H_PL, C_L, K_B = 6.626e-34, 2.998e8, 1.381e-23

# rotor (matches windscrew.py)
R_ROT = 0.150
A_ROT = math.pi * R_ROT**2
CT0, LAM0 = 0.20, 2.0      # -> Cp_max = Ct0 lam0 / 4 = 0.10 (independent Liam tests: ~0.1-0.2)
T_BEARING = 3e-3           # N m: two axle 608s and the drum's two 608s

# gearbox + governor (matches windscrew.py)
RATIO = 1 + 57 / 15
ETA_GEAR = 0.75            # printed PETG planetary, light load
T_DRAG_SUN = 1.5e-3        # N m at the sun: two 688s, planet pins, mesh friction
N_SHOES = 3
MU_PAD = 0.40              # cork on PETG
R_HUB = 0.016
R_SHOE = 0.027             # shoe centre of mass radius
M_SHOE = 0.0039            # printed shoe (3.1 cm3 PETG)
M_NUT = 0.0012             # M5 nut; up to 4 per shoe
NUTS_NOMINAL = 2           # O-ring force is chosen so this setting gives PHI_TARGET

# quartz (alpha, room temperature, directional averages)
E_Q, NU_Q = 95e9, 0.08
K_Q = 7.0
KAPPA_Q = K_Q / (2650 * 740)
H_Q = 10e9
MU_Q = 0.5
SIGMA_F = 70e6
R_STONE = 0.009
T_AMB = 288.0
T_ALPHA_BETA = 846.0
T_CEIL = 1950.0

# eye
PUPIL_D = 7e-3
VIEW_L = 0.30              # m, eye to the drum at night
T_PETG = 0.85              # transmission of one clear 2.4 mm PETG wall
HSP = (54, 148)

_SC = np.array([[400, 9.29e-3], [420, 9.66e-2], [440, 0.328], [460, 0.567], [480, 0.793],
                [500, 0.982], [507, 1.0], [520, 0.935], [540, 0.650], [560, 0.329],
                [580, 0.121], [600, 3.32e-2], [620, 7.37e-3], [640, 1.50e-3],
                [660, 3.13e-4], [680, 7.15e-5], [700, 1.78e-5]])

# drum torque vs fraction of critical, tabulated once from tumbler.py
_PHI = np.linspace(0.02, 1.0, 50)
_TD = np.array([TB.drum_torque(TB.DEFAULT_RIB, p) for p in _PHI])
W_C = TB.omega_c()

def drum_torque(phi):
    return float(np.interp(phi, _PHI, _TD))

def ct(lam):
    return CT0 * np.clip(1 - lam / LAM0, 0, None)

def rotor_torque(v, lam):
    return 0.5 * RHO_AIR * A_ROT * R_ROT * v**2 * ct(lam)

def t_sun(phi):
    """Torque the sun must deliver to hold the drum at phi."""
    return drum_torque(phi) / (RATIO * ETA_GEAR) + T_DRAG_SUN

def shoe_mass(nuts):
    return M_SHOE + nuts * M_NUT

def oring_force():
    """Inward force per shoe that puts the nominal setting at PHI_TARGET."""
    w = TB.PHI_TARGET * W_C * RATIO
    return shoe_mass(NUTS_NOMINAL) * w * w * R_SHOE + t_sun(TB.PHI_TARGET) / (N_SHOES * MU_PAD * R_HUB)

def governed_phi(nuts, torque_scale=1.0, F=None):
    """Drum speed (fraction of critical) the governor settles at; solves the load/speed loop."""
    F = oring_force() if F is None else F
    phi = TB.PHI_TARGET
    for _ in range(60):
        L = torque_scale * t_sun(phi) / (N_SHOES * MU_PAD * R_HUB)
        w2 = max(F - L, 0.0) / (shoe_mass(nuts) * R_SHOE)
        phi = 0.5 * phi + 0.5 * math.sqrt(w2) / (RATIO * W_C)
    return phi

def operating_point(v, nuts=NUTS_NOMINAL):
    """(rotor omega, drum phi, governed?) at wind v."""
    phi_g = governed_phi(nuts)
    lam = np.linspace(1e-3, LAM0, 3000)
    om = lam * v / R_ROT
    phi_d = np.minimum(om / (RATIO * W_C), phi_g)
    load = np.array([t_sun(p) for p in phi_d]) + T_BEARING
    net = rotor_torque(v, lam) - load
    i = np.where(np.diff(np.sign(net)) < 0)[0]
    if len(i) == 0:
        return 0.0, 0.0, False
    k = i[0]
    return om[k], phi_d[k], bool(om[k] / (RATIO * W_C) > phi_g)

def cut_in():
    """Rotor static torque must lift the resting charge through the gearbox."""
    charge = TB.N_STONES * TB.M_ROCK * G * (TB.R_W - TB.ROCK_A) * math.sin(TB.PSI_BED)
    need = charge / (RATIO * ETA_GEAR) + T_DRAG_SUN + T_BEARING
    return math.sqrt(need / (0.5 * RHO_AIR * A_ROT * R_ROT * CT0))

# axle (matches windscrew.py): world x of the bearings and loads, mm
X_NOSE, X_MOUTH_IN, X_MOUTH_OUT = -163.5, 163.5, 194.5
X_SUN, X_DRUM_HUB, DRUM_ARM = 233.0, 246.5, 41.0     # drum+lid+stones centroid sits DRUM_ARM beyond its hub
W_DRUM, W_GEAR, W_ROTOR, F_GUST = 9.0, 1.0, 5.0, 22.0
RODS = {"steel M8 threaded (root 6.47 mm)": (200e3, 6.47, 640.0, 250.0),
        "aluminium 8 mm 6061-T6 / 6082-T6": (69e3, 8.0, 275.0, 96.0)}   # E MPa, d, yield, fatigue (polished)

def shaft_beam(E, d, supports, gust=False, n=480):
    """Euler-Bernoulli FE beam of the axle.  Returns x, deflection (mm), bending stress (MPa)."""
    xs = np.linspace(-215.0, 262.0, n + 1)
    I = math.pi * d**4 / 64
    K = np.zeros((2 * len(xs), 2 * len(xs))); F = np.zeros(2 * len(xs))
    for e in range(n):
        L = xs[e + 1] - xs[e]
        k = E * I / L**3 * np.array([[12, 6 * L, -12, 6 * L], [6 * L, 4 * L * L, -6 * L, 2 * L * L],
                                     [-12, -6 * L, 12, -6 * L], [6 * L, 2 * L * L, -6 * L, 4 * L * L]])
        i = [2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3]
        K[np.ix_(i, i)] += k
    node = lambda x: int(np.argmin(abs(xs - x)))
    loads = [(X_DRUM_HUB, -W_DRUM), (X_SUN, -W_GEAR), (-60.0, -W_ROTOR)] + ([(-30.0, -F_GUST)] if gust else [])
    for x, P in loads:
        F[2 * node(x)] += P
    F[2 * node(X_DRUM_HUB) + 1] += -W_DRUM * DRUM_ARM
    fixed = [2 * node(x) for x in supports]
    free = [i for i in range(len(F)) if i not in fixed]
    u = np.zeros(len(F))
    u[free] = np.linalg.solve(K[np.ix_(free, free)], F[free])
    M = E * I * np.gradient(u[1::2], xs)
    return xs, u[0::2], M / (math.pi * d**3 / 32)

def hertz_impact(v_n, m=TB.M_ROCK, r1=R_STONE, r2=R_STONE):
    e_star = E_Q / (2 * (1 - NU_Q**2))
    r_star = r1 * r2 / (r1 + r2)
    d = (15 * m * v_n**2 / (16 * e_star * math.sqrt(r_star)))**0.4
    f = 4 / 3 * e_star * math.sqrt(r_star) * d**1.5
    a = math.sqrt(r_star * d)
    p0 = 3 * f / (2 * math.pi * a**2)
    return f, a, p0, (1 - 2 * NU_Q) / 3 * p0

def flash_rise(q, a, v):
    """Tian & Kennedy (1994), uniform circular source, flux q, Pe = v a / 2 kappa on the
    moving body; Blok partition against a body that sees the source stationary."""
    pe = v * a / (2 * KAPPA_Q)
    t_moving = 2 * q * a / (K_Q * math.sqrt(math.pi * (1.273 + pe)))
    t_still = 2 * q * a / (K_Q * math.sqrt(math.pi * 1.273))
    return 1 / (1 / t_moving + 1 / t_still)

def flash_T(q, a, v):
    return min(T_AMB + flash_rise(q, a, v), T_CEIL)

def photons_to_eye(T, area, tau, eps, dist=VIEW_L):
    lam_nm = np.linspace(400, 700, 301)
    lam = lam_nm * 1e-9
    mp = 2 * math.pi * C_L / lam**4 / np.expm1(H_PL * C_L / (lam * K_B * T))
    w = np.interp(lam_nm, _SC[:, 0], _SC[:, 1])
    n_scot = np.trapezoid(mp * w, lam) * area * tau * eps
    frac = (PUPIL_D**2 / 4) / dist**2 * T_PETG                  # Lambertian, on-axis
    return n_scot * frac

if __name__ == "__main__":
    F = oring_force()
    print("== governor ==")
    print(f"drum torque ({TB.DEFAULT_RIB}, {TB.N_STONES} stones) at {TB.PHI_TARGET:.0%}: "
          f"{drum_torque(TB.PHI_TARGET)*1e3:.1f} mN m -> sun {t_sun(TB.PHI_TARGET)*1e3:.1f} mN m")
    print(f"O-ring inward force per shoe for {NUTS_NOMINAL} nuts at {TB.PHI_TARGET:.0%}: {F:.2f} N "
          f"(band tension ~{F / (2 * math.sin(math.pi / 3)):.2f} N)")
    print(" nuts/shoe   drum rpm   phi    phi at -30 % / +30 % drum torque")
    for n in range(5):
        p = governed_phi(n)
        print(f"   {n}        {TB.rpm(p * W_C):6.0f}   {p:4.2f}     {governed_phi(n, 0.7):4.2f} / {governed_phi(n, 1.3):4.2f}")
    print(f"\ncut-in wind: {cut_in():.1f} m/s")
    print("\n wind  rotor rpm  drum rpm   phi  governed  strikes/s  quartz strike m/s  slip W")
    for v in (3, 4, 5, 6, 8, 10, 12, 15, 20, 25):
        om, phi, gov = operating_point(v)
        if om == 0:
            print(f"{v:5.0f}   (below cut-in)")
            continue
        slip = t_sun(phi) * max(om - phi * W_C * RATIO, 0)
        vq = TB.quartz_strike(TB.DEFAULT_RIB, phi) if phi > 0.3 else 0.0
        print(f"{v:5.0f} {om*30/math.pi:9.0f} {TB.rpm(phi*W_C):9.0f} {phi:6.2f}  {'yes' if gov else 'no ':8s}"
              f"{TB.strikes_per_s(phi):9.1f} {vq:15.2f} {slip:10.2f}")

    v_imp = TB.quartz_strike(TB.DEFAULT_RIB, TB.PHI_TARGET)
    f, a, p0, st = hertz_impact(v_imp)
    print(f"\n== a governed strike: {v_imp:.2f} m/s normal ==")
    print(f"Hertz peak force {f:.0f} N, contact radius {a*1e6:.0f} um, p0 {p0/1e9:.2f} GPa, "
          f"rim tension {st/SIGMA_F:.0f}x the surface strength -> cone cracks, fracture light")

    v_s = TB.PHI_TARGET * W_C * (TB.R_W - TB.ROCK_A)
    n = TB.M_ROCK * G * (1 + TB.PHI_TARGET**2)
    print(f"\n== sliding on the liner: {v_s:.2f} m/s under ~{n:.2f} N ==")
    for asp in (3e-6, 10e-6):
        print(f"  asperity a = {asp*1e6:.0f} um: flash T {flash_T(MU_Q * H_Q * v_s, asp, v_s):.0f} K")
    area, tau = n / H_Q, 0.01 / v_s
    T = flash_T(MU_Q * H_Q * v_s, 3e-6, v_s)
    eps_need = next((e for e in np.logspace(-8, 0, 800) if photons_to_eye(T, area, tau, e) > HSP[0]), None)
    print(f"  real contact {area*1e12:.0f} um2 at {T:.0f} K for {tau*1e3:.0f} ms; emissivity for the HSP floor: "
          f"{eps_need:.1e}" if eps_need else "  thermal channel below threshold at any emissivity")

    print("\n== axle: sag at the sun gear (backlash 0.14-0.20 mm) and bending stress ==")
    for name, (E, d, sy, sf) in RODS.items():
        for sup, lab in (([X_NOSE, X_MOUTH_IN], "one mouth bearing "), ([X_NOSE, X_MOUTH_IN, X_MOUTH_OUT], "two mouth bearings")):
            xs, w, s = shaft_beam(E, d, sup)
            _, _, sg = shaft_beam(E, d, sup, gust=True)
            i = int(np.argmin(abs(xs - X_SUN)))
            rot = abs(s[xs > 150]).max()
            print(f"  {name:34s} {lab}: sag {abs(w[i]):.3f} mm, rotating {rot:4.1f} MPa "
                  f"(fatigue {sf:.0f} polished), {F_GUST:.0f} N gust peak {abs(sg).max():4.1f} MPa (yield {sy:.0f})")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))
    winds = np.linspace(2, 25, 90)
    ops = [operating_point(v) for v in winds]
    ax[0, 0].plot(winds, [o[0] * 30 / math.pi for o in ops], color="#3a6ea5", label="rotor rpm")
    ax[0, 0].set_xlabel("wind m/s"); ax[0, 0].set_ylabel("rotor rpm", color="#3a6ea5")
    axb = ax[0, 0].twinx()
    axb.plot(winds, [o[1] for o in ops], color="#c0552b")
    axb.axhspan(0.70, 0.80, color="#f3d9c8", alpha=0.5)
    axb.set_ylabel("drum speed, fraction of critical", color="#c0552b"); axb.set_ylim(0, 1.05)
    ax[0, 0].axvline(cut_in(), color="0.6", ls="--")
    ax[0, 0].set_title("Governed drum: wind sets the rotor, not the drum")
    nuts = np.arange(0, 5)
    for sc, ls, lab in ((1.0, "-", "nominal drum torque"), (0.7, ":", "-30 %"), (1.3, "--", "+30 %")):
        ax[0, 1].plot(nuts, [governed_phi(n, sc) for n in nuts], ls, color="#1d3f66", marker="o", label=lab)
    ax[0, 1].axhspan(0.70, 0.80, color="#f3d9c8")
    ax[0, 1].set_xlabel("M5 nuts per shoe"); ax[0, 1].set_ylabel("governed fraction of critical")
    ax[0, 1].legend(fontsize=8); ax[0, 1].set_title("Tuning the governor")
    vv = np.logspace(-1, 1, 100)
    for asp, ls in ((1e-6, "--"), (3e-6, "-"), (10e-6, "-.")):
        ax[1, 0].plot(vv, [flash_T(MU_Q * H_Q * v, asp, v) for v in vv], ls, color="#c0552b",
                      label=f"asperity a={asp*1e6:.0f} um")
    ax[1, 0].axhline(T_ALPHA_BETA, color="0.5", ls=":"); ax[1, 0].text(0.11, 880, "alpha-beta 573 C", fontsize=8)
    ax[1, 0].axvspan(0.1, 0.3, color="0.9"); ax[1, 0].text(0.105, 1700, "hand\nrubbing", fontsize=8)
    ax[1, 0].axvline(v_s, color="#1d3f66", ls=":"); ax[1, 0].text(v_s * 1.05, 400, "drum\nsliding", fontsize=8)
    ax[1, 0].set_xscale("log"); ax[1, 0].set_xlabel("sliding speed m/s"); ax[1, 0].set_ylabel("flash T, K")
    ax[1, 0].legend(fontsize=8); ax[1, 0].set_title("Flash temperature (Tian-Kennedy, capped at softening)")
    vi = np.linspace(0.1, 2.0, 60)
    ax[1, 1].plot(vi, [hertz_impact(x)[3] / SIGMA_F for x in vi], color="#1d3f66")
    ax[1, 1].axhline(1, color="#c0552b", ls=":"); ax[1, 1].text(0.12, 1.15, "cone-crack threshold", fontsize=8)
    ax[1, 1].axvline(v_imp, color="0.5", ls="--"); ax[1, 1].text(v_imp + 0.03, 2, "governed strike", fontsize=8)
    ax[1, 1].set_yscale("log"); ax[1, 1].set_xlabel("normal impact speed m/s")
    ax[1, 1].set_ylabel("Hertz rim tension / surface strength")
    ax[1, 1].set_title("Every strike cracks quartz (fracture light)")
    fig.tight_layout()
    fig.savefig("physics.png", dpi=110)
    print("\nwrote physics.png")
