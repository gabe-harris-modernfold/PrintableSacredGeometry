"""Build every part, verify it, and print the slicing sheet."""

import math, os
import numpy as np
import trimesh
import params as P
import vessel, cascade, spiral

BED = P.BED

def panes():
    """Flat-printed viewing panes, laid out on the bed. Printing them separately and
    flat is what buys the clarity: two optically flat faces instead of the corrugation
    a tilted in-place wall would get."""
    made, seen = [], {}
    for mod in range(3):
        _, ports = vessel.shell(mod, verbose=False)
        for f, z0, z1 in ports:
            if f != 0: continue                       # one per band; x12 identical
            a = math.cos(math.pi / P.N_FACE)
            r0, r1 = float(P.vessel_r(z0)) * a, float(P.vessel_r(z1)) * a
            w0 = r0 * math.tan(math.pi / P.N_FACE) - vessel.FRAME_SIDE + 4
            w1 = r1 * math.tan(math.pi / P.N_FACE) - vessel.FRAME_SIDE + 4
            h = math.hypot(z1 - z0, r1 - r0)
            key = (round(w0, 1), round(w1, 1), round(h, 1))
            if key in seen: continue
            seen[key] = True
            V = np.array([[-w0, 0, 0], [w0, 0, 0], [w1, h, 0], [-w1, h, 0],
                          [-w0, 0, P.PANE_T], [w0, 0, P.PANE_T],
                          [w1, h, P.PANE_T], [-w1, h, P.PANE_T]])
            F = np.array([[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],
                          [1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]])
            m = trimesh.Trimesh(V, F, process=False)
            if m.volume < 0: m.invert()
            made.append((f"pane_m{mod}_b{len(made)}", m, 12, w0, w1, h))
    return made


def main():
    parts = []
    for k in range(3):
        m, _ = vessel.shell(k, verbose=False)
        m.export(f"vessel_shell_{k}.stl"); parts.append((f"vessel_shell_{k}", m, 1))
    r = vessel.reservoir(verbose=False); r.export("vessel_reservoir.stl")
    parts.append(("vessel_reservoir", r, 1))
    b = vessel.bed_tray(verbose=False); b.export("vessel_bedtray.stl")
    parts.append(("vessel_bedtray (+siphon)", b, 1))
    c = cascade.build(verbose=False)[0]; c.export("cascade_screen.stl")
    parts.append(("cascade_screen", c, 1))
    g = spiral.build(verbose=False, part="gutters"); g.export("drip_gutters.stl")
    parts.append(("drip_gutters (16 spirals)", g, 1))
    sp = spiral.build(verbose=False, part="splitter"); sp.export("drip_splitter.stl")
    parts.append(("drip_splitter (10->16)", sp, 1))
    for name, m, n, *_ in panes():
        m.export(f"{name}.stl"); parts.append((name, m, n))

    print()
    print("=" * 82)
    print(f"{'part':30s} {'x':>6s} {'y':>6s} {'z':>6s} {'cm3':>7s} {'qty':>4s}  "
          f"{'fits 320':>8s}  {'solid':>6s}")
    print("-" * 82)
    tot_v = 0.0
    for name, m, n, *_ in parts:
        e = m.extents
        fits = "yes" if max(e) <= BED and sorted(e)[-2] <= BED else "NO"
        wt = "yes" if m.is_watertight else "no"
        tot_v += m.volume / 1000 * n
        print(f"{name:30s} {e[0]:6.0f} {e[1]:6.0f} {e[2]:6.0f} "
              f"{m.volume/1000:7.0f} {n:4d}  {fits:>8s}  {wt:>6s}")
    print("-" * 82)
    print(f"{'TOTAL':30s} {'':6s} {'':6s} {'':6s} {tot_v:7.0f} cm3  "
          f"= {tot_v*1.27/1000:.2f} kg PETG")
    print("=" * 82)
    print(f"assembled: {2*P.R_BODY:.0f} across corners x {P.VESSEL_H:.0f} tall "
          f"(footprint limit 320, bed {BED:.0f})")
    print(f"reservoir {5.9:.1f} L to the waterline at z={P.Z_RES:.0f}; "
          f"living bed at z={P.Z_BED:.0f} on a 32 deg trickle slope")
    tf = math.sqrt(2 * P.STEP * 1e-3 / P.GRAV)
    dps = P.Q_TRICKLE / 3600 * 1e-3 / (P.DROP_V * 1e-9)
    casc_sites = P.LEVELS * 217
    print(f"drip sites: {casc_sites} on the cascade + 1680 on the spirals = "
          f"{casc_sites+1680}")
    print(f"drops: {P.DROP_V:.1f} uL / {P.DROP_D:.2f} mm off a {P.LIP_W} mm lip; "
          f"{dps:.0f}/s crossing each level")
    print(f"airborne at any instant: ~{dps*(P.LEVELS+28)*tf:.0f} at trickle, "
          f"~{dps*8*(P.LEVELS+28)*tf:.0f} during a siphon dump")

if __name__ == "__main__":
    main()
