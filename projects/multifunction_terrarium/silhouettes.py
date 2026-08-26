"""Ten candidate silhouettes. Alchemical-glassware family, tapering upward,
cascade tiers NESTED concentrically rather than stacked.

Constraints held in every one: max half-width 155 (310 across corners <= 320 bed
and footprint), overall height ~640, and enough belly to hold three nested cascade
cones plus a base reservoir."""

import numpy as np
from scipy.interpolate import PchipInterpolator

H, RMAX = 640.0, 155.0
CW, CH = 330, 500          # cell size, px
MARG_T, MARG_B = 46, 34
SC = (CH - MARG_T - MARG_B) / H          # px per mm

def prof(ctrl, smooth=True, n=260):
    z = np.array([c[0] for c in ctrl], float)
    r = np.array([c[1] for c in ctrl], float)
    zz = np.linspace(z.min(), z.max(), n)
    return zz, (PchipInterpolator(z, r)(zz) if smooth else np.interp(zz, z, r))

DESIGNS = [
 ("Alembic", "cucurbit belly, hard shoulder, still-head above", True, [
   (0,96),(14,124),(34,152),(120,155),(235,149),(298,122),(358,80),(418,58),
   (478,55),(518,74),(558,82),(600,52),(640,20)]),
 ("Cucurbit + cone", "straight body, conical head - the plainest still", False, [
   (0,110),(18,150),(40,155),(250,155),(268,146),(300,120),(340,92),(400,64),
   (470,42),(545,26),(610,12),(640,6)]),
 ("Double gourd", "two bellies, upper smaller - one chamber each", True, [
   (0,100),(20,136),(70,155),(150,150),(232,116),(282,80),(322,70),(372,88),
   (424,108),(472,111),(522,90),(572,56),(620,28),(640,16)]),
 ("Ogee / onion", "reflex curve, reliquary dome, no straight run anywhere", True, [
   (0,120),(20,150),(62,155),(142,149),(222,130),(292,103),(342,76),(382,57),
   (412,48),(442,55),(472,70),(502,79),(532,72),(562,53),(592,30),(618,13),(640,4)]),
 ("Nested bells", "one concave sweep - the nesting IS the silhouette", True, [
   (0,155),(58,151),(138,137),(220,114),(302,91),(382,69),(462,51),(540,35),
   (600,23),(640,13)]),
 ("Monstrance", "stepped foot, knop, radiant head - most ceremonial", False, [
   (0,155),(24,150),(38,120),(58,116),(74,86),(108,60),(138,48),(158,68),
   (184,50),(210,45),(300,58),(380,92),(458,114),(518,106),(568,78),(618,38),(640,16)]),
 ("Pagoda eaves", "flared drip eaves - drops fall clear into open air", False, [
   (0,140),(28,140),(28,155),(58,155),(58,126),(128,126),(128,150),(156,150),
   (156,114),(226,114),(226,136),(252,136),(252,98),(322,98),(322,118),(346,118),
   (346,80),(416,80),(416,98),(438,98),(438,60),(504,60),(504,76),(524,76),
   (524,40),(586,40),(610,22),(640,10)]),
 ("Pinnacle", "single concave taper, stalagmite / gothic spire", True, [
   (0,155),(80,139),(160,119),(240,98),(320,78),(400,59),(480,41),(550,25),
   (610,11),(640,4)]),
 ("Pelican", "belly with two return arms - the vessel that recirculates", True, [
   (0,104),(18,138),(60,155),(140,153),(215,132),(268,104),(300,86),(330,76),
   (378,70),(430,60),(492,46),(556,30),(610,14),(640,6)]),
 ("Beaked flask", "asymmetric retort - a front, a back, and a pouring side", True, [
   (0,100),(20,136),(66,155),(148,152),(226,128),(286,96),(330,72),(374,55),
   (430,44),(492,34),(552,24),(604,13),(640,5)]),
]

# asymmetric right-hand overrides (retort neck sweep, pelican arms)
ASYM = {
 "Beaked flask": [(0,100),(20,138),(66,155),(150,155),(232,140),(296,118),
                  (348,102),(400,92),(452,80),(508,62),(560,42),(606,20),(640,7)],
}
ARMS = {"Pelican": True}


def stair(z0, z1, r0, r1, n=9):
    """A terraced cascade cone in profile: narrow at top, widening as it falls."""
    pts = []
    for i in range(n + 1):
        f = i / n
        z = z1 + (z0 - z1) * f
        r = r1 + (r0 - r1) * f
        pts.append((z, r))
        if i < n:
            pts.append((z - (z1 - z0) / n * 0 - (z1 - z0) / n, r))
    out = []
    for i in range(n):
        f0, f1 = i / n, (i + 1) / n
        za, zb = z1 - (z1 - z0) * f0, z1 - (z1 - z0) * f1
        ra, rb = r1 + (r0 - r1) * f0, r1 + (r0 - r1) * f1
        out += [(za, ra), (zb, ra), (zb, rb)]
    return out


def svg():
    cols, rows = 5, 2
    W, Hh = cols * CW, rows * CH + 66
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{Hh}" '
         f'viewBox="0 0 {W} {Hh}" font-family="Helvetica,Arial,sans-serif">',
         f'<rect width="{W}" height="{Hh}" fill="#f7f5f1"/>',
         f'<text x="20" y="30" font-size="19" fill="#1d1d1b">Terrarium silhouettes '
         f'&#8212; alchemical vessel, tapering, nested cascade tiers</text>',
         f'<text x="20" y="50" font-size="12.5" fill="#6d6a65">every profile holds '
         f'310 across corners &#215; 640 tall &#8212; grey = nested cascade cones, '
         f'dashed = 320 mm bed / footprint limit</text>']
    for i, (name, note, smooth, ctrl) in enumerate(DESIGNS):
        cx = (i % cols) * CW + CW / 2
        cy = (i // cols) * CH + 66
        base = cy + CH - MARG_B
        zz, rr = prof(ctrl, smooth)
        if name in ASYM:
            z2, r2 = prof(ASYM[name], True)
            rr2 = np.interp(zz, z2, r2)
        else:
            rr2 = rr
        L = [(cx - r * SC, base - z * SC) for z, r in zip(zz, rr)]
        R = [(cx + r * SC, base - z * SC) for z, r in zip(zz, rr2)][::-1]
        path = "M" + " L".join(f"{x:.1f},{y:.1f}" for x, y in L + R) + " Z"

        o.append(f'<g><clipPath id="c{i}"><path d="{path}"/></clipPath>')
        o.append(f'<path d="{path}" fill="#20262b"/>')
        # nested cascade cones, clipped inside the vessel
        for (z0, z1, r0, r1) in ((0.02, 0.40, 0.92, 0.30), (0.30, 0.68, 0.66, 0.22),
                                 (0.58, 0.93, 0.44, 0.13)):
            pts = stair(z0 * H, z1 * H, r0 * RMAX, r1 * RMAX)
            for sgn in (-1, 1):
                d = "M" + " L".join(f"{cx+sgn*r*SC:.1f},{base-z*SC:.1f}" for z, r in pts)
                o.append(f'<path d="{d}" fill="none" stroke="#8fb7cc" '
                         f'stroke-width="1.5" clip-path="url(#c{i})" opacity=".95"/>')
        o.append('</g>')
        # footprint rule
        o.append(f'<line x1="{cx-RMAX*SC:.1f}" y1="{base+7}" x2="{cx+RMAX*SC:.1f}" '
                 f'y2="{base+7}" stroke="#b9b3aa" stroke-width="1" '
                 f'stroke-dasharray="4 3"/>')
        o.append(f'<text x="{cx}" y="{cy+18}" font-size="14.5" font-weight="600" '
                 f'text-anchor="middle" fill="#1d1d1b">{i+1}. {name}</text>')
        o.append(f'<text x="{cx}" y="{cy+34}" font-size="11" text-anchor="middle" '
                 f'fill="#6d6a65">{note}</text>')
    o.append('</svg>')
    return "\n".join(o)

open("silhouettes.svg", "w", encoding="utf-8").write(svg())
print("wrote silhouettes.svg")
