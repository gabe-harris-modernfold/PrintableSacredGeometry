"""Round two: named alchemical apparatus, drawn as real vessels with their
appendages. Asymmetry allowed. All hold 310 across corners x ~660 tall, with three
cascade cones NESTED concentrically (inner = higher + narrower, spilling outward)."""

import numpy as np
from scipy.interpolate import PchipInterpolator

H, RMAX = 660.0, 155.0

def prof(ctrl, n=300, zmax=H):
    z = np.array([c[0] for c in ctrl], float); r = np.array([c[1] for c in ctrl], float)
    zz = np.linspace(0, zmax, n)
    return zz, np.clip(PchipInterpolator(z, r, extrapolate=True)(zz), 0.5, None)

def tube(pts, w0, w1, n=60):
    """Swept round tube along a quadratic bezier through 3 control points."""
    p = np.array(pts, float); t = np.linspace(0, 1, n)[:, None]
    c = (1-t)**2*p[0] + 2*(1-t)*t*p[1] + t**2*p[2]
    d = np.gradient(c, axis=0); d /= np.linalg.norm(d, axis=1, keepdims=True)
    nrm = np.c_[-d[:, 1], d[:, 0]]
    w = (w0 + (w1-w0)*t)
    return np.vstack([c + nrm*w, (c - nrm*w)[::-1]])

# name, note, left ctrl, right ctrl (None = mirror), appendages
DESIGNS = [
 ("Pelican", "belly + two return arms - the vessel that feeds itself",
  [(0,102),(18,138),(62,155),(150,152),(228,128),(285,96),(330,70),(380,52),
   (440,44),(510,40),(570,36),(620,26),(660,10)], None,
  [("arm", [(-138,132),(-150,330),(-44,462)], 12, 8),
   ("arm", [( 138,132),( 150,330),( 44,462)], 12, 8)]),

 ("Retort", "swan neck sweeping off one side - a front and a back",
  [(0,96),(18,132),(60,152),(140,150),(215,124),(268,90),(305,62),(340,40),
   (380,26),(430,16),(500,10),(580,7),(660,5)],
  [(0,96),(18,134),(62,155),(150,155),(232,142),(292,124),(340,110),(386,98),
   (432,84),(486,64),(540,44),(600,24),(660,8)], []),

 ("Alembic + beak", "still-head with a condensing spout, drips to a receiver",
  [(0,98),(16,126),(36,152),(126,155),(238,148),(298,120),(356,78),(414,58),
   (474,56),(514,76),(556,84),(602,54),(660,18)], None,
  [("beak", [(70,556),(140,512),(150,392)], 10, 7)]),

 ("Tribikos", "three-armed still - one belly, three condensing limbs",
  [(0,104),(20,140),(66,155),(154,150),(230,120),(280,86),(320,64),(370,50),
   (430,42),(500,36),(570,30),(620,22),(660,10)], None,
  [("arm", [(-46,442),(-146,330),(-128,150)], 11, 8),
   ("arm", [( 46,442),( 146,330),( 128,150)], 11, 8)]),

 ("Kerotakis", "sealed tower with a condensing shelf - the reflux vessel",
  [(0,150),(26,155),(120,155),(200,150),(236,120),(262,124),(300,118),(330,96),
   (360,100),(392,92),(420,72),(452,76),(486,66),(520,50),(570,34),(620,18),(660,6)],
  None, []),

 ("Aludel", "sublimation pots nested mouth-to-mouth, each inside the last",
  [(0,155),(60,152),(96,120),(104,138),(150,134),(196,104),(206,120),(250,116),
   (296,88),(306,102),(348,98),(392,72),(400,86),(440,82),(482,58),(490,70),
   (528,66),(568,44),(610,22),(660,6)], None, []),

 ("Athanor", "furnace tower - heavy base, tapering stack, chimney crown",
  [(0,155),(90,155),(120,148),(140,120),(150,132),(230,126),(252,100),(262,110),
   (340,104),(360,80),(370,88),(446,82),(464,60),(472,66),(544,60),(560,40),
   (600,30),(630,34),(660,26)], None, []),

 ("Matrass", "bolt-head: one true sphere, one long neck. the purest still",
  [(0,60),(14,104),(40,140),(84,155),(150,152),(206,132),(250,102),(282,72),
   (306,52),(330,42),(380,38),(450,36),(530,34),(600,32),(660,30)], None, []),

 ("Cucurbit + cap", "squat gourd, domed cap, low and wide - most volume",
  [(0,120),(20,150),(70,155),(170,153),(250,140),(320,116),(374,88),(410,70),
   (426,86),(470,84),(514,72),(560,54),(604,34),(636,16),(660,6)], None, []),

 ("Ampulla", "double belly with a pinched waist - two chambers, one throat",
  [(0,96),(18,134),(66,155),(146,150),(222,116),(268,78),(300,62),(336,74),
   (382,98),(424,112),(462,110),(500,92),(540,66),(584,42),(624,20),(660,6)],
  None, []),
]

def nested(zz, rr):
    """Three cascade cones, clipped to whatever the profile allows."""
    out = []
    for (f0, f1, s0, s1) in ((0.03,0.40,0.94,0.34),(0.32,0.66,0.68,0.25),
                             (0.60,0.92,0.46,0.15)):
        z0, z1 = f0*H, f1*H
        pts, n = [], 11
        for i in range(n):
            za = z1 - (z1-z0)*i/n; zb = z1 - (z1-z0)*(i+1)/n
            ra = (s1 + (s0-s1)*(i/n)) * RMAX
            rb = (s1 + (s0-s1)*((i+1)/n)) * RMAX
            lim_a = 0.88*np.interp(za, zz, rr); lim_b = 0.88*np.interp(zb, zz, rr)
            pts += [(za, min(ra, lim_a)), (zb, min(ra, lim_b)), (zb, min(rb, lim_b))]
        out.append(pts)
    return out
