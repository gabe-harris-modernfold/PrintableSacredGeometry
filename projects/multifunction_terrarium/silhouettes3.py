"""Round three: seedpod and spherical. Nested cascade tiers, 320 footprint.

Each is annotated with its own drop budget, because form and physics trade directly
here: detachment stages = usable cascade height / STEP, and no amount of nesting
buys extra stages -- water only falls as far as the object is tall. Elongated pods
keep the drop count; true spheres roughly halve it."""

import numpy as np, math
from scipy.interpolate import PchipInterpolator
import params as P

RMAX = 155.0
def prof(ctrl, n=320):
    z = np.array([c[0] for c in ctrl], float); r = np.array([c[1] for c in ctrl], float)
    zz = np.linspace(z.min(), z.max(), n)
    return zz, np.clip(PchipInterpolator(z, r)(zz), 0.5, None)

def tube(pts, w0, w1, n=50):
    p = np.array(pts, float); t = np.linspace(0, 1, n)[:, None]
    c = (1-t)**2*p[0] + 2*(1-t)*t*p[1] + t**2*p[2]
    d = np.gradient(c, axis=0); d /= np.linalg.norm(d, axis=1, keepdims=True)
    nn = np.c_[-d[:,1], d[:,0]]
    return np.vstack([c + nn*(w0+(w1-w0)*t), (c - nn*(w0+(w1-w0)*t))[::-1]])

# name, note, ctrl, appendages, (cascade z_lo, z_hi fraction of H)
DESIGNS = [
 ("Poppy head", "sphere under a fluted crown, on a short foot",
  [(0,46),(18,62),(38,50),(66,58),(100,96),(142,132),(192,152),(246,155),(298,144),
   (340,118),(368,90),(386,68),(398,86),(410,80),(422,56),(430,28)], [], (.10,.86)),

 ("Sea urchin", "oblate test, flattened poles - widest for its height",
  [(0,58),(14,108),(38,140),(78,153),(150,155),(220,150),(258,133),(284,102),
   (298,58),(304,26)], [], (.08,.90)),

 ("Acorn", "nut with a cupule cap - the cap lifts off to prune",
  [(0,28),(24,78),(58,116),(108,141),(168,153),(238,155),(296,150),(330,140),
   (344,152),(392,146),(432,116),(462,72),(478,32)], [], (.06,.72)),

 ("Banksia", "elongated follicle cone - keeps almost all the drop count",
  [(0,52),(24,94),(60,124),(110,146),(180,155),(280,155),(380,148),(458,130),
   (518,104),(558,70),(582,36),(592,14)], [], (.06,.92)),

 ("Lotus pod", "inverted cone, flat perforated face - you look down into it",
  [(0,28),(38,34),(68,42),(98,60),(158,94),(228,127),(298,149),(346,155),
   (372,155),(380,148)], [], (.14,.96)),

 ("Gumnut", "squat cup with a heavy rim - lowest and widest",
  [(0,68),(18,108),(48,137),(98,152),(168,151),(216,143),(246,140),(258,155),
   (272,152),(282,138)], [], (.10,.88)),

 ("Pinecone", "bi-pointed ovoid of overlapping scales - ports become scales",
  [(0,22),(28,72),(68,113),(128,144),(198,155),(298,153),(398,140),(468,117),
   (528,87),(568,50),(596,18),(606,6)], [], (.05,.93)),

 ("Nigella", "inflated pod with radiating horns - horns are the return tubes",
  [(0,48),(24,94),(58,127),(108,150),(178,155),(258,152),(318,137),(368,111),
   (408,78),(436,44),(446,20)], 
  [("h",[(-30,430),(-108,530),(-138,612)],7,3),("h",[(30,430),(108,530),(138,612)],7,3),
   ("h",[(-14,442),(-46,556),(-58,640)],7,3),("h",[(14,442),(46,556),(58,640)],7,3),
   ("h",[(0,446),(0,540),(0,650)],7,3)], (.06,.86)),

 ("Dehiscing pod", "split along three seams - the split IS the viewing port",
  [(0,34),(26,84),(64,122),(118,146),(190,155),(290,154),(384,143),(452,122),
   (508,92),(546,56),(566,24),(574,8)], [], (.05,.90)),

 ("Geode", "sphere with a face cut away, terraced shells stepping inward",
  [(0,52),(16,102),(42,136),(84,152),(155,155),(226,151),(268,135),(298,106),
   (318,66),(328,28)], [], (.07,.90)),
]

def nested(zz, rr, frac, tiers=3):
    """Concentric cascade cones: inner one highest and narrowest, spilling outward."""
    H = zz.max(); z_lo, z_hi = frac[0]*H, frac[1]*H
    span = z_hi - z_lo
    out = []
    for i in range(tiers):
        a = i / tiers
        z0 = z_lo + span*a*0.62
        z1 = z_lo + span*(0.42 + 0.58*a)
        s0 = 0.94 - 0.26*i          # radius fraction at the wide (low) end
        s1 = 0.34 - 0.09*i
        pts, n = [], 11
        for k in range(n):
            za = z1 - (z1-z0)*k/n; zb = z1 - (z1-z0)*(k+1)/n
            ra = (s1 + (s0-s1)*(k/n)); rb = (s1 + (s0-s1)*((k+1)/n))
            la = 0.9*np.interp(za, zz, rr); lb = 0.9*np.interp(zb, zz, rr)
            pts += [(za, min(ra*RMAX, la)), (zb, min(ra*RMAX, lb)), (zb, min(rb*RMAX, lb))]
        out.append(pts)
    return out

def budget(zz, frac):
    """Stages and airborne drops for this profile's usable cascade height."""
    h = (frac[1]-frac[0]) * zz.max()
    stages = int(h / P.STEP)
    tf = math.sqrt(2*P.STEP*1e-3/P.GRAV)
    dps = P.Q_TRICKLE/3600*1e-3/(P.DROP_V*1e-9)
    return h, stages, dps*stages*tf
