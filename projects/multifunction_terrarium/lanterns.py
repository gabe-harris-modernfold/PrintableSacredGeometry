"""Round four: faceted gem-cut lanterns, per the reference.

Vertical prism of few large facets with narrow chamfers between, a crystal crown
chamfering in to a ridge/table/apex, a foot, and a waterline. The chevron groove
field lives on the facets -- see chevron.py for the wetting physics."""

import numpy as np, math
import params as P

RMAX, HB, HC, FOOT = 155.0, 400.0, 140.0, 26.0
H = HB + HC

def plan(n, alt=None, twist=0.0):
    """Plan polygon: n primary facets, optional narrow chamfers at radius alt*R."""
    if alt is None:
        a = np.arange(n) * 2*math.pi/n
        return a, np.full(n, RMAX)
    a = np.concatenate([[k*2*math.pi/n, (k+0.5)*2*math.pi/n] for k in range(n)])
    r = np.array([RMAX, alt*RMAX] * n)
    return a, r

def body_r(z, taper=0.045, barrel=0.0):
    """Radius scale down the body: slight taper, optional barrel swell."""
    t = z / HB
    return 1.0 - taper*t + barrel*math.sin(math.pi*t)

def crown_r(z, kind, table):
    """Radius scale through the crown."""
    t = (z - HB) / HC
    if kind == "pyramid": return (1-t)
    if kind == "ridge":   return (1-t)          # radius collapses; ridge handled in x
    if kind == "table":   return 1 - (1-table)*t
    if kind == "dome":    return math.sqrt(max(0.0, 1 - t*t)) * (1-table) + table*(1-t)
    return 1-t

VARIANTS = [
 ("Hex, ridge crown",        6, None, "ridge",   0.00, 0.00, 0.0),
 ("Hex + chamfers  (12)",    6, 0.92, "ridge",   0.00, 0.00, 0.0),
 ("Octagon, table crown",    8, None, "table",   0.52, 0.00, 0.0),
 ("Octagon + chamfers (16)", 8, 0.94, "pyramid", 0.00, 0.00, 0.0),
 ("Decagon, ridge",         10, None, "ridge",   0.00, 0.00, 0.0),
 ("Dodecagon, dome",        12, None, "dome",    0.30, 0.00, 0.0),
 ("Emerald cut  (4+4)",      4, 0.78, "table",   0.58, 0.00, 0.0),
 ("Trefoil  (3+3)",          3, 0.86, "pyramid", 0.00, 0.00, 0.0),
 ("Barrel, twisted facets",  8, 0.95, "ridge",   0.00, 0.13, 22.0),
 ("Phyllotactic  (11+11)",  11, 0.96, "dome",    0.26, 0.00, 0.0),
]

def stats():
    casc = HB - FOOT - 40           # reservoir below, crown above
    stages = int(casc / P.STEP)
    tf = math.sqrt(2*P.STEP*1e-3/P.GRAV)
    dps = P.Q_TRICKLE/3600*1e-3/(P.DROP_V*1e-9)
    return casc, stages, dps*stages*tf
