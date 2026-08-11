#!/usr/bin/env python3
"""Mirror ray tracing shared by the optical models.

facets() builds inward-oriented triangle normals using RADIAL terms only -- a
n.centroid test lets the z term dominate and flips orientation inconsistently.

shoot() appends a display stub past the last hit; never reverse a trace from its
end. Use reverse_from_last(): too short a stub and the reverse starts inside the
geometry, too long and (facets being two-sided here) the returning ray strikes the
outside first. Both break time reversal silently."""

import math

import numpy as np


def facets(ring_list):
    """Antiprism-topology triangle strip per band: an exact tiling, no louvres.
    Normals are oriented inward using the RADIAL terms only -- testing n.centroid
    lets the z term dominate high up the stack and flips bands inconsistently."""
    tris = []
    N = len(ring_list[0])
    for U, L in zip(ring_list[:-1], ring_list[1:]):        # U above, L below
        for m in range(N):
            m2 = (m + 1) % N
            tris.append([U[m], U[m2], L[m]])
            tris.append([L[m], L[m2], U[m2]])
    t = np.asarray(tris, float)
    n = np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0])
    n /= np.linalg.norm(n, axis=1, keepdims=True)
    c = t.mean(axis=1)
    n[np.einsum("ij,ij->i", n[:, :2], c[:, :2]) > 0] *= -1.0
    return t, n

def hit(p, d, tris, eps=1e-7):
    """Nearest forward ray/triangle intersection (Moller-Trumbore, vectorised)."""
    e1 = tris[:, 1] - tris[:, 0]
    e2 = tris[:, 2] - tris[:, 0]
    h = np.cross(d, e2)
    a = np.einsum("ij,ij->i", e1, h)
    ok = np.abs(a) > 1e-12
    f = np.where(ok, 1.0 / np.where(ok, a, 1.0), 0.0)
    s = p - tris[:, 0]
    u = f * np.einsum("ij,ij->i", s, h)
    q = np.cross(s, e1)
    v = f * (q @ d)
    t = f * np.einsum("ij,ij->i", e2, q)
    good = ok & (u >= -1e-9) & (v >= -1e-9) & (u + v <= 1 + 1e-9) & (t > eps)
    if not good.any():
        return None, None
    idx = np.where(good)[0]
    k = idx[np.argmin(t[idx])]
    return k, float(t[k])

def Lz(p, d):
    """Axial angular momentum r*v_phi = (p x d).z -- the invariant we must break."""
    return float(p[0] * d[1] - p[1] * d[0])

def shoot(p, d, tris, nrm, N, max_bounce=120):
    """Trace outward from the focus. Returns polyline + per-bounce log."""
    p = np.asarray(p, float)
    d = np.asarray(d, float)
    d = d / np.linalg.norm(d)
    pts, log = [p.copy()], []
    for _ in range(max_bounce):
        k, t = hit(p, d, tris)
        if k is None:
            break
        p = p + t * d
        n = nrm[k]
        d = d - 2.0 * (d @ n) * n
        pts.append(p.copy())
        r = math.hypot(p[0], p[1])
        log.append(dict(facet=int(k), band=int(k) // (2 * N), r=r, z=float(p[2]),
                        L=Lz(p, d), vphi=(Lz(p, d) / r if r > 1e-9 else 0.0)))
    # Display stub only. Do NOT reverse a trace from its end: too short and the
    # reverse starts inside the geometry, too long and -- because facets here are
    # two-sided -- the returning ray strikes the OUTSIDE of the tube first. Either
    # way time reversal breaks silently. Reverse from the last hit instead
    # (see entry_rays).
    pts.append(p + d * 45.0)
    return np.asarray(pts), log

def winding(pts):
    """(total turns swept, azimuthal direction reversals). Summing |da| rather
    than taking the net keeps a ray that winds +2 then -2 turns from scoring 0,
    and the reversal count is what separates a vortex from a wobble."""
    a = np.unwrap(np.arctan2(pts[:, 1], pts[:, 0]))
    da = np.diff(a)
    big = np.abs(da) > 1e-9
    sg = np.sign(da[big])
    return float(np.abs(da).sum() / (2.0 * np.pi)), int((np.diff(sg) != 0).sum())

def launch_from_focus(theta, az=0.0, z_f=0.0):
    """Start ON the axis at the base centre, so L = 0 and the focus is exact.
    Only ONE degree of freedom: at r = 0 every horizontal direction is
    equivalent, so an azimuthal-vs-radial mix angle would be redundant with az.
    All the swirl the ray acquires comes from the faceting, not the launch."""
    a = math.radians(az)
    th = math.radians(theta)
    return (np.array([0.0, 0.0, z_f]),
            np.array([math.sin(th) * math.cos(a), math.sin(th) * math.sin(a),
                      math.cos(th)]))

def reverse_from_last(pts, log, nudge=0.5):
    """Exact inbound ray: the outbound final leg traversed backwards. Starts just
    off the last hit facet so it retraces bounce for bounce."""
    p_last = pts[-2]
    d_out = pts[-1] - pts[-2]
    d_out = d_out / np.linalg.norm(d_out)
    return p_last + d_out * nudge, -d_out

def tilt(nrm, sigma_deg, rng):
    """Rotate each normal by a random small angle about a random in-plane axis."""
    if sigma_deg <= 0:
        return nrm.copy()
    n = nrm / np.linalg.norm(nrm, axis=1, keepdims=True)
    a = rng.normal(size=n.shape)
    a -= n * np.einsum("ij,ij->i", a, n)[:, None]          # project into the plane
    a /= np.maximum(np.linalg.norm(a, axis=1, keepdims=True), 1e-12)
    eps = np.radians(rng.normal(0.0, sigma_deg, len(n)))[:, None]
    p = n * np.cos(eps) + a * np.sin(eps)
    return p / np.linalg.norm(p, axis=1, keepdims=True)

def miss(p, d):
    """Distance from the origin (the design focus) to the ray p + t d, t >= 0."""
    t = float(-(p @ d))
    if t < 0:
        return float(np.linalg.norm(p))
    return float(np.linalg.norm(p + t * d))
