#!/usr/bin/env python3
"""
Three whole tetrahelices laid up as a three-ply rope, cut into a module that
STACKS: set the next one on top, give it a twist, and the braid carries on.

YOU WERE RIGHT ABOUT THE STACKING, AND IT COST NOTHING
------------------------------------------------------
A straight tetrahelix cannot repeat. Gluing regular tetrahedra face-to-face forces
a twist of exactly theta = arccos(-2/3) = 131.8103149 deg per cell, and theta/360 =
0.366139764... is irrational, so the column never comes back to the same azimuth --
that is the whole content of Fuller's 5.69 deg deficit, and it is why tetrahelix.py
has to fudge theta to 132 deg to get a stackable module out of it.

Laying three of them up as rope hands back a free parameter, and it lands exactly
on this problem. Carry each strand along a lay helix and the lay's TORSION rotates
the strand's cross-section as it goes. In the frame that follows the lay, the twist
per cell is no longer theta but

    psi = theta - tau*H          tau = lay torsion, H = rise per cell

and the module repeats after P cells as soon as psi is a rational part of a turn:

    P * psi = 360 * m        ->        psi = 360 m / P

So the lay's torsion absorbs the irrationality that the tetrahedron refuses to give
up. No cell is compromised to get it: theta stays exactly arccos(-2/3), and the
correction is spent entirely on the lay, which is a knob a single tetrahelix does
not have.

Which P are available? Only those where 360m/P is close enough to theta that the
lay can actually supply the difference (a unit-speed helix of radius A has torsion
at most 1/2A). Those are the continued-fraction convergents of theta/360:

    m/P     psi           tau*H needed     what it is
    1/3     120.000000    +11.8103 deg     tight rope. THE DEFAULT.
    3/8     135.000000     -3.1897 deg     gentle, opposite-handed lay
    4/11    130.909091     +0.9012 deg     nearly straight
    11/30   132.000000     -0.1897 deg     <- tetrahelix.py's 132 deg ten-cycle,
                                              recovered as the P=30 convergent
    15/41   131.707317     +0.1030 deg     nearly straight, other hand

The P=30 row is the check that this is the right frame for the problem: Fuller's
ten-cycle is not a special trick, it is one convergent among many, the one you are
forced into when the lay is not available to help.

FACE ON FACE, STILL
-------------------
Face-sharing is combinatorial, not metric: cell k is {k,k+1,k+2,k+3} and cell k+1 is
{k+1,k+2,k+3,k+4}, so they share the whole triangle {k+1,k+2,k+3} whatever the
vertices do. Bending a strand moves the vertices but cannot separate those cells --
_verify() asserts every consecutive pair still shares three vertices and that the
shared triangle is identical from both sides. What bending costs is regularity, not
adjacency: the faces stop being equilateral. That is measured and printed, not
waved away.

(One correction, since it matters for the design: the tetrahedron is not the only
shape that shares faces -- cubes do, octahedra and tetrahedra do together,
triangular prisms do. What is special is the opposite. The regular tetrahedron is
the one that CANNOT stack periodically: it does not tile space, and its face-to-face
chain has that irrational twist. The braid is a way of getting periodicity back
out of a shape that has none.)

HOW IT STACKS
-------------
The module is the infinite braid cut by two HORIZONTAL planes exactly one module
apart. Horizontal planes are preserved by the screw, so the cut piece tiles the
column: rotate by the module's lay advance, raise by its height, repeat. The cut
runs straight through beams, so both faces are flat sets of beam cross-sections --
flat on the bed for module one, flat glue faces at every seam, and once assembled
every tetrahedron is whole and every consecutive pair shares its face.

    ASSEMBLY at the defaults: each module 219 mm tall, 391 mm2 of flat mating
    face, and each one goes on +90.84 deg from the one below (the module's lay
    advance is 210.84 deg, and the braid is 3-fold, so 90.84 is the same thing).

That the faces really mate is asserted, not eyeballed: the top cut's triangles,
screwed back by one module, land on the bottom cut's to 1e-7 mm. Getting there needed
one non-obvious thing -- see Rope.basis. The vertices being screw-periodic is not
enough. Anything with an ORIENTATION (a beam's flats, a node ball's axes) has to be
built in a frame that commutes with the screw, and the rotation-minimising frame
that places the vertices does not: it is unwound by tau*s, so it accumulates an
extra tau*P*H = 35.43 deg per period. Orient solids with it and the vertices tile
while the solids do not, which surfaces as cut faces 1.16 mm out of register.

PHASE
-----
Because psi is rational the cells do not creep around the strand -- they land on
exactly P azimuthal rails. So there is a free choice the straight tetrahelix never
offers: roll all three strands about their own axes together, and aim the rails
between the neighbouring plies rather than at them. It is worth 3.6 mm of ply gap
unrolled against 12.8 mm at the 87.85 deg the search picks, and it costs nothing --
a common phase preserves the 3-fold symmetry and cannot touch k*psi. See best_phase.

Run:  python braided_tetrahelix.py [--edge 22] [--steps 39] [--period 3]
                                   [--clear 5] [--left] [--hard-lay]
                                   [--stack 2] [--base] [--weld]
"""

import argparse
import importlib.util
import math
import os

import numpy as np

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))

import mesh_kit as MK
import tetrahelix as TH

EDGE = 22.0             # tetrahedron edge of each strand. Big cells on purpose:
                        # they have to read as tetrahedra or the braid disappears.
STEPS = 39              # cells per strand per module; rounded up to a multiple of
                        # PERIOD. 39 = 13 periods -> a 61 x 61 x 219 mm module.
PERIOD = 3              # cells per screw period. See the convergent table above.
CLEAR = 5.0             # mm between strand vertex cylinders. 0 packs the plies into
                        # one clogged mass; 5 opens daylight between them.
BEAM_F = 0.145          # beam circumradius / edge. Fatter than tetrahelix.py's
NODE_F = 0.200          # 0.110 -- a see-through strand does not read as a ply.
                        # Node must exceed beam or the six beams there only kiss.
TIE_F = 1.05            # tie beam radius / beam radius
TIE_MAX = 0.85          # a mutual-nearest pair further apart than this * edge is
                        # not a contact, it is two strands passing
RING_W = 8.0            # optional base ring bar width, mm
PAD = 8                 # cells built beyond each cut plane before clipping. Must
                        # exceed the reach of anything that can cross a plane; the
                        # margin is ~34 mm against a worst case of ~12 mm, and
                        # _verify() rebuilds at PAD+4 and demands the same mesh.
EPS = 1e-9


def clip_convex(verts, faces, nrm, off):
    """Clip a CONVEX solid to the half-space n.p >= off. None if nothing survives.

    Convexity is what makes this short: the cut is a single convex polygon, so the
    cap needs no hole-finding -- collect the on-plane points, order them by angle
    about their centroid, fan them. Every solid here is convex (triangular-prism
    beams, octahedral nodes, rectangular base bars), and the clipper is checked
    against known volumes in _verify(): a split solid's two halves must sum to the
    whole, exactly."""
    V = np.asarray(verts, float)
    nrm = np.asarray(nrm, float)
    d = V @ nrm - off
    if (d >= -EPS).all():
        return [np.asarray(p, float) for p in verts], list(faces)
    if (d <= EPS).all():
        return None

    out_v, out_f, on = [], [], []

    def push(p):
        out_v.append(p)
        return len(out_v) - 1

    for tri in faces:
        poly = []
        n = len(tri)
        for i in range(n):
            a, b = int(tri[i]), int(tri[(i + 1) % n])
            da, db = d[a], d[b]
            if da >= -EPS:
                poly.append(V[a])
            if (da > EPS and db < -EPS) or (da < -EPS and db > EPS):
                poly.append(V[a] + (da / (da - db)) * (V[b] - V[a]))
        if len(poly) < 3:
            continue
        idx = [push(p) for p in poly]
        for i in range(1, len(poly) - 1):
            out_f.append((idx[0], idx[i], idx[i + 1]))
        on += [p for p in poly if abs(float(p @ nrm) - off) <= 1e-7]

    if len(on) >= 3:
        P = np.array(on)
        _, keep = np.unique(np.round(P, 6), axis=0, return_index=True)
        P = P[np.sort(keep)]
        if len(P) >= 3:
            c = P.mean(axis=0)
            u = P[0] - c
            if np.linalg.norm(u) > 1e-12:
                u = u / np.linalg.norm(u)
                w = np.cross(nrm, u)
                P = P[np.argsort(np.arctan2((P - c) @ w, (P - c) @ u))]
                ci = push(c)
                ring = [push(p) for p in P]
                for i in range(len(ring)):
                    j = (i + 1) % len(ring)
                    flip = np.cross(P[i] - c, P[j] - c) @ nrm > 0
                    out_f.append((ci, ring[j], ring[i]) if flip
                                 else (ci, ring[i], ring[j]))
    return out_v, out_f


def feasible_periods(edge, clear, pmax=41):
    """(P, m, psi_deg, tauH_deg, ok) for every period, ok where the lay can supply
    the torsion. Used for the report and for the error message on a bad --period."""
    theta = math.degrees(math.acos(-2.0 / 3.0))
    R = 3.0 * math.sqrt(3.0) / 10.0 * edge
    H = edge / math.sqrt(10.0)
    A = (2.0 * R + clear) / math.sqrt(3.0)
    tmax = math.degrees(H / (2.0 * A))          # tau*H ceiling, in degrees
    out = []
    for P in range(2, pmax + 1):
        m = max(1, round(P * theta / 360.0))
        psi = 360.0 * m / P
        out.append((P, m, psi, theta - psi, abs(theta - psi) <= tmax))
    return out, tmax


def solve_lay(period, edge, clear, hard=False):
    """Lay curvature kappa that makes the module repeat after `period` cells.

    tau = kappa*sqrt(1 - A^2 kappa^2) is quadratic in kappa^2, so a needed torsion
    has TWO lays that deliver it -- a gentle one and a tight one either side of the
    torsion maximum at A*kappa = 1/sqrt2. Default to the gentle root: same
    periodicity, less bend, rounder cells. --hard-lay takes the other."""
    theta = math.acos(-2.0 / 3.0)
    R = 3.0 * math.sqrt(3.0) / 10.0 * edge
    H = edge / math.sqrt(10.0)
    A = (2.0 * R + clear) / math.sqrt(3.0)
    m = max(1, round(period * math.degrees(theta) / 360.0))
    psi = 2.0 * math.pi * m / period
    tau = (theta - psi) / H
    if abs(tau) > 1.0 / (2.0 * A):
        good = [str(p) for p, _m, _p, _t, ok in feasible_periods(edge, clear)[0] if ok]
        raise SystemExit(
            f"--period {period} is not reachable: it needs "
            f"{math.degrees(theta - psi):+.4f} deg of torsion per cell but a lay of "
            f"radius {A:.2f} mm can only supply {math.degrees(H/(2*A)):.4f}. "
            f"Reachable periods at this edge/clear: {', '.join(good)}")
    disc = math.sqrt(max(0.0, 1.0 - 4.0 * A * A * tau * tau))
    k2 = (1.0 + disc if hard else 1.0 - disc) / (2.0 * A * A)
    kappa = math.copysign(math.sqrt(k2), tau) if tau else 0.0
    return kappa, m, psi, A


def octa_at(c, r, ax):
    """TH.octa, but oriented in a given right-handed triad instead of world axes.

    A world-axis octahedron is not invariant under the module's screw rotation
    (only under multiples of 90 deg), so node balls built with TH.octa break the
    stacking. Same 8 faces and winding, local axes."""
    e1, e2, t = ax
    c = np.asarray(c, float)
    verts = [c + r * d for d in (e1, e2, -e1, -e2, t, -t)]
    faces = [(0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4),
             (1, 0, 5), (2, 1, 5), (3, 2, 5), (0, 3, 5)]
    return verts, faces


def ply_gap(rope):
    """Closest approach between the vertex sets of two neighbouring plies."""
    return min(float(np.linalg.norm(rope.verts[j][:, None]
                                    - rope.verts[(j + 1) % 3][None], axis=2).min())
               for j in range(3))


def best_phase(edge, steps, period, clear, left=False, hard=False):
    """Roll all three strands about their own axes to open the most daylight.

    This knob only exists because psi is rational. With psi = 360m/P a strand's
    vertices do not creep around it -- they land on exactly P azimuthal rails -- so
    whether a rail points AT the neighbouring ply or between two of them is a real
    choice, and on the defaults it is worth 3.6 mm of ply gap at phase 0 against
    12.7 mm at 32 deg. Nothing else moves: a phase common to all three strands keeps
    the 3-fold symmetry, and it cannot touch the periodicity, which depends only on
    k*psi.

    The gap has period 360/P in phase, so a 1 deg sweep of one period plus a 0.05 deg
    refinement is exhaustive and costs a few hundred frame evaluations."""
    span = 360.0 / period

    def g(ph):
        return ply_gap(Rope(edge, 3 * period, period, clear, left, hard,
                            phase=float(ph), pad=2))

    coarse = max(np.arange(0.0, span, 1.0), key=g)
    return math.radians(float(max(np.arange(coarse - 1.0, coarse + 1.05, 0.05), key=g)))


class Lay:
    """The rope lay: three strand axes 120 deg apart on one unit-speed helix.

    Unit speed matters. Parameterised by z instead of arc length the strands would
    be stretched axially as well as bent, and the cells would grow taller towards
    the outside of the lay on top of everything else."""

    def __init__(self, A, kappa):
        self.A, self.k = A, kappa
        q = (A * kappa) ** 2
        assert q < 1.0, "lay too tight: A*kappa must be < 1"
        self.b = math.sqrt(1.0 - q)          # rise per unit arc length
        self.curv = abs(A * kappa * kappa)   # curvature of a strand axis
        self.tau = self.b * kappa            # torsion -- the whole point, see above

    def frame(self, j, s):
        """Axis point, tangent, and a ROTATION-MINIMISING pair of normals.

        Frenet would spin about the tangent at the torsion rate and add tau*H to
        every cell's twist. Here that spin is the mechanism, not an accident, so it
        has to be applied deliberately and once: the frame is unwound by tau*s, and
        the strand's own theta*k is laid off in it, which puts the Frenet-frame
        azimuth at exactly k*psi -- the quantity the whole periodicity argument is
        about. Using Frenet raw instead applies the torsion a second time and drives
        the step families apart by several percent in opposite directions, a bias
        rather than a spread. _verify() measures it, and report() prints it, so the
        difference cannot pass unnoticed."""
        a = self.k * s + 2.0 * math.pi * j / 3.0
        ca, sa = math.cos(a), math.sin(a)
        C = np.array([self.A * ca, self.A * sa, self.b * s])
        T = np.array([-self.A * self.k * sa, self.A * self.k * ca, self.b])
        e1, e2 = np.array([ca, sa, 0.0]), None
        e2 = np.cross(T, e1)
        w = -self.tau * s
        cw, sw = math.cos(w), math.sin(w)
        return C, T, cw * e1 + sw * e2, -sw * e1 + cw * e2


class Rope:
    """Three bent tetrahelices on a periodic lay, plus the ties between them."""

    def __init__(self, edge=EDGE, steps=STEPS, period=PERIOD, clear=CLEAR,
                 left=False, hard=False, phase=None, pad=PAD):
        kappa, m, psi, A = solve_lay(period, edge, clear, hard)
        self.theta, self.R, self.H = TH.geometry(edge, closed=False)
        if left:                                    # true mirror: both hands flip
            self.theta, kappa = -self.theta, -kappa
        self.edge, self.steps, self.period, self.clear = edge, steps, period, clear
        self.m, self.psi, self.left, self.hard = m, psi, left, hard
        self.lay = Lay(A, kappa)
        self.beam_r, self.node_r = BEAM_F * edge, NODE_F * edge
        self.arc = steps * self.H                   # arc of one module's axis
        self.height = self.lay.b * self.arc
        self.radius = A + self.R + self.node_r
        self.lay_turns = abs(kappa) * self.arc / (2.0 * math.pi)

        self.phase = (best_phase(edge, steps, period, clear, left, hard)
                      if phase is None else math.radians(phase))
        self.ks = np.arange(-pad, steps + pad + 1)
        self.verts = np.stack([self._strand(j) for j in range(3)])   # (3, K, 3)
        self.beams = [(i, i + d) for d in (1, 2, 3)
                      for i in range(len(self.ks) - d)]
        self.cells = [(i, i + 1, i + 2, i + 3) for i in range(len(self.ks) - 3)]
        self.ties = self._ties()

    # --- indexing helpers: self.verts is 0-based over self.ks -------------------
    def at(self, j, k):
        return self.verts[j][int(k) - int(self.ks[0])]

    def _strand(self, j):
        out = np.empty((len(self.ks), 3))
        for i, k in enumerate(self.ks):
            C, _T, e1, e2 = self.lay.frame(j, k * self.H)
            a = k * self.theta + self.phase
            out[i] = C + self.R * (math.cos(a) * e1 + math.sin(a) * e2)
        return out

    def basis(self, j, s):
        """A SCREW-COVARIANT triad at arc s: (outward radial, binormal, tangent).

        This is the lay's raw Frenet triad, deliberately NOT the rotation-minimising
        one. The RMF is right for placing vertices -- it is what makes psi come out
        at theta - tau*H -- but it is unwound by tau*s, so it does not commute with
        the screw: advance by one period and the RMF has rolled an extra tau*P*H
        (35.43 deg on the defaults). Orient a beam's flats or a node ball's axes with
        it and the SOLIDS stop being periodic even though their vertices still are,
        which shows up as cut faces that will not mate -- the module stops stacking.
        The Frenet triad rotates with the screw by construction, so use it for
        anything that has an orientation rather than just a position."""
        a = self.lay.k * s + 2.0 * math.pi * j / 3.0
        ca, sa = math.cos(a), math.sin(a)
        e1 = np.array([ca, sa, 0.0])
        T = np.array([-self.lay.A * self.lay.k * sa,
                      self.lay.A * self.lay.k * ca, self.lay.b])
        T = T / np.linalg.norm(T)
        return e1, np.cross(T, e1), T

    def local_out(self, j, s):
        """Direction to aim the beam flats, so a strand reads as faceted stock."""
        return self.basis(j, s)[0]

    def screw(self, n_cells):
        """The transform that maps the braid onto itself, n_cells further along."""
        return self.lay.k * n_cells * self.H, self.lay.b * n_cells * self.H

    def screw_matrix(self, n_cells):
        a, dz = self.screw(n_cells)
        c, s = math.cos(a), math.sin(a)
        return (np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]]),
                np.array([0.0, 0.0, dz]))

    def _ties(self):
        """One tie per mutual-nearest vertex pair, made exactly PERIOD-periodic.

        Mutual-nearest, not nearest: one-sided nearest fans a run of ties onto the
        same vertex wherever two strands are close over a stretch. And the pattern
        is generated from one period and replicated by the screw rather than taken
        as found, because a tie set that is only approximately periodic puts
        different ties either side of a seam, and then stacked modules do not agree
        about where the joint is."""
        base = []
        for j in range(3):
            P, Q = self.verts[j], self.verts[(j + 1) % 3]
            D = np.linalg.norm(P[:, None, :] - Q[None, :, :], axis=2)
            a, b = D.argmin(axis=1), D.argmin(axis=0)
            for i in range(len(P)):
                k = int(self.ks[i])
                if k < 0 or k >= self.period:
                    continue
                l = int(a[i])
                if int(b[l]) == i and D[i, l] < TIE_MAX * self.edge:
                    base.append((j, k, (j + 1) % 3, int(self.ks[l])))
        out = []
        lo, hi = int(self.ks[0]), int(self.ks[-1])
        span = (hi - lo) // self.period + 2
        for t in range(-span, span + 1):
            for j, k, j2, l in base:
                kk, ll = k + t * self.period, l + t * self.period
                if lo <= kk <= hi and lo <= ll <= hi:
                    d = float(np.linalg.norm(self.at(j2, ll) - self.at(j, kk)))
                    out.append((j, kk, j2, ll, d))
        return out

    def beam_lengths(self):
        """Measured beam lengths per step family, over one module only."""
        res = {}
        for d in (1, 2, 3):
            q = [np.linalg.norm(self.at(j, k + d) - self.at(j, k))
                 for j in range(3) for k in range(self.steps)]
            res[d] = np.array(q)
        return res

    def cell_volumes(self):
        v = np.array([[self.at(j, k + i) for i in range(4)]
                      for j in range(3) for k in range(self.steps)])
        return np.abs(np.einsum("ij,ij->i", v[:, 1] - v[:, 0],
                                np.cross(v[:, 2] - v[:, 0],
                                         v[:, 3] - v[:, 0]))) / 6.0


def ring_height(rope):
    """Tall enough to swallow a leg end -- see base_ring()."""
    return 2.0 * rope.beam_r + 1.2


def base_ring(rope, z_lo):
    """Optional open bar triangle on the bed, with legs up to the lowest vertices.

    Not needed for stacking -- the bottom cut is already flat -- but the module is
    3.5:1 on a 48 mm footprint, so it is there for anyone who would rather not
    trust a brim. Corners go at radius A+R, not at the strand axes: axis corners
    give a 28 mm triangle under a 63 mm object."""
    ring_h = ring_height(rope)
    rad = rope.lay.A + rope.R
    corners, legs = [], []
    for j in range(3):
        a = 2.0 * math.pi * j / 3.0
        corners.append(np.array([rad * math.cos(a), rad * math.sin(a),
                                 z_lo + 0.5 * ring_h]))
        # legs must land INSIDE the module, or they reach for geometry the slab
        # clip has already removed
        inside = [rope.at(j, k) for k in range(0, 8) if rope.at(j, k)[2] > 0.0]
        low = sorted(inside, key=lambda p: p[2])[:2]
        legs += [(corners[-1], p) for p in low]
    bars = [(corners[i], corners[(i + 1) % 3]) for i in range(3)]
    return corners, bars, legs, ring_h


def _solids(rope, base):
    """Every solid of the infinite braid that could touch the module slab, untrimmed.

    Built over the padded cell range and handed to the clipper; whatever falls
    outside the slab is dropped there. No seam bookkeeping: the module is just a
    slice of one continuous object."""
    out = []
    for j in range(3):
        P = rope.verts[j]
        for a, b in rope.beams:
            s_mid = 0.5 * (rope.ks[a] + rope.ks[b]) * rope.H
            out.append((TH.beam(P[a], P[b], rope.beam_r,
                                rope.local_out(j, s_mid)), f"strand{j}"))
        for i, k in enumerate(rope.ks):
            out.append((octa_at(P[i], rope.node_r, rope.basis(j, k * rope.H)),
                        f"node{j}"))
    for j, k, j2, l, _d in rope.ties:
        out.append((TH.beam(rope.at(j, k), rope.at(j2, l), rope.beam_r * TIE_F,
                            rope.local_out(j, k * rope.H)), "tie"))
    if base:
        corners, bars, legs, ring_h = base_ring(rope, 0.0)
        for p, q in bars:
            u = q[:2] - p[:2]
            u = u / np.linalg.norm(u)
            n = np.array([-u[1], u[0]]) * (RING_W / 2.0)
            out.append((MK.prism([p[:2] - n, q[:2] - n, q[:2] + n, p[:2] + n],
                                 -ring_h, 0.0), "base"))
        for p, q in legs:
            out.append((MK.tube(p - np.array([0.0, 0.0, ring_h]), q,
                                rope.beam_r, nseg=8), "base"))
        for p in corners:
            out.append((TH.octa(p - np.array([0.0, 0.0, ring_h]),
                                0.5 * ring_h), "base"))
    return out


def build(edge=EDGE, steps=STEPS, period=PERIOD, clear=CLEAR, left=False,
          hard=False, base=False, stack=1, phase=None, pad=PAD):
    if steps % period:
        steps += period - steps % period          # must be a whole number of periods
    rope = Rope(edge, steps, period, clear, left, hard, phase, pad)
    dz = rope.height

    m = MK.Mesh()
    cut = {"lo": 0, "hi": 0}
    kept = []
    for (verts, faces), tag in _solids(rope, base):
        z = np.asarray(verts, float)[:, 2]
        lo_cut, hi_cut = z.min() < 0.0 < z.max(), z.min() < dz < z.max()
        r = (verts, faces)
        if tag != "base":               # the base lives BELOW the slab on purpose,
            r = clip_convex(*r, (0.0, 0.0, 1.0), 0.0)   # so it is exempt from the
            if r is None:                               # bottom cut -- clipping it
                continue                                # deleted it outright
            cut["lo"] += bool(lo_cut)
        r = clip_convex(r[0], r[1], (0.0, 0.0, -1.0), -dz)
        if r is None:
            continue
        cut["hi"] += bool(hi_cut)
        kept.append((r, tag))
    # with a base the module no longer starts at z=0; drop it onto the bed
    drop = np.array([0.0, 0.0, ring_height(rope) if base else 0.0])

    for i in range(max(1, stack)):
        Rz, t = rope.screw_matrix(i * steps)
        for (verts, faces), tag in kept:
            V = np.asarray(verts, float) @ Rz.T + t + drop
            m.add_solid(V, faces, tag if i == 0 else f"{tag}#{i}")
    return m, rope, cut


class Welded:
    """Result of weld(): quacks like MK.Mesh for write_stl()."""

    def __init__(self, v, f):
        self._v, self._f = v, f

    def arrays(self):
        return self._v, self._f


def weld(m, log=lambda *a: print(*a, flush=True)):
    """CSG-union every shell into one 2-manifold solid, or None if unavailable.

    Slow: hundreds of pairwise unions on a growing accumulator, tens of minutes.
    Slicers union the raw soup correctly, so this is only for an STL that has to
    pass a strict check and report a volume that does not double-count overlaps.
    pymeshlab ships its MeshLab DLLs unregistered -- same workaround as
    validate_mesh.py."""
    spec = importlib.util.find_spec("pymeshlab")
    if spec is None:
        return None
    if spec.submodule_search_locations:
        os.add_dll_directory(spec.submodule_search_locations[0])
    import pymeshlab as ml

    V, F = m.arrays()
    ms = ml.MeshSet()
    acc = None
    for i, (base, f0, nf) in enumerate(m._solid_start):
        ff = F[f0:f0 + nf]
        idx = np.unique(ff)
        remap = np.zeros(idx.max() + 1, np.int32)
        remap[idx] = np.arange(len(idx))
        ms.add_mesh(ml.Mesh(V[idx], remap[ff].astype(np.int32)))
        nxt = ms.mesh_number() - 1
        if acc is None:
            acc = nxt
        else:
            ms.apply_filter("generate_boolean_union", first_mesh=acc, second_mesh=nxt)
            acc = ms.mesh_number() - 1
        if i and i % 100 == 0:
            log(f"    welded {i}/{len(m._solid_start)} shells "
                f"({ms.mesh(acc).face_number()} faces)")
    ms.set_current_mesh(acc)
    cm = ms.current_mesh()
    return Welded(np.asarray(cm.vertex_matrix(), float),
                  np.asarray(cm.face_matrix(), np.int64))


def shell_stats(v, f):
    """(components, (open edges, pinched edges), signed volume), BY POSITION.

    By position, not by vertex index, because that is how a slicer sees it -- and
    because the clipper deliberately emits duplicate vertices, so an index-space
    check would report nonsense here."""
    v, f = np.asarray(v, float), np.asarray(f)
    tri = v[f]
    vol = float(np.einsum("ij,ij->i", tri[:, 0],
                          np.cross(tri[:, 1], tri[:, 2])).sum() / 6.0)
    _, inv = np.unique(np.round(v, 5), axis=0, return_inverse=True)
    g = inv.reshape(-1)[f]
    parent = list(range(int(g.max()) + 1))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, b, c in g:
        for x, y in ((a, b), (b, c)):
            ra, rb = find(int(x)), find(int(y))
            if ra != rb:
                parent[ra] = rb
    comps = len({find(int(i)) for i in np.unique(g)})
    e = np.sort(np.concatenate([g[:, [0, 1]], g[:, [1, 2]], g[:, [2, 0]]]), axis=1)
    _, cnt = np.unique(e, axis=0, return_counts=True)
    return comps, (int((cnt == 1).sum()), int((cnt > 2).sum())), vol


def _vol(v, f):
    t = np.asarray(v, float)[np.asarray(f)]
    return float(np.einsum("ij,ij->i", t[:, 0],
                           np.cross(t[:, 1], t[:, 2])).sum() / 6.0)


def _verify(edge=EDGE, steps=STEPS, period=PERIOD, clear=CLEAR):
    """Every claim in the docstring, asserted."""
    # --- the clipper, against solids whose volume is known in closed form -------
    cv, cf = MK.prism([(0, 0), (10, 0), (10, 10), (0, 10)], 0.0, 10.0)
    assert abs(_vol(cv, cf) - 1000.0) < 1e-9
    for z in (2.5, 5.0, 7.5):
        hi = clip_convex(cv, cf, (0, 0, 1.0), z)
        lo = clip_convex(cv, cf, (0, 0, -1.0), -z)
        assert abs(_vol(*hi) - (1000 - 100 * z)) < 1e-7, z
        assert abs(_vol(*hi) + _vol(*lo) - 1000.0) < 1e-7, z
    ov, of = TH.octa((1, 2, 3), 4.6)                 # a node ball, split
    for z in (0.0, 3.0, 4.5):
        hi = clip_convex(ov, of, (0, 0, 1.0), z)
        lo = clip_convex(ov, of, (0, 0, -1.0), -z)
        assert abs(_vol(*hi) + _vol(*lo) - _vol(ov, of)) < 1e-7, z
    bv, bf = TH.beam((0, 0, -5), (0, 3, 12), 3.3, (1, 0, 0))   # a beam, split
    for z in (-4.0, 0.0, 6.0, 11.0):
        hi = clip_convex(bv, bf, (0, 0, 1.0), z)
        lo = clip_convex(bv, bf, (0, 0, -1.0), -z)
        assert abs(_vol(*hi) + _vol(*lo) - _vol(bv, bf)) < 1e-7, z

    for left in (False, True):
        r = Rope(edge, steps, period, clear, left)
        th = math.degrees(abs(r.theta))
        # --- the periodicity identity -----------------------------------------
        assert abs(period * math.degrees(r.psi) - 360.0 * r.m) < 1e-9
        tauH = th - math.degrees(r.psi)
        assert abs(abs(math.degrees(r.lay.tau) * r.H) - abs(tauH)) < 1e-9, tauH
        # the Frenet-frame azimuth must be exactly k*psi
        for k in (-3, 0, 7, 25):
            _C, T, e1, _e2 = r.lay.frame(0, k * r.H)
            T = T / np.linalg.norm(T)
            fren = np.array([math.cos(r.lay.k * k * r.H),
                             math.sin(r.lay.k * k * r.H), 0.0])
            # angle from the Frenet normal to the RMF normal, about the TANGENT
            got = math.atan2(float(np.cross(fren, e1) @ T), float(fren @ e1))
            want = -r.lay.tau * k * r.H
            assert abs((got - want + math.pi) % (2 * math.pi) - math.pi) < 1e-9, k

        # --- STACKING: the braid maps onto itself under the screw, exactly -----
        Rz, t = r.screw_matrix(period)
        klo, khi = int(r.ks[0]), int(r.ks[-1]) - period    # k+period must exist
        for j in range(3):
            for k in range(klo, khi + 1):
                assert np.allclose(r.at(j, k) @ Rz.T + t, r.at(j, k + period),
                                   atol=1e-9), (j, k)
        Rz, t = r.screw_matrix(steps)                       # and per whole module
        for j in range(3):
            assert np.allclose(r.at(j, 0) @ Rz.T + t, r.at(j, steps), atol=1e-9)
        # ties are periodic too, or seams would disagree about the joint
        tset = {(j, k % period, j2, l - k + k % period) for j, k, j2, l, _d in r.ties}
        assert len(tset) == len({(j, k, j2, l - k) for j, k, j2, l, _d in r.ties
                                 if 0 <= k < period})

        # --- FACE ON FACE: consecutive cells share a whole triangle ------------
        for j in range(3):
            for k in range(steps):
                a, b = set(range(k, k + 4)), set(range(k + 1, k + 5))
                assert len(a & b) == 3, "consecutive cells must share 3 vertices"
                fa = [r.at(j, i) for i in sorted(a & b)]
                fb = [r.at(j, i) for i in sorted(b & a)]
                assert np.allclose(fa, fb, atol=1e-12), "shared face must be one face"

        # --- what bending costs: spread, not bias -----------------------------
        L = r.beam_lengths()
        pred = r.R * r.lay.curv
        for d, q in L.items():
            assert abs(q.mean() - edge) < 0.03 * edge, f"m={d} mean drifted {q.mean()}"
            assert np.ptp(q) / edge < 2.5 * pred, f"m={d} spread {np.ptp(q)/edge}"
        assert (r.cell_volumes() > 0).all(), "no cell may collapse or invert"
        # and the Frenet-frame bias the RMF avoids, measured
        P = np.empty((steps + 1, 3))
        for k in range(steps + 1):
            C, _T, e1, e2 = r.lay.frame(0, k * r.H)
            w = r.lay.tau * k * r.H
            cw, sw = math.cos(w), math.sin(w)
            a = k * r.theta
            P[k] = C + r.R * (math.cos(a) * (cw * e1 + sw * e2)
                              + math.sin(a) * (-sw * e1 + cw * e2))
        bias = [np.linalg.norm(P[d:] - P[:-d], axis=1).mean() for d in (1, 2)]
        # Driven APART; which family goes which way depends on sign(tau). Only
        # checkable when there is torsion worth double-applying -- at period 30
        # tauH is 0.19 deg and Frenet vs RMF is not measurably different, which is
        # exactly why tetrahelix.py never had to care about the distinction.
        if abs(math.degrees(r.lay.tau) * r.H) >= 3.0:
            assert max(bias) > 1.03 * edge and min(bias) < 0.97 * edge, bias
        spread = lambda q: max(q) - min(q)
        assert spread(bias) >= spread([v.mean() for v in L.values()]) - 1e-12

        # --- connectivity ------------------------------------------------------
        pairs = {(j, j2) for j, _k, j2, _l, _d in r.ties}
        assert pairs == {(0, 1), (1, 2), (2, 0)}, pairs
        gap = ply_gap(r)
        # A ply must not swallow its neighbour's vertices. Full daylight between
        # node balls (gap > 2*node_r) is what the phase search buys at the default
        # period, but it is not available at every period -- at P=30 the balls
        # overlap 0.4 mm -- so the invariant is the weaker, real one.
        assert gap > r.node_r, f"plies interpenetrating: {gap:.2f} mm"
        assert gap < TIE_MAX * edge, "plies too far apart to tie"
        # the auto phase must beat phase 0, or the search is doing nothing
        zero = ply_gap(Rope(edge, 3 * period, period, clear, left, phase=0.0, pad=2))
        assert gap > zero, (gap, zero)

    # --- PAD really is enough: more padding must change nothing ---------------
    a = build(edge, steps, period, clear, pad=PAD)[0].arrays()
    b = build(edge, steps, period, clear, pad=PAD + 4)[0].arrays()
    assert len(a[1]) == len(b[1]), (len(a[1]), len(b[1]))
    assert abs(_vol(*a) - _vol(*b)) < 1e-6

    # --- the cut faces mate: top cut screwed back must land on the bottom cut --
    m, r, _cut = build(edge, steps, period, clear)
    v, f = m.arrays()
    tri = v[f]
    Rz, t = r.screw_matrix(steps)
    top = tri[np.abs(tri[:, :, 2] - r.height).max(axis=1) < 1e-9].mean(axis=1)
    bot = tri[np.abs(tri[:, :, 2]).max(axis=1) < 1e-9].mean(axis=1)
    assert len(top) and len(top) == len(bot), (len(top), len(bot))
    back = (top - t) @ Rz
    d = np.linalg.norm(back[:, None, :] - bot[None, :, :], axis=2)
    assert d.min(axis=1).max() < 1e-7, d.min(axis=1).max()


def report(a):
    _verify(a.edge, a.steps, a.period, a.clear)
    m, r, cut = build(a.edge, a.steps, a.period, a.clear, a.left, a.hard_lay,
                      a.base, a.stack)
    v, f = m.arrays()
    tag = np.asarray(m.tag)
    lo, hi = v.min(axis=0), v.max(axis=0)
    ext = hi - lo
    th = math.degrees(abs(r.theta))
    psi = math.degrees(r.psi)
    tauH = math.degrees(abs(r.lay.tau) * r.H)

    print(f"Braided rope of three tetrahelices -- STACKABLE MODULE"
          f"{'  (LEFT-hand mirror)' if a.left else ''}"
          f"{'  hard lay root' if a.hard_lay else ''}")
    print(f"  STRAND: true Boerdijk-Coxeter, edge {a.edge:.2f} mm, "
          f"twist {th:.7f} deg/cell (= arccos(-2/3), unfudged)")
    print(f"    {r.steps} cells per strand per module, {3*r.steps} cells total, "
          f"strand radius {r.R:.3f} mm, rise {r.H:.4f} mm/cell")

    print(f"  PERIODICITY -- this is what makes it stack:")
    print(f"    theta/360 = {th/360:.9f} is irrational, so a STRAIGHT tetrahelix "
          f"never repeats")
    print(f"    lay torsion {tauH:.4f} deg/cell takes the twist from {th:.6f} to "
          f"psi = {psi:.6f} deg in the lay-following frame")
    print(f"    psi = 360 * {r.m}/{r.period} exactly -> the braid maps onto itself "
          f"every {r.period} cells (asserted to 1e-9)")
    per_a, per_z = r.screw(r.period)
    mod_a, mod_z = r.screw(r.steps)
    print(f"    screw per period: {math.degrees(per_a):+.4f} deg, "
          f"{per_z:.4f} mm rise")
    print(f"    SCREW PER MODULE: {math.degrees(mod_a) % 360:+.3f} deg, "
          f"{mod_z:.3f} mm rise; the braid is 3-fold so twist each module by "
          f"{math.degrees(mod_a) % 120:+.3f} deg")
    tbl, tmax = feasible_periods(a.edge, a.clear)
    good = [(P, mm, ps, tt) for P, mm, ps, tt, ok in tbl
            if ok and P <= 41 and math.gcd(mm, P) == 1]     # primitive periods only
    print(f"    other reachable periods (lay can give at most {tmax:.4f} deg/cell; "
          f"multiples of a period are the same lay, so only primitives are listed):")
    for P, mm, ps, tt in good:
        note = ""
        if abs(ps - 132.0) < 1e-9:
            note = "  <- tetrahelix.py's 132 deg ten-cycle"
        elif P == r.period:
            note = "  <- this build"
        print(f"        P={P:3d}  m={mm:3d}  psi={ps:11.6f}  "
              f"tauH={tt:+8.4f}{note}")

    print(f"  LAY: radius {r.lay.A:.3f} mm (= (2R + clear {a.clear:.1f})/sqrt3), "
          f"{r.lay_turns:.4f} turns per module")
    print(f"    rise {r.lay.b:.4f} mm per mm of strand -> "
          f"{math.degrees(math.acos(r.lay.b)):.1f} deg off the axis")
    print(f"    axis curvature 1/{1/max(r.lay.curv, 1e-12):.1f} mm -> predicted "
          f"beam spread +-R/Rc = +-{100*r.R*r.lay.curv:.1f}%")

    L = r.beam_lengths()
    print("  BEAMS, measured (the price of bending a tetrahelix -- see docstring):")
    for d, q in L.items():
        print(f"      k->k+{d}:  {q.min():6.3f} .. {q.max():6.3f} mm   "
              f"mean {q.mean():6.3f}  (+-{100*np.ptp(q)/2/a.edge:4.1f}% of edge)")
    cv = r.cell_volumes()
    reg = a.edge ** 3 / (6 * math.sqrt(2))
    print(f"      cell volume {cv.min():.1f} .. {cv.max():.1f} mm3 vs {reg:.1f} "
          f"regular ({100*cv.min()/reg:.0f}% .. {100*cv.max()/reg:.0f}%), "
          f"none inverted")
    print(f"      every consecutive pair still shares a WHOLE face (asserted); "
          f"what bending costs is equilateral faces, not adjacency")
    print(f"      section: 3-sided, side {r.beam_r*math.sqrt(3):.2f} mm, "
          f"node ball r {r.node_r:.2f} mm")

    inside = [t for t in r.ties if 0 <= t[1] < r.steps]
    td = np.array([d for *_x, d in inside])
    print(f"  TIES: {len(inside)} per module, r {r.beam_r*TIE_F:.2f} mm, "
          f"{len(inside)//r.steps if r.steps else 0} per cell-step, "
          f"exactly {r.period}-periodic")
    print(f"      length {td.min():.2f} .. {td.max():.2f} mm; each end lands ON a "
          f"node centre, so it is buried {r.node_r:.2f} mm deep at both ends")
    day = td.min() - 2.0 * r.node_r
    print(f"      plies " + (f"clear each other by {day:.2f} mm of daylight, node "
                             f"ball to node ball, so the ties ARE the joint"
                             if day > 0 else
                             f"also fuse directly: node balls overlap "
                             f"{-day:.2f} mm at the closest approach"))
    print(f"  PHASE: all three strands rolled {math.degrees(r.phase):.2f} deg about "
          f"their own axes -- ply gap {ply_gap(r):.2f} mm against "
          f"{ply_gap(Rope(a.edge, 3*r.period, r.period, a.clear, a.left, a.hard_lay, phase=0.0, pad=2)):.2f} mm "
          f"unrolled. Only possible because psi is rational: the cells sit on "
          f"{r.period} azimuthal rails, so aiming them between the neighbours "
          f"instead of at them is a free choice (see best_phase)")

    print(f"  CUT / STACK: {cut['lo']} solids cut by the bottom plane, "
          f"{cut['hi']} by the top")
    tri = v[f]
    z_off = ring_height(r) if a.base else 0.0
    for name, z in (("bottom", z_off), ("top", r.height + z_off)):
        sel = np.abs(tri[:, :, 2] - z).max(axis=1) < 1e-9
        ar = 0.5 * np.linalg.norm(np.cross(tri[sel][:, 1] - tri[sel][:, 0],
                                           tri[sel][:, 2] - tri[sel][:, 0]),
                                  axis=1).sum() if sel.any() else 0.0
        print(f"      {name} face: {int(sel.sum())} coplanar triangles, "
              f"{ar:.0f} mm2")
    print(f"      the top cut screwed back by one module lands on the bottom cut "
          f"to 1e-7 mm (asserted) -> modules mate exactly")

    print("  triangles by feature:", "  ".join(
        f"{q}={int((tag == q).sum())}" for q in sorted(set(tag)) if "#" not in q))
    if a.stack > 1:
        print(f"  STACKED {a.stack} modules in this file")
    print(f"  bbox x {lo[0]:7.1f}..{hi[0]:6.1f}  y {lo[1]:7.1f}..{hi[1]:6.1f}"
          f"  z {lo[2]:7.1f}..{hi[2]:6.1f}")
    fits = "FITS" if (ext <= 320.0 + 1e-9).all() else "OVER BED"
    print(f"  envelope {ext[0]:.1f} x {ext[1]:.1f} x {ext[2]:.1f} mm   "
          f"(bed 320 cubed: {fits}),  {ext[2]/max(ext[0], ext[1]):.1f}:1 aspect")
    print(f"  solids {len(m._solid_start)}  triangles {len(f)}")

    inc, span = [], []
    for j in range(3):
        for k in range(r.steps):
            for d in (1, 2, 3):
                q = r.at(j, k + d) - r.at(j, k)
                n = np.linalg.norm(q)
                inc.append(math.degrees(math.asin(abs(q[2]) / n)))
                span.append(math.hypot(q[0], q[1]))
    inc, span = np.array(inc), np.array(span)
    sh = inc < 20.0
    print(f"  PRINT: beam inclination {inc.min():.1f} .. {inc.max():.1f} deg from "
          f"horizontal, median {np.median(inc):.1f}")
    print(f"      {int(sh.sum())}/{len(inc)} under 20 deg, longest unsupported "
          f"horizontal run {span[sh].max():.0f} mm -> TREE SUPPORTS")
    print(f"      module stands on its own flat bottom cut"
          + (" plus the base ring" if a.base else "")
          + f"; {ext[2]/max(ext[0], ext[1]):.1f}:1 aspect -- add a brim")
    return m, r


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--edge", type=float, default=EDGE,
                    help="tetrahedron edge of each strand, mm (default 22)")
    ap.add_argument("--steps", type=int, default=STEPS,
                    help="cells per strand per module; rounded up to a whole "
                         "number of periods (default 39)")
    ap.add_argument("--period", type=int, default=PERIOD,
                    help="cells per screw period (default 3). Only the convergents "
                         "of theta/360 are reachable -- the report lists them")
    ap.add_argument("--clear", type=float, default=CLEAR,
                    help="mm between strand vertex cylinders (default 5)")
    ap.add_argument("--left", action="store_true",
                    help="mirror image: strand twist and lay both flip")
    ap.add_argument("--hard-lay", action="store_true",
                    help="the other lay that gives the same period -- tighter "
                         "braid, more cell distortion")
    ap.add_argument("--stack", type=int, default=1,
                    help="write N modules already screwed into place (default 1)")
    ap.add_argument("--base", action="store_true",
                    help="add the optional bar-triangle base under module one")
    ap.add_argument("--weld", action="store_true",
                    help="CSG-union every shell into one 2-manifold solid. Tens of "
                         "minutes -- slicers union the raw soup themselves")
    ap.add_argument("--out", default="braided_tetrahelix.stl")
    a = ap.parse_args()
    m, r = report(a)

    out = m
    if a.weld:
        print(f"\n  welding {len(m._solid_start)} shells (pymeshlab CSG union)...")
        out = weld(m)
        if out is None:
            print("  pymeshlab unavailable -- writing the raw soup instead")
            out = m
    v, f = out.arrays()
    comps, (open_e, pinch_e), vol = shell_stats(v, f)
    print(f"\n  OUTPUT ({'welded' if out is not m else 'raw soup'}), by position: "
          f"{len(f)} triangles, {comps} connected "
          f"{'body' if comps == 1 else 'bodies'}, {open_e} open edges, "
          f"{pinch_e} pinched edges")
    if out is not m:
        if comps == 1:
            print("  -> ONE body: three tetrahelices and every tie, a single closed "
                  "solid, measured")
        if pinch_e:
            print(f"     {pinch_e} pinches are CSG seam points with coincident "
                  f"vertices -- closed and correctly wound either side, so slicers "
                  f"fill it")
    elif a.stack > 1 and pinch_e:
        print(f"  -> {comps} bodies from {len(m._solid_start)} overlapping closed "
              f"shells, which slicers union on slice")
        print(f"     the {pinch_e} pinched edges are the {a.stack - 1} seams: two "
              f"modules' cut faces landing on each other exactly, which is the "
              f"stacking working, not a defect")
    else:
        print(f"  -> {comps} bodies is expected: {len(m._solid_start)} overlapping "
              f"closed shells, which slicers union on slice; --weld does it up front")
    print(f"  volume {vol/1000:.1f} cm3 -> {vol*1.27/1000:.0f} g PETG"
          + ("" if out is not m else "  (soup -- overlaps double-counted)"))
    n = MK.write_stl(a.out, out)
    print(f"\nwrote {a.out} ({n} triangles)")


if __name__ == "__main__":
    main()
