#!/usr/bin/env python3
"""Self-contained orthographic z-buffer rasteriser (numpy + PIL).

Flat per-triangle Lambert + specular, 2x supersampled, two-sided shading so
section views light their interior faces."""

import math

import numpy as np
from PIL import Image

SS = 2
BG = (242, 242, 240)
MAT = {'frame': {'col': (108, 112, 120), 'amb': 0.3, 'spec': 0.1, 'shin': 14.0}, 'glass': {'col': (176, 198, 214), 'amb': 0.34, 'spec': 0.85, 'shin': 48.0}, 'beam': {'col': (86, 226, 122), 'amb': 0.72, 'spec': 0.2, 'shin': 20.0}, 'pad': {'col': (196, 148, 96), 'amb': 0.34, 'spec': 0.12, 'shin': 14.0}}


def basis(az, el):
    a, e = math.radians(az), math.radians(el)
    d = np.array([math.cos(e) * math.cos(a), math.cos(e) * math.sin(a), math.sin(e)])
    up_w = np.array([0.0, 0.0, 1.0]) if abs(d[2]) < 0.985 else np.array([0.0, 1.0, 0.0])
    r = np.cross(up_w, d)
    r /= np.linalg.norm(r)
    u = np.cross(d, r)
    return r, u, d

def render(layers, az, el, w=1100, h=1100, target=None, span=None,
           cull=None, out="out.png"):
    """layers: list of (verts, faces, material-name)."""
    r, u, d = basis(az, el)
    key = 0.40 * r + 0.58 * u + 0.72 * d
    key /= np.linalg.norm(key)
    fill = -0.62 * r + 0.18 * u + 0.55 * d
    fill /= np.linalg.norm(fill)
    half = key + d
    half /= np.linalg.norm(half)

    prepared = []
    for layer in layers:
        verts, faces, mat = layer[:3]
        cullable = layer[3] if len(layer) > 3 else True
        tri = verts[faces]
        if cull is not None and cullable:
            keep = ~cull(tri.mean(axis=1))
            tri = tri[keep]
        if not len(tri):
            continue
        nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        ln = np.linalg.norm(nrm, axis=1, keepdims=True)
        good = ln[:, 0] > 1e-9
        tri, nrm, ln = tri[good], nrm[good], ln[good]
        nrm = nrm / ln
        nrm[nrm @ d < 0] *= -1.0                       # two-sided shading
        m = MAT[mat]
        shade = (m["amb"] + 0.78 * np.clip(nrm @ key, 0, None)
                 + 0.30 * np.clip(nrm @ fill, 0, None))
        spec = m["spec"] * np.clip(nrm @ half, 0, None) ** m["shin"]
        col = np.clip(np.asarray(m["col"], float) * shade[:, None]
                      + 255.0 * spec[:, None], 0, 255)
        prepared.append((tri, col))

    allpts = np.concatenate([t.reshape(-1, 3) for t, _ in prepared])
    if target is None:
        target = 0.5 * (allpts.min(axis=0) + allpts.max(axis=0))
    target = np.asarray(target, float)
    proj = np.stack([(allpts - target) @ r, (allpts - target) @ u], axis=1)
    if span is None:
        ext = np.abs(proj).max(axis=0)
        span = 2.10 * max(ext[0], ext[1] * w / h)
    W, Hh = w * SS, h * SS
    s = W / span
    cx, cy = W / 2.0, Hh / 2.0

    img = np.zeros((Hh, W, 3), np.float32)
    img[:] = BG
    zbuf = np.full((Hh, W), -1e18, np.float32)

    for tri, col in prepared:
        q = tri - target
        sx = cx + (q @ r) * s
        sy = cy - (q @ u) * s
        dz = q @ d
        x0 = np.clip(np.floor(sx.min(axis=1)).astype(int), 0, W - 1)
        x1 = np.clip(np.ceil(sx.max(axis=1)).astype(int), 0, W - 1)
        y0 = np.clip(np.floor(sy.min(axis=1)).astype(int), 0, Hh - 1)
        y1 = np.clip(np.ceil(sy.max(axis=1)).astype(int), 0, Hh - 1)
        for i in range(len(tri)):
            if x1[i] < x0[i] or y1[i] < y0[i]:
                continue
            ax, ay = sx[i, 0], sy[i, 0]
            bx, by = sx[i, 1], sy[i, 1]
            gx, gy = sx[i, 2], sy[i, 2]
            area = (bx - ax) * (gy - ay) - (by - ay) * (gx - ax)
            if abs(area) < 1e-9:
                continue
            px = np.arange(x0[i], x1[i] + 1) + 0.5
            py = np.arange(y0[i], y1[i] + 1) + 0.5
            PX, PY = np.meshgrid(px, py)
            w0 = ((bx - ax) * (PY - ay) - (by - ay) * (PX - ax)) / area
            w1 = ((gx - bx) * (PY - by) - (gy - by) * (PX - bx)) / area
            inside = (w0 >= -1e-6) & (w1 >= -1e-6) & (w0 + w1 <= 1 + 1e-6)
            if not inside.any():
                continue
            lam2, lam0 = w0, w1
            lam1 = 1.0 - lam0 - lam2
            z = lam0 * dz[i, 0] + lam1 * dz[i, 1] + lam2 * dz[i, 2]
            sub = zbuf[y0[i]:y1[i] + 1, x0[i]:x1[i] + 1]
            hit = inside & (z > sub)
            if not hit.any():
                continue
            sub[hit] = z[hit]
            img[y0[i]:y1[i] + 1, x0[i]:x1[i] + 1][hit] = col[i]

    im = Image.fromarray(img.astype(np.uint8)).resize((w, h), Image.LANCZOS)
    im.save(out)
    print("wrote", out)
