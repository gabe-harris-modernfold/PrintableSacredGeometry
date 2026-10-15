"""numpy re-implementation of web/volumetric-sound-forms.html's pipeline.

Same frame placement, periodic Hann, band map, dB/gamma/contrast, layouts, blur,
erosion, foot, solid/hollow cylinder topology and binary STL as the page, so the
two can be diffed field-by-field and vertex-by-vertex.

    python sound_form_reference.py --write-test-wav test_tone.wav
    python sound_form_reference.py --wav test_tone.wav [--layout mirror|plain|transpose] [--hollow] [--nt 192] [--nf 96]
    python sound_form_reference.py --compare        # against run_pure.mjs's demo_field.*, demo_cavity.* and js_*.stl

Rounding conventions pinned to the JS: Math.round -> floor(x + 0.5); Math.ceil(x - 1e-9);
values stored to Float32Array are cast to float32 at the same points."""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TAU = 2.0 * np.pi
F32 = np.float32


def js_round(x):
    return int(np.floor(x + 0.5))


def defaults():
    return dict(layout="mirror", fft=2048, nt=256, nf=96, fmin=40, fmax=12000, scale="log",
                db=-60, gamma=1, contrast=1, bz=1, bt=1,
                r0=45, depth=12, h=160, foot=6, fade=0.02, hollow=0, wall=2, flr=3,
                fit=1, bed=300, t0=0, t1=0, slope=0,
                form="chart", cr0=50, cr1=40, ch=120, cwall=2, cflr=3, ctop="open", ctune=60, cq=30,
                cfill=15, csmooth=1, cvox=1.2, csrc=0.62, cmode="clip", ct=0, chalf=0)


# ─── signal ───────────────────────────────────────────────────────────────────
DEMO_HZ = 100.0


def synth_demo(sr, seconds):
    """the page's demo: a pure 100 Hz sine at 0.5, with 20 ms raised-cosine fades at both ends"""
    n = int(np.floor(sr * seconds))
    i = np.arange(n, dtype=np.float64)
    nf = min(int(np.floor(0.02 * sr)), n >> 1)
    g = np.full(n, 0.5)
    head = i < nf
    tail = (i >= n - nf) & ~head
    g[head] *= 0.5 - 0.5 * np.cos(np.pi * i[head] / nf)
    g[tail] *= 0.5 - 0.5 * np.cos(np.pi * (n - 1 - i[tail]) / nf)
    return (g * np.sin(TAU * DEMO_HZ * i / sr)).astype(F32)


def span(n, sr, t0, t1):
    s0 = int(np.floor(max(0.0, t0 or 0.0) * sr))
    s1 = int(np.floor(t1 * sr)) if (t1 and t1 > 0) else n
    s1 = min(n, max(s1, 0))
    s0 = min(s0, max(0, s1 - 2))
    if s1 - s0 < 2:
        s0, s1 = 0, n
    return s0, s1


def effective_rows(n, sr, P):
    s0, s1 = span(n, sr, P["t0"], P["t1"])
    return max(2, min(int(P["nt"]), (s1 - s0) // 16))


def hann(N):
    return (0.5 - 0.5 * np.cos(TAU * np.arange(N) / N)).astype(F32)


def stft(samples, sr, P):
    N = int(P["fft"])
    nT = effective_rows(len(samples), sr, P)
    nBins = N // 2 + 1
    s0, s1 = span(len(samples), sr, P["t0"], P["t1"])
    hop = (s1 - s0) / nT
    K = 1 if hop <= N else min(8, int(np.ceil(hop / N)))
    w = hann(N).astype(np.float64)
    x = samples.astype(np.float64)
    spec = np.zeros((nT, nBins), dtype=F32)
    i = np.arange(nT)
    c0 = s0 + (i + 0.5) * hop
    for k in range(K):
        c = c0 + (k - (K - 1) / 2) * hop / K
        start = np.floor(c - N // 2 + 0.5).astype(np.int64)
        idx = start[:, None] + np.arange(N)[None, :]
        valid = (idx >= 0) & (idx < len(x))
        frames = np.where(valid, x[np.clip(idx, 0, len(x) - 1)], 0.0) * w[None, :]
        X = np.fft.rfft(frames, axis=1)
        spec += ((X.real ** 2 + X.imag ** 2) * (1.0 / K)).astype(F32)
    return dict(spec=spec, nT=nT, nBins=nBins, N=N, hop=hop, K=K, s0=s0, s1=s1)


def rms_rows(samples, sr, t0, t1, nT):
    s0, s1 = span(len(samples), sr, t0, t1)
    hop = (s1 - s0) / nT
    x = samples.astype(np.float64)
    env = np.zeros(nT)
    for i in range(nT):
        a = int(np.floor(s0 + i * hop))
        b = max(a + 1, int(np.floor(s0 + (i + 1) * hop)))
        env[i] = np.sqrt(np.mean(x[a:b] ** 2))
    if env.max() > 0:
        env /= env.max()
    return env.astype(F32)


# ─── bands and field ──────────────────────────────────────────────────────────
def band_edges(fmin, fmax, nF, scale):
    k = np.arange(nF + 1, dtype=np.float64)
    if scale == "lin":
        return fmin + (fmax - fmin) * k / nF
    return fmin * np.power(fmax / fmin, k / nF)


def band_map(edges, sr, N):
    nBins = N // 2 + 1
    nF = len(edges) - 1
    b2b = np.full(nBins, -1, dtype=np.int32)
    count = np.zeros(nF, dtype=np.int32)
    k = 0
    for b in range(nBins):
        f = b * sr / N
        if f < edges[0]:
            continue
        while k < nF and f >= edges[k + 1]:
            k += 1
        if k >= nF:
            break
        b2b[b] = k
        count[k] += 1
    lo = np.full(nF, -1, dtype=np.int32)
    hi = np.full(nF, -1, dtype=np.int32)
    w = np.zeros(nF, dtype=F32)
    prev = -1
    for i in range(nF):
        if count[i] > 0:
            prev = i
        lo[i] = prev
    nxt = -1
    for i in range(nF - 1, -1, -1):
        if count[i] > 0:
            nxt = i
        hi[i] = nxt
    for i in range(nF):
        if count[i] > 0:
            lo[i] = hi[i] = i
            w[i] = 0
            continue
        if lo[i] < 0:
            lo[i] = hi[i]
        if hi[i] < 0:
            hi[i] = lo[i]
        w[i] = 0 if hi[i] == lo[i] else F32((i - lo[i]) / (hi[i] - lo[i]))
    return dict(binToBand=b2b, count=count, lo=lo, hi=hi, w=w, nF=nF, nBins=nBins)


def bin_bands(spec, m):
    nT = spec.shape[0]
    nF = m["nF"]
    M = np.zeros((m["nBins"], nF))
    used = m["binToBand"] >= 0
    M[np.nonzero(used)[0], m["binToBand"][used]] = 1.0
    power = spec.astype(np.float64) @ M
    power = power.astype(F32)
    for k in range(nF):
        if m["count"][k] > 0 or m["lo"][k] < 0:
            continue
        wk = np.float64(m["w"][k])
        power[:, k] = (power[:, m["lo"][k]].astype(np.float64) * (1 - wk)
                       + power[:, m["hi"][k]].astype(np.float64) * wk).astype(F32)
    return power


def to_field(power, floor_db, gamma, contrast):
    pmax = float(power.max())
    if not pmax > 0:
        return np.zeros_like(power, dtype=F32)
    p = power.astype(np.float64)
    db = 10.0 * np.log10((p + 1e-20) * (1.0 / pmax))
    f = np.clip((db - floor_db) * (1.0 / (-floor_db)), 0, 1)
    f = np.power(f, gamma)
    f = np.clip(0.5 + (f - 0.5) * contrast, 0, 1)
    return f.astype(F32)


def layout_surface(field, layout, fade):
    if layout == "transpose":
        S = field.T.astype(np.float64).copy()
        nR, nC = S.shape
        nf = js_round((fade or 0) * nC)
        j = 0
        while j < nf and j < nC - 1 - j:
            g = 0.5 - 0.5 * np.cos(np.pi * (j + 0.5) / nf)
            S[:, j] = (S[:, j].astype(F32) * g).astype(F32)
            S[:, nC - 1 - j] = (S[:, nC - 1 - j].astype(F32) * g).astype(F32)
            j += 1
        return S.astype(F32)
    if layout == "plain":
        return field.astype(F32).copy()
    return np.concatenate([field, field[:, ::-1]], axis=1).astype(F32)


def _wrap_sum(S64, r):
    acc = np.zeros_like(S64)
    for d in range(-r, r + 1):
        acc += np.roll(S64, -d, axis=1)
    return acc


def _clamp_sum(S64, r):
    nR = S64.shape[0]
    i = np.arange(nR)
    acc = np.zeros_like(S64)
    for d in range(-r, r + 1):
        acc += S64[np.clip(i + d, 0, nR - 1)]
    return acc


def blur_surface(S, rz, rt, passes=2):
    nR, nC = S.shape
    rz = max(0, js_round(rz))
    rt = max(0, min(js_round(rt), (nC - 1) >> 1))
    if not rz and not rt:
        return S
    for _ in range(passes):
        if rt > 0:
            S = (_wrap_sum(S.astype(np.float64), rt) * (1.0 / (2 * rt + 1))).astype(F32)
        if rz > 0:
            S = (_clamp_sum(S.astype(np.float64), rz) * (1.0 / (2 * rz + 1))).astype(F32)
    return S


def js_ceil(x):
    return int(np.ceil(x - 1e-9))


def erode_surface(S, rz, rt):
    nR, nC = S.shape
    rz = max(0, js_ceil(rz))
    rt = max(0, min(js_ceil(rt), (nC - 1) >> 1))
    if not rz and not rt:
        return S
    if rt > 0:
        m = S.copy()
        for d in range(-rt, rt + 1):
            m = np.minimum(m, np.roll(S, -d, axis=1))
        S = m
    if rz > 0:
        i = np.arange(nR)
        m = S.copy()
        for d in range(-rz, rz + 1):
            m = np.minimum(m, S[np.clip(i + d, 0, nR - 1)])
        S = m
    return S


def pad_foot(S, n_foot):
    if n_foot <= 0:
        return S
    return np.concatenate([np.zeros((n_foot, S.shape[1]), dtype=F32), S], axis=0)


# ─── the solid ────────────────────────────────────────────────────────────────
def grid_spec(P, nT):
    nF = int(P["nf"])
    layout = P["layout"]
    nR0 = nF if layout == "transpose" else nT
    nC = 2 * nF if layout == "mirror" else nF if layout == "plain" else nT
    dz = P["h"] / (nR0 - 1)
    n_foot = max(0, js_round(P["foot"] / dz))
    nR = nR0 + n_foot
    nRi = 0
    if int(P["hollow"]):
        H = (nR - 1) * dz
        floor_mm = min(max(1.2, float(P["flr"])), H * 0.5)
        k0 = int(np.ceil(floor_mm / dz - 1e-9))
        if k0 > nR - 2:
            k0 = nR - 2
        nRi = nR - k0
    F = 2 * nC * (nR + nRi) if int(P["hollow"]) else 2 * nR * nC
    return dict(nR0=nR0, nR=nR, nC=nC, nRi=nRi, nFoot=n_foot, dz=dz, F=F)


def _sheet(base, n_rows, nC, flip):
    i = np.arange(n_rows - 1)[:, None]
    j = np.arange(nC)[None, :]
    j1 = (j + 1) % nC
    a = base + i * nC + j
    b = base + i * nC + j1
    c = base + (i + 1) * nC + j1
    d = base + (i + 1) * nC + j
    if flip:
        t = np.stack([np.stack([a, c, b], -1), np.stack([a, d, c], -1)], -2)
    else:
        t = np.stack([np.stack([a, b, c], -1), np.stack([a, c, d], -1)], -2)
    return t.reshape(-1, 3)


def _fan(apex, ring, nC, up):
    j = np.arange(nC)
    j1 = (j + 1) % nC
    if up:
        return np.stack([np.full(nC, apex), ring + j, ring + j1], -1)
    return np.stack([np.full(nC, apex), ring + j1, ring + j], -1)


def _annulus(outer, inner, nC, up):
    j = np.arange(nC)
    j1 = (j + 1) % nC
    a, b, c, d = outer + j, outer + j1, inner + j1, inner + j
    if up:
        t = np.stack([np.stack([a, b, c], -1), np.stack([a, c, d], -1)], -2)
    else:
        t = np.stack([np.stack([a, c, b], -1), np.stack([a, d, c], -1)], -2)
    return t.reshape(-1, 3)


def mesh_cylinder(S, P, dz):
    nR, nC = S.shape
    R0 = float(P["r0"])
    depth = max(0.0, float(P["depth"]))
    H = (nR - 1) * dz
    hollow = bool(int(P["hollow"]))
    r_min = max(1.0, 0.05 * R0)
    Rout = np.maximum(r_min, R0 + depth * S.astype(np.float64)).astype(F32)
    slope = max(0.0, min(80.0, float(P.get("slope", 0) or 0)))
    tanS = np.tan(slope * np.pi / 180) if slope > 0 else 0.0
    if slope > 0:
        lim = dz * tanS
        for i in range(nR - 2, -1, -1):
            Rout[i] = np.maximum(Rout[i].astype(np.float64), Rout[i + 1].astype(np.float64) - lim).astype(F32)
    nRi, Rin, floor_mm, wall = 0, None, 0.0, 0.0
    if hollow:
        wall = max(0.8, float(P["wall"]))
        floor_mm = min(max(1.2, float(P["flr"])), H * 0.5)
        k0 = int(np.ceil(floor_mm / dz - 1e-9))
        if k0 > nR - 2:
            k0 = nR - 2
        nRi = nR - k0
        Rer = erode_surface(Rout.copy(), wall / dz, wall / (R0 * TAU / nC))
        r_in_min = max(1.0, 0.3 * R0)
        Rin = np.zeros((nRi, nC), dtype=F32)
        for k in range(nRi):
            z = floor_mm + (H - floor_mm) * k / (nRi - 1)
            fi = z / dz
            i0 = min(int(np.floor(fi)), nR - 1)
            i1 = min(i0 + 1, nR - 1)
            t = min(1.0, max(0.0, fi - i0))
            r = Rer[i0].astype(np.float64) * (1 - t) + Rer[i1].astype(np.float64) * t
            Rin[k] = np.maximum(r_in_min, r - wall).astype(F32)
        if slope > 0 and nRi > 1:
            lim_in = (H - floor_mm) / (nRi - 1) * tanS
            for k in range(nRi - 2, -1, -1):
                Rin[k] = np.minimum(Rin[k].astype(np.float64), Rin[k + 1].astype(np.float64) + lim_in).astype(F32)
    nV = nR * nC + nRi * nC + 2
    th = TAU * np.arange(nC) / nC
    cs, sn = np.cos(th), np.sin(th)
    pos = np.zeros((nV, 3), dtype=F32)
    z = (np.arange(nR) * dz)[:, None]
    outer = np.stack([Rout.astype(np.float64) * cs[None, :], Rout.astype(np.float64) * sn[None, :], np.broadcast_to(z, (nR, nC))], -1)
    pos[: nR * nC] = outer.reshape(-1, 3).astype(F32)
    inner_base = nR * nC
    if hollow:
        zi = (floor_mm + (H - floor_mm) * np.arange(nRi) / (nRi - 1))[:, None]
        inner = np.stack([Rin.astype(np.float64) * cs[None, :], Rin.astype(np.float64) * sn[None, :], np.broadcast_to(zi, (nRi, nC))], -1)
        pos[inner_base: inner_base + nRi * nC] = inner.reshape(-1, 3).astype(F32)
    B = inner_base + nRi * nC
    T = B + 1
    pos[B] = (0, 0, 0)
    pos[T] = (0, 0, floor_mm if hollow else H)
    parts = [_sheet(0, nR, nC, False), _fan(B, 0, nC, False)]
    if not hollow:
        parts.append(_fan(T, (nR - 1) * nC, nC, True))
    else:
        parts.append(_sheet(inner_base, nRi, nC, True))
        parts.append(_fan(T, inner_base, nC, True))
        parts.append(_annulus((nR - 1) * nC, inner_base + (nRi - 1) * nC, nC, True))
    idx = np.concatenate(parts, 0).astype(np.int64)
    return dict(pos=pos, idx=idx, nV=nV, nF=len(idx), nR=nR, nC=nC, nRi=nRi, hollow=hollow)


def fit_to_bed(pos, bed_max):
    size = pos.max(0) - pos.min(0)
    big = float(size.max())
    s = bed_max / big if big > bed_max else 1.0
    if s < 1:
        pos[:] = (pos.astype(np.float64) * s).astype(F32)
    return s


def solid_check(pos, idx):
    nV = len(pos)
    a, b, c = idx[:, 0], idx[:, 1], idx[:, 2]
    und = np.concatenate([np.minimum(a, b) * nV + np.maximum(a, b), np.minimum(b, c) * nV + np.maximum(b, c),
                          np.minimum(c, a) * nV + np.maximum(c, a)])
    dire = np.concatenate([a * nV + b, b * nV + c, c * nV + a])
    _, cnt = np.unique(und, return_counts=True)
    _, dcnt = np.unique(dire, return_counts=True)
    p = pos.astype(np.float64)
    vol = float(np.sum(np.einsum("ij,ij->i", p[a], np.cross(p[b], p[c]))) / 6.0)
    return dict(V=nV, E=len(cnt), F=len(idx), euler=nV - len(cnt) + len(idx), vol=vol,
                badEdges=int(np.sum(cnt != 2)), badDir=int(np.sum(dcnt != 1)))


STL_DTYPE = np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)), ("attr", "<u2")])


def write_stl(path, pos, idx, header="Promethean Theoria - volumetric sound form (numpy reference)"):
    p = pos.astype(np.float64)
    tri = p[idx]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    L = np.linalg.norm(n, axis=1)
    n = np.where(L[:, None] > 0, n / np.where(L[:, None] > 0, L[:, None], 1), 0)
    rec = np.zeros(len(idx), dtype=STL_DTYPE)
    rec["n"] = n.astype(F32)
    rec["v"] = pos[idx]
    hdr = header.encode("ascii", "replace")[:80].ljust(80, b"\0")
    with open(path, "wb") as fh:
        fh.write(hdr)
        fh.write(np.uint32(len(idx)).tobytes())
        fh.write(rec.tobytes())
    return len(idx)


def read_stl(path):
    with open(path, "rb") as fh:
        fh.read(80)
        n = int(np.frombuffer(fh.read(4), dtype="<u4")[0])
        rec = np.frombuffer(fh.read(n * 50), dtype=STL_DTYPE)
    return rec["v"]


# ─── the pipeline ─────────────────────────────────────────────────────────────
def pipeline(samples, sr, P):
    an = stft(samples, sr, P)
    m = band_map(band_edges(P["fmin"], P["fmax"], int(P["nf"]), P["scale"]), sr, an["N"])
    power = bin_bands(an["spec"], m)
    field = to_field(power, P["db"], P["gamma"], P["contrast"])
    S = layout_surface(field, P["layout"], P["fade"])
    S = blur_surface(S, P["bz"], P["bt"], 2)
    g = grid_spec(P, an["nT"])
    S = pad_foot(S, g["nFoot"])
    mesh = mesh_cylinder(S, P, g["dz"])
    scale = fit_to_bed(mesh["pos"], float(P["bed"])) if int(P["fit"]) else 1.0
    stats = solid_check(mesh["pos"], mesh["idx"])
    stats["scale"] = scale
    stats["size"] = (mesh["pos"].max(0) - mesh["pos"].min(0)).tolist()
    return dict(an=an, map=m, power=power, field=field, S=S, mesh=mesh, stats=stats, grid=g)


def load_wav(path):
    from scipy.io import wavfile
    sr, data = wavfile.read(path)
    if data.dtype == np.int16:
        x = data.astype(np.float64) / 32768.0
    elif data.dtype == np.int32:
        x = data.astype(np.float64) / 2147483648.0
    elif data.dtype == np.uint8:
        x = (data.astype(np.float64) - 128.0) / 128.0
    else:
        x = data.astype(np.float64)
    if x.ndim == 2:
        x = x.mean(axis=1)
    return int(sr), x.astype(F32)


def write_test_wav(path, seconds=8, sr=44100):
    from scipy.io import wavfile
    x = synth_demo(sr, seconds).astype(np.float64) * 0.7
    pcm = np.round(x * 32767).astype(np.int16)
    wavfile.write(path, sr, np.stack([pcm, pcm], axis=1))
    print(f"wrote {path}: {seconds} s stereo 16-bit at {sr} Hz")


def compare(tol_field=5e-3, tol_mean=5e-4, tol_mm=1e-3):
    with open(os.path.join(HERE, "demo_field.json")) as fh:
        js = json.load(fh)
    P = dict(defaults())
    P.update(js["P"])
    sr = js["sr"]
    x = synth_demo(sr, js["seconds"])
    ref = pipeline(x, sr, P)
    fails = 0

    def ok(cond, msg):
        nonlocal fails
        print(("  ok   " if cond else "  FAIL ") + msg)
        if not cond:
            fails += 1

    jf = js["field"]
    ok(ref["an"]["nT"] == jf["nT"] and ref["an"]["N"] == jf["N"], f"rows {ref['an']['nT']} == {jf['nT']}, N {ref['an']['N']} == {jf['N']}")
    ok(abs(ref["an"]["hop"] - jf["hop"]) < 1e-9 and ref["an"]["K"] == jf["K"], f"hop {ref['an']['hop']:.4f} K {ref['an']['K']}")
    ok(np.array_equal(ref["map"]["binToBand"], np.array(jf["binToBand"], dtype=np.int32)), "binToBand identical")
    ok(np.array_equal(ref["map"]["count"], np.array(jf["count"], dtype=np.int32)), "band bin counts identical")
    js_field = np.fromfile(os.path.join(HERE, "demo_field.f32"), dtype="<f4").reshape(jf["nT"], jf["nF"])
    d = np.abs(ref["field"].astype(np.float64) - js_field)
    ok(d.max() <= tol_field and d.mean() <= tol_mean, f"field max|d| {d.max():.2e} mean {d.mean():.2e}")
    js_S = np.fromfile(os.path.join(HERE, "demo_surface.f32"), dtype="<f4").reshape(js["surface"]["nR"], js["surface"]["nC"])
    dS = np.abs(ref["S"].astype(np.float64) - js_S)
    ok(ref["S"].shape == js_S.shape and dS.max() <= tol_field, f"surface {ref['S'].shape} max|d| {dS.max():.2e}")
    for run in js["runs"]:
        Q = dict(P)
        Q["layout"] = run["layout"]
        Q["hollow"] = run["hollow"]
        Q["slope"] = run.get("slope", 0)
        r = pipeline(x, sr, Q)
        s = r["stats"]
        ok(s["V"] == run["V"] and s["F"] == run["F"] and s["euler"] == 2 and s["badEdges"] == 0 and s["badDir"] == 0,
           f"{run['name']}: V {s['V']} F {s['F']} Euler {s['euler']} odd {s['badEdges']} rev {s['badDir']}")
        ok(abs(s["vol"] - run["vol"]) <= 1e-3 * abs(run["vol"]) + 1e-6, f"{run['name']}: volume {s['vol'] / 1000:.3f} vs js {run['vol'] / 1000:.3f} cm3")
        js_tri = read_stl(os.path.join(HERE, run["name"])).astype(np.float64)
        py_tri = r["mesh"]["pos"][r["mesh"]["idx"]].astype(np.float64)
        if js_tri.shape == py_tri.shape:
            dv = np.abs(js_tri - py_tri).max()
            ok(dv <= tol_mm, f"{run['name']}: vertex max|d| {dv:.2e} mm over {len(py_tri)} triangles in order")
        else:
            ok(False, f"{run['name']}: triangle count {js_tri.shape[0]} vs {py_tri.shape[0]}")
        write_stl(os.path.join(HERE, "py_" + run["name"][3:]), r["mesh"]["pos"], r["mesh"]["idx"])
    print("compare:", "PASS" if not fails else f"{fails} FAILURE(S)")
    return 1 if fails else 0


# ─── the cavity form ──────────────────────────────────────────────────────────
# Independent of the page's numerics where it matters: zeros of J_m' from scipy.special.jnp_zeros, axial
# eigenpairs from LAPACK (scipy.linalg.eigh_tridiagonal) on the same finite-volume operator, the field from
# scipy.special.jv on every voxel (the page uses radial tables), smoothing by scipy.ndimage.gaussian_filter1d,
# islands by scipy.ndimage.label, and the surface by VTK's marching cubes (the page uses marching tetrahedra).
KMAX = 22.0
K11 = 1.8411837813406593


def cavity_geom(P):
    cr0, cr1, ch = max(1.0, float(P["cr0"])), max(1.0, float(P["cr1"])), max(1.0, float(P["ch"]))
    a0 = max(cr0, cr1)
    return dict(cr0=cr0, cr1=cr1, ch=ch, a0=a0, Hn=ch / a0, open=P.get("ctop", "open") != "closed")


def cavity_modes(P):
    from scipy.special import jnp_zeros, jv
    from scipy.linalg import eigh_tridiagonal
    G = cavity_geom(P)
    cr0, cr1, a0, Hn = G["cr0"], G["cr1"], G["a0"], G["Hn"]
    nz = int(min(1600, max(240, np.ceil(4 * KMAX * Hn))))
    h = Hn / nz
    ih2 = 1.0 / (h * h)
    i = np.arange(nz)
    aC = (cr0 + (cr1 - cr0) * (i + 0.5) / nz) / a0
    aF2 = ((cr0 + (cr1 - cr0) * np.arange(nz + 1) / nz) / a0) ** 2
    cl = np.where(i > 0, aF2[i], 0.0)
    cr = np.where(i < nz - 1, aF2[np.minimum(i + 1, nz)], 0.0)
    e = -aF2[1:nz] * ih2 / (aC[:-1] * aC[1:])
    groups, modes = [], []
    m = 0
    while True:
        zs = ([0.0] if m == 0 else []) + [float(z) for z in jnp_zeros(m, 40)]
        zs = [z for z in zs if z <= KMAX]
        if not zs:
            break
        for n, kap in enumerate(zs):
            d = ((cl + cr) * ih2 + kap * kap) / (aC * aC)
            if G["open"]:
                d[-1] += 2 * aF2[nz] * ih2 / (aC[-1] * aC[-1])
            vals, vecs = eigh_tridiagonal(d, e, select="v", select_range=(1e-6, KMAX * KMAX))
            if not len(vals):
                continue
            Jw = 1.0 if kap == 0 else float(jv(m, kap))
            C = np.pi if kap == 0 else (np.pi if m else 2 * np.pi) * 0.5 * (1 - m * m / (kap * kap)) * Jw * Jw
            lst = []
            for q in range(len(vals)):
                lst.append(len(modes))
                modes.append(dict(m=m, g=len(groups), k=float(np.sqrt(vals[q])), N=C * h, Z=(vecs[:, q] / aC).astype(F32)))
            groups.append(dict(m=m, n=n, kappa=kap, C=C, Jw=Jw, modes=lst))
        m += 1
    return dict(G=G, a0=a0, Hn=Hn, nz=nz, h=h, groups=groups, modes=modes, mMax=max(g["m"] for g in groups))


def resonance_G(Q):
    al = np.sqrt(max(1e-9, 4 - 1 / (Q * Q)))
    ha, be = al / 2, 1 / (2 * Q)
    be2 = be * be
    return lambda r: (np.log(((r + ha) ** 2 + be2) / ((r - ha) ** 2 + be2)) / (4 * al)
                      + (np.arctan((r + ha) / be) + np.arctan((r - ha) / be)) / (4 * be))


def cavity_drive(an, M, P, sr):
    s = K11 / max(1.0, float(P["ctune"]))
    Q = max(1.0, float(P["cq"]))
    N, nT, nB = an["N"], an["nT"], an["nBins"]
    df = sr / N
    spec = an["spec"].astype(np.float64)
    fTop = KMAX / s
    bmax = int(max(1, min(nB - 1, np.floor(1.6 * fTop / df))))
    band = spec[:, 1:bmax + 1]
    row = -1
    if P.get("cmode") == "moment":
        row = int(np.floor((float(P["ct"]) * sr - an["s0"]) / an["hop"] - 0.5 + 0.5))
        row = max(0, min(nT - 1, row))
        Pb = band[row]
    else:
        Pb = band.mean(0)
    G = resonance_G(Q)
    u = max(0.0, min(1.0, float(P["csrc"]))) * M["Hn"] / M["h"] - 0.5
    i0, tt, nz = int(np.floor(u)), 0.0, M["nz"]
    tt = u - i0
    if i0 < 0:
        i0, tt = 0, 0.0
    elif i0 > nz - 2:
        i0, tt = nz - 2, 1.0
    edges = np.arange(bmax + 1) + 0.5
    nM = len(M["modes"])
    f, D, phis2, w, ub = np.zeros(nM), np.zeros(nM), np.zeros(nM), np.zeros(nM), np.zeros(bmax)
    for j, md in enumerate(M["modes"]):
        fj = md["k"] / s
        f[j] = fj
        Z = md["Z"].astype(np.float64)
        zs = Z[i0] + (Z[i0 + 1] - Z[i0]) * tt
        Jw = M["groups"][md["g"]]["Jw"]
        phis2[j] = Jw * Jw * zs * zs
        cj = phis2[j] / (md["N"] ** 2 * md["k"] ** 4)
        wb = np.diff(G(edges * df / fj)) * (fj / df)
        D[j] = float((Pb * wb).sum())
        w[j] = D[j] * cj
        ub += wb * cj
    ring = band @ ub
    rm = ring.mean()
    gain = (ring[row] / rm if rm > 0 else 0.0) if row >= 0 else 1.0
    return dict(s=s, Q=Q, f=f, D=D, phis2=phis2, w=w, ring=ring, gain=gain, row=row, bmax=bmax, window=(f.min(), fTop))


def cavity_field(M, drive, g, P):
    """exact p² on the page's voxel grid, smoothed and normalised the same way; shape (nz, ny, nx)"""
    from scipy.special import jv
    from scipy.ndimage import gaussian_filter1d
    nx, ny, nz, vx, vz = g["nx"], g["ny"], g["nz"], g["vx"], g["vz"]
    G = M["G"]
    cr0, cr1, ch, a0, flr = G["cr0"], G["cr1"], G["ch"], M["a0"], g["flr"]
    X, Y = np.meshgrid(g["x0"] + np.arange(nx) * vx, g["y0"] + np.arange(ny) * vx)
    R, TH = np.hypot(X, Y), np.arctan2(Y, X)
    cs = max(0.0, float(P["csmooth"]))
    margin = g["tw"] + 3 * cs + 2 * vx
    cos2 = [np.cos(m * TH) ** 2 for m in range(M["mMax"] + 1)]
    Zs = np.stack([md["Z"].astype(np.float64) for md in M["modes"]])
    gOf = np.array([md["g"] for md in M["modes"]])
    w = drive["w"]
    p2 = np.zeros((nz, ny, nx))
    for k in range(nz):
        zz = min(ch, max(0.0, g["z0"] + k * vz - flr))
        a = cr0 + (cr1 - cr0) * zz / ch
        u = zz / a0 / M["h"] - 0.5
        i0 = int(np.floor(u))
        tz = u - i0
        if i0 < 0:
            i0, tz = 0, 0.0
        elif i0 > M["nz"] - 2:
            i0, tz = M["nz"] - 2, 1.0
        zv = Zs[:, i0] + (Zs[:, i0 + 1] - Zs[:, i0]) * tz
        zp = np.bincount(gOf, weights=w * zv * zv, minlength=len(M["groups"]))
        rho = np.minimum(1.0, R / a)
        acc = np.zeros((ny, nx))
        for gi, gr in enumerate(M["groups"]):
            if zp[gi] == 0:
                continue
            J2 = 1.0 if gr["kappa"] == 0 else jv(gr["m"], gr["kappa"] * rho) ** 2
            acc += zp[gi] * J2 * cos2[gr["m"]]
        p2[k] = np.where(R <= a + margin, acc, 0.0)
    if cs > 0:
        for axis, sig in ((2, cs / vx), (1, cs / vx), (0, cs / vz)):
            p2 = gaussian_filter1d(p2, sig, axis=axis, mode="nearest", truncate=3.0)
    zz = g["z0"] + np.arange(nz) * vz - flr
    a = cr0 + (cr1 - cr0) * zz / ch
    cav = (R[None, :, :] <= a[:, None, None]) & ((zz >= 0) & (zz <= ch))[:, None, None]
    mx = p2[cav].max()
    return p2 / mx, cav


def cavity_implicit(p2n, cav, g, P, level, half=False):
    """F > 0 is material, the page's formula; islands removed with scipy.ndimage.label (6-connected)"""
    from scipy.ndimage import label
    nx, ny, nz, vx, vz = g["nx"], g["ny"], g["nz"], g["vx"], g["vz"]
    G = cavity_geom(P)
    cr0, cr1, ch, flr, tw, zRim, zTop = G["cr0"], G["cr1"], G["ch"], g["flr"], g["tw"], g["zRim"], g["zTop"]
    closed = P.get("ctop") == "closed"
    Ls = max(0.5, float(P["cwall"]))
    x = g["x0"] + np.arange(nx) * vx
    y = g["y0"] + np.arange(ny) * vx
    z = g["z0"] + np.arange(nz) * vz
    Zg, Yg, Xg = np.meshgrid(z, y, x, indexing="ij")
    R = np.hypot(Xg, Yg).astype(F32).astype(np.float64)
    zc = np.clip(Zg - flr, 0, ch)
    a = cr0 + (cr1 - cr0) * zc / ch
    aO = a + tw
    zA = np.minimum(Zg, zTop - Zg)
    sh = np.minimum(np.minimum(R - a, aO - R), zA)
    sh = np.maximum(sh, np.minimum(np.minimum(Zg, flr - Zg), cr0 + tw - R))
    if closed:
        sh = np.maximum(sh, np.minimum(np.minimum(Zg - zRim, zTop - Zg), cr1 + tw - R))
    lobe = np.minimum((p2n / level - 1) * Ls, np.minimum(aO - R, zA))
    F = np.maximum(lobe, sh)
    if half:
        F = np.minimum(F, Yg)[:, g["c"]:, :]
    lab, n = label(F > 0)
    seed = lab[g["below"], 1 if half else g["c"], g["c"]]
    F = np.where((F > 0) & (lab != seed), -0.5 * vx, F)
    return F


def iso_volume(F, g, half=False):
    import pyvista as pv
    nz, ny, nx = F.shape
    y0 = g["y0"] + (g["c"] * g["vx"] if half else 0.0)
    img = pv.ImageData(dimensions=(nx, ny, nz), spacing=(g["vx"], g["vx"], g["vz"]), origin=(g["x0"], y0, g["z0"]))
    img.point_data["F"] = F.astype(np.float64).ravel(order="C")
    surf = img.contour([0.0], scalars="F")
    return float(surf.volume), surf.n_cells


def stl_volume(path):
    tri = read_stl(path).astype(np.float64)
    return float(np.einsum("ij,ij->i", tri[:, 0], np.cross(tri[:, 1], tri[:, 2])).sum() / 6.0)


def compare_cavity():
    path = os.path.join(HERE, "demo_cavity.json")
    if not os.path.exists(path):
        print("  (no demo_cavity.json: run node run_pure.mjs first)")
        return 1
    with open(path) as fh:
        js = json.load(fh)
    P = dict(defaults())
    P.update(js["P"])
    sr = js["sr"]
    fails = 0

    def ok(cond, msg):
        nonlocal fails
        print(("  ok   " if cond else "  FAIL ") + msg)
        if not cond:
            fails += 1

    x = synth_demo(sr, js["seconds"])
    an = stft(x, sr, P)
    M = cavity_modes(P)
    jm = js["modes"]
    ok(len(M["groups"]) == len(jm["groups"]) and len(M["modes"]) == jm["count"],
       f"cavity: {len(M['groups'])} groups, {len(M['modes'])} modes (js {len(jm['groups'])}, {jm['count']})")
    dk = max(abs(a["kappa"] - b["kappa"]) for a, b in zip(M["groups"], jm["groups"]))
    ok(dk <= 1e-9, f"cavity: zeros of J_m' vs scipy.special.jnp_zeros, max|d| {dk:.1e}")
    kk = np.array([md["k"] for md in M["modes"]])
    dr = np.abs(kk / np.array(jm["k"]) - 1).max() if len(kk) == jm["count"] else np.inf
    ok(dr <= 1e-6, f"cavity: axial eigenvalues vs LAPACK eigh_tridiagonal, max rel {dr:.1e}")
    drive = cavity_drive(an, M, P, sr)
    jd = js["drive"]
    wr = np.abs(drive["w"] - np.array(jd["w"])).max() / np.abs(drive["w"]).max()
    ok(wr <= 1e-5, f"cavity: modal weights (bin-integrated resonance, 1/w^4 speaker drive), max|d|/max {wr:.1e}")
    g = js["grid"]
    p2n, cav = cavity_field(M, drive, g, P)
    jp = np.fromfile(os.path.join(HERE, "demo_cavity_field.f32"), dtype="<f4").reshape(g["nz"], g["ny"], g["nx"])
    dp = np.abs(p2n - jp).max()
    ok(dp <= 1e-3, f"cavity: p^2 on {g['nx']}x{g['ny']}x{g['nz']} voxels, exact Bessel vs page tables, max|d| {dp:.1e} of 1")
    for run in js["runs"]:
        Q = dict(P)
        Q.update(run["P"])
        st = run["stats"]
        gq = run["grid"]
        if Q.get("cmode") == "moment" or Q.get("ctop") == "closed":
            Mq = cavity_modes(Q) if Q.get("ctop") == "closed" else M
            dq = cavity_drive(an, Mq, Q, sr)
            pq, cq = cavity_field(Mq, dq, gq, Q)
        else:
            dq, pq, cq = drive, p2n, cav
        fill = min(0.6, float(Q["cfill"]) / 100 * dq["gain"])
        thr = float(np.quantile(pq[cq], 1 - fill))
        ok(abs(thr - st["level"]) <= 2e-3, f"{run['name']}: fill {fill * 100:.1f} % threshold {thr:.4f} vs page {st['level']:.4f}")
        half = bool(int(Q.get("chalf", 0)))
        F = cavity_implicit(pq, cq, gq, Q, st["level"], half)
        v, nc = iso_volume(F, gq, half)
        v *= 2 if half else 1
        vj = stl_volume(os.path.join(HERE, run["name"]))
        ok(abs(v / vj - 1) <= 0.02, f"{run['name']}: page (marching tetrahedra) {vj / 1000:.2f} cm3 vs VTK marching cubes {v / 1000:.2f} cm3 ({(v / vj - 1) * 100:+.2f} %)")
    b = js["bare"]
    ok(abs(b["vol"] / b["exact"] - 1) <= 0.02, f"bare vessel {b['vol'] / 1000:.2f} vs closed form {b['exact'] / 1000:.2f} cm3")
    print("compare cavity:", "PASS" if not fails else f"{fails} FAILURE(S)")
    return 1 if fails else 0



def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--write-test-wav", metavar="PATH")
    ap.add_argument("--wav", metavar="PATH")
    ap.add_argument("--layout", default="mirror", choices=["mirror", "plain", "transpose"])
    ap.add_argument("--hollow", action="store_true")
    ap.add_argument("--nt", type=int, default=192)
    ap.add_argument("--nf", type=int, default=96)
    ap.add_argument("--compare", action="store_true")
    a = ap.parse_args()
    rc = 0
    if a.write_test_wav:
        write_test_wav(a.write_test_wav)
    if a.wav:
        sr, x = load_wav(a.wav)
        P = defaults()
        P.update(layout=a.layout, hollow=1 if a.hollow else 0, nt=a.nt, nf=a.nf)
        r = pipeline(x, sr, P)
        s = r["stats"]
        out = os.path.join(HERE, "py_" + os.path.splitext(os.path.basename(a.wav))[0] + f"_{a.layout}{'_hollow' if a.hollow else ''}.stl")
        write_stl(out, r["mesh"]["pos"], r["mesh"]["idx"])
        print(f"{os.path.basename(out)}: {sr} Hz {len(x) / sr:.1f} s -> V {s['V']} F {s['F']} Euler {s['euler']} odd {s['badEdges']} rev {s['badDir']} vol {s['vol'] / 1000:.1f} cm3 size {' x '.join(f'{v:.1f}' for v in s['size'])}")
        if not (s["euler"] == 2 and s["badEdges"] == 0 and s["badDir"] == 0 and s["vol"] > 0):
            rc = 1
    if a.compare:
        rc = max(rc, compare())
        rc = max(rc, compare_cavity())
    if not (a.write_test_wav or a.wav or a.compare):
        ap.print_help()
    return rc


if __name__ == "__main__":
    sys.exit(main())
