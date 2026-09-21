"""Turn the thermal scan into (a) an RGBA body sprite and (b) a lattice collision mask.

    python build_body_asset.py

The exhibit draws a real anterior thermal photograph instead of a drawn silhouette, so
the solver's solid mask has to come from the *same* outline or smoke will visibly clip
through the body. This script derives both from one source image.

What it does NOT do is read temperature out of the photograph. The scan's palette is
uncalibrated false colour for one individual; the physics keeps the measured
population values in body_model.SEGMENTS. The image supplies shape and appearance
only, and the page says so.

Note on the source: isolating the body as the largest connected component also drops
the "MEDICAL INFRARED" watermark that sits on the black background. Confirm the image
is cleared for publication before deploying the site.
"""

import base64
import io
import json
import os
import sys

import numpy as np
from PIL import Image
from scipy import ndimage

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import body_model as bm

SRC = "Thermal-Infrared-Full-Body-Scan.jpg"
SPRITE_OUT = "body_thermal.png"
SITE_SPRITE = os.path.join("..", "..", "site", "assets", "plume", "body_thermal.png")
WEB_SNIPPET = "body_mask.js"

# Lattice the exhibit runs on. Must match the constants in the HTML.
NX, NY = 130, 228
DOMAIN_W = 1.90
DXM = DOMAIN_W / NX
BODY_Y = 2
LUMA_MIN = 30           # background is pure black; anything above this is subject
SPRITE_EMBED_W = 300    # px; the canvas never draws it wider than about 260

# Segment lookup by name, so the ids we emit index body_model.SEGMENTS.
SEG_IX = {s["name"]: i for i, s in enumerate(bm.SEGMENTS)}


def isolate_body(img):
    """Largest connected non-black component, holes filled."""
    a = np.asarray(img.convert("RGB")).astype(np.int32)
    luma = (0.299 * a[..., 0] + 0.587 * a[..., 1] + 0.114 * a[..., 2])
    raw = luma > LUMA_MIN
    raw = ndimage.binary_opening(raw, np.ones((3, 3)))
    lab, n = ndimage.label(raw)
    if n == 0:
        raise RuntimeError("no subject found in the scan")
    sizes = ndimage.sum(raw, lab, range(1, n + 1))
    keep = int(np.argmax(sizes)) + 1
    mask = ndimage.binary_fill_holes(lab == keep)
    dropped = int(raw.sum() - mask.sum())
    return mask, dropped, n


def classify(mask_lat, l_body, body_x):
    """Assign a measured segment to every solid lattice cell.

    Height alone will not do it: in this pose the hands hang at about half body
    height, level with the thighs, so a purely vertical rule would label them thigh.
    Each row is split into contiguous runs instead, and a run far from the midline is
    an arm or a hand whatever its height.
    """
    seg = np.full(mask_lat.shape, -1, dtype=np.int16)
    ny, nx = mask_lat.shape
    arm_offset_cells = 0.20 / DXM

    for j in range(ny):
        row = mask_lat[j]
        if not row.any():
            continue
        hf = (j - BODY_Y) / l_body          # 0 at the floor, 1 at the crown
        idx = np.nonzero(row)[0]
        splits = np.nonzero(np.diff(idx) > 1)[0]
        runs = np.split(idx, splits + 1)
        for run in runs:
            centre = run.mean()
            off = abs(centre - body_x)
            if off > arm_offset_cells:
                # Measured off the scan: fingertips sit at hf 0.45-0.57 and the
                # forearm runs 0.55-0.66. Above the armpit (hf ~0.68) the arm is
                # contiguous with the torso, so no separate upper-arm run exists in
                # this pose and those cells take the chest value.
                name = ("hand" if hf < 0.54 else
                        "forearm" if hf < 0.66 else "upper_arm")
            elif hf < 0.050:
                name = "foot"
            elif hf < 0.300:
                name = "lower_leg"
            elif hf < 0.545:
                name = "thigh"
            elif hf < 0.595:
                name = "pelvis"
            elif hf < 0.647:
                name = "abdomen"
            elif hf < 0.841:
                name = "chest"
            elif hf < 0.909:
                name = "neck"
            else:
                name = "head"
            seg[j, run] = SEG_IX[name]
    return seg


def thermal_ramp(rgb, mask, n_stops=9):
    """Recover the scan's own false-colour ramp, coldest to hottest.

    The palette is a rainbow map, so hue runs monotonically from blue through
    cyan, green and yellow to red as temperature rises, with the very hottest
    pixels desaturating toward white. Ordering the body's pixels by that path and
    sampling it at even quantiles reproduces the photograph's colouring rather
    than a generic rainbow that merely resembles it.
    """
    px = rgb[mask].astype(np.float32) / 255.0
    mx = px.max(axis=1)
    mn = px.min(axis=1)
    delta = np.maximum(mx - mn, 1e-6)

    r, g, b = px[:, 0], px[:, 1], px[:, 2]
    hue = np.zeros(len(px), dtype=np.float32)
    m = (mx == r)
    hue[m] = ((g[m] - b[m]) / delta[m]) % 6
    m = (mx == g) & (mx != r)
    hue[m] = (b[m] - r[m]) / delta[m] + 2
    m = (mx == b) & (mx != r) & (mx != g)
    hue[m] = (r[m] - g[m]) / delta[m] + 4
    hue *= 60.0                                   # degrees

    # Rank cold -> hot: blue (240) down through green (120) to red (0), then the
    # washed-out hot spots, which have lost saturation and would otherwise sort
    # back to the cold end.
    sat = delta / np.maximum(mx, 1e-6)
    score = np.where(sat < 0.25, 400.0 + mx * 100.0, 300.0 - hue)
    keep = mx > 0.12
    score, px = score[keep], px[keep]
    order = np.argsort(score)
    px = px[order]

    stops = []
    for q in np.linspace(0.02, 0.98, n_stops):
        lo = int(max(0, q * len(px) - 0.005 * len(px)))
        hi = int(min(len(px), q * len(px) + 0.005 * len(px))) or 1
        stops.append([round(float(c), 4) for c in px[lo:hi].mean(axis=0)])
    return stops


def segment_anchors(seg, l_body, body_x):
    """A label anchor per segment: its vertical centre and the outer edge there."""
    out = {}
    for k in range(len(bm.SEGMENTS)):
        jj, ii = np.nonzero(seg == k)
        if jj.size == 0:
            continue
        jc = int(np.median(jj))
        row = np.nonzero(seg[jc] == k)[0]
        if row.size == 0:
            continue
        left = int(row.min())
        right = int(row.max())
        # Put the label on whichever side that segment actually reaches further.
        side = 1 if abs(right - body_x) >= abs(left - body_x) else -1
        out[bm.SEGMENTS[k]["name"]] = dict(
            row=jc, col=(right if side > 0 else left), side=side,
            c=bm.SEGMENTS[k]["T_skin"],
            f=round(bm.SEGMENTS[k]["T_skin"] * 9 / 5 + 32, 1))
    return out


def main():
    img = Image.open(SRC)
    mask, dropped, ncomp = isolate_body(img)
    ys, xs = np.nonzero(mask)
    y0, y1, x0, x1 = ys.min(), ys.max() + 1, xs.min(), xs.max() + 1
    print(f"source {img.size}, {ncomp} components, kept the largest; "
          f"dropped {dropped} px outside it (watermark and speckle)")
    print(f"body bounding box x[{x0}:{x1}] y[{y0}:{y1}] -> {x1 - x0}x{y1 - y0}")

    rgb = np.asarray(img.convert("RGB"))[y0:y1, x0:x1]
    m = mask[y0:y1, x0:x1]
    rgba = np.dstack([rgb, (m * 255).astype(np.uint8)])
    sprite = Image.fromarray(rgba, "RGBA")
    sprite.save(SPRITE_OUT)
    os.makedirs(os.path.dirname(SITE_SPRITE), exist_ok=True)
    sprite.save(SITE_SPRITE)

    ih, iw = m.shape
    aspect = iw / ih
    body_w_m = bm.BODY_HEIGHT * aspect
    l_body = bm.BODY_HEIGHT / DXM
    w_cells = body_w_m / DXM
    body_x = NX // 2
    print(f"aspect {aspect:.3f} -> body spans {body_w_m:.3f} m "
          f"({w_cells:.1f} cells) for a {bm.BODY_HEIGHT} m figure")
    print(f"lattice {NX}x{NY}, dx {DXM * 1000:.2f} mm, body {l_body:.1f} cells tall")
    clear = (0.5 * DOMAIN_W - 0.5 * body_w_m)
    print(f"side clearance beyond the fingertips: {clear:.3f} m")

    # Resample the mask onto the lattice. Area-average then threshold, so a cell
    # counts as solid when the body covers at least a third of it.
    col0 = body_x - w_cells / 2.0
    lat = np.zeros((NY, NX), dtype=bool)
    for j in range(NY):
        fy = (j - BODY_Y) / l_body
        if not (0.0 <= fy < 1.0):
            continue
        sy0 = int((1.0 - fy - 1.0 / l_body) * ih)
        sy1 = max(sy0 + 1, int((1.0 - fy) * ih))
        for i in range(NX):
            fx = (i - col0) / w_cells
            if not (0.0 <= fx < 1.0):
                continue
            sx0 = int(fx * iw)
            sx1 = max(sx0 + 1, int((fx + 1.0 / w_cells) * iw))
            patch = m[sy0:sy1, sx0:sx1]
            if patch.size and patch.mean() >= 0.34:
                lat[j, i] = True

    seg = classify(lat, l_body, body_x)
    lat = seg >= 0
    print(f"solid lattice cells: {lat.sum()}")

    jj, ii = np.nonzero(lat)
    bj0, bj1, bi0, bi1 = jj.min(), jj.max() + 1, ii.min(), ii.max() + 1
    box = (seg[bj0:bj1, bi0:bi1] + 1).astype(np.uint8)   # 0 = fluid, 1.. = segment+1
    payload = base64.b64encode(box.tobytes()).decode("ascii")
    print(f"mask box {bi1 - bi0}x{bj1 - bj0} = {box.size} bytes "
          f"-> {len(payload)} base64 chars")

    counts = {bm.SEGMENTS[k]["name"]: int((seg == k).sum())
              for k in range(len(bm.SEGMENTS)) if (seg == k).any()}
    print("cells per segment:", json.dumps(counts))

    theta = [round((s["T_skin"] - bm.T_AMBIENT) / bm.DELTA_T_REF, 5)
             for s in bm.SEGMENTS]

    # Embed the sprite in the page rather than linking it. The visualisers are meant
    # to be self-contained, and a linked asset silently falls back to flat shading
    # whenever the page is opened from anywhere but the publish directory.
    sprite_small = sprite.resize(
        (SPRITE_EMBED_W, max(1, round(SPRITE_EMBED_W * ih / iw))), Image.LANCZOS)
    buf = io.BytesIO()
    sprite_small.save(buf, "WEBP", quality=86, method=6, exact=True)
    data_uri = "data:image/webp;base64," + base64.b64encode(buf.getvalue()).decode()
    print(f"embedded sprite {sprite_small.size} "
          f"= {len(buf.getvalue()) / 1024:.0f} kB webp -> {len(data_uri) / 1024:.0f} kB base64")

    ramp = thermal_ramp(rgb, m)
    anchors = segment_anchors(seg, l_body, body_x)
    print(f"palette: {len(ramp)} stops, coldest {ramp[0]} -> hottest {ramp[-1]}")
    print(f"label anchors for {len(anchors)} segments")

    with open(WEB_SNIPPET, "w") as fh:
        fh.write("// Generated by build_body_asset.py -- do not edit by hand.\n")
        fh.write("// Body outline traced from the anterior thermal scan, so the solver's\n")
        fh.write("// collision mask and the drawn sprite share one silhouette.\n")
        fh.write(f"const MASK_X0 = {bi0}, MASK_Y0 = {bj0};\n")
        fh.write(f"const MASK_W = {bi1 - bi0}, MASK_H = {bj1 - bj0};\n")
        fh.write(f"const SPRITE_X0 = {col0:.4f}, SPRITE_W = {w_cells:.4f};\n")
        fh.write(f"const SPRITE_Y0 = {BODY_Y}, SPRITE_H = {l_body:.4f};\n")
        fh.write(f"const SEG_THETA = {json.dumps(theta)};\n")
        fh.write(f"const SEG_NAMES = {json.dumps([s['name'] for s in bm.SEGMENTS])};\n")
        fh.write("// The scan's own false-colour ramp, coldest first, as 0-1 RGB.\n")
        fh.write(f"const THERMAL_RAMP = {json.dumps(ramp)};\n")
        fh.write("// Where to hang each segment's temperature label, in lattice cells.\n")
        fh.write(f"const SEG_ANCHORS = {json.dumps(anchors)};\n")
        fh.write(f"const T_AMBIENT_C = {bm.T_AMBIENT}, DELTA_T_REF_K = {bm.DELTA_T_REF};\n")
        fh.write("// The scan itself, embedded so the page needs no external asset and\n")
        fh.write("// renders the same wherever it is opened from.\n")
        fh.write(f'const SPRITE_DATA = "{data_uri}";\n')
        fh.write(f'const MASK_B64 = "{payload}";\n')

    print(f"wrote {SPRITE_OUT}, {SITE_SPRITE}, {WEB_SNIPPET}")


if __name__ == "__main__":
    main()
