#!/usr/bin/env python3
"""2D layout traced from the reference photo (reference_photo.png, 360 x 340 px).

The silhouette comes from thresholding the photo: the white sculpture against a dark ground.
The features inside it are polylines in photo pixels, mirrored about the photo's axis (x = 179),
since the sculpture is symmetric. The hollow's upper part and the grooves are hand-traced off
6-8x zoomed crops. The hollow's lower walls, the canal and the lips' rounded ends are taken
row by row (and column by column) at the photo's steepest brightness edge, the canal's two
sides averaged:
    hollow  - the cavity of the body, from under the fundus band down to the internal orifice
    canal   - the slit between the lips: it pinches below the orifice, swells, pinches again at
              the lips' knees and opens round their rounded ends into the cup
    cup     - the U-shaped recess around the lips' rounded ends at the base
    groove  - the tube lumen, from the cavity's top corner out through each tube tab

Run alone to write photo_trace_check.png, the traced masks laid over the photo.
"""

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage as nd

AXIS = 179.0            # photo x of the symmetry axis (candle, tabs and neck all centre here)
BASE = 342              # photo y of the flat base (the photo is cropped at 340)
SIZE = (360, 342)       # mask size: photo width, rows down to the base

# Left halves; each is mirrored about AXIS.
HOLLOW_LEFT = [(179, 95), (160, 94), (140, 92), (120, 90), (100, 89), (92, 91), (90, 97),
               (100, 101), (110, 108), (118, 117), (124, 128), (128, 140), (132, 152),
               (135, 164), (139, 176), (143, 188), (147, 200), (150, 210), (155, 220),
               (158, 226), (161, 230), (163, 234), (165, 238), (168, 242), (169, 246),
               (170, 250), (174, 253), (179, 256)]
# (the hollow's tip, below y = 250, lies inside the canal, which takes over there)
# canal, left side, axis to axis: a narrow straight run up inside the hollow (hidden by the
# candle), then the lip's edge down to its rounded end, then straight down through the cup
CANAL_LEFT = [(179, 205), (170, 205), (170, 250), (171, 254), (172, 258), (172.5, 262),
              (172.5, 270), (171.5, 274), (171, 277), (170.5, 280), (170, 283), (169.5, 287),
              (169, 291), (169, 296), (169.5, 299), (170, 302), (170.5, 305), (171, 308),
              (171.5, 310), (172.5, 313), (173.5, 317), (173.5, 322), (173, 324), (172, 327),
              (170, 329.5), (167, 331.5), (164, 333), (164, BASE + 2), (179, BASE + 2)]
# cup, left half: outer rim from its top corner down to the base, then the lip's rounded end
# back up from the axis to the same corner
CUP_RIM = [(133, 311), (130, 316), (130, 322), (132, 326), (134, 331), (136, 336),
           (138, BASE + 2)]
CUP_LIP = [(179, 334), (164, 333.5), (160, 333), (156, 332), (152, 330), (149, 328),
           (146, 326), (144, 323.5), (142, 321), (140, 318.5), (138, 316), (136, 314)]
GROOVE_LEFT = [(40, 93), (60, 93), (80, 93), (95, 94), (110, 97)]   # tab end -> cavity corner
CANDLE_TOP = 222        # photo y of the candle's top (wick base)
ORIFICE = 266           # photo y where the lips close narrowest, below the hollow: the squeeze
CANAL_FLOOR = 318       # photo y the canal's candle channel ends (at the lips' knees)


def mirror(pts):
    return [(2 * AXIS - x, y) for x, y in pts]


def closed(left):
    """Left half (axis point first and last) + its mirror, as one closed outline."""
    right = mirror(left)[::-1]
    return left + right


def poly_mask(pts, size=SIZE, scale=4):
    """Rasterise a polygon at 4x and box-filter down, for sub-pixel edges (0..1 coverage)."""
    W, H = size
    img = Image.new("L", (W * scale, H * scale), 0)
    ImageDraw.Draw(img).polygon([(x * scale, y * scale) for x, y in pts], fill=255)
    a = np.asarray(img, float).reshape(H, scale, W, scale).mean(axis=(1, 3)) / 255
    return a


def silhouette():
    g = np.asarray(Image.open("reference_photo.png").convert("L"), float)
    sil = nd.binary_fill_holes(g > 100)
    sil = nd.binary_opening(sil, iterations=2)
    lab, n = nd.label(sil)
    sil = lab == 1 + np.argmax(nd.sum(sil, lab, range(1, n + 1)))
    sil = np.vstack([sil[:338], np.repeat(sil[337:338], SIZE[1] - 338, axis=0)])
    # symmetrise: a pixel is in when the average of it and its mirror is
    cols = np.arange(SIZE[0])
    mx = np.clip(np.round(2 * AXIS - cols).astype(int), 0, SIZE[0] - 1)
    return (0.5 * (sil + sil[:, mx]) >= 0.5)


def smooth_outline(pts, step=0.5, sigma=3.0):
    """A closed polyline resampled every `step` px and smoothed along its length (sigma in
    samples, so corners round off at about 1.5 px), for a distance field free of facets."""
    p = np.array(pts, float)
    p = p[np.r_[True, np.any(np.diff(p, axis=0) != 0, axis=1)]]
    p = np.vstack([p, p[:1]])
    s = np.r_[0, np.cumsum(np.hypot(*np.diff(p, axis=0).T))]
    t = np.arange(0, s[-1], step)
    q = np.c_[np.interp(t, s, p[:, 0]), np.interp(t, s, p[:, 1])]
    return nd.gaussian_filter1d(q, sigma, axis=0, mode="wrap")


def masks():
    """(silhouette, cavity) boolean masks in photo pixels, plus the lumen groove polylines."""
    sil = silhouette()
    hollow = poly_mask(closed(HOLLOW_LEFT)) > 0.5
    cup_outline = CUP_RIM + mirror(CUP_RIM)[::-1] + mirror(CUP_LIP)[::-1] + CUP_LIP[1:]
    cup = poly_mask(cup_outline) > 0.5
    cav = (hollow | cup) & sil
    grooves = [np.array(GROOVE_LEFT, float), np.array(mirror(GROOVE_LEFT), float)]
    return sil, cav, grooves


def canal_outline():
    """The canal as a smooth closed outline in photo pixels. It is kept apart from the cavity
    and not rasterised: it is cut in 3D, back to the candle, from its exact distance field."""
    return smooth_outline(closed(CANAL_LEFT))


if __name__ == "__main__":
    sil, cav, grooves = masks()
    im = Image.open("reference_photo.png").convert("RGB").crop((0, 0) + SIZE)
    a = np.asarray(im, float)
    edge = sil ^ nd.binary_erosion(sil)
    cedge = cav ^ nd.binary_erosion(cav)
    a[edge] = (0, 255, 0)
    a[cedge] = (255, 0, 255)
    out = Image.fromarray(a.astype(np.uint8))
    d = ImageDraw.Draw(out)
    d.line([tuple(p) for p in canal_outline()], fill=(255, 0, 255), width=1)
    for g in grooves:
        d.line([tuple(p) for p in g], fill=(0, 200, 255), width=1)
    out.resize((SIZE[0] * 3, SIZE[1] * 3), Image.NEAREST).save("photo_trace_check.png")
    ys, xs = np.nonzero(sil)
    print("silhouette x", xs.min(), xs.max(), "y", ys.min(), ys.max())
