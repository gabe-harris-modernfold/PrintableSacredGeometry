#!/usr/bin/env python3
"""Sync the web page and the staged site files with the current design.

1. Re-solves the voices (voices.py) and rewrites voices.json, keeping the layout block
   that chandelier.py last wrote (run chandelier.py first after any geometry change).
2. Injects a trimmed copy of voices.json into web/hymn-of-the-ogdoad.html.
3. Stages the site: page copy, preview webps in site/assets/ogdoad/, the STL zip in
   site/stl/, and the menu-card preview in site/previews/.

Run from this directory:  python sync_page.py
"""
import json
import os
import re
import shutil
import zipfile

from PIL import Image

import voices as vc

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
PAGE = os.path.join(ROOT, "web", "hymn-of-the-ogdoad.html")
SITE = os.path.join(ROOT, "site")


def autocrop(im, margin=36):
    """Trim the render's flat background to the model plus a margin."""
    import numpy as np
    a = np.asarray(im).astype(int)
    bg = a[2, 2]
    mask = np.abs(a - bg).sum(axis=2) > 12
    ys, xs = np.where(mask)
    if not len(xs):
        return im
    x0, x1 = max(0, xs.min() - margin), min(im.width, xs.max() + margin)
    y0, y1 = max(0, ys.min() - margin), min(im.height, ys.max() + margin)
    return im.crop((x0, y0, x1, y1))


def main():
    with open(os.path.join(HERE, "voices.json")) as fh:
        layout = json.load(fh).get("layout")
    assert layout, "voices.json has no layout: run chandelier.py first"
    vc.export_json(vc.design(), layout)
    with open(os.path.join(HERE, "voices.json")) as fh:
        data = json.load(fh)
    for v in data["voices"]:
        o = v["outer"]
        v["outer"] = [[round(r, 1), round(z, 1)] for r, z in o[::2]] + ([o[-1]] if len(o) % 2 == 0 else [])
    blob = json.dumps(data, separators=(",", ":"))
    with open(PAGE, encoding="utf-8") as fh:
        page = fh.read()
    page, n = re.subn(r'(<script id="voices-data" type="application/json">).*?(</script>)',
                      lambda m: m.group(1) + blob + m.group(2), page, flags=re.S)
    assert n == 1
    with open(PAGE, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(page)
    print(f"injected {len(blob) / 1024:.1f} kB of data into {os.path.relpath(PAGE, ROOT)}")

    shutil.copy(PAGE, os.path.join(SITE, "hymn-of-the-ogdoad.html"))
    ad = os.path.join(SITE, "assets", "ogdoad")
    os.makedirs(ad, exist_ok=True)
    for src, dst, w in (("preview_three_quarter.png", "three-quarter.webp", 1400),
                        ("preview_front.png", "front.webp", 900),
                        ("preview_side.png", "side.webp", 900),
                        ("preview_head_section.png", "head-section.webp", 1200)):
        im = autocrop(Image.open(os.path.join(HERE, src)).convert("RGB"))
        im.thumbnail((w, w))
        im.save(os.path.join(ad, dst), "WEBP", quality=86)
    sd = os.path.join(SITE, "stl")
    os.makedirs(sd, exist_ok=True)
    zp = os.path.join(sd, "ogdoad_wind_chandelier.zip")
    with zipfile.ZipFile(zp, "w", zipfile.ZIP_DEFLATED) as z:
        for f in sorted(os.listdir(os.path.join(HERE, "stl"))):
            z.write(os.path.join(HERE, "stl", f), f)
        z.write(os.path.join(HERE, "README.md"), "README.md") if os.path.exists(os.path.join(HERE, "README.md")) else None
    print(f"staged page, 4 webps and {os.path.relpath(zp, ROOT)} ({os.path.getsize(zp) / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
