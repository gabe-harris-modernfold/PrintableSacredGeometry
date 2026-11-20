"""Check The Shadow Solids page and write its Markdown reference.

Run from this directory:  python check_solids.py

The page web/the-shadow-solids.html carries all its content in one JSON block
(<script id="ss-data">). This script
  1. recomputes every geometric number the page states (faces, vertices, edges,
     duals, share of the circumsphere filled, sphericity) and checks the data,
  2. checks the ten hex colours and reports each two-state shift in HSL,
  3. scans every string for the excluded vocabulary (older correspondence
     systems, planets, classical elements, metals, god names, classical incense),
  4. checks the 2006-2026 window on the cultural objects,
  5. writes the-shadow-solids.md here, the repo copy of the reference.
Exit code is non-zero if anything fails.
"""
import colorsys
import json
import math
import os
import re
import sys

import numpy as np
from scipy.spatial import ConvexHull

HERE = os.path.dirname(os.path.abspath(__file__))
PAGE = os.path.join(HERE, "..", "..", "web", "the-shadow-solids.html")
OUT = os.path.join(HERE, "the-shadow-solids.md")

failures = []


def fail(msg):
    failures.append(msg)
    print("  FAIL", msg)


def load_data():
    html = open(PAGE, encoding="utf-8").read()
    m = re.search(r'<script type="application/json" id="ss-data">(.*?)</script>', html, re.S)
    if not m:
        sys.exit("no ss-data block found")
    return json.loads(m.group(1))


# ---------------------------------------------------------------- geometry
PHI = (1 + 5 ** 0.5) / 2


def verts(kind):
    pm = (-1, 1)
    if kind == "tetra":
        v = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
    elif kind == "octa":
        v = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
    elif kind == "cube":
        v = [(a, b, c) for a in pm for b in pm for c in pm]
    elif kind == "icosa":
        v = []
        for a in pm:
            for b in pm:
                v += [(0, a, b * PHI), (a, b * PHI, 0), (b * PHI, 0, a)]
    elif kind == "dodeca":
        v = [(a, b, c) for a in pm for b in pm for c in pm]
        for a in pm:
            for b in pm:
                v += [(0, a / PHI, b * PHI), (a / PHI, b * PHI, 0), (b * PHI, 0, a / PHI)]
    v = np.array(v, float)
    return v / np.linalg.norm(v, axis=1)[:, None]  # unit circumradius


def solid_numbers(kind):
    v = verts(kind)
    d = np.linalg.norm(v[:, None] - v[None], axis=2)
    dmin = d[d > 1e-9].min()
    edges = int((np.abs(d - dmin) < 1e-6).sum() // 2)
    hull = ConvexHull(v)
    # merge coplanar hull triangles into faces
    planes = {tuple(np.round(eq, 6)) for eq in hull.equations}
    faces = len(planes)
    vol, area = hull.volume, hull.area
    return {
        "vertices": len(v),
        "edges": edges,
        "faces": faces,
        "sphereFill": vol / (4 / 3 * math.pi),
        "sphericity": math.pi ** (1 / 3) * (6 * vol) ** (2 / 3) / area,
    }


DUALS = {"tetra": "itself", "cube": "octahedron", "octa": "cube", "icosa": "dodecahedron", "dodeca": "icosahedron"}
FACE_SIDES = {"triangle": 3, "square": 4, "pentagon": 5}


def check_geometry(solids):
    print("Geometry")
    edges_seq = []
    for s in solids:
        n = solid_numbers(s["shape"])
        edges_seq.append(n["edges"])
        for key in ("vertices", "edges", "faces"):
            if s[key] != n[key]:
                fail(f"{s['solid']}: {key} {s[key]} != computed {n[key]}")
        for key in ("sphereFill", "sphericity"):
            if abs(s[key] - n[key]) > 0.0015:
                fail(f"{s['solid']}: {key} {s[key]} != computed {n[key]:.4f}")
        if s["dual"] != DUALS[s["shape"]]:
            fail(f"{s['solid']}: dual {s['dual']} should be {DUALS[s['shape']]}")
        # Euler and face-shape consistency: F * sides = 2E
        if n["vertices"] - n["edges"] + n["faces"] != 2:
            fail(f"{s['solid']}: Euler characteristic")
        if s["faces"] * FACE_SIDES[s["face"]] != 2 * s["edges"]:
            fail(f"{s['solid']}: {s['faces']} {s['face']}s cannot give {s['edges']} edges")
        print(f"  {s['solid']:<20} V{n['vertices']:>3} E{n['edges']:>3} F{n['faces']:>3}"
              f"  fills {n['sphereFill']:.4f}  sphericity {n['sphericity']:.4f}")
    if edges_seq != sorted(edges_seq):
        fail(f"edge counts do not climb: {edges_seq}")
    # the claims made in prose on the page
    by = {s["shape"]: solid_numbers(s["shape"]) for s in solids}
    if max(by, key=lambda k: by[k]["sphericity"]) != "icosa":
        fail("icosahedron is not the most spherical")
    if max(by, key=lambda k: by[k]["sphereFill"]) != "dodeca":
        fail("dodecahedron does not fill the most of its sphere")
    if not by["octa"]["sphereFill"] < by["cube"]["sphereFill"]:
        fail("no dip from cube to octahedron")
    # pentagon diagonal / side = phi
    ang = [2 * math.pi * k / 5 for k in range(5)]
    p = [(math.cos(a), math.sin(a)) for a in ang]
    ratio = math.dist(p[0], p[2]) / math.dist(p[0], p[1])
    if abs(ratio - PHI) > 1e-12:
        fail("pentagon diagonal/side is not phi")


# ---------------------------------------------------------------- colour
def hsl(hexs):
    r, g, b = (int(hexs[i:i + 2], 16) / 255 for i in (1, 3, 5))
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    return h * 360, s * 100, l * 100


def check_colours(solids):
    print("Colour shifts (H deg, S %, L %)")
    seen = set()
    for s in solids:
        for key in ("colorInt", "colorSha"):
            hx = s[key]["hex"]
            if not re.fullmatch(r"#[0-9A-F]{6}", hx):
                fail(f"{s['solid']}: bad hex {hx}")
            if hx in seen:
                fail(f"duplicate colour {hx}")
            seen.add(hx)
        a, b = hsl(s["colorInt"]["hex"]), hsl(s["colorSha"]["hex"])
        print(f"  {s['solid']:<20} {s['colorInt']['hex']} H{a[0]:5.0f} S{a[1]:5.1f} L{a[2]:5.1f}"
              f"  ->  {s['colorSha']['hex']} H{b[0]:5.0f} S{b[1]:5.1f} L{b[2]:5.1f}")
    return {s["id"]: (hsl(s["colorInt"]["hex"]), hsl(s["colorSha"]["hex"])) for s in solids}


def check_shift_claims(h):
    """The colourShift sentences make quantitative claims; hold them to it."""
    (_, s0, l0), (_, s1, l1) = h["tetrahedron"]
    if not (s1 < s0 / 3 * 1.15 and l1 > l0):
        fail("tetrahedron: 'saturation falls by two-thirds and lightness rises' is not true")
    (_, s0, l0), (_, s1, l1) = h["cube"]
    if not (s1 < 8 and abs(l1 - l0) < 15):
        fail("cube: 'saturation drained almost to zero at nearly the same lightness' is not true")
    (h0, s0, _), (h1, s1, _) = h["octahedron"]
    if not (s1 > 99 and s1 > s0 and 300 < h1 < 345):
        fail("octahedron: 'saturation to maximum, hue toward magenta' is not true")
    (h0, _, l0), (h1, _, l1) = h["icosahedron"]
    if not (l1 > 85 and 195 < h1 < 230 and 170 < h0 < 185):
        fail("icosahedron: 'lightness near white, hue cooled from teal to blue' is not true")
    (_, _, l0), (h1, s1, l1) = h["dodecahedron"]
    if not (l1 < 20 and 240 < h1 < 290):
        fail("dodecahedron: 'near black with a trace of violet' is not true")


# ---------------------------------------------------------------- exclusions
# Whole words that belong to the excluded systems. Each hit is listed; a word in
# ALLOWED is permitted only inside the given phrase (a modern, literal use).
EXCLUDED = [
    # older correspondence systems
    "kabbalah", "qabalah", "cabala", "kabbalistic", "sephirot", "sefirot", "sephiroth", "sephira",
    "golden dawn", "hermetic", "hermes", "trismegistus", "tree of life", "tarot", "zodiac",
    # Hebrew letters
    "aleph", "beth", "gimel", "daleth", "vav", "zayin", "cheth", "teth", "yod", "kaph", "lamed",
    "mem", "nun", "samekh", "ayin", "tzaddi", "qoph", "resh", "shin", "tau",
    # classical planets and luminaries as attributions
    "saturn", "jupiter", "mars", "venus", "mercury", "luna", "sol", "planet", "planetary",
    # classical elements
    "fire", "water", "air", "earth", "aether", "ether", "quintessence", "elemental",
    # planetary metals
    "lead", "tin", "iron", "gold", "copper", "quicksilver", "silver",
    # god names
    "zeus", "apollo", "athena", "aphrodite", "ares", "hera", "artemis", "dionysus", "hades",
    "poseidon", "kronos", "cronus", "thoth", "isis", "osiris", "horus", "odin", "loki",
    # classical / renaissance incense and materia
    "frankincense", "myrrh", "olibanum", "kyphi", "benzoin", "dragon's blood", "storax", "galbanum",
]
ALLOWED = {
    "air": ["sea air", "filtered air", "radiator air", "air handling", "air freshener"],
}


def strings(obj):
    if isinstance(obj, dict):
        for v in obj.values():
            yield from strings(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from strings(v)
    elif isinstance(obj, str):
        yield obj


def check_exclusions(data, html_text):
    print("Exclusions")
    texts = list(strings(data))
    # prose outside the data block too (header, closing note)
    texts.append(re.sub(r"<script.*?</script>|<style.*?</style>|<[^>]+>", " ", html_text, flags=re.S))
    hits = 0
    for t in texts:
        low = t.lower()
        for w in EXCLUDED:
            for m in re.finditer(r"(?<![a-z])" + re.escape(w) + r"(?![a-z])", low):
                ctx = low[max(0, m.start() - 20): m.end() + 20]
                if any(p in ctx for p in ALLOWED.get(w, [])):
                    continue
                hits += 1
                fail(f"excluded term '{w}' in: ...{ctx}...")
    if not hits:
        print("  none of the excluded vocabulary appears")


def check_objects(data):
    print("Cultural objects")
    for s in data["solids"]:
        for o in s["objects"]:
            years = [int(y) for y in re.findall(r"\b(19\d\d|20\d\d)\b", o["year"])]
            if not years or min(years) < 2006 or max(years) > 2026:
                fail(f"{o['name']}: year '{o['year']}' outside 2006-2026")
            print(f"  {o['name']:<45} {o['year']}")


def check_shape(data):
    print("Structure")
    names = [s["shape"] for s in data["solids"]]
    if names != ["tetra", "cube", "octa", "icosa", "dodeca"]:
        fail(f"order is {names}")
    for s in data["solids"]:
        if not 2 <= len(s["integrated"]) <= 3 or not 2 <= len(s["shadowTraits"]) <= 3:
            fail(f"{s['solid']}: trait counts")
        if set(s["senses"]) != {"sight", "sound", "smell", "taste", "touch"}:
            fail(f"{s['solid']}: senses {sorted(s['senses'])}")
        if not 1 <= len(s["objects"]) <= 2:
            fail(f"{s['solid']}: object count")
        sentences = len(re.findall(r"[.!?](\s|$)", s["transmutation"]))
        if sentences > 2:
            fail(f"{s['solid']}: transmutation has {sentences} sentences")
    print("  order, trait counts, five senses, 1-2 objects, 1-2 sentence transmutation")


# ---------------------------------------------------------------- markdown
TOK = re.compile(r"\{\{([^|}]+)\|([^}]+)\}\}")
SENSES = [("sight", "Sight"), ("sound", "Sound"), ("smell", "Smell"), ("taste", "Taste"), ("touch", "Touch")]


def plain(s):
    return TOK.sub(lambda m: f"{m.group(1)} [{m.group(2)}]", s)


def cell(s):
    return plain(s).replace("|", "\\|")


def dual_txt(s):
    return "self-dual" if s["dual"] == "itself" else f"dual of the {s['dual']}"


def arch_full(s):
    return f"{s['archetype']} ({s['archetypeTerm']})" if s.get("archetypeTerm") else s["archetype"]


def sha_full(s):
    return f"{s['shadow']} ({s['shadowTerm']})" if s.get("shadowTerm") else s["shadow"]


def write_markdown(data):
    sol = data["solids"]
    head = ["Platonic Solid & Geometry", "Jungian Archetype & Polarity", "Positive Traits (Integrated)",
            "Negative Traits (Shadow)", "Core Drive & Fear", "Somatic Cue", "Modern Color Theory (2-State Shift)",
            "5-Sense Modern Correspondences", "Cultural Object (2006–2026)", "Transmutation / Magical Function"]
    lines = [
        "# The Shadow Solids — a modern table of correspondences",
        "",
        "*Generated by `check_solids.py` from the data block in `web/the-shadow-solids.html`. Edit the page, not this file.*",
        "",
        "Five Platonic solids, five Jungian archetypes, each with an integrated face and a shadow face. "
        "Built without inherited correspondence tables: the qualities are read from the geometry, the anchors "
        "come from 2006–2026.",
        "",
        "## Master table",
        "",
        "| " + " | ".join(head) + " |",
        "|" + "|".join("---" for _ in head) + "|",
    ]
    for s in sol:
        row = [
            f"**{s['solid']}**<br>{s['faces']} {s['face']}s · {s['vertices']} vertices · {s['edges']} edges · {dual_txt(s)}<br>*{s['quality']}*",
            f"**{arch_full(s)}**<br>↔ *{sha_full(s)}*",
            "<br>".join(s["integrated"]),
            "<br>".join(s["shadowTraits"]),
            f"**Drive:** {s['drive']}<br>**Fear:** {s['fear']}",
            f"**Integrated:** {s['somaticInt']}<br>**Shadow:** {s['somaticSha']}",
            f"**Balanced** `{s['colorInt']['hex']}` {s['colorInt']['name']}<br>**Shadow** `{s['colorSha']['hex']}` "
            f"{s['colorSha']['name']}<br>{s['colorEffect']}",
            "<br>".join(f"**{l}:** {s['senses'][k]['k']}" for k, l in SENSES),
            "<br>".join(f"{o['name']} ({o['year']})" for o in s["objects"]),
            s["transmutation"],
        ]
        lines.append("| " + " | ".join(cell(c) for c in row) + " |")
    lines += ["", "## Sensory & ritual breakdown", ""]
    for s in sol:
        lines += [
            f"### {s['solid']} — {plain(arch_full(s))} ↔ {plain(sha_full(s))}",
            "",
            f"*{s['faces']} {s['face']}s · {s['vertices']} vertices · {s['edges']} edges · {dual_txt(s)} · "
            f"fills {round(s['sphereFill'] * 100)}% of its circumsphere.* **{s['quality']}.** {s['qualityText']}",
            "",
            f"- **Integrated — {plain(s['archetype'])}:** {s['archetypeGloss']}",
            f"- **Shadow — {s['shadow']}:** {s['shadowGloss']}",
            f"- **Drive:** {s['drive']} **Fear:** {s['fear']}",
            f"- **Somatic, integrated:** {s['somaticInt']}",
            f"- **Somatic, shadow:** {s['somaticSha']}",
            f"- **Colour:** `{s['colorInt']['hex']}` {s['colorInt']['name']} → `{s['colorSha']['hex']}` "
            f"{s['colorSha']['name']}. {s['colorShift']} {s['colorEffect']}",
            "",
            "| Sense | Anchor (invoke) | Shadow tell (recognise) |",
            "|---|---|---|",
        ]
        lines += [f"| {l} | {cell(s['senses'][k]['a'])} | {cell(s['senses'][k]['t'])} |" for k, l in SENSES]
        lines += ["", "**Cultural objects:**"]
        lines += [f"- *{o['name']}* ({o['year']}). Light: {o['light']} Shadow: {o['shadow']}" for o in s["objects"]]
        lines += ["", f"**Transmutation:** {s['transmutation']}", "", f"**Sigil:** {s['sigil']}", ""]
    lines += [
        "## How the geometry mirrors the growth",
        "",
        "- **Connections only grow.** Edges climb 6, 12, 12, 30, 30: development as more relationships between "
        "the parts, not more parts.",
        "- **One seed, then two mirrors.** The tetrahedron is its own dual (the Child meets only itself). Cube and "
        "octahedron are duals (the persona's outer face and the inner other Jung called complementary to it). "
        "Icosahedron and dodecahedron are duals (the Sage knows *about* wholeness from twenty viewpoints; the Self "
        "*is* it, and the inflated guru is its counterfeit).",
        "- **Jung's own sequence.** *The Relations between the Ego and the Unconscious* (1928) runs persona, "
        "anima/animus, mana-personality, self; with the Child (1940) as the seed, the five line up one to one.",
        "- **Faces change kind twice.** Triangles for spark, mirror and knowledge; squares for structure; pentagons, "
        "carrying the golden ratio, only at the end: wholeness built on a proportion that never closes.",
        "- **Round is not whole.** The icosahedron is the most spherical (0.94 vs 0.91); the dodecahedron fills the "
        "most of its sphere (66% vs 61%). One dip, cube to octahedron (37% to 32%): turning to the inner other "
        "first feels like losing ground.",
        "",
        "*These are correspondences: useful patterns to work with, not causes. A reflective practice, not therapy.*",
        "",
    ]
    with open(OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines))
    print(f"Wrote {os.path.relpath(OUT, HERE)} ({sum(len(l) + 1 for l in lines)} bytes)")


if __name__ == "__main__":
    data = load_data()
    html = open(PAGE, encoding="utf-8").read()
    check_shape(data)
    check_geometry(data["solids"])
    check_shift_claims(check_colours(data["solids"]))
    check_exclusions(data, html)
    check_objects(data)
    write_markdown(data)
    print(f"\n{len(failures)} failure(s)" if failures else "\nAll checks passed.")
    sys.exit(1 if failures else 0)
