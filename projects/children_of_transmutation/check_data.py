"""Integrity check for the data blocks in web/the-children-of-transmutation.html.

Run from this directory:  python check_data.py
Parses every <script type="application/json" id="d-..."> block and checks that ids are unique,
every cross-reference resolves (sources, chronicle events, octagon pairs, ladders, chains, cards,
element arrows), registers and lanes are known, and chain efficiencies are sane.
Exits non-zero on any failure.
"""
import json
import os
import re
import sys

PAGE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "web", "the-children-of-transmutation.html")
REGS = {"dem", "mea", "pra", "rep", "cla", "ref"}
TREGS = {"inherit", "bridge", "rhyme", "contrast"}
LANES = {"energy", "matter", "art", "inner"}


def blocks(html):
    out = {}
    for m in re.finditer(r'<script type="application/json" id="(d-[a-z]+)">(.*?)</script>', html, re.S):
        out[m.group(1)] = json.loads(m.group(2))
    return out


def main():
    html = open(PAGE, encoding="utf-8").read()
    d = blocks(html)
    errs = []
    need = ["d-sources", "d-body", "d-events", "d-forms", "d-conv", "d-chains", "d-elements", "d-ladders", "d-touch"]
    for k in need:
        if k not in d:
            errs.append(f"missing block {k}")
    if errs:
        print("\n".join(errs)); sys.exit(1)

    src, ev, forms, conv = d["d-sources"], d["d-events"], d["d-forms"], d["d-conv"]
    chains, els, lads, touch, body = d["d-chains"], d["d-elements"], d["d-ladders"], d["d-touch"], d["d-body"]
    form_ids = {f["id"] for f in forms}
    ev_ids = [e["id"] for e in ev]
    pairs = {c["f"] + ">" + c["to"] for c in conv}
    arrow_ids = {a["id"] for a in els["arrows"]}
    lad_ids = {l["id"] for l in lads}
    touch_ids = {t["id"] for t in touch}
    chain_ids = {c["id"] for c in chains}
    used_src = set()

    def chk_src(ids, where):
        for s in ids or []:
            used_src.add(s)
            if s not in src:
                errs.append(f"{where}: unknown source '{s}'")

    def chk_go(go, where):
        if not go:
            return
        (k, v), = go.items()
        ok = {"oct": lambda v: v in pairs or v in form_ids, "el": lambda v: v in arrow_ids or v in {f"z{p['z']}" for p in els["synth"]},
              "lad": lambda v: v in lad_ids, "joule": lambda v: v in chain_ids, "touch": lambda v: v in touch_ids,
              "ev": lambda v: v in ev_ids}.get(k)
        if ok is None or not ok(v):
            errs.append(f"{where}: bad go {go}")

    dup = {i for i in ev_ids if ev_ids.count(i) > 1}
    if dup:
        errs.append(f"duplicate event ids: {sorted(dup)}")
    for e in ev:
        w = f"event {e['id']}"
        for k in ("y", "d", "lane", "reg", "t", "x"):
            if k not in e:
                errs.append(f"{w}: missing {k}")
        if e.get("lane") not in LANES:
            errs.append(f"{w}: bad lane")
        if e.get("reg") not in REGS:
            errs.append(f"{w}: bad reg")
        if not -1.2e6 <= e.get("y", 0) <= 2026:
            errs.append(f"{w}: year out of range")
        chk_src(e.get("s"), w); chk_go(e.get("go"), w)

    for c in conv:
        w = f"conv {c['f']}>{c['to']} '{c['t']}'"
        if c["f"] not in form_ids or c["to"] not in form_ids:
            errs.append(f"{w}: unknown form")
        if c.get("via") and c["via"] not in form_ids:
            errs.append(f"{w}: unknown via")
        if c.get("ev") and c["ev"] not in ev_ids:
            errs.append(f"{w}: unknown event {c['ev']}")
        chk_src(c.get("s"), w)

    for ch in chains:
        w = f"chain {ch['id']}"
        starts = [s for s in ch["stages"] if s.get("start")]
        if len(starts) != 1:
            errs.append(f"{w}: needs exactly one start stage")
        seen_start = False
        for s in ch["stages"]:
            if s["form"] not in form_ids:
                errs.append(f"{w}: unknown form {s['form']}")
            if seen_start and not s.get("lamp") and not (0 < s.get("eff", -1) <= 1):
                errs.append(f"{w}: stage '{s['label']}' needs 0 < eff <= 1")
            seen_start = seen_start or s.get("start", False)
        if not ch["stages"][-1].get("lamp"):
            errs.append(f"{w}: last stage must be the lamp")
        chk_src(ch.get("s"), w)

    zs = [p["z"] for p in els["synth"]]
    if len(zs) != len(set(zs)):
        errs.append("duplicate synthesised element")
    for p in els["synth"]:  # here "s" is the element symbol; sources live in "src"
        chk_src(p.get("src"), f"element {p['z']}")
    for a in els["arrows"]:
        w = f"arrow {a['id']}"
        if a.get("ev") and a["ev"] not in ev_ids:
            errs.append(f"{w}: unknown event")
        if not all(1 <= z <= 118 for z in [a["f"], *a["to"]]):
            errs.append(f"{w}: Z out of range")
        chk_src(a.get("s"), w)

    for l in lads:
        w = f"ladder {l['id']}"
        if len(l["rungs"]) < 3:
            errs.append(f"{w}: fewer than 3 rungs")
        if l["reg"] not in REGS:
            errs.append(f"{w}: bad reg")
        chk_src(l.get("s"), w)

    for t in touch:
        w = f"card {t['id']}"
        if t["reg"] not in TREGS:
            errs.append(f"{w}: bad reg")
        if not set(t["lanes"]) <= LANES:
            errs.append(f"{w}: bad lane")
        chk_src(t.get("s"), w); chk_go(t.get("go"), w)

    tot = sum(o["p"] for o in body)
    if abs(tot - 100) > 0.05:
        errs.append(f"body origins sum to {tot}, not 100")

    # relative page links must exist in web/ (index.html lives in site/ only); #src- anchors must resolve
    for m in re.finditer(r'href=\\?"([^"\\#:]+\.html)', html):
        name = m.group(1)
        if name != "index.html" and not os.path.exists(os.path.join(os.path.dirname(PAGE), name)):
            errs.append(f"link to missing page {name}")
    for m in re.finditer(r'href=\\?"#src-([a-z0-9_]+)', html):
        used_src.add(m.group(1))
        if m.group(1) not in src:
            errs.append(f"anchor to unknown source {m.group(1)}")
    for m in re.finditer(r'href=\\?"#t-([a-z0-9_]+)', html):
        if m.group(1) not in touch_ids:
            errs.append(f"anchor to unknown card {m.group(1)}")

    unused = sorted(set(src) - used_src)
    n_dem = sum(1 for e in ev if e["reg"] == "dem")
    print(f"{len(ev)} events ({n_dem} demonstrated), {len(conv)} conversions over {len(pairs)} pairs, "
          f"{len(chains)} chains, {len(els['synth'])} made elements, {len(els['arrows'])} arrows, "
          f"{len(lads)} ladders, {len(touch)} cards, {len(src)} sources")
    for lane in sorted(LANES):
        print(f"  {lane:7s} {sum(1 for e in ev if e['lane'] == lane)} events")
    if unused:
        print("note: sources cited nowhere in data:", ", ".join(unused))
    if errs:
        print("\nFAIL:\n  " + "\n  ".join(errs)); sys.exit(1)
    print("OK")


if __name__ == "__main__":
    main()
