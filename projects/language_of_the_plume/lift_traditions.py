"""Lift the smoke-reading traditions from the research writeup into the page.

docs/smoke-reading-traditions.md ends in a JSON block (lenses, forms, timeline) and carries a Sources
table. This copies both, verbatim apart from markdown emphasis, into web/language-of-the-plume.html
between the TRADITIONS markers, so the page shows exactly what the writeup sources. Run it after the
writeup changes, then core.py and check_plume.js:

    cd projects/language_of_the_plume && python lift_traditions.py && python core.py && node check_plume.js
"""
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
DOC = os.path.join(HERE, "..", "..", "docs", "smoke-reading-traditions.md")
PAGE = os.path.join(HERE, "..", "..", "web", "language-of-the-plume.html")
BEGIN, END = "/* TRADITIONS:BEGIN */", "/* TRADITIONS:END */"

doc = open(DOC, encoding="utf-8").read()
data = json.loads(re.findall(r"```json\s*(.*?)```", doc, re.S)[-1])
for k in ("lenses", "forms", "timeline"):
    data.setdefault(k, [])

sources = {}
for line in doc.splitlines():
    cells = [c.strip() for c in line.strip().strip("|").split("|")]
    if len(cells) >= 6 and re.fullmatch(r"S\d+", cells[0]):
        url = cells[5] if cells[5].startswith("http") else None
        sources[cells[0]] = {"author": cells[1], "title": re.sub(r"[*_]", "", cells[2]),
                             "year": cells[3], "kind": cells[4], "url": url}
data["sources"] = sources

# every src must resolve, or the page would cite nothing
refs = [x.get("src") for k in ("lenses", "forms", "timeline") for x in data[k]]
refs += [e.get("src") for l in data["lenses"] for e in l["entries"]]
missing = sorted({r for r in refs if r and r not in sources})
if missing:
    raise SystemExit(f"sources missing from the table: {missing}")
data["timeline"].sort(key=lambda t: t["year"])

blob = json.dumps(data, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")
html = open(PAGE, encoding="utf-8").read()
a, b = html.index(BEGIN) + len(BEGIN), html.index(END)
html = html[:a] + "\nconst TRADITIONS=" + blob + ";\n" + html[b:]
open(PAGE, "w", encoding="utf-8", newline="\n").write(html)
print(f"lifted {len(data['lenses'])} lenses, {len(data['forms'])} forms, {len(data['timeline'])} timeline entries, "
      f"{len(sources)} sources")
