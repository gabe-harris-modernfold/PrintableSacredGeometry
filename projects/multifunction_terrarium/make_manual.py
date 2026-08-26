"""Build TERRARIUM_MANUAL.pdf -- what the object is, what it does, and how it is used.

Text is sourced from BRIEF.md / IDEA.md and the built geometry (params.py, verify.py).
Figures: cropped panels of terrarium_final.png plus drawn circuit / clock diagrams.
"""
import os
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.platypus import (BaseDocTemplate, PageTemplate, Frame, Paragraph,
                                Spacer, Table, TableStyle, Image as RLImage,
                                PageBreak)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "TERRARIUM_MANUAL.pdf")

import trimesh as _tm
import params as _P

#: Measured off the built part rather than typed in. This number already moved once,
#: when the tap was seated against the dome (94 -> 23 mm), and a manual that disagrees
#: with the mesh is worse than one that leaves the number out.
_STEM_PROUD = float(_tm.load(os.path.join(HERE, "drip_splitter.stl")).bounds[1][2])               - _P.VESSEL_H

INK = colors.HexColor("#1d2b33")
MUTED = colors.HexColor("#5c6b74")
ACC = colors.HexColor("#22648f")
RULE = colors.HexColor("#c8d4db")
BAND = colors.HexColor("#eef4f8")


# ------------------------------------------------------------------ figures
def crop_panels():
    """terrarium_final.png is a 2x2 of assembled / section / eye-level / crown."""
    im = Image.open(os.path.join(HERE, "terrarium_final.png"))
    w, h = im.size
    names = {"fig_assembled": (0, 0), "fig_section": (1, 0),
             "fig_eye": (0, 1), "fig_crown": (1, 1)}
    out = {}
    for name, (cx, cy) in names.items():
        box = (cx * w // 2, cy * h // 2, (cx + 1) * w // 2, (cy + 1) * h // 2)
        p = os.path.join(HERE, name + ".png")
        im.crop(box).save(p)
        out[name] = p
    return out


def circuit_diagram():
    """The water loop, drawn as it actually runs."""
    fig, ax = plt.subplots(figsize=(7.4, 5.4), dpi=220)
    ax.set_xlim(0, 10.4)
    ax.set_ylim(0.6, 11.4)
    ax.axis("off")
    fig.patch.set_facecolor("white")

    L, C, R = 2.0, 5.3, 8.6
    boxes = [
        ("sun", R, 10.4, "SUN + clear PETG shell\ngreenhouse dT 3-8 K, daily",
         "#f7f0da", "#8a7b45"),
        ("fog", R, 8.0, "crown condensate\nfog -> gutters -> drip",
         "#eaf3ea", "#4d7a54"),
        ("pump", L, 6.6, "solar pump ~0.5 W\n6 L/h via restrictor",
         "#e6eef4", "#22648f"),
        ("split", C, 10.4, "SPLITTER\n10 mm in, 4 levels", "#e6eef4", "#22648f"),
        ("spiral", C, 8.4, "16 SPIRALS 2.5 mm bore\n1680 notches",
         "#e6eef4", "#22648f"),
        ("casc", C, 6.2, "CASCADE SCREEN\n14 terraces x 10 mm\n~3000 drip sites",
         "#dcecf6", "#22648f"),
        ("bed", C, 4.0, "LIVING BED  32 deg trickle\nmoss, springtails, isopods,\n"
         "worms, snails", "#e3f0e3", "#4d7a54"),
        ("siph", C, 2.0, "BELL SIPHON\nfill -> trip -> dump -> reset",
         "#f6e8e8", "#8a4b4b"),
        ("res", L, 2.0, "RESERVOIR ~5.9 L\nwaterline z = 92 mm", "#e6eef4", "#22648f"),
    ]
    pos = {}
    for key, x, y, t, fc, ec in boxes:
        n = t.count("\n") + 1
        w, h = 2.9, 0.34 * n + 0.34
        ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                    boxstyle="round,pad=0.06,rounding_size=0.12",
                                    fc=fc, ec=ec, lw=1.1))
        ax.text(x, y, t, ha="center", va="center", fontsize=7.4, color="#1d2b33",
                linespacing=1.35)
        pos[key] = (x, y, w, h)

    def arrow(a, b, label=None, dx=0.0, dy=0.0):
        x0, y0, w0, h0 = pos[a]
        x1, y1, w1, h1 = pos[b]
        vx, vy = x1 - x0, y1 - y0
        n = max((vx ** 2 + vy ** 2) ** 0.5, 1e-9)

        def trim(w, h):
            tx = (w / 2) / abs(vx / n) if abs(vx / n) > 1e-6 else 1e9
            ty = (h / 2) / abs(vy / n) if abs(vy / n) > 1e-6 else 1e9
            return min(tx, ty) + 0.10

        s0, s1 = trim(w0, h0), trim(w1, h1)
        p0 = (x0 + vx / n * s0, y0 + vy / n * s0)
        p1 = (x1 - vx / n * s1, y1 - vy / n * s1)
        ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=11,
                                     lw=1.3, color="#3f9ad0"))
        if label:
            ax.text((p0[0] + p1[0]) / 2 + dx, (p0[1] + p1[1]) / 2 + dy, label,
                    fontsize=6.8, color="#5c6b74", ha="center", va="center")

    arrow("res", "pump", "lift", dx=-0.55)
    arrow("pump", "split")
    arrow("split", "spiral")
    arrow("spiral", "casc")
    arrow("sun", "fog")
    arrow("fog", "casc")
    arrow("casc", "bed")
    arrow("bed", "siph")
    arrow("siph", "res", "dump", dy=-0.78)
    ax.text(5.2, 0.85, "one loop, four clocks:   drop ~0.15 s   |   drip ~1 s   |   "
                       "siphon ~10 min   |   solar day 24 h",
            ha="center", fontsize=7.6, color="#22648f")
    p = os.path.join(HERE, "fig_circuit.png")
    fig.savefig(p, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


def clock_diagram():
    """The same accumulate->threshold->discharge->reset law at four timescales."""
    fig, ax = plt.subplots(figsize=(7.4, 2.4), dpi=220)
    t = np.linspace(0, 4, 2000)
    rows = [("drop at a notch", "~0.15 s", "#3f9ad0"),
            ("terrace fill / drip", "~1 s", "#22648f"),
            ("bell siphon", "~10 min", "#8a4b4b"),
            ("solar day", "24 h", "#8a7b45")]
    for i, (name, per, c) in enumerate(rows):
        u = t % 1.0
        y = np.where(u > 0.86, 1 - (u - 0.86) / 0.14, u / 0.86)
        ax.plot(t, y * 0.72 + (3 - i) * 1.0, color=c, lw=1.5)
        ax.text(-0.12, (3 - i) * 1.0 + 0.36, name, ha="right", va="center",
                fontsize=7.6, color="#1d2b33")
        ax.text(4.12, (3 - i) * 1.0 + 0.36, per, ha="left", va="center",
                fontsize=7.6, color="#5c6b74")
    ax.set_xlim(-2.3, 5.2)
    ax.set_ylim(-0.45, 4.1)
    ax.axis("off")
    ax.text(2.0, -0.32, "accumulate  ->  threshold  ->  discharge  ->  reset",
            ha="center", fontsize=7.8, color="#22648f")
    p = os.path.join(HERE, "fig_clocks.png")
    fig.savefig(p, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


# ------------------------------------------------------------------ styles
ss = getSampleStyleSheet()
S = {}
S["title"] = ParagraphStyle("title", parent=ss["Title"], fontName="Helvetica-Bold",
                            fontSize=27, leading=31, textColor=INK, alignment=0,
                            spaceAfter=2)
S["sub"] = ParagraphStyle("sub", parent=ss["Normal"], fontName="Helvetica",
                          fontSize=11.5, leading=16, textColor=ACC, spaceAfter=10)
S["h1"] = ParagraphStyle("h1", parent=ss["Heading1"], fontName="Helvetica-Bold",
                         fontSize=14.5, leading=18, textColor=INK,
                         spaceBefore=14, spaceAfter=5)
S["h2"] = ParagraphStyle("h2", parent=ss["Heading2"], fontName="Helvetica-Bold",
                         fontSize=10.8, leading=14, textColor=ACC,
                         spaceBefore=10, spaceAfter=3)
S["body"] = ParagraphStyle("body", parent=ss["BodyText"], fontName="Helvetica",
                           fontSize=9.5, leading=13.6, textColor=INK,
                           alignment=TA_JUSTIFY, spaceAfter=6)
S["cap"] = ParagraphStyle("cap", parent=ss["Normal"], fontName="Helvetica-Oblique",
                          fontSize=8, leading=11, textColor=MUTED, spaceBefore=3,
                          spaceAfter=8)
S["bullet"] = ParagraphStyle("bullet", parent=S["body"], leftIndent=11,
                             bulletIndent=2, spaceAfter=3, alignment=0)
S["cell"] = ParagraphStyle("cell", parent=ss["Normal"], fontName="Helvetica",
                           fontSize=8.3, leading=11, textColor=INK)
S["cellb"] = ParagraphStyle("cellb", parent=S["cell"], fontName="Helvetica-Bold",
                            textColor=colors.white)


def P_(t, s="body"):
    return Paragraph(t, S[s])


def B(t):
    return Paragraph(t, S["bullet"], bulletText="•")


def table(data, widths, head=True):
    rows = [[Paragraph(c, S["cellb"] if (head and i == 0) else S["cell"]) for c in r]
            for i, r in enumerate(data)]
    t = Table(rows, colWidths=widths, hAlign="LEFT", repeatRows=1 if head else 0)
    cmds = [("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("TOPPADDING", (0, 0), (-1, -1), 4),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
            ("LEFTPADDING", (0, 0), (-1, -1), 6),
            ("RIGHTPADDING", (0, 0), (-1, -1), 6),
            ("LINEBELOW", (0, 0), (-1, -2), 0.4, RULE)]
    if head:
        cmds.append(("BACKGROUND", (0, 0), (-1, 0), ACC))
    start = 1 if head else 0
    for i in range(start, len(rows)):
        if (i - start) % 2 == 1:
            cmds.append(("BACKGROUND", (0, i), (-1, i), BAND))
    t.setStyle(TableStyle(cmds))
    return t


def note_box(t):
    tb = Table([[Paragraph(t, S["cell"])]], colWidths=[168 * mm], hAlign="LEFT")
    tb.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, -1), BAND),
                            ("LINEBEFORE", (0, 0), (0, -1), 2.2, ACC),
                            ("LEFTPADDING", (0, 0), (-1, -1), 8),
                            ("RIGHTPADDING", (0, 0), (-1, -1), 8),
                            ("TOPPADDING", (0, 0), (-1, -1), 7),
                            ("BOTTOMPADDING", (0, 0), (-1, -1), 7)]))
    return tb


def img(path, w_mm, cap=None):
    im = Image.open(path)
    w = w_mm * mm
    out = [RLImage(path, width=w, height=w * im.size[1] / im.size[0])]
    if cap:
        out.append(P_(cap, "cap"))
    return out


# ------------------------------------------------------------------ page frame
def decorate(canvas, doc):
    canvas.saveState()
    W, H = A4
    canvas.setStrokeColor(RULE)
    canvas.setLineWidth(0.5)
    canvas.line(21 * mm, H - 15 * mm, W - 21 * mm, H - 15 * mm)
    canvas.setFont("Helvetica", 7.5)
    canvas.setFillColor(MUTED)
    canvas.drawString(21 * mm, H - 13 * mm, "TERRARIUM  —  operating manual")
    canvas.drawRightString(W - 21 * mm, H - 13 * mm,
                           "PrintableSacredGeometry / terrarium")
    canvas.line(21 * mm, 15 * mm, W - 21 * mm, 15 * mm)
    canvas.drawRightString(W - 21 * mm, 10.5 * mm, str(doc.page))
    canvas.drawString(21 * mm, 10.5 * mm,
                      "clear PETG · 304 mm across corners · 540 mm tall "
                      "· 13 printed parts")
    canvas.restoreState()


def build():
    figs = crop_panels()
    circ = circuit_diagram()
    clocks = clock_diagram()

    doc = BaseDocTemplate(OUT, pagesize=A4,
                          leftMargin=21 * mm, rightMargin=21 * mm,
                          topMargin=21 * mm, bottomMargin=20 * mm,
                          title="Terrarium — operating manual",
                          author="PrintableSacredGeometry")
    frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id="f")
    doc.addPageTemplates([PageTemplate(id="main", frames=[frame], onPage=decorate)])

    E = []
    A = E.append

    # ---------------------------------------------------------- cover
    A(P_("Terrarium", "title"))
    A(P_("An instrument for making ambient gradients visible — "
         "a sealed, self-watering living column in clear PETG", "sub"))
    A(note_box(
        "<b>What it does.</b> Sunlight on a clear PETG shell makes a small temperature "
        "difference between the warm headspace and the cool wet bed. A half-watt solar "
        "pump lifts water to the crown; from there gravity, surface tension and one "
        "bell siphon do everything else. The water is split into sixteen spiral "
        "gutters, shed as thousands of separate drops down a fourteen-terrace cascade, "
        "trickled through a living moss-and-invertebrate bed, collected underneath, and "
        "dumped back to the reservoir when the siphon trips. Nothing is timed, switched "
        "or engraved: every visible rhythm is the same "
        "<i>accumulate → threshold → discharge → reset</i> law running at "
        "a different scale, and every dimension in the object is set by a physical "
        "constant rather than by taste."))
    A(Spacer(1, 6))
    E.extend(img(figs["fig_assembled"], 105))
    A(P_("Assembled: a twelve-sided tapered body under an ovoid crown, 304 mm across "
         "corners and 540 mm tall, standing on a 750 mm table.", "cap"))

    A(table([["At a glance", ""],
             ["Overall", "304 mm across corners × 540 mm tall; footprint inside the "
              "320 mm print bed"],
             ["Material", "Clear PETG throughout; 0.6 mm flat-printed viewing panes"],
             ["Parts", "13 printed parts — 3 shell modules, reservoir, bed tray "
              "with siphon, cascade screen, spiral gutters, splitter, 5 pane types "
              "× 12"],
             ["Water", "~5.9 L reservoir; 6 L/h trickle; ~48 L/h equivalent during a "
              "siphon dump"],
             ["Power", "One ~0.5 W solar pump. No heater, no timer, no electronics"],
             ["Livestock", "Moss, springtails, isopods, worms, snails. No fish"],
             ["Viewing", "60 ports; facets tilted so the worst-case sightline meets "
              "glass at 14°"]],
            [30 * mm, 138 * mm]))

    A(PageBreak())

    # ---------------------------------------------------------- 1
    A(P_("1 · What it is", "h1"))
    A(P_("This is not an ornament with plants in it, and it is not an aquarium. It is a "
         "<b>demonstration instrument</b>. Its claim is that a force too small to feel "
         "— a few kelvin of greenhouse warming, a surface-tension gradient, the "
         "head inside a 2.5 mm bore — can be made plainly visible if you put it "
         "somewhere it has no choice but to become a pattern. Geometry is the only "
         "actuator here. It never pushes; it sets boundary conditions and waits."))
    A(P_("There are exactly three ways shape does that, and every feature of the object "
         "is one of them:"))
    A(table([["Amplifier", "What the geometry does", "Where you see it"],
             ["Threshold", "Poise the system at a tipping point so an imperceptible "
              "input produces a discrete, visible event.",
              "Every drip notch; the bell siphon"],
             ["Accumulation", "Collect a trickle over minutes into one thing that falls.",
              "Terrace fill; the siphon flood volume"],
             ["Selection", "A dimension picks one mode out of broadband noise — "
              "depth sets cell size, slot width sets rise height.",
              "Drop spacing; convection cell size"]],
            [24 * mm, 82 * mm, 62 * mm]))
    A(Spacer(1, 4))
    A(P_("Because the shell is clear PETG it passes sunlight but blocks thermal "
         "infrared — so it is a greenhouse, and in any sunlit room it generates a "
         "3–8 K gradient between headspace and wet bed: free, every day, dying at "
         "night and returning with the sun. That gradient is the prime mover behind "
         "everything the pump does not do."))
    A(P_("The same material sets the look. PETG is only genuinely <i>clear</i> where the "
         "wall is one extrusion wide, so “thin wall, mostly viewing ports” is "
         "enforced by physics rather than chosen as a style: the frames carry the "
         "structure, and the panes are 0.6 mm, printed flat so both faces come out "
         "optically flat."))
    A(Spacer(1, 2))
    E.extend(img(figs["fig_section"], 168,
           "Section: reservoir at the bottom, living bed and bell siphon above it, the "
           "cascade screen through the mid volume, and the drip spirals riding the "
           "inside of the shell up into the crown."))

    A(PageBreak())

    # ---------------------------------------------------------- 2
    A(P_("2 · How the water goes round", "h1"))
    E.extend(img(circ, 146, "The single loop. The pump is the only powered step; everything "
                      "after it is gravity, surface tension and one siphon."))
    A(P_("<b>Lift.</b> A ~0.5 W solar pump sits in the reservoir bay and pushes about "
         "6 L/h up a 10 mm intake to the crown. That is the entire energy budget — "
         "enough for a pump, nowhere near enough for a heater, which is why the "
         "livestock list contains nothing that needs one."))
    A(P_("<b>Split.</b> A four-level splitter divides the 10 mm intake into sixteen "
         "2.5 mm bores, each level conserving cross-sectional area so no branch starves. "
         "Sixteen is not decorative: below about nine branches, per-branch flow during a "
         "siphon dump exceeds the jet limit and the drip sites stop dripping and start "
         "squirting."))
    A(P_("<b>Tap.</b> The splitter is fed twice. The pump supplies the bores; the dome "
         "supplies the outside. A cone hangs under the apex with a flat seat that "
         "<i>touches</i> the crown's inner face on a 24 mm circle, and that contact is "
         "the whole mechanism: film creeping down the ceiling is drawn into a wetted "
         "contact line instead of hanging up and nucleating pendant drops every 17 mm, "
         "which is what a ceiling with nothing touching it does. Inside the seat the "
         "gap to the apex is 1.3 mm — under the 2.7 mm capillary length, so it fills "
         "and feeds the same contact rather than dripping. Condensate runs down the "
         "cone flank, sheds at its rim and continues on the branch exteriors to the "
         "same sixteen mouths the pump feeds from within."))
    A(P_("<b>Spiral.</b> The sixteen gutters descend the inside of the shell as open "
         "C-sections, not sealed pipes. At 2.5 mm a closed pipe would need about 576 Pa "
         "to break the meniscus at a hole in its underside while the bore supplies only "
         "about 25 Pa of head; it would hold every drop and discharge at the far end "
         "alone. An open gutter has a free surface, sheds cleanly from a notched lip, "
         "and catches the condensate running down the crown as well. Each spiral is "
         "stepped rather than smoothly helical, so every level run fills and spills "
         "along its whole notched lip instead of racing to the first notch."))
    A(P_("<b>Cascade.</b> Water reaching the mid volume falls onto a staircase cone: "
         "fourteen terraces at a 10 mm riser, each lip sitting directly over the next "
         "tread so a drop lands and <i>re-forms</i> rather than falling clear. One fall "
         "becomes twenty-one detachments. The lips are scalloped into 36 lobes purely to "
         "buy arc length — more drip sites without moving the cone out toward the "
         "viewing ports."))
    A(P_("<b>Bed.</b> The living bed is a tray on a 32° slope over an air gap, so "
         "it <i>trickles</i> and never ponds. This is the one non-negotiable rule of the "
         "whole object: a litter bed that ponds goes anaerobic and poisons everything "
         "below it."))
    A(P_("<b>Siphon and reset.</b> Under the bed sits a bell siphon — standpipe, "
         "bell, snorkel. Water accumulates until it reaches the standpipe crown, the "
         "bell primes, and the whole flood volume dumps to the reservoir in seconds at "
         "roughly 48 L/h equivalent; then air breaks the seal through the snorkel and it "
         "resets. No timer, no float switch, no electronics — a genuine relaxation "
         "oscillator, which is the honest way to demonstrate vibration."))

    A(PageBreak())

    # ---------------------------------------------------------- 3
    A(P_("3 · What you are meant to watch", "h1"))
    A(P_("The object runs four clocks at once, all obeying the same law. That "
         "correspondence is the point, and it is the part you can put a stopwatch on "
         "rather than assert."))
    E.extend(img(clocks, 168, "Four spans, one law: fill slowly, cross a threshold, dump "
                        "fast, reset."))
    A(table([["Timescale", "What happens", "How to see it"],
             ["~0.15 s", "A drop grows at a 0.4 mm notch until surface tension loses to "
              "gravity, then falls.", "Watch one notch against a dark card"],
             ["~1 s", "A terrace fills and spills along its whole lip.",
              "Count drops crossing one theatre gap"],
             ["~10 min", "The bell siphon floods, trips, drains and resets.",
              "The loud clock — audible across the room"],
             ["24 h", "The greenhouse gradient rises with the sun and dies at night; the "
              "pump follows it.", "Log headspace vs bed temperature"]],
            [22 * mm, 86 * mm, 60 * mm]))
    A(Spacer(1, 5))
    A(P_("Numbers that make the drops countable", "h2"))
    A(table([["Quantity", "Value", "Why it is that number"],
             ["Drop", "5.54 µL, 2.19 mm",
              "Tate's law off a 0.4 mm lip — the smallest drop the printer can make"],
             ["Drip sites", "~3,000 cascade + 1,680 spiral",
              "Notch pitch 6 mm and 12 mm: wider than two drop diameters, so drops never "
              "coalesce"],
             ["Terrace riser", "10 mm",
              "Floor is 8.8 mm (4 × drop diameter), below which a drop bridges and "
              "never detaches"],
             ["Terrace tread", "2.86 mm", "At or above the 2.71 mm capillary length"],
             ["Jet limit", "19.3 drops/s per site",
              "Above this a site stops dripping and becomes a jet"],
             ["Airborne at once", "tens at trickle, hundreds during a dump",
              "Fall time across a 10 mm step"]],
            [30 * mm, 42 * mm, 96 * mm]))
    A(Spacer(1, 5))
    A(note_box("<b>Two registers, one object.</b> The siphon and the cascade are the "
               "loud clock — bulk water, audible, unmissable. Riding on top of it "
               "and running on nothing is the quiet layer: hexagonal convection cells in "
               "a shallow dish, dew spacing itself along a wetted fibre, films meeting "
               "at exactly 120°. Same law, five orders of magnitude apart."))

    A(PageBreak())

    # ---------------------------------------------------------- 4
    A(P_("4 · Setting it up", "h1"))
    A(P_("Printing", "h2"))
    A(table([["Part", "Qty", "Notes"],
             ["Shell modules 0 / 1 / 2", "1 each",
              "Up to 306 × 306 × 181 mm; they stack"],
             ["Viewing panes", "5 types × 12 = 60",
              "0.6 mm, printed <b>flat</b> — this is what buys the clarity"],
             ["Reservoir pan", "1", "1.6 mm minimum wall below any waterline"],
             ["Bed tray + siphon", "1",
              "One body; the 32° slope is structural, do not flatten it"],
             ["Cascade screen", "1", "218 × 192 × 140 mm, 0.9 mm walls"],
             ["Drip gutters", "1",
              "All 16 spirals as one part, 275 × 275 × 289 mm"],
             ["Drip splitter", "1",
              "10 mm in, 16 × 2.5 mm out, condenser cone on top"]],
            [42 * mm, 24 * mm, 102 * mm]))
    A(Spacer(1, 4))
    A(B("Clear PETG, 0.4 mm nozzle, 0.2 mm layer. Wet walls 1.6 mm; lattice ribs and "
        "screen walls 0.9 mm."))
    A(B("Overhangs are held within 45° of vertical so nothing needs support inside "
        "the water path — support scars in a gutter pin droplets and kill drip "
        "sites."))
    A(B("Do not sand or polish the panes. Print them flat and leave them: the "
        "as-printed faces are flatter than anything achieved by hand."))
    A(Spacer(1, 3))
    A(P_("Not printed — you supply", "h2"))
    A(B("Solar panel, ~0.5 W submersible pump, silicone tubing and a flow restrictor set "
        "to 6 L/h."))
    A(B("Pane gaskets and retention clips."))
    A(B("Substrate, moss, leaf litter and livestock."))
    A(Spacer(1, 4))
    A(P_("Assembly order", "h2"))
    steps = [
        "Stand the reservoir pan on a level surface. Fit the pump in its bay and route "
        "the 10 mm feed line up through the centre.",
        "Drop shell module 0 over the reservoir; glaze its ports with the flat panes and "
        "clips before the next module buries them.",
        "Set the bed tray with its bell siphon on the module 0 rim. Check by eye that "
        "the bell cap sits clear <i>below</i> the tray rim — if it does not, the "
        "tray overflows before the siphon ever trips.",
        "Seat the cascade screen over the bed, then module 1, glazing as you go.",
        "Lower in the spiral gutter assembly, then close with module 2 and the crown "
        "panes. The splitter’s stem passes up through the boss at the apex and "
        f"stands {_STEM_PROUD:.0f} mm proud of it; the feed line connects there. "
        "Press the assembly up "
        "until the cone’s seat meets the inside of the dome — it is meant to "
        "touch. A gap there costs you the condensate, which drips off the ceiling at "
        "random instead of being collected.",
        "Fill the reservoir to the waterline at z = 92 mm (~5.9 L). Use rainwater, "
        "distilled or RO water — never softened or heavily mineralised tap water.",
        "Run the pump from a bench supply for an hour before trusting the solar panel, "
        "and walk the whole path looking for water going anywhere it should not.",
    ]
    for i, t in enumerate(steps):
        A(B(f"<b>{i + 1}.</b> {t}"))

    A(PageBreak())

    # ---------------------------------------------------------- 5
    A(P_("5 · Living in it", "h1"))
    A(P_("Planting and stocking", "h2"))
    A(P_("The bed is a trickle bed, not a pot. Build it shallow: a thin drainage layer, "
         "then leaf litter and moss laid directly on the slope. Anything that dams the "
         "32° fall of the tray — a deep peat layer, a flat stone laid across "
         "the slope — defeats the drainage and creates exactly the anaerobic pocket "
         "the geometry exists to prevent."))
    A(table([["Species", "Role", "Notes"],
             ["Moss (sheet or cushion)",
              "Ground cover, humidity buffer, the visible health indicator",
              "Wants light and trickle, never standing water"],
             ["Springtails", "Mould control, first-line cleanup",
              "Seed heavily; they are why mould never gets established"],
             ["Isopods (dwarf species)", "Litter breakdown",
              "Keep to small species; large ones tunnel and slump the bed"],
             ["Composting worms", "Litter to castings, aeration",
              "A few only — the bed is thin, this is not a wormery"],
             ["Snails (small)", "Algae grazing on panes and gutters",
              "Least forgiving of neglect; watch them first"],
             ["No fish", "—",
              "A 320 mm footprint cannot hold enough water for one without wrecking the "
              "form"]],
            [40 * mm, 52 * mm, 76 * mm]))
    A(Spacer(1, 4))
    A(P_("Where to put it", "h2"))
    A(B("Bright indirect light, or a few hours of direct sun. The greenhouse gradient is "
        "the prime mover — in a dark corner the object still runs on the pump, but "
        "the quiet layer stops."))
    A(B("On a table at about 750 mm. The facets are tilted to split the difference "
        "between a standing and a seated eye 700 mm back; at that geometry the worst "
        "sightline meets glass at 14°, well inside the 45° where reflection "
        "takes over and a pane turns into a mirror."))
    A(B("Away from draughts and heat sources. Both flatten the very gradient the object "
        "is built to display."))
    A(Spacer(1, 4))
    A(P_("Routine", "h2"))
    A(table([["When", "Do"],
             ["Daily (10 s)", "Look. Is water dripping at the crown? Is the siphon still "
              "cycling? Is the bed damp rather than glossy-wet?"],
             ["Weekly", "Top up evaporation loss with rain, distilled or RO water. Check "
              "the top spiral's notches for blockage."],
             ["Monthly", "Lift the crown module; prune moss back off the panes and clear "
              "litter from the gutter lips."],
             ["Twice a year", "Strip to modules, flush the reservoir and gutters, re-seat "
              "the panes. No detergent anywhere in the water path."]],
            [30 * mm, 138 * mm]))
    A(Spacer(1, 4))
    A(note_box("<b>Never put detergent, soap or fertiliser in the water.</b> Beyond the "
               "obvious harm to the livestock, organic surfactants flatten "
               "surface-tension-driven flow dead — the quiet half of the "
               "demonstration stops and does not come back until every trace is flushed "
               "out."))

    A(PageBreak())

    # ---------------------------------------------------------- 6
    A(P_("6 · When it misbehaves", "h1"))
    A(table([["Symptom", "Most likely cause", "Fix"],
             ["Siphon never trips", "Flow below the priming rate, or a leak at the "
              "standpipe joint",
              "Open the restrictor toward 6 L/h; check the bell seats flat on its feet"],
             ["Siphon never <i>stops</i>",
              "Snorkel blocked, so air cannot break the seal",
              "Clear the snorkel bore; it must stay open to the headspace"],
             ["Bed tray overflows", "Water arriving faster than the siphon drains, or "
              "the bell cap sitting at or above the tray rim",
              "Reduce flow; verify the cap sits clear below the rim"],
             ["Drips stop, water sheets instead",
              "Per-site rate above the jet limit — usually a blocked branch pushing "
              "flow into fewer notches", "Flush the splitter; restore all 16 branches"],
             ["One spiral runs dry", "Air trapped in a 2.5 mm bore",
              "Briefly raise the flow to purge, then return to trickle"],
             ["Panes fog and stay blind",
              "Condensate pinning on layer lines instead of running off",
              "Confirm the crown gutters are clear; a flame or plasma pass on the inner "
              "face restores run-off"],
             ["Bed smells sour",
              "It is ponding — the one failure that kills everything",
              "Clear the drain path at once; thin the substrate; restore the air gap "
              "under the tray"],
             ["Algae on the panes", "Too much direct sun plus nutrients",
              "Move to indirect light; add snails; never dose fertiliser"],
             ["Pump runs only in bright sun", "Correct behaviour",
              "Solar-only was chosen; the daily cycle is part of the demonstration"]],
            [36 * mm, 62 * mm, 70 * mm]))

    A(Spacer(1, 6))
    A(P_("7 · What this object does not do", "h1"))
    A(P_("Stated rather than quietly omitted:"))
    A(B("<b>No heating.</b> A 0.5 W solar budget runs a pump; a heater needs ~25 W. "
        "Stock accordingly."))
    A(B("<b>No fish, and no closed-loop nitrogen claim.</b> This is a moss-and-"
        "invertebrate system, not an aquaponic one."))
    A(B("<b>Only one bell siphon is built.</b> The design intent is three, with flood "
        "volumes at 1:2:3 so their periods beat on a common 6T — three "
        "oscillators, one law, a harmonic ratio. This vessel has a single bed level, so "
        "it runs one."))
    A(B("<b>The electrokinetic option was analysed and rejected</b>, not forgotten. "
        "Streaming potential in terrarium water is millivolts, and the Kelvin-dropper "
        "induction that does reach kilovolts would turn 2,000 drip notches into 2,000 "
        "corona points — ozone inside a sealed living volume. It stays out."))
    A(B("<b>No engraved symbolism.</b> Where the object is hermetic it is hermetic "
        "through working parts: hexagons that arrive because a dish is 2 mm deep, "
        "120° because that is the only angle three films will meet at, drop spacing "
        "because a fibre has a radius."))

    A(Spacer(1, 8))
    A(P_("Every dimension, and the constant behind it", "h2"))
    A(table([["Dimension", "Value", "Constant it comes from"],
             ["Cell size floor", "3.0 mm",
              "Capillary length 2.709 mm — below it one drop spans the cell"],
             ["Rib height", "3.5 mm", "Maximum sessile puddle height, 3.42 mm"],
             ["Drip spacing", "≥ 6 mm",
              "Rayleigh–Taylor spacing 17.0 mm: closer sites are ours, wider and "
              "the physics chooses them"],
             ["Terrace riser", "10 mm",
              "4 × drop diameter = the 8.8 mm bridging floor"],
             ["Headspace", "~50 mm",
              "Rayleigh number 12,400 against a critical 1,708 — it convects on its "
              "own"],
             ["Marangoni dish", "2 mm deep",
              "Supercritical at ΔT ≈ 0.04 K; cells come out 4–6 mm"],
             ["Capillary slot", "0.4 mm",
              "Jurin's law at PETG's 70° contact angle → ~13 mm lift"],
             ["Facet tilt", "per band",
              "Fresnel: clear below 45° incidence, a mirror by 70°"],
             ["Footprint", "306 mm", "The 320 mm print bed"]],
            [34 * mm, 26 * mm, 108 * mm]))
    A(Spacer(1, 5))
    A(P_("Sources: BRIEF.md (requirements), IDEA.md (physics), params.py (every number "
         "above), verify.py (the built geometry re-measured from the exported STLs).",
         "cap"))

    A(Spacer(1, 10))
    A(P_("Two more views", "h2"))
    def _im(path, w_mm):
        im = Image.open(path)
        w = w_mm * mm
        return RLImage(path, width=w, height=w * im.size[1] / im.size[0])
    pair = Table([[_im(figs["fig_eye"], 82), _im(figs["fig_crown"], 82)]],
                 colWidths=[84 * mm, 84 * mm], hAlign="LEFT")
    pair.setStyle(TableStyle([("LEFTPADDING", (0, 0), (-1, -1), 0),
                              ("RIGHTPADDING", (0, 0), (-1, -1), 2),
                              ("TOPPADDING", (0, 0), (-1, -1), 0),
                              ("BOTTOMPADDING", (0, 0), (-1, -1), 0)]))
    A(pair)
    A(P_("Left: what a standing eye sees from 700 mm back, with the object on a 750 mm "
         "table. Right: the crown — the condenser tap, sixteen spirals, and the "
         "ovoid head it taps the condensate from.", "cap"))

    doc.build(E)
    print("wrote", OUT)


if __name__ == "__main__":
    build()
