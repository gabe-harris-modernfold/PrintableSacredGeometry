# -*- coding: utf-8 -*-
"""Build 'The Spiral Instruments' treatise PDF: 50% fluid mechanics,
50% hermetic axioms, documenting the three vortex devices."""
import os, math
from PIL import Image as PILImage
from reportlab.lib.pagesizes import LETTER
from reportlab.lib.units import inch
from reportlab.lib.colors import HexColor, Color
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer,
                                PageBreak, Image, Table, TableStyle,
                                HRFlowable, KeepTogether)
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

pdfmetrics.registerFont(TTFont('ArialU', r'C:\Windows\Fonts\arial.ttf'))
pdfmetrics.registerFont(TTFont('SegoeSym', r'C:\Windows\Fonts\seguisym.ttf'))

HERE = os.path.dirname(os.path.abspath(__file__))
OUT  = os.path.join(HERE, "The_Spiral_Instruments.pdf")

INK    = HexColor("#2b2b33")
INDIGO = HexColor("#33334f")
GOLD   = HexColor("#8f6f24")
GOLD_L = HexColor("#c9b070")
HERM   = HexColor("#4a4160")
GRAY   = HexColor("#6f6f78")
PARCH  = HexColor("#f4efdf")

# Greek and operators via embedded Unicode TTFs (base-14 fonts lack them)
OMEGA = '<font face="ArialU">ω</font>'
RHO   = '<font face="ArialU">ρ</font>'
DELTA = '<font face="ArialU">Δ</font>'
NABLA = '<font face="SegoeSym">∇</font>'

S = {}
S['series']   = ParagraphStyle('series', fontName='Helvetica', fontSize=9.5,
                textColor=GOLD, alignment=TA_CENTER, leading=14)
S['title']    = ParagraphStyle('title', fontName='Times-Bold', fontSize=30,
                textColor=INDIGO, alignment=TA_CENTER, leading=36)
S['subtitle'] = ParagraphStyle('subtitle', fontName='Times-Italic', fontSize=12.5,
                textColor=INK, alignment=TA_CENTER, leading=17)
S['epigraph'] = ParagraphStyle('epigraph', fontName='Times-Italic', fontSize=12,
                textColor=HERM, alignment=TA_CENTER, leading=17,
                leftIndent=54, rightIndent=54)
S['attrib']   = ParagraphStyle('attrib', fontName='Times-Roman', fontSize=9,
                textColor=GRAY, alignment=TA_CENTER, leading=12)
S['h1']       = ParagraphStyle('h1', fontName='Times-Bold', fontSize=16,
                textColor=INDIGO, alignment=TA_CENTER, leading=20,
                spaceBefore=18, spaceAfter=4)
S['h1sub']    = ParagraphStyle('h1sub', fontName='Times-Italic', fontSize=11,
                textColor=GOLD, alignment=TA_CENTER, leading=14, spaceAfter=10)
S['voice']    = ParagraphStyle('voice', fontName='Helvetica-Bold', fontSize=8.5,
                textColor=GOLD, alignment=TA_LEFT, leading=12,
                spaceBefore=7, spaceAfter=3)
S['sci']      = ParagraphStyle('sci', fontName='Helvetica', fontSize=9.6,
                textColor=INK, alignment=TA_JUSTIFY, leading=14, spaceAfter=7)
S['herm']     = ParagraphStyle('herm', fontName='Times-Italic', fontSize=10.8,
                textColor=HERM, alignment=TA_JUSTIFY, leading=15.5, spaceAfter=7)
S['axiom']    = ParagraphStyle('axiom', fontName='Times-Italic', fontSize=11.5,
                textColor=INDIGO, alignment=TA_CENTER, leading=16,
                leftIndent=40, rightIndent=40, spaceBefore=6, spaceAfter=2)
S['caption']  = ParagraphStyle('caption', fontName='Helvetica-Oblique', fontSize=8,
                textColor=GRAY, alignment=TA_CENTER, leading=11,
                spaceBefore=3, spaceAfter=7)
S['cell']     = ParagraphStyle('cell', fontName='Helvetica', fontSize=8.4,
                textColor=INK, alignment=TA_LEFT, leading=11.4)
S['cellb']    = ParagraphStyle('cellb', fontName='Helvetica-Bold', fontSize=8.4,
                textColor=INDIGO, alignment=TA_LEFT, leading=11.4)
S['cellq']    = ParagraphStyle('cellq', fontName='Times-Italic', fontSize=8.8,
                textColor=HERM, alignment=TA_LEFT, leading=11.6)
S['plain']    = ParagraphStyle('plain', fontName='Helvetica', fontSize=9.6,
                textColor=INK, alignment=TA_JUSTIFY, leading=14, spaceAfter=6)
S['colo']     = ParagraphStyle('colo', fontName='Helvetica', fontSize=8.6,
                textColor=GRAY, alignment=TA_LEFT, leading=12.6, spaceAfter=5)

def rule(width=2.2*inch, color=GOLD_L, th=0.7, before=4, after=6):
    return HRFlowable(width=width, thickness=th, color=color,
                      spaceBefore=before, spaceAfter=after, hAlign='CENTER')

def img(name, width_in):
    path = os.path.join(HERE, name)
    w, h = PILImage.open(path).size
    return Image(path, width=width_in*inch, height=width_in*inch*h/w,
                 hAlign='CENTER')

def axiom(text, source):
    return KeepTogether([rule(1.4*inch), Paragraph(text, S['axiom']),
                         Paragraph(source, S['attrib']), rule(1.4*inch, after=8)])

def voice(label):
    return Paragraph(label, S['voice'])

# ---------------- title page canvas: golden log-spiral ----------------
def title_page(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(Color(0.63, 0.49, 0.16, alpha=0.28))
    canvas.setLineWidth(1.0)
    cx, cy = LETTER[0]/2.0, 1.85*inch
    pts, t = [], 0.0
    while t <= 5.7*math.pi:
        r = 3.4*math.exp(0.150*t)
        pts.append((cx + r*math.cos(t), cy + r*math.sin(t)))
        t += 0.05
    p = canvas.beginPath(); p.moveTo(*pts[0])
    for x, y in pts[1:]:
        p.lineTo(x, y)
    canvas.drawPath(p)
    canvas.setFillColor(Color(0.63, 0.49, 0.16, alpha=0.5))
    canvas.circle(cx, cy, 1.6, stroke=0, fill=1)
    canvas.restoreState()

def later_pages(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(GOLD_L); canvas.setLineWidth(0.5)
    canvas.line(0.9*inch, 0.62*inch, LETTER[0]-0.9*inch, 0.62*inch)
    canvas.setFont('Helvetica', 7.2); canvas.setFillColor(GRAY)
    canvas.drawString(0.9*inch, 0.48*inch,
        "THE SPIRAL INSTRUMENTS  ·  PRINTABLE SACRED GEOMETRY")
    canvas.drawRightString(LETTER[0]-0.9*inch, 0.48*inch, "%d" % doc.page)
    canvas.restoreState()

doc = SimpleDocTemplate(OUT, pagesize=LETTER,
        leftMargin=0.95*inch, rightMargin=0.95*inch,
        topMargin=0.85*inch, bottomMargin=0.95*inch,
        title="The Spiral Instruments",
        author="G. Harris · Printable Sacred Geometry",
        subject="Three printable vortex devices, in fluid mechanics and hermetic axioms")

E = []  # the story

# ============================ TITLE ============================
E.append(Spacer(1, 40))
E.append(Paragraph("P R I N T A B L E   S A C R E D   G E O M E T R Y", S['series']))
E.append(Spacer(1, 6)); E.append(rule(3.2*inch))
E.append(Spacer(1, 14))
E.append(Paragraph("THE SPIRAL<br/>INSTRUMENTS", S['title']))
E.append(Spacer(1, 12))
E.append(Paragraph("A Treatise on Three Printable Vortex Devices, rendered in the "
    "twin languages of Fluid Mechanics and Hermetic Philosophy", S['subtitle']))
E.append(Spacer(1, 26))
E.append(Paragraph("“The Sun is its father, the Moon its mother, the Wind hath "
    "carried it in its belly, the Earth is its nurse.”", S['epigraph']))
E.append(Spacer(1, 4))
E.append(Paragraph("— The Emerald Tablet of Hermes Trismegistus, "
    "in the translation of Isaac Newton", S['attrib']))
E.append(Spacer(1, 24))
tt = Table([
    [Paragraph("<i>The Sun is its father</i>", S['cellq']),
     Paragraph("INSTRUMENT I · The Hearth Tornado", S['cellb'])],
    [Paragraph("<i>The Wind hath carried it in its belly</i>", S['cellq']),
     Paragraph("INSTRUMENT II · The Wind Tornado", S['cellb'])],
    [Paragraph("<i>The Earth is its nurse</i>", S['cellq']),
     Paragraph("INSTRUMENT III · The Dust-Devil Plate", S['cellb'])],
], colWidths=[2.5*inch, 2.6*inch])
tt.setStyle(TableStyle([
    ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
    ('TOPPADDING', (0,0), (-1,-1), 4), ('BOTTOMPADDING', (0,0), (-1,-1), 4),
    ('LINEBELOW', (0,0), (-1,1), 0.4, GOLD_L),
]))
tt.hAlign = 'CENTER'
E.append(tt)
E.append(Spacer(1, 150))
E.append(Paragraph("Set down in the month of July, 2026 · the incense_vortex folio",
    S['attrib']))
E.append(PageBreak())

# ======================= PREAMBLE =======================
E.append(Paragraph("THE DOCTRINE OF TWO LENSES", S['h1']))
E.append(rule())
E.append(Paragraph(
    "This treatise documents three desktop instruments, printed in PETG and "
    "raised around a common standard — the clear tube of a hundred millimetres — whose "
    "single purpose is to persuade air to reveal itself. The first stands a "
    "helix of incense smoke by heat alone; the second does the same work with "
    "harvested wind; the third teaches a spoonful of powdered earth to rise up "
    "and dance. All three are children of one geometry, the tangent circle, "
    "and of one law: that converging flow must spin faster, because angular "
    "momentum, once granted, is conserved.", S['plain']))
E.append(Paragraph(
    "Each instrument is described twice. <b>The Measured Account</b> speaks "
    "the language of fluid mechanics — pressures in pascals, angles in "
    "degrees, areas in square millimetres, every figure taken from the "
    "parametric models as built and verified. <b>The Hermetic Reading</b> "
    "speaks the older language of the seven principles of the Kybalion and of "
    "the Emerald Tablet, in which the same motions are read as doctrine. The "
    "second voice makes no claim upon the first: equations do not require "
    "blessing, and axioms do not require proof. They are two lanterns held to "
    "one whirlwind, and the reader may walk by either light, or by both.",
    S['plain']))

# ======================= OF VORTICES =======================
E.append(Paragraph("OF VORTICES IN GENERAL", S['h1']))
E.append(rule())
E.append(voice("THE MEASURED ACCOUNT"))
E.append(Paragraph(
    "A vortex is flow organised about an axis. Its measure is vorticity, "
    "%s = %s × <i>u</i> — the curl of the velocity field, the local rate at "
    "which the fluid turns about itself. Two ideal forms bracket every real "
    "whirl. In the <i>forced</i> vortex the fluid rotates as a rigid body, "
    "v = %sr, and must be continuously driven. In the <i>free</i> vortex, "
    "conservation of angular momentum sets v·r = constant, so fluid drifting "
    "inward is compelled to accelerate: a rim breeze of half a metre per "
    "second, carried from a 76 mm wall to a 5 mm core, arrives spinning some "
    "fifteen times faster. Every instrument in this folio is a machine for "
    "arranging that drift. Geometry grants the air its first angular "
    "momentum; convergence multiplies it; and because friction forever taxes "
    "the spin, each device must also supply a steady revenue of new momentum "
    "— from heat, from wind, or from breath." % (OMEGA, NABLA, OMEGA), S['sci']))
E.append(voice("THE HERMETIC READING"))
E.append(axiom("“As above, so below; as below, so above.”",
               "— The Kybalion, on the Principle of Correspondence (1908)"))
E.append(Paragraph(
    "Fluid mechanics states the second principle with unusual candour. The "
    "equations that govern a smoke-helix in a tube of 100 millimetres govern "
    "likewise the dust devil on the plate, the waterspout, the hurricane, "
    "and, in their gravitational translation, the spiral galaxy. Scale falls "
    "away in the dimensionless numbers — Reynolds, swirl — which are the "
    "algebra of correspondence: shrink the heavens to a desk and the form "
    "survives. Whoever has watched the small tornado stand in its tube has "
    "seen the shape of larger weathers; the pattern below answers the "
    "pattern above, and instructs in it.", S['herm']))

# ======================= INSTRUMENT I =======================
E.append(PageBreak())
E.append(Paragraph("INSTRUMENT I · THE HEARTH TORNADO", S['h1']))
E.append(Paragraph("“The Sun is its father”", S['h1sub']))
E.append(img("preview_assembly.png", 1.85))
E.append(Paragraph("PLATE I — The Hearth Tornado assembled: slotted drum, "
    "clear tube, and exit collar. Twelve chord-slots visible in the base.",
    S['caption']))
E.append(voice("THE MEASURED ACCOUNT"))
E.append(Paragraph(
    "The hearth tornado is a chimney taught to spin. A standard incense cone "
    "burns on a heat-shielded pedestal inside a 176 mm drum; the clear tube "
    "stands over it. The ember's few watts warm the column and the stack "
    "effect supplies the draft, %sP = %s g H (%sT/T) — roughly six-tenths of a "
    "pascal for a 394 mm tube. Feeble, but free. Make-up air can reach the "
    "plume only through twelve slots cut not on radii but on chords, tangent "
    "to an inner circle, so every entering filament arrives at 68 degrees "
    "from radial. Slot area is held at 0.76 of the bore area, resistance low "
    "enough that the whisper of draft survives. Converging from the 45 mm "
    "chamber wall to the 6 mm smoke core, the swirl multiplies about sevenfold, "
    "and the plume is wrapped into a standing clockwise helix. A 9 mm "
    "mouthpiece bore, coaxial with one slot, accepts the operator's breath as "
    "a tangential injection of angular momentum — the column visibly tightens "
    "at a puff. The pedestal is hollow, a sealed air gap against conduction, "
    "and its crown recess takes a metal disc, for the ember is hotter than "
    "PETG's comfort of about 80 °C." % (DELTA, RHO, DELTA), S['sci']))
E.append(voice("THE HERMETIC READING"))
E.append(axiom("“Nothing rests; everything moves; everything vibrates.”",
               "— The Kybalion, on the Principle of Vibration"))
E.append(Paragraph(
    "Fire is the father here in the Tablet's exact sense: the ember performs "
    "the first transmutation, dividing dense resin into volatile spirit, and "
    "the warmth of that division is itself the engine that carries the spirit "
    "upward. Nothing in the standing column is still. What appears a fixed "
    "and patient form — a rope of smoke holding its place in the tube — is "
    "motion all the way through; the shape persists only because movement is "
    "continually renewed, as the third principle insists of all seemingly "
    "solid things. And the instrument honours the first principle before the "
    "third: it existed complete upon the mental plane — a page of pure "
    "relations, radii and angles in a parametric file — before any matter "
    "condensed about the idea. The printer is merely the scribe.", S['herm']))

# ======================= INSTRUMENT II =======================
E.append(PageBreak())
E.append(Paragraph("INSTRUMENT II · THE WIND TORNADO", S['h1']))
E.append(Paragraph("“The Wind hath carried it in its belly”", S['h1sub']))
E.append(img("preview_wind_assembly.png", 1.85))
E.append(Paragraph("PLATE II — The Wind Tornado: pinwheel intake fins below, "
    "stacked-disc venturi cowl above, staked skirt at the foot.", S['caption']))
E.append(voice("THE MEASURED ACCOUNT"))
E.append(Paragraph(
    "Outdoors the budget improves by an order of magnitude. The ram pressure "
    "of moving air is q = ½ %s v² — about 5.4 Pa in a 3 m/s breeze, some "
    "nine times the hearth draft — and the wind edition spends that wealth "
    "twice. At the base, twelve tall guide fins radiate from the slotted drum "
    "as a fixed stator: from any compass direction the windward fins form "
    "converging funnels of roughly four-to-one contraction, accelerating the "
    "breeze and discharging it through the same 68-degree chords, so that all "
    "weather is translated into the one clockwise currency. At the crown, a "
    "stacked-disc venturi cowl squeezes crosswind through a narrowing gap "
    "above the exit throat; by Bernoulli's accounting fast air is cheap air, "
    "and the pressure minimum draws the column upward. Push below, pull above "
    "— one wind applied at both poles of the tube. The gap's perimeter area "
    "equals the throat area, so in a calm the cowl does not choke the chimney "
    "and the device degrades gracefully into its indoor cousin. In gusts the "
    "swirl briefly outruns the axial flow and the core shatters — vortex "
    "breakdown — then gathers and stands again. Three countersunk holes "
    "anchor the skirt; a tower of 450 mm keeps real company with the wind."
    % RHO, S['sci']))
E.append(voice("THE HERMETIC READING"))
E.append(axiom("“Everything flows, out and in; all things rise and fall; "
    "the pendulum-swing manifests in everything.”",
    "— The Kybalion, on the Principle of Rhythm"))
E.append(Paragraph(
    "This is the instrument of the fourth and fifth principles together. "
    "<b>Polarity:</b> ram and suction — the crowding of air below, the "
    "thinning of air above — are not two winds but one, identical in nature "
    "and differing only in degree; and the column stands precisely between "
    "the poles, as all manifest things do. <b>Rhythm:</b> the gust and the "
    "lull, the bursting and regathering of the core, are the pendulum made "
    "visible in smoke; the vortex does not fail when it breaks — it keeps "
    "time. And the Tablet's third clause is enacted without metaphor: the "
    "wind receives the smoke into its belly at the fins, gestates it in the "
    "spiral, and delivers it skyward through the cowl. What the operator once "
    "supplied with breath, the atmosphere now supplies with weather.",
    S['herm']))

# ======================= INSTRUMENT III =======================
E.append(PageBreak())
E.append(Paragraph("INSTRUMENT III · THE DUST-DEVIL PLATE", S['h1']))
E.append(Paragraph("“The Earth is its nurse”", S['h1sub']))
E.append(img("preview_dust_plate_persp.png", 2.65))
E.append(Paragraph("PLATE III — The Dust-Devil Plate: twelve-fin crown, "
    "removable roof ring, dished arena with the powder saucer at centre.",
    S['caption']))
E.append(voice("THE MEASURED ACCOUNT"))
E.append(Paragraph(
    "The third instrument removes the tube and asks the vortex to stand "
    "unclothed. A flat arena, 200 mm across and 79 mm tall assembled, is "
    "walled by a twelve-fin crown; wind from any direction is pressed through "
    "9 mm throats at 58 degrees and enters as tangential jets. The floor is "
    "dished — high at the rim, descending to a 40 mm spherical saucer at dead "
    "centre — and in the saucer waits the working powder: Borozin, micronised "
    "zinc stearate, an airflow tracer manufactured to float like smoke. The "
    "core's pressure minimum forms directly over the saucer; converging "
    "inflow scrubs the dish, entrains the powder, and renders the invisible "
    "whirl white. Centrifuged grains climb the dished floor, stall, and slide "
    "home to be lifted again. A removable roof ring keys onto the fin tops, "
    "keeping the wind committed to the channels; its 136 mm opening is the "
    "devil's door. Unconfined, the devil is intermittent by nature — it "
    "stands, wanders, collapses, and is reborn with the gusts.", S['sci']))
E.append(Paragraph(
    "Size is its ally, and the scaling laws are worth stating plainly: "
    "captured wind power grows as the square of diameter, circulation grows "
    "linearly, and resistance to gusts — the vortex's stored angular momentum "
    "— grows as the <i>fourth</i> power, which is why the devils of the "
    "desert outlive those of the desk. Peak velocity alone does not scale: it "
    "remains the wind's gift. A dark-coloured print in sunshine adds the "
    "convective updraft that feeds the wild kind. The folio accordingly "
    "holds the plate at two stations: the desk plate, 200 mm across and "
    "79 mm assembled, and the XL280 — 280 mm across with an 88 mm crown, "
    "99 mm assembled, the honest ceiling of the flat-form envelope, whose "
    "13 mm throats tax the entering wind noticeably less.", S['sci']))
E.append(img("preview_dust_plate_section.png", 2.9))
E.append(Paragraph("PLATE IV — Section through the arena: dished floor "
    "descending to the saucer; fins carrying the roof ring.", S['caption']))
E.append(voice("THE HERMETIC READING"))
E.append(axiom("“Every Cause has its Effect; every Effect has its Cause; "
    "Chance is but a name for Law not recognised.”",
    "— The Kybalion, on the Principle of Cause and Effect"))
E.append(Paragraph(
    "Earth is the nurse, says the Tablet, and here earth at its most finely "
    "divided is nursed into the air and weaned back again. The plate performs "
    "the Tablet's central circulation in miniature — <i>“It ascends from the "
    "earth to the heaven, and again it descends to the earth, and receives "
    "the strength of things above and below”</i> — for the powder rises white "
    "in the core, receives the strength of the wind above and the sun-warmed "
    "floor below, and settles to the saucer only to ascend again. The devil's "
    "wandering looks like caprice; the sixth principle corrects us. Every "
    "stagger of the vortex is written jointly in the geometry and in the "
    "gust, and what we call chance is turbulence we declined to compute. The "
    "seventh principle attends quietly: a projective stream and a receptive "
    "vessel, jet and dish — and generation only in their meeting. Neither "
    "alone has ever raised a single grain.", S['herm']))

# =================== THE SEVEN PRINCIPLES, MAPPED ===================
E.append(PageBreak())
E.append(Paragraph("THE SEVEN PRINCIPLES, MAPPED TO THE PHYSICS", S['h1']))
E.append(rule())
E.append(Paragraph(
    "The Kybalion's seven principles, each set beside its shadow in the "
    "measured behaviour of the three instruments — half doctrine, half "
    "datasheet, in equal standing.", S['plain']))
rows = [
    [Paragraph("<b>Principle</b>", S['cellb']),
     Paragraph("<b>The Axiom (Kybalion, 1908)</b>", S['cellb']),
     Paragraph("<b>Its shadow in the physics</b>", S['cellb'])],
    ["Mentalism", "“THE ALL is MIND; the Universe is Mental.”",
     "Each device existed first as pure relation — a parametric file of radii "
     "and angles — before matter condensed about the idea."],
    ["Correspondence", "“As above, so below; as below, so above.”",
     "Vortex self-similarity: one set of equations serves dust devil, "
     "waterspout, hurricane and galaxy; dimensionless numbers carry the form "
     "across scales."],
    ["Vibration", "“Nothing rests; everything moves; everything vibrates.”",
     "Vorticity is the curl of velocity: the standing column is organised "
     "motion, its stillness a persistence and not a rest."],
    ["Polarity", "“Everything is Dual; everything has poles.”",
     "All flow is driven by difference of pressure: ram against suction, rim "
     "against core. No gradient, no motion."],
    ["Rhythm", "“Everything flows, out and in; the pendulum-swing manifests "
     "in everything.”",
     "Gust and lull; vortex breakdown and re-formation; the devil that dies "
     "and is reborn keeps the measure."],
    ["Cause and Effect", "“Every Cause has its Effect; Chance is but a name "
     "for Law not recognised.”",
     "The fluid equations are deterministic: every wander of the whirl is "
     "written in geometry and boundary condition."],
    ["Gender", "“Gender is in everything; everything has its Masculine and "
     "Feminine Principles.”",
     "Projective and receptive: jet and chamber, updraft and dish — the "
     "helix is generated only in their union."],
]
data = [rows[0]] + [[Paragraph(r[0], S['cellb']), Paragraph(r[1], S['cellq']),
                     Paragraph(r[2], S['cell'])] for r in rows[1:]]
t7 = Table(data, colWidths=[1.18*inch, 2.32*inch, 2.9*inch], repeatRows=1)
t7.setStyle(TableStyle([
    ('GRID', (0,0), (-1,-1), 0.4, GOLD_L),
    ('BACKGROUND', (0,0), (-1,0), PARCH),
    ('VALIGN', (0,0), (-1,-1), 'TOP'),
    ('TOPPADDING', (0,0), (-1,-1), 4), ('BOTTOMPADDING', (0,0), (-1,-1), 4),
    ('LEFTPADDING', (0,0), (-1,-1), 5), ('RIGHTPADDING', (0,0), (-1,-1), 5),
]))
E.append(t7)

# =================== COMPARISON ===================
E.append(Paragraph("THE THREE INSTRUMENTS COMPARED", S['h1']))
E.append(rule())
comp = [
    ["", "I · Hearth Tornado", "II · Wind Tornado", "III · Dust-Devil Plate"],
    ["Clause of the Tablet", "The Sun is its father",
     "The Wind hath carried it in its belly", "The Earth is its nurse"],
    ["Driving power", "Ember heat (stack effect)",
     "Wind ram + venturi suction", "Wind ram (+ solar floor)"],
    ["Pressure budget", "~0.6 Pa", "~5 Pa at 3 m/s", "~5 Pa at 3 m/s"],
    ["Swirl generator", "12 chord slots at 68°",
     "12-fin stator, 4:1 funnels, 68°",
     "12-fin crown at 58°; throats 9 mm (XL: 13 mm)"],
    ["Visible working body", "Incense smoke", "Incense smoke",
     "Zinc stearate powder"],
    ["Confinement", "100 mm clear tube", "100 mm clear tube",
     "None — a free devil"],
    ["Stature, assembled", "~460 mm", "~490 mm", "79 mm · XL280: 99 mm"],
    ["Fair weather", "A still room", "Steady breeze, 1–4 m/s",
     "Sun and steady breeze, 2–5 m/s"],
]
cdata = [[Paragraph("<b>%s</b>" % c if i == 0 or j == 0 else c,
          S['cellb'] if i == 0 or j == 0 else S['cell'])
          for j, c in enumerate(row)] for i, row in enumerate(comp)]
tc = Table(cdata, colWidths=[1.35*inch, 1.68*inch, 1.68*inch, 1.68*inch],
           repeatRows=1)
tc.setStyle(TableStyle([
    ('GRID', (0,0), (-1,-1), 0.4, GOLD_L),
    ('BACKGROUND', (0,0), (-1,0), PARCH),
    ('BACKGROUND', (0,1), (0,-1), PARCH),
    ('VALIGN', (0,0), (-1,-1), 'TOP'),
    ('TOPPADDING', (0,0), (-1,-1), 4), ('BOTTOMPADDING', (0,0), (-1,-1), 4),
    ('LEFTPADDING', (0,0), (-1,-1), 5), ('RIGHTPADDING', (0,0), (-1,-1), 5),
]))
E.append(tc)

# =================== CAUTIONS ===================
E.append(KeepTogether([
    Paragraph("THE CAUTIONS", S['h1']), rule(),
    Paragraph("Here the treatise sets aside the second voice and speaks "
              "plainly, once.", S['plain'])]))
for b in [
    "<b>The ember is hotter than the vessel.</b> PETG softens near 80 °C. "
    "The metal disc or foil in the pedestal recess is not optional; the "
    "hollow pedestal is a heat break, not a heat shield.",
    "<b>Never leave a burning cone unattended.</b> Light it with the tube "
    "lifted away; lower the tube afterwards.",
    "<b>Dispersed zinc stearate is a combustible dust.</b> The dust plate "
    "and the flame-bearing instruments must never share a table. Keep the "
    "powder cloud away from every ignition source.",
    "<b>Do not breathe the cloud</b> — neither concentrated smoke nor "
    "airborne powder.",
    "<b>Anchor the outdoor instruments.</b> The stake holes are provided "
    "because a tower in wind, and a plate beside a flame, both deserve "
    "respect for the Principle of Cause and Effect.",
]:
    E.append(Paragraph("· " + b, S['plain']))

# =================== COLOPHON ===================
E.append(Paragraph("COLOPHON", S['h1']))
E.append(rule())
E.append(Paragraph(
    "Instrument I — vortex_base.stl, vortex_top_collar.stl "
    "(source: incense_vortex.scad; tube: 100 mm OD x 2 mm acrylic, 393.7 mm). "
    "Instrument II — wind_base.stl, wind_venturi_head.stl, "
    "wind_venturi_hat.stl (source: incense_vortex_wind.scad). "
    "Instrument III — dust_plate_arena.stl, dust_plate_roof.stl; XL280 "
    "variant — dust_plate_arena_XL280.stl, dust_plate_roof_XL280.stl "
    "(source: vortex_dust_plate.scad, rendered at two parameter stations).",
    S['colo']))
E.append(Paragraph(
    "All solids were modelled parametrically in OpenSCAD 2021.01 (CGAL "
    "kernel) and every exported mesh verified watertight and two-manifold "
    "with pymeshlab. Every vortex in the family turns clockwise when viewed "
    "from above. Printed matter: PETG or ABS, three perimeters, no supports "
    "required.", S['colo']))
E.append(Paragraph(
    "The axioms are quoted from public-domain sources: The Kybalion (Yogi "
    "Publication Society, 1908) and the Emerald Tablet in Isaac Newton's "
    "translation. The physics is quoted from the instruments themselves.",
    S['colo']))
E.append(Spacer(1, 10)); E.append(rule(1.4*inch))
E.append(Paragraph("As above, so below — and at 68 degrees between.",
    S['attrib']))

doc.build(E, onFirstPage=title_page, onLaterPages=later_pages)
print("built:", OUT)
