# The Orifice Grimoire

*New mouths for the Incense Vortex Ring, read through Vega's "Dynamics of Vortex Filaments
with Corners" (Sapporo, 2014) and de la Hoz & Vega, "Vortex Filament Equation for a Regular
Polygon" (2013). Written in the voice of a smoke sorcerer who would rather be right than
mystical, and keeps the two apart with labels.*

Every claim below carries one of three marks:

- **[F]** fluid dynamics: published, measured, or follows from the binormal law.
- **[S]** speculation: a prediction of mine that nobody has tested on this instrument.
- **[C]** capnomancy: the reading. Pseudoscience, kept here on purpose and never allowed to
  leak into [F].

The instrument is unchanged: Ø110 membrane, 12 mm stroke, 86 mL slug, 100 mm × 2 mm acrylic
tube (96 mm bore, 394 mm tall), drop-in orifice inserts with a Ø56 body under the lid. Only
the **shape of the mouth** changes. See [README.md](README.md) for the build.

---

## 1. What the paper lets us predict [F]

Vega's model is the **binormal flow**, X_t = X_s ∧ X_ss: every point of a thin vortex
filament moves along its binormal at a speed equal to its curvature. A smoke ring is such a
filament. A circle has constant curvature, so the whole ring moves together and the circle
persists. A **corner** has infinite curvature at a point, and the paper's whole point is what
a corner does next:

1. **A corner does not stay a corner.** It becomes a self-similar hairpin, X = √t·G(s/√t),
   whose tip runs at speed a/√t: fastest at birth, then slowing. In a real fluid viscosity
   rounds it off at once, which only hurries the same process.
2. **A regular M-gon is a Talbot carpet.** Start the filament as a planar regular polygon
   with M corners. In the model's units the motion is periodic with period **T = 2π/M²**, and
   at every rational time t = T·(p/q) the curve is again a polygon, now a *skew* one with
   **Mq sides** (q odd) or **Mq/2 sides** (q even). The corner angles come from quadratic
   Gauss sums. At q = 2, half the Gauss sums vanish, so the curve has only M sides again:
   **at T/2 the original polygon reappears, rotated about the axis by π/M.** Triangle: 60°
   (point-up becomes point-down). Square: 45°. Pentagon: 36°. Hexagon: 30°.
3. **The whole thing climbs.** The polygon's centre of mass rises at a constant speed c_M
   that grows with M and tends to 1 (the circle's speed) as M → ∞. A triangle ring is the
   slowest ring of all; a many-sided one is nearly a circle.
4. **Higher polygons blur faster.** The period scales as 1/M², so a hexagon runs through its
   whole cycle four times faster than a triangle. In smoke, where every cycle costs
   definition, the triangle is the most legible shape and the one that holds its character
   longest.
5. **A corner draws a fractal.** The path of one corner, projected on a plane, resembles
   Riemann's non-differentiable function ∑ e^{−2πi t k²}/k², proved multifractal by Jaffard.
   The more sides, the closer the resemblance.
6. **This is seen in real fluids.** Axis switching of non-circular jets is documented for
   elliptic, rectangular, square and triangular nozzles (Gutmark & Grinstein's review; Toyoda
   1992 for square/triangle/rectangle). For a square orifice the first switch happens at
   about **0.7 orifice diameters** of travel. Vega's own cardboard smoke cannon with a
   triangular hole showed the triangle flipping. He also reports that smooth polygons with
   *more* angles than the nozzle have been photographed for triangular jets (the Talbot
   skew polygons), and calls their study worthwhile. That is an open invitation.

What the model does **not** do: it has no viscosity, no core thickness, no tube wall. The
tube wall is an image vortex that slows the ring (~11 % for a Ø36 ring in the 96 mm bore)
and will bias any corner that comes within a few mm of the acrylic. Expect the sequence of
shapes to be right and every number to be off.

### The sorcerer's three laws [S, distilled from F]

- **Convex corners run ahead; concave corners run behind.** Curvature times binormal: where
  the outline bends the other way, the velocity flips sign relative to the ring's travel. A
  star-shaped ring is therefore a tug of war between its points and its notches.
- **Regular returns; irregular never returns.** The Talbot recurrence needs the M-fold
  symmetry. Break it (unequal sides, unequal angles) and the Gauss-sum cancellations fail:
  the shape evolves without ever reassembling itself.
- **The fewer the corners, the longer the memory.** Period ∝ 1/M².

---

## 2. The orifice catalogue

All openings are sized to the **same area as a Ø30 circle (707 mm²)** unless stated, so the
86 mL slug keeps the same length and the ring carries the same impulse. Fit limit: the insert
body is Ø56, so the opening's circumscribed circle must stay ≤ Ø46 to leave a 5 mm wall. The
ring born from a shape is ~1.2× its size, and the README's rule still holds: **keep the ring
under about half the 96 mm bore**, or the wall image eats it.

The slug length rule changes with shape. Use the hydraulic diameter D_h = 4A/P: a cornered
opening has more perimeter per area, so the full 12 mm press gives a *longer* slug in D_h
terms. The right stroke for L/D_h ≈ 4 is in the table.

| Mouth | Dimensions (equal area to Ø30) | Circum. Ø | D_h | Full-press L/D_h | Press to |
|---|---|---|---|---|---|
| Circle (reference) | Ø30 | 30.0 | 30.0 | 4.0 | full |
| Equilateral triangle | side 40.4 | 46.7 | 23.3 | 5.2 | ~75 % |
| Square | side 26.6 | 37.6 | 26.6 | 4.6 | ~85 % |
| Pentagon | side 20.3 | 34.5 | 27.9 | 4.4 | ~90 % |
| Hexagon | side 16.5 | 33.0 | 28.6 | 4.3 | full |
| Ellipse AR 2 | 42.4 × 21.2 | 42.4 | 27.5 | 4.4 | ~90 % |
| Ellipse AR 3 (Ø25 area) | 43.3 × 14.4 (Ø30 area would be 52 wide: too big) | 43.3 | 20.4 | 8.6 | ~45 % |
| Slot AR 5 (Ø25 area) | 50.6 × 10.1, rounded ends | 50.6 (marginal) | 17.4 | 10.1 | ~40 % |

The triangle's circumscribed Ø46.7 is right at the fit limit and its ring (~Ø56 across the
points) is over half the bore. **Print the triangle at Ø25-equivalent area** (side 33.7,
circum. Ø38.9, D_h 19.4, full-press L/D_h 9.0, so press ~45 %) as the first one; the
Ø30-equivalent only if the small one is too faint. The Ø25-area shapes all want a half
press: the 86 mL slug is simply too long for a small cornered mouth, and an over-long slug
pinches off a trailing jet that smears the ring.

### Tier I: the Vega polygons [F predictions, S for this instrument]

**△ The Triangle: "the Flip."** The slowest, most legible ring, with three hairpin corners.
Prediction: the triangle rises point-forward, passes through a skew hexagon (t = T/4: q = 4
→ Mq/2 = 6 sides) and a skew 9-gon (t = T/3), and at T/2 comes back as a triangle turned 60°,
i.e. upside down. Then the sequence runs in reverse. Vega photographed the flip and could not
confirm the hexagon; the hexagon is the prize.

**□ The Square: "the Diamond."** Flips 45° at T/2 (square → diamond → square), passing an
octagon at T/4. Axis switch measured at ~0.7 D downstream in square jets, so with D_h 27 mm
the first diamond should appear ~20 mm above the lid, within the first 5 % of the tube. The
"eightfold pattern" reported for square jets is the T/4 octagon.

**⬠ The Pentagon** and **⬡ the Hexagon.** Included for completeness and for the reading; at
M = 5 and 6 the flips are 36° and 30° and the periods are 2.8× and 4× shorter than the
triangle's. Expect them to look like slightly lumpy circles that wobble. If you only print
three inserts, skip these.

### Tier II: ellipses [F]

**⬭ The Ellipse AR 2: "the Breath."** The classic axis-switching ring (Dhanak & de
Bernardinis 1981): the long axis becomes the short axis and back again, a visible oscillation
rather than a flip. Below a critical aspect ratio the ring just oscillates, with a period that
lengthens as AR → 1.

**The Slot AR 5: "the Cleaving."** Above the critical aspect ratio (roughly 3.5–4 for
isolated rings; synthetic-jet work at AR ≥ 13 saw only one switch before breakup) the ring
switches once, collides with itself, and **bifurcates into two or three rings** that drift
apart. One press, two rings. Fit is marginal at Ø25-equivalent area; shorten the slot to
AR 4 if it rubs.

### Tier III: cornered and cusped [S]

**✦ The Pentagram outline (a star, five points, five notches): "the Ring that Eats Itself."**
Five convex corners run ahead, five concave corners run behind (first law). The filament
folds across itself within a diameter or two and reconnects. I expect the result to be either
one ragged ring with five trailing hairpins, or a central ring plus five small satellite rings
spat backward. Not a stable sigil, which is the point of trying it.

**☾ The Crescent / Vesica.** Two cusps and one concave arc. The convex side leads, the
concave side lags: the crescent should fold toward its own horns. The **vesica piscis** (two
cusps, both sides convex) is a lens: a 2-gon in Vega's language, and the simplest cornered
shape that is symmetric enough for the recurrence to hold. Untested anywhere that I know of.

**◎ The Annulus (a ring-shaped slot around a central plug).** Two coaxial vortex rings,
the inner one born with opposite sense. Known from annular-jet work to produce a leading ring
chased by a counter-rotating one that either merges or stalls the first. Needs the plug held
by three thin spokes; spokes are themselves three small corners.

**∞ Twin circles (two Ø20 holes, 26 mm apart): "the Lovers."** Two rings side by side
attract, touch, reconnect into one peanut-shaped ring, which then axis-switches like an
ellipse. The single best-documented reconnection you can make with a hole; the README's
leapfrogging is the coaxial version of the same thing.

**⟋ The Oblique Mouth (a Ø30 bore whose exit face is cut at 20°).** Not a shape but a
tilt: the ring is born inclined and curves its path toward the tube wall instead of rising
straight. In the tube it will scrape. Mentioned so you don't print it by accident.

### What a hole cannot do [F]

A hole can only birth a filament that is **planar** at t = 0. Knotted or linked rings (the
Kleckner & Irvine trefoil) need a *hydrofoil* swept through the fluid, whose trailing edge
traces the knot in 3D, and in water the trefoil unties itself by reconnection within a few
diameters anyway. The nearest the ring instrument can come is the skew polygons that the
Talbot effect generates for free: the T/4 hexagon of a triangle is a genuinely non-planar
curve, "tied" by the dynamics rather than by the mouth.

---

## 3. Can you form a personal sigil?

Three different answers.

**As a shape in the air [F]:** yes, for about a third of an orifice diameter. The vortex
sheet rolls up along the whole edge of any opening, so the ring at birth *is* the outline,
rounded by the core radius (a few mm). After ~0.3 D the corners have become hairpins; after
~0.7 D the first axis switch has happened. For a 30 mm mouth that is 10 mm and 20 mm of
rise. The sigil exists for a few hundredths of a second as a readable figure, inside the lid
collar, where you cannot see it. A 240 fps phone looking straight down the tube can.

**As an irregular polygon [S]:** this is the real sorcerer's craft. The paper says the
future of a cornered ring is fixed entirely by its **corner angles and their spacing** (the
θ_m structure), not by any drawn line. A sigil that is a simple closed polygon with 4–7
corners, mostly convex, is a set of initial conditions; the ring then "computes" it. By the
second law an irregular sigil **never returns** to its initial shape, while a regular one
returns inverted. The reading writes itself: *a question put as a regular figure comes back
reversed; a question put as a personal figure is answered once and dissolved.*

**As Spare's glyph [S + C]:** Austin Osman Spare's method was to write the desire as a
sentence, strike the repeated letters, fuse the survivors into one glyph that no longer
suggests the sentence, charge it, and forget it. The sigil's job is to be *forgotten by the
conscious mind*. Smoke is the perfect medium: the glyph is formed, flies, and un-forms itself
within a second, with no trace left to remember. The orifice insert is the only durable copy,
and it is a negative.

### A sigil that the air will carry: constraints for the mouth [F]

1. One closed, **non-self-intersecting** outline. (A self-crossing sigil is two openings that
   touch; the air treats them as Twin Circles.)
2. **≤ 7 corners**, and no more than 2 of them concave (reflex). Each concave corner is a
   piece of ring that runs backward; three or more and the sigil tears on the lid.
3. Smallest inscribed width **≥ 8 mm** everywhere, or that limb necks off into its own ring.
4. Area between the Ø25 and Ø35 circles (491–962 mm²); circumscribed circle ≤ Ø46.
5. No interior islands unless you accept three spokes (and their three corners).
6. Corners print sharp enough: a 0.4 mm nozzle leaves ~0.3 mm radius, far below the core
   radius. The air will never know the difference.

Pipeline: Spare the sentence → draw the glyph → reduce it to a polygon that satisfies 1–5 →
scale to the Ø30-equivalent area → generate the insert. The generator is the follow-up: an
`orifice_sigil.py` that takes a polygon (list of points or an SVG path), checks the six
constraints, scales to area, and emits an insert with the standard Ø70 flange, a 45° inlet
chamfer in place of the bellmouth (a bellmouth needs a revolve; a chamfer works on any
outline), and the square exit edge the README insists on.

---

## 4. Capnomancy: the old readings

The practice is older than the word. What the sources actually say, with their weight:

- **Babylon.** Smoke divination (*libanomancy*, from Greek *libanos*, frankincense) is the
  earliest attested form. Shulgi of Ur (c. 2094–2047 BCE) boasted of mastering it; the
  specialist was the *bārû ša qutrenni*, the "diviner of the smoke." The surviving manuals
  are Old Babylonian (c. 2000–1600 BCE); the best edition of a new tablet is Finkel, "A New
  Piece of Libanomancy," *AfO* 29–30 (1983–84) 50–55, with Cingolo's 2016 remarks on the
  Old Babylonian texts. The Akkadian *muššakkum* is "incense used for libanomancy" (CT 15,
  5f.). The omens are protasis-apodosis like all Mesopotamian omens: *if the smoke* does X,
  *then* Y. The popularly quoted readings (smoke drifting right = victory over the enemy,
  smoke splitting = a divided outcome) come through a modern occult compilation, not from the
  tablets I could check: treat them as the folk layer.
- **Greece.** Sacrificial smoke was read at the altar. The rule the 1797 *Britannica* still
  repeats: thin, light smoke rising straight is good; the contrary is ill. The great literary
  case is Tiresias in Sophocles' *Antigone*: the fire will not blaze, the fat spits and
  smokes instead of burning, and he reads it as the gods refusing the offering because of
  Creon's crimes. A ring that stalls at the lid is this omen.
- **The Hebrew Bible.** The smoke of incense is a veil, not a message: on the Day of
  Atonement the priest must fill the inner sanctuary with incense smoke so that its cloud
  covers the mercy seat and he does not die (Leviticus 16:13). In Revelation 8:4 the smoke of
  the incense rises with the prayers of the saints. Smoke that rises straight is a prayer
  received; smoke that covers is a mercy.
- **The Celts** are said to have read the smoke of oak and mistletoe (*dendromancy*). Thin
  evidence; included because it is the only tradition that reads smoke by what was burnt.

None of these traditions read **shape**. They read direction, density, straightness,
division, and whether the fire took. That is exactly what the instrument can vary with a
mouth: straightness (the ring's wobble), division (bifurcation), and whether it takes
(formation at all). The old readings transfer better than they have any right to.

### The reading table [C]

A divination layer on deterministic fluid dynamics. The honest part: Gutmark & Grinstein
note that axis switching is *sensitive to the initial geometric conditions*, and the press is
a human hand. The oracle's randomness is real; its meaning is yours.

| What the smoke does | Fluid cause [F] | Reading [C] |
|---|---|---|
| Ring rises straight, tight, and leaves the tube whole | good slug, L/D ≈ 4 | The straight smoke of the altar: favourable, the question is well put |
| Ring stalls at the lid, smears | sharp jab, or weak smoke | Tiresias' altar: the offering is refused; ask again later |
| Triangle returns inverted | T/2 axis switch | Reversal: the answer is the opposite of what was asked |
| The hexagon appears between triangles | T/4 skew polygon | A hidden structure doubles the visible one |
| Square becomes diamond | 45° switch | The same thing seen from another side; reframe |
| Ellipse breathes | oscillatory axis switching | Alternation; neither outcome holds |
| Slot cleaves into two rings | bifurcation above critical AR | Two paths, and they diverge |
| Twin rings become one | reconnection | Union |
| Second ring passes through the first | leapfrogging | The later overtakes the earlier |
| Star folds into itself, sheds satellites | concave corners run backward | The question is ill-formed; the wish fights itself |
| Irregular sigil dissolves and does not return | no Talbot recurrence | Answered once; do not ask twice |
| Regular sigil returns | Talbot recurrence | The question will come back |
| Ring touches the wall and dies | image vortex | The asker's circumstances constrain the answer |

---

## 5. Observing it: the part that will actually decide whether any of this is seen

- **Look down the tube, not across it.** Every event above is a rotation or deformation
  *about the axis*. A side view of a triangle ring and its 60°-flipped twin look identical.
  A phone on a stand above the tube, looking straight down the bore, or a mirror at 45° on
  the rim, is the whole experiment. Side lighting, black cloth behind.
- **Slow the ring.** The first switch is at ~0.7 D ≈ 20 mm above the lid, inside the first
  5 % of the tube. A slow press (0.4 m/s rise) gives ~0.05 s per diameter; 240 fps is 12
  frames per diameter. A jab gives 1–2 frames.
- **Fire into the free air first.** The tube wall will blur the result; a bare lid on a
  still day, with the tube off, is the control. Then add the tube and see what the image
  vortex does to the sequence.
- **Score blind.** Decide in advance what counts as "the triangle came back inverted," and
  have the person pressing the key not be the person scoring the frames. The oracle is
  allowed to be unfalsifiable; the fluid dynamics is not.

## 6. Print order

1. **Triangle, Ø25-equivalent** (side 33.7 mm). The flip and the possible hexagon.
2. **Square, Ø30-equivalent** (side 26.6 mm). The fastest, best-documented switch.
3. **Ellipse AR 2** (42.4 × 21.2 mm). The breath.
4. **Slot AR 4–5, Ø25-equivalent.** The cleaving.
5. **Twin Ø20.** The lovers.
6. Then the sigil generator, and whatever you write with it.

All are the same Ø70 flange, Ø56 body, 12 mm re-entrant collar and square exit edge as the
round inserts, printed flange-down with no supports. Four of them cost one evening of
printing.

---

## Sources

Physics: Vega, *The dynamics of vortex filaments with corners* (Sapporo, Aug 2014; the PDF
you gave me); [de la Hoz & Vega, *Vortex Filament Equation for a Regular Polygon*,
arXiv:1304.5521](https://arxiv.org/pdf/1304.5521); [*Multifractality and polygonal vortex
filaments*, arXiv:2412.04926](https://arxiv.org/pdf/2412.04926); [axis switching of
pulsating jets from rectangular nozzles, arXiv:2307.11266](https://arxiv.org/pdf/2307.11266)
(cites Toyoda 1992 for square/triangle/rectangle and Grinstein & DeVore 1996 for the 0.7 D
switch); [*Vortex rings from notched nozzles*, arXiv:1708.05639](https://arxiv.org/pdf/1708.05639)
(Dhanak & de Bernardinis, elliptic ring axis switching, crosslinks, bifurcation);
[Straccia & Farnsworth, Phys. Rev. Fluids 6, 054702](https://link.aps.org/accepted/10.1103/PhysRevFluids.6.054702)
(rectangular synthetic jets, AR 2–19, bifurcation after one switch);
[elliptical ring holography, arXiv:2310.19730](https://arxiv.org/pdf/2310.19730);
[Kleckner & Irvine, knotted vortices (ar5iv)](https://ar5iv.labs.arxiv.org/html/1310.3321).

Capnomancy: [Libanomancy (Wikipedia)](https://en.wikipedia.org/wiki/Libanomancy);
[Capnomancy (Wikipedia)](https://en.wikipedia.org/wiki/Capnomancy);
[Divination in Mesopotamia (Shulgi, bārû ša qutrenni)](https://uomus.edu.iq/CollegeEn/Details/87596);
[Hallo, *Before Tea Leaves*, BAR 31:2 (2005)](https://library.biblicalarchaeology.org/?p=71182);
[Finkel, *A New Piece of Libanomancy*, AfO 29–30 (KeiBi record)](https://keibi.uni-tuebingen.de/Record/KEI00064283?lng=en);
[Cingolo, *Some Remarks about the Old Babylonian Libanomancy Texts* (2016)](https://ocb.uni-tuebingen.de/Record/KEI00129049?lng=en);
[SEAL: *muššakkum*](https://seal.huji.ac.il/taxonomy/term/5996);
[*Britannica* 3rd ed. (1797), "Capnomancy"](https://jimclifford.ca/early_encyclopedia_britannica/articles/1797/eb_3rd_1797_005246.html);
[Sophocles, *Antigone* 991–1048 (Perseus)](https://cts.perseids.org/read/greekLit/tlg0011/tlg002/perseus-eng2/991-1048).

Sigils: [Spare's method, quoted at sacred-texts (TOPY)](https://sacred-texts.com/eso/topy/black3.htm);
[Frater U∴D∴, *Practical Sigil Magic* (review)](https://paralibrum.com/reviews/practical-sigil-magic).
