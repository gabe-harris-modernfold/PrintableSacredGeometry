# Crystal systems, Bravais lattices and Wigner-Seitz cells

Printable models of the 7 crystal systems, the 14 Bravais lattices and the three cubic
Wigner-Seitz cells, dimensioned from measured lattice constants of real minerals rather
than from idealised ratios.

| set | what | scale | parts |
|---|---|---|---|
| `stl/supercell/` | **each lattice as a 2×2×2 block of cells** | cell edge 50 mm | 14 |
| `stl/solids/` | 14 solid unit cells, Pearson symbol embossed | longest cell edge 60 mm | 14 |
| `stl/assembled/` | one cell each, one-piece ball-and-stick | longest cell edge 60 mm | 14 |
| `stl/kit/` | the same frames as separate nodes and rods | longest cell edge 90 mm | 387 |
| `stl/wigner_seitz/` | cube, truncated octahedron, rhombic dodecahedron | a = 60 mm | 3 |

**Start with `supercell/`.** A single unit cell cannot show periodicity, and periodicity is
the entire content of the word *lattice*. See below for what a block shows that one cell
cannot. `assembled/` is the same object as one cell, useful as the reference the block
repeats. `kit/` is the take-apart version — 387 pieces, and the only thing it buys is that
cells detach and rejoin into other arrangements (three hP cells make the hexagonal prism).

Previews are in `previews/`; the full part list, staged into build plates, is in
[PARTS.md](PARTS.md).

```bash
python build.py && python verify.py && python preview.py
```

## Why a block and not one cell

One cell shows the *shape* of the repeating unit. It cannot show the thing that makes a
lattice a lattice, and three facts only become visible once cells sit next to each other:

**Shared corners are why counting works.** A block of 8 cubic P cells has 27 hubs, not 64.
Merging is not an optimisation — it is the geometry. An interior corner is one lattice point
belonging to eight cells at once, and 8 corners × ⅛ = 1 is why a primitive cell contains a
single lattice point despite having eight corners. The same counting gives cI two points per
conventional cell (8 × ⅛ corners + 1 body centre, shared with nobody) and cF four
(8 × ⅛ + 6 × ½, each face centre shared by two cells). Those are the numbers in the hub
column of [PARTS.md](PARTS.md), arrived at by construction rather than asserted.

**A centring point is not a special kind of point.** In the cI block the body-centring hubs
and the corner hubs end up with identical surroundings — each has eight nearest neighbours at
the same distance. Which hubs you call "corners" is a choice of where to draw the cell, not a
property of the lattice. This is the same fact that makes monoclinic I and monoclinic C one
lattice, and it is invisible in a single cell, where the centre looks like an interloper.

**The pattern continues.** A lone cell reads as an object; a block reads as a fragment of
something unbounded, which is what a crystal is.

The cell at the origin is printed **1.35× heavier** in both hub and strut diameter, so the
repeating unit is identifiable inside the block by eye and by touch. At 2×2×2 it always
lands on a corner of the block, so it is visible from outside.

`SUPER_N` in `build.py` sets the block size. It is 2 by default because 3 is *almost*
possible: at 3×3×3 thirteen of the fourteen still fit, and only corundum busts the bed, at
330 mm against 320 — its rhombohedron is elongated 2.5 : 1 along the 3-fold axis, so its
block is by far the longest. If you want 3×3×3, set `SUPER_N = 3` and either drop
`SUPER_LONGEST` to about 48 mm or accept that hR is the one you cannot print; the build
reports the overrun rather than silently writing an unprintable file.

### Hub size

Hubs are Ø9.8 on the frames and Ø7.7 on the blocks, both 30 % down from where they started.
Two things move when you change them, and `build.py` prints the first at startup:

- **A strut's end cap has to finish inside its hub.** The cap sits at the strut inset from
  the node centre but its *rim* is `hypot(inset, strut_radius)` out, so it reaches further
  than the inset suggests. Shrink past that and the flat disc pokes through the ball as a
  visible ring. `_hub_hides_cap` raises rather than letting it through; current clearance is
  0.70 mm on the frames and 0.45 mm on the blocks. At these strut diameters the floor is
  about **−43 % for the frames and −39 % for the blocks** — 30 % is comfortable, another
  10 % would need thinner struts too.
- **Bridges get longer.** The free span is hub-centre spacing less the two hubs it anchors
  into, so smaller hubs mean more air: the frames went from 46 mm to 50.2 mm and the blocks
  from 39 mm to 42.3 mm. Still fine in PETG at Ø5–6, just worth knowing that slimming the
  nodes is not free.

## Reference card

Pearson symbol, the mineral each cell is measured from, and the conventional cell. The
centring letter is the modern IUCr **S** for base-centred, historically written **C**.

| Pearson | system | centring | mineral | a, b, c (Å) | angles | space group |
|---|---|---|---|---|---|---|
| aP | triclinic | primitive | turquoise | 7.48, 9.95, 7.68 | 111.65°, 115.38°, 69.43° | P1̄ |
| mP | monoclinic | primitive | azurite | 5.008, 5.844, 10.336 | β = 92.333° | P2₁/a |
| mS | monoclinic | base | gypsum | 6.29, 15.18, 5.68 | β = 113.833° | C2/c |
| oP | orthorhombic | primitive | topaz | 4.6499, 8.7968, 8.3909 | 90° | Pbnm |
| oS | orthorhombic | base | cordierite | 17.13, 9.80, 9.35 | 90° | Cccm |
| oI | orthorhombic | body | hemimorphite | 8.37, 10.719, 5.12 | 90° | Imm2 |
| oF | orthorhombic | face | natrolite | 18.27, 18.587, 6.56 | 90° | Fdd2 |
| tP | tetragonal | primitive | rutile | 4.594, —, 2.958 | 90° | P4₂/mnm |
| tI | tetragonal | body | zircon | 6.604, —, 5.979 | 90° | I4₁/amd |
| hP | hexagonal | primitive | beryl | 9.215, —, 9.192 | γ = 120° | P6/mcc |
| hR | rhombohedral | primitive | corundum | 5.120 (rhombohedral) | α = 55.286° | R3̄c |
| cP | cubic | primitive | pyrite | 5.417 | 90° | Pa3̄ |
| cI | cubic | body | almandine | 11.526 | 90° | Ia3̄d |
| cF | cubic | face | diamond | 3.5668 | 90° | Fd3̄m |

### Two things the usual table gets wrong

**Trigonal is not a lattice type.** The familiar row "trigonal: a = b = c, α = β = γ ≠ 90°,
primitive" describes the *rhombohedral* setting of the R lattice, and only some trigonal
minerals use it. Trigonal crystals sit on either the hexagonal lattice (**hP**) or the
rhombohedral one (**hR**). Quartz is trigonal — space group P3₁21 — but its lattice is
hexagonal primitive, the same lattice as beryl. Corundum, tourmaline and calcite are the
ones that are really rhombohedral. Getting this right is what makes the count come to 14:

    1 aP + 2 mP mS + 4 oP oS oI oF + 2 tP tI + 1 hP + 1 hR + 3 cP cI cF = 14

**Monoclinic I is monoclinic C.** Gypsum is the demonstration. Its literature refinement is
I2/a, webmineral lists it as A2/a, and the conventional setting is C2/c — three labels, one
Bravais lattice, because a body-centred monoclinic cell is a base-centred one in a
different basis. That is exactly why the monoclinic system has two lattices and not three,
and it is the same fact that makes a "centred" point in cI or cF an ordinary corner once
you change basis.

Also: the `≠` signs in the textbook table are the *absence* of an imposed constraint, not a
requirement. A triclinic cell with a = b by coincidence is still triclinic.

## Sources

Cell constants come from `webmineral.com/data/<Mineral>.shtml`, verified 2026-08-13, with
three deliberate departures:

- **topaz** — webmineral lists a = 4.35 Å, a typo. Ribbe & Gibbs, *The crystal structure of
  topaz and its relation to physical properties*, Am. Mineral. **56** (1971) 24, gives
  a = 4.6499, b = 8.7968, c = 8.3909; those are used.
- **gypsum** — webmineral's A2/a cell has a and c swapped relative to the conventional
  C-centred setting. Swapping them puts the centring vector at (½, ½, 0) like every other
  mS in the set. Cross-checked against the literature I2/a cell (a = 5.679, b = 15.202,
  c = 6.522, β = 118.43°), which is the same lattice.
- **corundum** — quoted in the hexagonal setting (a = 4.751, c = 12.97). `lattices.py`
  converts it to the primitive rhombohedron the model shows, via
  `a_rh = √(a²/3 + c²/9)` and `sin(α/2) = 3 / (2√(3 + (c/a)²))`, giving a = 5.120 Å,
  α = 55.286°. The check that this is right: the rhombohedron's body diagonal comes out at
  12.97 Å, i.e. exactly the hexagonal c.

webmineral's space-group column has typos that do not affect the lattice type — rutile
`P4/mnm` for P4₂/mnm, beryl `P6/mmc` for P6/mcc, turquoise `P1` alongside the
centrosymmetric pinacoidal class. The symbols in the reference card are the corrected ones.
Turquoise also has a second cell in circulation (a = 7.424, b = 7.629, c = 9.910,
α = 68.61°, β = 69.71°, γ = 65.05°): same lattice, different axis choice.

## Printing

PETG on a 320 × 320 × 320 mm bed. Largest single part is the corundum solid at
128 × 67 × 46 mm.

**Solids and Wigner-Seitz cells** — the STLs are solid, so hollow them in the slicer:
3 walls, 10 % gyroid, no support needed for any of the 17. The shallowest downward surface
among them is the truncated octahedron's at 55° (its {111} faces, arccos(1/√3) = 54.74°),
and the shallowest cell is gypsum's at 66° — both well clear of the 45° limit. Printed
solid these 17 parts would be 2.8 kg; at 3 walls and 10 % infill the 14 cells come to about
605 g and the three Wigner-Seitz cells about 125 g.

All five sets together are roughly 2.9 kg. `verify.py` prints the breakdown; its two infill
factors (0.75 for the strut sets, 0.26 for the bulky ones) are estimates, not slicer output.

**Repeating blocks** — print solid. Ø7.7 hubs and Ø5 struts, slimmer than a single cell or the
block becomes a lump. Largest is corundum at 223 × 121 × 85 mm; the cubic ones are ~108 mm
cubes, and the set is about 1115 g at 3 walls + 10 % infill, roughly 80 g each.
Orientation is chosen the same way as for the single frames, and the numbers are
larger for the same reasons: 24 bridges for most lattices (every horizontal strut above the
bottom level), 56 for the face-centred pair. **oF natrolite and cF diamond have 8 islands
each and need support**; the other twelve have none. The islands are the face-centring hubs
lying in horizontal planes above the bed, reached only by struts within their own plane —
one set per level, which is why 8 and not 4.

**Fused frames** — Ø6 struts and Ø9.8 hubs, so three walls leaves them mostly solid and infill
saves little: about 340 g for all 14, roughly 24 g each. Largest is corundum at
138 × 77 × 54 mm.

`build.py` tries all six ways of standing each cell on one of its own faces and picks the
best, reporting two different things per lattice in [PARTS.md](PARTS.md):

- **bridges** — horizontal struts with air under them. Every frame has four, the top face's
  edges, spanning at most 46 mm of Ø6 rod between hubs. PETG will sag a little across that
  and it will not fail; the bottom-level horizontals are excluded because they sag onto the
  bed where it does not show.
- **islands** — hubs that no strut reaches from below, which start in mid-air and *must*
  have support. Only **oF natrolite and cF diamond** have one: the face-centring hub on the
  top face, reached solely by struts in its own horizontal plane. Being face-centred, every
  face has a centring hub, so one is always on top and no orientation avoids it.

The other twelve have zero islands, so supports are optional. Worth noting that the
orientation search earned that: gypsum, cordierite and hemimorphite would each have had an
island resting on the obvious a-b face, and the search found faces that put their centring
hubs on the *sides*, where struts reach them from below.

Standing a cubic frame on a body diagonal instead would put all twelve edges at 35.3° and
need no bridging at all, but it balances a 104 mm tower on a single hub. A print that falls
over costs more than a support that has to be snipped, so the orientations are all face-down.

**Kit** — about 679 g for all 387 parts. Struts lie flat along the bed so their
layers run along the axis rather than across the joint. Nodes are oriented automatically to
lift their sleeves as far above horizontal as possible, and truncated at the bottom to give
a ~8.7 mm bed pad instead of a point contact.

Three node types want a little support: the body-centre nodes of cI, tI and oI carry eight
sleeves on the body diagonals, which come in antipodal pairs, so *no* orientation gets them
all above horizontal — the best available leaves the worst sleeve at −36° (cI), −35° (tI)
and −23° (oI). The face-centre and base-centre nodes of cF, oF, mS and oS sit near −2°,
which is a short horizontal cantilever. Everything else is at +21° or better and needs
nothing.

### Fits

| feature | dimension |
|---|---|
| hub | Ø14 mm |
| sleeve | Ø9.2 mm outside, Ø6.30 mm bore, 6 mm socket depth |
| strut | Ø6.00 mm rod, length = cell edge − 14 mm |

The Ø6.30 bore on a Ø6.00 rod is 0.30 mm of clearance, which is a friction fit in PETG on a
well-tuned printer. If yours runs tight, raise `BORE_R` in `build.py`; if the struts fall
out, lower it. Everything else follows from the parameters at the top of that file.

### Assembly (kit only — `assembled/` needs none)

Every node of a given type is interchangeable with every other of that type — the strut
spiders that hold the centring points use *all* the diagonals rather than the minimum
needed for rigidity, precisely so that this stays true. Filenames carry the socket count
and the quantity (`node_A_6way_x8.stl`), so match a node to a site by counting sockets.

Two things worth doing once assembled:

- **Three hP cells make the hexagonal prism.** The Bravais lattice is the primitive
  120° rhombic cell; the familiar hexagonal prism is three of them.
- **Stack a cI or cF cell** on a second copy and the centring point of one lands on the
  corner of the next, which is the whole content of the phrase "centring translation".

## Wigner-Seitz cells

The region closer to one lattice point than to any other. For the three cubic lattices,
these are three well-known solids, and their volumes are the primitive-cell volumes:

| lattice | Wigner-Seitz cell | volume | Brillouin zone |
|---|---|---|---|
| cP | cube | a³ | cube |
| cI | truncated octahedron | a³/2 | rhombic dodecahedron |
| cF | rhombic dodecahedron | a³/4 | truncated octahedron |

`build.py` asserts those volumes exactly, which is a sharp check on the vertex sets — the
first version of the truncated octahedron assigned a/4 to the lower-indexed free axis every
time, yielding 12 vertices instead of 24 and a completely different solid, and the volume
assertion is what caught it.

cI and cF **swap solids** between real and reciprocal space: the reciprocal lattice of bcc
is fcc and vice versa, so each one's Brillouin zone is the other's Wigner-Seitz cell. All
three space-fill, and at a common a = 60 mm they interlock.

Two links to the rest of this repo: the 12-fold coordination shell of cF is the
**cuboctahedron**, and diamond's cF lattice is two interpenetrating face-centred lattices
with **tetrahedral** coordination.

## Files

| file | role |
|---|---|
| `lattices.py` | the 14 lattices, cell-vector maths, hexagonal→rhombohedral conversion, node/strut graph, supercell merging, congruence dedup |
| `glyphs.py` | eleven-glyph stroke font (`a m o t h c` + `P S I F R`) for the embossed Pearson symbols |
| `solids.py` | mesh primitives: parallelepipeds, hulls, spheres, sleeved hubs, print-orientation search, overhang analysis |
| `build.py` | writes every STL plus `PARTS.md` |
| `verify.py` | reloads every STL in trimesh and checks each connected body independently |
| `preview.py` | the three preview figures |

### A note on watertightness

These parts are unions of *overlapping* watertight solids — a hub plus its sleeves, a cell
plus its embossed glyph prisms — because this environment has no mesh-boolean backend
(no manifold3d, fcl or rtree) to fuse them with. That is the same contract `mesh_kit.py`
uses elsewhere in the repo, and every slicer unions overlapping bodies on import.

So `trimesh.load(path).is_watertight` is `False` for a labelled part *by construction*, and
checking it would be checking the wrong thing. `verify.py` splits each file and requires
every connected body to be watertight, consistently wound and of positive volume; `build.py`
additionally checks that every solid is edge-manifold. Both pass for all 96 STLs.

One consequence worth knowing: solids that overlap are fine, but solids that are *exactly
coincident* are not — merging duplicate vertices on import welds them into a doubled shell
with zero volume, or leaves an edge used four times. Two places in this build had to be
written around it:

- the glyph code emits one join dot per distinct vertex, not per occurrence, because `P`'s
  stem meets its bowl at a shared point and a closed ring repeats its first point;
- fused struts stop 4.2 mm short of the hub centre, because in every centred lattice the
  run corner → centring point → opposite corner is *collinear*, so struts drawn all the way
  to the centre would put two end caps exactly on top of each other.

`verify.py` splits components with `scipy.sparse.csgraph` rather than `trimesh.split`, which
routes through `fill_holes` and needs networkx — absent here, so a genuinely broken body
crashed the checker with `ModuleNotFoundError` instead of being reported. That is how the
collinear-strut defect surfaced.
