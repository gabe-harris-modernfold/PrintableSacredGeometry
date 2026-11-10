# Raking-light book cover — brief

Project start plan, 2026-11-10. Nothing is built or printed yet. Every claim about how the
surface will look is reasoning, not measurement, until the coupon in "First build step" is printed.

## What it is

A printed plate bonded permanently over the front board of a hardcover book. In room light it
reads as a densely detailed ornamental cover. Under a flashlight held at a low angle it shows
three hidden things: a word, a symbol and an image.

The hiding principle: **encode in surface orientation, not in height or colour.** Diffuse light
averages over orientation, so the surface looks uniform. A directional raking beam does not, so
regions whose top-layer extrusion lines run in different directions separate into light and dark.

## Confirmed requirements

- Hardcover book; the plate is bonded permanently over the original cover.
- FDM, on the repo printer: 320 × 320 × 320 mm bed, PETG (see `CLAUDE.md`).
- Three hidden contents: a word, a symbol, an image.
- To normal observation the plate looks very detailed.

## Settled

### Visible design

- **A dense tessellation** (scales, star-and-polygon lattice, knotwork cells; style still open).
  Every cell needs a border anyway, so the seams between hatch directions sit under the visible
  linework, and varied hatching inside cells reads as ornament.
- **Engraved, never embossed.** A groove keeps its shadow inside itself. Raised ornament throws
  long shadows across the reveal (a 0.4 mm step casts about 4.6 mm at 5°) and outlines itself in
  room light.
- **Second colour from a layer swap.** Lower layers in a bright colour, top two layers in the
  field colour, grooves cut two layers deep so they show as bright lines. One manual filament
  change, no multi-material unit.

### Hidden layer

Each cell's top-layer hatch direction is the data. A cell lit across its lines goes dark and
textured (response about |sin| of the angle between light and lines); lit along them it stays
bright. Direction carries exactly two independent channels, much like polarisation:

| Cell class | Hatch directions | Revealed by light from | Under the other key |
|---|---|---|---|
| A | 0° or 90° (parallel / perpendicular to the edges) | a side edge; the adjacent edge shows the negative | uniform mid-tone (0.71 for both states) |
| B | 45° or 135° (the two diagonals) | a corner, along the diagonal | uniform mid-tone |

- The two classes are interleaved evenly across the cover (checkerboard where the tessellation
  allows, otherwise any balanced scatter), so each key shows its own content only.
- **Edge light:** the word in a border band and the symbol in a centre medallion, together.
- **Corner light:** the image across the whole field.
- Proposed 2026-11-10, not yet confirmed by Gabe: two keys for three contents, with word and
  symbol sharing one.

### Plate

- About 1.0 mm, fully solid: 5 layers at 0.2 mm, layers 1–3 bright, 4–5 field colour.
- Grooves 0.4 mm deep, 0.6–0.8 mm wide.
- Top fill: monotonic lines, slow, one perimeter per cell on the top layer.
- Stop the plate 8–10 mm short of the spine groove so the hinge still works, and inset it
  2–3 mm from the other three edges.
- Bond with a sheet transfer adhesive (3M 467 type), not liquid glue.

## Rejected — do not propose again

1. **Plain raised or recessed message.** The slicer draws perimeters round it; visible in any light.
2. **On-edge sawtooth plaque.** A thin cover plate has to print flat.
3. **Raised title, raised border, bosses.** They shadow the raking beam.
4. **A third independent light key from hatch direction.** There are only two. The one other
   mechanism considered, ironed cells that glint, shows in room light.
5. **Wipe-on paint or gilding wax to pick out the grooves.** It would fill the gaps between
   beads, which are the hidden channel.
6. **Greyscale by intermediate hatch angles.** It breaks the independence of the two keys.
   Greyscale comes from dithering one-bit cells.

## Limits to design within

- **Image resolution.** Each key owns half the cells. At 4–6 mm cells on a roughly 310 mm plate
  that is about a 36 × 36 to 55 × 55 one-bit picture. A bold silhouette or stencil-style portrait
  works; a photo does not.
- **Lettering.** A stroke needs about 10 beads (4–5 mm) of width before its direction reads, so
  letters are 25 mm tall or more.
- **Light distance.** A light held at the corner of the cover sends rays across it at anything
  from 0° to 90°, which mixes the two keys. Design for a light a metre or more away at table
  height, in a dark room. Alternative: compute the hatch for one marked light position and put a
  sun motif there in the visible art as the clue. With two near positions the keys only stay
  independent where their rays cross at about 45°.
- **Shadow length** = height ÷ tan(light elevation): 0.2 mm gives 2.3 mm at 5° and 5.7 mm at 2°.

## Risks

- **PETG.** It is glossier than matte PLA, so the hatch direction may show as sheen in room
  light, and it strings more, which matters across several thousand small top-layer islands.
  A matte or carbon-filled PETG may fix the sheen (carbon-filled needs a hardened nozzle).
  Untested; the coupon decides.
- **Field colour.** Expected best: opaque matte mid-tone. Expected worst: black (only sheen
  shows) and translucent (fine shadows wash out).
- **Clutter.** The engraved lines stay visible during the reveal, like a screen in front of the
  hidden content.
- **Uniformity.** Extrusion inconsistency, travel scars and blobs on the top layer are revealed
  by the same light as the message.
- **Wear.** A cover is handled constantly; burnishing and skin oil may expose or erase the
  message over time. Effect of a matte clear coat on the ridges is unknown.
- **Slicer.** The whole approach needs the slicer to honour a per-modifier fill angle on the
  top solid layer, with the same layer parity for all four modifiers.

## Open

1. **Cover dimensions.** "320 × 320" was given, and that is the bed. A plate up to about
   310 mm a side prints in one piece; a larger cover is tiled along cell borders.
2. Two light keys accepted, or is a third wanted badly enough to take a weaker mechanism?
3. Tessellation style, and whether it ties to this repo's `lib/sacred_grid.py` (R = 15 mm).
4. The actual word, symbol and image.
5. Distant light, or one marked near position with a clue in the art?
6. Filament: which PETG, which two colours.

## First build step

A 60 mm coupon from the same generator that will make the plate. It answers, in order:

1. Does the slicer honour the four fill angles on the top layer?
2. Is the hatch direction hidden in room light with the chosen PETG?
3. Which cell size reads under raking light: 4, 5 or 6 mm?
4. How much do the grooves clutter the reveal?
5. Bead texture strength: standard 0.42 mm lines, versus 0.6 mm lines at slightly reduced flow.

Pick the weakest texture that still reads.

## Planned files

Per repo convention, scripts run from this directory and write bare filenames.

- `cover.py` — three bitmaps in (word, symbol, image) plus tessellation parameters; writes
  `cover_plate.stl` (engraved plate) and `mod_000.stl`, `mod_045.stl`, `mod_090.stl`,
  `mod_135.stl` (one modifier mesh per hatch angle, each the union of that angle's cells through
  the top layer).
- `coupon.py` — the 60 mm test piece from the same code.
- `preview_cover.py` — simulated reveal for a given light azimuth, using the |sin| response, to
  check the content reads before printing.
- `verify.py` — `trimesh` extents within the bed, `is_watertight`, and the two-key independence
  check (each class is uniform under the other key).
- `README.md` — print settings and the slicer steps, once the coupon has settled them.

Bitmaps can be read with Matplotlib, which is already installed; read `docs/PYTHON_TOOLING.md`
before adding any dependency.
