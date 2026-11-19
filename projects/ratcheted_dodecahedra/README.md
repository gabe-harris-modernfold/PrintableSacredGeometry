# ratcheted_dodecahedra

Dan Winter's "DNA helix is ratcheted dodecahedra" ([credit](#credit)) as a printable strut lattice: ten dodecahedra,
each turned 36° and lifted φ from the last, making one full turn of a double helix. The two helix
strands are printed in heavier struts.

```bash
python ratcheted_dodecahedra.py
```

| Output | |
| --- | --- |
| `ratcheted_dodecahedra.stl` | one turn, 129 × 123 × 308 mm, no base, one watertight body |
| `ratcheted_dodecahedra_views.png` | side, oblique and axial views, the two strands in red and blue |

## The geometry, and why the graphic's numbers are exact

Stand a dodecahedron on a face. Its 20 corners lie in four rings of five: a small face ring at the
top, a wide zig-zag girdle round the waist, and the mirror-image pair below. Choose the edge as
2φ² sin 36° and the graphic's labels are all exact:

| Label in the graphic | Here | Value |
| --- | --- | --- |
| radius of turn | girdle radius | φ³ |
| horizontal increment | girdle side seen down the axis (2φ³ sin 18°) | φ² |
| vertical increment | rise of each girdle zig-zag edge | φ |
| | face-ring radius / height | φ² / φ² + φ/2 |

**The ratchet.** Turn a copy 36° and lift it φ. Its lower girdle lands exactly on the original's
upper girdle, so the two share five corners and lock together, one 36° tooth per step, ten steps to
the turn. A girdle edge that climbs +36° continues straight into the next copy's, so the chain of
them is a helix: radius φ³, pitch 10φ.

**Why it is a double helix.** Five right-handed strands thread the column, 72° apart. Five
left-handed ones mirror them, because turning a dodecahedron 36° is the same as flipping it end for
end. The print highlights two right-handed strands 144° apart. That leaves a narrow and a wide
groove (144° and 216°), the way DNA has a minor and a major groove.

`verify()` asserts each of these numbers before anything is built.

**Size.** Seen down the axis, the girdle decagon has radius 4R = 60 mm, which is the outer circle of
the [matched set](../from_one_circle) (`lib/sacred_grid.py`). One step is φ·14.16 = 22.9 mm, the edge
is 43.6 mm, and a turn is 229 mm.

## Printing (PETG, no slicer supports)

- **Diamond struts.** Every strut has a diamond cross-section with one diagonal vertical. The
  shallow girdle edges (31.7°) have 53° faces underneath, and the level pentagon edges are 45°
  ridges that bridge 43.6 mm from corner to corner. Under 30 mm² of the 2,600 cm² surface hangs
  steeper than 45°.
- **Posts.** Every face-ring corner of every copy lies on one of ten vertical lines (radius φ²,
  one every 36°), with copies two steps apart on the same line. Posts along those lines carry the
  face pentagons that would otherwise float, so nothing starts in mid-air. Together they look like
  the column's inner decagonal cage. Set `POSTS = False` to get the bare dodecahedra, then print
  with tree supports.
- Struts are 6 mm (diagonal), strands 9 mm and posts 4.5 mm. Joints are the convex hull of the
  strut ends that meet there.
- **No base.** The column stands on its bottom pentagon. The ridges of those struts are trimmed
  1.2 mm flat (`BED_FLAT`), and five post feet come down beside them, giving 576 mm² on the bed. A
  308 mm column on that footprint needs a slicer brim (8–10 mm), which peels off afterwards.
- 0.2 mm layers, 3 walls. The struts are thin, so they print as almost solid walls (about 260 cm³).

## Credit

The concept and the graphic "DNA Helix is Wratcheted Dodecahedra" are by **Dan Winter**:
[goldenmean.info/dnamanifesto](https://www.goldenmean.info/dnamanifesto). This project builds his
picture as an exact, printable solid. The derivation of his φ, φ² and φ³ labels above is ours.
