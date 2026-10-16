# Negative Space — a 1 m ribbon around an invisible sphere

A single spiral ribbon, 1 m in diameter, wound around a sphere that is not there.
The ribbon's inner face lies exactly on the surface of that sphere, so the viewer's
visual system completes it: the sphere exists only in the viewer's head. The data
lives in the *holes* in the ribbon, not in the plastic. Light from the center of
the absent sphere passes through both absences and paints the year onto the walls
of the room.

Status: **design discussion, no generator or STLs yet.** Everything below is a
proposal with concrete numbers so it can be turned into a script.

## Perceptual basis

- The seed image was the "Simplicity Metric, based on symmetry" figure: a 3D object
  (a ribbon wound on a cylinder) projected onto two planes, a messy gray silhouette
  on one and two clean mirror-image blobs on the other. One 2D image is consistent
  with infinitely many 3D shapes (the inverse projection problem). The visual
  system picks the interpretation that scores best on a simplicity metric.
- Names for the principle: the Gestalt **law of Prägnanz** (law of simplicity /
  good figure), its formal version the **minimum principle** (Hochberg & McAlister,
  Leeuwenberg's structural information theory), and the specific symmetry flavor,
  Zygmunt Pizlo's **symmetry prior** for 3D shape recovery. It is one side of the
  classic simplicity-vs-likelihood (Gestalt vs Helmholtz) debate.
- Why the sphere works: every turn of the ribbon lies on one common surface.
  Gestalt closure / amodal completion fills it in, and the sphere is the maximally
  symmetric candidate, so the symmetry prior selects it.
- Two layers of absence: the sphere (implied by the ribbon) and the data (holes in
  the ribbon). The sculpture is itself a demonstration of the principle: from the
  side it is a tangle, from one projection it is a legible calendar.
- Optional metaphor: the sphere is Earth. Indiana sits at 39.8° N, so the
  temperature data describes a point on an object that is not there.

## Geometry (1 m diameter sphere)

| Parameter | Value | Note |
|---|---|---|
| Sphere radius | 500 mm | ribbon's inner face sits exactly on it |
| Path | spherical spiral, latitude linear in day-of-year, pole to pole | Jan 1 at the top pole, Dec 31 at the bottom: the year descends |
| Turns | 12.37 (one per synodic month) | every full moon lands on the same meridian: a visible line of discs down the sphere; new moons on the opposite meridian |
| Meridional pitch | ≈ 127 mm per turn | 60–80 mm band leaves ~50 mm of "sphere" visible between turns |
| Total path length | ≈ 24 m | average 66 mm per day; ~100 mm at the equator, ~14 mm near the poles |
| Taper | width ∝ cos(latitude) near the poles | ribbon narrows to a point at each pole like an apple peel; also avoids impossible tight turns |
| Thickness | 3 mm PETG | ≈ 6 kg total at 80 mm width, ≈ 3.5 kg at 50 mm |

Spiral in sphere coordinates, t = (day − 1) / 365:

```
lat(t) = -90 deg + 180 deg * t        (top pole = Jan 1, bottom pole = Dec 31)
lon(t) = 360 deg * N * t              (N = 12.37 turns, one per lunation)
```

Turn-count feel: fewer than ~8 turns reads as separate rings, more than ~16 reads
as a solid ball. About 12 keeps the sphere implied but absent.

## The three encodings (all negative space)

1. **Days of the year** — one cell per day along the path; a 2 mm notch on the
   outer edge each day, a 10 mm full-width slit at each month boundary.
2. **Phases of the moon** — one aperture per day, cut to the *illuminated* shape
   of the moon (nothing at new, crescent, gibbous, full disc). Aperture diameter
   about 0.6 x the local band width. With one turn per lunation the full moons
   stack on one meridian, so the phase data becomes visible *structure*, not a
   code to read.
3. **Average Indiana temperature** — band width on top of the polar taper:
   40 mm at the January normal (about 29 °F) up to 90 mm in July (about 76 °F).
   The ribbon's silhouette is the climate curve wrapped around the column.

Indianapolis 1991–2020 monthly normals, approximate (°F): Jan 29, Feb 33, Mar 43,
Apr 54, May 64, Jun 73, Jul 76, Aug 75, Sep 68, Oct 56, Nov 44, Dec 34. Use the
real daily normals, not these.

### Data sources

- Temperature: NOAA NCEI 1991–2020 **daily** normals, Indianapolis International
  Airport, station `USW00093819`.
- Moon: Meeus synodic formula (29.530589 d from a reference new moon, e.g.
  2000-01-06 18:14 UTC) is within a day; Skyfield for exact. Phases shift ~11
  days per year, so a specific year is required. **2027 assumed**; stamp it on
  the base.

## Readout modes

- **Shadow lamp (recommended):** an LED filament at the sphere's center (inside
  the axis tube) throws every aperture radially onto the room. The room is the
  projection plane from the seed image: tick marks for days, a band of moon
  shapes at the height of each turn, the climate curve as the shadow envelope.
- **Naked eye:** the full-moon meridian and the width swell in summer are
  readable without a key; day notches need a close look.
- **Backlit:** an opaque sphere-shaped diffuser inside would make the sphere
  present, which defeats the point. Avoid.

## Structure, since nothing is holding it up

The ribbon surface is a patch of sphere, so it is **not developable**. It cannot
be bent from flat stock; it must be printed or thermoformed.

- **Recommended skeleton:** a 20 mm axis tube through the poles that doubles as
  the stand post, plus one 2 mm stainless wire or clear acrylic spoke per half
  turn (about 25) from the tube to hidden bosses on the ribbon's inner face. The
  spokes are radii pointing to a center that is not there, which fits the idea.
  The LED filament for projection lives in the tube.
- **Stronger story, more work:** assemble the segments on a 1 m foam or
  inflatable form, bond them, then remove the form. The sphere literally existed
  and was taken away.
- **Avoid a clear acrylic sphere.** It makes the sphere present, and 1 m domes
  cost real money.
- Weight at 80 mm width is ~6 kg; the skeleton has to be real, not decorative.

## Print plan (320 x 320 x 320 mm bed, PETG)

- Segments of up to ~300 mm along the path, ~80 pieces total; dovetail joints
  with a hidden spline on the inner face.
- Each piece prints **standing on its edge**. Across an 80 mm width the sagitta
  is only 80^2 / (8 x 500) = 1.6 mm, so the lean is negligible; along a 300 mm
  length the sagitta is ~22 mm, so the footprint is roughly 300 x 25 mm.
- Near the poles the turn radius drops (about 130 mm at ±75°), so polar pieces
  are shorter arcs; the taper keeps them narrow anyway.
- Opaque or matte filament so only the holes pass light.
- Verify every part with `trimesh.load(path).extents` and `is_watertight`
  before calling it done (repo rule).

## Modeling pipeline (proposed)

- Model directly in sphere coordinates. For each day t: compute (lat, lon),
  offset across the surface for the local width, extrude radially outward 3 mm,
  project the hole outlines radially and subtract. Python with NumPy + manifold3d
  (or build123d, both already installed here). OpenSCAD works but 365 booleans is
  slow.
- Alternative: Blender, flat ribbon with holes -> curve modifier along the spiral
  -> shrinkwrap to the sphere -> solidify 3 mm -> split with boolean cutters.
  Python is more repeatable.
- Emit per-segment STLs plus a flat SVG proof sheet of the unrolled year (the
  unrolled strip is itself a readable 1-D timeline, even though the real part is
  not developable).

## Earlier iteration (superseded): band helix on a cylinder

Before the sphere, the design was a band helix on a 1 m tall cylinder: 12 turns
(one per month), 83 mm pitch, 100 mm radius, 70 mm band, ~21 mm per day cell,
printed as 12 one-turn parts on a 25 mm core tube. Its advantage was that a band
helix *is* developable (design flat, wrap exactly). It was dropped because the
cylinder is present, not implied; the sphere makes the absence the subject.

## Open decisions

- One turn per lunation (moon alignment becomes structure, months scatter) versus
  one turn per calendar month (months become structure, moons scatter).
  Recommendation: lunation; the stacked full-moon meridian is the strongest
  visual.
- Year for the moon data (2027 assumed).
- Band width at the summer peak: 90 mm (stronger silhouette, 6 kg) or 60 mm
  (lighter, more sphere visible).
- Skeleton: axis-and-spokes vs. lost-form build.

## Next steps

1. Pull NCEI daily normals and compute 2027 moon illumination per day.
2. Write `negative_space_sphere.py` generating the full ribbon, then the
   segmented STLs and the SVG proof sheet.
3. Print one equatorial and one near-polar segment to check the on-edge print and
   the dovetail fit before committing to ~80 pieces.
