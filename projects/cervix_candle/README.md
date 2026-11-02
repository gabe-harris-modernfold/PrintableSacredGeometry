# Uterus candle holder

A uterus candle holder modelled on the user's reference photo (`reference_photo.png`) of a white
sculpture after a coronal-section diagram (`reference_diagram.png`). The uterine tubes are cut
off at two flat tabs. The front is flat, like a section through the organ, with softly rounded
edges. The back is domed. The cavity is hollowed into the face like the inside of a vessel cut in
half. A wax birthday candle is held by light pressure where the walls close round it at the
internal orifice. Its flame burns in the hollow under the closed fundus band, and heat leaves
along the tube-lumen grooves that run out through the tabs. The sculpture stands in a socket in
the middle of a round dish that catches the wax.

| file | what |
|---|---|
| `photo_trace.py` | 2D layout: silhouette thresholded from the photo. The hollow, canal, bottom cup and lumen grooves are traced in photo pixels and mirrored. The canal and lips are taken at the photo's brightness edges, row by row. Run alone → `photo_trace_check.png` overlay |
| `cervix_candle.py` | sculpture → `cervix_candle.stl` (print pose, face down) and `cervix_candle_upright.stl` |
| `dish.py` | dish with the footprint socket → `dish.stl` (run after the sculpture) |
| `preview_cervix_candle.py` | `cervix_candle_views.png`, `cervix_candle_vs_photo.png` (photo beside a candle-lit render), `cervix_candle_canal_vs_photo.png` (the lips of both, enlarged) and `cervix_candle_closeup.png` (the canal) |

## Numbers

- Sculpture 83.7 wide × 100.2 tall × 30.3 deep mm, 88 cm³, watertight. The face edges have a
  2.5 mm roundover, and every cut's rim a 3 mm bullnose (the photo's soft edges).
- Back dome: the outline inflated as the union of its inscribed spheres, at 0.85 of a round
  section. The hollow is inflated the same way, fully round. The thinnest back wall behind the
  hollow is 10.3 mm.
- Canal: the slit between the lips follows the photo's S-curve. It is 4.2 mm wide where the
  lips close just below the hollow (z = 23–26). It swells to 6.4 mm (z = 15–17), pinches to
  3.4 mm at the lips' knees (z = 8) and opens round the lips' rounded ends into the cup.
- Candle (0.23 in = 5.84 mm): it stands down in the canal, recessed 6 mm behind the face, as in
  the photo. Behind the lips runs a smooth hourglass channel: 6.2 mm wide down the canal (the
  candle slides), narrowing to a 5.5 mm ring at z = 24.2, behind the canal's upper pinch.
  Pressing the candle down through the ring gives the light-pressure grip. The channel ends in
  a rounded floor, and the candle's foot rests on it at z = 9.7. **Pre-cut the candle to about
  28.5 mm** so the wick sits where the photo's does (z = 38.2). If the ring prints too tight or
  too loose, change `RING_R` (radius, mm) in `cervix_candle.py`.
- The candle can't come out forward. The lips are narrower than it at the ring (4.2 mm at
  z = 23–27) and in front of its foot (3.8–5.4 mm at z = 9–12). Through the swell between
  (z = 13–19, up to 6.4 mm) it is in view, as the photo's canal is.
- No ledges: every cut (hollow, canal, channel, cup, grooves) is blended into its neighbours
  with fillets and rounded into the face. The only flat faces are the front face and the base,
  and the base sits in the dish socket.
- Flame: the hollow's ceiling over the wick is at z = 78.8. A generous 25 mm flame leaves 15.6 mm
  to the fundus band, and a typical 15–20 mm birthday flame leaves 21–26 mm. Use the Milliput
  coat on the hollow's upper part. Keep the canal, the candle channel and the 3 mm lumen grooves
  free of Milliput (mask them), or the candle won't fit.
- Dish 155 × 155 × 18 mm. A wax well 22 mm wide runs round a 10 mm boss. The boss has a 5 mm-deep
  socket (0.4 mm clearance) shaped like the sculpture's base, so the sculpture seats one way
  only.

## Printing (PETG)

- Sculpture: flat face down on the bed (2,660 mm² of contact). The domed back prints on top
  without support. The hollow is a bowl over the bed: enable supports, blocked everywhere except
  inside the hollow (tree supports are easiest to pull out of it). The Milliput coat hides the
  scars. The canal, candle channel and lumen grooves bridge without support.
- Dish: flat on the bed. The outer wall leans out at about 30 degrees, and the boss and inner
  wall slope at 45 degrees or less.
