#!/usr/bin/env python3
"""
Cymatics petri dish for a 12" driver -- OVERHANG version.

At Ø277 the dish is wider than any 12" cone mouth (a 305 mm frame with a
~287 mm bolt circle cannot carry a cone wider than about Ø270), so the dish can
no longer sit down inside the cone. Instead it bridges over the surround and the
skirt reaches UP to meet it:

  * DISH    -- flat floor + vertical rim, Ø277. Nothing on it touches the
               speaker. The floor is flat and of uniform thickness because
               shallow-water wavelength goes as sqrt(g*h); any dome or cone in
               the floor smears every pattern into a radial gradient.
               Prints floor-down, no supports.

  * SKIRT   -- three sections, one continuous 1 mm shell:
                 FLANGE, a horizontal annulus glued to the floor underside,
                         full Ø277 so it is flush with the dish and supports the
                         floor's overhang out to the rim,
                 COLLAR, a vertical wall dropping from the flange to the cone,
                         tall enough to hold the flange clear of the surround,
                 BAND,   the conical part whose outer face lies exactly on the
                         speaker cone, Ø264 down to Ø150. This is the only part
                         that touches the driver.
               Prints flange-down: every section steps inward going up, so no
               supports.

  * RING    -- Ø130 tube, floor underside to the cone. Shortens the floor's
               unsupported span and drives the middle of the floor directly.

  * BUMPS   -- 4 spheres on the collar, just above the contact circle, sitting a
               hair clear of the cone. They stop the assembly walking sideways
               without bearing on the paper; if they bore, they would lift the
               band off its seat and become the only contact.

GEOMETRY CHAIN
--------------
The cone model fixes one number, the slope K, and the rest follows.

    K            = depth / (r_mouth - r_dustcap)          [dz/dr]
    contact top  = (R_CON, -FLOOR_T - COLLAR_H)           [on the cone]
    band outer   : z(r) = that point extended at slope K
    COLLAR_H     : solved so the flange underside clears the surround crown

Measure your own driver and edit the CONE_* and SURROUND_* constants; the dish
is unaffected, only the skirt is cut to that angle.
"""

import numpy as np
import trimesh

# ---- speaker cone model (typical 12" PA driver -- MEASURE AND EDIT) -------
CONE_MOUTH_D   = 268.0     # cone dia where the paper meets the surround (mm)
CONE_DUSTCAP_D = 75.0      # dust cap dia (mm)
CONE_DEPTH     = 70.0      # vertical drop, cone mouth -> dust cap (mm)
SURROUND_H     = 9.0       # how high the surround crown stands above the cone
                           # mouth plane (mm) -- the dish has to bridge this
SURROUND_CLR   = 3.0       # air gap left under the flange over the surround (mm)

# ---- dish ----------------------------------------------------------------
OD        = 277.0          # outer diameter (mm)
RIM_T     = 1.0            # rim wall thickness (mm)
FLOOR_T   = 0.39375        # floor thickness (mm) -- 3 exact layers at 0.13125
RIM_H     = 11.25          # rim height above the floor's water side (mm)
FILLET    = 0.5            # floor->wall inner fillet radius (mm)
RIM_ROUND = 0.5            # rim top round-over radius (mm) = RIM_T/2, a true
                           # bullnose: any larger would undercut the wall
WATER     = 3.0            # nominal water depth used for the report (mm)

# ---- skirt ---------------------------------------------------------------
CONTACT_D   = 264.0        # top of the contact band (mm) -- must clear the cone
                           # mouth, and is what the dish overhangs
SKIRT_T     = 1.0          # shell thickness (mm)
SKIRT_D_BOT = 150.0        # small-end diameter (mm) -- stays clear of the cap
SKIRT_TOE   = 6.0          # radial length of the feathered toe (mm): the band
SKIRT_T_TOE = 0.6          # thins to this at its inner end so it does not
                           # terminate as a square edge digging into the paper
N_VENT      = 6            # vent notches across the flange's glue face
VENT_W      = 4.0          # notch width (mm)
VENT_D      = 0.5          # notch depth (mm) -- shallower than SKIRT_T so it
                           # channels air without slotting through the flange

# ---- locating bumps on the collar, 8 at 45 degrees -----------------------
N_BUMP     = 8
BUMP_D     = 3.0           # bump diameter (mm). Proud = BUMP_D - BUMP_EMBED, so
BUMP_EMBED = 0.3           # growing the ball is how a sphere stands further out
BUMP_CLR   = 0.2           # gap to the cone, measured on the cone's normal (mm)

# ---- floor support ring (a separate part) --------------------------------
SUPPORT_D   = 130.0        # ring diameter (mm) -- inboard of the band's toe, so
SUPPORT_T   = 1.0          # it lands straight on the cone
SUPPORT_GAP = 0.2          # clearance to whatever it lands on, for glue (mm)

SECTIONS = 480             # angular facets (~1.8 mm chords at Ø277)
NSKIRT   = 24              # samples along the band's slant
DENSITY  = 1.27e-3         # PETG, g/mm^3

OUT_DISH  = "cymatics_dish.stl"
OUT_SKIRT = "cymatics_skirt.stl"          # assembly position, for checking fit
OUT_SKIRT_PRINT = "cymatics_skirt_print.stl"   # flipped flange down, on the bed
OUT_SUPPORT = "cymatics_support_ring.stl"

R_OUT = OD / 2.0
R_IN  = R_OUT - RIM_T
R_CON = CONTACT_D / 2.0
R_MTH = CONE_MOUTH_D / 2.0
K     = CONE_DEPTH / (R_MTH - CONE_DUSTCAP_D / 2.0)
NRM   = np.array([-K, 1.0]) / np.hypot(K, 1.0)     # cone normal, up-and-inward
# collar tall enough that the flange underside clears the surround crown
COLLAR_H = SKIRT_T + (R_MTH - R_CON) * K + SURROUND_H + SURROUND_CLR

assert R_OUT > R_CON, "contact circle must be inboard of the dish rim"
assert CONTACT_D < CONE_MOUTH_D, "contact band would start outside the cone"
assert RIM_ROUND <= RIM_T / 2, "rim round-over would undercut the wall"
assert BUMP_D / 2 > BUMP_EMBED, "bump would not stand proud of the collar"
# --------------------------------------------------------------------------


def arc(cr, cz, rad, a0, a1, n=24):
    """Polyline arc in the (r, z) half-plane, angles in degrees."""
    a = np.radians(np.linspace(a0, a1, n))
    return np.column_stack([cr + rad * np.cos(a), cz + rad * np.sin(a)])


def cone_z(r):
    """Height of the speaker cone surface at radius r, in dish coordinates.

    Anchored on the top of the contact band. Material is clear of the cone when
    its z is ABOVE this (inside the funnel)."""
    return -FLOOR_T - COLLAR_H + (r - R_CON) * K


def solid(profile):
    """Revolve a closed (r, z) profile into a watertight solid."""
    if not np.allclose(profile[0], profile[-1]):
        profile = np.vstack([profile, profile[:1]])    # revolve wants it closed
    m = trimesh.creation.revolve(profile, sections=SECTIONS)
    m.merge_vertices()
    m = max(m.split(only_watertight=False), key=lambda c: c.area)
    m.remove_unreferenced_vertices()
    m.fix_normals()
    return m


def dish_profile():
    """Axis -> floor -> fillet -> rim wall -> rim round -> outside -> underside.

    The outer bottom edge is a free edge now: it hangs in air over the surround,
    so there is nothing for it to bed against."""
    return np.vstack([
        [[0.0, 0.0], [R_IN - FILLET, 0.0]],
        arc(R_IN - FILLET, FILLET, FILLET, -90.0, 0.0),          # inner fillet
        [[R_IN, RIM_H - RIM_ROUND]],
        arc(R_OUT - RIM_ROUND, RIM_H - RIM_ROUND, RIM_ROUND, 180.0, 0.0),
        [[R_OUT, -FLOOR_T], [0.0, -FLOOR_T]],
    ])


def skirt_profile():
    """Flange (flush to Ø277) -> collar (down to the cone) -> band (on the cone,
    Ø264 to Ø150), as one 1 mm shell."""
    r_bot = SKIRT_D_BOT / 2.0
    r_o = np.unique(np.append(np.linspace(R_CON, r_bot, NSKIRT),
                              r_bot + SKIRT_TOE))[::-1]
    z_o = cone_z(r_o)
    t   = np.interp(r_o, [r_bot, r_bot + SKIRT_TOE, R_CON],
                         [SKIRT_T_TOE, SKIRT_T, SKIRT_T])
    r_i = r_o + t * NRM[0]
    z_i = z_o + t * NRM[1]

    # the band's inner face runs out past the collar's inner wall; clip it there
    r_wall = R_CON - SKIRT_T
    keep = r_i <= r_wall
    j = int(np.argmax(keep))
    assert 0 < j < len(r_o), "band inner face never reaches the collar wall"
    f = (r_wall - r_i[j - 1]) / (r_i[j] - r_i[j - 1])
    z_cut = z_i[j - 1] + f * (z_i[j] - z_i[j - 1])
    inner = np.vstack([[[r_wall, z_cut]], np.column_stack([r_i, z_i])[keep]])

    return np.vstack([
        [[R_OUT, -FLOOR_T],                       # flange top, outer (flush)
         [R_OUT, -FLOOR_T - SKIRT_T],             # flange outer edge
         [R_CON, -FLOOR_T - SKIRT_T]],            # flange underside -> collar
        np.column_stack([r_o, z_o]),              # collar outer, then the band
        inner[::-1],                              # end cap, then back up inside
        [[r_wall, -FLOOR_T]],                     # collar inner -> flange top
    ])


def bump_centre():
    """(r, z) of a locating bump's centre on the collar.

    r is set by how far the sphere sinks into the collar; z is then solved so
    the perpendicular distance from the centre to the cone is exactly the
    sphere radius plus BUMP_CLR."""
    rb = BUMP_D / 2.0
    rc = R_CON - BUMP_EMBED + rb
    zc = cone_z(rc) + (rb + BUMP_CLR) * np.hypot(K, 1.0)
    return rc, zc


def bumps():
    """The four spheres, as one solid to union onto the skirt."""
    rc, zc = bump_centre()
    out = []
    for i in range(N_BUMP):
        s = trimesh.creation.icosphere(subdivisions=3, radius=BUMP_D / 2.0)
        th = 2 * np.pi * i / N_BUMP
        s.apply_translation([rc * np.cos(th), rc * np.sin(th), zc])
        out.append(s)
    return trimesh.util.concatenate(out)


def skirt_inner_z(r):
    """Height of the skirt's inner face, or of the bare cone inboard of the
    band's toe -- whichever the support ring actually lands on."""
    inside_toe = np.asarray(r) < SKIRT_D_BOT / 2.0
    return cone_z(r) + np.where(inside_toe, 0.0, SKIRT_T * np.hypot(K, 1.0))


def support_profile():
    """Plain tube: flat top glues to the floor underside, bottom cut on the cone
    angle to seat on the cone."""
    ro, ri = SUPPORT_D / 2.0, SUPPORT_D / 2.0 - SUPPORT_T
    return np.array([[ri, -FLOOR_T], [ro, -FLOOR_T],
                     [ro, float(skirt_inner_z(ro)) + SUPPORT_GAP],
                     [ri, float(skirt_inner_z(ri)) + SUPPORT_GAP]])


def vent_cutter():
    """N_VENT radial channels across the flange's glue face, as one solid."""
    boxes = []
    for i in range(N_VENT):
        b = trimesh.creation.box(extents=[R_OUT - R_CON + 8.0, VENT_W, 4 * VENT_D])
        b.apply_translation([(R_OUT + R_CON) / 2, 0.0, -FLOOR_T - VENT_D])
        b.apply_transform(trimesh.transformations.rotation_matrix(
            2 * np.pi * i / N_VENT, [0, 0, 1]))
        boxes.append(b)
    return trimesh.util.concatenate(boxes)


def report(name, mesh):
    e = mesh.extents
    print(f"  {name:<6} Ø{max(e[0], e[1]):.1f} x {e[2]:.1f} mm   "
          f"{mesh.volume * DENSITY:6.1f} g   "
          f"watertight={mesh.is_watertight}  winding={mesh.is_winding_consistent}  "
          f"faces={len(mesh.faces)}")
    return mesh.volume * DENSITY


def main():
    dish = solid(dish_profile())
    skirt = solid(skirt_profile()).difference(vent_cutter()).union(bumps())
    skirt.fix_normals()
    support = solid(support_profile())

    dish.export(OUT_DISH)
    skirt.export(OUT_SKIRT)

    # print-ready skirt: flipped so the flange's glue face is the bed face and
    # every section steps inward going up -- self-supporting, no supports needed
    flat = skirt.copy()
    flat.apply_transform(trimesh.transformations.rotation_matrix(np.pi, [1, 0, 0]))
    flat.apply_translation([0, 0, -flat.bounds[0][2]])
    flat.export(OUT_SKIRT_PRINT)
    support.export(OUT_SUPPORT)

    r_bot = SKIRT_D_BOT / 2.0
    water = np.pi * R_IN**2 * WATER * 1e-3          # g (== ml)
    rc, zc = bump_centre()

    print(f"cone model   : Ø{CONE_MOUTH_D:.0f} mouth, Ø{CONE_DUSTCAP_D:.0f} cap, "
          f"{CONE_DEPTH:.0f} deep -> slope {np.degrees(np.arctan(K)):.1f}° "
          f"from horizontal")
    print(f"overhang     : dish Ø{OD:.0f} bridges out {R_OUT - R_CON:.1f} mm past "
          f"the Ø{CONTACT_D:.0f} contact circle, over the surround")
    print(f"contact band : Ø{CONTACT_D:.0f} down to Ø{SKIRT_D_BOT:.0f}, "
          f"{np.hypot(R_CON - r_bot, cone_z(r_bot) - cone_z(R_CON)):.1f} mm slant")
    print(f"  toe        : feathers to {SKIRT_T_TOE:.1f} mm over the last "
          f"{SKIRT_TOE:.0f} mm of radius")
    print(f"collar       : {COLLAR_H:.2f} mm tall, holding the flange "
          f"{SURROUND_CLR:.1f} mm clear of a {SURROUND_H:.0f} mm surround crown")
    print(f"flange       : Ø{OD:.0f} glue face {R_OUT - (R_CON - SKIRT_T):.1f} mm "
          f"wide, {N_VENT} x {VENT_W:.0f} mm vent notches")
    print(f"bumps        : {N_BUMP} × Ø{BUMP_D:.0f} on the collar at r={rc:.2f} "
          f"z={zc:.2f}, {rc + BUMP_D / 2 - R_CON:.2f} mm proud, "
          f"{BUMP_CLR:.1f} mm off the cone")
    print("parts:")
    m = report("dish", dish) + report("skirt", skirt)
    ms = report("ring", support)
    print(f"assembly     : {m:.0f} g plastic (+{ms:.0f} g with the ring) "
          f"+ {water:.0f} g water ({WATER:.0f} mm deep) = {m + water:.0f} g")
    print("levels (dish floor water-side = 0):")
    print(f"  rim top          z = +{RIM_H:.1f}")
    print(f"  flange underside z = {-FLOOR_T - SKIRT_T:+.1f}")
    print(f"  surround crown   z = {cone_z(R_MTH) + SURROUND_H:+.1f}  "
          f"(assumed {SURROUND_H:.0f} mm above the Ø{CONE_MOUTH_D:.0f} mouth)")
    print(f"  contact top      z = {cone_z(R_CON):+.1f}  (Ø{CONTACT_D:.0f})")
    print(f"  band toe         z = {cone_z(r_bot):+.1f}  (Ø{SKIRT_D_BOT:.0f})")
    print(f"  dust cap edge    z = {cone_z(CONE_DUSTCAP_D / 2.0):+.1f}  "
          f"(Ø{CONE_DUSTCAP_D:.0f})")
    print(f"exported {OUT_DISH}, {OUT_SKIRT}, {OUT_SKIRT_PRINT}, {OUT_SUPPORT}")


if __name__ == "__main__":
    main()
