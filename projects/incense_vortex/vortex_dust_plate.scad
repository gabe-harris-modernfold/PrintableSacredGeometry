/* ============================================================
   VORTEX DUST PLATE — wind-driven dust-devil arena
   ============================================================
   A flat, open-topped swirl arena that lifts a tracer powder
   (Borozin / zinc stearate, or any micronized powder) into a
   visible free-standing vortex using only ambient wind.

   How it works: twelve guide fins form a converging crown that
   turns wind from ANY azimuth into same-handed tangential jets
   (~4:1 area contraction through ~9 mm throats at 58 deg from
   radial). Inside the arena the flow spirals inward, spinning up
   as angular momentum is conserved; the core's pressure minimum
   sits over a shallow saucer dimple at dead center where the
   powder pools. The dished floor returns centrifuged powder to
   the middle for re-entrainment; a removable roof ring keeps the
   wind in the channels and leaves a 136 mm sky opening for the
   devil to rise through.

   Expectations: this makes dancing, intermittent dust whirls in
   a steady 2-5 m/s breeze, not a permanent standing column - a
   free vortex has no tube to confine it. A dark-colored print in
   sunshine adds thermal updraft (real dust-devil physics) and
   noticeably improves it.

   Swirl is CLOCKWISE viewed from above, matching the family.

   Parts (set `part`):
     0 = assembly preview
     1 = arena plate  (print upright, no supports)
     2 = roof ring    (print flat, no supports)
     9 = cross-section (preview/PNG only)
   ============================================================ */

part = 0;

/* [Overall] */
plate_or   = 100;    // 200 mm plate
skirt_th   = 6;      // base slab
total_note = 79;     // fins 76 + roof 3 (<= 100 mm budget)

/* [Arena] */
arena_r    = 60;     // swirl arena radius
ped_r      = 63;     // arena pedestal (fins embed into it)
ped_h      = 12;     // arena floor plane
dish_in_r  = 22;     // dished floor: rim z12 -> z8 at this radius
dimple_r   = 20;     // powder saucer half-width (40 mm dish)
dimple_d   = 3;      // saucer depth

/* [Fin crown] */
n_fins     = 12;
fin_angle  = 58;     // jet entry angle from radial
fin_th     = 2.4;
fin_h      = 70;
fin_clip_r = 92;     // fin outer trim
fin_x0     = -77;    // fin outer end (pre-clip)
pin_d      = 5;      // roof locating pins on fin tops
pin_h      = 3;
pin_r      = 81;     // radius of the pin circle (on the fin chords)

/* [Roof ring] */
roof_ir    = 68;     // sky opening r (136 mm)
roof_th    = 3;
rib_r      = 86;     // stiffening rib inner radius

/* [Anchoring] */
stake_r    = 96;

/* derived */
chord_p    = arena_r*sin(fin_angle);          // 50.88
fin_x1     = -sqrt(arena_r*arena_r - chord_p*chord_p); // tip at arena rim, embeds pedestal ring
pin_x      = -sqrt(pin_r*pin_r - chord_p*chord_p);     // pin on the chord
dimple_R   = (dimple_r*dimple_r + dimple_d*dimple_d) / (2*dimple_d);

$fa = 3;
$fs = 0.5;

/* ---------------- fin crown ------------------ */
module fins() {
    intersection() {
        union() {
            for (i = [0 : n_fins-1])
                rotate([0, 0, i*360/n_fins])
                    translate([fin_x0, chord_p, skirt_th])
                        cube([fin_x1 - fin_x0 + 0.01, fin_th, fin_h]);
        }
        cylinder(r=fin_clip_r, h=skirt_th + fin_h + 1);
    }
}

module roof_pins() {
    for (i = [0 : n_fins-1])
        rotate([0, 0, i*360/n_fins])
            translate([pin_x, chord_p + fin_th/2, skirt_th + fin_h])
                cylinder(d=pin_d, h=pin_h);
}

/* ---------------- arena plate ------------------ */
module plate() {
    difference() {
        union() {
            cylinder(r=plate_or, h=skirt_th);          // base slab
            cylinder(r=ped_r, h=ped_h);                // arena pedestal
            fins();
            roof_pins();
        }
        // dished floor: rim z12 descending to z8 at r22, flat inside
        rotate_extrude()
            polygon([[0, ped_h - 4], [dish_in_r, ped_h - 4],
                     [arena_r, ped_h + 0.01], [0, ped_h + 0.01]]);
        // powder saucer dimple (40 mm dia x 3 deep, spherical)
        translate([0, 0, (ped_h - 4) + dimple_R - dimple_d])
            sphere(r=dimple_R);
        // 3 countersunk stake / screw holes
        for (a = [5, 125, 245]) rotate([0, 0, a]) {
            translate([stake_r, 0, -1]) cylinder(d=4.5, h=skirt_th + 2);
            translate([stake_r, 0, skirt_th - 1.6]) cylinder(d=9, h=2);
        }
    }
}

/* ---------------- roof ring ------------------ */
module roof() {
    difference() {
        union() {
            difference() {
                cylinder(r=plate_or, h=roof_th);
                translate([0, 0, -1]) cylinder(r=roof_ir, h=roof_th + 2);
            }
            // stiffening rib
            rotate_extrude()
                polygon([[rib_r, roof_th], [rib_r + 3, roof_th],
                         [rib_r + 3, roof_th + 2], [rib_r, roof_th + 2]]);
        }
        // pin sockets: through-holes so the ring seats on the fin tops
        for (i = [0 : n_fins-1])
            rotate([0, 0, i*360/n_fins])
                translate([pin_x, chord_p + fin_th/2, -1])
                    cylinder(d=pin_d + 0.8, h=roof_th + 4);
    }
}

/* ---------------- assembly ------------------ */
module assembly() {
    plate();
    translate([0, 0, skirt_th + fin_h]) roof();
}

if (part == 0) assembly();
if (part == 1) plate();
if (part == 2) roof();
if (part == 9) difference() {
    assembly();
    translate([-plate_or-5, -plate_or-5, -1])
        cube([2*plate_or+10, plate_or+5, 120]);
}
