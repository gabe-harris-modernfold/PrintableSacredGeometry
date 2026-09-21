/* ============================================================
   WIND-POWERED INCENSE TORNADO — omnidirectional wind harvester
   ============================================================
   Variant of the passive incense tornado optimised to CATCH THE
   WIND and convert it into swirl.

   INTAKE PINWHEEL (base): 12 tall guide fins radiate from the
   drum, each extending one tangential slot's guide wall. From any
   wind azimuth, the windward fins form converging funnels (~4:1
   area contraction) that accelerate the breeze and inject it as
   same-handed tangential jets. Ram pressure of a 2-3 m/s breeze
   (~5 Pa) is ~8x the stack-effect draft of the 394 mm tube.

   The tube exit runs fully OPEN. A stacked-disc venturi cowl was
   tried on the crown and removed after field testing: it did not
   enhance the helix. (Consequence: no rain / downdraft hat — this
   is now a fair-weather instrument.)

   In dead calm it still works like the indoor version (heat-only,
   gentler). Swirl is CLOCKWISE viewed from above.

   Parts (set `part`):
     0 = assembly preview (tube ghosted)
     1 = wind base: slotted drum + pinwheel fins + staked skirt
   Tube: 100 mm OD acrylic, 2 mm wall, 15.5" (393.7 mm) long.
   ============================================================ */

part = 0;

/* [Tube fit] */
tube_od    = 100;
tube_wall  = 2;
tube_len   = 393.7;
fit_gap    = 0.4;

/* [Swirl geometry] */
chamber_r  = 45.5;
vane_angle = 68;
n_slots    = 12;
slot_w     = 10;
slot_h     = 47;     // taller than indoor version: more catch area
drum_or    = 78;

/* [Wind intake] */
fin_len    = 50;     // guide fin length beyond the drum
fin_th     = 3.2;
fin_clip_r = 119;    // fins trimmed to this radius
skirt_r    = 121;    // stability skirt
stake_r    = 112;    // bolt-down / tent-stake holes

/* [Base build] */
floor_th   = 3;
ceil_th    = 3;
collar_h   = 10;
ped_or     = 17;
ped_wall   = 2.4;
recess_r   = 15;
recess_d   = 2;

/* derived */
tube_or   = tube_od/2;
tube_ir   = tube_or - tube_wall;
groove_or = tube_or + fit_gap;           // 50.4
groove_ir = tube_ir - fit_gap;           // 47.6
chord_p   = chamber_r*sin(vane_angle);   // 42.19 -> 68 deg entry
z_slot0   = floor_th;                    // 3
z_slot1   = floor_th + slot_h;           // 50
drum_h    = z_slot1 + ceil_th;           // 53 (groove floor)
base_h    = drum_h + collar_h;           // 63

$fa = 3;
$fs = 0.5;

/* ---------------- swirl slot (one chord channel) ------------- */
module slot_cut() {
    // drum wall sits at x=-65.6 on the chord; chamber wall at
    // -17.0; cut spans -78..-8, stopping short of the far wall
    translate([-43, chord_p, z_slot0 + slot_h/2])
        cube([70, slot_w, slot_h], center=true);
}

/* ---------------- pinwheel guide fins ------------------ */
module fins() {
    // each fin extends its slot's +y guide wall outward; 0.2 mm
    // overlap into the channel avoids coplanar CSG ambiguity
    intersection() {
        union() {
            for (i = [0 : n_slots-1])
                rotate([0, 0, i*360/n_slots])
                    translate([-62 - fin_len, chord_p + slot_w/2 - 0.2, 0])
                        cube([fin_len + 0.2, fin_th, drum_h]);
        }
        cylinder(r=fin_clip_r, h=drum_h);
    }
}

/* ---------------- wind base ------------------ */
module wind_base() {
    union() {
        difference() {
            union() {
                cylinder(r=drum_or, h=drum_h);           // slotted drum
                cylinder(r=skirt_r, h=floor_th);         // stability skirt
                fins();
                // outer groove collar with sloped shoulder + entry flare
                rotate_extrude()
                    polygon([[groove_or, drum_h], [drum_or, drum_h],
                             [58, base_h], [52, base_h],
                             [groove_or, base_h-2]]);
                // inner groove collar, tapered tube lead-in
                rotate_extrude()
                    polygon([[chamber_r, drum_h], [groove_ir, drum_h],
                             [groove_ir, base_h-2], [chamber_r, base_h]]);
            }
            // swirl chamber + exit bore
            translate([0, 0, floor_th])
                cylinder(r=chamber_r, h=base_h + 5);
            // 12 tangential slots
            for (i = [0 : n_slots-1])
                rotate([0, 0, i*360/n_slots]) slot_cut();
            // 3 countersunk stake / screw holes in the skirt
            for (a = [10, 130, 250]) rotate([0, 0, a]) {
                translate([stake_r, 0, -1]) cylinder(d=4.5, h=5);
                translate([stake_r, 0, 1.4]) cylinder(d=9, h=2);
            }
        }
        pedestal();
    }
}

module pedestal() {
    difference() {
        translate([0, 0, floor_th])
            cylinder(r=ped_or, h=29);                 // top at z=32
        translate([0, 0, floor_th+3])
            cylinder(r=ped_or-ped_wall, h=21);        // sealed air gap
        translate([0, 0, 30])
            cylinder(r=recess_r, h=5);                // cap/foil recess
    }
}

/* ---------------- printable tube (same as indoor version) ---- */
module ptube() {
    difference() {
        cylinder(r=tube_or, h=tube_len);
        translate([0, 0, -1]) cylinder(r=tube_or-1.6, h=tube_len+2);
    }
}

/* ---------------- assembly ------------------ */
module assembly() {
    wind_base();
    %translate([0, 0, drum_h]) ptube();
}

if (part == 0) assembly();
if (part == 1) wind_base();
