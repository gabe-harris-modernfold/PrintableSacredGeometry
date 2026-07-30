/* ============================================================
   WIND-POWERED INCENSE TORNADO — omnidirectional wind harvester
   ============================================================
   Variant of the passive incense tornado optimised to CATCH THE
   WIND and convert it into swirl. Two harvesting stages:

   INTAKE PINWHEEL (base): 12 tall guide fins radiate from the
   drum, each extending one tangential slot's guide wall. From any
   wind azimuth, the windward fins form converging funnels (~4:1
   area contraction) that accelerate the breeze and inject it as
   same-handed tangential jets. Ram pressure of a 2-3 m/s breeze
   (~5 Pa) is ~15x the cone's stack-effect draft.

   VENTURI COWL (top): stacked-disc chimney-cowl. Crosswind is
   squeezed through the converging gap between a 45-degree skirt
   and a flat hat, dropping static pressure right over the exit
   throat (Bernoulli) and sucking the column upward from any wind
   direction. Also a rain / downdraft hat.

   In dead calm it still works like the indoor version (heat-only,
   gentler). Swirl is CLOCKWISE viewed from above.

   Parts (set `part`):
     0 = assembly preview (tube ghosted)
     1 = wind base: slotted drum + pinwheel fins + staked skirt
     2 = venturi head: slips over tube top (print as oriented)
     3 = venturi hat: press-fits onto the head's three posts
   Tube: 3" (76.2 mm OD) acrylic, or reuse vortex_tube_printable.stl.
   ============================================================ */

part = 0;

/* [Tube fit] */
tube_od    = 76.2;
tube_wall  = 3.175;
tube_len   = 220;
fit_gap    = 0.4;

/* [Swirl geometry] */
chamber_r  = 32.5;
vane_angle = 68;
n_slots    = 12;
slot_w     = 7;
slot_h     = 34;     // taller than indoor version: more catch area
drum_or    = 55;

/* [Wind intake] */
fin_len    = 35;     // guide fin length beyond the drum
fin_th     = 2.47;
fin_clip_r = 84;     // fins trimmed to this radius
skirt_r    = 86;     // stability skirt
stake_r    = 79;     // bolt-down / tent-stake holes

/* [Base build] */
floor_th   = 3;
ceil_th    = 3;
collar_h   = 10;
ped_or     = 17;
ped_wall   = 2.4;
recess_r   = 15;
recess_d   = 2;

/* [Venturi cowl] */
throat_r   = 23;     // exit throat radius
post_r_pos = 30;     // post circle radius
post_d     = 6;
post_len   = 17.5;   // 13 mm gap + 4.5 mm socket
hat_r      = 55;

/* derived */
tube_or   = tube_od/2;
tube_ir   = tube_or - tube_wall;
groove_or = tube_or + fit_gap;           // 38.5
groove_ir = tube_ir - fit_gap;           // 34.5
chord_p   = chamber_r*sin(vane_angle);   // 30.13 -> 68 deg entry
z_slot0   = floor_th;                    // 3
z_slot1   = floor_th + slot_h;           // 37
drum_h    = z_slot1 + ceil_th;           // 40 (groove floor)
base_h    = drum_h + collar_h;           // 50

$fa = 3;
$fs = 0.5;

/* ---------------- swirl slot (one chord channel) ------------- */
module slot_cut() {
    translate([-32.75, chord_p, z_slot0 + slot_h/2])
        cube([48.5, slot_w, slot_h], center=true);
}

/* ---------------- pinwheel guide fins ------------------ */
module fins() {
    // each fin extends its slot's +y guide wall outward; 0.2 mm
    // overlap into the channel avoids coplanar CSG ambiguity
    intersection() {
        union() {
            for (i = [0 : n_slots-1])
                rotate([0, 0, i*360/n_slots])
                    translate([-42 - fin_len, chord_p + slot_w/2 - 0.2, 0])
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
                             [46, base_h], [40, base_h],
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

/* ---------------- venturi head ------------------ */
module venturi_head() {
    union() {
        rotate_extrude()
            polygon([
                [34, 8],       // ledge inner, underside
                [38.5, 8],
                [38.5, 0],     // sleeve inner wall
                [41, 0],
                [41, 10],      // sleeve outer top
                [30, 22],      // dome outer -> skirt launch
                [48, 4],       // skirt underside, 45 deg out-down
                [50, 5.5],     // skirt rim
                [32.5, 23],    // skirt topside = venturi lower disc
                [28, 23],      // flat post band
                [25.5, 24.5],  // throat lip chamfer
                [23, 24.5],    // throat lip inner edge
                [25, 21],      // bore below lip
                [36, 9.6],     // interior dome, 44 deg
                [34, 9.6]
            ]);
        for (a = [0, 120, 240]) rotate([0, 0, a])
            translate([post_r_pos, 0, 22])
                cylinder(d=post_d, h=post_len + 1);   // seats into band
    }
}

/* ---------------- venturi hat ------------------ */
module venturi_hat() {
    difference() {
        union() {
            cylinder(r=hat_r, h=3);
            translate([0, 0, 3]) cylinder(r1=hat_r, r2=0.5, h=7);
        }
        for (a = [0, 120, 240]) rotate([0, 0, a])
            translate([post_r_pos, 0, -0.5])
                cylinder(d=post_d + 0.6, h=5);        // press-fit sockets
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
    translate([0, 0, drum_h + tube_len - 8]) venturi_head();
    translate([0, 0, drum_h + tube_len - 8 + 36]) venturi_hat();
}

if (part == 0) assembly();
if (part == 1) wind_base();
if (part == 2) venturi_head();
if (part == 3) venturi_hat();
