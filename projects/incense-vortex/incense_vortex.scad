/* ============================================================
   PASSIVE INCENSE TORNADO — heat-driven helical smoke vortex
   ============================================================
   No fans, no electronics. The burning cone's own buoyant plume
   creates a stack-effect draft inside the clear tube; makeup air
   is forced through 12 tangential slots in the base drum, which
   impose angular momentum on the inflow (slotted-chimney fire-
   tornado principle). Conservation of angular momentum tightens
   the rotation toward the axis, wrapping the smoke into a helix.
   One slot carries an external mouthpiece nozzle: a gentle breath
   injects extra swirl to spin the column up on demand.

   Parts (set `part`):
     0 = assembly preview (tube ghosted)
     1 = base drum        (print upright, no supports needed*)
     2 = top collar       (print as oriented, flange down)
     3 = printable tube   (only if you don't buy 3" acrylic;
                           clear PETG, 0.8 nozzle, spiral or 2 walls)
   *mouthpiece boss: 3 mm of bridging, fine on any printer.

   Swirl is CLOCKWISE viewed from above. Tube: 3" (76.2 mm OD)
   acrylic, 1/8" wall — or part 3. PETG or ABS recommended.
   ============================================================ */

part = 0;

/* [Tube fit] */
tube_od    = 76.2;   // 3 inch acrylic
tube_wall  = 3.175;  // 1/8 inch acrylic wall
tube_len   = 220;    // for printable tube / preview
fit_gap    = 0.4;    // radial slip-fit clearance

/* [Swirl geometry] */
chamber_r  = 32.5;   // swirl chamber radius (= ceiling exit hole)
vane_angle = 68;     // entry angle from radial; tan = swirl ratio
n_slots    = 12;
slot_w     = 7;      // slot channel width
slot_h     = 30;     // slot height
base_or    = 60;     // base drum outer radius (120 mm plinth)

/* [Base build] */
floor_th   = 3;
ceil_th    = 3;
collar_h   = 10;     // groove wall height above ceiling
ped_or     = 17;     // cone pedestal radius
ped_wall   = 2.4;    // pedestal shell (hollow = heat break)
recess_r   = 15;     // fits 30 mm metal cap / foil disc
recess_d   = 2;

/* [Top collar] */
exit_r     = 25;     // nozzle exit radius (50 mm opening)

/* [Printable tube] */
ptube_wall = 1.6;

/* derived */
tube_or   = tube_od/2;            // 38.1
tube_ir   = tube_or - tube_wall;  // 34.925
groove_or = tube_or + fit_gap;    // 38.5 outer groove wall (inner face)
groove_ir = tube_ir - fit_gap;    // 34.5 inner collar outer face
chord_p   = chamber_r*sin(vane_angle);   // chord offset -> 68 deg entry
z_slot0   = floor_th;                    // 3
z_slot1   = floor_th + slot_h;           // 33
drum_h    = z_slot1 + ceil_th;           // 36  (groove floor)
base_h    = drum_h + collar_h;           // 46

$fa = 3;
$fs = 0.5;

/* ---------------- swirl slot (one chord channel) ------------- */
module slot_cut() {
    // one-sided chord cut: pierces outer wall and chamber wall on
    // the -x side only, so every slot injects the SAME handedness
    translate([-32.75, chord_p, z_slot0 + slot_h/2])
        cube([48.5, slot_w, slot_h], center=true);
}

/* ---------------- mouthpiece bore (slot #0) ------------------ */
module breath_bore() {
    // 9 mm bore coaxial with slot 0's chord; blow here to add swirl
    translate([-60, chord_p, 12]) rotate([0, 90, 0])
        cylinder(d=9, h=45);
}

module breath_boss() {
    translate([-58, chord_p, 12]) rotate([0, 90, 0])
        cylinder(d=18, h=18);
}

/* ---------------- base ------------------ */
module base() {
    union() {
        difference() {
            union() {
                cylinder(r=base_or, h=drum_h);            // main drum
                breath_boss();
                // outer groove collar with sloped shoulder + entry flare
                rotate_extrude()
                    polygon([[groove_or, drum_h], [base_or, drum_h],
                             [46, base_h], [40, base_h],
                             [groove_or, base_h-2]]);
                // inner groove collar, top tapered for tube lead-in
                rotate_extrude()
                    polygon([[chamber_r, drum_h], [groove_ir, drum_h],
                             [groove_ir, base_h-2], [chamber_r, base_h]]);
            }
            // swirl chamber + ceiling exit hole
            translate([0, 0, floor_th])
                cylinder(r=chamber_r, h=base_h);
            // 12 tangential slots
            for (i = [0 : n_slots-1])
                rotate([0, 0, i*360/n_slots]) slot_cut();
            breath_bore();
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

/* ---------------- top collar ------------------ */
module top_collar() {
    lip = [ for (a = [0:15:180])
            [26.25 + 1.25*cos(a), 24 + 1.25*sin(a)] ];
    rotate_extrude()
        polygon(concat(
            [[34, 8], [38.5, 8], [38.5, 0], [41, 0], [41, 10],
             [27.5, 24]],
            lip,
            [[25, 21], [34.6, 9.6], [34, 9.6]]
        ));
}

/* ---------------- printable clear tube ------------------ */
module ptube() {
    difference() {
        cylinder(r=tube_or, h=tube_len);
        translate([0, 0, -1])
            cylinder(r=tube_or - ptube_wall, h=tube_len+2);
    }
}

/* ---------------- assembly ------------------ */
module assembly() {
    base();
    %translate([0, 0, drum_h]) ptube();               // ghost tube
    translate([0, 0, drum_h + tube_len - 8]) top_collar();
}

if (part == 0) assembly();
if (part == 1) base();
if (part == 2) top_collar();
if (part == 3) ptube();
