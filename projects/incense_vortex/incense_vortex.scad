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
     3 = printable tube   (reference only at this size: 394 mm
                           exceeds the 320 mm bed — use the acrylic)
   *mouthpiece boss: 3 mm of bridging, fine on any printer.

   Swirl is CLOCKWISE viewed from above. Tube: 100 mm OD acrylic,
   2 mm wall, 15.5" (393.7 mm) long. PETG or ABS recommended.
   ============================================================ */

part = 0;

/* [Tube fit] */
tube_od    = 100;    // 100 mm OD acrylic
tube_wall  = 2;      // 2 mm acrylic wall
tube_len   = 393.7;  // 15.5 inch tube
fit_gap    = 0.4;    // radial slip-fit clearance

/* [Swirl geometry] */
chamber_r  = 45.5;   // swirl chamber radius (= ceiling exit hole)
vane_angle = 68;     // entry angle from radial; tan = swirl ratio
n_slots    = 12;
slot_w     = 10;     // slot channel width
slot_h     = 41;     // slot height (12*10*41 / pi*45.5^2 = 0.76 bore)
base_or    = 88;     // base drum outer radius (176 mm plinth)

/* [Base build] */
floor_th   = 3;
ceil_th    = 3;
collar_h   = 10;     // groove wall height above ceiling
ped_or     = 17;     // cone pedestal radius
ped_wall   = 2.4;    // pedestal shell (hollow = heat break)
recess_r   = 15;     // fits 30 mm metal cap / foil disc
recess_d   = 2;

/* [Top collar] */
exit_r     = 35;     // nozzle exit radius (70 mm opening)

/* [Printable tube] */
ptube_wall = 1.6;

/* derived */
tube_or   = tube_od/2;            // 50
tube_ir   = tube_or - tube_wall;  // 48
groove_or = tube_or + fit_gap;    // 50.4 outer groove wall (inner face)
groove_ir = tube_ir - fit_gap;    // 47.6 inner collar outer face
chord_p   = chamber_r*sin(vane_angle);   // 42.19 -> 68 deg entry
z_slot0   = floor_th;                    // 3
z_slot1   = floor_th + slot_h;           // 44
drum_h    = z_slot1 + ceil_th;           // 47  (groove floor)
base_h    = drum_h + collar_h;           // 57

$fa = 3;
$fs = 0.5;

/* ---------------- swirl slot (one chord channel) ------------- */
module slot_cut() {
    // one-sided chord cut: pierces outer wall and chamber wall on
    // the -x side only, so every slot injects the SAME handedness
    // (outer wall sits at x=-77.2 on the chord; chamber wall at
    //  -17.0; cut spans -84..-8, stopping short of the far wall)
    translate([-46, chord_p, z_slot0 + slot_h/2])
        cube([76, slot_w, slot_h], center=true);
}

/* ---------------- mouthpiece bore (slot #0) ------------------ */
module breath_bore() {
    // 9 mm bore coaxial with slot 0's chord; blow here to add swirl
    translate([-89, chord_p, 12]) rotate([0, 90, 0])
        cylinder(d=9, h=55);
}

module breath_boss() {
    translate([-86, chord_p, 12]) rotate([0, 90, 0])
        cylinder(d=18, h=22);
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
                             [58, base_h], [52, base_h],
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
            [36.25 + 1.25*cos(a), 27 + 1.25*sin(a)] ];
    rotate_extrude()
        polygon(concat(
            [[47, 8], [50.4, 8], [50.4, 0], [53, 0], [53, 10],
             [37.5, 27]],
            lip,
            [[35, 24], [46.6, 12.4], [47, 12.4]]
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
