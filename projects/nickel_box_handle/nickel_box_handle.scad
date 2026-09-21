// Nickel-box carry handle -- slip-on end cap for a bank box of nickels.
// Print two: one on each end of the box gives a two-handled ~10 kg workout weight.
//
// The part is an extruded end-view profile (rectangular sleeve + D-handle on top)
// closed by an end plate. Print with the end plate on the bed (as modelled):
// every feature grows straight up from the bed, no supports, no bridges.
// Loads in use are in the layer plane, the strongest direction.

// ---- box (nominal cardboard) ----------------------------------------------
BOX_W   = 88;    // width of the end face (Y)
BOX_H   = 100;   // height of the end face (Z)
BOX_L   = 240;   // length (X) -- for the preview only

// ---- cap ------------------------------------------------------------------
CLEAR   = 1.5;   // per-side clearance to the cardboard
WALL    = 3.6;   // sleeve wall
PLATE   = 4;     // end plate thickness
DEPTH   = 70;    // how far the sleeve reaches along the box (incl. plate)
LEAD_IN = 1.5;    // 45-degree chamfer at the mouth for easy slip-on
CORNER  = 4;     // outer vertical-edge rounding

// ---- handle ---------------------------------------------------------------
HANDLE_T   = 30;   // thickness along the box (X), starts flush with end face
GRIP_H     = 24;   // grip bar height
FINGERS    = 4;    // scalloped finger grooves under the grip bar
GROOVE_P   = 22;   // groove pitch
GROOVE_R   = 10;   // groove radius
GROOVE_D   = 3;    // groove depth
EAR_R      = 12;   // paracord ears at the lower corners of the hoop (holes run along the box); must exceed EAR_HOLE/2 + R_EDGE
EAR_HOLE   = 7;    // hole diameter (4 mm paracord, room for a doubled pass or knot)
OPEN_H     = 44;   // finger clearance under the grip bar
OPEN_W_TOP = 112;  // opening width just under the grip bar
UPRIGHT    = 14;   // width of the flared uprights
GRIP_R     = 8;    // grip bar edge rounding
GUSSET     = 25;   // 45-ish degree fillet running from the handle back along the top wall
BAND       = 6;    // the handle section continues as a hoop this much thicker than the sleeve, all the way round
R_EDGE     = 3;    // rounding of every outside edge (PETG-friendly, no sharp corners)

$fn = 48;

IN_W = BOX_W + 2*CLEAR;      // cavity
IN_H = BOX_H + 2*CLEAR;
OUT_W = IN_W + 2*WALL;       // sleeve outside
OUT_H = IN_H + 2*WALL;
TOP   = OUT_H/2;             // z of the top wall outer face (cap centred on origin in Y,Z)

module rrect(w, h, r) { offset(r) offset(-r) square([w, h], center=true); }

// Extrude a 2D profile along +X with every outside edge rounded by R_EDGE.
// Stepped fillet: R_STEPS thin slices at each end, each an offset of the profile that
// follows a quarter circle. With 0.5 mm steps it prints identically to a true round,
// and unlike minkowski() with a sphere it renders in seconds on OpenSCAD 2021/CGAL.
R_STEPS = 6;
module rounded_extrude_x(h) {
    rotate([0, 90, 0]) rotate([0, 0, 90]) {
        for (i = [1 : R_STEPS]) {
            a0 = 90 * (i - 1) / R_STEPS;   a1 = 90 * i / R_STEPS;
            h0 = R_EDGE * (1 - cos(a0));   h1 = R_EDGE * (1 - cos(a1));
            shrink = R_EDGE * (1 - sin(a0));                 // inside the true round
            for (z0 = [h0, h - h1])
                translate([0, 0, z0]) linear_extrude(h1 - h0 + 0.01)
                    offset(r = -shrink) children();
        }
        translate([0, 0, R_EDGE]) linear_extrude(h - 2*R_EDGE) children();
    }
}

// End-view profile of the handle section: a hoop BAND thicker than the sleeve running
// all the way round the box, with the D-loop growing out of its top.
module handle_profile() {
    hw_base = OUT_W/2 + BAND;              // half-width where the loop leaves the hoop
    hw_top  = OPEN_W_TOP/2 + UPRIGHT;      // half-width at grip level
    z0      = TOP + BAND;                  // hoop outer face on top = floor of the finger opening
    z_grip0 = z0 + OPEN_H;
    z_grip1 = z_grip0 + GRIP_H;
    difference() {
        union() {
            rrect(OUT_W + 2*BAND, OUT_H + 2*BAND, CORNER + BAND);
            offset(GRIP_R) offset(-GRIP_R)
                polygon([[-hw_base, TOP - 8], [ hw_base, TOP - 8],
                         [ hw_top, z_grip0], [ hw_top, z_grip1],
                         [-hw_top, z_grip1], [-hw_top, z_grip0]]);
        }
        // finger opening: flared trapezoid down to the hoop's top face
        offset(GRIP_R/2) offset(-GRIP_R/2)
            polygon([[-(hw_base - UPRIGHT), z0], [ (hw_base - UPRIGHT), z0],
                     [ OPEN_W_TOP/2, z_grip0], [-OPEN_W_TOP/2, z_grip0]]);
        // finger grooves scalloped into the underside of the grip bar
        for (i = [0 : FINGERS - 1])
            translate([(i - (FINGERS - 1)/2) * GROOVE_P, z_grip0 - GROOVE_R + GROOVE_D])
                circle(GROOVE_R);
        // paracord holes
        for (sx = [-1, 1]) translate(ear_centre(sx)) circle(d=EAR_HOLE);
    }
    // ears are added outside the difference so their material is never cut by the opening
}

function ear_centre(sx) = [sx * (OUT_W/2 + BAND + EAR_R - 2), -(TOP + BAND) + EAR_R];

module handle_profile_with_ears() {
    difference() {
        union() {
            handle_profile();
            for (sx = [-1, 1]) hull() {
                translate(ear_centre(sx)) circle(EAR_R);
                // blend into the hoop side wall
                translate([sx * (OUT_W/2 + BAND - 4), -(TOP + BAND) + EAR_R]) square([8, 2*EAR_R], center=true);
            }
        }
        for (sx = [-1, 1]) translate(ear_centre(sx)) circle(d=EAR_HOLE);
    }
}

module cap() {
    difference() {
        union() {
            // sleeve + end plate
            rounded_extrude_x(DEPTH) rrect(OUT_W, OUT_H, CORNER);
            // handle loop, flush with the end face
            rounded_extrude_x(HANDLE_T) handle_profile_with_ears();
            // gusset: the loop's footprint on the top wall, tapered back along X
            hull() {
                rounded_extrude_x(HANDLE_T) intersection() {
                    handle_profile();
                    translate([-OUT_W/2 - BAND, TOP - 8]) square([OUT_W + 2*BAND, GUSSET + 8 + BAND]);
                }
                translate([HANDLE_T + GUSSET - R_EDGE, 0, TOP]) rotate([90, 0, 0])
                    cylinder(r=R_EDGE, h=OUT_W - 2*CORNER - 2*R_EDGE, center=true, $fn=16);
            }
        }
        // cavity
        translate([PLATE, -IN_W/2, -IN_H/2]) cube([DEPTH, IN_W, IN_H]);
        // mouth lead-in
        translate([DEPTH, 0, 0]) rotate([0, -90, 0]) rotate([0, 0, 90])
            linear_extrude(LEAD_IN, scale=[IN_W/(IN_W + 2*LEAD_IN), IN_H/(IN_H + 2*LEAD_IN)])
                square([IN_W + 2*LEAD_IN, IN_H + 2*LEAD_IN], center=true);
    }
}

// PART: end plate on the bed
module part() { rotate([0, -90, 0]) cap(); }

// PREVIEW: box with a cap on each end
module assembly() {
    color("BurlyWood") translate([0, -BOX_W/2, -BOX_H/2]) cube([BOX_L, BOX_W, BOX_H]);
    color("SteelBlue") translate([-PLATE - 0.01, 0, 0]) cap();
    color("SteelBlue") translate([BOX_L + PLATE + 0.01, 0, 0]) mirror([1, 0, 0]) cap();
}

MODE = "part";   // overridden from the command line with -D MODE='"assembly"'
if (MODE == "assembly") assembly(); else part();
