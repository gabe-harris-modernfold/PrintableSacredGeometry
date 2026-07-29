// ============================================================
// Rainbow Marble Tray - two-layer "sandwich" version
// Marbles are locked between a cupped base plate and a
// countersunk cap plate. Cap grabs each marble just above its
// equator so ~90% of the dome stays exposed to the light.
//
// part = "base"     -> printable bottom plate (cups + bosses)
// part = "cap"      -> printable top plate (cone holes), prints flat
// part = "assembly" -> both plates + marbles, assembled
// part = "cutaway"  -> assembly sliced in half to show the lock
//
// Assembly: drop marbles in cups, set cap on the 6 corner
// bosses, drive M3 x 8 self-tapping screws into the bosses.
// ============================================================

part = "cutaway"; // ["base","cap","assembly","cutaway"]

/* [Marble] */
marble_d   = 10.0;  // marble diameter (0.39" = 10 mm)
clearance  = 0.3;   // cup oversize for marble tolerance
seat_depth = 3.5;   // marble depth below base surface

/* [Tile] */
pitch       = 11.0; // marble center-to-center spacing
tile_af     = 150;  // hexagon across flats
base_t      = 6;    // base plate thickness
cap_t       = 2.0;  // cap plate thickness
standoff_h  = 2.5;  // gap between plates (sets cap at marble shoulder)
edge_margin = 6;    // cups stay this far inside the edge

/* [Cap holes] */
hole_top    = 9.2;  // cap opening at top surface (< marble_d locks it in)
hole_bottom = 10.4; // cap opening at underside (clears the equator)

/* [Fasteners] */
boss_r_pos  = 76;   // bosses sit at hex corners, this far from center
boss_d      = 8;    // boss outer diameter (doubles as cap standoff)
pilot_d     = 2.8;  // M3 self-tap pilot
screw_d     = 3.4;  // cap through-hole
csk_d       = 6.4;  // countersink for screw head

$fn = 48;

// ------------------------------------------------------------
seat_d    = marble_d + clearance;
apothem   = tile_af / 2;
hex_r     = apothem / cos(30);
row_h     = pitch * sin(60);
a_in      = apothem - edge_margin;
marble_z  = base_t + marble_d/2 - seat_depth;      // marble center height
cap_z     = base_t + standoff_h;                   // cap underside height
boss_clear = 10;                                   // no cups this close to a boss

function in_hex(x, y, a) =
    abs(y) <= a &&
    abs(0.866025*x + 0.5*y) <= a &&
    abs(0.866025*x - 0.5*y) <= a;

boss_pos = [for (a = [0:60:300]) [boss_r_pos*cos(a), boss_r_pos*sin(a)]];

function d2(x, y, b) = (x-b[0])*(x-b[0]) + (y-b[1])*(y-b[1]);
function clear_of_bosses(x, y) =
    min([for (b = boss_pos) d2(x, y, b)]) > boss_clear*boss_clear;

pts = [
    for (iy = [-14:14], ix = [-14:14])
        let (x = ix*pitch + (abs(iy) % 2) * pitch/2,
             y = iy*row_h)
        if (in_hex(x, y, a_in) && clear_of_bosses(x, y)) [x, y]
];

echo(str("Marbles per tile: ", len(pts)));

// ------------------------------------------------------------
module base_plate() {
    difference() {
        union() {
            cylinder(h = base_t, r = hex_r, $fn = 6);
            for (b = boss_pos)
                translate([b[0], b[1], 0])
                    cylinder(h = base_t + standoff_h, d = boss_d);
        }
        for (p = pts)
            translate([p[0], p[1], base_t + seat_d/2 - seat_depth])
                sphere(d = seat_d);
        for (b = boss_pos)
            translate([b[0], b[1], base_t + standoff_h - 7])
                cylinder(h = 7.1, d = pilot_d);
    }
}

// cap is modeled flat on z=0 in print orientation
module cap_plate() {
    difference() {
        cylinder(h = cap_t, r = hex_r, $fn = 6);
        for (p = pts)
            translate([p[0], p[1], -0.05])
                cylinder(h = cap_t + 0.1, d1 = hole_bottom, d2 = hole_top);
        for (b = boss_pos) {
            translate([b[0], b[1], -0.05])
                cylinder(h = cap_t + 0.1, d = screw_d);
            translate([b[0], b[1], cap_t - 1.5])
                cylinder(h = 1.55, d1 = screw_d, d2 = csk_d);
        }
    }
}

module marbles() {
    for (p = pts)
        translate([p[0], p[1], marble_z])
            sphere(d = marble_d, $fn = 32);
}

module assembly() {
    color([0.14, 0.14, 0.15]) base_plate();
    color([0.22, 0.22, 0.24]) translate([0, 0, cap_z]) cap_plate();
    color([0.55, 0.72, 0.85, 0.75]) marbles();
}

// ------------------------------------------------------------
if (part == "base")     base_plate();
if (part == "cap")      cap_plate();
if (part == "assembly") assembly();
module cutbox(n) { translate([-500 - n, -500, -1 - n]) cube([1000 + 2*n, 500, 50 + 2*n]); }

if (part == "cutaway") {
    color([0.16, 0.16, 0.18]) difference() { base_plate(); cutbox(0); }
    color([0.38, 0.38, 0.42]) difference() { translate([0, 0, cap_z]) cap_plate(); cutbox(1); }
    color([0.45, 0.68, 0.88]) difference() { marbles(); cutbox(2); }
}
