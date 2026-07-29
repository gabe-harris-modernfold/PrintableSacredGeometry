// ============================================================
// Rainbow Marble Tray - printable hex tile
// One tile of the tabletop glass-bead-rainbow display.
// Print in matte black PETG, seat clear 10 mm (0.39") glass
// marbles in the cups. Tiles butt together edge-to-edge.
// ============================================================

/* [Marble] */
marble_d   = 10.0;   // marble diameter, mm (0.39" = 9.9-10 mm)
clearance  = 0.3;    // extra cup diameter for marble tolerance + print shrink
seat_depth = 3.5;    // how deep the marble sits below the tile surface

/* [Tile] */
pitch       = 11.0;  // center-to-center marble spacing (marble_d + wall)
tile_af     = 150;   // hexagon size across flats, mm (fits 180+ mm beds)
base_t      = 6;     // plate thickness
edge_margin = 6;     // keep cups this far inside the tile edge

/* [Quality] */
$fn = 48;

// ------------------------------------------------------------
seat_d   = marble_d + clearance;
apothem  = tile_af / 2;
hex_r    = apothem / cos(30);          // corner radius of the hexagon
row_h    = pitch * sin(60);            // hex-grid row spacing
a_in     = apothem - edge_margin;      // cups stay inside this apothem

// point-in-hexagon (flat-top hexagon, flats facing +/-Y)
function in_hex(x, y, a) =
    abs(y) <= a &&
    abs(0.866025*x + 0.5*y) <= a &&
    abs(0.866025*x - 0.5*y) <= a;

// all cup centers on a hex lattice, clipped to the inset hexagon
pts = [
    for (iy = [-14:14], ix = [-14:14])
        let (x = ix*pitch + (abs(iy) % 2) * pitch/2,
             y = iy*row_h)
        if (in_hex(x, y, a_in)) [x, y]
];

echo(str("Seats on this tile: ", len(pts)));

difference() {
    cylinder(h = base_t, r = hex_r, $fn = 6);   // flat-top hex plate
    for (p = pts)
        translate([p[0], p[1], base_t + seat_d/2 - seat_depth])
            sphere(d = seat_d);
}
