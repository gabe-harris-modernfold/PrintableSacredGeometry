// ============================================================
// Rainbow Marble Tray - full tabletop scene (visualization only)
// 0.6 m field of hex tiles + 10 mm clear marbles.
// Marbles inside the rainbow band (light held ~600 mm above,
// eye next to the light) are colored by the actual glass
// rainbow angles: violet ~20.2 deg -> red ~22.0 deg.
// ============================================================

marble_d = 10.0;
pitch    = 11.0;
row_h    = pitch * sin(60);
base_t   = 6;
seat_depth = 3.5;

field_r  = 295;          // marble field radius (0.6 m tabletop)
tile_af  = 150;
hex_r    = (tile_af/2) / cos(30);

light_h  = 600;                      // light height above tray, mm
r_violet = light_h * tan(20.24);     // n = 1.53 (blue/violet in glass)
r_red    = light_h * tan(22.0);      // n = 1.51 (red in glass)

$fn = 14;                            // low-poly spheres, ~2600 of them

// hue (0-300) -> rgb, s=v=1
function h2rgb(h) =
    let (hh = h/60, i = floor(hh), f = hh - i)
    i == 0 ? [1, f, 0] :
    i == 1 ? [1-f, 1, 0] :
    i == 2 ? [0, 1, f] :
    i == 3 ? [0, 1-f, 1] :
             [f, 0, 1];

// ---- tiles: flat-top hexes tiled to cover the field ----
tile_dx = 1.5 * hex_r;
tile_dy = tile_af;
color([0.13, 0.13, 0.14])
for (col = [-3:3], row = [-3:3]) {
    tx = col * tile_dx;
    ty = row * tile_dy + (abs(col) % 2) * tile_af/2;
    if (sqrt(tx*tx + ty*ty) < field_r + 60)
        translate([tx, ty, 0])
            cylinder(h = base_t, r = hex_r - 0.6, $fn = 6);
}

// ---- marbles on one continuous hex lattice ----
for (iy = [-31:31], ix = [-31:31]) {
    x = ix*pitch + (abs(iy) % 2) * pitch/2;
    y = iy*row_h;
    r = sqrt(x*x + y*y);
    if (r <= field_r) {
        in_band = (r >= r_violet - 4 && r <= r_red + 4);
        t = min(1, max(0, (r - r_violet) / (r_red - r_violet)));
        c = in_band ? h2rgb(270 * (1 - t)) : [0.72, 0.76, 0.80];
        a = in_band ? 1.0 : 0.45;
        color([c[0], c[1], c[2], a])
            translate([x, y, base_t + marble_d/2 - seat_depth])
                sphere(d = marble_d);
    }
}
