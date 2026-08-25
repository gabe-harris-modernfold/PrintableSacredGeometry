// Spindle torus — hollow shell, sliced vertically in half.
// Spindle torus: tube radius r > ring radius R, so the tube overlaps
// itself through the center. Hollow 2mm wall, cut on a plane through
// the spin axis, laid flat on the cut face for printing.
//
// Footprint on bed (cut face down): 2*(R+r) x 2*r, height R+r.
// Defaults: 300 x 200 mm footprint, 150 mm tall -> fits a 340x320 bed.

R    = 50;    // ring radius: axis -> tube center
r    = 100;   // tube radius (must be > R for a spindle)
wall = 2;     // shell wall thickness (mm)

lay_flat = true;  // true: rest on cut face for printing; false: upright

// Smoothness
seg_circle = 160;  // facets around the tube cross-section
seg_revolve = 240; // facets around the revolution

// --- hollow tube cross-section in the X>=0 half-plane ---------------
// Model the torus TUBE as a hollow pipe (annulus: outer r, inner r-wall)
// and revolve it. Because r > R the pipe reaches past the axis, so the
// revolution self-intersects: the part of the pipe that crosses to -X
// reappears as the mirrored annulus on +X. Their union, revolved, gives
// both the outer wall and the inner spindle wall through the center.
// Clip to X>=0 so rotate_extrude is valid.
module annulus(tube) {
    difference() {
        circle(r = tube,        $fn = seg_circle);
        circle(r = tube - wall, $fn = seg_circle);
    }
}

module half_tube_profile() {
    intersection() {
        union() {
            translate([ R, 0]) annulus(r);
            translate([-R, 0]) annulus(r);
        }
        translate([0, -(R + r) - 1]) square([R + r + 1, 2 * (R + r) + 2]);
    }
}

module spindle_shell() {
    rotate_extrude($fn = seg_revolve) half_tube_profile();
}

// Keep the half with y >= 0; the cut face lies on the XZ plane.
module half_shell() {
    big = 2 * (R + r) + 10;
    intersection() {
        spindle_shell();
        translate([-big/2, 0, -big/2]) cube([big, big, big]);
    }
}

if (lay_flat)
    // rotate cut face (y=0) down onto the bed (z=0)
    rotate([90, 0, 0]) half_shell();
else
    half_shell();
