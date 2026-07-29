
// Sacred Geometry Generator (Metric, aligned to the same grid)
// All shapes share the same base radius R, allowing them to perfectly overlap.
// To generate STLs, uncomment the module you want at the bottom, press F6 to render, and File -> Export -> Export as STL.

R = 10;           // Base radius of the fundamental circles
T = 2;            // Thickness for 2D extruded shapes
$fn = 64;         // Resolution of curves

// Helper: 2D Ring
module ring(r, thickness=1) {
    difference() {
        circle(r=r);
        circle(r=r-thickness);
    }
}

// 1. A Single Circle
module single_circle() {
    linear_extrude(height=T, center=true)
    ring(R);
}

// 2. Two Intersecting Circles (Vesica Piscis)
module two_circles() {
    linear_extrude(height=T, center=true) {
        translate([-R/2, 0, 0]) ring(R);
        translate([R/2, 0, 0]) ring(R);
    }
}

// 3. Seed of Life (7 overlapping circles)
module seed_of_life() {
    linear_extrude(height=T, center=true) {
        ring(R);
        for (i = [0:5]) {
            translate([R * cos(i * 60), R * sin(i * 60), 0])
            ring(R);
        }
    }
}

// 4. Egg of Life (8 Spheres)
// The Egg of Life is 8 spheres arranged like the corners of a cube.
module egg_of_life() {
    d = R; 
    for (x = [-d/2, d/2]) {
        for (y = [-d/2, d/2]) {
            for (z = [-d/2, d/2]) {
                translate([x, y, z]) sphere(r=R/2);
            }
        }
    }
}

// 5. Fruit of Life (13 Circles: Center + 2 in each of the 6 directions)
module fruit_of_life() {
    linear_extrude(height=T, center=true) {
        ring(R);
        for (i = [0:5]) {
            translate([2 * R * cos(i * 60), 2 * R * sin(i * 60), 0]) ring(R);
            translate([4 * R * cos(i * 60), 4 * R * sin(i * 60), 0]) ring(R);
        }
    }
}

// Helper: Centers for Metatron's Cube
function metatron_centers() = [
    [0, 0, 0],
    for (i = [0:5]) [2 * R * cos(i * 60), 2 * R * sin(i * 60), 0],
    for (i = [0:5]) [4 * R * cos(i * 60), 4 * R * sin(i * 60), 0]
];

// 6. Metatron's Cube 
module metatrons_cube() {
    fruit_of_life();
    pts = metatron_centers();
    linear_extrude(height=T/2, center=true) {
        for (i = [0:12]) {
            for (j = [i+1:12]) {
                hull() {
                    translate([pts[i][0], pts[i][1], 0]) circle(r=T/4);
                    translate([pts[j][0], pts[j][1], 0]) circle(r=T/4);
                }
            }
        }
    }
}

// 7. Tree of Life (Aligned to the Seed of Life grid)
module tree_of_life() {
    // Standard coordinates scaling to match the Seed of Life intersections
    nodes = [
        [0, 2*R, 0],                 // Kether
        [R*cos(30), R*sin(30)+R, 0], // Chokhmah
        [-R*cos(30), R*sin(30)+R, 0],// Binah
        [R*cos(30), R*sin(30)-R, 0], // Chesed
        [-R*cos(30), R*sin(30)-R, 0],// Gevurah
        [0, 0, 0],                   // Tiferet
        [R*cos(30), R*sin(30)-3*R, 0], // Netzach
        [-R*cos(30), R*sin(30)-3*R, 0],// Hod
        [0, -2*R, 0],                // Yesod
        [0, -4*R, 0]                 // Malkuth
    ];

    linear_extrude(height=T, center=true) {
        for (p = nodes) {
            translate([p[0], p[1], 0]) circle(r=R/3);
        }
    }

    connections = [
        [0,1], [0,2], [0,5], [1,2], [1,3], [1,5], [2,4], [2,5],
        [3,4], [3,5], [3,6], [4,5], [4,7], [5,6], [5,7], [5,8],
        [6,7], [6,8], [7,8], [8,9]
    ];

    linear_extrude(height=T/2, center=true) {
        for (c = connections) {
            hull() {
                translate([nodes[c[0]][0], nodes[c[0]][1], 0]) circle(r=T/4);
                translate([nodes[c[1]][0], nodes[c[1]][1], 0]) circle(r=T/4);
            }
        }
    }
}

// 8. Platonic Solids (Basic Wireframes scaled to the same grid)
module tetrahedron() {
    cylinder(h=4*R, r1=4*R, r2=0, $fn=3, center=true);
}
module hexahedron() {
    cube(4*R, center=true);
}
module octahedron() {
    intersection() {
        cylinder(h=4*R, r1=4*R, r2=0, $fn=4, center=true);
        rotate([180,0,0]) cylinder(h=4*R, r1=4*R, r2=0, $fn=4, center=true);
    }
}

// Render ALL aligned together (so you can see they fit the same bounding box)
module render_all() {
    single_circle();
    color("red") two_circles();
    color("green") seed_of_life();
    color("blue") egg_of_life();
    color("orange") fruit_of_life();
    color("purple") metatrons_cube();
    color("cyan") tree_of_life();
}

// Uncomment the one you want to export as STL:

// single_circle();
// two_circles();
// seed_of_life();
// egg_of_life();
// fruit_of_life();
// metatrons_cube();
// tree_of_life();

// render_all();
