// Horn torus — hollow 2mm shell, sliced in half, with a twisted vortex.
//
// Vortex = N flutes on the surface + a twist about the vertical axis
// proportional to height z. z flips sign across the equator, so the top
// spirals one way and the bottom the other (CW top / CCW bottom).
//
// Built as a single parametric polyhedron (cheap memory, smooth & high
// res). The model spans u = 0..180 deg, i.e. it is ALREADY the half cut
// on the vertical (XZ) plane. R is set just above r so the tube pinches
// to a ~1mm pinhole at the center instead of self-intersecting — a clean
// manifold that still reads as a horn torus.

r    = 75;    // tube radius
wall = 2;     // shell wall thickness
amp  = 2.5;   // flute depth (mm)
N    = 16;    // number of flutes (vortex arms)
R    = r + amp + 0.5;  // ring radius -> outer surface grazes center (~1mm hole)
twist_per_mm = 1.0;    // deg of twist per mm of height
cw_on_top    = true;   // true: CW vortex on top, CCW on bottom

Nu = 144;  // segments around the half revolution (twist smoothness)
Nv = 96;   // segments around the tube
lay_flat = true;

sgn = cw_on_top ? -1 : 1;

function uu(i) = 180 * i / Nu;
function vv(j) = 360 * j / Nv;
function flute(i, j) = amp * cos(N * (uu(i) + sgn * twist_per_mm * (r * sin(vv(j)))));
function vert(i, j, base) =
    let (rh = base + flute(i, j), v = vv(j), u = uu(i), RR = R + rh * cos(v))
        [RR * cos(u), RR * sin(u), rh * sin(v)];

outer = [for (i = [0:Nu]) for (j = [0:Nv-1]) vert(i, j, r)];
inner = [for (i = [0:Nu]) for (j = [0:Nv-1]) vert(i, j, r - wall)];
verts = concat(outer, inner);

no = (Nu + 1) * Nv;                  // index offset to the inner sheet
function oi(i, j) = i * Nv + (j % Nv);
function ii(i, j) = no + i * Nv + (j % Nv);

// outer skin
of = [for (i = [0:Nu-1]) for (j = [0:Nv-1]) each [
    [oi(i,j), oi(i+1,j), oi(i+1,j+1)],
    [oi(i,j), oi(i+1,j+1), oi(i,j+1)] ]];
// inner skin (reversed winding -> normals face the cavity)
inf = [for (i = [0:Nu-1]) for (j = [0:Nv-1]) each [
    [ii(i,j), ii(i+1,j+1), ii(i+1,j)],
    [ii(i,j), ii(i,j+1), ii(i+1,j+1)] ]];
// cap the open ring at u=0
cap0 = [for (j = [0:Nv-1]) each [
    [oi(0,j), oi(0,j+1), ii(0,j+1)],
    [oi(0,j), ii(0,j+1), ii(0,j)] ]];
// cap the open ring at u=180 (reversed)
capN = [for (j = [0:Nv-1]) each [
    [oi(Nu,j), ii(Nu,j+1), oi(Nu,j+1)],
    [oi(Nu,j), ii(Nu,j), ii(Nu,j+1)] ]];

module vortex() polyhedron(points = verts, faces = concat(of, inf, cap0, capN), convexity = 12);

if (lay_flat)
    rotate([90, 0, 0]) vortex();  // drop the y=0 cut face onto the bed
else
    vortex();
