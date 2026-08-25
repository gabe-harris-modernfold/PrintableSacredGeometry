#!/usr/bin/env python3
"""Generate printable STLs for the 7 crystal systems, the 14 Bravais lattices and the
three cubic Wigner-Seitz cells.

Five output sets
----------------
supercell/      Each lattice as a block of cells rather than one cell, because periodicity
                is the whole content of the word "lattice" and a single cell cannot show it.
                Everything shared between neighbours is merged exactly once, so a block of
                8 cubic P cells is 27 hubs and 54 struts rather than 64 and 96 -- an interior
                corner is one hub belonging to eight cells at a time, which is precisely why
                a primitive cell holds one lattice point and not eight. The unit cell at the
                origin is drawn heavier so you can see what is repeating.

solids/         14 solid unit cells, longest edge 60 mm, resting on the a-b face with the
                Pearson symbol embossed on the a-c face. This is the comparison set: it
                shows the axial ratios and interaxial angles, which is what the crystal
                system actually is. All 14 fit one 320 mm plate.

assembled/      One cell each as a one-piece ball-and-stick frame, longest edge 60 mm --
                the reference object the supercell repeats. `orientations` and
                `frame_printability` pick which cell face each one stands on.

kit/            The same 14 cells as ball-and-stick frames at 90 mm, but supplied as
                sleeved hub nodes and plain rods. More work (387 parts) and it buys
                something the fused version cannot do: the cells come apart and join into
                supercells -- three hP cells make the familiar hexagonal prism, and
                stacking a cI or cF cell shows the centring translation actually tiling.

wigner_seitz/   Cube, truncated octahedron and rhombic dodecahedron at a shared 60 mm
                lattice constant. These are the primitive cells of cP, cI and cF drawn as
                the region closer to one lattice point than to any other, so their volumes
                are exactly a^3, a^3/2 and a^3/4 -- and cI's cell is cF's Brillouin zone
                and vice versa. All three space-fill.

Why the nodes are sleeves rather than drilled sockets
----------------------------------------------------
There is no mesh-boolean backend in this environment, so a socket cannot be subtracted
from a sphere. Instead each socket is an open-bore tube standing off the hub, starting
just inside the hub surface: the hub closes the bottom of the bore, giving a blind socket
6 mm deep with no CSG at all. It also removes the reason to enlarge centring nodes -- a
drilled body-centre node would have had its core eaten by eight intersecting holes, but a
sleeved one does not, so every node in the set is the same Ø14 hub. Lattice points in a
Bravais lattice are all equivalent, and now the parts are too.

Outputs
-------
stl/supercell/*.stl, stl/solids/*.stl, stl/assembled/*.stl, stl/kit/<pearson>/*.stl,
stl/wigner_seitz/*.stl
PARTS.md        per-lattice part list, staged into two build plates
"""

import itertools
import math
import os
import sys

import numpy as np
from scipy.spatial import ConvexHull

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))
from mesh_kit import Mesh, tube                                        # noqa: E402

from glyphs import text_polys, check_coverage                          # noqa: E402
from lattices import (LATTICES, dedup_lengths, dedup_nodes,             # noqa: E402
                      frame, node_sockets, supercell)
import solids as S                                                     # noqa: E402

# ---------------------------------------------------------------- parameters

SOLID_LONGEST = 60.0        # longest cell edge of the solid comparison set, mm
KIT_LONGEST = 90.0          # longest cell edge of the ball-and-stick kit, mm
MONO_LONGEST = 60.0         # longest cell edge of the one-piece fused frames, mm
WS_A = 60.0                 # cubic lattice constant for the Wigner-Seitz cells, mm

MONO_HUB_R = 4.90           # Ø9.8 hubs, Ø6 struts (hubs 30% down from the original Ø14)
MONO_STRUT_R = 3.00
#: How far short of the hub centre a fused strut stops. Running struts all the way to the
#: centre looks equivalent but is not: in every centred lattice the pair corner -> centring
#: point -> opposite corner is *collinear*, so the two struts' end caps land exactly
#: coincident at the shared hub. Merging duplicate vertices then leaves each cap edge used
#: four times -- non-manifold, and trimesh reports the body as not watertight. Stopping
#: 4.2 mm short keeps 2.8 mm of overlap inside the Ø14 hub, so the union is still solid and
#: collinear struts still read as one continuous rod through the ball.
MONO_STRUT_INSET = 0.6 * MONO_HUB_R

SUPER_N = 2                 # cells per side of the repeating block; 3 busts the bed for hR
SUPER_LONGEST = 50.0        # longest edge of ONE cell inside the block, mm
SUPER_HUB_R = 3.85          # Ø7.7 hubs, Ø5 struts (hubs 30% down from the original Ø11)
SUPER_STRUT_R = 2.5
SUPER_HIGHLIGHT = 1.35      # the origin cell is drawn this much heavier than the rest


def _hub_hides_cap(hub_r, strut_r, inset, name):
    """A strut's end cap has to finish *inside* its hub.

    The cap sits `inset` from the node centre with radius `strut_r`, so its rim is
    hypot(inset, strut_r) from the centre -- not `inset`. Shrinking a hub without checking
    that leaves the flat disc poking through the ball as a visible ring, and the two radii
    are independent parameters, so it is easy to do by accident. Enforced rather than
    eyeballed because it survives future edits to these numbers."""
    rim = math.hypot(inset, strut_r)
    if rim > hub_r - 0.15:
        raise ValueError(
            f"{name}: strut end cap rim reaches {rim:.2f} mm but the hub radius is only "
            f"{hub_r:.2f} mm -- the cap would show as a ring. Raise the hub radius or "
            f"lower the strut radius / inset.")
    return hub_r - rim

HUB_R = 7.0                 # node sphere radius -> Ø14 nodes
SLEEVE_R0 = 6.0             # sleeve starts inside the hub, so the hub caps the bore
SLEEVE_R1 = 13.0            # ... giving a socket 13 - 7 = 6 mm deep
BORE_R = 3.15               # Ø6.30 bore for a Ø6.00 rod: 0.30 mm PETG clearance
SLEEVE_OUTER_R = 4.60       # Ø9.20 sleeve, 1.45 mm wall
STRUT_R = 3.00              # Ø6.00 rods

CAP_MAX = 10.0              # embossed Pearson symbol cap height, mm
EMBOSS = 1.2                # emboss relief, mm
LABEL_MARGIN = 3.0          # clear space around the label on its face, mm

BED = 320.0                 # printer envelope, mm (PETG, 320 x 320 x 320)

CORE_SEVEN = ["aP", "mP", "oP", "tP", "hP", "hR", "cP"]     # one per crystal system
CENTRED_SEVEN = ["mS", "oS", "oI", "oF", "tI", "cI", "cF"]  # completes the 14

HERE = os.path.dirname(os.path.abspath(__file__))
STL = os.path.join(HERE, "stl")


# ---------------------------------------------------------------- labelling

def emboss(mesh, sym, origin, ex, ey, outward, outline, tag="label"):
    """Emboss `sym` on a convex planar face, as large as it will fit.

    `outline` is the face boundary in local (ex, ey) coordinates relative to `origin`. The
    label is centred on the face centroid and the cap height is found by shrinking until
    every glyph vertex clears the face boundary by LABEL_MARGIN. Testing the real boundary
    rather than a bounding box matters for the sheared cells -- turquoise's a-c face leans
    about 24 degrees, so a box-fitted label would hang off the low corner -- and it lets the
    rhombic and hexagonal Wigner-Seitz faces use the same code path."""
    outline = np.asarray(outline, float)
    centre_local = outline.mean(axis=0)
    eqs = ConvexHull(outline).equations               # a . x + b <= 0 inside

    cap = None
    for trial in np.arange(CAP_MAX, 3.75, -0.25):
        pts = np.array([p for poly in text_polys(sym, float(trial)) for p in poly])
        pts = pts + centre_local
        if float((pts @ eqs[:, :2].T + eqs[:, 2]).max()) <= -LABEL_MARGIN:
            cap = float(trial)
            break
    if cap is None:
        return None

    origin = np.asarray(origin, float)
    at = origin + centre_local[0] * np.asarray(ex) + centre_local[1] * np.asarray(ey)
    for poly in text_polys(sym, cap):
        mesh.add_solid(*S.prism_on_plane(poly, at, ex, ey, outward, EMBOSS), tag)
    return cap


# ---------------------------------------------------------------- solid cells

def build_solid(latt):
    """One solid unit cell, standing on its a-b face, label on the a-c face."""
    m = latt.vectors(longest=SOLID_LONGEST)
    A, B, C = m[0], m[1], m[2]

    mesh = Mesh()
    mesh.add_solid(*S.parallelepiped(A, B, C), "cell")

    # The a-c face at b = 0, outward normal pointing away from +B.
    n = np.cross(C, A)
    n /= np.linalg.norm(n)
    if float(n @ B) > 0.0:
        n = -n
    ex = A / np.linalg.norm(A)
    perp = C - float(C @ ex) * ex
    ey = perp / np.linalg.norm(perp)
    la, cx, cy = float(np.linalg.norm(A)), float(C @ ex), float(np.linalg.norm(perp))
    outline = [(0.0, 0.0), (la, 0.0), (la + cx, cy), (cx, cy)]
    cap = emboss(mesh, latt.pearson, np.zeros(3), ex, ey, n, outline)

    v, f = mesh.arrays()
    shift = np.array([-(v[:, 0].min() + v[:, 0].max()) / 2.0,
                      -(v[:, 1].min() + v[:, 1].max()) / 2.0, -v[:, 2].min()])
    mesh.v = [tuple(np.asarray(p) + shift) for p in mesh.v]
    return mesh, cap


# ---------------------------------------------------------------- kit parts

def build_node(dirs):
    """One hub with a sleeve per socket, rotated into its best print orientation."""
    up = S.best_up(dirs)
    R = S.align(up, np.array([0.0, 0.0, 1.0]))
    d = np.asarray(dirs, float) @ R.T

    mesh = Mesh()
    mesh.add_solid(*S.hub(HUB_R), "hub")
    for k, di in enumerate(d):
        mesh.add_solid(*S.sleeve(di, SLEEVE_R0, SLEEVE_R1, BORE_R, SLEEVE_OUTER_R),
                       f"sleeve{k}")

    v, _ = mesh.arrays()
    dz = -v[:, 2].min()
    mesh.v = [(x, y, z + dz) for x, y, z in mesh.v]
    elev = float(np.degrees(np.arcsin(np.clip(d[:, 2], -1.0, 1.0))).min())
    return mesh, elev


def build_strut(length):
    """A plain rod, lying on the bed along +x so it prints without support and its layers
    run along the axis rather than across the joint."""
    mesh = Mesh()
    mesh.add_solid(*tube((0.0, 0.0, STRUT_R), (length, 0.0, STRUT_R), STRUT_R, nseg=24),
                   "strut")
    return mesh


def build_kit(latt):
    """Deduped node and strut parts for one lattice."""
    f, pos, dirs = node_sockets(latt, KIT_LONGEST)
    classes, _ = dedup_nodes(dirs)

    lengths = [float(np.linalg.norm(pos[j] - pos[i])) - 2.0 * HUB_R for i, j in f.bonds]
    values, label = dedup_lengths(lengths)
    counts = [label.count(k) for k in range(len(values))]

    nodes = []
    for c, members in enumerate(classes):
        mesh, elev = build_node(dirs[members[0]])
        nodes.append({"name": chr(ord("A") + c), "count": len(members),
                      "sockets": len(dirs[members[0]]), "mesh": mesh, "min_elev": elev,
                      "kind": f.kind[members[0]]})
    struts = [{"length": values[k], "count": counts[k], "mesh": build_strut(values[k])}
              for k in range(len(values))]
    return nodes, struts


# ---------------------------------------------------------------- one-piece frames

def orientations(latt, longest):
    """The six ways to stand a cell on one of its own faces, as rotation matrices.

    Restricting to face-down orientations guarantees a stable print. The tempting
    alternative for a cubic frame is to stand it on a body diagonal, which puts all twelve
    edges at 35.3 degrees and needs no bridging at all -- but it also balances a 104 mm tower
    on one hub, and a print that falls over costs more than a support that has to be
    snipped."""
    m = latt.vectors(longest=longest)
    down = np.array([0.0, 0.0, -1.0])
    out = []
    for i, j, k in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):
        n = np.cross(m[i], m[j])
        n /= np.linalg.norm(n)
        side = 1.0 if float(n @ m[k]) > 0.0 else -1.0
        names = "abc"
        for t, outward in ((0, -side * n), (1, side * n)):
            out.append((f"{names[i]}{names[j]}@{t}", S.align(outward, down)))
    return out


def frame_printability(pos, bonds, bed_tol=0.6, flat_deg=6.0):
    """How a one-piece frame prints in a given orientation.

    Two failure modes matter and they are not the same thing:

      bridges  horizontal struts with air underneath. Harmless at the bottom level, where
               they sag a fraction of a millimetre onto the bed; a real bridge higher up.
      islands  hubs that no strut reaches from below. These start in mid-air and *must*
               have support -- every face-centred and base-centred lattice has one, the
               centring hub on the top face, reached only by struts in its own horizontal
               plane."""
    z = pos[:, 2]
    zmin = z.min()
    bridges, span = 0, 0.0
    for i, j in bonds:
        d = pos[j] - pos[i]
        L = float(np.linalg.norm(d))
        elev = math.degrees(math.asin(abs(float(d[2])) / L))
        if elev < flat_deg and min(z[i], z[j]) > zmin + bed_tol:
            bridges += 1
            # The free span is hub-centre spacing less the two hubs the bridge anchors
            # into -- that gap is what the slicer actually has to cross.
            span = max(span, L - 2.0 * MONO_HUB_R)

    lower = {k: False for k in range(len(pos))}
    for i, j in bonds:
        if z[i] < z[j] - bed_tol:
            lower[j] = True
        elif z[j] < z[i] - bed_tol:
            lower[i] = True
    islands = sum(1 for k in range(len(pos))
                  if not lower[k] and z[k] > zmin + bed_tol)

    footprint = float((z < zmin + bed_tol).sum())
    return {"bridges": bridges, "islands": islands, "span": span, "feet": footprint}


def build_assembled(latt):
    """One printed piece per lattice: a hub at every lattice point, a strut along every
    bond, fused by overlap. No sleeves, no bores, no assembly.

    Scaled to MONO_LONGEST rather than the kit's 90 mm because span is what limits this
    version: at 90 mm a cubic cell asks PETG to bridge 76 mm of Ø6 rod in mid-air, at 60 mm
    it asks for 46 mm. It also puts these at the same scale as the solid comparison set, so
    a frame and its solid sit side by side."""
    f, pos0, _ = node_sockets(latt, MONO_LONGEST)
    best = None
    for name, R in orientations(latt, MONO_LONGEST):
        p = np.asarray(pos0) @ R.T
        rep = frame_printability(p, f.bonds)
        key = (rep["islands"], rep["bridges"], -rep["feet"], rep["span"])
        if best is None or key < best[0]:
            best = (key, name, p, rep)
    _, face, p, rep = best

    p = p - np.array([(p[:, 0].min() + p[:, 0].max()) / 2.0,
                      (p[:, 1].min() + p[:, 1].max()) / 2.0,
                      p[:, 2].min() - MONO_HUB_R])
    zmin = p[:, 2].min()

    mesh = Mesh()
    for k, q in enumerate(p):
        if q[2] < zmin + 0.6:
            # On the bed: truncated hub, so the frame stands on flats instead of points.
            v, fc = S.hub(MONO_HUB_R)
            v = np.asarray(v) + q
            mesh.add_solid(v, fc, "hub")
        else:
            v, fc = S.sphere(MONO_HUB_R)
            mesh.add_solid(np.asarray(v) + q, fc, "hub")
    for i, j in f.bonds:
        d = p[j] - p[i]
        d = d / np.linalg.norm(d)
        mesh.add_solid(*tube(p[i] + MONO_STRUT_INSET * d, p[j] - MONO_STRUT_INSET * d,
                             MONO_STRUT_R, nseg=16), "strut")

    v, _ = mesh.arrays()
    dz = -v[:, 2].min()
    mesh.v = [(x, y, z + dz) for x, y, z in mesh.v]
    rep["face"] = face
    return mesh, rep


# ---------------------------------------------------------------- supercells

def build_supercell(latt, n=None):
    """An n x n x n block of fused cells: the lattice actually repeating.

    The unit cell at the origin is drawn SUPER_HIGHLIGHT times heavier than the rest, so the
    repeating unit is visible inside the block. Without it the object is a thicket; with it
    you can see one cell and then see it recur, which is the point of the model.

    A single cell cannot show what this does. Here an interior corner is one hub belonging to
    eight cells simultaneously, a face-centring hub is shared by two, and a body-centring hub
    is shared by none -- which is the geometric content of "1, 2 and 4 lattice points per
    conventional cell" for P, I and F. In a body-centred block the centring hubs and the
    corner hubs also end up with identical surroundings, which is the sense in which a
    centred lattice's "extra" point is not extra at all."""
    n = SUPER_N if n is None else n
    f, origin_nodes, origin_bonds = supercell(latt, n)
    m = latt.vectors(longest=SUPER_LONGEST)
    pos0 = np.array([np.asarray(fr, float) @ m for fr in f.frac])

    best = None
    for name, R in orientations(latt, SUPER_LONGEST):
        p = pos0 @ R.T
        rep = frame_printability(p, f.bonds, bed_tol=0.6 * SUPER_HUB_R)
        key = (rep["islands"], rep["bridges"], -rep["feet"], rep["span"])
        if best is None or key < best[0]:
            best = (key, name, p, rep)
    _, face, p, rep = best

    p = p - np.array([(p[:, 0].min() + p[:, 0].max()) / 2.0,
                      (p[:, 1].min() + p[:, 1].max()) / 2.0,
                      p[:, 2].min() - SUPER_HUB_R * SUPER_HIGHLIGHT])
    zfloor = min(p[k, 2] - (SUPER_HUB_R * (SUPER_HIGHLIGHT if k in origin_nodes else 1.0))
                 for k in range(len(p)))

    mesh = Mesh()
    for k, q in enumerate(p):
        r = SUPER_HUB_R * (SUPER_HIGHLIGHT if k in origin_nodes else 1.0)
        tag = "cell" if k in origin_nodes else "hub"
        if q[2] - r < zfloor + 0.6:
            mesh.add_solid(np.asarray(S.hub(r)[0]) + q, S.hub(r)[1], tag)
        else:
            v, fc = S.sphere(r, nu=20, nv=12)
            mesh.add_solid(np.asarray(v) + q, fc, tag)
    for i, j in f.bonds:
        heavy = (i, j) in origin_bonds
        r = SUPER_STRUT_R * (SUPER_HIGHLIGHT if heavy else 1.0)
        d = p[j] - p[i]
        d = d / np.linalg.norm(d)
        inset = 0.6 * SUPER_HUB_R
        mesh.add_solid(*tube(p[i] + inset * d, p[j] - inset * d, r, nseg=14),
                       "cell" if heavy else "strut")

    v, _ = mesh.arrays()
    dz = -v[:, 2].min()
    mesh.v = [(x, y, z + dz) for x, y, z in mesh.v]
    rep.update(face=face, cells=n ** 3, nodes=len(f.frac), struts=len(f.bonds))
    return mesh, rep


# ---------------------------------------------------------------- Wigner-Seitz

def wigner_seitz(a=WS_A):
    """Cube, truncated octahedron and rhombic dodecahedron -- the Wigner-Seitz cells of
    cP, cI and cF, all at cubic lattice constant `a`.

    The volume of each is the primitive-cell volume of its lattice, so they come out as
    a^3, a^3/2 and a^3/4; `main` asserts that, which is a sharp check on the vertex sets.
    Reciprocal space swaps the last two: cI's Brillouin zone is the rhombic dodecahedron
    and cF's is the truncated octahedron."""
    h = a / 2.0
    q = a / 4.0

    cube = [(x, y, z) for x in (-h, h) for y in (-h, h) for z in (-h, h)]

    # All 24 *permutations* of (0, +-a/4, +-a/2). Assigning a/4 to the lower-indexed free
    # axis every time yields only 12 points, which hull to a different solid entirely --
    # the volume assertion below is what catches that.
    trunc_oct = set()
    for perm in itertools.permutations(range(3)):
        for s1 in (-1, 1):
            for s2 in (-1, 1):
                p = [0.0, 0.0, 0.0]
                p[perm[1]] = s1 * q
                p[perm[2]] = s2 * h
                trunc_oct.add(tuple(p))
    trunc_oct = sorted(trunc_oct)
    assert len(trunc_oct) == 24, len(trunc_oct)

    rhombic_dodec = [(s * h, 0.0, 0.0) for s in (-1, 1)] \
        + [(0.0, s * h, 0.0) for s in (-1, 1)] \
        + [(0.0, 0.0, s * h) for s in (-1, 1)] \
        + [(x, y, z) for x in (-q, q) for y in (-q, q) for z in (-q, q)]

    return [
        ("cP", "cube", cube, a ** 3, None),
        ("cI", "truncated_octahedron", trunc_oct, a ** 3 / 2.0, None),
        # Rest the rhombic dodecahedron on a face: its face normals are the <110>
        # directions, so bring (1,1,0) down to -z.
        ("cF", "rhombic_dodecahedron", rhombic_dodec, a ** 3 / 4.0,
         S.align(np.array([1.0, 1.0, 0.0]) / math.sqrt(2.0), np.array([0.0, 0.0, -1.0]))),
    ]


def build_ws(pearson, points, rotation):
    v, f = S.hull_solid(points)
    if rotation is not None:
        v = v @ rotation.T
        v, f = S.orient_outward(v, f)
    v = v - np.array([(v[:, 0].min() + v[:, 0].max()) / 2.0,
                      (v[:, 1].min() + v[:, 1].max()) / 2.0, v[:, 2].min()])

    mesh = Mesh()
    mesh.add_solid(v, f, "cell")
    centre, ex, ey, normal, outline = S.face_frame(v, f, target=(0.0, 0.0, 1.0))
    cap = emboss(mesh, pearson, centre, ex, ey, normal, outline)
    return mesh, cap


# ---------------------------------------------------------------- reporting

def audit(mesh, name, problems):
    """Per-solid edge-manifoldness plus the whole-part print metrics.

    Overhang is measured on the body only. Every embossed glyph contributes horizontal
    downward ledges of stroke-width x emboss-depth -- a 1.6 x 1.2 mm shelf -- which any
    slicer bridges without a thought, exactly as it does for embossed text on any print.
    Counting those as overhangs reported 0 degrees and ~60 mm2 of 'support needed' for
    every cell in the set and buried the two cells that genuinely need it."""
    bad = mesh.validate()
    if bad:
        problems.append(f"{name}: {len(bad)} non-manifold solid(s)")
    v, f = mesh.arrays()
    ext = v.max(axis=0) - v.min(axis=0)
    if ext.max() > BED:
        problems.append(f"{name}: {ext.max():.1f} mm exceeds the {BED:.0f} mm bed")
    vb, fb = mesh.pick(drop_tags=["label"])
    slope, area = S.overhang_slope(vb, fb)
    return {"name": name, "tris": len(f), "extents": ext, "slope": slope,
            "support_area": area}


def clean_output():
    """Delete STLs from a previous run before writing new ones.

    Part *identity* changes when the geometry changes: the first build of this set gave
    natrolite 14 node types, the current spider gives it 5, and the files for types F..N
    stayed on disk describing a design that no longer exists. A stale node in the output
    directory is a part you can print by mistake and then fail to fit, so the build starts
    from an empty tree rather than overwriting selectively."""
    removed = 0
    for root, dirs, files in os.walk(STL, topdown=False):
        for f in files:
            if f.endswith(".stl"):
                os.remove(os.path.join(root, f))
                removed += 1
        for d in dirs:
            p = os.path.join(root, d)
            try:
                if not os.listdir(p):
                    os.rmdir(p)
            except OSError:
                pass          # a viewer or indexer holding the directory is not our problem
    return removed


def main():
    check_coverage([l.pearson for l in LATTICES])
    print(f"hub clearance over strut cap: "
          f"frames {_hub_hides_cap(MONO_HUB_R, MONO_STRUT_R, MONO_STRUT_INSET, 'frames'):.2f} mm, "
          f"blocks {_hub_hides_cap(SUPER_HUB_R, SUPER_STRUT_R, 0.6 * SUPER_HUB_R, 'blocks'):.2f} mm")
    stale = clean_output()
    if stale:
        print(f"removed {stale} STLs from the previous run\n")
    for d in ("solids", "kit", "assembled", "supercell", "wigner_seitz"):
        os.makedirs(os.path.join(STL, d), exist_ok=True)
    problems = []

    # ---- solid comparison set
    print(f"{'sym':4} {'mineral':14} {'extents (mm)':>26} {'vol':>8} "
          f"{'overhang':>9} {'support':>9} {'cap':>5}")
    solid_rows = []
    for l in LATTICES:
        mesh, cap = build_solid(l)
        path = os.path.join(STL, "solids", f"{l.pearson}_{l.mineral}_solid_60mm.stl")
        S.write_stl(path, mesh, f"{l.pearson} {l.mineral} unit cell - solid 60mm")
        r = audit(mesh, f"solids/{l.pearson}", problems)
        v, f = mesh.arrays()
        vol = abs(float(np.linalg.det(l.vectors(longest=SOLID_LONGEST)))) / 1000.0
        e = r["extents"]
        print(f"{l.pearson:4} {l.mineral:14} {e[0]:7.1f} x{e[1]:7.1f} x{e[2]:7.1f} "
              f"{vol:8.1f} {r['slope']:8.1f}d {r['support_area']:8.0f}mm2 "
              f"{(cap or 0):5.1f}")
        solid_rows.append((l, r, vol, cap))

    # ---- ball-and-stick kit
    print()
    kit = {}
    tot_nodes = tot_struts = 0
    for l in LATTICES:
        nodes, struts = build_kit(l)
        d = os.path.join(STL, "kit", f"{l.pearson}_{l.mineral}")
        os.makedirs(d, exist_ok=True)
        for nd in nodes:
            S.write_stl(os.path.join(d, f"node_{nd['name']}_{nd['sockets']}way_x{nd['count']}.stl"),
                        nd["mesh"], f"{l.pearson} node {nd['name']} - {nd['sockets']} sockets")
            audit(nd["mesh"], f"kit/{l.pearson}/node_{nd['name']}", problems)
        for st in struts:
            S.write_stl(os.path.join(d, f"strut_{st['length']:.1f}mm_x{st['count']}.stl"),
                        st["mesh"], f"{l.pearson} strut {st['length']:.1f} mm")
            audit(st["mesh"], f"kit/{l.pearson}/strut_{st['length']:.1f}", problems)
        n_parts = sum(n["count"] for n in nodes)
        s_parts = sum(s["count"] for s in struts)
        tot_nodes += n_parts
        tot_struts += s_parts
        worst = min(n["min_elev"] for n in nodes)
        kit[l.pearson] = (nodes, struts)
        print(f"{l.pearson:4} {l.mineral:14} {len(nodes)} node types -> {n_parts:3} parts, "
              f"{len(struts)} strut lengths -> {s_parts:3} parts   "
              f"worst sleeve elevation {worst:+6.1f}d")
    print(f"     {'TOTAL':14} {tot_nodes} nodes + {tot_struts} struts = "
          f"{tot_nodes + tot_struts} parts")

    # ---- one-piece fused frames
    print()
    print(f"{'sym':4} {'mineral':14} {'extents (mm)':>26} {'face':>7} "
          f"{'feet':>5} {'bridges':>8} {'span':>7} {'islands':>8}")
    mono_rows = []
    for l in LATTICES:
        mesh, rep = build_assembled(l)
        path = os.path.join(STL, "assembled",
                            f"{l.pearson}_{l.mineral}_frame_{MONO_LONGEST:.0f}mm.stl")
        S.write_stl(path, mesh, f"{l.pearson} {l.mineral} one-piece cell frame")
        r = audit(mesh, f"assembled/{l.pearson}", problems)
        e = r["extents"]
        print(f"{l.pearson:4} {l.mineral:14} {e[0]:7.1f} x{e[1]:7.1f} x{e[2]:7.1f} "
              f"{rep['face']:>7} {rep['feet']:5.0f} {rep['bridges']:8} "
              f"{rep['span']:6.1f}mm {rep['islands']:8}")
        mono_rows.append((l, r, rep))

    # ---- repeating blocks
    print()
    print(f"{'sym':4} {'mineral':14} {'extents (mm)':>26} {'cells':>6} {'hubs':>5} "
          f"{'struts':>7} {'face':>7} {'bridges':>8} {'islands':>8}")
    super_rows = []
    for l in LATTICES:
        mesh, rep = build_supercell(l)
        path = os.path.join(STL, "supercell",
                            f"{l.pearson}_{l.mineral}_{SUPER_N}x{SUPER_N}x{SUPER_N}"
                            f"_{SUPER_LONGEST:.0f}mm.stl")
        S.write_stl(path, mesh, f"{l.pearson} {l.mineral} {SUPER_N}x{SUPER_N}x{SUPER_N} block")
        r = audit(mesh, f"supercell/{l.pearson}", problems)
        e = r["extents"]
        print(f"{l.pearson:4} {l.mineral:14} {e[0]:7.1f} x{e[1]:7.1f} x{e[2]:7.1f} "
              f"{rep['cells']:6} {rep['nodes']:5} {rep['struts']:7} {rep['face']:>7} "
              f"{rep['bridges']:8} {rep['islands']:8}")
        super_rows.append((l, r, rep))

    # ---- Wigner-Seitz cells
    print()
    ws_rows = []
    for pearson, name, pts, expect_vol, rot in wigner_seitz():
        mesh, cap = build_ws(pearson, pts, rot)
        v, f = mesh.pick(tags=["cell"])
        got = _mesh_volume(v, f)
        assert abs(got - expect_vol) / expect_vol < 1e-6, \
            f"{name}: volume {got:.1f} != {expect_vol:.1f}"
        path = os.path.join(STL, "wigner_seitz",
                            f"{pearson}_wigner_seitz_{name}_{WS_A:.0f}mm.stl")
        S.write_stl(path, mesh, f"{pearson} Wigner-Seitz cell - {name}")
        r = audit(mesh, f"wigner_seitz/{name}", problems)
        e = r["extents"]
        print(f"{pearson:4} {name:24} {e[0]:6.1f} x{e[1]:6.1f} x{e[2]:6.1f}  "
              f"vol {got / 1000.0:7.1f} cm3 = a^3/{a_ratio(expect_vol):.0f}  "
              f"overhang {r['slope']:5.1f}d  cap {cap:4.1f}")
        ws_rows.append((pearson, name, r, got, cap))

    write_manifest(solid_rows, kit, mono_rows, super_rows, ws_rows)

    print()
    if problems:
        print(f"PROBLEMS ({len(problems)}):")
        for p in problems:
            print("  " + p)
    else:
        print("all solids edge-manifold, all parts inside the 320 mm bed")
    return 1 if problems else 0


def a_ratio(vol):
    return WS_A ** 3 / vol


def _mesh_volume(v, f):
    tri = v[f]
    return abs(float(np.einsum("ij,ij->", tri[:, 0],
                              np.cross(tri[:, 1], tri[:, 2])) / 6.0))


def write_manifest(solid_rows, kit, mono_rows, super_rows, ws_rows):
    lines = ["# Crystal lattice parts", "",
             "Generated by `build.py`. Cell constants and their sources are documented in",
             "`lattices.py`; print settings and assembly are in `README.md`.", ""]

    lines += ["## Solid comparison set (one 320 mm plate)", "",
              "Longest cell edge 60 mm, resting on the a-b face, Pearson symbol embossed",
              f"{EMBOSS} mm on the a-c face.", "",
              "| Pearson | mineral | bounding box (mm) | cell volume (cm3) | "
              "shallowest overhang | support |", "|---|---|---|---|---|---|"]
    for l, r, vol, cap in solid_rows:
        e = r["extents"]
        need = "none" if r["slope"] >= 45.0 else f"{r['support_area']:.0f} mm2"
        lines.append(f"| {l.pearson} | {l.mineral} | {e[0]:.1f} x {e[1]:.1f} x {e[2]:.1f} "
                     f"| {vol:.1f} | {r['slope']:.1f}d | {need} |")

    lines += ["", "## One-piece fused frames (no assembly)", "",
              f"Longest cell edge {MONO_LONGEST:.0f} mm, hubs and struts fused into a single",
              "print. `face` is which cell face stands on the bed and `feet` how many hubs",
              "touch it. `bridges` counts horizontal struts with air underneath, above the",
              "bottom level; `islands` counts hubs no strut reaches from below, which are the",
              "ones that genuinely require support.", "",
              "| Pearson | mineral | bounding box (mm) | face | feet | bridges | longest span | islands |",
              "|---|---|---|---|---|---|---|---|"]
    for l, r, rep in mono_rows:
        e = r["extents"]
        lines.append(f"| {l.pearson} | {l.mineral} | {e[0]:.1f} x {e[1]:.1f} x {e[2]:.1f} "
                     f"| {rep['face']} | {rep['feet']:.0f} | {rep['bridges']} "
                     f"| {rep['span']:.1f} mm | {rep['islands']} |")

    lines += ["", f"## Repeating blocks ({SUPER_N}x{SUPER_N}x{SUPER_N} cells, no assembly)",
              "",
              f"One cell edge is {SUPER_LONGEST:.0f} mm, so each block is {SUPER_N} cells across.",
              "Corners and struts shared between neighbouring cells are merged, so `hubs` and",
              f"`struts` are far fewer than {SUPER_N ** 3} copies of one cell. The cell at the",
              f"origin is drawn {SUPER_HIGHLIGHT:g}x heavier so the repeating unit is visible",
              "inside the block.", "",
              "| Pearson | mineral | bounding box (mm) | cells | hubs | struts | face | bridges | islands |",
              "|---|---|---|---|---|---|---|---|---|"]
    for l, r, rep in super_rows:
        e = r["extents"]
        lines.append(f"| {l.pearson} | {l.mineral} | {e[0]:.1f} x {e[1]:.1f} x {e[2]:.1f} "
                     f"| {rep['cells']} | {rep['nodes']} | {rep['struts']} | {rep['face']} "
                     f"| {rep['bridges']} | {rep['islands']} |")

    lines += ["", "## Wigner-Seitz cells", "",
              f"Shared cubic lattice constant a = {WS_A:.0f} mm. The volumes are the",
              "primitive-cell volumes, so they are exactly a^3, a^3/2 and a^3/4.", "",
              "| lattice | solid | volume (cm3) | as a fraction of a^3 |",
              "|---|---|---|---|"]
    for pearson, name, r, vol, cap in ws_rows:
        lines.append(f"| {pearson} | {name.replace('_', ' ')} | {vol / 1000.0:.1f} "
                     f"| 1/{a_ratio(vol):.0f} |")

    for title, group in (("Plate 1 - the seven crystal systems", CORE_SEVEN),
                         ("Plate 2 - the seven centred lattices", CENTRED_SEVEN)):
        n = sum(sum(x["count"] for x in kit[p][0]) for p in group)
        s = sum(sum(x["count"] for x in kit[p][1]) for p in group)
        lines += ["", f"## {title}", "", f"{n} nodes + {s} struts = {n + s} parts.", ""]
        for p in group:
            nodes, struts = kit[p]
            latt = next(l for l in LATTICES if l.pearson == p)
            lines += [f"### {p} - {latt.mineral}", "",
                      "| part | sockets / length | quantity |", "|---|---|---|"]
            for nd in nodes:
                lines.append(f"| node {nd['name']} | {nd['sockets']} sockets "
                             f"| {nd['count']} |")
            for st in struts:
                lines.append(f"| strut | {st['length']:.1f} mm | {st['count']} |")
            lines.append("")

    with open(os.path.join(HERE, "PARTS.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
