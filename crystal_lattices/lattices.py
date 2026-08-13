#!/usr/bin/env python3
"""The 14 Bravais lattices, parameterised from verified mineral lattice constants.

Parameters
----------
Each of the 7 crystal systems is a set of constraints on the cell edges a, b, c and the
interaxial angles alpha, beta, gamma. Adding a centring type -- primitive (P), body (I),
face (F) or base (S, historically C) -- to the systems that admit it gives 14 distinct
lattices, one per row of `LATTICES`:

    1 aP + 2 mP mS + 4 oP oS oI oF + 2 tP tI + 1 hP + 1 hR + 3 cP cI cF = 14

Two counts that look like they should be larger are not, and both are worth knowing
because they explain the table rather than decorate it:

  * Monoclinic has 2 lattices, not 3, because monoclinic I and monoclinic C are the
    same lattice in different bases. Gypsum below is the demonstration: the literature
    refinement is I2/a, webmineral lists A2/a, and the conventional setting is C2/c --
    one lattice, three labels.
  * Trigonal is not a lattice type. Trigonal crystals sit on either the hexagonal
    lattice (hP -- quartz, despite the axial ratios usually quoted for "trigonal") or
    the rhombohedral one (hR -- corundum, tourmaline, calcite). Only hR has the
    a = b = c, alpha = beta = gamma /= 90 cell.

The "/=" signs in the usual textbook table are the *absence* of an imposed constraint,
not a requirement. A triclinic cell with a = b by coincidence is still triclinic.

Sources
-------
Cell constants are the conventional cell for each lattice, from
webmineral.com/data/<Mineral>.shtml, verified 2026-08-13. Three deliberate departures
from the raw source:

  topaz     webmineral gives a = 4.35 A, which is a typo -- a 7% error in the shortest
            axis. Ribbe & Gibbs, "The crystal structure of topaz and its relation to
            physical properties", Am. Mineral. 56 (1971) 24, gives a = 4.6499,
            b = 8.7968, c = 8.3909; those are used here.
  gypsum    webmineral gives the A2/a setting (a = 5.68, b = 15.18, c = 6.29). a and c
            are swapped below to reach the conventional C-centred setting, so gypsum's
            centring vector is (1/2, 1/2, 0) like every other mS in the table. The
            literature I2/a cell (a = 5.679, b = 15.202, c = 6.522, beta = 118.43) is
            the same lattice -- see the monoclinic note above.
  corundum  quoted in the hexagonal setting (a = 4.751, c = 12.97). `Lattice.cell()`
            converts it to the primitive rhombohedral cell, which is what the model
            shows; see `rhombohedral_from_hex`.

webmineral's space-group column carries typos that do not affect the lattice type
(rutile "P4/mnm" -> P4_2/mnm, beryl "P6/mmc" -> P6/mcc, turquoise "P1" -> P-1 given
that the same page lists the centrosymmetric pinacoidal class). The symbols below are
the corrected ones. Turquoise additionally has two cells in circulation; the reduced
form is a = 7.424, b = 7.629, c = 9.910, alpha = 68.61, beta = 69.71, gamma = 65.05,
the same lattice on a different axis choice.
"""

import itertools
import math
from dataclasses import dataclass, field

import numpy as np

# ---------------------------------------------------------------- parameter table

#: Fractional coordinates added to each corner by a centring type.
CENTRING = {
    "P": [],
    "I": [(0.5, 0.5, 0.5)],
    "S": [(0.5, 0.5, 0.0), (0.5, 0.5, 1.0)],
    "F": [(0.5, 0.5, 0.0), (0.5, 0.5, 1.0), (0.5, 0.0, 0.5),
          (0.5, 1.0, 0.5), (0.0, 0.5, 0.5), (1.0, 0.5, 0.5)],
}

SYSTEM_ORDER = ["triclinic", "monoclinic", "orthorhombic",
                "tetragonal", "hexagonal", "rhombohedral", "cubic"]


@dataclass(frozen=True)
class Lattice:
    """One Bravais lattice, carrying the mineral it is measured from."""

    pearson: str            # aP, mP, mS, oP, oS, oI, oF, tP, tI, hP, hR, cP, cI, cF
    system: str
    centring: str           # P, I, F, S
    mineral: str
    space_group: str
    a: float
    b: float
    c: float
    alpha: float = 90.0
    beta: float = 90.0
    gamma: float = 90.0
    hex_setting: bool = False   # True => (a, c) are hexagonal, convert to rhombohedral
    note: str = ""

    def cell(self):
        """Conventional cell constants actually used for the model, in Angstroms and
        degrees. For hR this is the primitive rhombohedral cell, not the hexagonal one."""
        if self.hex_setting:
            a_rh, alpha_rh = rhombohedral_from_hex(self.a, self.c)
            return a_rh, a_rh, a_rh, alpha_rh, alpha_rh, alpha_rh
        return self.a, self.b, self.c, self.alpha, self.beta, self.gamma

    def vectors(self, longest=None):
        """Cell vectors as rows of a 3x3 array. `longest` rescales so that the longest
        cell *edge* is that many mm -- the bounding box will exceed it when the cell
        is sheared."""
        a, b, c, al, be, ga = self.cell()
        m = basis(a, b, c, al, be, ga)
        if longest is not None:
            m = m * (longest / max(a, b, c))
        return m

    @property
    def is_centred(self):
        return self.centring != "P"


LATTICES = [
    Lattice("aP", "triclinic", "P", "turquoise", "P-1",
            7.48, 9.95, 7.68, 111.65, 115.383, 69.433,
            note="webmineral setting; a reduced cell with all-acute angles also circulates"),
    Lattice("mP", "monoclinic", "P", "azurite", "P2_1/a",
            5.008, 5.844, 10.336, beta=92.333),
    Lattice("mS", "monoclinic", "S", "gypsum", "C2/c",
            6.29, 15.18, 5.68, beta=113.833,
            note="A2/a on webmineral, I2/a in the literature -- one lattice, three settings"),
    Lattice("oP", "orthorhombic", "P", "topaz", "Pbnm",
            4.6499, 8.7968, 8.3909,
            note="a from Am. Mineral. 56 (1971) 24; webmineral's 4.35 is a typo"),
    Lattice("oS", "orthorhombic", "S", "cordierite", "Cccm",
            17.13, 9.80, 9.35),
    Lattice("oI", "orthorhombic", "I", "hemimorphite", "Imm2",
            8.37, 10.719, 5.12),
    Lattice("oF", "orthorhombic", "F", "natrolite", "Fdd2",
            18.27, 18.587, 6.56,
            note="face-centred orthorhombic is rare in minerals; natrolite is the clean case"),
    Lattice("tP", "tetragonal", "P", "rutile", "P4_2/mnm",
            4.594, 4.594, 2.958),
    Lattice("tI", "tetragonal", "I", "zircon", "I4_1/amd",
            6.604, 6.604, 5.979),
    Lattice("hP", "hexagonal", "P", "beryl", "P6/mcc",
            9.215, 9.215, 9.192, gamma=120.0,
            note="quartz shares this lattice; three of these cells assemble the hexagonal prism"),
    Lattice("hR", "rhombohedral", "P", "corundum", "R-3c",
            4.751, 4.751, 12.97, hex_setting=True,
            note="hexagonal setting converted to the primitive rhombohedral cell"),
    Lattice("cP", "cubic", "P", "pyrite", "Pa-3",
            5.417, 5.417, 5.417),
    Lattice("cI", "cubic", "I", "almandine", "Ia-3d",
            11.526, 11.526, 11.526),
    Lattice("cF", "cubic", "F", "diamond", "Fd-3m",
            3.5668, 3.5668, 3.5668),
]

assert len(LATTICES) == 14
assert len({l.pearson for l in LATTICES}) == 14


# ---------------------------------------------------------------- cell geometry

def basis(a, b, c, alpha, beta, gamma):
    """Cell vectors from the six constants, in the standard crystallographic setting:
    A along +x, B in the xy plane, C completing a right-handed set.

    The z component of C is c * V / sin(gamma) where V is the dimensionless cell-volume
    factor; it goes imaginary for angle triples that cannot close a parallelepiped, so
    this doubles as a validity check on the table."""
    al, be, ga = map(math.radians, (alpha, beta, gamma))
    ca, cb, cg, sg = math.cos(al), math.cos(be), math.cos(ga), math.sin(ga)
    vf = 1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg
    if vf <= 0.0:
        raise ValueError(f"angles {alpha}, {beta}, {gamma} do not close a cell")
    return np.array([
        [a, 0.0, 0.0],
        [b * cg, b * sg, 0.0],
        [c * cb, c * (ca - cb * cg) / sg, c * math.sqrt(vf) / sg],
    ])


def rhombohedral_from_hex(a_hex, c_hex):
    """Primitive rhombohedral cell (edge, angle in degrees) from the hexagonal setting
    of an R lattice. The hexagonal cell holds 3 rhombohedral lattice points, so the
    rhombohedron is the primitive cell of the same lattice.

        a_rh = sqrt(a^2/3 + c^2/9)
        sin(alpha_rh / 2) = 3 / (2 * sqrt(3 + (c/a)^2))
    """
    a_rh = math.sqrt(a_hex ** 2 / 3.0 + c_hex ** 2 / 9.0)
    alpha = 2.0 * math.degrees(math.asin(3.0 / (2.0 * math.sqrt(3.0 + (c_hex / a_hex) ** 2))))
    return a_rh, alpha


def cell_volume(vectors):
    return abs(float(np.linalg.det(vectors)))


# ---------------------------------------------------------------- node / strut graph

CORNERS = [(i, j, k) for i in (0, 1) for j in (0, 1) for k in (0, 1)]


@dataclass
class Frame:
    """The ball-and-stick graph of one unit cell, in fractional coordinates."""

    frac: list = field(default_factory=list)    # (x, y, z) per node
    kind: list = field(default_factory=list)    # 'corner' | 'centre'
    bonds: list = field(default_factory=list)   # (i, j) index pairs


def frame(latt):
    """Nodes and struts for one cell: the 8 corners and 12 cell edges, plus the centring
    points and a spider of diagonals holding each one.

    The spider uses *every* diagonal -- all 8 body diagonals to a body centre, all 4 face
    corners to a face or base centre -- rather than the minimum needed for rigidity. A
    tetrahedral 4-strut body spider and a 2-strut face spider would save 44 struts across
    the set, but they pick out particular corners and so break the cell's own symmetry.
    The cost is severe and not obvious up front: with a minimal spider, natrolite's 14
    lattice points need 14 *different* printed nodes, because corners stop being
    equivalent once only some of them carry a diagonal. With the full spider natrolite
    needs 3, and every corner node is interchangeable with every other.

    So: 256 struts instead of 212, and roughly a third the number of distinct node types.
    Struts are identical repeated rods, node types are things you have to tell apart
    during assembly, so this is the right way round.

    Every lattice point in a Bravais lattice is equivalent. The struts are a printing
    scaffold, not a claim that centring points differ from corners -- and the symmetric
    spider is what keeps the parts honest about that."""
    f = Frame()
    index = {}
    for cor in CORNERS:
        index[cor] = len(f.frac)
        f.frac.append(tuple(float(v) for v in cor))
        f.kind.append("corner")

    for p, q in itertools.combinations(CORNERS, 2):
        if sum(1 for u, v in zip(p, q) if u != v) == 1:
            f.bonds.append((index[p], index[q]))

    for cen in CENTRING[latt.centring]:
        ci = len(f.frac)
        f.frac.append(tuple(float(v) for v in cen))
        f.kind.append("centre")
        for nb in _spider(cen):
            f.bonds.append((ci, index[nb]))
    return f


def _spider(cen):
    """Corners a centring node is strutted to: every corner of the cell for a body centre,
    every corner of the containing face for a face or base centre. Symmetric by
    construction, which is what keeps corner nodes interchangeable."""
    fixed = [i for i, v in enumerate(cen) if abs(v - 0.5) >= 1e-9]
    return [c for c in CORNERS
            if all(abs(c[i] - cen[i]) < 1e-9 for i in fixed)]


def supercell(latt, n):
    """An n x n x n block of cells as one lattice graph, with everything shared between
    neighbours merged exactly once.

    Merging is the whole job. Tiling `frame` n^3 times and concatenating would emit a
    duplicate hub at every shared corner and a duplicate strut along every shared edge --
    coincident geometry, which welds into non-manifold shells the moment a slicer merges
    vertices. Deduplicating by position and by node-index pair is also what makes the model
    *say* something: an interior corner ends up as one node belonging to eight cells at
    once, which is why a primitive cell contains one lattice point and not eight.

    Returns (Frame, origin_nodes, origin_bonds) -- the last two identify the single cell at
    the origin, so the caller can draw it heavier and show which unit is repeating."""
    base = frame(latt)
    index, frac, kind = {}, [], []
    bonds, origin_nodes, origin_bonds = set(), set(), set()

    for i in range(n):
        for j in range(n):
            for k in range(n):
                local = []
                for fr, kd in zip(base.frac, base.kind):
                    p = (fr[0] + i, fr[1] + j, fr[2] + k)
                    key = tuple(round(v, 6) for v in p)
                    if key not in index:
                        index[key] = len(frac)
                        frac.append(p)
                        kind.append(kd)
                    local.append(index[key])
                for a, b in base.bonds:
                    e = (min(local[a], local[b]), max(local[a], local[b]))
                    bonds.add(e)
                    if (i, j, k) == (0, 0, 0):
                        origin_bonds.add(e)
                if (i, j, k) == (0, 0, 0):
                    origin_nodes.update(local)

    return Frame(frac=frac, kind=kind, bonds=sorted(bonds)), origin_nodes, origin_bonds


def node_sockets(latt, longest):
    """Per node: cartesian position and the unit socket directions leaving it."""
    m = latt.vectors(longest=longest)
    f = frame(latt)
    pos = [np.asarray(fr, float) @ m for fr in f.frac]
    dirs = [[] for _ in pos]
    for i, j in f.bonds:
        d = pos[j] - pos[i]
        d /= np.linalg.norm(d)
        dirs[i].append(d)
        dirs[j].append(-d)
    return f, [np.asarray(p, float) for p in pos], [np.asarray(d, float) for d in dirs]


# ---------------------------------------------------------------- congruence / dedup

def congruent(d1, d2, tol=1e-6):
    """True if two sets of unit socket directions are the same up to a rotation, i.e.
    the two nodes are the same printed part.

    Rather than search all permutations, build every candidate rotation from a pair of
    directions in d1 mapped to a pair in d2 -- two non-parallel vectors fix a frame --
    then test the whole set. O(n^2) candidates instead of O(n!)."""
    if len(d1) != len(d2):
        return False
    if len(d1) == 1:
        return True
    for p in range(len(d2)):
        for q in range(len(d2)):
            if p == q:
                continue
            r = _frame_rotation(d1[0], d1[1], d2[p], d2[q], tol)
            if r is None:
                continue
            moved = d1 @ r.T
            if _same_set(moved, d2, tol):
                return True
    return False


def _frame_rotation(u1, u2, v1, v2, tol):
    """Rotation carrying (u1, u2) onto (v1, v2), or None if the pairs are not congruent."""
    if abs(float(u1 @ u2) - float(v1 @ v2)) > 1e-4:
        return None
    fu, fv = _orthoframe(u1, u2), _orthoframe(v1, v2)
    if fu is None or fv is None:
        return None
    r = fv.T @ fu
    if abs(float(np.linalg.det(r)) - 1.0) > 1e-4:
        return None
    if np.linalg.norm(r @ u1 - v1) > 1e-4 or np.linalg.norm(r @ u2 - v2) > 1e-4:
        return None
    return r


def _orthoframe(u1, u2):
    e1 = u1 / np.linalg.norm(u1)
    w = u2 - float(u2 @ e1) * e1
    n = np.linalg.norm(w)
    if n < 1e-9:
        return None
    e2 = w / n
    return np.array([e1, e2, np.cross(e1, e2)])


def _same_set(x, y, tol):
    used = [False] * len(y)
    for a in x:
        for k, b in enumerate(y):
            if not used[k] and np.linalg.norm(a - b) < 1e-4:
                used[k] = True
                break
        else:
            return False
    return all(used)


def dedup_nodes(dirs):
    """Group node indices into congruence classes. Returns (classes, label_per_node)
    where classes is a list of representative index lists."""
    classes, label = [], [None] * len(dirs)
    for i, d in enumerate(dirs):
        for c, members in enumerate(classes):
            if congruent(dirs[members[0]], d):
                members.append(i)
                label[i] = c
                break
        else:
            label[i] = len(classes)
            classes.append([i])
    return classes, label


def dedup_lengths(lengths, tol=0.1):
    """Group strut lengths that a printer cannot tell apart. Returns (values, label)."""
    order = sorted(range(len(lengths)), key=lambda i: lengths[i])
    values, label = [], [None] * len(lengths)
    for i in order:
        if values and abs(lengths[i] - values[-1]) <= tol:
            label[i] = len(values) - 1
        else:
            values.append(lengths[i])
            label[i] = len(values) - 1
    return values, label


if __name__ == "__main__":
    print(f"{'sym':4} {'system':13} {'cen':4} {'mineral':14} "
          f"{'a':>7} {'b':>7} {'c':>7} {'alpha':>8} {'beta':>8} {'gamma':>8}  nodes struts")
    tot_n = tot_s = 0
    for l in LATTICES:
        a, b, c, al, be, ga = l.cell()
        f = frame(l)
        tot_n += len(f.frac)
        tot_s += len(f.bonds)
        print(f"{l.pearson:4} {l.system:13} {l.centring:4} {l.mineral:14} "
              f"{a:7.3f} {b:7.3f} {c:7.3f} {al:8.3f} {be:8.3f} {ga:8.3f}"
              f"  {len(f.frac):5} {len(f.bonds):6}")
    print(f"\ntotal kit parts: {tot_n} nodes + {tot_s} struts = {tot_n + tot_s}")
