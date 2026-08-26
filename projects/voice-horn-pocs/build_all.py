#!/usr/bin/env python3
"""Build every printable part of the three voice-horn POCs and report the set.

    python build_all.py              plain smooth bodies (the product)
    python build_all.py --skin       add the ornament pass: claws, grip hand
                                     prints, braille (kept for the doc's Skin
                                     rules; rejected on looks for the print set)

Each row is measured on the EXPORTED file, not the in-memory mesh: footprint is
the minimum-area rectangle round the XY hull (parts may rotate on the plate),
and the mass estimate assumes the F1 sandwich (two 0.8 mm perimeters over 40%
Grid infill = 59% of the solid wall volume) in PETG at 1.27 g/cm3.
"""
import argparse

import clarion
import halo
import horn_lib as H
import volute

PETG, SANDWICH = 1.27, 0.592


def main(skin=True):
    rows = []
    print("=" * 78)
    print("POC-A CLARION")
    print("=" * 78)
    rows += clarion.build_projector(skin=skin)
    rows += clarion.build_squillo(skin=skin)
    rows.append(H.report("fleece_coupon", clarion.fleece_coupon(),
                         clarion.OUT + "fleece_coupon.stl"))
    print()
    print("=" * 78)
    print("POC-C HALO")
    print("=" * 78)
    rows += halo.build(skin=skin)
    print()
    print("=" * 78)
    print("POC-B VOLUTE")
    print("=" * 78)
    rows += volute.build(skin=skin)

    print()
    H.print_table(rows)
    vol = sum(r["vol_cm3"] for r in rows)
    print(f"\n{len(rows)} parts, {sum(r['faces'] for r in rows):,} triangles, "
          f"{vol:.0f} cm3 of solid wall")
    print(f"estimated PETG: {vol * PETG * SANDWICH / 1000:.2f} kg sliced as the F1 "
          f"sandwich ({vol * PETG / 1000:.2f} kg if printed solid)")
    for tag, name in (("clarion_p", "  Clarion projector"), ("clarion_s", "  Clarion squillo"),
                      ("halo", "  Halo"), ("volute", "  Volute")):
        v = sum(r["vol_cm3"] for r in rows if r["part"].startswith(tag))
        print(f"{name:<22} {v:7.0f} cm3 -> {v * PETG * SANDWICH / 1000:.2f} kg")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--skin", action="store_true",
                    help="add the ornament pass (claws, grip prints, braille); "
                         "default is the plain smooth bodies")
    main(skin=ap.parse_args().skin)
