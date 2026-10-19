#!/usr/bin/env python3
"""Check Hercules topography counts, spacing, base and perimeter relief (NumPy)."""

# USAGE (Python 3 + NumPy; run this example from the repository root)
#   python3 script/validate_topography.py inputfiles/topography/topography.in \
#       --lx 40000 --ly 40000 --base 1250
#   python3 script/validate_topography.py --help
#
# Required: input filename, --lx (north extent), --ly (east extent), and --base
# (the run's topographybase_zcoord), all in metres. Relative paths use the
# working directory. File format: positive uniform spacing followed by nodal
# elevations, positive up, X outer/Y inner; spacing must divide both lengths.
# Checks: finite values, exact point count, positive lengths/base/spacing, and
# a base above the highest elevation. --base is a datum inside the domain;
# it does not add depth to region_depth_deep_m.
# Output: grid dimensions, elevation range, and perimeter elevation range.
# No files are written; all elevations are loaded into memory.
# Non-flat boundaries are reported, not rejected. A flat perimeter alone does
# not prove the entire CPML collar is flat. This does not validate seabed
# consistency, material coverage, or source clearance. See script/README.md.

import argparse
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", type=Path, help="topography.in")
    parser.add_argument("--lx", type=float, required=True, help="North extent, metres")
    parser.add_argument("--ly", type=float, required=True, help="East extent, metres")
    parser.add_argument("--base", type=float, required=True, help="topographybase_zcoord, metres")
    args = parser.parse_args()
    if not np.isfinite([args.lx, args.ly, args.base]).all() or min(args.lx, args.ly, args.base) <= 0:
        parser.error("domain lengths and base must be finite and positive")
    raw = np.array(args.file.read_text().split(), dtype=float)
    if raw.size < 2 or not np.isfinite(raw).all() or raw[0] <= 0:
        parser.error("need a positive spacing followed by finite elevations")
    spacing, values = raw[0], raw[1:]
    nx, ny = (round(length / spacing) + 1 for length in (args.lx, args.ly))
    if (nx < 2 or ny < 2
            or not np.allclose([(nx - 1) * spacing, (ny - 1) * spacing],
                               [args.lx, args.ly], rtol=1e-12, atol=1e-8)):
        parser.error("topographic spacing does not divide the domain")
    if values.size != nx * ny:
        parser.error(f"expected {nx * ny} elevations, found {values.size}")
    elevation = values.reshape(nx, ny)
    if args.base <= elevation.max():
        parser.error("topographybase_zcoord must exceed the highest elevation")
    perimeter = np.r_[elevation[0, :], elevation[-1, :], elevation[:, 0], elevation[:, -1]]
    print(f"grid={nx}x{ny}, range={elevation.min():.3f}..{elevation.max():.3f} m")
    print(f"perimeter range={np.ptp(perimeter):.6g} m")


if __name__ == "__main__":
    main()
