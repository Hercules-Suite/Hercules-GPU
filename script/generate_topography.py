#!/usr/bin/env python3
"""Generate the documentation's Gaussian hill in solver X-major order (NumPy)."""

# USAGE (Python 3 + NumPy; run this example from the repository root)
#   python3 script/generate_topography.py \
#       --lx 40000 --ly 40000 --spacing 100 --height 800 --sigma 3000 \
#       --out inputfiles/topography/topography.in
#   python3 script/generate_topography.py --help
#
# Required: --lx (north extent), --ly (east extent), --spacing (both axes),
# --height (peak elevation), --sigma (Gaussian width), and --out. All dimensions
# and elevations are metres. Lengths/spacing/sigma must be positive and spacing
# must divide both lengths. Heights are positive UP relative to the datum.
# Optional: --center-x and --center-y default to the domain midpoint; --force
# permits overwriting an existing output (otherwise it is rejected).
# Output: topography.in format, one spacing value followed by all nodal heights
# in X-outer/Y-inner order (east fastest), written to six decimal places.
# Parent directories are created; relative paths use the working directory.
# Arrays are held in memory. This creates only an illustrative Gaussian hill;
# it does not taper the boundary flat, create seabed/material files, or choose
# topographybase_zcoord. Configure those separately and regenerate the mesh.
# See script/README.md for validation and boundary-taper examples.

import argparse
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lx", type=float, required=True, help="North extent, metres")
    parser.add_argument("--ly", type=float, required=True, help="East extent, metres")
    parser.add_argument("--spacing", type=float, required=True, help="Common X/Y spacing, metres")
    parser.add_argument("--height", type=float, required=True, help="Peak elevation, metres")
    parser.add_argument("--sigma", type=float, required=True, help="Gaussian width, metres")
    parser.add_argument("--center-x", type=float, help="North coordinate; default domain midpoint")
    parser.add_argument("--center-y", type=float, help="East coordinate; default domain midpoint")
    parser.add_argument("--out", type=Path, required=True, help="Output topography.in")
    parser.add_argument("--force", action="store_true", help="Replace an existing output file")
    args = parser.parse_args()
    for name in ("lx", "ly", "spacing", "sigma"):
        value = getattr(args, name)
        if not np.isfinite(value) or value <= 0:
            parser.error(f"--{name} must be finite and positive")
    center_x = args.lx / 2 if args.center_x is None else args.center_x
    center_y = args.ly / 2 if args.center_y is None else args.center_y
    if not np.isfinite([args.height, center_x, center_y]).all():
        parser.error("height and centre coordinates must be finite")
    nx, ny = (round(length / args.spacing) + 1 for length in (args.lx, args.ly))
    if (nx < 2 or ny < 2
            or not np.allclose([(nx - 1) * args.spacing, (ny - 1) * args.spacing],
                               [args.lx, args.ly], rtol=1e-12, atol=1e-8)):
        parser.error("spacing must divide both domain lengths exactly")
    x = np.arange(nx) * args.spacing
    y = np.arange(ny) * args.spacing
    north, east = np.meshgrid(x, y, indexing="ij")
    elevation = args.height * np.exp(
        -((north - center_x) ** 2 + (east - center_y) ** 2) / (2 * args.sigma ** 2))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w" if args.force else "x") as output:
        output.write(f"{args.spacing:g}\n")
        np.savetxt(output, elevation.reshape(-1, order="C"), fmt="%.6f")
    print(f"wrote {nx} x {ny} = {elevation.size} elevations to {args.out}")


if __name__ == "__main__":
    main()
