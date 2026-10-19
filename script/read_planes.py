#!/usr/bin/env python3
"""Inspect native-float64 plane displacements with manual grid/timing (NumPy)."""

# USAGE (Python 3 + NumPy; run this example from the repository root)
#   python3 script/read_planes.py outputfiles/planes/planedisplacements.0 \
#       --nx 361 --ny 361 --dt 0.1
#   python3 script/read_planes.py --help
#
# Required: completed binary filename, --nx (along-strike POINT count), --ny
# (down-dip POINT count), and --dt (time between SAVED snapshots, seconds).
# Usually dt = simulation_delta_time_sec * output_planes_print_rate; it is not
# the solver timestep alone. --first-time defaults to 0 s; set it for restarts.
# Headerless native-float64 layout: (time, nx, ny, 3), C order, XYZ components
# contiguous, down-dip index fastest. Displacements are metres, Z positive down.
# Output: array shape, time span, and peak absolute GLOBAL Z displacement.
# All snapshots are scanned for nonfinite XYZ values through a read-only memory
# map. No output files are created; relative paths use the working directory.
# Counts/timing are not read from parameters.in. A wrong shape with the same
# product, or missing complete frames, cannot be identified from size alone.
# See script/README.md; animations_mpi.py additionally renders horizontal planes.

import argparse
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", type=Path, help="planedisplacements.N")
    parser.add_argument("--nx", type=int, required=True, help="Along-strike point count")
    parser.add_argument("--ny", type=int, required=True, help="Down-dip point count")
    parser.add_argument("--dt", type=float, required=True, help="SAVED snapshot interval, seconds")
    parser.add_argument("--first-time", type=float, default=0.0, help="First saved time, seconds")
    args = parser.parse_args()
    if args.nx < 1 or args.ny < 1 or not np.isfinite(args.dt) or args.dt <= 0:
        parser.error("grid counts and snapshot interval must be positive and finite")
    if not np.isfinite(args.first_time):
        parser.error("first-time must be finite")
    frame_bytes = args.nx * args.ny * 3 * np.dtype(np.float64).itemsize
    nt, remainder = divmod(args.file.stat().st_size, frame_bytes)
    if remainder or nt == 0:
        parser.error(f"file must contain complete {frame_bytes}-byte snapshots")
    displacement = np.memmap(args.file, dtype=np.float64, mode="r",
                             shape=(nt, args.nx, args.ny, 3), order="C")
    peak = 0.0
    for index, frame in enumerate(displacement):
        if not np.isfinite(frame).all():
            parser.error(f"snapshot {index} contains non-finite displacements")
        peak = max(peak, float(np.abs(frame[..., 2]).max()))
    print(f"layout={displacement.shape}: time, strike point, dip point, XYZ")
    print(f"time range={args.first_time:g}..{args.first_time + (nt - 1) * args.dt:g} s")
    print(f"peak |uz| = {peak:.6g} m (all snapshots scanned)")


if __name__ == "__main__":
    main()
