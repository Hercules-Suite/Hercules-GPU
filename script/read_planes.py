#!/usr/bin/env python3
"""Inspect native-float64 plane displacements with manual grid/timing (NumPy)."""

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
