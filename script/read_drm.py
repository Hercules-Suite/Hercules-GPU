#!/usr/bin/env python3
"""Inspect a flat native-float64 DRM displacement file (Python standard library)."""

import argparse
from array import array
import math
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", type=Path)
    parser.add_argument("--points", type=int, required=True, help="Number of point rows per snapshot")
    parser.add_argument("--dt", type=float, required=True, help="Saved snapshot interval in seconds")
    parser.add_argument("--first-time", type=float, default=0.0)
    args = parser.parse_args()
    if args.points <= 0 or not math.isfinite(args.dt) or args.dt <= 0:
        parser.error("points and dt must be positive and finite")
    if not math.isfinite(args.first_time):
        parser.error("first-time must be finite")
    frame_bytes = args.points * 3 * 8
    frames, remainder = divmod(args.file.stat().st_size, frame_bytes)
    if remainder or not frames:
        parser.error(f"file must contain complete {frame_bytes}-byte snapshots")
    peak = [0.0, 0.0, 0.0]
    with args.file.open("rb") as source:
        for frame in range(frames):
            values = array("d")
            try:
                values.fromfile(source, args.points * 3)
            except EOFError:
                parser.error("file changed or was truncated while reading")
            for i, value in enumerate(values):
                if not math.isfinite(value):
                    parser.error(f"snapshot {frame}, point {i // 3} contains non-finite displacement")
                peak[i % 3] = max(peak[i % 3], abs(value))
    print(f"layout=({frames}, {args.points}, 3): time, point, XYZ")
    print(f"time range={args.first_time:g}..{args.first_time + (frames - 1) * args.dt:g} s")
    print("peak absolute XYZ displacement (m):", *peak)


if __name__ == "__main__":
    main()
