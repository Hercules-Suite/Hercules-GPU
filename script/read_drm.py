#!/usr/bin/env python3
"""Inspect a flat native-float64 station or DRM displacement file (Python standard library)."""

import argparse
from array import array
import math
from pathlib import Path


def box_layout(length_x, length_y, depth, spacing):
    """Return (name, rows, columns, point_offset) for all ten complete faces."""
    if not all(math.isfinite(v) for v in (length_x, length_y, depth, spacing)) or spacing <= 0:
        raise ValueError("box dimensions must be finite and spacing positive")
    counts = []
    for length in (length_x, length_y, depth):
        count = length / spacing
        if not math.isfinite(count):
            raise ValueError("box dimensions exceed count limits")
        rounded = round(count)
        if rounded < 1 or rounded > 2147483644 or abs(count-rounded) > 1e-8 + 16*math.ulp(1.0)*abs(count):
            raise ValueError("box dimensions must be positive integer multiples of spacing")
        counts.append(rounded)
    nx, ny, nz = counts
    if nx < 3 or ny < 3 or nz < 2:
        raise ValueError("exterior lengths must exceed 2*spacing and depth must exceed spacing")
    layout, offset = [], 0
    for side, label in enumerate(("exterior", "interior")):
        x, y, z = nx - 2*side, ny - 2*side, nz - side
        for face, rows, columns in (("front", y+1, z+1), ("back", y+1, z+1),
                                    ("left", x+1, z+1), ("right", x+1, z+1),
                                    ("bottom", x+1, y+1)):
            layout.append((f"{label}_{face}", rows, columns, offset))
            offset += rows * columns
    if offset > 2147483647 // 3:
        raise ValueError("box faces exceed MPI point-count limits")
    return layout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", type=Path)
    geometry = parser.add_mutually_exclusive_group(required=True)
    geometry.add_argument("--points", type=int, help="Number of output_stations rows per snapshot")
    geometry.add_argument("--box", type=float, nargs=4, metavar=("LX", "LY", "DEPTH", "SPACING"),
                          help="Exterior box dimensions and spacing in metres")
    parser.add_argument("--dt", type=float, required=True, help="Saved snapshot interval in seconds")
    parser.add_argument("--first-time", type=float, default=0.0)
    args = parser.parse_args()
    layout = []
    if args.box:
        try:
            layout = box_layout(*args.box)
        except ValueError as error:
            parser.error(str(error))
        args.points = sum(rows * columns for _, rows, columns, _ in layout)
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
    for name, rows, columns, offset in layout:
        print(f"{name}: shape=({frames}, {rows}, {columns}, 3), point_offset={offset}")
    print(f"time range={args.first_time:g}..{args.first_time + (frames - 1) * args.dt:g} s")
    print("peak absolute XYZ displacement (m):", *peak)


if __name__ == "__main__":
    main()
