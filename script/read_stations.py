#!/usr/bin/env python3
"""Summarize station time and displacement; inspect its header for optional columns."""

import argparse
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", type=Path, help="ASCII station.N file; requires NumPy")
    args = parser.parse_args()
    # Echo the solver's header so optional velocity/pressure/stress columns are visible.
    with args.file.open() as source:
        for line in source:
            if line.lstrip().startswith("#"):
                print(line.rstrip())
            elif line.strip():
                break
    station = np.loadtxt(args.file, comments="#", ndmin=2)
    if station.shape[0] == 0 or station.shape[1] < 4:
        parser.error("need time and three displacement columns")
    time, displacement = station[:, 0], station[:, 1:4]
    if not np.isfinite(station[:, :4]).all() or np.any(np.diff(time) <= 0):
        parser.error("time/displacement must be finite and time strictly increasing")
    if time.size > 1:
        steps = np.diff(time)
        print(f"{time.size} samples, dt_out range={steps.min():.6g}..{steps.max():.6g} s")
    else:
        print("1 sample; output interval unavailable")
    print(f"time range={time[0]:.6g}..{time[-1]:.6g} s; columns={station.shape[1]}")
    print(f"peak |u| = {np.linalg.norm(displacement, axis=1).max():.6g} m")


if __name__ == "__main__":
    main()
