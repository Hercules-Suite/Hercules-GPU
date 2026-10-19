#!/usr/bin/env python3
"""Summarize station time and displacement; inspect its header for optional columns."""

# USAGE (Python 3 + NumPy; run this example from the repository root)
#   python3 script/read_stations.py outputfiles/stations/station.0
#   python3 script/read_stations.py --help
#
# Required: one completed ASCII station.N file. Relative paths use the working
# directory. The first four columns must be time [s], ux/uy/uz [m], with finite
# values and strictly increasing time. Leading # comment headers are echoed.
# Output: sample count, saved-time interval range, time span, column count,
# and peak displacement VECTOR magnitude sqrt(ux^2+uy^2+uz^2).
# Optional velocity/acceleration/pressure/stress columns are not checked or
# summarized; inspect the echoed header for their meaning. This is not a PGA
# calculation. The full table is loaded into memory; no files are written.
# For flat binary station_displacements output, use read_drm.py --points.
# See script/README.md for output-format details.

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
