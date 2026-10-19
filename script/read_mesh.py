#!/usr/bin/env python3
"""Read a rank's packed native-endian mesh coordinate/material files (NumPy)."""

# USAGE (Python 3 + NumPy; run this example from the repository root)
#   python3 script/read_mesh.py \
#       outputfiles/mesh_coordinates.0 outputfiles/mesh_data.0
#   python3 script/read_mesh.py --help
#
# Required: coordinate file then material file from the SAME run and MPI rank;
# use their actual output directory. Relative paths use the working directory.
# Packed native-endian formats: coordinates = int64 ID + float64 X/Y/Z (32
# bytes per corner); materials = int64 ID + float32 Vs/Vp/rho (20 bytes/element).
# There must be eight consecutive coordinate records per material record with
# matching IDs. Coordinates are metres, X north/Y east/Z down; speeds m/s;
# density kg/m^3 (unlike tonne/m^3 in material-input rho.in).
# Output: element/corner counts and coordinate ranges. Both files are loaded
# into memory; incomplete records, ID mismatches and nonfinite values fail.
# No files are written. Finite sentinel material values are not rejected just
# because they are negative. This is NOT a mesh.e/archive reader or plotter.
# See script/README.md for the binary layouts.

import argparse
from pathlib import Path

import numpy as np

COORD_DTYPE = np.dtype([("geid", np.int64), ("x", np.float64),
                        ("y", np.float64), ("z", np.float64)])
MATERIAL_DTYPE = np.dtype([("geid", np.int64), ("Vs", np.float32),
                           ("Vp", np.float32), ("rho", np.float32)])


def read_records(path, dtype):
    if path.stat().st_size % dtype.itemsize:
        raise ValueError(f"{path}: incomplete {dtype.itemsize}-byte record")
    return np.fromfile(path, dtype=dtype)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("coordinates", type=Path, help="mesh_coordinates.RANK")
    parser.add_argument("materials", type=Path, help="mesh_data.RANK for the same rank")
    args = parser.parse_args()
    coords = read_records(args.coordinates, COORD_DTYPE)
    materials = read_records(args.materials, MATERIAL_DTYPE)
    if coords.size != 8 * materials.size:
        parser.error("need eight coordinate records for each material record")
    if not np.array_equal(coords["geid"], np.repeat(materials["geid"], 8)):
        parser.error("coordinate and material element IDs do not match")
    for records in (coords, materials):
        for name in records.dtype.names[1:]:
            if not np.isfinite(records[name]).all():
                parser.error(f"non-finite {name} values")
    print(f"{materials.size} elements, {coords.size} corner records")
    if coords.size:
        for name in ("x", "y", "z"):
            print(f"{name} range={coords[name].min():g}..{coords[name].max():g} m")


if __name__ == "__main__":
    main()
