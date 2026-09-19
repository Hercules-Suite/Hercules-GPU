#!/usr/bin/env python3
"""Read a rank's packed native-endian mesh coordinate/material files (NumPy)."""

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
