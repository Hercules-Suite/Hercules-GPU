#!/usr/bin/env python3
"""Convert a row-wise CSV model to Hercules count-prefixed 3-D model files.

Required CSV columns and units:
    x [m north], y [m east], z [m positive downward],
    vs [m/s], vp [m/s], rho [tonne/m^3]

Optional columns qs and qp must either both be present or both be absent.
"""

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


REQUIRED_COLUMNS = ("x", "y", "z", "vs", "vp", "rho")
QUALITY_COLUMNS = ("qs", "qp")


def require_uniform_axis(values, name):
    if len(values) < 2:
        raise ValueError(f"{name}: need at least two unique coordinates")
    spacing = values[1] - values[0]
    if spacing <= 0.0:
        raise ValueError(f"{name}: coordinates must increase")
    for left, right in zip(values, values[1:]):
        if not math.isclose(right - left, spacing, rel_tol=1.0e-9,
                            abs_tol=1.0e-8):
            raise ValueError(f"{name}: coordinates are not uniformly spaced")


def write_counted(path, values):
    with path.open("w", encoding="utf-8") as output:
        output.write(f"{len(values)}\n")
        for value in values:
            if isinstance(value, int):
                output.write(f"{value}\n")
            else:
                output.write(f"{value:.17g}\n")


def read_csv_rows(path):
    with path.open(newline="", encoding="utf-8-sig") as source:
        reader = csv.DictReader(source)
        if reader.fieldnames is None:
            raise ValueError("CSV is missing a header row")
        fieldnames = [name.strip().lower() for name in reader.fieldnames]
        if len(fieldnames) != len(set(fieldnames)):
            raise ValueError("CSV header contains duplicate column names")
        reader.fieldnames = fieldnames

        missing = [name for name in REQUIRED_COLUMNS if name not in fieldnames]
        if missing:
            raise ValueError(f"CSV is missing columns: {', '.join(missing)}")
        has_qs = "qs" in fieldnames
        has_qp = "qp" in fieldnames
        if has_qs != has_qp:
            raise ValueError("CSV must contain both qs and qp, or neither")
        property_columns = REQUIRED_COLUMNS + (QUALITY_COLUMNS if has_qs else ())

        rows = []
        for line_number, raw in enumerate(reader, start=2):
            try:
                row = {name: float(raw[name]) for name in property_columns}
            except (KeyError, TypeError, ValueError) as error:
                raise ValueError(
                    f"line {line_number}: every model field must be numeric"
                ) from error
            if not all(math.isfinite(value) for value in row.values()):
                raise ValueError(f"line {line_number}: values must be finite")
            if (row["z"] < 0.0 or row["vs"] < 0.0 or
                    row["vp"] <= 0.0 or row["rho"] <= 0.0):
                raise ValueError(
                    f"line {line_number}: require z >= 0, Vs >= 0, "
                    "Vp > 0, rho > 0"
                )
            if has_qs and (row["qs"] <= 0.0 or row["qp"] <= 0.0):
                raise ValueError(f"line {line_number}: Qs and Qp must be > 0")
            rows.append(row)

    if not rows:
        raise ValueError("CSV contains no model rows")
    return rows, property_columns


def convert(csv_path, output_dir):
    rows, property_columns = read_csv_rows(csv_path)
    x_values = sorted({row["x"] for row in rows})
    y_values = sorted({row["y"] for row in rows})
    require_uniform_axis(x_values, "x")
    require_uniform_axis(y_values, "y")

    profiles = defaultdict(list)
    for row in rows:
        profiles[(row["y"], row["x"])].append(row)

    expected_profiles = {
        (y_value, x_value)
        for y_value in y_values
        for x_value in x_values
    }
    missing_profiles = expected_profiles.difference(profiles)
    if missing_profiles:
        missing_text = ", ".join(
            f"(x={x:g}, y={y:g})" for y, x in sorted(missing_profiles)
        )
        raise ValueError(f"rectangular grid is missing profiles: {missing_text}")

    starts = []
    ordered_rows = []
    for y_value in y_values:          # Y is the outer loop.
        for x_value in x_values:      # X is the inner loop.
            profile = sorted(
                profiles[(y_value, x_value)], key=lambda row: row["z"]
            )
            depths = [row["z"] for row in profile]
            if any(right <= left for left, right in zip(depths, depths[1:])):
                raise ValueError(
                    f"profile (x={x_value:g}, y={y_value:g}) has duplicate depths"
                )
            starts.append(len(ordered_rows))
            ordered_rows.extend(profile)

    output_dir.mkdir(parents=True, exist_ok=True)
    write_counted(output_dir / "x.in", x_values)
    write_counted(output_dir / "y.in", y_values)
    write_counted(output_dir / "index.in", starts)
    for name in property_columns:
        if name not in ("x", "y"):
            write_counted(output_dir / f"{name}.in",
                          [row[name] for row in ordered_rows])

    print(
        f"wrote {len(x_values)} x {len(y_values)} profiles and "
        f"{len(ordered_rows)} samples to {output_dir}"
    )


def main():
    parser = argparse.ArgumentParser(
        description="Convert CSV to Hercules 3-D velocity-model .in files"
    )
    parser.add_argument("csv_file", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    convert(args.csv_file, args.output_dir)


if __name__ == "__main__":
    main()
