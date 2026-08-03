#!/usr/bin/env python3
"""Check CPU/GPU nonlinear station parity and the production speedup gate."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys


def numeric_rows(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if not fields or line.lstrip().startswith("#"):
            continue
        try:
            row = [float(field) for field in fields]
        except ValueError:
            continue
        rows.append(row)
    return rows


def compare_file(cpu_path: Path, gpu_path: Path, rtol: float,
                 atol: float) -> tuple[float, float, int]:
    cpu = numeric_rows(cpu_path)
    gpu = numeric_rows(gpu_path)
    if len(cpu) != len(gpu):
        raise ValueError(
            f"row count differs for {cpu_path.name}: CPU={len(cpu)} GPU={len(gpu)}"
        )

    max_abs = 0.0
    max_rel = 0.0
    compared = 0
    for row_index, (cpu_row, gpu_row) in enumerate(zip(cpu, gpu), start=1):
        if len(cpu_row) != len(gpu_row):
            raise ValueError(
                f"column count differs for {cpu_path.name} row {row_index}: "
                f"CPU={len(cpu_row)} GPU={len(gpu_row)}"
            )
        for column_index, (cpu_value, gpu_value) in enumerate(
            zip(cpu_row, gpu_row), start=1
        ):
            if not (math.isfinite(cpu_value) and math.isfinite(gpu_value)):
                if math.isnan(cpu_value) and math.isnan(gpu_value):
                    continue
                raise ValueError(
                    f"non-finite mismatch in {cpu_path.name} row {row_index} "
                    f"column {column_index}: CPU={cpu_value} GPU={gpu_value}"
                )
            absolute = abs(cpu_value - gpu_value)
            relative = absolute / max(abs(cpu_value), abs(gpu_value), 1.0e-30)
            max_abs = max(max_abs, absolute)
            max_rel = max(max_rel, relative)
            compared += 1
            if absolute > atol and relative > rtol:
                raise ValueError(
                    f"tolerance failure in {cpu_path.name} row {row_index} "
                    f"column {column_index}: CPU={cpu_value:.17e} "
                    f"GPU={gpu_value:.17e} abs={absolute:.3e} rel={relative:.3e}"
                )
    return max_abs, max_rel, compared


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--cpu-dir", required=True, type=Path)
    result.add_argument("--gpu-dir", required=True, type=Path)
    result.add_argument("--station-glob", default="station.*")
    result.add_argument("--cpu-seconds", required=True, type=float)
    result.add_argument("--gpu-seconds", required=True, type=float)
    result.add_argument("--nonlinear-elements", required=True, type=int)
    result.add_argument("--total-elements", required=True, type=int)
    result.add_argument("--rtol", default=1.0e-6, type=float)
    result.add_argument("--atol", default=1.0e-9, type=float)
    result.add_argument("--minimum-speedup", default=2.0, type=float)
    return result


def main() -> int:
    args = parser().parse_args()
    if args.total_elements <= 0 or args.nonlinear_elements < 0:
        raise ValueError("element counts must be nonnegative and total must be positive")
    nonlinear_fraction = args.nonlinear_elements / args.total_elements
    if nonlinear_fraction < 0.25:
        raise ValueError(
            f"production gate requires at least 25% nonlinear elements; "
            f"observed {nonlinear_fraction:.2%}"
        )
    if args.cpu_seconds <= 0.0 or args.gpu_seconds <= 0.0:
        raise ValueError("wall-clock measurements must be positive")

    cpu_files = sorted(args.cpu_dir.glob(args.station_glob))
    if not cpu_files:
        raise ValueError(
            f"no CPU station files match {args.station_glob!r} in {args.cpu_dir}"
        )

    maximum_absolute = 0.0
    maximum_relative = 0.0
    values_compared = 0
    for cpu_file in cpu_files:
        gpu_file = args.gpu_dir / cpu_file.name
        if not gpu_file.is_file():
            raise ValueError(f"missing GPU station file: {gpu_file}")
        absolute, relative, count = compare_file(
            cpu_file, gpu_file, args.rtol, args.atol
        )
        maximum_absolute = max(maximum_absolute, absolute)
        maximum_relative = max(maximum_relative, relative)
        values_compared += count

    speedup = args.cpu_seconds / args.gpu_seconds
    print(f"nonlinear fraction: {nonlinear_fraction:.2%}")
    print(f"station files: {len(cpu_files)}; values compared: {values_compared}")
    print(f"maximum absolute difference: {maximum_absolute:.6e}")
    print(f"maximum relative difference: {maximum_relative:.6e}")
    print(f"end-to-end speedup: {speedup:.3f}x")
    if speedup < args.minimum_speedup:
        print(
            f"FAIL: speedup is below {args.minimum_speedup:.3f}x",
            file=sys.stderr,
        )
        return 1
    print("PASS: parity tolerances and performance gate satisfied")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValueError as error:
        print(f"FAIL: {error}", file=sys.stderr)
        raise SystemExit(2)
