#!/usr/bin/env python3
"""Convert a canonical OpenSees MD strain path to a Hercules point input."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


STRAIN_FIELDS = ("eps_xx", "eps_yy", "eps_zz", "eps_xy", "eps_yz", "eps_xz")
STRESS_FIELDS = tuple(f"sigma_{component}_kpa"
                      for component in ("xx", "yy", "zz", "xy", "yz", "xz"))
STATE_TENSOR_FIELDS = {
    "alpha": tuple(f"alpha_{component}"
                   for component in ("xx", "yy", "zz", "xy", "yz", "xz")),
    "alpha_in": tuple(f"alpha_in_{component}"
                      for component in ("xx", "yy", "zz", "xy", "yz", "xz")),
    "fabric": tuple(f"fabric_{component}"
                    for component in ("xx", "yy", "zz", "xy", "yz", "xz")),
}
KNOWN_TESTS = {
    "md_drained_triaxial": ("5", "complete"),
    "md_undrained_cyclic_shear": ("1", "cyclic"),
}


def read_path(path: Path) -> tuple[str, list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise SystemExit(f"OpenSees path is empty: {path}")
    required = (*STRAIN_FIELDS, *STRESS_FIELDS,
                *(field for fields in STATE_TENSOR_FIELDS.values()
                  for field in fields),
                "test_id", "driver_version", "phase")
    missing = [name for name in required
               if name not in rows[0]]
    if missing:
        raise SystemExit("OpenSees CSV is missing fields: " + ", ".join(missing))
    test_ids = {row["test_id"] for row in rows}
    if len(test_ids) != 1 or next(iter(test_ids)) not in KNOWN_TESTS:
        raise SystemExit(f"unsupported or inconsistent test_id values: {test_ids}")
    test_id = next(iter(test_ids))
    expected_version, required_phase = KNOWN_TESTS[test_id]
    if any(row["driver_version"] != expected_version for row in rows):
        raise SystemExit(
            f"OpenSees CSV was not produced by {test_id} driver version "
            f"{expected_version}"
        )
    if not any(row.get("phase") == required_phase for row in rows):
        raise SystemExit(
            f"incomplete OpenSees CSV: no {required_phase!r} states were recorded"
        )
    for index, row in enumerate(rows):
        if int(row["step"]) != index:
            raise SystemExit(f"nonsequential OpenSees step at row {index}")
        values = [float(row[field]) for field in STRAIN_FIELDS]
        if not all(math.isfinite(value) for value in values):
            raise SystemExit(f"nonfinite strain at OpenSees step {index}")
    if test_id == "md_drained_triaxial":
        consolidated = [row for row in rows if row.get("phase") == "consolidation"]
        completed = [row for row in rows if row.get("phase") == "complete"]
        if not consolidated or len(completed) != 1 or rows[-1] is not completed[0]:
            raise SystemExit("incomplete drained path: missing final completion state")
        axial_increment = (
            float(completed[0]["eps_zz"]) - float(consolidated[-1]["eps_zz"])
        )
        if not math.isclose(axial_increment, -0.15, rel_tol=0.0, abs_tol=1.0e-8):
            raise SystemExit(
                "incomplete drained path: final axial increment is "
                f"{axial_increment:.17e}, expected -1.5e-1"
            )
    return test_id, rows


def replay_rows(test_id: str, rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Start at the committed state where OpenSees enables plasticity."""
    if test_id == "md_drained_triaxial":
        base = max(index for index, row in enumerate(rows)
                   if row["phase"] == "consolidation")
    else:
        base = next(index for index, row in enumerate(rows)
                    if row["phase"] == "cyclic")
    return rows[base:]


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a Hercules --nonlinear-material-point input that replays "
            "the exact strain path solved by a canonical OpenSees driver."
        )
    )
    parser.add_argument("--opensees-csv", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--output-prefix")
    parser.add_argument("--backbone-errtol", type=float, default=1.0e-4)
    parser.add_argument("--yield-tolerance-pa", type=float, default=1.0e-1)
    parser.add_argument("--maximum-substeps", type=int, default=100)
    args = parser.parse_args()

    test_id, source_rows = read_path(args.opensees_csv)
    rows = replay_rows(test_id, source_rows)
    initial = rows[0]
    prefix = args.output_prefix or f"hercules_{test_id}"
    if args.backbone_errtol <= 0.0:
        raise SystemExit("--backbone-errtol must be positive")
    if args.yield_tolerance_pa <= 0.0:
        raise SystemExit("--yield-tolerance-pa must be positive")
    if args.maximum_substeps <= 0:
        raise SystemExit("--maximum-substeps must be positive")

    lines = [
        "# Generated from a canonical OpenSees Manzari--Dafalias strain path.",
        f"# Source: {args.opensees_csv}",
        f"# Test: {test_id}",
        "# Constants are numerical characterization inputs, not a soil calibration.",
        "material_model = ManzariDafalias",
        "material_point_delta_time_sec = 1.0",
        f"material_point_steps = {len(rows)}",
        f"no_substeps = {args.maximum_substeps}",
        f"error_tolerance = {args.yield_tolerance_pa:.17e}",
        f"backbone_errTol = {args.backbone_errtol:.17e}",
        "tension_cutoff = no",
        # Long adaptive histories accumulate small host/device libm differences.
        # Tensor components use their tensor magnitude as the relative scale;
        # constitutive convergence is checked separately against error_tolerance.
        "material_point_compare_rtol = 1.0e-6",
        "material_point_compare_atol = 1.0e-8",
        f"material_point_output_prefix = {prefix}",
        "",
        "# G0 nu e_init Mc c lambda_c e0 ksi P_atm m h0 ch nb A0 nd z_max cz",
        "material_point_constants =",
        "125.0 0.05 0.80 1.25 0.712 0.019 0.934 0.70 101300.0 0.01 7.05 0.968 1.10 0.704 3.50 4.0 600.0",
        "",
        "# Committed OpenSees state at the elastic-to-plastic stage boundary.",
        "# Stress is tension-positive and expressed in Pa for Hercules.",
        "material_point_initial_strain =",
        " ".join(f"{float(initial[field]):.17e}" for field in STRAIN_FIELDS),
        "material_point_initial_stress =",
        " ".join(f"{1000.0*float(initial[field]):.17e}"
                 for field in STRESS_FIELDS),
    ]
    for state_name, fields in STATE_TENSOR_FIELDS.items():
        lines.extend([
            f"material_point_initial_{state_name} =",
            " ".join(f"{float(initial[field]):.17e}" for field in fields),
        ])
    lines.extend([
        "",
        "# Hercules tensor strain: exx eyy ezz exy eyz exz",
        "material_point_strain_history =",
    ])
    for row in rows:
        lines.append(" ".join(f"{float(row[field]):.17e}" for field in STRAIN_FIELDS))
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"Canonical test: {test_id}")
    print(f"OpenSees source states: {len(source_rows)}")
    print(f"Stage-1 replay states: {len(rows)}")
    print(f"Initial OpenSees step: {initial['step']}")
    print(f"Hercules input: {args.output.resolve()}")
    print(f"Material-point output prefix: {prefix}")
    print(f"Expected output names: {prefix}.cpu.csv and {prefix}.gpu.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
