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
    "md_drained_triaxial": ("9", "complete"),
    "md_undrained_cyclic_shear": ("1", "cyclic"),
    "md_sspbrickup_undrained_triaxial": ("3", "complete"),
}


def normalize_isotropic_seed(
    initial: dict[str, str],
) -> tuple[tuple[float, ...], dict[str, tuple[float, ...]]]:
    """Remove roundoff-only directional state after isotropic consolidation."""
    stress_kpa = tuple(float(initial[field]) for field in STRESS_FIELDS)
    mean_normal_kpa = sum(stress_kpa[:3]) / 3.0
    stress_scale_kpa = max(abs(mean_normal_kpa), 1.0)
    maximum_stress_departure = max(
        *(abs(value - mean_normal_kpa) for value in stress_kpa[:3]),
        *(abs(value) for value in stress_kpa[3:]),
    )
    if maximum_stress_departure > 1.0e-10 * stress_scale_kpa:
        raise SystemExit(
            "OpenSees stage-boundary stress is not isotropic within the "
            f"roundoff-cleanup tolerance: {maximum_stress_departure:.3e} kPa"
        )

    for state_name, fields in STATE_TENSOR_FIELDS.items():
        maximum_state = max(abs(float(initial[field])) for field in fields)
        if maximum_state > 1.0e-10:
            raise SystemExit(
                f"OpenSees stage-boundary {state_name} is not zero within "
                f"the roundoff-cleanup tolerance: {maximum_state:.3e}"
            )

    stress_pa = (
        1000.0 * mean_normal_kpa,
        1000.0 * mean_normal_kpa,
        1000.0 * mean_normal_kpa,
        0.0,
        0.0,
        0.0,
    )
    zero_state = (0.0,) * 6
    return stress_pa, {
        state_name: zero_state for state_name in STATE_TENSOR_FIELDS
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
        if "confinement_kpa" not in rows[0]:
            raise SystemExit("drained OpenSees CSV is missing confinement_kpa")
        consolidated = [row for row in rows if row.get("phase") == "consolidation"]
        monotonic = [row for row in rows if row.get("phase") == "monotonic"]
        completed = [row for row in rows if row.get("phase") == "complete"]
        if (not consolidated or not monotonic or len(completed) != 1
                or rows[-1] is not completed[0]):
            raise SystemExit("incomplete drained path: missing final completion state")
        confinement_kpa = float(monotonic[0]["confinement_kpa"])
        maximum_confinement_error_kpa = max(
            abs(-float(row[field]) - confinement_kpa)
            for row in monotonic
            for field in ("sigma_xx_kpa", "sigma_yy_kpa")
        )
        if maximum_confinement_error_kpa > 1.0e-3:
            raise SystemExit(
                "invalid drained path: OpenSees radial material stress misses "
                "the confinement target by "
                f"{maximum_confinement_error_kpa:.6e} kPa"
            )
        axial_increment = (
            float(completed[0]["eps_zz"]) - float(consolidated[-1]["eps_zz"])
        )
        if not math.isclose(axial_increment, -0.15, rel_tol=0.0, abs_tol=1.0e-8):
            raise SystemExit(
                "incomplete drained path: final axial increment is "
                f"{axial_increment:.17e}, expected -1.5e-1"
            )
    elif test_id == "md_sspbrickup_undrained_triaxial":
        coupled_fields = (
            "pore_pressure_kpa", "total_sigma_xx_kpa",
            "total_sigma_yy_kpa", "total_sigma_zz_kpa",
            "fluid_bulk_modulus_kpa", "porosity", "biot_alpha",
            "confinement_kpa",
        )
        missing = [field for field in coupled_fields if field not in rows[0]]
        if missing:
            raise SystemExit(
                "OpenSees SSPbrickUP CSV is missing fields: "
                + ", ".join(missing)
            )
        consolidated = [row for row in rows
                        if row.get("phase") == "consolidation"]
        coupled = [row for row in rows if row.get("phase") == "coupled"]
        completed = [row for row in rows if row.get("phase") == "complete"]
        if (not consolidated or not coupled or len(completed) != 1
                or rows[-1] is not completed[0]):
            raise SystemExit(
                "incomplete SSPbrickUP path: missing coupled or completion state"
            )
        initial = consolidated[-1]
        initial_pressure = float(initial["pore_pressure_kpa"])
        if abs(initial_pressure) > 1.0e-8:
            raise SystemExit(
                "SSPbrickUP stage-boundary pore pressure is not zero: "
                f"{initial_pressure:.6e} kPa"
            )
        fluid_bulk = float(initial["fluid_bulk_modulus_kpa"])
        porosity = float(initial["porosity"])
        biot_alpha = float(initial["biot_alpha"])
        if fluid_bulk <= 0.0 or not 0.0 < porosity <= 1.0:
            raise SystemExit("invalid SSPbrickUP fluid bulk modulus or porosity")
        if not math.isclose(biot_alpha, 1.0, rel_tol=0.0, abs_tol=1.0e-12):
            raise SystemExit(
                "the SSPbrickUP/Hercules benchmark currently requires alpha=1"
            )
        storage_modulus = fluid_bulk / porosity
        initial_ev = sum(float(initial[field]) for field in STRAIN_FIELDS[:3])
        maximum_storage_error = max(
            abs(
                float(row["pore_pressure_kpa"])
                - initial_pressure
                + storage_modulus * (
                    sum(float(row[field]) for field in STRAIN_FIELDS[:3])
                    - initial_ev
                )
            )
            for row in coupled
        )
        if maximum_storage_error > 1.0e-3:
            raise SystemExit(
                "invalid SSPbrickUP sealed-storage path: maximum residual is "
                f"{maximum_storage_error:.6e} kPa"
            )
        confinement = float(initial["confinement_kpa"])
        maximum_confinement_error = max(
            abs(-float(row[field]) - confinement)
            for row in coupled
            for field in ("total_sigma_xx_kpa", "total_sigma_yy_kpa")
        )
        if maximum_confinement_error > 1.0e-2:
            raise SystemExit(
                "invalid SSPbrickUP total radial confinement path: maximum "
                f"error is {maximum_confinement_error:.6e} kPa"
            )
        axial_increment = (
            float(completed[0]["eps_zz"]) - float(initial["eps_zz"])
        )
        if not math.isclose(axial_increment, -0.05,
                            rel_tol=0.0, abs_tol=1.0e-8):
            raise SystemExit(
                "incomplete SSPbrickUP path: final axial increment is "
                f"{axial_increment:.17e}, expected -5.0e-2"
            )
    return test_id, rows


def replay_rows(test_id: str, rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Start at the committed state where OpenSees enables plasticity."""
    if test_id in ("md_drained_triaxial",
                    "md_sspbrickup_undrained_triaxial"):
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
    parser.add_argument("--backbone-errtol", type=float, default=1.0e-5)
    parser.add_argument("--yield-tolerance-pa", type=float, default=1.0e-1)
    parser.add_argument("--maximum-substeps", type=int, default=400)
    args = parser.parse_args()

    test_id, source_rows = read_path(args.opensees_csv)
    rows = replay_rows(test_id, source_rows)
    initial = rows[0]
    initial_stress_pa, initial_state_tensors = normalize_isotropic_seed(initial)
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
    ]
    if test_id == "md_sspbrickup_undrained_triaxial":
        fluid_bulk_pa = 1000.0 * float(initial["fluid_bulk_modulus_kpa"])
        lines.extend([
            "",
            "# Sealed uniform BrickUP coupling; alpha=1 matches SSPbrickUP.",
            "material_point_brickup = yes",
            ("material_point_brickup_porosity = "
             f"{float(initial['porosity']):.17e}"),
            ("material_point_brickup_fluid_bulk_modulus_pa = "
             f"{fluid_bulk_pa:.17e}"),
            ("material_point_brickup_biot_alpha = "
             f"{float(initial['biot_alpha']):.17e}"),
            ("material_point_initial_pore_pressure_pa = "
             f"{1000.0*float(initial['pore_pressure_kpa']):.17e}"),
        ])
    lines.extend([
        "",
        "# G0 nu e_init Mc c lambda_c e0 ksi P_atm m h0 ch nb A0 nd z_max cz",
        "material_point_constants =",
        "125.0 0.05 0.80 1.25 0.712 0.019 0.934 0.70 101300.0 0.01 7.05 0.968 1.10 0.704 3.50 4.0 600.0",
        "",
        "# Committed OpenSees state at the elastic-to-plastic stage boundary.",
        "# The verified isotropic seed is normalized to remove roundoff-only",
        "# directional state. Stress is tension-positive and expressed in Pa.",
        "material_point_initial_strain =",
        " ".join(f"{float(initial[field]):.17e}" for field in STRAIN_FIELDS),
        "material_point_initial_stress =",
        " ".join(f"{value:.17e}" for value in initial_stress_pa),
    ])
    for state_name in STATE_TENSOR_FIELDS:
        lines.extend([
            f"material_point_initial_{state_name} =",
            " ".join(f"{value:.17e}"
                     for value in initial_state_tensors[state_name]),
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
    print("Seed normalization: verified isotropic stress; zeroed "
          "roundoff-only directional state")
    if test_id == "md_sspbrickup_undrained_triaxial":
        print("BrickUP replay: sealed zero-permeability pressure update; "
              "production CUDA pressure kernels enabled")
    print(f"Hercules input: {args.output.resolve()}")
    print(f"Material-point output prefix: {prefix}")
    print(f"Expected output names: {prefix}.cpu.csv and {prefix}.gpu.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
