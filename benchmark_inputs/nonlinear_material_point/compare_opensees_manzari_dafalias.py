#!/usr/bin/env python3
"""Characterize Hercules Manzari--Dafalias against the OpenSees fixture."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


STRESS_PAIRS = (
    ("sig_xx", "sigma_xx_kpa"),
    ("sig_yy", "sigma_yy_kpa"),
    ("sig_zz", "sigma_zz_kpa"),
    ("sig_xy", "sigma_xy_kpa"),
    ("sig_yz", "sigma_yz_kpa"),
    ("sig_xz", "sigma_xz_kpa"),
)
STATE_PAIRS = tuple(
    [(f"eta_{component}", f"alpha_{component}")
     for component in ("xx", "yy", "zz", "xy", "yz", "xz")]
    + [(f"fabric_{component}", f"fabric_{component}")
       for component in ("xx", "yy", "zz", "xy", "yz", "xz")]
    + [(f"alpha_in_{component}", f"alpha_in_{component}")
       for component in ("xx", "yy", "zz", "xy", "yz", "xz")]
    + [("void_ratio", "void_ratio")]
)


def read_rows(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return [
            {key: float(value) for key, value in row.items() if key != "status"}
            for row in csv.DictReader(stream)
        ]


def normalized_rmse(errors: list[float], reference: list[float]) -> float:
    scale = max((abs(value) for value in reference), default=0.0)
    if scale == 0.0:
        return 0.0 if not any(errors) else math.inf
    return math.sqrt(sum(value * value for value in errors) / len(errors)) / scale


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hercules-csv", type=Path, required=True)
    parser.add_argument("--opensees-csv", type=Path, required=True)
    parser.add_argument("--integration-scheme",
                        choices=("modified_euler", "forward_euler"),
                        default="modified_euler")
    parser.add_argument("--output-dir", type=Path,
                        default=Path("opensees_manzari_dafalias_comparison"))
    args = parser.parse_args()

    hercules = read_rows(args.hercules_csv)
    opensees = read_rows(args.opensees_csv)
    if len(hercules) != len(opensees):
        raise SystemExit(
            f"history length mismatch: Hercules={len(hercules)} "
            f"OpenSees={len(opensees)}"
        )
    if not hercules:
        raise SystemExit("comparison histories are empty")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    detailed_path = args.output_dir / "comparison.csv"
    stress_errors: list[float] = []
    stress_references: list[float] = []
    state_errors: list[float] = []
    max_strain_error = 0.0
    detailed: list[dict[str, float]] = []
    expected_scheme_id = 1.0 if args.integration_scheme == "modified_euler" else 5.0

    for step, (hrow, orow) in enumerate(zip(hercules, opensees)):
        driver_version = orow.get("driver_version")
        if driver_version not in (1.0, 2.0):
            raise SystemExit("OpenSees CSV has the wrong driver version")
        actual_scheme_id = (1.0 if driver_version == 1.0 else
                            orow.get("integration_scheme_id"))
        if actual_scheme_id != expected_scheme_id:
            raise SystemExit(
                "OpenSees CSV integration scheme does not match "
                f"--integration-scheme={args.integration_scheme}"
            )
        expected_time = float(step + 1)
        if not math.isclose(orow.get("domain_time", math.nan), expected_time,
                            rel_tol=0.0, abs_tol=1.0e-12):
            raise SystemExit(f"OpenSees domain time mismatch at step {step}")
        for hfield, ofield in (("eps_xx", "eps_xx"),
                               ("eps_yy", "eps_yy"),
                               ("eps_zz", "eps_zz"),
                               ("eps_xy", "eps_xy")):
            max_strain_error = max(
                max_strain_error, abs(hrow[hfield] - orow[ofield]))

        output = {"step": float(step), "eps_xy": hrow["eps_xy"]}
        for hfield, ofield in STRESS_PAIRS:
            reference_pa = 1000.0 * orow[ofield]
            error = hrow[hfield] - reference_pa
            stress_errors.append(error)
            stress_references.append(reference_pa)
            output[f"hercules_{hfield}"] = hrow[hfield]
            output[f"opensees_{hfield}_pa"] = reference_pa
            output[f"delta_{hfield}_pa"] = error
        for hfield, ofield in STATE_PAIRS:
            error = hrow[hfield] - orow[ofield]
            state_errors.append(error)
            output[f"delta_{hfield}"] = error
        detailed.append(output)

    with detailed_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(detailed[0]))
        writer.writeheader()
        writer.writerows(detailed)

    summary = {
        "steps": len(hercules),
        "maximum_strain_difference": max_strain_error,
        "maximum_stress_difference_pa": max(abs(value) for value in stress_errors),
        "normalized_stress_rmse": normalized_rmse(stress_errors, stress_references),
        "maximum_dimensionless_state_difference": max(
            abs(value) for value in state_errors
        ),
        "integration_scheme": args.integration_scheme,
        "note": (
            "Characterization only: the Hercules and OpenSees implementations "
            f"use related {args.integration_scheme.replace('_', ' ').title()} "
            "algorithms but are not expected to be bitwise identical."
        ),
    }
    summary_path = args.output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(f"Hercules/OpenSees steps: {summary['steps']}")
    print(f"Maximum strain difference: {max_strain_error:.6e}")
    print("Maximum stress difference: "
          f"{summary['maximum_stress_difference_pa']:.6e} Pa")
    print(f"Normalized stress RMSE: {summary['normalized_stress_rmse']:.6e}")
    print("Maximum alpha/fabric/void-ratio difference: "
          f"{summary['maximum_dimensionless_state_difference']:.6e}")
    print(f"Detailed comparison: {detailed_path.resolve()}")
    print(f"Summary: {summary_path.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
