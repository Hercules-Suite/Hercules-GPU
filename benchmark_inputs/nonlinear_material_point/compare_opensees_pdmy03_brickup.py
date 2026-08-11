#!/usr/bin/env python3
"""Compare coupled OpenSees SSPbrickUP and Hercules PDMY03 BrickUP histories."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


SUPPORTED_TESTS = {
    "pdmy03_sspbrickup_undrained_triaxial": {
        "confinement_components": (0, 1),
        "confinement_label": "radial",
        "control_field": "axial_compression",
    },
    "pdmy03_sspbrickup_undrained_cyclic_shear": {
        "confinement_components": (2,),
        "confinement_label": "vertical",
        "control_field": "engineering_gamma_xz",
    },
}
DRIVER_VERSION = "2"
COMPONENTS = ("xx", "yy", "zz", "xy", "yz", "xz")
STRAINS = tuple(f"eps_{component}" for component in COMPONENTS)
PDMY03_DIAGNOSTIC_FIELDS = (
    "pdmy03_active_surface",
    "pdmy03_requested_substeps",
    "pdmy03_crossed_surfaces",
    "pdmy03_substeps_clamped",
    "pdmy03_pressure_d",
    "pdmy03_center_norm",
)


def read_hercules(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.reader(stream)
        try:
            header = next(reader)
        except StopIteration:
            raise SystemExit("Hercules PDMY03 BrickUP history is empty")
        raw_rows = list(reader)
    if not raw_rows:
        raise SystemExit("Hercules PDMY03 BrickUP history is empty")

    # PDMY03 initially appended six diagnostics to each material-point row
    # without appending their names to the shared CSV header. Repair those
    # already-produced histories in memory so users do not need to rebuild and
    # rerun a valid constitutive test merely to compare it. New Hercules output
    # contains these names directly and takes the normal path below.
    legacy_header = not any(field in header
                            for field in PDMY03_DIAGNOSTIC_FIELDS)
    if legacy_header and all(
            len(row) == len(header) + len(PDMY03_DIAGNOSTIC_FIELDS)
            for row in raw_rows):
        try:
            insertion = header.index("accepted_substeps") + 1
        except ValueError:
            raise SystemExit(
                "legacy Hercules CSV has no accepted_substeps column"
            )
        header[insertion:insertion] = PDMY03_DIAGNOSTIC_FIELDS

    for index, row in enumerate(raw_rows):
        if len(row) != len(header):
            raise SystemExit(
                f"Hercules CSV column mismatch at row {index}: "
                f"header={len(header)} data={len(row)}"
            )
    source = [dict(zip(header, row)) for row in raw_rows]
    rows: list[dict[str, float]] = []
    for index, row in enumerate(source):
        if int(float(row["status"])) != 0:
            raise SystemExit(
                f"Hercules constitutive update failed at replay step {index}"
            )
        rows.append({key: float(value) for key, value in row.items()
                     if key != "status"})
    required = ("pore_pressure_pa", *(f"total_sig_{c}" for c in COMPONENTS))
    missing = [field for field in required if field not in rows[0]]
    if missing:
        raise SystemExit("Hercules CSV is not a BrickUP history: " + ", ".join(missing))
    return rows


def read_opensees(path: Path) -> tuple[str, list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as stream:
        source = list(csv.DictReader(stream))
    if not source:
        raise SystemExit("OpenSees PDMY03 SSPbrickUP history is empty")
    test_id = source[0].get("test_id", "")
    if test_id not in SUPPORTED_TESTS:
        raise SystemExit(
            "unsupported OpenSees PDMY03 SSPbrickUP test_id: " + test_id
        )
    for index, row in enumerate(source):
        if row.get("test_id") != test_id or row.get("driver_version") != DRIVER_VERSION:
            raise SystemExit(
                f"OpenSees CSV is not {test_id} driver version {DRIVER_VERSION} "
                f"at row {index}"
            )
    if source[-1].get("phase") != "complete":
        raise SystemExit("OpenSees PDMY03 SSPbrickUP history is incomplete")
    base = max(index for index, row in enumerate(source)
               if row.get("phase") == "consolidation")
    return test_id, source[base:]


def normalized_rmse(errors: list[float], references: list[float]) -> float:
    scale = max((abs(value) for value in references), default=0.0)
    if scale == 0.0:
        return 0.0 if not any(errors) else math.inf
    return math.sqrt(sum(error*error for error in errors) / len(errors)) / scale


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hercules-csv", type=Path, required=True)
    parser.add_argument("--opensees-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path,
                        default=Path("opensees_pdmy03_sspbrickup_comparison"))
    parser.add_argument("--max-effective-stress-nrmse", type=float)
    parser.add_argument("--max-pore-pressure-difference-kpa", type=float)
    parser.add_argument("--max-total-confinement-error-kpa", type=float)
    args = parser.parse_args()

    hercules = read_hercules(args.hercules_csv)
    test_id, opensees = read_opensees(args.opensees_csv)
    test = SUPPORTED_TESTS[test_id]
    if len(hercules) != len(opensees):
        raise SystemExit(
            f"history length mismatch: Hercules={len(hercules)} "
            f"OpenSees={len(opensees)}"
        )

    effective_errors: list[float] = []
    effective_references: list[float] = []
    total_errors: list[float] = []
    total_references: list[float] = []
    maximum_strain_difference = 0.0
    maximum_stress_difference = 0.0
    maximum_pressure_difference = 0.0
    maximum_ru_difference = 0.0
    maximum_hercules_confinement_error = 0.0
    maximum_opensees_confinement_error = 0.0
    maximum_hercules_ru = 0.0
    maximum_opensees_ru = 0.0
    comparison: list[dict[str, float | int | str]] = []

    for index, (hrow, orow) in enumerate(zip(hercules, opensees)):
        for field in STRAINS:
            maximum_strain_difference = max(
                maximum_strain_difference,
                abs(hrow[field] - float(orow[field])),
            )
        h_eff = [hrow[f"sig_{component}"] for component in COMPONENTS]
        o_eff = [1000.0*float(orow[f"sigma_{component}_kpa"])
                 for component in COMPONENTS]
        h_total = [hrow[f"total_sig_{component}"] for component in COMPONENTS]
        o_total = [1000.0*float(orow[f"total_sigma_{component}_kpa"])
                   for component in COMPONENTS]
        for hvalue, ovalue in zip(h_eff, o_eff):
            error = hvalue - ovalue
            effective_errors.append(error)
            effective_references.append(ovalue)
            maximum_stress_difference = max(maximum_stress_difference, abs(error))
        for hvalue, ovalue in zip(h_total, o_total):
            total_errors.append(hvalue - ovalue)
            total_references.append(ovalue)

        hp = hrow["pore_pressure_pa"]
        op = 1000.0*float(orow["pore_pressure_kpa"])
        maximum_pressure_difference = max(maximum_pressure_difference,
                                          abs(hp-op))
        confinement = 1000.0*float(orow["confinement_kpa"])
        h_ru = hp/confinement
        o_ru = op/confinement
        maximum_hercules_ru = max(maximum_hercules_ru, h_ru)
        maximum_opensees_ru = max(maximum_opensees_ru, o_ru)
        maximum_ru_difference = max(maximum_ru_difference, abs(h_ru-o_ru))
        h_confinement = max(
            abs(-h_total[component]-confinement)
            for component in test["confinement_components"]
        )
        o_confinement = max(
            abs(-o_total[component]-confinement)
            for component in test["confinement_components"]
        )
        maximum_hercules_confinement_error = max(
            maximum_hercules_confinement_error, h_confinement
        )
        maximum_opensees_confinement_error = max(
            maximum_opensees_confinement_error, o_confinement
        )
        control_value = (
            -float(orow["eps_zz"])
            if test_id == "pdmy03_sspbrickup_undrained_triaxial"
            else 2.0*float(orow["eps_xz"])
        )
        comparison.append({
            "local_step": index,
            "opensees_step": int(orow["step"]),
            "phase": orow["phase"],
            test["control_field"]: control_value,
            "hercules_pore_pressure_kpa": hp/1000.0,
            "opensees_pore_pressure_kpa": op/1000.0,
            "hercules_ru": h_ru,
            "opensees_ru": o_ru,
            "hercules_total_confinement_error_kpa": h_confinement/1000.0,
            "opensees_total_confinement_error_kpa": o_confinement/1000.0,
        })

    effective_nrmse = normalized_rmse(effective_errors, effective_references)
    total_nrmse = normalized_rmse(total_errors, total_references)
    summary = {
        "test_id": test_id,
        "states": len(hercules),
        "maximum_strain_difference": maximum_strain_difference,
        "maximum_effective_stress_difference_pa": maximum_stress_difference,
        "effective_stress_nrmse": effective_nrmse,
        "total_stress_nrmse": total_nrmse,
        "maximum_pore_pressure_difference_kpa": maximum_pressure_difference/1000.0,
        "maximum_ru_hercules": maximum_hercules_ru,
        "maximum_ru_opensees": maximum_opensees_ru,
        "maximum_ru_difference": maximum_ru_difference,
        "confinement_label": test["confinement_label"],
        "maximum_hercules_total_confinement_error_kpa":
            maximum_hercules_confinement_error/1000.0,
        "maximum_opensees_total_confinement_error_kpa":
            maximum_opensees_confinement_error/1000.0,
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir/"comparison.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=comparison[0].keys())
        writer.writeheader()
        writer.writerows(comparison)
    (args.output_dir/"summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    print(f"Coupled canonical test: {test_id}")
    print(f"Hercules/OpenSees states: {len(hercules)}")
    print(f"Maximum strain difference: {maximum_strain_difference:.6e}")
    print(f"Maximum effective-stress difference: {maximum_stress_difference:.6e} Pa")
    print(f"Effective-stress NRMSE: {effective_nrmse:.6e}")
    print(f"Total-stress NRMSE: {total_nrmse:.6e}")
    print(f"Maximum pore-pressure difference: {maximum_pressure_difference/1000.0:.6e} kPa")
    print("Maximum ru (Hercules/OpenSees): "
          f"{maximum_hercules_ru:.6e} / {maximum_opensees_ru:.6e}")
    print(f"Maximum ru difference: {maximum_ru_difference:.6e}")
    print(f"Maximum total {test['confinement_label']} error "
          "(Hercules/OpenSees): "
          f"{maximum_hercules_confinement_error/1000.0:.6e} / "
          f"{maximum_opensees_confinement_error/1000.0:.6e} kPa")
    print(f"Detailed comparison: {(args.output_dir/'comparison.csv').resolve()}")
    print(f"Summary: {(args.output_dir/'summary.json').resolve()}")

    failures: list[str] = []
    if (args.max_effective_stress_nrmse is not None and
            effective_nrmse > args.max_effective_stress_nrmse):
        failures.append("effective-stress NRMSE")
    if (args.max_pore_pressure_difference_kpa is not None and
            maximum_pressure_difference/1000.0 >
            args.max_pore_pressure_difference_kpa):
        failures.append("pore-pressure difference")
    if (args.max_total_confinement_error_kpa is not None and
            maximum_hercules_confinement_error/1000.0 >
            args.max_total_confinement_error_kpa):
        failures.append("Hercules total confinement")
    if failures:
        raise SystemExit("FAIL: " + ", ".join(failures) + " exceeds threshold")
    if any(value is not None for value in (
            args.max_effective_stress_nrmse,
            args.max_pore_pressure_difference_kpa,
            args.max_total_confinement_error_kpa)):
        print("PASS: coupled PDMY03 BrickUP characterization agrees")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
