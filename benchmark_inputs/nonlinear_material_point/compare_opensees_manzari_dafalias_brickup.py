#!/usr/bin/env python3
"""Compare the coupled SSPbrickUP/Hercules BrickUP one-element history."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


STRAIN_FIELDS = ("eps_xx", "eps_yy", "eps_zz", "eps_xy", "eps_yz", "eps_xz")
COMPONENTS = ("xx", "yy", "zz", "xy", "yz", "xz")
STATE_PAIRS = tuple(
    [(f"eta_{component}", f"alpha_{component}") for component in COMPONENTS]
    + [(f"fabric_{component}", f"fabric_{component}")
       for component in COMPONENTS]
    + [(f"alpha_in_{component}", f"alpha_in_{component}")
       for component in COMPONENTS]
    + [("void_ratio", "void_ratio")]
)
TEST_ID = "md_sspbrickup_undrained_triaxial"
DRIVER_VERSION = "2"


def read_hercules(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows: list[dict[str, float]] = []
        for index, row in enumerate(csv.DictReader(stream)):
            if int(row["status"]) != 0:
                raise SystemExit(
                    "Hercules coupled history is invalid: constitutive update "
                    f"failed at local replay step {index}"
                )
            rows.append({key: float(value) for key, value in row.items()
                         if key != "status"})
    if not rows:
        raise SystemExit("Hercules coupled history is empty")
    required = ("pore_pressure_pa", *(f"total_sig_{c}" for c in COMPONENTS))
    missing = [field for field in required if field not in rows[0]]
    if missing:
        raise SystemExit(
            "Hercules CSV was not generated with material_point_brickup=yes: "
            + ", ".join(missing)
        )
    return rows


def read_opensees(path: Path) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as stream:
        source = list(csv.DictReader(stream))
    if not source:
        raise SystemExit("OpenSees SSPbrickUP history is empty")
    for index, row in enumerate(source):
        if (row.get("test_id") != TEST_ID
                or row.get("driver_version") != DRIVER_VERSION):
            raise SystemExit(
                f"OpenSees CSV is not {TEST_ID} driver version "
                f"{DRIVER_VERSION} at row {index}"
            )
        if int(row["step"]) != index:
            raise SystemExit(f"nonsequential OpenSees step at row {index}")
    if source[-1].get("phase") != "complete":
        raise SystemExit("OpenSees SSPbrickUP history is incomplete")
    base = max(index for index, row in enumerate(source)
               if row.get("phase") == "consolidation")
    replay = source[base:]
    if not any(row.get("phase") == "coupled" for row in replay):
        raise SystemExit("OpenSees SSPbrickUP history has no coupled states")
    return source, replay


def normalized_rmse(errors: list[float], reference: list[float]) -> float:
    scale = max((abs(value) for value in reference), default=0.0)
    if scale == 0.0:
        return 0.0 if not any(errors) else math.inf
    return math.sqrt(sum(value * value for value in errors) / len(errors)) / scale


def invariants(stress: tuple[float, ...]) -> tuple[float, float]:
    sxx, syy, szz, sxy, syz, sxz = stress
    pressure = -(sxx + syy + szz) / 3.0
    dxx = sxx + pressure
    dyy = syy + pressure
    dzz = szz + pressure
    q = math.sqrt(1.5 * (
        dxx*dxx + dyy*dyy + dzz*dzz
        + 2.0*(sxy*sxy + syz*syz + sxz*sxz)
    ))
    return pressure, q


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hercules-csv", type=Path, required=True)
    parser.add_argument("--opensees-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path,
                        default=Path("opensees_md_sspbrickup_comparison"))
    parser.add_argument("--max-effective-stress-nrmse", type=float)
    parser.add_argument("--max-pore-pressure-difference-kpa", type=float)
    parser.add_argument("--max-total-radial-error-kpa", type=float)
    args = parser.parse_args()

    hercules = read_hercules(args.hercules_csv)
    source, opensees = read_opensees(args.opensees_csv)
    if len(hercules) != len(opensees):
        raise SystemExit(
            f"history length mismatch: Hercules={len(hercules)} "
            f"OpenSees={len(opensees)}"
        )

    phase_indices = [index for index, row in enumerate(opensees)
                     if row["phase"] == "coupled"]
    confinement_pa = 1000.0 * float(opensees[phase_indices[0]]["confinement_kpa"])
    fluid_bulk_pa = 1000.0 * float(opensees[0]["fluid_bulk_modulus_kpa"])
    porosity = float(opensees[0]["porosity"])
    storage_modulus_pa = fluid_bulk_pa / porosity
    initial_ev = sum(float(opensees[0][field]) for field in STRAIN_FIELDS[:3])
    initial_op = 1000.0 * float(opensees[0]["pore_pressure_kpa"])
    initial_hp = hercules[0]["pore_pressure_pa"]

    maximum_strain_difference = 0.0
    effective_errors: list[float] = []
    effective_references: list[float] = []
    total_errors: list[float] = []
    total_references: list[float] = []
    pressure_errors: list[float] = []
    pressure_references: list[float] = []
    state_errors: list[float] = []
    detailed: list[dict[str, float | str]] = []

    for index, (hrow, orow) in enumerate(zip(hercules, opensees)):
        if int(hrow["step"]) != index:
            raise SystemExit(f"nonsequential Hercules step at row {index}")
        output: dict[str, float | str] = {
            "step": float(index),
            "opensees_step": float(orow["step"]),
            "phase": orow["phase"],
        }
        for field in STRAIN_FIELDS:
            difference = hrow[field] - float(orow[field])
            maximum_strain_difference = max(
                maximum_strain_difference, abs(difference))
            output[field] = hrow[field]
            output[f"delta_{field}"] = difference

        heffective = tuple(hrow[f"sig_{component}"] for component in COMPONENTS)
        oeffective = tuple(
            1000.0 * float(orow[f"sigma_{component}_kpa"])
            for component in COMPONENTS
        )
        htotal = tuple(
            hrow[f"total_sig_{component}"] for component in COMPONENTS
        )
        ototal = tuple(
            1000.0 * float(orow[f"total_sigma_{component}_kpa"])
            for component in COMPONENTS
        )
        hp_eff, hq = invariants(heffective)
        op_eff, oq = invariants(oeffective)
        output.update({
            "hercules_p_eff_pa": hp_eff,
            "opensees_p_eff_pa": op_eff,
            "hercules_q_pa": hq,
            "opensees_q_pa": oq,
        })
        for component, hv, ov, htv, otv in zip(
                COMPONENTS, heffective, oeffective, htotal, ototal):
            effective_errors.append(hv - ov)
            effective_references.append(ov)
            total_errors.append(htv - otv)
            total_references.append(otv)
            output[f"hercules_sig_{component}_pa"] = hv
            output[f"opensees_sig_{component}_pa"] = ov
            output[f"delta_sig_{component}_pa"] = hv - ov
            output[f"hercules_total_sig_{component}_pa"] = htv
            output[f"opensees_total_sig_{component}_pa"] = otv
            output[f"delta_total_sig_{component}_pa"] = htv - otv

        hpore = hrow["pore_pressure_pa"]
        opore = 1000.0 * float(orow["pore_pressure_kpa"])
        pressure_errors.append(hpore - opore)
        pressure_references.append(opore)
        output["hercules_pore_pressure_pa"] = hpore
        output["opensees_pore_pressure_pa"] = opore
        output["delta_pore_pressure_pa"] = hpore - opore
        output["hercules_ru"] = hpore / confinement_pa
        output["opensees_ru"] = opore / confinement_pa

        for hfield, ofield in STATE_PAIRS:
            difference = hrow[hfield] - float(orow[ofield])
            state_errors.append(difference)
            output[f"delta_{hfield}"] = difference
        detailed.append(output)

    coupled_h_radial = [
        -hercules[index][f"total_sig_{component}"]
        for index in phase_indices for component in ("xx", "yy")
    ]
    coupled_o_radial = [
        -1000.0 * float(opensees[index][f"total_sigma_{component}_kpa"])
        for index in phase_indices for component in ("xx", "yy")
    ]
    h_storage_residual = []
    o_storage_residual = []
    for index in phase_indices:
        hev = sum(hercules[index][field] for field in STRAIN_FIELDS[:3])
        oev = sum(float(opensees[index][field]) for field in STRAIN_FIELDS[:3])
        h_storage_residual.append(
            hercules[index]["pore_pressure_pa"] - initial_hp
            + storage_modulus_pa*(hev - initial_ev)
        )
        o_storage_residual.append(
            1000.0*float(opensees[index]["pore_pressure_kpa"]) - initial_op
            + storage_modulus_pa*(oev - initial_ev)
        )

    effective_nrmse = normalized_rmse(effective_errors, effective_references)
    pressure_nrmse = normalized_rmse(pressure_errors, pressure_references)
    total_nrmse = normalized_rmse(total_errors, total_references)
    maximum_pressure_difference_kpa = max(
        abs(value) for value in pressure_errors
    ) / 1000.0
    maximum_h_radial_error_kpa = max(
        abs(value - confinement_pa) for value in coupled_h_radial
    ) / 1000.0
    maximum_o_radial_error_kpa = max(
        abs(value - confinement_pa) for value in coupled_o_radial
    ) / 1000.0
    maximum_radial_difference_kpa = max(
        abs(left - right)
        for left, right in zip(coupled_h_radial, coupled_o_radial)
    ) / 1000.0
    maximum_ru_difference = max(
        abs(hercules[index]["pore_pressure_pa"]
            - 1000.0*float(opensees[index]["pore_pressure_kpa"]))
        / confinement_pa
        for index in phase_indices
    )

    summary: dict[str, object] = {
        "test_id": TEST_ID,
        "steps": len(hercules),
        "opensees_source_steps": len(source),
        "opensees_initial_replay_step": int(opensees[0]["step"]),
        "coupled_steps": len(phase_indices),
        "maximum_strain_difference": maximum_strain_difference,
        "maximum_effective_stress_difference_pa": max(
            abs(value) for value in effective_errors
        ),
        "normalized_effective_stress_rmse": effective_nrmse,
        "maximum_total_stress_difference_pa": max(
            abs(value) for value in total_errors
        ),
        "normalized_total_stress_rmse": total_nrmse,
        "maximum_pore_pressure_difference_kpa": maximum_pressure_difference_kpa,
        "normalized_pore_pressure_rmse": pressure_nrmse,
        "maximum_ru_difference": maximum_ru_difference,
        "maximum_alpha_fabric_void_ratio_difference": max(
            abs(value) for value in state_errors
        ),
        "target_total_radial_confinement_kpa": confinement_pa / 1000.0,
        "maximum_total_radial_error_hercules_kpa": maximum_h_radial_error_kpa,
        "maximum_total_radial_error_opensees_kpa": maximum_o_radial_error_kpa,
        "maximum_hercules_opensees_total_radial_difference_kpa": (
            maximum_radial_difference_kpa
        ),
        "maximum_sealed_storage_residual_hercules_pa": max(
            abs(value) for value in h_storage_residual
        ),
        "maximum_sealed_storage_residual_opensees_pa": max(
            abs(value) for value in o_storage_residual
        ),
        "maximum_ru_hercules": max(
            hercules[index]["pore_pressure_pa"] / confinement_pa
            for index in phase_indices
        ),
        "maximum_ru_opensees": max(
            1000.0*float(opensees[index]["pore_pressure_kpa"])
            / confinement_pa for index in phase_indices
        ),
        "interpretation": (
            "The OpenSees SSPbrickUP path and Hercules BrickUP replay use "
            "zero permeability, unit Biot coupling, row-summed pressure "
            "storage, and the same affine strain history. This validates the "
            "uniform sealed "
            "storage/Q^T du mode plus the effective-stress constitutive law; "
            "it does not validate Darcy diffusion or a spatial pressure mode."
        ),
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    detailed_path = args.output_dir / "comparison.csv"
    fieldnames: list[str] = []
    for row in detailed:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with detailed_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(detailed)
    summary_path = args.output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n",
                            encoding="utf-8")

    failures: list[str] = []
    if (args.max_effective_stress_nrmse is not None
            and effective_nrmse > args.max_effective_stress_nrmse):
        failures.append(
            f"effective-stress NRMSE {effective_nrmse:.6e} exceeds "
            f"{args.max_effective_stress_nrmse:.6e}"
        )
    if (args.max_pore_pressure_difference_kpa is not None
            and maximum_pressure_difference_kpa
            > args.max_pore_pressure_difference_kpa):
        failures.append(
            f"pore-pressure difference {maximum_pressure_difference_kpa:.6e} "
            f"kPa exceeds {args.max_pore_pressure_difference_kpa:.6e} kPa"
        )
    if (args.max_total_radial_error_kpa is not None
            and maximum_h_radial_error_kpa > args.max_total_radial_error_kpa):
        failures.append(
            f"Hercules total radial error {maximum_h_radial_error_kpa:.6e} "
            f"kPa exceeds {args.max_total_radial_error_kpa:.6e} kPa"
        )

    print(f"Coupled canonical test: {TEST_ID}")
    print(f"Hercules/OpenSees states: {len(hercules)}")
    print(f"Maximum strain difference: {maximum_strain_difference:.6e}")
    print(f"Effective-stress NRMSE: {effective_nrmse:.6e}")
    print(f"Total-stress NRMSE: {total_nrmse:.6e}")
    print("Maximum pore-pressure difference: "
          f"{maximum_pressure_difference_kpa:.6e} kPa")
    print("Maximum ru (Hercules/OpenSees): "
          f"{summary['maximum_ru_hercules']:.6e} / "
          f"{summary['maximum_ru_opensees']:.6e}")
    print("Maximum total radial error (Hercules/OpenSees): "
          f"{maximum_h_radial_error_kpa:.6e} / "
          f"{maximum_o_radial_error_kpa:.6e} kPa")
    print("Maximum Hercules/OpenSees total radial difference: "
          f"{maximum_radial_difference_kpa:.6e} kPa")
    print(f"Detailed comparison: {detailed_path.resolve()}")
    print(f"Summary: {summary_path.resolve()}")
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
