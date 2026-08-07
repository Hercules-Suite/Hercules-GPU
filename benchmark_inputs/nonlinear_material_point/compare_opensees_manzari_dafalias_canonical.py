#!/usr/bin/env python3
"""Compare canonical Hercules/OpenSees Manzari--Dafalias histories."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


STRAIN_FIELDS = ("eps_xx", "eps_yy", "eps_zz", "eps_xy", "eps_yz", "eps_xz")
STRESS_PAIRS = tuple(
    (f"sig_{component}", f"sigma_{component}_kpa")
    for component in ("xx", "yy", "zz", "xy", "yz", "xz")
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
TEST_MODES = {
    "drained": ("md_drained_triaxial", "monotonic", "2"),
    "undrained": ("md_undrained_cyclic_shear", "cyclic", "1"),
}


def read_hercules(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return [
            {key: float(value) for key, value in row.items() if key != "status"}
            for row in csv.DictReader(stream)
        ]


def read_opensees(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


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
        dxx * dxx + dyy * dyy + dzz * dzz
        + 2.0 * (sxy * sxy + syz * syz + sxz * sxz)
    ))
    return pressure, q


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def reversal_count(values: list[float]) -> int:
    scale = max((abs(value) for value in values), default=0.0)
    cutoff = max(scale * 1.0e-10, 1.0e-12)
    signs = [1 if value > 0.0 else -1 for value in values if abs(value) > cutoff]
    return sum(left != right for left, right in zip(signs, signs[1:]))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=sorted(TEST_MODES), required=True)
    parser.add_argument("--hercules-csv", type=Path, required=True)
    parser.add_argument("--opensees-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--max-stress-nrmse", type=float)
    parser.add_argument("--max-radial-confinement-error-kpa", type=float)
    parser.add_argument("--max-ru-difference", type=float)
    args = parser.parse_args()

    test_id, comparison_phase, driver_version = TEST_MODES[args.mode]
    output_dir = args.output_dir or Path(f"opensees_{test_id}_comparison")
    hercules = read_hercules(args.hercules_csv)
    opensees = read_opensees(args.opensees_csv)
    if len(hercules) != len(opensees):
        raise SystemExit(
            f"history length mismatch: Hercules={len(hercules)} "
            f"OpenSees={len(opensees)}"
        )
    if not hercules:
        raise SystemExit("comparison histories are empty")

    output_dir.mkdir(parents=True, exist_ok=True)
    stress_errors: list[float] = []
    stress_references: list[float] = []
    state_errors: list[float] = []
    max_strain_error = 0.0
    detailed: list[dict[str, float | str]] = []
    phase_indices: list[int] = []

    for index, (hrow, orow) in enumerate(zip(hercules, opensees)):
        if (orow.get("test_id") != test_id
                or orow.get("driver_version") != driver_version):
            raise SystemExit(
                f"OpenSees CSV is not {test_id} driver version "
                f"{driver_version} at row {index}"
            )
        if int(orow["step"]) != index or int(hrow["step"]) != index:
            raise SystemExit(f"nonsequential step at row {index}")
        if orow["phase"] == comparison_phase:
            phase_indices.append(index)

        output: dict[str, float | str] = {
            "step": float(index),
            "phase": orow["phase"],
        }
        for field in STRAIN_FIELDS:
            error = hrow[field] - float(orow[field])
            max_strain_error = max(max_strain_error, abs(error))
            output[field] = hrow[field]
            output[f"delta_{field}"] = error

        hstress = tuple(hrow[hfield] for hfield, _ in STRESS_PAIRS)
        ostress = tuple(1000.0 * float(orow[ofield]) for _, ofield in STRESS_PAIRS)
        hp, hq = invariants(hstress)
        op, oq = invariants(ostress)
        output.update({
            "hercules_p_eff_pa": hp,
            "opensees_p_eff_pa": op,
            "delta_p_eff_pa": hp - op,
            "hercules_q_pa": hq,
            "opensees_q_pa": oq,
            "delta_q_pa": hq - oq,
        })
        for (hfield, _), hvalue, ovalue in zip(STRESS_PAIRS, hstress, ostress):
            error = hvalue - ovalue
            stress_errors.append(error)
            stress_references.append(ovalue)
            output[f"hercules_{hfield}_pa"] = hvalue
            output[f"opensees_{hfield}_pa"] = ovalue
            output[f"delta_{hfield}_pa"] = error
        for hfield, ofield in STATE_PAIRS:
            error = hrow[hfield] - float(orow[ofield])
            state_errors.append(error)
            output[f"delta_{hfield}"] = error
        detailed.append(output)

    if not phase_indices:
        raise SystemExit(f"OpenSees CSV has no {comparison_phase!r} states")

    stress_nrmse = normalized_rmse(stress_errors, stress_references)
    summary: dict[str, object] = {
        "test_id": test_id,
        "steps": len(hercules),
        "comparison_phase": comparison_phase,
        "comparison_phase_steps": len(phase_indices),
        "maximum_strain_difference": max_strain_error,
        "maximum_stress_difference_pa": max(abs(value) for value in stress_errors),
        "normalized_stress_rmse": stress_nrmse,
        "maximum_alpha_fabric_void_ratio_difference": max(
            abs(value) for value in state_errors
        ),
    }

    if args.mode == "drained":
        confinement_pa = 1000.0 * float(opensees[phase_indices[0]]["confinement_kpa"])
        hradial: list[float] = []
        oradial: list[float] = []
        heta: list[float] = []
        oeta: list[float] = []
        for index in phase_indices:
            hrow = hercules[index]
            orow = opensees[index]
            hradial.append(-(hrow["sig_xx"] + hrow["sig_yy"]) / 2.0)
            oradial.append(-500.0 * (
                float(orow["sigma_xx_kpa"]) + float(orow["sigma_yy_kpa"])
            ))
            hp = float(detailed[index]["hercules_p_eff_pa"])
            hq = float(detailed[index]["hercules_q_pa"])
            op = float(detailed[index]["opensees_p_eff_pa"])
            oq = float(detailed[index]["opensees_q_pa"])
            heta.append(hq / hp if hp > 0.0 else math.nan)
            oeta.append(oq / op if op > 0.0 else math.nan)
        tail_count = max(5, len(phase_indices) // 10)
        base_index = phase_indices[0] - 1 if phase_indices[0] > 0 else phase_indices[0]
        base_ev = sum(hercules[base_index][field] for field in ("eps_xx", "eps_yy", "eps_zz"))
        delta_ev = [
            sum(hercules[index][field] for field in ("eps_xx", "eps_yy", "eps_zz"))
            - base_ev
            for index in phase_indices
        ]
        radial_error_kpa = max(abs(value - confinement_pa) for value in hradial) / 1000.0
        summary["drained_metrics"] = {
            "target_radial_confinement_kpa": confinement_pa / 1000.0,
            "maximum_opensees_radial_confinement_error_kpa": (
                max(abs(value - confinement_pa) for value in oradial) / 1000.0
            ),
            "maximum_hercules_radial_confinement_error_kpa": radial_error_kpa,
            "peak_q_over_p_hercules": max(heta),
            "peak_q_over_p_opensees": max(oeta),
            "tail_mean_q_over_p_hercules": mean(heta[-tail_count:]),
            "tail_mean_q_over_p_opensees": mean(oeta[-tail_count:]),
            "final_q_over_p_difference": heta[-1] - oeta[-1],
            "maximum_incremental_volumetric_contraction": max(0.0, -min(delta_ev)),
            "maximum_incremental_volumetric_dilation": max(0.0, max(delta_ev)),
            "interpretation": (
                "OpenSees solves the constant-confinement lateral strain. "
                "Hercules replays that strain; its radial-stress drift therefore "
                "measures disagreement in contraction/dilation tendency."
            ),
        }
    else:
        first = phase_indices[0]
        hp0 = float(detailed[first]["hercules_p_eff_pa"])
        op0 = float(detailed[first]["opensees_p_eff_pa"])
        if hp0 <= 0.0 or op0 <= 0.0:
            raise SystemExit("initial cyclic effective pressure must be positive")
        hru: list[float] = []
        oru: list[float] = []
        hshear: list[float] = []
        oshear: list[float] = []
        base_ev = sum(
            hercules[first][field] for field in ("eps_xx", "eps_yy", "eps_zz")
        )
        volume_drift: list[float] = []
        max_reported_ru_error = 0.0
        for index in phase_indices:
            hp = float(detailed[index]["hercules_p_eff_pa"])
            op = float(detailed[index]["opensees_p_eff_pa"])
            hru.append((hp0 - hp) / hp0)
            oru.append((op0 - op) / op0)
            hshear.append(hercules[index]["sig_xy"])
            oshear.append(1000.0 * float(opensees[index]["sigma_xy_kpa"]))
            volume_drift.append(abs(
                sum(hercules[index][field]
                    for field in ("eps_xx", "eps_yy", "eps_zz")) - base_ev
            ))
            max_reported_ru_error = max(
                max_reported_ru_error,
                abs(oru[-1] - float(opensees[index]["ru"])),
            )
            detailed[index]["hercules_excess_pore_pressure_pa"] = hp0 - hp
            detailed[index]["opensees_excess_pore_pressure_pa"] = op0 - op
            detailed[index]["hercules_ru"] = hru[-1]
            detailed[index]["opensees_ru"] = oru[-1]
            detailed[index]["delta_ru"] = hru[-1] - oru[-1]
        ru_difference = max(abs(left - right) for left, right in zip(hru, oru))
        summary["undrained_metrics"] = {
            "initial_effective_pressure_hercules_kpa": hp0 / 1000.0,
            "initial_effective_pressure_opensees_kpa": op0 / 1000.0,
            "maximum_ru_hercules": max(hru),
            "maximum_ru_opensees": max(oru),
            "maximum_ru_difference": ru_difference,
            "maximum_cyclic_volumetric_strain_drift": max(volume_drift),
            "minimum_effective_pressure_ratio_hercules": 1.0 - max(hru),
            "minimum_effective_pressure_ratio_opensees": 1.0 - max(oru),
            "hercules_shear_stress_reversals": reversal_count(hshear),
            "opensees_shear_stress_reversals": reversal_count(oshear),
            "maximum_opensees_reported_ru_consistency_error": max_reported_ru_error,
            "interpretation": (
                "For the zero-volume material path, excess pore pressure is "
                "inferred as p'_0-p'. This characterizes the effective-stress "
                "law; it does not exercise Hercules pore-pressure DOFs, Darcy "
                "flow, or solid-fluid coupling."
            ),
        }

    detailed_path = output_dir / "comparison.csv"
    fieldnames: list[str] = []
    for row in detailed:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with detailed_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(detailed)

    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    failures: list[str] = []
    if args.max_stress_nrmse is not None and stress_nrmse > args.max_stress_nrmse:
        failures.append(
            f"stress NRMSE {stress_nrmse:.6e} exceeds {args.max_stress_nrmse:.6e}"
        )
    if args.mode == "drained" and args.max_radial_confinement_error_kpa is not None:
        measured = float(summary["drained_metrics"]["maximum_hercules_radial_confinement_error_kpa"])  # type: ignore[index]
        if measured > args.max_radial_confinement_error_kpa:
            failures.append(
                f"radial-confinement error {measured:.6e} kPa exceeds "
                f"{args.max_radial_confinement_error_kpa:.6e} kPa"
            )
    if args.mode == "undrained" and args.max_ru_difference is not None:
        measured = float(summary["undrained_metrics"]["maximum_ru_difference"])  # type: ignore[index]
        if measured > args.max_ru_difference:
            failures.append(
                f"maximum ru difference {measured:.6e} exceeds "
                f"{args.max_ru_difference:.6e}"
            )

    print(f"Canonical test: {test_id}")
    print(f"Hercules/OpenSees states: {len(hercules)}")
    print(f"Maximum strain difference: {max_strain_error:.6e}")
    print(f"Maximum stress difference: {summary['maximum_stress_difference_pa']:.6e} Pa")
    print(f"Normalized stress RMSE: {stress_nrmse:.6e}")
    if args.mode == "drained":
        metrics = summary["drained_metrics"]
        print("Hercules maximum radial-confinement error: "
              f"{metrics['maximum_hercules_radial_confinement_error_kpa']:.6e} kPa")  # type: ignore[index]
        print("Tail mean q/p (Hercules/OpenSees): "
              f"{metrics['tail_mean_q_over_p_hercules']:.6e} / "  # type: ignore[index]
              f"{metrics['tail_mean_q_over_p_opensees']:.6e}")  # type: ignore[index]
    else:
        metrics = summary["undrained_metrics"]
        print("Maximum ru (Hercules/OpenSees): "
              f"{metrics['maximum_ru_hercules']:.6e} / "  # type: ignore[index]
              f"{metrics['maximum_ru_opensees']:.6e}")  # type: ignore[index]
        print(f"Maximum ru difference: {metrics['maximum_ru_difference']:.6e}")  # type: ignore[index]
    print(f"Detailed comparison: {detailed_path.resolve()}")
    print(f"Summary: {summary_path.resolve()}")
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
