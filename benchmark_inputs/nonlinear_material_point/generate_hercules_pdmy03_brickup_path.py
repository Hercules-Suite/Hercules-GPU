#!/usr/bin/env python3
"""Generate a Hercules PDMY03 BrickUP material-point replay from OpenSees."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


SUPPORTED_TESTS = {
    "pdmy03_sspbrickup_undrained_triaxial": "pdmy03_brickup_undrained_triaxial",
    "pdmy03_sspbrickup_undrained_cyclic_shear":
        "pdmy03_brickup_undrained_cyclic_shear",
}
DRIVER_VERSION = "2"
STRAIN_FIELDS = ("eps_xx", "eps_yy", "eps_zz", "eps_xy", "eps_yz", "eps_xz")
STRESS_FIELDS = (
    "sigma_xx_kpa", "sigma_yy_kpa", "sigma_zz_kpa",
    "sigma_xy_kpa", "sigma_yz_kpa", "sigma_xz_kpa",
)


def read_rows(
    path: Path,
) -> tuple[str, list[dict[str, str]], list[dict[str, str]]]:
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
        if int(row["step"]) != index:
            raise SystemExit(f"nonsequential OpenSees step at row {index}")
    if source[-1].get("phase") != "complete":
        raise SystemExit("OpenSees PDMY03 SSPbrickUP history is incomplete")
    base = max(index for index, row in enumerate(source)
               if row.get("phase") == "consolidation")
    replay = source[base:]
    if not any(row.get("phase") == "coupled" for row in replay):
        raise SystemExit("OpenSees PDMY03 SSPbrickUP history has no coupled states")
    return test_id, source, replay


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--opensees-csv", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--output-prefix")
    parser.add_argument("--maximum-substeps", type=int, default=400)
    args = parser.parse_args()

    test_id, source, rows = read_rows(args.opensees_csv)
    initial = rows[0]
    if args.maximum_substeps <= 0:
        raise SystemExit("--maximum-substeps must be positive")
    radial = [1000.0 * float(initial[field]) for field in STRESS_FIELDS[:3]]
    if max(radial) - min(radial) > 1.0e-4:
        raise SystemExit("OpenSees stage-boundary effective stress is not isotropic")

    prefix = args.output_prefix or SUPPORTED_TESTS[test_id]
    fluid_bulk_pa = 1000.0 * float(initial["fluid_bulk_modulus_kpa"])
    zeros = "0.0 0.0 0.0 0.0 0.0 0.0"
    lines = [
        "# Generated coupled PDMY03/BrickUP replay from OpenSees SSPbrickUP.",
        f"# Source: {args.opensees_csv}",
        "# Constants are numerical characterization inputs, not a calibration.",
        "material_model = PressureDependMultiYield03",
        "material_point_delta_time_sec = 1.0",
        f"material_point_steps = {len(rows)}",
        f"no_substeps = {args.maximum_substeps}",
        "error_tolerance = 1.0e-1",
        "backbone_errTol = 1.0e-4",
        "tension_cutoff = no",
        "material_point_compare_rtol = 1.0e-6",
        "material_point_compare_atol = 1.0e-6",
        f"material_point_output_prefix = {prefix}",
        "material_point_brickup = yes",
        f"material_point_brickup_porosity = {float(initial['porosity']):.17e}",
        f"material_point_brickup_fluid_bulk_modulus_pa = {fluid_bulk_pa:.17e}",
        f"material_point_brickup_biot_alpha = {float(initial['biot_alpha']):.17e}",
        ("material_point_initial_pore_pressure_pa = "
         f"{1000.0*float(initial['pore_pressure_kpa']):.17e}"),
        "",
        "# Gref Kref phi peak_gamma ref_pressure pressure_exponent PT_angle mType",
        "# ca cb cc cd ce da db dc surfaces liquefaction1 liquefaction2 Patm cohesion",
        "material_point_constants =",
        "9.0e7 2.2e8 32.0 0.10 1.013e5 0.50 26.0 0",
        "0.067 0.50 0.23 1.0 0.0 0.06 3.0 0.27",
        "20 1.0e4 0.01 1.013e5 1.73e3",
        "",
        "material_point_initial_strain =",
        " ".join(f"{float(initial[field]):.17e}" for field in STRAIN_FIELDS),
        "material_point_initial_stress =",
        " ".join(f"{1000.0*float(initial[field]):.17e}" for field in STRESS_FIELDS),
        "material_point_initial_alpha =",
        zeros,
        "material_point_initial_alpha_in =",
        zeros,
        "material_point_initial_fabric =",
        zeros,
        "",
        "# Hercules tensor strain: exx eyy ezz exy eyz exz",
        "material_point_strain_history =",
    ]
    for row in rows:
        lines.append(" ".join(f"{float(row[field]):.17e}" for field in STRAIN_FIELDS))
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"Coupled canonical test: {test_id}")
    print(f"OpenSees source states: {len(source)}")
    print(f"BrickUP replay states: {len(rows)}")
    print(f"Initial OpenSees step: {initial['step']}")
    print(f"Hercules input: {args.output.resolve()}")
    print(f"Expected output names: {prefix}.cpu.csv and {prefix}.gpu.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
