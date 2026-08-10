#!/usr/bin/env python3
"""Compare the synthetic Hercules PDMY03 history with OpenSees 3.8."""

import argparse
import csv
import json
import math
from pathlib import Path


STRESS_FIELDS = ("sigma_xx", "sigma_yy", "sigma_zz", "sigma_xy", "sigma_yz", "sigma_xz")
HERCULES_STRESS_FIELDS = ("sig_xx", "sig_yy", "sig_zz", "sig_xy", "sig_yz", "sig_xz")


def rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hercules-csv", type=Path, required=True)
    parser.add_argument("--opensees-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("opensees_pdmy03_comparison"))
    parser.add_argument("--max-stress-nrmse", type=float, default=1.0e-8)
    parser.add_argument("--max-shear-strain-difference", type=float, default=1.0e-12)
    args = parser.parse_args()

    hercules = rows(args.hercules_csv)
    opensees = rows(args.opensees_csv)
    if len(hercules) != len(opensees) or not hercules:
        raise SystemExit("Hercules and OpenSees histories have different nonzero lengths")
    if any(row.get("test_id") != "pdmy03_cyclic_shear" or row.get("driver_version") != "1" for row in opensees):
        raise SystemExit("OpenSees CSV was not produced by the PDMY03 cyclic-shear driver version 1")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    maximum_strain = 0.0
    maximum_stress = 0.0
    stress_square = 0.0
    stress_scale = 0.0
    detailed = []
    for index, (hrow, orow) in enumerate(zip(hercules, opensees)):
        hstrain = float(hrow["eps_xy"])
        ostrain = float(orow["epsilon_xy"])
        maximum_strain = max(maximum_strain, abs(hstrain-ostrain))
        output = {"step": index, "epsilon_xy": hstrain, "gamma_xy": 2.0*hstrain}
        for hfield, ofield in zip(HERCULES_STRESS_FIELDS, STRESS_FIELDS):
            hvalue = float(hrow[hfield])
            ovalue = float(orow[ofield])
            difference = hvalue-ovalue
            output[f"hercules_{ofield}"] = hvalue
            output[f"opensees_{ofield}"] = ovalue
            output[f"difference_{ofield}"] = difference
            maximum_stress = max(maximum_stress, abs(difference))
            stress_square += difference*difference
            stress_scale = max(stress_scale, abs(hvalue), abs(ovalue))
        detailed.append(output)

    stress_nrmse = math.sqrt(stress_square/(len(detailed)*len(STRESS_FIELDS)))/max(stress_scale, 1.0e-300)
    detail_path = args.output_dir/"comparison.csv"
    with detail_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=detailed[0].keys())
        writer.writeheader()
        writer.writerows(detailed)
    summary = {
        "test_id": "pdmy03_cyclic_shear",
        "states": len(detailed),
        "maximum_shear_strain_difference": maximum_strain,
        "maximum_stress_difference_pa": maximum_stress,
        "stress_nrmse": stress_nrmse,
        "maximum_stress_nrmse": args.max_stress_nrmse,
        "maximum_allowed_shear_strain_difference": args.max_shear_strain_difference,
    }
    summary_path = args.output_dir/"summary.json"
    summary_path.write_text(json.dumps(summary, indent=2)+"\n")

    print("PDMY03 Hercules/OpenSees states:", len(detailed))
    print(f"Maximum shear-strain difference: {maximum_strain:.6e}")
    print(f"Maximum stress difference: {maximum_stress:.6e} Pa")
    print(f"Normalized stress RMSE: {stress_nrmse:.6e}")
    print("Detailed comparison:", detail_path.resolve())
    print("Summary:", summary_path.resolve())
    if maximum_strain > args.max_shear_strain_difference or stress_nrmse > args.max_stress_nrmse:
        raise SystemExit("FAIL: PDMY03 OpenSees characterization exceeds its configured threshold")
    print("PASS: PDMY03 Hercules/OpenSees characterization agrees")


if __name__ == "__main__":
    main()
