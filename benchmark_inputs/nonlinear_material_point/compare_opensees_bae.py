#!/usr/bin/env python3
"""Run and compare Hercules/OpenSees Borja--Amies material histories."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


STRESS_FIELDS = ("sig_xx", "sig_yy", "sig_zz", "sig_xy", "sig_yz", "sig_xz")
OPEN_STRESS_FIELDS = (
    "sigma_xx",
    "sigma_yy",
    "sigma_zz",
    "sigma_xy",
    "sigma_yz",
    "sigma_xz",
)
ASSIGNMENT = re.compile(r"(?m)^\s*([A-Za-z][A-Za-z0-9_]*)\s*=\s*")


@dataclass(frozen=True)
class Fixture:
    path: Path
    steps: int
    output_prefix: str
    strains: tuple[tuple[float, float, float, float, float, float], ...]
    bulk: float
    shear: float
    su: float
    psi0: float
    exponent: float


def parse_assignments(path: Path) -> dict[str, str]:
    text = re.sub(r"#.*", "", path.read_text(encoding="utf-8"))
    matches = list(ASSIGNMENT.finditer(text))
    values: dict[str, str] = {}
    for index, match in enumerate(matches):
        end = matches[index + 1].start() if index + 1 < len(matches) else len(text)
        values[match.group(1)] = text[match.end():end].strip()
    return values


def parse_float_list(value: str) -> list[float]:
    return [float(token) for token in value.split()]


def load_fixture(path: Path) -> Fixture:
    values = parse_assignments(path)
    if values.get("material_model", "").lower() != "vonmises_bae":
        raise ValueError("OpenSees comparison requires material_model=vonMises_baE")
    constants = parse_float_list(values["material_point_constants"])
    if len(constants) != 24:
        raise ValueError(f"expected 24 material constants, found {len(constants)}")
    steps = int(values["material_point_steps"])
    strain_values = parse_float_list(values["material_point_strain_history"])
    if len(strain_values) != 6 * steps:
        raise ValueError(
            f"expected {6 * steps} strain values for {steps} steps, "
            f"found {len(strain_values)}"
        )
    strains = tuple(
        tuple(strain_values[index:index + 6])
        for index in range(0, len(strain_values), 6)
    )
    for step, strain in enumerate(strains):
        if any(value != 0.0 for index, value in enumerate(strain) if index != 3):
            raise ValueError(
                f"step {step} is not pure xy simple shear; the OpenSees brick "
                "driver currently supports only epsilon_xy histories"
            )
    lam, shear = constants[0], constants[1]
    return Fixture(
        path=path,
        steps=steps,
        output_prefix=values.get("material_point_output_prefix", path.stem),
        strains=strains,
        bulk=lam + 2.0 * shear / 3.0,
        shear=shear,
        su=constants[5],
        psi0=constants[10],
        exponent=constants[11],
    )


def tcl_float(value: float) -> str:
    return f"{value:.17e}"


def tcl_path(path: Path) -> str:
    return "{" + str(path).replace("\\", "/").replace("}", "\\}") + "}"


def write_opensees_driver(
    fixture: Fixture,
    path: Path,
    csv_path: Path,
    beta: float,
    k0: float,
) -> None:
    hardening_scale = fixture.psi0 * fixture.shear
    lines = [
        "wipe",
        "model BasicBuilder -ndm 3 -ndf 3",
        (
            "nDMaterial MultiaxialCyclicPlasticity 1 "
            f"0.0 {tcl_float(fixture.bulk)} {tcl_float(fixture.shear)} "
            f"{tcl_float(fixture.su)} 0.0 {tcl_float(hardening_scale)} "
            f"{tcl_float(fixture.exponent)} {tcl_float(beta)} {tcl_float(k0)}"
        ),
        "node 1 0.0 0.0 0.0",
        "node 2 1.0 0.0 0.0",
        "node 3 1.0 1.0 0.0",
        "node 4 0.0 1.0 0.0",
        "node 5 0.0 0.0 1.0",
        "node 6 1.0 0.0 1.0",
        "node 7 1.0 1.0 1.0",
        "node 8 0.0 1.0 1.0",
        "fix 1 1 1 1",
        "fix 2 1 1 1",
        "fix 5 1 1 1",
        "fix 6 1 1 1",
        "fix 3 0 1 1",
        "fix 4 0 1 1",
        "fix 7 0 1 1",
        "fix 8 0 1 1",
        "equalDOF 3 4 1",
        "equalDOF 3 7 1",
        "equalDOF 3 8 1",
        "element stdBrick 1 1 2 3 4 5 6 7 8 1 0.0 0.0 0.0",
        (
            "timeSeries Path 1 -time {"
            + " ".join(str(step) for step in range(fixture.steps + 1))
            + "} -values {0.0 "
            + " ".join(tcl_float(2.0 * strain[3]) for strain in fixture.strains)
            + "} -useLast"
        ),
        "pattern Plain 1 1 { sp 3 1 1.0 }",
        "constraints Penalty 1.0e18 1.0e18",
        "numberer RCM",
        "system BandGeneral",
        "test NormDispIncr 1.0e-12 50 0",
        "algorithm Newton",
        "integrator LoadControl 1.0",
        "analysis Static",
        f"set output [open {tcl_path(csv_path)} w]",
        (
            "puts $output \"step,epsilon_xy,gamma_xy,sigma_xx,sigma_yy,"
            "sigma_zz,sigma_xy,sigma_yz,sigma_xz,driver_version,domain_time\""
        ),
    ]
    for step, strain in enumerate(fixture.strains):
        gamma_xy = 2.0 * strain[3]
        lines.extend(
            [
                (
                    f"if {{[analyze 1] != 0}} {{ puts stderr \"OpenSees analysis "
                    f"failed at material step {step}\"; close $output; exit 42 }}"
                ),
                "set actualGamma [nodeDisp 3 1]",
                (
                    f"if {{abs($actualGamma - {tcl_float(gamma_xy)}) > 1.0e-12}} {{ "
                    f"puts stderr \"OpenSees displacement mismatch at material step {step}: "
                    f"target={tcl_float(gamma_xy)} actual=$actualGamma\"; "
                    "close $output; exit 44 }"
                ),
                "set stress [eleResponse 1 material 1 stress]",
                (
                    "if {[llength $stress] != 6} { puts stderr \"OpenSees brick "
                    "did not return six material stresses\"; close $output; exit 43 }"
                ),
                (
                    "puts $output [format \"%d,%.17e,%.17e,%.17e,%.17e,%.17e,"
                    "%.17e,%.17e,%.17e,%d,%.17e\" "
                    f"{step} [expr {{0.5 * $actualGamma}}] $actualGamma "
                    "[lindex $stress 0] [lindex $stress 1] [lindex $stress 2] "
                    "[lindex $stress 3] [lindex $stress 4] [lindex $stress 5] "
                    "2 [getTime]]"
                ),
            ]
        )
    lines.extend(["close $output", "wipe", "exit 0", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def run_command(command: list[str], cwd: Path, stdout_path: Path, stderr_path: Path) -> None:
    completed = subprocess.run(command, cwd=cwd, text=True, capture_output=True, check=False)
    stdout_path.write_text(completed.stdout, encoding="utf-8")
    stderr_path.write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"command failed with status {completed.returncode}: {' '.join(command)}; "
            f"see {stderr_path}"
        )


def read_csv(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = []
        for row in csv.DictReader(stream):
            rows.append({key: float(value) for key, value in row.items() if key != "status"})
    return rows


def loop_work(rows: Iterable[dict[str, float]], gamma_field: str, stress_field: str) -> float:
    gamma_previous = 0.0
    stress_previous = 0.0
    work = 0.0
    for row in rows:
        gamma = row[gamma_field]
        stress = row[stress_field]
        work += 0.5 * (stress_previous + stress) * (gamma - gamma_previous)
        gamma_previous, stress_previous = gamma, stress
    return work


def compare_histories(
    fixture: Fixture,
    hercules_rows: list[dict[str, float]],
    opensees_rows: list[dict[str, float]],
    comparison_path: Path,
) -> dict[str, object]:
    if len(hercules_rows) != fixture.steps or len(opensees_rows) != fixture.steps:
        raise ValueError(
            f"history length mismatch: fixture={fixture.steps}, "
            f"Hercules={len(hercules_rows)}, OpenSees={len(opensees_rows)}"
        )
    for step, row in enumerate(opensees_rows):
        if row.get("driver_version") != 2.0:
            raise ValueError(
                "OpenSees CSV was not produced by prescribed-displacement driver "
                "version 2; rerun opensees_bae_synthetic.tcl"
            )
        if not math.isclose(
            row.get("domain_time", math.nan),
            float(step + 1),
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ):
            raise ValueError(f"OpenSees domain time mismatch at step {step}")
    field_metrics: dict[str, dict[str, float]] = {}
    with comparison_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ["step", "epsilon_xy", "gamma_xy"]
            + [
                item
                for field in STRESS_FIELDS
                for item in (
                    f"hercules_{field}",
                    f"opensees_{field}",
                    f"difference_{field}",
                )
            ]
        )
        for step, (hercules, opensees) in enumerate(zip(hercules_rows, opensees_rows)):
            expected_epsilon = fixture.strains[step][3]
            if not math.isclose(hercules["eps_xy"], expected_epsilon, rel_tol=0.0, abs_tol=1.0e-15):
                raise ValueError(f"Hercules epsilon_xy mismatch at step {step}")
            if not math.isclose(opensees["epsilon_xy"], expected_epsilon, rel_tol=0.0, abs_tol=1.0e-15):
                raise ValueError(f"OpenSees epsilon_xy mismatch at step {step}")
            if not math.isclose(
                opensees["gamma_xy"],
                2.0 * expected_epsilon,
                rel_tol=0.0,
                abs_tol=1.0e-15,
            ):
                raise ValueError(f"OpenSees gamma_xy mismatch at step {step}")
            output = [step, expected_epsilon, 2.0 * expected_epsilon]
            for hercules_field, opensees_field in zip(STRESS_FIELDS, OPEN_STRESS_FIELDS):
                left, right = hercules[hercules_field], opensees[opensees_field]
                output.extend((left, right, right - left))
            writer.writerow(output)

    for hercules_field, opensees_field in zip(STRESS_FIELDS, OPEN_STRESS_FIELDS):
        differences = [
            opensees[opensees_field] - hercules[hercules_field]
            for hercules, opensees in zip(hercules_rows, opensees_rows)
        ]
        scale = max(
            fixture.su,
            max(abs(row[hercules_field]) for row in hercules_rows),
            max(abs(row[opensees_field]) for row in opensees_rows),
            sys.float_info.min,
        )
        rms = math.sqrt(sum(value * value for value in differences) / len(differences))
        field_metrics[hercules_field] = {
            "maximum_absolute_difference": max(abs(value) for value in differences),
            "root_mean_square_difference": rms,
            "normalized_root_mean_square_difference": rms / scale,
            "normalization_scale": scale,
        }

    hercules_work = loop_work(hercules_rows, "eps_xy", "sig_xy") * 2.0
    opensees_work = loop_work(opensees_rows, "gamma_xy", "sigma_xy")
    work_scale = max(abs(hercules_work), abs(opensees_work), sys.float_info.min)
    return {
        "steps": fixture.steps,
        "stress_metrics": field_metrics,
        "shear_work": {
            "hercules": hercules_work,
            "opensees": opensees_work,
            "absolute_difference": abs(opensees_work - hercules_work),
            "relative_difference": abs(opensees_work - hercules_work) / work_scale,
        },
    }


def run_selftest(fixture_path: Path) -> int:
    fixture = load_fixture(fixture_path)
    assert math.isclose(fixture.bulk, 2.1333333333333333e8)
    assert math.isclose(fixture.psi0 * fixture.shear, 7.4648e7)
    assert any(fixture.strains[index] == fixture.strains[index - 1] for index in range(1, fixture.steps))
    standalone = fixture_path.with_name("opensees_bae_synthetic.tcl").read_text(
        encoding="utf-8"
    )
    gamma_match = re.search(r"set gammaHistory\s*\{(.*?)\}", standalone, re.DOTALL)
    assert gamma_match is not None
    standalone_gamma = tuple(float(value) for value in gamma_match.group(1).split())
    assert standalone_gamma == tuple(2.0 * strain[3] for strain in fixture.strains)
    for name, expected in {
        "bulkModulus": fixture.bulk,
        "shearModulus": fixture.shear,
        "undrainedShearStrength": fixture.su,
        "hardeningScale": fixture.psi0 * fixture.shear,
        "hardeningExponent": fixture.exponent,
    }.items():
        value_match = re.search(rf"(?m)^set {name}\s+(\S+)$", standalone)
        assert value_match is not None
        assert math.isclose(float(value_match.group(1)), expected)
    assert "set beta 1.0" in standalone
    assert "sigma_xy,sigma_yz,sigma_xz" in standalone
    assert "driver_version,domain_time" in standalone
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        driver = root / "driver.tcl"
        write_opensees_driver(fixture, driver, root / "output.csv", 1.0, 0.5)
        generated = driver.read_text(encoding="utf-8")
        assert "timeSeries Path 1" in generated
        assert "pattern Plain 1 1 { sp 3 1 1.0 }" in generated
        assert "constraints Penalty" in generated
        assert "integrator LoadControl 1.0" in generated
        assert "set actualGamma [nodeDisp 3 1]" in generated
        assert "integrator LoadControl 0.0" not in generated
        assert "integrator DisplacementControl" not in generated
        assert "element stdBrick" in generated
        assert "7.46480000000000000e+07" in generated
        hercules_rows = []
        opensees_rows = []
        for step, strain in enumerate(fixture.strains):
            epsilon_xy = strain[3]
            sigma_xy = 2.0 * fixture.shear * epsilon_xy
            hercules_rows.append(
                {
                    "step": float(step),
                    "eps_xy": epsilon_xy,
                    "sig_xx": 0.0,
                    "sig_yy": 0.0,
                    "sig_zz": 0.0,
                    "sig_xy": sigma_xy,
                    "sig_yz": 0.0,
                    "sig_xz": 0.0,
                }
            )
            opensees_rows.append(
                {
                    "step": float(step),
                    "epsilon_xy": epsilon_xy,
                    "gamma_xy": 2.0 * epsilon_xy,
                    "sigma_xx": 0.0,
                    "sigma_yy": 0.0,
                    "sigma_zz": 0.0,
                    "sigma_xy": sigma_xy,
                    "sigma_yz": 0.0,
                    "sigma_xz": 0.0,
                    "driver_version": 2.0,
                    "domain_time": float(step + 1),
                }
            )
        stale_rows = [dict(row) for row in opensees_rows]
        for row in stale_rows:
            row.pop("driver_version")
        try:
            compare_histories(
                fixture, hercules_rows, stale_rows, root / "stale-comparison.csv"
            )
        except ValueError as error:
            assert "prescribed-displacement driver version 2" in str(error)
        else:
            raise AssertionError("stale OpenSees CSV was not rejected")
        comparison = compare_histories(
            fixture, hercules_rows, opensees_rows, root / "comparison.csv"
        )
        assert comparison["stress_metrics"]["sig_xy"]["normalized_root_mean_square_difference"] == 0.0
        assert comparison["shear_work"]["relative_difference"] == 0.0
    print("PASS: OpenSees BAE comparison parser and driver-generation self-test")
    return 0


def main() -> int:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fixture",
        type=Path,
        default=script_dir / "vonmises_bae_opensees_synthetic.in",
        help="Hercules BAE material-point fixture",
    )
    opensees_source = parser.add_mutually_exclusive_group()
    opensees_source.add_argument(
        "--opensees", type=Path, help="OpenSees executable to run"
    )
    opensees_source.add_argument(
        "--opensees-csv",
        type=Path,
        help="existing OpenSees CSV produced by opensees_bae_synthetic.tcl",
    )
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--psolve", type=Path, help="Hercules psolve executable to run")
    source.add_argument("--hercules-csv", type=Path, help="existing Hercules CPU material-point CSV")
    parser.add_argument("--mpiexec", default="mpirun", help="MPI launcher used with --psolve")
    parser.add_argument(
        "--beta",
        type=float,
        default=1.0,
        help="OpenSees integration beta; 1 matches the Hercules end-state rule",
    )
    parser.add_argument("--k0", type=float, default=0.5, help="OpenSees K0 coefficient; inactive in material stage 2")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("opensees_bae_comparison_output"),
    )
    parser.add_argument(
        "--max-shear-nrmse",
        type=float,
        help="optional failure threshold for normalized sigma_xy RMS error",
    )
    parser.add_argument(
        "--max-work-relative-error",
        type=float,
        help="optional failure threshold for relative shear-work error",
    )
    parser.add_argument("--selftest", action="store_true")
    args = parser.parse_args()

    fixture_path = args.fixture.resolve()
    if args.selftest:
        return run_selftest(fixture_path)
    if args.opensees is None and args.opensees_csv is None:
        parser.error("provide either --opensees or --opensees-csv")
    if args.psolve is None and args.hercules_csv is None:
        parser.error("provide either --psolve or --hercules-csv")
    if not 0.0 <= args.beta <= 1.0:
        parser.error("--beta must be between 0 and 1")

    fixture = load_fixture(fixture_path)
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    if args.psolve is not None:
        run_command(
            [
                args.mpiexec,
                "-np",
                "1",
                str(args.psolve.resolve()),
                "--nonlinear-material-point",
                str(fixture_path),
            ],
            output_dir,
            output_dir / "hercules.stdout.log",
            output_dir / "hercules.stderr.log",
        )
        hercules_csv = output_dir / f"{fixture.output_prefix}.cpu.csv"
    else:
        hercules_csv = args.hercules_csv.resolve()

    if args.opensees_csv is not None:
        opensees_csv = args.opensees_csv.resolve()
    else:
        tcl_pathname = output_dir / "opensees_bae_driver.tcl"
        opensees_csv = output_dir / "opensees_bae.csv"
        write_opensees_driver(fixture, tcl_pathname, opensees_csv, args.beta, args.k0)
        run_command(
            [str(args.opensees.resolve()), str(tcl_pathname)],
            output_dir,
            output_dir / "opensees.stdout.log",
            output_dir / "opensees.stderr.log",
        )
    try:
        summary = compare_histories(
            fixture,
            read_csv(hercules_csv),
            read_csv(opensees_csv),
            output_dir / "comparison.csv",
        )
    except ValueError as error:
        print(f"Comparison input error: {error}", file=sys.stderr)
        return 2
    summary.update(
        {
            "fixture": str(fixture_path),
            "hercules_csv": str(hercules_csv),
            "opensees_csv": str(opensees_csv),
            "parameter_mapping": {
                "K": fixture.bulk,
                "G": fixture.shear,
                "Su": fixture.su,
                "Ho": 0.0,
                "h": fixture.psi0 * fixture.shear,
                "m": fixture.exponent,
                "beta": args.beta,
                "K0": args.k0,
            },
            "interpretation": (
                "Characterization comparison, not an equivalence proof. OpenSees and Hercules "
                "use different nonlinear integration and zero-increment history rules."
            ),
        }
    )
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    shear = summary["stress_metrics"]["sig_xy"]
    work = summary["shear_work"]
    print(f"Hercules/OpenSees steps: {fixture.steps}")
    print(f"Maximum |delta sigma_xy|: {shear['maximum_absolute_difference']:.6e} Pa")
    print(f"Normalized sigma_xy RMSE: {shear['normalized_root_mean_square_difference']:.6e}")
    print(f"Relative shear-work difference: {work['relative_difference']:.6e}")
    print(f"Detailed comparison: {output_dir / 'comparison.csv'}")
    print(f"Summary: {summary_path}")

    failed = False
    if (
        args.max_shear_nrmse is not None
        and shear["normalized_root_mean_square_difference"] > args.max_shear_nrmse
    ):
        print("FAIL: normalized sigma_xy RMSE exceeds configured threshold", file=sys.stderr)
        failed = True
    if (
        args.max_work_relative_error is not None
        and work["relative_difference"] > args.max_work_relative_error
    ):
        print("FAIL: relative shear-work error exceeds configured threshold", file=sys.stderr)
        failed = True
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
