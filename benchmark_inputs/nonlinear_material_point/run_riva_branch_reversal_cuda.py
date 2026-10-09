"""Generate or execute smooth-path RIVASand02/Custom02 CUDA parity fixtures.

Without --psolve this only writes inputs and reports that CUDA was not run.
The psolve material-point runner compares all 138 state values internally.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import shlex
import subprocess


def fixture(steps, custom, prefix, reversal_type=2):
    material = "RIVASandCustom02" if custom else "RIVASand02"
    row = "0.90 900"
    if custom:
        row += " 1.25 1.125 122.44207260468994 0.945 0.025 0.78 0.51 10 1.5 0.65"
    history = [[0.0]*6]
    for i in range(1, 2*steps+1):
        t = i/steps
        history.append([0, 0, 0.0001*t, 0, 0, 0.0015*math.sin(2*math.pi*t)])
    history.append(history[-1].copy())  # Exact zero-increment state check.
    text = f"""# Loading-branch research; two smooth shear cycles with vertical strain.
material_model = {material}
rivasand02_reversaltype = {reversal_type}
material_point_delta_time_sec = {1/steps:.17g}
material_point_steps = {len(history)}
no_substeps = 4
error_tolerance = 1.0e-1
backbone_errTol = 1.0e-5
tension_cutoff = no
material_point_compare_rtol = 1.0e-8
material_point_compare_atol = 1.0e-10
material_point_output_prefix = {prefix}
material_point_constants =
  {row}
material_point_initial_strain =
  0 0 0 0 0 0
material_point_initial_stress =
  -100000 -100000 -100000 0 0 0
material_point_initial_alpha =
  0 0 0 0 0 0
material_point_initial_alpha_in =
  0 0 0 0 0 0
material_point_initial_fabric =
  0 0 0 0 0 0
material_point_strain_history =
"""
    return text + "".join("  " + " ".join(f"{v:.17g}" for v in row) + "\n" for row in history)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--psolve", type=Path)
    parser.add_argument("--launcher", default="", help="Optional launcher, e.g. 'mpirun -np 1'")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reversal-type", type=int, choices=(1, 2, 3), default=2)
    args = parser.parse_args()
    checkpoint = {1: 8, 2: 7, 3: 9}[args.reversal_type]
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    results = []
    for steps in (128, 512, 2048):
        prefixes = []
        for custom in (False, True):
            name = f"branch_{'custom02' if custom else '02'}_n{steps}"
            prefixes.append(name)
            path = output / (name + ".in")
            path.write_text(fixture(steps, custom, name, args.reversal_type))
            if args.psolve:
                command = shlex.split(args.launcher) + [str(args.psolve.resolve()),
                           "--nonlinear-material-point", str(path)]
                result = subprocess.run(command, cwd=output, stdout=subprocess.PIPE,
                                        stderr=subprocess.STDOUT, text=True, timeout=600)
                (output / (name + ".log")).write_text(result.stdout)
                for marker in ("PASS: CPU/GPU material-point histories agree",
                               f"Validated CPU/GPU registered reversals: {4 if args.reversal_type in (2,3) or steps==128 else 0}",
                               f"kernel 4, checkpoint {checkpoint}",
                               f"rivasand02_reversaltype={args.reversal_type}"):
                    if result.returncode != 0 or marker not in result.stdout:
                        raise RuntimeError(f"{name} failed or used an old executable:\n{result.stdout[-4000:]}")
                for backend in ("cpu", "gpu"):
                    with (output / f"{name}.{backend}.csv").open() as stream:
                        rows = list(csv.DictReader(stream))
                    assert len(rows) == 2*steps+2
                    assert all(math.isfinite(float(v)) for row in rows for v in row.values())
                    assert all(int(row["status"]) == 0 for row in rows)
            results.append({"case": name, "status": "passed" if args.psolve else "generated; not run"})
        if args.psolve:
            for backend in ("cpu", "gpu"):
                assert (output / f"{prefixes[0]}.{backend}.csv").read_bytes() == (
                    output / f"{prefixes[1]}.{backend}.csv").read_bytes(), "Standard02/Custom02 mismatch"
    (output / "summary.json").write_text(json.dumps(results, indent=2) + "\n")
    print("PASS: CUDA branch-reference parity and event counts" if args.psolve else
          f"Generated six type-{args.reversal_type} CUDA fixtures; CUDA execution remains pending")


if __name__ == "__main__":
    main()
