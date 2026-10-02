#!/usr/bin/env python3
"""Compare completed acoustic station pressure and motion histories (NumPy).

Read columns by their header, verify ownership/coordinates, and retain raw
signals. Metrics describe sensitivity and integration consistency; they do
not establish physical accuracy. Optional plots require Matplotlib.
"""
import argparse
import csv
import json
from pathlib import Path
import re

import numpy as np


FIELDS = {
    "u": ("X|(m)", "Y-(m)", "Z.(m)"),
    "v": ("X|(m/s)", "Y-(m/s)", "Z.(m/s)"),
    "a": ("X|(m/s2)", "Y-(m/s2)", "Z.(m/s2)"),
    "p": ("Pressure(Pa)",),
}


def signature(path):
    stat = path.stat()
    return stat.st_size, stat.st_mtime_ns


def read_station(path, expected_samples=None):
    before = signature(path)
    with path.open() as stream:
        header = next((line[1:].split() for line in stream
                       if line.startswith("#") and "Time(s)" in line), None)
    if not header or len(header) != len(set(header)):
        raise ValueError(f"{path}: missing or ambiguous column header")
    data = np.loadtxt(path, comments="#", ndmin=2)
    if signature(path) != before:
        raise ValueError(f"{path}: file changed while reading; use completed runs")
    if data.shape[1] != len(header) or data.shape[0] < 3:
        raise ValueError(f"{path}: incomplete history or column count mismatch")
    if expected_samples is not None and len(data) != expected_samples:
        raise ValueError(f"{path}: expected {expected_samples} samples, got {len(data)}")
    try:
        time = data[:, header.index("Time(s)")]
        fields = {k: data[:, [header.index(c) for c in columns]] for k, columns in FIELDS.items()}
    except ValueError as error:
        raise ValueError(f"{path}: enable station velocity, acceleration and pressure output") from error
    if not np.isfinite(time).all() or np.any(np.diff(time) <= 0):
        raise ValueError(f"{path}: invalid/non-increasing time axis")
    if not all(np.isfinite(values).all() for values in fields.values()):
        raise ValueError(f"{path}: nonfinite acoustic motion/pressure")
    return time, fields


def read_mediums(path):
    with path.open() as stream:
        header = stream.readline().lstrip("# ").split()
        rows = [dict(zip(header, line.split())) for line in stream
                if line.strip() and not line.startswith("#")]
    result = {int(row["station_id"]): row for row in rows}
    if len(result) != len(rows):
        raise ValueError(f"{path}: duplicate station metadata")
    return result


def integrate(time, values):
    increments = .5 * (values[:-1] + values[1:]) * np.diff(time)[:, None]
    return np.vstack((np.zeros(values.shape[1]), np.cumsum(increments, axis=0)))


def metrics(time, fields, late_start, rho, gravity_surface=False):
    late = time >= late_start
    if np.count_nonzero(late) < 3:
        raise ValueError("Late window contains fewer than three samples")
    result = {}
    for key, values in fields.items():
        result[key + "_peak_abs"] = np.max(np.abs(values), axis=0).tolist()
        result[key + "_final"] = values[-1].tolist()
        result[key + "_late_mean"] = np.mean(values[late], axis=0).tolist()
        result[key + "_late_rms"] = np.sqrt(np.mean(values[late]**2, axis=0)).tolist()
    result["u_late_slope_m_s"] = np.polyfit(time[late] - time[late][0], fields["u"][late], 1)[0].tolist()
    # These use SAVED samples; coarsening or text rounding can create residuals.
    result["u_integral_v_max_residual_m"] = np.max(np.abs(
        fields["u"] - fields["u"][0] - integrate(time, fields["v"])), axis=0).tolist()
    result["v_integral_a_max_residual_m_s"] = np.max(np.abs(
        fields["v"] - fields["v"][0] - integrate(time, fields["a"])), axis=0).tolist()
    if gravity_surface:
        residual = fields["u"][:, 2] + fields["p"][:, 0] / (rho * 9.81)
        result["surface_relation_rms_m"] = float(np.sqrt(np.mean(residual**2)))
        result["surface_relation_late_rms_m"] = float(np.sqrt(np.mean(residual[late]**2)))
    return result


def check_probe(probe, medium):
    if medium["owner_medium"] != "acoustic" or medium["pressure_defined"] != "yes":
        raise ValueError(f"station {probe['station']}: expected an acoustic owner with pressure")
    for key in ("x_m", "y_m", "z_m"):
        if key in probe and abs(float(medium[key]) - float(probe[key])) > .002:
            raise ValueError(f"station {probe['station']}: {key} disagrees with probe manifest")
    rho = float(medium["rho_kgm3"])
    if not np.isfinite(rho) or rho <= 0:
        raise ValueError(f"station {probe['station']}: invalid density")
    return rho


def plot_probe(output, probe, curves, gravity_cases):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(4, 3, figsize=(13, 10), sharex=True)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for case_index, (name, (time, fields, rho)) in enumerate(curves.items()):
        color = colors[case_index % len(colors)]
        for row, (key, unit) in enumerate((("u", "m"), ("v", "m/s"), ("a", "m/s²"))):
            for component in range(3):
                axes[row, component].plot(time, fields[key][:, component], label=name, color=color, lw=1)
                axes[row, component].set_ylabel(f"{key}{'XYZ'[component]} ({unit})")
        axes[3, 0].plot(time, fields["p"][:, 0], label=name, color=color, lw=1)
        axes[3, 1].plot(time, fields["u"][:, 2], label=f"{name}: uZ", color=color, lw=1)
        if name in gravity_cases and probe.get("surface") and probe.get("cpml") == 0:
            pressure_depth = -fields["p"][:, 0] / (rho * 9.81)
            axes[3, 1].plot(time, pressure_depth, label=f"{name}: −p/(ρg)", color=color, ls="--", lw=1)
            axes[3, 2].plot(time, fields["u"][:, 2] - pressure_depth, color=color, lw=1)
    axes[3, 0].set_ylabel("Pressure (Pa)")
    axes[3, 1].set_ylabel("Surface relation (m)")
    axes[3, 2].set_ylabel("uZ + p/(ρg) (m)")
    for axis in axes.flat:
        axis.grid(alpha=.2)
    for axis in axes[-1]:
        axis.set_xlabel("Time (s)")
    axes[0, 0].legend(fontsize=8)
    axes[3, 1].legend(fontsize=7)
    figure.suptitle(f"Station {probe['station']} · {probe['label']} · CPML={probe.get('cpml', 'unknown')}\n"
                   "Raw histories; surface relation is a diagnostic, not an acceptance test")
    figure.tight_layout()
    figure.savefig(output / f"station.{probe['station']}.png", dpi=140)
    plt.close(figure)


def run(args):
    probes = json.loads(args.probes.read_text())
    ids = [int(p["station"]) for p in probes]
    if not probes or len(ids) != len(set(ids)) or any(i < 0 for i in ids):
        raise ValueError("Manifest must contain unique nonnegative station IDs")
    cases = {}
    for item in args.cases:
        name, sep, path = item.partition("=")
        if not sep or not re.fullmatch(r"[A-Za-z0-9_-]+", name) or name in cases or not path:
            raise ValueError("Cases must be unique NAME=CASE_DIRECTORY entries")
        cases[name] = Path(path)
    if set(args.gravity_cases) - set(cases):
        raise ValueError("Every gravity case must be present in --cases")
    if not np.isfinite(args.late_start):
        raise ValueError("late-start must be finite")
    if args.history_stride is not None and args.history_stride < 1:
        raise ValueError("history-stride must be positive")
    args.output_dir.mkdir(parents=True, exist_ok=False)
    metadata = {name: read_mediums(path / "output/stations/station_mediums.txt")
                for name, path in cases.items()}
    report = {"cases": {name: str(path.resolve()) for name, path in cases.items()},
              "late_start_s": args.late_start, "stations": [],
              "interpretation": "Raw signals, no detrending/filtering. Integral residuals use saved samples. "
              "Surface relation applies only to the declared flat gravity surface outside CPML; "
              "mass lumping, gradient recovery and time staggering limit pointwise agreement. "
              "Shared station/plane reconstruction is not an independent physics validation."}
    csv_rows = []
    for probe in probes:
        station = int(probe["station"])
        curves = {}
        record = {"probe": probe, "cases": {}, "differences_from_first_case": {}}
        baseline = None
        for name, path in cases.items():
            if station not in metadata[name]:
                raise ValueError(f"{name}: missing station {station} ownership")
            rho = check_probe(probe, metadata[name][station])
            time, fields = read_station(path / f"output/stations/station.{station}", args.expected_samples)
            if baseline is None:
                baseline = time, fields, rho
            else:
                if len(time) != len(baseline[0]) or not np.allclose(time, baseline[0], rtol=0, atol=1e-8):
                    raise ValueError(f"{name}: station time axes differ")
                if not np.isclose(rho, baseline[2], rtol=1e-8):
                    raise ValueError(f"{name}: probe density changed")
                record["differences_from_first_case"][name] = {
                    key + "_max_abs": np.max(np.abs(values - baseline[1][key]), axis=0).tolist()
                    for key, values in fields.items()}
            values = metrics(time, fields, args.late_start, rho,
                             name in args.gravity_cases and probe.get("surface") and probe.get("cpml") == 0)
            values.update(samples=len(time), first_time_s=float(time[0]), last_time_s=float(time[-1]),
                          rho_kgm3=rho, owner=metadata[name][station])
            record["cases"][name] = values
            curves[name] = time, fields, rho
            for axis, label in enumerate("XYZ"):
                csv_rows.append([name, station, probe["label"], probe.get("cpml", "unknown"), label,
                                 values["u_peak_abs"][axis], values["u_final"][axis],
                                 values["u_late_slope_m_s"][axis], values["v_late_mean"][axis],
                                 values["v_peak_abs"][axis], values["a_peak_abs"][axis],
                                 values["p_peak_abs"][0], values["p_final"][0]])
        if args.plot:
            plot_probe(args.output_dir, probe, curves, args.gravity_cases)
        if args.history_stride:
            archive = {}
            for name, (time, fields, rho) in curves.items():
                selected = np.unique(np.r_[np.arange(0, len(time), args.history_stride), len(time)-1])
                archive[name + "_time"] = time[selected]
                archive[name + "_rho"] = rho
                archive.update({name + "_" + key: values[selected] for key, values in fields.items()})
            np.savez_compressed(args.output_dir / f"station.{station}.npz", **archive)
        report["stations"].append(record)
    with (args.output_dir / "metrics.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["case", "station", "label", "cpml", "axis", "peak_u_m", "final_u_m",
                         "late_u_slope_m_s", "late_mean_v_m_s", "peak_v_m_s", "peak_a_m_s2",
                         "peak_pressure_Pa", "final_pressure_Pa"])
        writer.writerows(csv_rows)
    (args.output_dir / "summary.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(f"Compared {len(cases)} cases at {len(probes)} acoustic stations; wrote {args.output_dir}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", nargs="+", required=True, metavar="NAME=CASE_DIRECTORY")
    parser.add_argument("--probes", type=Path, required=True, help="JSON list: station, label, optional XYZ, surface, cpml")
    parser.add_argument("--expected-samples", type=int)
    parser.add_argument("--late-start", type=float, default=80., help="Start of late-time window in seconds (default 80)")
    parser.add_argument("--gravity-cases", nargs="*", default=[], help="Cases with gravity enabled for surface-relation diagnostics")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--plot", action="store_true", help="Write raw-history plots; requires Matplotlib")
    parser.add_argument("--history-stride", type=int, help="Export every Nth raw sample, plus the last, to compact NPZ files; metrics still use all samples")
    args = parser.parse_args()
    try:
        run(args)
    except (OSError, ValueError, KeyError, ImportError) as error:
        parser.exit(1, f"Error: {error}\n")


if __name__ == "__main__":
    main()
