#!/usr/bin/env python3
"""Locate plane displacement peaks and export compact time histories (NumPy).

Scan completed headerless float64 XYZ files twice, one snapshot at a time.
Optional geometry describes a horizontal, unrotated plane in solver X/Y metres.
Distance-to-boundary flags are geometric candidates, not element classification.
Exact classification requires plane_point_metadata.RANK.csv from the same mesh
and plane configuration, exported with HERCULES_PLANE_METADATA=1.
"""

import argparse
import csv
import json
import os
from pathlib import Path
import sys

import numpy as np

from compare_plane_outputs import file_signature, snapshot_count

AXES = "XYZ"


def metadata(paths, plane, nx, ny):
    """Keep unknown ownership explicit; reject duplicate and out-of-grid rows."""
    records = {}
    for path in paths:
        with path.open(newline="") as stream:
            reader = csv.DictReader(stream)
            required = {"plane", "i", "j", "medium", "cpml", "x_m", "y_m", "z_m",
                        "valid", "acoustic_reconstruction", "pressure_defined"}
            if not required.issubset(reader.fieldnames or []):
                raise ValueError(f"{path}: missing plane metadata columns")
            for row in reader:
                if int(row["plane"]) != plane:
                    continue
                i, j = int(row["i"]), int(row["j"])
                if not (0 <= i < nx and 0 <= j < ny):
                    raise ValueError(f"{path}: point ({i},{j}) outside requested grid")
                k = i * ny + j
                if k in records:
                    raise ValueError(f"Duplicate metadata for point ({i},{j})")
                if row["medium"] not in ("acoustic", "solid", "inactive", "unknown"):
                    raise ValueError(f"{path}: unknown medium {row['medium']!r}")
                if row["cpml"] not in ("-1", "0", "1"):
                    raise ValueError(f"{path}: invalid CPML flag")
                for key in ("x_m", "y_m", "z_m"):
                    if not np.isfinite(float(row[key])):
                        raise ValueError(f"{path}: nonfinite coordinate")
                records[k] = row
    if paths and not records:
        raise ValueError(f"Metadata contains no rows for plane {plane}")
    return records


def point_groups(args, records):
    n = args.nx * args.ny
    groups = {"all": np.ones(n, dtype=bool)}
    xy = None
    boundary_distance = None
    if args.spacing is not None:
        i, j = np.divmod(np.arange(n), args.ny)
        xy = np.column_stack((args.xy_origin[0] + i * args.spacing[0],
                              args.xy_origin[1] + j * args.spacing[1]))
        # Use the actual mesh query coordinates when metadata is available.
        for k, row in records.items():
            actual = (float(row["x_m"]), float(row["y_m"]))
            if not np.allclose(xy[k], actual, rtol=0, atol=1e-4):
                raise ValueError("Supplied horizontal-plane geometry disagrees with metadata")
            xy[k] = actual
        boundary_distance = np.minimum.reduce((xy[:, 0], xy[:, 1],
                                               args.domain[0] - xy[:, 0],
                                               args.domain[1] - xy[:, 1]))
        if np.any(boundary_distance < -1e-4):
            raise ValueError("Plane points extend beyond the supplied domain")
        groups["lateral_boundary_band"] = boundary_distance < args.cpml_thickness
        groups["geometric_interior"] = ~groups["lateral_boundary_band"]
    if records:
        medium = np.full(n, "unknown", dtype="U16")
        cpml = np.full(n, -1, dtype=np.int8)
        for k, row in records.items():
            if row["valid"] == "1":
                medium[k], cpml[k] = row["medium"], int(row["cpml"])
        for name in ("acoustic", "solid", "inactive", "unknown"):
            groups[name] = medium == name
        groups["cpml"] = cpml == 1
        groups["non_cpml"] = cpml == 0
        for name in ("acoustic", "solid"):
            groups[name + "_cpml"] = (medium == name) & (cpml == 1)
            groups[name + "_non_cpml"] = (medium == name) & (cpml == 0)
    return {name: mask for name, mask in groups.items() if mask.any()}, xy, boundary_distance


def frames(stream, nt, n, dtype):
    for step in range(nt):
        data = np.fromfile(stream, dtype=dtype, count=n * 3)
        if data.size != n * 3:
            raise ValueError(f"Unexpected end of file at snapshot {step}")
        data = data.reshape(n, 3)
        bad = np.argwhere(~np.isfinite(data))
        if bad.size:
            point, component = bad[0]
            raise ValueError(f"Nonfinite displacement at snapshot {step}, flat point "
                             f"{point}, component {AXES[component]}")
        yield step, data


def diagnose(args):
    path = args.file
    if args.output_dir.exists():
        raise FileExistsError(f"Use a new diagnostic directory: {args.output_dir}")
    signature = file_signature(path.stat())
    nt = snapshot_count(path, args.nx, args.ny)
    if args.expected_snapshots is not None and nt != args.expected_snapshots:
        raise ValueError(f"Expected {args.expected_snapshots} snapshots; found {nt}")
    records = metadata(args.metadata, args.plane_id, args.nx, args.ny)
    groups, xy, distance = point_groups(args, records)
    n = args.nx * args.ny
    peak = np.full((n, 3), -1.0)
    signed = np.zeros((n, 3))
    when = np.zeros((n, 3), dtype=np.int64)
    envelope = np.zeros((nt, len(groups) * 3))
    dtype = np.dtype({"native": "=f8", "little": "<f8", "big": ">f8"}[args.endian])
    with path.open("rb") as stream:
        if file_signature(os.fstat(stream.fileno())) != signature:
            raise ValueError("File changed before scanning")
        for step, data in frames(stream, nt, n, dtype):
            absolute = np.abs(data)
            update = absolute > peak
            peak[update], signed[update], when[update] = absolute[update], data[update], step
            for g, mask in enumerate(groups.values()):
                envelope[step, 3*g:3*g+3] = absolute[mask].max(axis=0)
        if file_signature(os.fstat(stream.fileno())) != signature:
            raise ValueError("File changed while scanning; use completed output")

    selected = {i * args.ny + j for i, j in args.points}
    peak_rows, group_stats = [], {}
    for name, mask in groups.items():
        ids = np.flatnonzero(mask)
        group_stats[name] = {"points": int(mask.sum()),
                             "peak_abs_m": peak[mask].max(axis=0).tolist(),
                             "points_exceeding_threshold_by_component":
                                 (peak[mask] > args.threshold).sum(axis=0).tolist()}
        for c, component in enumerate(AXES):
            # Global top-N distinct points, plus each group's strongest point.
            count = min(args.top if name == "all" else 1, ids.size)
            ranked = ids[np.argsort(-peak[ids, c], kind="stable")[:count]]
            for rank, k in enumerate(ranked, 1):
                selected.add(int(k))
                i, j = divmod(int(k), args.ny)
                row = records.get(int(k), {})
                coords = [row.get(key, "") for key in ("x_m", "y_m", "z_m")]
                if not row and xy is not None:
                    coords[:2] = xy[k].tolist()
                peak_rows.append([name, component, rank, i, j, int(when[k, c]),
                                  args.first_time + when[k, c] * args.dt,
                                  signed[k, c], peak[k, c], *coords,
                                  "" if distance is None else distance[k],
                                  row.get("medium", "unknown"), row.get("cpml", "-1"),
                                  row.get("acoustic_reconstruction", ""),
                                  row.get("pressure_defined", ""),
                                  *[row.get(key, "") for key in
                                    ("vp_m_s", "vs_m_s", "rho_kg_m3", "h_m",
                                     "rank", "element", "topostrip")]])

    # Refuse to overwrite a previous diagnostic. Summary is written only after
    # both scans succeed; a failed second pass can leave partial CSV files.
    args.output_dir.mkdir(parents=True, exist_ok=False)
    with (args.output_dir / "peak_points.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["group", "component", "rank", "i", "j", "snapshot", "time_s",
                         "value_m", "abs_value_m", "x_m", "y_m", "z_m",
                         "lateral_boundary_distance_m", "medium", "cpml",
                         "acoustic_reconstruction", "pressure_defined", "vp_m_s", "vs_m_s",
                         "rho_kg_m3", "h_m", "owner_rank", "owner_element", "topostrip"])
        writer.writerows(peak_rows)
    with (args.output_dir / "frame_peaks.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["snapshot", "time_s"] +
                        [f"{name}_{c}_abs_m" for name in groups for c in AXES])
        for step in range(nt):
            writer.writerow([step, args.first_time + step * args.dt, *envelope[step]])
    selected = sorted(selected)
    with path.open("rb") as source, (args.output_dir / "histories.csv").open("w", newline="") as out:
        if file_signature(os.fstat(source.fileno())) != signature:
            raise ValueError("File changed before extracting histories")
        writer = csv.writer(out)
        writer.writerow(["snapshot", "time_s", "i", "j", "ux_m", "uy_m", "uz_m"])
        for step, data in frames(source, nt, n, dtype):
            t = args.first_time + step * args.dt
            writer.writerows([step, t, k // args.ny, k % args.ny, *data[k]] for k in selected)
        if (file_signature(os.fstat(source.fileno())) != signature or
                file_signature(path.stat()) != signature):
            raise ValueError("File changed while extracting histories")
    summary = {"file": str(path.resolve()), "nx": args.nx, "ny": args.ny,
               "snapshots": nt, "saved_dt_s": args.dt, "first_time_s": args.first_time,
               "last_time_s": args.first_time + (nt - 1) * args.dt,
               "threshold_m": args.threshold, "metadata_points": len(records),
               "history_points": len(selected), "groups": group_stats,
               "notes": ["Indices are zero-based; components are solver north/east/down.",
                         "Boundary-band geometry is not actual element CPML membership.",
                         "Unknown medium/CPML requires metadata from the same mesh and planes.",
                         "Threshold exceedance is a diagnostic, not a physical acceptance test."]}
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Scanned {nt} snapshots, t={args.first_time:g}..{summary['last_time_s']:g} s")
    for row in peak_rows:
        if row[0] == "all" and row[2] == 1:
            print(f"{row[1]} peak: {row[7]:.9g} m at ({row[3]},{row[4]}), "
                  f"snapshot={row[5]}, t={row[6]:g} s; medium={row[13]}, cpml={row[14]}")
    for name, stats in group_stats.items():
        print(f"{name}: points={stats['points']}, XYZ peaks(m)={stats['peak_abs_m']}")
    print(f"Wrote {len(selected)} complete point histories and reports to {args.output_dir}")
    return summary


def point(value):
    try:
        i, j = map(int, value.split(","))
        return i, j
    except ValueError as exc:
        raise argparse.ArgumentTypeError("point must be I,J (zero-based)") from exc


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", type=Path)
    parser.add_argument("--nx", type=int, required=True)
    parser.add_argument("--ny", type=int, required=True)
    parser.add_argument("--dt", type=float, required=True, help="SAVED snapshot interval in seconds")
    parser.add_argument("--first-time", type=float, default=0.0)
    parser.add_argument("--expected-snapshots", type=int)
    parser.add_argument("--top", type=int, default=5, help="Top distinct points per component")
    parser.add_argument("--threshold", type=float, default=1.0, help="Report point counts exceeding this many metres")
    parser.add_argument("--points", nargs="*", type=point, default=[], help="Additional I,J histories")
    parser.add_argument("--spacing", nargs=2, type=float, metavar=("DX", "DY"))
    parser.add_argument("--xy-origin", nargs=2, type=float, default=(0.0, 0.0))
    parser.add_argument("--domain", nargs=2, type=float, metavar=("LX", "LY"))
    parser.add_argument("--cpml-thickness", type=float)
    parser.add_argument("--metadata", nargs="+", type=Path, default=[])
    parser.add_argument("--plane-id", type=int, default=0, help="Plane ID to select from metadata")
    parser.add_argument("--endian", choices=("native", "little", "big"), default="native")
    parser.add_argument("--output-dir", type=Path, required=True, help="New directory; must not exist")
    args = parser.parse_args(argv)
    if args.nx < 1 or args.ny < 1 or args.top < 1 or args.plane_id < 0:
        parser.error("grid counts and top must be positive; plane-id must be nonnegative")
    if not np.isfinite([args.dt, args.first_time, args.threshold, *args.xy_origin]).all() or args.dt <= 0 or args.threshold < 0:
        parser.error("times, origin, and threshold must be finite; dt > 0 and threshold >= 0")
    if args.expected_snapshots is not None and args.expected_snapshots < 1:
        parser.error("expected-snapshots must be positive")
    geometry = (args.spacing is not None, args.domain is not None, args.cpml_thickness is not None)
    if any(geometry):
        if not all(geometry):
            parser.error("spacing, domain, and cpml-thickness must be supplied together")
        values = [*args.spacing, *args.domain, args.cpml_thickness]
        if not np.isfinite(values).all() or min(values) <= 0:
            parser.error("spacing, domain, and cpml-thickness must be positive and finite")
        if 2 * args.cpml_thickness >= min(args.domain):
            parser.error("lateral CPML bands must leave an interior")
    if any(not (0 <= i < args.nx and 0 <= j < args.ny) for i, j in args.points):
        parser.error("requested point outside grid")
    try:
        diagnose(args)
    except (OSError, ValueError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
