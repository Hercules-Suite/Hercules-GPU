#!/usr/bin/env python3
"""Render Hercules horizontal plane displacements as vertical-velocity PNGs.

Select one binary file and supply its grid size and snapshot interval manually.
No parameters.in or other configuration file is read. Relative paths resolve
from the current working directory. --dt is the time between saved snapshots,
in seconds, not the solver time step.
Requires NumPy; rendering also requires PyGMT/GMT, xarray, and joblib.

    python animations_mpi.py --plane-file planes_fullDomain/planedisplacements.0 \\
        --nx 201 --ny 126 --dx 500 --dy 325 --dt 0.02 --check-only
    mpirun -np 2 python animations_mpi.py \\
        --plane-file planes_fullDomain/planedisplacements.0 \\
        --nx 201 --ny 126 --dx 500 --dy 325 --dt 0.02 --workers 4

Storage: [snapshot, along-strike point, down-dip point, Cartesian component].
This renderer assumes strike=0 (modulo 360), dip=0. For these planes,
rows advance in solver X and columns in solver Y. Nothing is transposed or
mirrored. Plot coordinates are distances from the plane origin, not lat/lon.
"""

import argparse
import math
import os
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class Plane:
    path: Path
    nx: int
    ny: int
    dx: float
    dy: float
    dt: float
    nt: int

    def open(self):
        # Native float64 matches the writer on the producing architecture.
        return np.memmap(self.path, dtype=np.float64, mode="r",
                         shape=(self.nt, self.nx, self.ny, 3), order="C")

    def coordinates(self):
        # xarray's x/y names are PLOT axes: x=solver Y, y=solver X.
        return np.arange(self.ny) * self.dy / 1000, np.arange(self.nx) * self.dx / 1000


def positive(value, name):
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def load_plane(plane_file, nx, ny, dt, dx=500.0, dy=325.0):
    """Validate the manually specified layout and map one plane's binary file."""
    if nx < 2 or ny < 2:
        raise ValueError("Rendering requires at least two points on each plane axis")
    positive(dt, "snapshot interval --dt")
    positive(dx, "along-strike spacing")
    positive(dy, "down-dip spacing")
    path = Path(plane_file).resolve()
    frame_bytes = nx * ny * 3 * np.dtype(np.float64).itemsize
    nt, remainder = divmod(path.stat().st_size, frame_bytes)
    if remainder:
        raise ValueError(f"{path}: size is not a multiple of {frame_bytes} bytes "
                         f"per snapshot (trailing bytes: {remainder})")
    if nt < 2:
        raise ValueError("At least two displacement snapshots are needed for velocity")
    return Plane(path, nx, ny, dx, dy, dt, nt)


def vertical_velocity(data, frame, dt):
    """Return unclipped, downward-positive cm/s in the original (X, Y) order."""
    positive(dt, "snapshot interval")
    nt = data.shape[0]
    if nt < 2 or not 0 <= frame < nt:
        raise ValueError("Velocity requires two snapshots and a valid frame index")
    # Central differences inside the record; one-sided differences at its ends.
    left, right = max(0, frame - 1), min(nt - 1, frame + 1)
    velocity = (data[right, :, :, 2] - data[left, :, :, 2]) * (
        100.0 / ((right - left) * dt))
    if not np.isfinite(velocity).all():
        raise ValueError(f"Non-finite velocity in frame {frame}; check input data")
    return velocity


def detect_rank_info(environ=None):
    environ = os.environ if environ is None else environ
    for rank_name, size_name in (
        ("OMPI_COMM_WORLD_RANK", "OMPI_COMM_WORLD_SIZE"),
        ("PMI_RANK", "PMI_SIZE"),
        ("PMIX_RANK", "PMIX_SIZE"),
        ("SLURM_PROCID", "SLURM_NTASKS"),
    ):
        if rank_name in environ and size_name in environ:
            rank, size = int(environ[rank_name]), int(environ[size_name])
            if size < 1 or not 0 <= rank < size:
                raise ValueError(f"Invalid rank/size: {rank}/{size}")
            return rank, size
    if any(key in environ for key in ("OMPI_COMM_WORLD_RANK", "PMI_RANK",
                                      "PMIX_RANK", "SLURM_PROCID")):
        raise ValueError("MPI rank is set but no matching world size is available")
    return 0, 1


def partition(total, rank, size):
    return total * rank // size, total * (rank + 1) // size


def render_chunk(plane, args, start, stop, rank):
    """One worker owns a mapping, GMT directory, and CPT for its whole chunk."""
    local_tmp = os.environ.get("SLURM_TMPDIR") or os.environ.get("TMPDIR") or "/tmp"
    with tempfile.TemporaryDirectory(prefix=f"hercules_gmt_rank{rank}_", dir=local_tmp) as scratch:
        # Set these before importing PyGMT in THIS worker process.
        os.environ["GMT_USERDIR"] = scratch
        os.environ["GMT_TMPDIR"] = scratch
        import pygmt
        import xarray as xr

        cpt = str(Path(scratch) / "vertical_velocity.cpt")
        # GMT B/F colors saturate the display without modifying the data.
        pygmt.makecpt(cmap="vik", series=[-args.velocity_limit, args.velocity_limit],
                      background=True, output=cpt)
        data = plane.open()
        horizontal, vertical = plane.coordinates()
        width, height = horizontal[-1], vertical[-1]
        region = [0, width, 0, height]
        projection = f"X{20 * width / height:g}c/20c"
        for frame in range(start, stop):
            velocity = vertical_velocity(data, frame, plane.dt)
            if args.progress_every and (frame == start or frame % args.progress_every == 0):
                saturated = np.count_nonzero(np.abs(velocity) > args.velocity_limit)
                print(f"Rank {rank}: frame {frame + 1}/{plane.nt}, "
                      f"vz=[{velocity.min():.6g}, {velocity.max():.6g}] cm/s, "
                      f"{100 * saturated / velocity.size:.2f}% beyond color limits", flush=True)
            grid = xr.DataArray(velocity, coords={"y": vertical, "x": horizontal},
                                dims=("y", "x"))
            fig = pygmt.Figure()
            with pygmt.config(FONT_LABEL="18p", FONT_ANNOT_PRIMARY="14p",
                              FONT_TITLE="20p", MAP_FRAME_TYPE="plain"):
                fig.grdimage(grid=grid, region=region, projection=projection, cmap=cpt,
                             interpolation="c+c", cores=1,
                             frame=["xaf+lY offset from plane origin (km)",
                                    "yaf+lX offset from plane origin (km)",
                                    f"+t{args.title}"])
                fig.colorbar(cmap=cpt,
                             frame=["a", "x+lVertical velocity (cm/s, positive down)"],
                             position="JMR+o1.5c/0c+w12c/0.8c")
                fig.text(text=f"Time = {args.first_time + frame * plane.dt:.3f} s",
                         x=0.025 * width, y=0.98 * height, font="16p,Helvetica,black",
                         justify="TL", fill="white@20", pen="0.5p,black", clearance="0.15c")
                target = args.output_dir / f"frame_{frame:05d}.png"
                temporary = target.with_name(f".{target.stem}.{os.getpid()}.png")
                try:
                    fig.savefig(str(temporary), dpi=args.dpi)
                    os.replace(temporary, target)
                finally:
                    temporary.unlink(missing_ok=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--plane-file", type=Path, required=True,
                        help="Binary displacement file for the selected plane")
    parser.add_argument("--nx", type=int, required=True,
                        help="Along-strike point count (solver X for strike=0, dip=0)")
    parser.add_argument("--ny", type=int, required=True,
                        help="Down-dip point count (solver Y for strike=0, dip=0)")
    parser.add_argument("--dt", type=float, required=True,
                        help="Time between SAVED snapshots, in seconds; not the solver time step")
    parser.add_argument("--dx", type=float, default=500.0,
                        help="Along-strike spacing in metres (default: 500)")
    parser.add_argument("--dy", type=float, default=325.0,
                        help="Down-dip spacing in metres (default: 325)")
    parser.add_argument("--output-dir", type=Path, default=Path("Animations_velocity"))
    parser.add_argument("--workers", type=int,
                        default=int(os.environ.get("FRAME_WORKERS_PER_RANK",
                                    os.environ.get("FRAME_WORKERS_PER_NODE", "1"))),
                        help="Worker processes PER MPI RANK (default: 1, or worker env setting)")
    parser.add_argument("--velocity-limit", type=float,
                        default=float(os.environ.get("VELOCITY_LIMIT_CM_S", "10")),
                        help="Symmetric color limits in cm/s (default: 10)")
    parser.add_argument("--progress-every", type=int,
                        default=int(os.environ.get("FRAME_PROGRESS_EVERY", "100")))
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--title", default="Hercules vertical velocity")
    parser.add_argument("--first-time", type=float, default=0.0,
                        help="Time of the first stored snapshot in seconds (default: 0); set for restarted output")
    parser.add_argument("--check-only", action="store_true",
                        help="Check layout and sample velocities without importing plotting libraries")
    args = parser.parse_args(argv)
    positive(args.workers, "workers")
    positive(args.velocity_limit, "velocity limit")
    positive(args.dpi, "DPI")
    if args.progress_every < 0:
        raise ValueError("progress-every must be nonnegative (0 disables progress)")
    rank, ranks = detect_rank_info()
    plane = load_plane(args.plane_file, args.nx, args.ny, args.dt, args.dx, args.dy)
    if not math.isfinite(args.first_time):
        raise ValueError("first-time must be finite")
    if rank == 0:
        horizontal, vertical = plane.coordinates()
        print(f"Plane: {plane.path}\n"
              "Manual geometry: assumes strike=0, dip=0\n"
              f"Layout: ({plane.nt}, {plane.nx}, {plane.ny}, 3), native float64, C order\n"
              f"Rows: solver X, spacing {plane.dx:g} m; columns: solver Y, spacing {plane.dy:g} m\n"
              f"Plot: Y rightward 0..{horizontal[-1]:g} km; X upward 0..{vertical[-1]:g} km\n"
              f"No transpose or mirror; coordinates relative to plane origin\n"
              f"Snapshot interval: {plane.dt:g} s; first/last time: "
              f"{args.first_time:g}/{args.first_time + (plane.nt - 1) * plane.dt:g} s\n"
              f"Velocity: cm/s, positive down; color limits +/-{args.velocity_limit:g}\n"
              "Endpoint velocities use one-sided differences; interiors use central differences.", flush=True)
    if args.check_only:
        if rank == 0:
            data = plane.open()
            for frame in sorted({0, plane.nt // 2, plane.nt - 1}):
                velocity = vertical_velocity(data, frame, plane.dt)
                print(f"Frame {frame}: vz=[{velocity.min():.6g}, {velocity.max():.6g}] cm/s")
            print("Layout checks passed; sampled velocities are finite. "
                  "This does not verify the source location or all snapshots.")
        return

    start, stop = partition(plane.nt, rank, ranks)
    jobs = min(args.workers, stop - start)
    print(f"Rank {rank}/{ranks}: frames [{start}, {stop}), {jobs} workers", flush=True)
    if jobs == 0:
        return
    args.output_dir = args.output_dir.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    # Use a separate directory for each run to avoid mixing old frame sequences.
    if jobs == 1:
        render_chunk(plane, args, start, stop, rank)
    else:
        from joblib import Parallel, delayed, parallel_backend
        # Each task is a whole chunk, so process-local state is initialized once
        # without relying on cloudpickle retaining mutable __main__ globals.
        with parallel_backend("loky", inner_max_num_threads=1):
            Parallel(n_jobs=jobs)(
                delayed(render_chunk)(plane, args, start + left, start + right, rank)
                for left, right in (partition(stop - start, job, jobs) for job in range(jobs))
            )
    print(f"Rank {rank}: local frames completed", flush=True)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, ImportError) as error:
        sys.exit(f"Animation error: {error}")
