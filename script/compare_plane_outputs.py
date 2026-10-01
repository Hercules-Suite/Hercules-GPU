#!/usr/bin/env python3
r"""Compare completed Hercules planedisplacements.N files using bounded memory.

Requires Python 3 and NumPy. Files contain headerless native-endian float64 in
[snapshot, along-strike point, down-dip point, XYZ component] order, in metres.
Grid shapes, sample times, plane geometry, and byte order are not stored in the
files: compare runs with matching inputs and supply the correct layout.

Example for seven planes (plane IDs are 0 through 6):
    python3 compare_plane_outputs.py planes planes_mpiio \
        --shapes 361x361 71x71 71x71 71x71 71x71 71x71 71x71 \
        --expected-snapshots 1000

The default is exact numerical equality. --atol (metres) and --rtol optionally
apply abs(candidate-reference) <= atol + rtol*abs(reference) to EVERY value.
Exit status: 0 = all finite values pass; 1 = differences/non-finite values;
2 = invalid input, incomplete files, changed files, or inconsistent layouts.
Byte identity is reported separately; signed zeros compare numerically equal.
"""

import argparse
from dataclasses import dataclass
import math
import os
from pathlib import Path
import re
import sys

import numpy as np


@dataclass
class ComponentStats:
    reference_peak: float = 0.0
    candidate_peak: float = 0.0
    max_abs_error: float = 0.0
    finite_differences: int = 0
    outside_tolerance: int = 0
    nonfinite_reference: int = 0
    nonfinite_candidate: int = 0
    max_location: object = None
    first_failure: object = None

    @property
    def passed(self):
        return not (self.outside_tolerance or self.nonfinite_reference
                    or self.nonfinite_candidate)


def location(flat_point, nx, ny, reference, candidate):
    frame, point = divmod(int(flat_point), nx * ny)
    i, j = divmod(point, ny)
    return (frame, i, j, float(reference), float(candidate))


def file_signature(stat):
    return stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns


def snapshot_count(path, nx, ny):
    if not path.is_file():
        raise ValueError(f"Not a regular file: {path}")
    size = path.stat().st_size
    frame_bytes = nx * ny * 3 * 8
    count, remainder = divmod(size, frame_bytes)
    if not count or remainder:
        raise ValueError(f"{path}: {size} bytes is not a positive whole number "
                         f"of {frame_bytes}-byte snapshots (remainder {remainder})")
    return count


def compare_plane(reference, candidate, nx, ny, snapshots, *,
                  chunk_points=262144, atol=0.0, rtol=0.0, endian="native"):
    dtype = np.dtype({"native": "=f8", "little": "<f8", "big": ">f8"}[endian])
    stats = [ComponentStats() for _ in range(3)]
    bitwise_identical = True
    total_points = snapshots * nx * ny
    expected_size = total_points * 3 * 8
    with reference.open("rb") as ref_file, candidate.open("rb") as test_file:
        signatures = [file_signature(os.fstat(f.fileno()))
                      for f in (ref_file, test_file)]
        if any(s[0] != expected_size for s in signatures):
            raise ValueError("File size changed after layout validation")
        for offset in range(0, total_points, chunk_points):
            count = min(chunk_points, total_points - offset)
            ref = np.fromfile(ref_file, dtype=dtype, count=count * 3)
            test = np.fromfile(test_file, dtype=dtype, count=count * 3)
            if ref.size != count * 3 or test.size != count * 3:
                raise ValueError("Unexpected end of file while reading "
                                 f"{reference.name}; use completed outputs")
            if not np.array_equal(ref.view(np.uint8), test.view(np.uint8)):
                bitwise_identical = False
            ref, test = ref.reshape(-1, 3), test.reshape(-1, 3)
            for c, stat in enumerate(stats):
                a, b = ref[:, c], test[:, c]
                finite_a, finite_b = np.isfinite(a), np.isfinite(b)
                valid = finite_a & finite_b
                stat.nonfinite_reference += int(np.count_nonzero(~finite_a))
                stat.nonfinite_candidate += int(np.count_nonzero(~finite_b))
                stat.reference_peak = max(stat.reference_peak, float(
                    np.max(np.abs(a[finite_a]), initial=0.0)))
                stat.candidate_peak = max(stat.candidate_peak, float(
                    np.max(np.abs(b[finite_b]), initial=0.0)))
                stat.finite_differences += int(np.count_nonzero(valid & (a != b)))
                with np.errstate(over="ignore", invalid="ignore"):
                    error = np.abs(b - a)
                    tolerance = atol + rtol * np.abs(a)
                outside = valid & ((error > tolerance) | ~np.isfinite(error))
                stat.outside_tolerance += int(np.count_nonzero(outside))
                failed = ~valid | outside
                if stat.first_failure is None and np.any(failed):
                    k = int(np.flatnonzero(failed)[0])
                    stat.first_failure = location(offset + k, nx, ny, a[k], b[k])
                # Invalid input values are counted above and always fail. They
                # must not hide the largest difference among the finite pairs.
                error[~valid] = 0.0
                k = int(np.argmax(error))
                if error[k] > stat.max_abs_error:
                    stat.max_abs_error = float(error[k])
                    stat.max_location = location(offset + k, nx, ny, a[k], b[k])
        for f, before in zip((ref_file, test_file), signatures):
            if file_signature(os.fstat(f.fileno())) != before:
                raise ValueError(f"File changed during comparison: {f.name}")
    return bitwise_identical, stats


def positive_int(text):
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return value


def nonnegative_float(text):
    value = float(text)
    if not math.isfinite(value) or value < 0:
        raise argparse.ArgumentTypeError("must be finite and nonnegative")
    return value


def grid_shape(text):
    match = re.fullmatch(r"([1-9][0-9]*)[xX]([1-9][0-9]*)", text)
    if match is None:
        raise argparse.ArgumentTypeError("expected NXxNY, for example 361x361")
    return tuple(map(int, match.groups()))


def validate_layout(reference, candidate, shapes, expected_snapshots):
    expected_names = {f"planedisplacements.{i}" for i in range(len(shapes))}
    for directory in (reference, candidate):
        if not directory.is_dir():
            raise ValueError(f"Not a directory: {directory}")
        names = {p.name for p in directory.iterdir()
                 if re.fullmatch(r"planedisplacements\.[0-9]+", p.name)}
        if names != expected_names:
            raise ValueError(f"{directory}: missing {sorted(expected_names - names)}; "
                             f"unexpected {sorted(names - expected_names)}")
    common_count = expected_snapshots
    for i, (nx, ny) in enumerate(shapes):
        for directory in (reference, candidate):
            path = directory / f"planedisplacements.{i}"
            count = snapshot_count(path, nx, ny)
            if common_count is None:
                common_count = count
            if count != common_count:
                raise ValueError(f"{path}: {count} snapshots; expected {common_count}")
    return common_count


def describe_location(value):
    frame, i, j, reference, candidate = value
    return (f"snapshot={frame}, point=({i},{j}), "
            f"reference={reference:.17g}, candidate={candidate:.17g}")


def print_result(name, snapshots, identical, stats):
    passed = all(stat.passed for stat in stats)
    exact = passed and not any(stat.finite_differences for stat in stats)
    print(f"{name}: {'PASS' if passed else 'FAIL'}; snapshots={snapshots}; "
          f"byte-identical={'yes' if identical else 'no'}; "
          f"exact finite values={'yes' if exact else 'no'}")
    print("  Comp   peak_ref(m)  peak_test(m)   max_abs(m)  err/ref_peak   "
          "finite_diff  outside_tol  nonfinite_ref/test")
    for axis, stat in zip("XYZ", stats):
        relative = (stat.max_abs_error / stat.reference_peak
                    if stat.reference_peak else
                    (0.0 if stat.max_abs_error == 0 else math.inf))
        print(f"  {axis}    {stat.reference_peak:12.5e} {stat.candidate_peak:12.5e} "
              f"{stat.max_abs_error:12.5e} {relative:12.5e} "
              f"{stat.finite_differences:12d} {stat.outside_tolerance:12d} "
              f"{stat.nonfinite_reference}/{stat.nonfinite_candidate}")
        if stat.max_location is not None:
            print(f"    max absolute error: {describe_location(stat.max_location)}")
        if stat.first_failure is not None:
            print(f"    first failure: {describe_location(stat.first_failure)}")
    return passed


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("reference", type=Path, help="legacy/reference output directory")
    parser.add_argument("candidate", type=Path, help="MPI-IO/candidate output directory")
    parser.add_argument("--shapes", nargs="+", required=True, type=grid_shape,
                        help="NXxNY for every plane, in plane-ID order starting at 0")
    parser.add_argument("--expected-snapshots", type=positive_int,
                        help="reject equally truncated files as well as count mismatches")
    parser.add_argument("--atol", type=nonnegative_float, default=0.0,
                        help="absolute displacement tolerance in metres (default: 0)")
    parser.add_argument("--rtol", type=nonnegative_float, default=0.0,
                        help="relative tolerance against each reference value (default: 0)")
    parser.add_argument("--chunk-points", type=positive_int, default=262144,
                        help="XYZ points per read, independent of file size (default: 262144)")
    parser.add_argument("--endian", choices=("native", "little", "big"), default="native",
                        help="byte order of both input sets (default: native)")
    args = parser.parse_args(argv)
    try:
        reference, candidate = args.reference.resolve(), args.candidate.resolve()
        if reference == candidate:
            raise ValueError("Reference and candidate directories must be different")
        snapshots = validate_layout(reference, candidate, args.shapes,
                                    args.expected_snapshots)
        print(f"Reference: {reference}\nCandidate: {candidate}")
        print(f"Criterion: abs(candidate-reference) <= {args.atol:g} m "
              f"+ {args.rtol:g} * abs(reference); NaN/Inf always fail.")
        print("Indices are zero-based: snapshot, along-strike point, down-dip point.\n"
              "err/ref_peak is a diagnostic; acceptance uses EACH reference value.")
        passed = True
        for i, (nx, ny) in enumerate(args.shapes):
            name = f"planedisplacements.{i}"
            print(f"\nComparing {name} ({nx}x{ny}, {snapshots} snapshots)...", flush=True)
            identical, stats = compare_plane(
                reference / name, candidate / name, nx, ny, snapshots,
                chunk_points=args.chunk_points, atol=args.atol, rtol=args.rtol,
                endian=args.endian)
            passed = print_result(name, snapshots, identical, stats) and passed
        print(f"\nOVERALL: {'PASS' if passed else 'FAIL'} for {len(args.shapes)} planes "
              "under the stated criterion.")
        return 0 if passed else 1
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
