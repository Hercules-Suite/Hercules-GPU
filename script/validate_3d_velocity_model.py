#!/usr/bin/env python3
"""Validate count-prefixed Hercules 3-D model files; requires NumPy."""

import argparse
from pathlib import Path

import numpy as np


def counted(root, name):
    raw = np.array((root / name).read_text().split(), dtype=float)
    if raw.size == 0 or not np.isfinite(raw).all():
        raise ValueError(f"{name}: empty file or non-finite values")
    if raw[0] < 0 or raw[0] != int(raw[0]) or raw.size - 1 != raw[0]:
        raise ValueError(f"{name}: invalid count header or mismatched value count")
    return raw[1:]


def validate(root):
    root = Path(root)
    x, y, z, starts, vs, vp, rho = (
        counted(root, name + ".in") for name in ("x", "y", "z", "index", "vs", "vp", "rho")
    )
    for name, axis in (("x", x), ("y", y)):
        steps = np.diff(axis)
        if axis.size < 2 or np.any(steps <= 0):
            raise ValueError(f"{name}.in: need at least two increasing coordinates")
        if not np.allclose(steps, steps[0], rtol=1e-9, atol=1e-8):
            raise ValueError(f"{name}.in: spacing is not uniform")
    if not (z.size == vs.size == vp.size == rho.size) or z.size == 0:
        raise ValueError("z, Vs, Vp, and density counts must match and be nonzero")
    if (starts.size != x.size * y.size or starts[0] != 0
            or np.any(starts != np.floor(starts)) or np.any(np.diff(starts) <= 0)
            or starts[-1] >= z.size):
        raise ValueError("index.in does not define one ordered profile per x-y column")
    starts = starts.astype(np.int64)
    for column, (begin, end) in enumerate(zip(starts, np.r_[starts[1:], z.size])):
        if np.any(np.diff(z[begin:end]) <= 0):
            raise ValueError(f"vertical profile {column} is not strictly increasing")
    if np.any(z < 0) or np.any(vs < 0) or np.any(vp <= 0) or np.any(rho <= 0):
        raise ValueError("require z >= 0, Vs >= 0, Vp > 0, and density > 0")
    has_qs, has_qp = (root / "qs.in").exists(), (root / "qp.in").exists()
    if has_qs != has_qp:
        raise ValueError("provide both qs.in and qp.in, or neither")
    if has_qs:
        for name in ("qs.in", "qp.in"):
            quality = counted(root, name)
            if quality.size != z.size or np.any(quality <= 0):
                raise ValueError(f"{name}: need one positive Q value per material sample")
    print(f"valid model: {x.size} x {y.size} columns, {z.size} material samples")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="Directory containing x.in, y.in, etc.")
    args = parser.parse_args()
    validate(args.directory)


if __name__ == "__main__":
    main()
