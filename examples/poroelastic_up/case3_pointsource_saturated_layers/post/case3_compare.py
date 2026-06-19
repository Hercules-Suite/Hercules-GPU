#!/usr/bin/env python3
"""Case 3 - point source under saturated layers: poroelastic vs undrained-equivalent elastic.

In the undrained limit (permeability = 0) the saturated poroelastic layers are
exactly elastic layers with the undrained Vp used in ../elastic_twin. So the
surface ground motion from the poroelastic run must match the elastic twin to
within discretization error. This script computes the per-offset normalized RMS
misfit of the surface velocity components and reports the pore-pressure traces.

Station columns (print_station_velocities=yes, accelerations=no, pressures=yes):
    t  ux uy uz  vx vy vz  [p]        (x=north, y=east, z=down)
The elastic twin has pressures off -> no trailing p column.

Usage:
    python3 case3_compare.py [poro_dir] [twin_dir]
defaults: ../output/stations  and  ../elastic_twin/output/stations
"""
import os
import sys
import numpy as np

N_STATIONS = 5
OFFSETS_M = [0, 200, 400, 600, 800]
TCOL = 0
VX, VY, VZ = 4, 5, 6     # velocity component columns (0-based)
P_COL = 7
NRMS_TOL = 0.10          # body-wave window tolerance (staggered scheme vs monolithic twin)
T_EARLY = 2.5            # s: direct body-wave window; coda/surface waves come later


def read_station(path):
    return np.loadtxt(path, comments="#", ndmin=2)


def nrms(a, b, mask=None):
    n = min(len(a), len(b))
    a, b = a[:n], b[:n]
    if mask is not None:
        a, b = a[mask[:n]], b[mask[:n]]
    if a.size == 0:
        return np.nan
    ref = np.sqrt(np.mean(b ** 2))
    if ref <= 0:
        return np.nan
    return np.sqrt(np.mean((a - b) ** 2)) / ref


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    poro = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "output", "stations")
    twin = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "..", "elastic_twin", "output", "stations")

    print("Poroelastic vs undrained-equivalent elastic twin. NRMS per component,")
    print(f"split into the direct body-wave window (t < {T_EARLY}s) and the full trace")
    print("(full includes the multiply-scattered coda / surface waves):\n")
    print("  offset(m) |   body-wave NRMS (vx vy vz)   |     full-trace NRMS (vx vy vz)   | peak|p|(Pa)")
    body_ok = True
    for k, off in enumerate(OFFSETS_M):
        pp = os.path.join(poro, f"station.{k}")
        tp = os.path.join(twin, f"station.{k}")
        if not (os.path.isfile(pp) and os.path.isfile(tp)):
            print(f"  {off:7d}    [missing station files]")
            body_ok = False
            continue
        a, b = read_station(pp), read_station(tp)
        em = a[:, TCOL] < T_EARLY                       # body-wave window mask
        ex, ey, ez = (nrms(a[:, c], b[:, c], em) for c in (VX, VY, VZ))
        fx, fy, fz = (nrms(a[:, c], b[:, c]) for c in (VX, VY, VZ))
        pmax = np.max(np.abs(a[:, P_COL])) if a.shape[1] > P_COL else float("nan")
        worst_body = np.nanmax([ex, ey, ez])
        ok = np.isfinite(worst_body) and worst_body < NRMS_TOL
        body_ok = body_ok and ok
        print(f"  {off:7d} | {ex:6.3f} {ey:6.3f} {ez:6.3f}  [{'ok' if ok else 'HI'}] "
              f"| {fx:6.3f} {fy:6.3f} {fz:6.3f}        | {pmax:9.3g}")

    print(f"\n-> Body-wave window: {'PASS' if body_ok else 'CHECK'} "
          f"(all components NRMS < {100*NRMS_TOL:.0f}% for t < {T_EARLY}s)")
    print("Interpretation:")
    print(" - body-wave NRMS small + full-trace large  -> coupling reproduces the undrained")
    print("   medium for direct waves; coda divergence is staggered-scheme dispersion (expected).")
    print(" - body-wave NRMS also large                -> genuine coupling discrepancy; investigate.")
    print("Reminder: the elastic twin needs the_threshold_Vp_over_Vs raised above the saturated")
    print("layers' undrained Vp/Vs, or it is clamped and is not the undrained-equivalent.")
    return 0 if body_ok else 1


if __name__ == "__main__":
    sys.exit(main())
