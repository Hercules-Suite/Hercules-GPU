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
VX, VY, VZ = 4, 5, 6     # velocity component columns (0-based)
P_COL = 7
NRMS_TOL = 0.03          # 3 %


def read_station(path):
    return np.loadtxt(path, comments="#", ndmin=2)


def nrms(a, b):
    n = min(len(a), len(b))
    ref = np.sqrt(np.mean(b[:n] ** 2))
    if ref <= 0:
        return np.nan
    return np.sqrt(np.mean((a[:n] - b[:n]) ** 2)) / ref


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    poro = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "output", "stations")
    twin = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "..", "elastic_twin", "output", "stations")

    print("Per-offset NRMS, poroelastic vs undrained-equivalent elastic twin:")
    print("  offset(m)   NRMS(vx)   NRMS(vy)   NRMS(vz)   peak|p|(Pa)")
    all_ok = True
    for k, off in enumerate(OFFSETS_M):
        pp = os.path.join(poro, f"station.{k}")
        tp = os.path.join(twin, f"station.{k}")
        if not (os.path.isfile(pp) and os.path.isfile(tp)):
            print(f"  {off:7d}    [missing station files]")
            all_ok = False
            continue
        a, b = read_station(pp), read_station(tp)
        nx, ny, nz = nrms(a[:, VX], b[:, VX]), nrms(a[:, VY], b[:, VY]), nrms(a[:, VZ], b[:, VZ])
        pmax = np.max(np.abs(a[:, P_COL])) if a.shape[1] > P_COL else float("nan")
        worst = np.nanmax([nx, ny, nz])
        ok = np.isfinite(worst) and worst < NRMS_TOL
        all_ok = all_ok and ok
        flag = "ok" if ok else "HIGH"
        print(f"  {off:7d}    {nx:7.4f}    {ny:7.4f}    {nz:7.4f}    {pmax:10.4g}  [{flag}]")

    print(f"\n-> Case 3 {'PASS' if all_ok else 'FAIL'} "
          f"(all components NRMS < {100*NRMS_TOL:.0f}% at every offset)")
    print("Note: the poro-vs-twin comparison is insensitive to absorbing-boundary "
          "reflections (both runs share them). For the optional CPS reflectivity "
          "cross-check, use a larger domain / longer record to keep the boundary "
          "out of the comparison window.")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
