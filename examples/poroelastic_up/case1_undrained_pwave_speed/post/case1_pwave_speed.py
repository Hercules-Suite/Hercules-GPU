#!/usr/bin/env python3
"""Case 1 - undrained P-wave speed (Gassmann limit) check.

Reads the Hercules station traces from a poroelastic run (and, if present, its
elastic twin), picks the P first arrival down the vertical receiver line, and
fits the propagation speed. In the undrained limit (permeability = 0) the
measured speed must equal

    Vp_u = sqrt(Vp_drained^2 + alpha^2 * M / rho),   M = K_f / porosity.

For the bundled material (Vp_d=400, rho=2000, K_f=2.2e9, n=0.4, alpha=1):
    Vp_u = sqrt(400**2 + 5.5e9/2000) = 1705.9 m/s.

Station columns (print_station_velocities/accelerations/pressures = yes):
    t  ux uy uz  vx vy vz  ax ay az  [p]      (x=north, y=east, z=down)

Usage:
    python3 case1_pwave_speed.py [poro_stations_dir] [twin_stations_dir]
defaults: ../output/stations  and  ../elastic_twin/output/stations
"""
import os
import sys
import numpy as np

# Geometry / material (keep in sync with ../parameters.in)
SOURCE_DEPTH_M = 32.0
STATION_DEPTHS_M = [0.0, 2.0, 6.0, 10.0, 14.0, 18.0, 22.0, 27.0, 31.5]
VP_DRAINED, RHO, KF, POROSITY, ALPHA = 400.0, 2000.0, 2.2e9, 0.4, 1.0

M = KF / POROSITY
VP_U_EXPECTED = np.sqrt(VP_DRAINED**2 + ALPHA**2 * M / RHO)

VZ_COL = 6          # 0-based: t,ux,uy,uz,vx,vy,vz,...
P_COL = 10
PICK_FRACTION = 0.05  # first-arrival threshold as a fraction of the trace peak
SPEED_TOL = 0.02      # 2 %
NRMS_TOL = 0.02       # 2 %


def read_station(path):
    a = np.loadtxt(path, comments="#", ndmin=2)
    return a


def pick_arrival(t, sig, frac=PICK_FRACTION):
    """First time |sig| exceeds frac * max|sig|."""
    peak = np.max(np.abs(sig))
    if peak <= 0:
        return np.nan
    idx = np.argmax(np.abs(sig) >= frac * peak)
    return t[idx]


def collect(stations_dir):
    t_arr, dist = [], []
    traces = []
    for k, z in enumerate(STATION_DEPTHS_M):
        path = os.path.join(stations_dir, f"station.{k}")
        if not os.path.isfile(path):
            print(f"  [warn] missing {path}")
            traces.append(None)
            continue
        a = read_station(path)
        t, vz = a[:, 0], a[:, VZ_COL]
        traces.append((t, vz, a))
        ta = pick_arrival(t, vz)
        t_arr.append(ta)
        dist.append(SOURCE_DEPTH_M - z)
    return np.array(dist), np.array(t_arr), traces


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    poro_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "output", "stations")
    twin_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "..", "elastic_twin", "output", "stations")

    print(f"Expected undrained Vp_u = {VP_U_EXPECTED:.2f} m/s "
          f"(Vp_drained={VP_DRAINED}, M=K_f/n={M:.3g} Pa)\n")

    dist, t_arr, traces = collect(poro_dir)
    good = np.isfinite(t_arr) & (t_arr > 0)
    if good.sum() < 2:
        print("ERROR: too few first-arrival picks; check the run output.")
        return 1
    # distance = Vp * traveltime  ->  slope of dist vs t_arr is Vp
    slope, intercept = np.polyfit(t_arr[good], dist[good], 1)
    vp_measured = slope
    err = abs(vp_measured - VP_U_EXPECTED) / VP_U_EXPECTED
    print("Depth(m)  travel_dist(m)  t_arrival(s)")
    for z, d, ta in zip(STATION_DEPTHS_M, dist, t_arr):
        print(f"  {z:6.1f}      {d:7.2f}       {ta:.6f}")
    print(f"\nMeasured Vp (slope fit) = {vp_measured:.2f} m/s "
          f"(expected {VP_U_EXPECTED:.2f}, error {100*err:.2f}%)")
    speed_ok = err < SPEED_TOL
    print(f"  -> P-wave speed {'PASS' if speed_ok else 'FAIL'} "
          f"(tol {100*SPEED_TOL:.0f}%)")

    # Self-consistency vs the elastic twin (same vz time base, same sampling).
    nrms_ok = True
    if os.path.isdir(twin_dir):
        print("\nNRMS(vz) poroelastic vs elastic twin, per station:")
        td, tt, twin_traces = collect(twin_dir)
        for k, z in enumerate(STATION_DEPTHS_M):
            if traces[k] is None or twin_traces[k] is None:
                continue
            _, vz_p, _ = traces[k]
            tt_t, vz_t, _ = twin_traces[k]
            n = min(len(vz_p), len(vz_t))
            ref = np.sqrt(np.mean(vz_t[:n] ** 2))
            if ref <= 0:
                continue
            nrms = np.sqrt(np.mean((vz_p[:n] - vz_t[:n]) ** 2)) / ref
            flag = "ok" if nrms < NRMS_TOL else "HIGH"
            nrms_ok = nrms_ok and (nrms < NRMS_TOL)
            print(f"  z={z:6.1f} m   NRMS={nrms:8.4f}  [{flag}]")
        print(f"  -> twin self-consistency {'PASS' if nrms_ok else 'FAIL'} "
              f"(tol {100*NRMS_TOL:.0f}%)")
    else:
        print(f"\n[info] elastic twin dir not found ({twin_dir}); skipping NRMS.")

    # Pore-pressure sign diagnostic (vertical compression -> positive excess p).
    k_mid = len(STATION_DEPTHS_M) // 2
    if traces[k_mid] is not None and traces[k_mid][2].shape[1] > P_COL:
        p = traces[k_mid][2][:, P_COL]
        ip = int(np.argmax(np.abs(p)))
        print(f"\nPore-pressure sign check at z={STATION_DEPTHS_M[k_mid]} m: "
              f"peak p = {p[ip]:.4g} Pa (expect compression -> sign consistent "
              f"with vertical loading convention)")

    return 0 if (speed_ok and nrms_ok) else 1


if __name__ == "__main__":
    sys.exit(main())
