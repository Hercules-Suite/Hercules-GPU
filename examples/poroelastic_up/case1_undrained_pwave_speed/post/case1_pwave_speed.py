#!/usr/bin/env python3
"""Case 1 - undrained P-wave speed (Gassmann limit) check.

Vertical P pulse up a homogeneous saturated column (perm=0). In the undrained
limit the propagating speed must be

    Vp_u = sqrt(Vp_drained^2 + alpha^2 * M / rho),   M = K_f / porosity
         = sqrt(400^2 + 5.5e9/2000) = 1705.9 m/s   (for the bundled material).

The column is 256 m = ~6 undrained wavelengths at the 40 Hz source, so the wave
genuinely propagates (a too-short / too-low-frequency column moves quasi-static
and the speed is unmeasurable).

Two reference runs disambiguate the result:
  * elastic_twin  -> Vp = Vp_u = 1705.872 m/s  (undrained-equivalent)
  * drained_twin  -> Vp = 400 m/s              (skeleton, coupling absent)
poro should match the UNDRAINED twin. If it matches the DRAINED twin instead,
the Biot coupling is not delivering the undrained stiffening.

Station columns (velocities/accelerations/pressures = yes):
    t  ux uy uz  vx vy vz  ax ay az  [p]      (x=north, y=east, z=down)

Usage: python3 case1_pwave_speed.py [poro_dir] [undrained_twin_dir] [drained_twin_dir]
"""
import os
import sys
import numpy as np

SOURCE_DEPTH_M = 1024.0
STATION_DEPTHS_M = [0.0, 128.0, 256.0, 384.0, 512.0, 640.0, 768.0, 896.0, 1020.0]
VP_DRAINED, RHO, KF, POROSITY, ALPHA = 400.0, 2000.0, 2.2e9, 0.4, 1.0
M = KF / POROSITY
VP_U = np.sqrt(VP_DRAINED**2 + ALPHA**2 * M / RHO)

VZ_COL = 6
P_COL = 10
PICK_FRACTION = 0.05
SPEED_TOL = 0.03      # 3 %
NRMS_TOL = 0.05       # 5 %


def read_station(path):
    return np.loadtxt(path, comments="#", ndmin=2)


def pick_arrival(t, sig, frac=PICK_FRACTION):
    peak = np.max(np.abs(sig))
    if peak <= 0:
        return np.nan
    return t[int(np.argmax(np.abs(sig) >= frac * peak))]


def load_vz(stations_dir):
    """Return list of (t, vz, full_array) per station (None if missing)."""
    out = []
    for k in range(len(STATION_DEPTHS_M)):
        p = os.path.join(stations_dir, f"station.{k}")
        if not os.path.isfile(p):
            out.append(None)
            continue
        a = read_station(p)
        out.append((a[:, 0], a[:, VZ_COL], a))
    return out


def nrms_vs(poro, ref):
    vals = []
    for k in range(len(STATION_DEPTHS_M)):
        if poro[k] is None or ref[k] is None:
            continue
        vp, vr = poro[k][1], ref[k][1]
        n = min(len(vp), len(vr))
        r = np.sqrt(np.mean(vr[:n] ** 2))
        if r > 0:
            vals.append(np.sqrt(np.mean((vp[:n] - vr[:n]) ** 2)) / r)
    return np.array(vals)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    poro_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "output", "stations")
    und_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "..", "elastic_twin", "output", "stations")
    dr_dir = sys.argv[3] if len(sys.argv) > 3 else os.path.join(here, "..", "drained_twin", "output", "stations")

    print(f"Expected undrained Vp_u = {VP_U:.2f} m/s; drained Vp = {VP_DRAINED} m/s "
          f"(M=K_f/n={M:.3g} Pa)\n")

    poro = load_vz(poro_dir)
    # ---- first-arrival speed (poro run) ----
    dist, tarr = [], []
    print("Depth(m)  travel_dist(m)  t_arrival(s)")
    for k, z in enumerate(STATION_DEPTHS_M):
        if poro[k] is None:
            print(f"  {z:6.1f}      [missing]")
            continue
        t, vz, _ = poro[k]
        ta = pick_arrival(t, vz)
        dist.append(SOURCE_DEPTH_M - z)
        tarr.append(ta)
        print(f"  {z:6.1f}      {SOURCE_DEPTH_M - z:7.1f}       {ta:.6f}")
    dist, tarr = np.array(dist), np.array(tarr)
    good = np.isfinite(tarr) & (tarr > 0)
    speed_ok = False
    if good.sum() >= 2:
        vp_meas = np.polyfit(tarr[good], dist[good], 1)[0]
        err = abs(vp_meas - VP_U) / VP_U
        speed_ok = err < SPEED_TOL
        print(f"\nMeasured Vp (slope fit) = {vp_meas:.1f} m/s "
              f"(expected {VP_U:.1f}, error {100*err:.1f}%) -> "
              f"{'PASS' if speed_ok else 'FAIL'} (tol {100*SPEED_TOL:.0f}%)")

    # ---- which medium does poro match? ----
    nu = nrms_vs(poro, load_vz(und_dir))
    nd = nrms_vs(poro, load_vz(dr_dir))
    mean_u = float(np.mean(nu)) if nu.size else float("nan")
    mean_d = float(np.mean(nd)) if nd.size else float("nan")
    print(f"\nMean NRMS(vz) vs UNDRAINED twin (Vp={VP_U:.0f}) = {mean_u:.4f}")
    print(f"Mean NRMS(vz) vs DRAINED   twin (Vp={VP_DRAINED:.0f}) = {mean_d:.4f}")
    match_ok = np.isfinite(mean_u) and mean_u < NRMS_TOL
    if match_ok:
        verdict = "matches UNDRAINED twin -> coupling delivers alpha^2*M correctly. PASS"
    elif np.isfinite(mean_d) and mean_d < NRMS_TOL:
        verdict = "matches DRAINED twin -> Biot coupling NOT adding undrained stiffening. FAIL"
    elif np.isfinite(mean_u) and np.isfinite(mean_d) and mean_u < mean_d:
        verdict = "closer to undrained but not within tol -> partial coupling / dispersion. CHECK"
    else:
        verdict = "matches NEITHER cleanly -> investigate (setup, BC, or coupling). CHECK"
    print(f"  -> {verdict}")

    # pore-pressure sign diagnostic
    km = len(STATION_DEPTHS_M) // 2
    if poro[km] is not None and poro[km][2].shape[1] > P_COL:
        p = poro[km][2][:, P_COL]
        print(f"\nPore-pressure peak at z={STATION_DEPTHS_M[km]} m: "
              f"{p[int(np.argmax(np.abs(p)))]:.4g} Pa")
    return 0 if (speed_ok and match_ok) else 1


if __name__ == "__main__":
    sys.exit(main())
