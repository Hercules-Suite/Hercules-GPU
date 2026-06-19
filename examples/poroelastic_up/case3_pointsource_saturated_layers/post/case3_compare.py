#!/usr/bin/env python3
"""Case 3 - point source under saturated layers: poroelastic vs undrained-equivalent elastic.

In the undrained limit (permeability = 0) the saturated poroelastic layers are exactly elastic
layers with the undrained Vp used in ../elastic_twin, so the surface ground motion should match.

Headline metric: the CROSS-CORRELATION LAG of poro vs twin over the direct body-wave window.
Unlike a relative-threshold first-break pick, the cross-correlation lag is amplitude-independent
(it finds the true time-shift regardless of a ~10% peak-amplitude difference) and sub-sample
(parabolic peak refinement). Kernel analysis predicts the discrete coupling reproduces the
undrained modulus to machine precision, with only ~0.05% (~0.1 ms over 400 m) dispersive lag at
the run's resolution -- so a correct implementation should show sub-ms lags here, NOT the +3..11 ms
the threshold picker reported. NRMS is kept for context.

Station columns (print_station_velocities=yes, accelerations=no, pressures=yes):
    t  ux uy uz  vx vy vz  [p]        (x=north, y=east, z=down)

Usage:
    python3 case3_compare.py [poro_dir] [twin_dir]
    python3 case3_compare.py --selftest      # verify the lag estimator (sign, sub-sample, amplitude)
"""
import os
import sys
import numpy as np

OFFSETS_M = [0, 200, 400, 600, 800]
TCOL = 0
VX, VY, VZ = 4, 5, 6
P_COL = 7
T_EARLY = 2.5            # s: direct body-wave window
CORR_MIN = 0.9          # peak correlation below this -> lag not trustworthy (node/small signal)
LAG_TOL_MS = 1.0        # sub-ms expected if the coupling matches the undrained twin


def read_station(path):
    return np.loadtxt(path, comments="#", ndmin=2)


def nrms(a, b, mask=None):
    n = min(len(a), len(b)); a, b = a[:n], b[:n]
    if mask is not None:
        a, b = a[mask[:n]], b[mask[:n]]
    ref = np.sqrt(np.mean(b ** 2)) if a.size else 0.0
    return np.sqrt(np.mean((a - b) ** 2)) / ref if ref > 0 else np.nan


def xcorr_lag_ms(poro, twin, dt):
    """Sub-sample cross-correlation lag of poro relative to twin (POSITIVE = poro arrives LATER),
    plus the normalized peak correlation. Amplitude-independent."""
    n = min(len(poro), len(twin))
    a = poro[:n] - np.mean(poro[:n])
    b = twin[:n] - np.mean(twin[:n])
    na, nb = np.sqrt(np.sum(a * a)), np.sqrt(np.sum(b * b))
    if na == 0 or nb == 0:
        return np.nan, np.nan
    c = np.correlate(a, b, mode="full")          # index n-1 == zero lag
    k = int(np.argmax(c))
    dp = 0.0
    if 0 < k < len(c) - 1:                        # parabolic sub-sample refinement
        denom = c[k - 1] - 2 * c[k] + c[k + 1]
        if denom != 0:
            dp = 0.5 * (c[k - 1] - c[k + 1]) / denom
    lag_samples = (k + dp) - (n - 1)             # poro delayed by +lag => peak at index n-1+lag
    return lag_samples * dt * 1000.0, c[k] / (na * nb)


def selftest():
    dt = 1e-3
    t = np.arange(0, 4.0, dt)
    base = np.sin(2 * np.pi * 2 * t) * np.exp(-((t - 1.5) / 0.3) ** 2)
    ok = True
    for true_ms in (-2.3, 0.0, 0.4, 5.0):
        shifted = np.interp(t - true_ms / 1000.0, t, base, left=0, right=0)  # delayed by true_ms
        poro = 1.12 * shifted + 0.02 * np.sin(2 * np.pi * 5 * t) * np.exp(-((t - 1.5) / 0.3) ** 2)
        lag, corr = xcorr_lag_ms(poro, base, dt)
        good = abs(lag - true_ms) < 0.15
        ok = ok and good
        print(f"  true={true_ms:+5.2f} ms  recovered={lag:+6.3f} ms  corr={corr:.3f}  "
              f"[{'ok' if good else 'FAIL'}]  (amp x1.12, +5Hz shape perturbation)")
    print(f"-> lag estimator self-test {'PASS' if ok else 'FAIL'} "
          f"(sign: +=poro later; amplitude-independent; sub-sample)")
    return 0 if ok else 1


def main():
    if "--selftest" in sys.argv:
        return selftest()
    here = os.path.dirname(os.path.abspath(__file__))
    poro = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "output", "stations")
    twin = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "..", "elastic_twin", "output", "stations")

    print("Poroelastic vs undrained-equivalent elastic twin, body-wave window (t < "
          f"{T_EARLY}s).")
    print("Cross-correlation lag (POSITIVE = poro later); corr<0.9 = unreliable (node/small signal).")
    print("Discrete-operator prediction: |lag| ~ 0.1 ms. The threshold picker reported +3..11 ms.\n")
    print("  offset(m) |     vx (lag, corr)     |     vy (lag, corr)     |     vz (lag, corr)")
    reliable = []
    for k, off in enumerate(OFFSETS_M):
        pp = os.path.join(poro, f"station.{k}")
        tp = os.path.join(twin, f"station.{k}")
        if not (os.path.isfile(pp) and os.path.isfile(tp)):
            print(f"  {off:7d}    [missing station files]")
            continue
        a, b = read_station(pp), read_station(tp)
        n = min(len(a), len(b))
        dt = float(np.median(np.diff(a[:n, TCOL])))
        m = a[:n, TCOL] < T_EARLY
        cells = []
        for c in (VX, VY, VZ):
            lag, corr = xcorr_lag_ms(a[m, c], b[m, c], dt)
            cells.append(f"{lag:+6.2f}ms {corr:4.2f}")
            if np.isfinite(corr) and corr >= CORR_MIN:
                reliable.append(abs(lag))
        print(f"  {off:7d} | {cells[0]:>20} | {cells[1]:>20} | {cells[2]:>20}")

    if reliable:
        med = float(np.median(reliable))
        print(f"\nMedian |lag| over reliable (corr>={CORR_MIN}) components = {med:.2f} ms")
        if med < LAG_TOL_MS:
            print(f"-> PASS: sub-ms lag confirms the +3..11 ms first-break diffs were picker artifacts;")
            print("   the poroelastic coupling matches the undrained-equivalent twin (no Vp deficit).")
            ok = True
        else:
            print(f"-> CHECK: lag exceeds {LAG_TOL_MS} ms on real signals -> a genuine residual to investigate.")
            ok = False
    else:
        print("\n[info] no reliable (high-correlation) components found.")
        ok = False

    # NRMS for context (body-wave window)
    print("\nContext - body-wave NRMS (vx vy vz):")
    for k, off in enumerate(OFFSETS_M):
        pp = os.path.join(poro, f"station.{k}"); tp = os.path.join(twin, f"station.{k}")
        if not (os.path.isfile(pp) and os.path.isfile(tp)):
            continue
        a, b = read_station(pp), read_station(tp)
        m = a[:, TCOL] < T_EARLY
        print(f"  {off:7d} | " + " ".join(f"{nrms(a[:, c], b[:, c], m):6.3f}" for c in (VX, VY, VZ)))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
