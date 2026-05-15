#!/usr/bin/env python3
"""PML run self-check: is the new binary's output actually different?

Usage:
    python3 pml_fix_selfcheck.py <new_station_dir> [<baseline_station_dir>]

  <new_station_dir>       station.* output you just produced
  <baseline_station_dir>  a known OLD full_gauss run to diff against
                          (default: Output/stations_pml_codex_v9)

Decision logic:
  * Compares the new run's station.24 against the baseline run.
  * Computes the direct-arrival-window (10-22 s) RMS misfit vs the CPS
    analytic (which is ABC-independent there). A genuine quadrature /
    stiffness-mode change moves this by >=5%; floating-point round-off
    or a stale binary moves it <0.1%.
  * Prints "FIX LIVE" or "STALE / UNCHANGED" with the numbers.

No third-party deps beyond numpy. CPS reference path is hard-coded to
the 5-offset analytic run produced by
run_elastic_halfspace_no_atten_dipslip_5offsets.sh.
"""

import os
import sys
import struct
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CPS = os.path.join(
    HERE,
    "runs/elastic_halfspace_no_atten_dipslip_5offsets/displacement")
# station.24 = offset 16 km (00016000), depth 1500 m (001500)
CPS_R = os.path.join(CPS, "hr_1.5000", "00016000_005000_001500.R")
CPS_Z = os.path.join(CPS, "hr_1.5000", "00016000_005000_001500.Z")

DEF_BASELINE = "/Users/wzhang/hercules/Output/stations_pml_codex_v9"


def read_sac(path):
    raw = open(path, "rb").read()
    d = struct.unpack("f", raw[0:4])[0]
    b = struct.unpack("f", raw[20:24])[0]
    n = struct.unpack("i", raw[316:320])[0]
    return (b + np.arange(n) * d,
            np.frombuffer(raw[632:632 + n * 4], dtype="<f4").copy() * 0.01)


def herc24(d):
    a = np.loadtxt(os.path.join(d, "station.24"))
    # t, R(=Y east), Z(=-Z down)
    return a[:, 0], a[:, 2], -a[:, 3]


def direct_misfit(group_dir):
    tc, R = read_sac(CPS_R)
    _, Z = read_sac(CPS_Z)
    th, Rh, Zh = herc24(group_dir)
    Ri = np.interp(tc, th, Rh)
    Zi = np.interp(tc, th, Zh)
    w = (tc >= 10.0) & (tc <= 22.0)
    return float(np.sum((Ri[w] - R[w]) ** 2) + np.sum((Zi[w] - Z[w]) ** 2))


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    new_dir = sys.argv[1]
    base_dir = sys.argv[2] if len(sys.argv) > 2 else DEF_BASELINE

    if not os.path.isdir(new_dir):
        print(f"ERROR: {new_dir} not found")
        sys.exit(2)

    # Raw waveform difference vs the baseline run
    an = np.loadtxt(os.path.join(new_dir, "station.24"))
    have_base = os.path.isdir(base_dir)
    if have_base:
        ab = np.loadtxt(os.path.join(base_dir, "station.24"))
        n = min(len(an), len(ab))
        diff = np.sqrt(((an[:n, 1:4] - ab[:n, 1:4]) ** 2).sum(1))
        peak = float(np.sqrt((an[:, 1:4] ** 2).sum(1)).max())
        rel = diff.max() / peak if peak > 0 else 0.0
    else:
        rel = None

    dm_new = direct_misfit(new_dir)
    dm_base = direct_misfit(base_dir) if have_base else None

    print("=" * 64)
    print(f"new run     : {new_dir}")
    print(f"baseline    : {base_dir}"
          + ("" if have_base else "  (NOT FOUND)"))
    print("-" * 64)
    print(f"direct-window (10-22 s) sq-misfit vs CPS:")
    print(f"  new      = {dm_new:.6e}")
    if dm_base is not None:
        print(f"  baseline = {dm_base:.6e}")
        chg = 100.0 * (dm_new / dm_base - 1.0) if dm_base > 0 else 0.0
        print(f"  change   = {chg:+.2f}%")
    if rel is not None:
        print(f"max |new - baseline| / peak at station.24 = {rel:.2e}")
    print("-" * 64)

    # Verdict
    stale = False
    if rel is not None and rel < 5e-3:
        stale = True
    if dm_base is not None and abs(dm_new / dm_base - 1.0) < 0.01:
        stale = True

    if stale:
        print("VERDICT: STALE / UNCHANGED")
        print("  The new output is within round-off of the baseline run.")
        print("  The intended stiffness/quadrature change did NOT take")
        print("  effect. Likely a stale psolve binary or the run wrote")
        print("  to a different output_stations_directory.")
        print("  -> make clean && make in quake/forward, confirm psolve")
        print("     mtime is newer than the source edit, confirm the")
        print("     two-phase PSOLVE_BIN points at the rebuilt binary,")
        print("     and that this dir is where the run actually wrote.")
    else:
        print("VERDICT: FIX LIVE")
        print("  The new output differs from the baseline by more than")
        print("  round-off in the ABC-independent direct window, so the")
        print("  stiffness-mode change is genuinely active. Proceed to")
        print("  the full CPS comparison.")
    print("=" * 64)


if __name__ == "__main__":
    main()
