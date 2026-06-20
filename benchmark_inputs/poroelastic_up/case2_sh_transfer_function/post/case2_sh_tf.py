#!/usr/bin/env python3
"""Case 2 - 1D SH site-response transfer function (diagnostic).

Reports the surface/input transfer function for BOTH the poroelastic and the
elastic runs, so we can separate a TEST problem from a SOLVER problem:

  * elastic run resonates at f0 but poro does not  -> poroelastic coupling is
    altering the (fluid-independent) shear response  -> SOLVER issue.
  * NEITHER run resonates (flat TF, no peak near f0) -> the 1D-column setup is
    wrong (e.g. absorbing side boundaries damp the laterally-uniform shear, or
    the base force isn't driving the column) -> switch to the incident
    plane-wave / DRM injector. NOT a poroelastic-specific problem.
  * both resonate at f0 and agree -> PASS.

Analytical SH layer-over-halfspace (Kramer/Roesset), surface / rock outcrop:
    TF(f) = 1 / |cos(k1 H) + i*alpha_z*sin(k1 H)|,  k1=2*pi*f/Vs1,
    alpha_z = (rho1*Vs1)/(rho2*Vs2);  f0 = Vs1/(4H);  peak ~ 1/alpha_z.

Station columns: t ux uy uz vx vy vz [...]; shear motion is the x (north) comp.

Usage: python3 case2_sh_tf.py [poro_dir] [elastic_dir] [stf_file]
"""
import os
import sys
import numpy as np

H, VS1, RHO1 = 20.0, 150.0, 1900.0
VS2, RHO2 = 450.0, 2200.0
ALPHA_Z = (RHO1 * VS1) / (RHO2 * VS2)
F0 = VS1 / (4.0 * H)
PEAK_EXPECTED = 1.0 / ALPHA_Z

VX_COL = 4
FBAND = (0.5, 12.0)
F0_TOL = 0.05
RESONANCE_AMP = 1.5          # a real site peak amplifies > this
PORO_ELASTIC_TOL = 0.03


def read_station(path):
    return np.loadtxt(path, comments="#", ndmin=2)


def read_stf(path):
    with open(path) as f:
        nt = int(f.readline().split()[0])
        f.readline()
        dt = float(f.readline().split()[0])
        vals = [float(x) for x in f.read().split()][:nt]
    return dt, np.asarray(vals)


def amp_spectrum(sig, dt, pad=8):
    sig = (sig - np.mean(sig)) * np.hanning(len(sig))
    nfft = int(2 ** np.ceil(np.log2(max(len(sig) * pad, 2))))
    return np.fft.rfftfreq(nfft, dt), np.abs(np.fft.rfft(sig, n=nfft))


def analytic_tf(f):
    k1 = 2 * np.pi * f / VS1
    return 1.0 / np.abs(np.cos(k1 * H) + 1j * ALPHA_Z * np.sin(k1 * H))


def tf_of(stations_dir, fi, Ai):
    """surface(station.0) vx / input spectrum -> (fs, tf, surf_rms, within_rms)."""
    s = read_station(os.path.join(stations_dir, "station.0"))
    t, vx = s[:, 0], s[:, VX_COL]
    dt = np.median(np.diff(t)) if len(t) > 1 else 5e-4
    fs, As = amp_spectrum(vx, dt)
    Ai_i = np.interp(fs, fi, Ai)
    tf = np.where(Ai_i > 0, As / Ai_i, 0.0)
    surf_rms = float(np.sqrt(np.mean(vx ** 2)))
    within_rms = np.nan
    wp = os.path.join(stations_dir, "station.1")
    if os.path.isfile(wp):
        w = read_station(wp)
        within_rms = float(np.sqrt(np.mean(w[:, VX_COL] ** 2)))
    return fs, tf, surf_rms, within_rms


def peak_in_band(fs, tf):
    b = (fs >= FBAND[0]) & (fs <= FBAND[1])
    fpk = fs[b][np.argmax(tf[b])]
    return fpk, float(np.max(tf[b]))


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    poro = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "poroelastic", "output", "stations")
    elas = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "..", "elastic", "output", "stations")
    stf = sys.argv[3] if len(sys.argv) > 3 else os.path.join(here, "..", "sourcefiles", "baseforcefunction.in")

    print(f"Analytical SH: f0 = {F0:.3f} Hz, alpha_z = {ALPHA_Z:.3f}, peak ~ {PEAK_EXPECTED:.2f}\n")
    dt_i, vi = read_stf(stf)
    fi, Ai = amp_spectrum(vi, dt_i)

    results = {}
    for tag, d in (("ELASTIC", elas), ("PORO", poro)):
        if not os.path.isfile(os.path.join(d, "station.0")):
            print(f"{tag}: station.0 missing ({d})")
            continue
        fs, tf, srms, wrms = tf_of(d, fi, Ai)
        fpk, apk = peak_in_band(fs, tf)
        amp_ratio = (srms / wrms) if (wrms and np.isfinite(wrms) and wrms > 0) else float("nan")
        resonates = (abs(fpk - F0) / F0 < F0_TOL) and (apk > RESONANCE_AMP)
        results[tag] = dict(fpk=fpk, apk=apk, srms=srms, wrms=wrms,
                            amp_ratio=amp_ratio, resonates=resonates, fs=fs, tf=tf)
        print(f"{tag:8s}: TF peak f={fpk:6.3f} Hz (f0={F0:.3f}), TF amp={apk:7.3f}, "
              f"surf/within RMS={amp_ratio:6.2f}, resonates={resonates}")

    if "ELASTIC" in results and "PORO" in results:
        e, p = results["ELASTIC"], results["PORO"]
        # poro vs elastic surface time series
        es = read_station(os.path.join(elas, "station.0"))[:, VX_COL]
        ps = read_station(os.path.join(poro, "station.0"))[:, VX_COL]
        n = min(len(es), len(ps)); ref = np.sqrt(np.mean(es[:n] ** 2))
        nrms = np.sqrt(np.mean((ps[:n] - es[:n]) ** 2)) / ref if ref > 0 else np.inf
        print(f"\nPoro vs elastic surface vx NRMS = {nrms:.4f} "
              f"({'PASS' if nrms < PORO_ELASTIC_TOL else 'FAIL'}, tol {100*PORO_ELASTIC_TOL:.0f}%)")
        print("\nVERDICT:")
        if e["resonates"] and p["resonates"] and nrms < PORO_ELASTIC_TOL:
            print("  Both runs resonate at f0 and agree -> PASS.")
        elif e["resonates"] and not p["resonates"]:
            print("  Elastic resonates but poroelastic does NOT -> the poroelastic "
                  "coupling is altering the shear response. SOLVER issue to investigate.")
        elif not e["resonates"]:
            print("  Elastic run does NOT resonate (flat TF) -> the 1D-column TEST is "
                  "wrong: absorbing side boundaries are damping the laterally-uniform "
                  "shear, or the base force isn't driving the column. Switch to the "
                  "incident plane-wave/DRM injector (a wide domain won't help a thin "
                  "column). This is NOT a poroelastic-specific failure.")
        else:
            print("  Both resonate near f0 but differ -> check amplitude/phase; "
                  "likely a moderate poro-vs-elastic discrepancy worth a closer look.")
    print(f"\n(Plot tf vs analytic_tf over {FBAND[0]}-{FBAND[1]} Hz for the full curve.)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
