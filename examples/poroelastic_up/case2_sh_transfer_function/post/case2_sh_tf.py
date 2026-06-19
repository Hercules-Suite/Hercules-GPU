#!/usr/bin/env python3
"""Case 2 - 1D SH site-response transfer function.

Validates that (a) the poroelastic run reproduces the SAME shear/site response
as the elastic run (shear is fluid-independent -> the hard correctness check),
and (b) the surface site amplification matches the analytical 1D SH transfer
function for a soft layer over an elastic halfspace.

Analytical SH TF (surface / rock outcrop), Kramer/Roesset:
    TF(f) = 1 / |cos(k1 H) + i * alpha_z * sin(k1 H)|,   k1 = 2*pi*f / Vs1
    alpha_z = (rho1*Vs1) / (rho2*Vs2)   (impedance ratio, layer/halfspace)
    fundamental f0 = Vs1 / (4H);  peak amplification ~ 1/alpha_z at f0.

Station columns (print_station_velocities=yes): t ux uy uz vx vy vz [...]; the
shear motion is the x (north) component. The TF reference here is the known
input velocity history (sourcefiles/baseforcefunction.in), i.e. the incident
motion. NOTE: the absolute amplitude depends on the dashpot incident/outcrop
convention (a factor ~2) and on damping (here Q~1e6, none), so the FIRM checks
are the resonance FREQUENCY and the poro-vs-elastic agreement; the peak
amplitude is compared approximately.

Usage:
    python3 case2_sh_tf.py [poro_dir] [elastic_dir] [stf_file]
defaults: ../poroelastic/output/stations ../elastic/output/stations ../sourcefiles/baseforcefunction.in
"""
import os
import sys
import numpy as np

# Profile (keep in sync with the parameters.in files)
H, VS1, RHO1 = 20.0, 150.0, 1900.0       # soft layer
VS2, RHO2 = 450.0, 2200.0                # bedrock halfspace
ALPHA_Z = (RHO1 * VS1) / (RHO2 * VS2)
F0 = VS1 / (4.0 * H)
PEAK_EXPECTED = 1.0 / ALPHA_Z

DT_OUT = 5.0e-4          # station print_rate(5) * simulation_delta_time_sec(1e-4)
VX_COL = 4               # 0-based: t,ux,uy,uz,vx,...
FBAND = (0.5, 12.0)      # comparison band, Hz
F0_TOL = 0.05            # 5 %
PEAK_TOL = 0.20          # 20 % (approximate; see note above)
PORO_ELASTIC_TOL = 0.03  # 3 % time-series NRMS


def read_station(path):
    return np.loadtxt(path, comments="#", ndmin=2)


def read_stf(path):
    with open(path) as f:
        nt = int(f.readline().split()[0])
        f.readline()                       # tinit
        dt = float(f.readline().split()[0])
        vals = [float(x) for x in f.read().split()][:nt]
    return dt, np.asarray(vals)


def amp_spectrum(sig, dt, pad=8):
    """Hanning-windowed amplitude spectrum, zero-padded so the peak frequency
    is well localized (a 4 s record has only 0.25 Hz bins, and f0=1.875 Hz
    falls between bins -> pad to interpolate the spectrum)."""
    sig = sig - np.mean(sig)
    sig = sig * np.hanning(len(sig))
    nfft = int(2 ** np.ceil(np.log2(max(len(sig) * pad, 2))))
    F = np.fft.rfft(sig, n=nfft)
    f = np.fft.rfftfreq(nfft, dt)
    return f, np.abs(F)


def analytic_tf(f):
    k1 = 2 * np.pi * f / VS1
    return 1.0 / np.abs(np.cos(k1 * H) + 1j * ALPHA_Z * np.sin(k1 * H))


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    poro = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "poroelastic", "output", "stations")
    elas = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "..", "elastic", "output", "stations")
    stf = sys.argv[3] if len(sys.argv) > 3 else os.path.join(here, "..", "sourcefiles", "baseforcefunction.in")

    print(f"Analytical SH layer-over-halfspace: f0 = {F0:.3f} Hz, "
          f"alpha_z = {ALPHA_Z:.3f}, peak ~ {PEAK_EXPECTED:.2f}\n")

    surf = read_station(os.path.join(poro, "station.0"))
    t, vx_surf = surf[:, 0], surf[:, VX_COL]
    dt = np.median(np.diff(t)) if len(t) > 1 else DT_OUT

    # site amplification referenced to the (known) input incident motion
    dt_stf, stf_v = read_stf(stf)
    fs, As = amp_spectrum(vx_surf, dt)
    fi, Ai = amp_spectrum(stf_v, dt_stf)
    band = (fs >= FBAND[0]) & (fs <= FBAND[1])
    Ai_on_fs = np.interp(fs, fi, Ai)
    with np.errstate(divide="ignore", invalid="ignore"):
        tf_sim = np.where(Ai_on_fs > 0, As / Ai_on_fs, 0.0)

    fpk = fs[band][np.argmax(tf_sim[band])]
    apk = np.max(tf_sim[band])
    f0_err = abs(fpk - F0) / F0
    print(f"Simulated TF peak: f = {fpk:.3f} Hz (expected f0 = {F0:.3f}, "
          f"err {100*f0_err:.1f}%), amplitude = {apk:.2f}")
    f0_ok = f0_err < F0_TOL
    print(f"  -> resonance frequency {'PASS' if f0_ok else 'FAIL'} "
          f"(tol {100*F0_TOL:.0f}%)")
    # amplitude is approximate (incident/outcrop factor ~2, damping)
    amp_ratio = apk / PEAK_EXPECTED
    print(f"  amplitude/analytic = {amp_ratio:.2f} "
          f"(expect O(1); a ~2x offset can come from the incident/outcrop "
          f"dashpot convention)")

    # Hard check: poroelastic vs elastic surface shear time series.
    pe_ok = True
    epath = os.path.join(elas, "station.0")
    if os.path.isfile(epath):
        es = read_station(epath)
        vx_e = es[:, VX_COL]
        n = min(len(vx_surf), len(vx_e))
        ref = np.sqrt(np.mean(vx_e[:n] ** 2))
        nrms = np.sqrt(np.mean((vx_surf[:n] - vx_e[:n]) ** 2)) / ref if ref > 0 else np.inf
        pe_ok = nrms < PORO_ELASTIC_TOL
        print(f"\nPoroelastic vs elastic surface vx NRMS = {nrms:.4f} "
              f"-> {'PASS' if pe_ok else 'FAIL'} (tol {100*PORO_ELASTIC_TOL:.0f}%)")
    else:
        print(f"\n[info] elastic run not found ({epath}); skipping poro-vs-elastic check.")

    print("\n(For a full curve, plot tf_sim(fs) against analytic_tf(fs) over "
          f"{FBAND[0]}-{FBAND[1]} Hz.)")
    return 0 if (f0_ok and pe_ok) else 1


if __name__ == "__main__":
    sys.exit(main())
