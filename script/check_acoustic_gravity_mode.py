#!/usr/bin/env python3
"""Check the production Hex8 gravity operator against a flat-water mode.

Requires NumPy and a host C++ compiler (CXX is honored); no CUDA or MPI.
Extracts the production stiffness and mass expressions, assembles a periodic
horizontal Fourier mode above a rigid bottom, and solves its vertical
generalized eigenproblem. This isolates uniform-grid, lossless, non-CPML
surface dispersion and gradient recovery, not the full coupled GPU solver.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile

import numpy as np


def function(source, signature):
    start = source.index(signature)
    return source[start:source.index("\n}\n", start) + 3]


def production_element(root):
    kernel = (root / "quake/forward/kernel.cu").read_text()
    solver = (root / "quake/forward/psolve.cc").read_text()
    header = (root / "quake/forward/psolve.h").read_text()
    gravity_define = re.search(r"^#define ACOUSTIC_GRAVITY_ACCELERATION_M_S2\s+[^\n]+",
                               header, flags=re.M).group(0)
    functions = "\n".join(function(kernel, signature) for signature in (
        "__host__ __device__ static inline double acoustic_node_sign(",
        "__device__ static inline void acoustic_scalar_stiffness_sums(",
        "__device__ static inline double acoustic_scalar_stiffness_from_sums("))
    start = solver.index("      const double bulk = edata->rho * edata->Vp * edata->Vp;",
                         solver.index("static void solver_init() {"))
    mass = solver[start:solver.index("      ep->c1 =", start)]
    start = solver.index("      acousticPressureCoeffArray[i] =")
    prefactor = solver[start:solver.index(";", start) + 1].replace(
        "acousticPressureCoeffArray[i]", "const double prefactor")
    cpp = r'''
#include <stdio.h>
#define __host__
#define __device__
#define YES 1
''' + gravity_define + r'''
typedef double solver_float;
struct Material { double rho, Vp, edgesize; };
static struct { int theAcousticFreeSurfaceGravity; double theDeltaTSquared; }
    Param = {YES, 0.04};
static int acoustic_element_has_free_surface(void*, double) { return 1; }
''' + functions + r'''
int main() {
    Material material = {1025.0, 1500.0, 2.0};
    Material* edata = &material;
    void* elemp = nullptr;
    const double acousticSurfaceZ = 0;
''' + mass + prefactor + r'''
    const double h = edata->edgesize;
    printf("%.17g %.17g\n", nodalPressureMass*edata->rho*edata->Vp*edata->Vp/(h*h*h),
           gravitySurfaceMass*edata->rho/(h*h));
    for (int j = 0; j < 8; ++j) {
        double p[8] = {}, sums[7]; p[j] = 1;
        acoustic_scalar_stiffness_sums(p, sums);
        for (int i = 0; i < 8; ++i)
            printf("%.17g%c", prefactor * acoustic_scalar_stiffness_from_sums(i, sums)
                   * edata->rho/(Param.theDeltaTSquared*h),
                   i == 7 ? '\n' : ' ');
    }
}
'''
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory)
        (path / "element.cc").write_text(cpp)
        subprocess.run(shlex.split(os.environ.get("CXX", "c++")) +
                       ["-std=c++11", "-O2", "-Wall", "-Wextra", "-Werror",
                        str(path / "element.cc"), "-o", str(path / "element")], check=True)
        output = subprocess.check_output([str(path / "element")], text=True).splitlines()
    volume_mass, surface_mass = map(float, output[0].split())
    stiffness = np.array([[float(v) for v in line.split()] for line in output[1:]]).T
    # Independent exact tensor-product Hex8 integral on a unit cube.
    m = np.array([[2, 1], [1, 2]]) / 6
    k = np.array([[1, -1], [-1, 1]])
    exact = np.kron(m, np.kron(m, k)) + np.kron(m, np.kron(k, m)) + np.kron(k, np.kron(m, m))
    if not np.allclose(stiffness, exact, rtol=0, atol=2e-15):
        raise ValueError("Production acoustic stiffness differs from exact Hex8 quadrature")
    if not np.isclose(volume_mass, 1/8, rtol=0, atol=1e-15) or not np.isclose(
            surface_mass, 1/(4*9.81), rtol=0, atol=1e-15):
        raise ValueError("Production volume/surface mass normalization differs")
    return stiffness, volume_mass, surface_mass


def analytical_omega2(k, depth, c, gravity=9.81):
    # For phi ~ cosh(mu*(H-z)), mu²=k²-omega²/c². Combining the
    # rigid-bottom condition and phi_tt=g*phi_z gives omega²=g*mu*tanh(mu*H).
    lo, hi = 0., min(gravity*k*math.tanh(k*depth), (c*k)**2)
    for _ in range(90):
        value = .5*(lo+hi)
        mu = math.sqrt(max(0., k*k-value/(c*c)))
        if value > gravity*mu*math.tanh(mu*depth):
            hi = value
        else:
            lo = value
    return .5*(lo+hi)


def mode(element, volume_mass, surface_mass, wavelength, depth, layers, c=1500.):
    h = depth/layers
    k = 2*math.pi/wavelength
    if k*h >= math.pi:
        raise ValueError("Horizontal wavelength must exceed two element widths")
    # One horizontal cell, Bloch periodic in X and constant in Y. Identifying
    # its corner DOFs assembles all sharing contributions without a full grid.
    projection = np.zeros((8, 2), dtype=complex)
    for n in range(8):
        projection[n, n >> 2] = np.exp(1j*k*h*(n & 1))
    local_k = projection.conj().T @ (h*element) @ projection
    local_m = np.sum(np.abs(projection)**2, axis=0)*volume_mass*h**3/c**2
    stiffness = np.zeros((layers+1, layers+1), dtype=complex)
    mass = np.zeros(layers+1)
    for i in range(layers):
        stiffness[i:i+2, i:i+2] += local_k
        mass[i:i+2] += local_m
    mass[0] += 4*surface_mass*h*h
    inv_sqrt = 1/np.sqrt(mass)
    eigenvalues, vectors = np.linalg.eigh(stiffness*inv_sqrt[:, None]*inv_sqrt[None, :])
    omega2 = float(eigenvalues[0])
    if omega2 <= 0:
        raise ValueError("Non-positive gravity eigenvalue")
    phi = vectors[:, 0]*inv_sqrt
    phi /= phi[0]
    if np.max(np.abs(phi.imag)) > 1e-9:
        raise ValueError("Unexpected complex vertical eigenfunction")
    phi = phi.real
    exact = analytical_omega2(k, depth, c)
    # At an X grid node, the one-sided Hex8 recovered vertical derivative.
    # For a pure mode uZ / [-p/(rho*g)] = -g*(dphi/dz)/(omega²*phi_surface).
    recovered_ratio = -9.81*(phi[1]-phi[0])/(h*omega2)
    residual = stiffness@phi - omega2*mass*phi
    relative_residual = float(np.linalg.norm(residual)/np.linalg.norm(stiffness@phi))
    # Time-integrated potential after a unit sinusoidal permanent bottom
    # displacement. A consistent bottom face load gives K*psi = B*u_bottom.
    # Continuum surface displacement is zero at static equilibrium, but
    # horizontal particle displacement need not vanish. The first-element
    # vertical derivative can retain a nonzero recovery error.
    bottom_load = np.zeros(layers+1)
    bottom_load[-1] = h*h*(2+math.cos(k*h))/3
    integrated_phi = np.linalg.solve(stiffness, bottom_load).real
    static_z = (integrated_phi[1]-integrated_phi[0])/h
    static_x = 2*math.sin(k*h/2)/h*integrated_phi[0]
    return {"layers": layers, "h_m": h, "wavelength_m": wavelength,
            "cells_per_wavelength": wavelength/h, "depth_m": depth,
            "frequency_hz": math.sqrt(omega2)/(2*math.pi),
            "analytic_frequency_hz": math.sqrt(exact)/(2*math.pi),
            "frequency_relative_error": math.sqrt(omega2/exact)-1,
            "surface_gradient_motion_ratio": float(recovered_ratio),
            "unit_bottom_static_surface_recovered_uz_m": float(static_z),
            "unit_bottom_static_surface_recovered_ux_amplitude_m": float(static_x),
            "unit_bottom_static_surface_exact_ux_amplitude_m": 1/math.sinh(k*depth),
            "eigen_residual_relative": relative_residual}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depth", type=float, default=175.78125)
    parser.add_argument("--wavelengths", type=float, nargs="+", default=[250., 500., 1000., 2000.])
    parser.add_argument("--layers", type=int, nargs="+", default=[8, 16, 32, 64])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not math.isfinite(args.depth) or args.depth <= 0 or any(n < 1 for n in args.layers) or any(
            not math.isfinite(w) or w <= 0 for w in args.wavelengths):
        parser.error("Depth, wavelengths and layer counts must be finite and positive")
    if args.output.exists():
        parser.error("Output already exists; choose a new report path")
    root = Path(__file__).resolve().parents[1]
    element, volume_mass, surface_mass = production_element(root)
    rows = [mode(element, volume_mass, surface_mass, w, args.depth, n)
            for w in args.wavelengths for n in args.layers]
    report = {"scope": "Production Hex8 stiffness and mass, uniform flat water, periodic X, constant Y, "
              "rigid bottom, lossless, no CPML or solid coupling; not a full solver validation.",
              "production_source_sha256": {str(p): hashlib.sha256((root/p).read_bytes()).hexdigest()
                  for p in (Path("quake/forward/kernel.cu"), Path("quake/forward/psolve.cc"),
                            Path("quake/forward/psolve.h"))},
              "numpy_version": np.__version__, "compiler": os.environ.get("CXX", "c++"),
              "reference": "omega^2 = g*mu*tanh(mu*H), mu^2 = k^2-omega^2/c^2; c=1500 m/s",
              "static_reference": "Unit sinusoidal permanent bottom displacement: equilibrium surface uZ=0, "
              "surface horizontal displacement amplitude=1/sinh(k*H). Recovery uses the first Hex8 element.",
              "cases": rows}
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    for row in rows:
        print(f"lambda={row['wavelength_m']:g} m h={row['h_m']:.6g} m "
              f"f={row['frequency_hz']:.7g} Hz error={100*row['frequency_relative_error']:+.4g}% "
              f"gradient/surface motion={row['surface_gradient_motion_ratio']:.6g} "
              f"static uZ error/unit bottom={row['unit_bottom_static_surface_recovered_uz_m']:.6g}")


if __name__ == "__main__":
    main()
