"""Exercise production plane/station reconstruction without CUDA or MPI.

A compact velocity-potential pulse has an analytic final displacement and
zero final velocity. Pressure mode also retains its constant-acceleration
solution. Both tests compare the actual plane and station implementations.
"""

import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest


def function(source, signature):
    start = source.index(signature)
    return source[start:source.index("\n}\n", start) + 3]


class AcousticPlaneTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(__file__).resolve().parents[2]
        planes = (root / "quake/forward/io_planes.cc").read_text()
        solver = (root / "quake/forward/psolve.cc").read_text()
        start = planes.index("static int plane_pressure_gradient(")
        end = planes.index("static void plane_interpolate_point(", start)
        station = function(solver, "static int station_pressure_gradient(int32_t iStation, double grad[3]) {")
        station += function(solver, "static void solver_update_acoustic_station_kinematics(int32_t step) {")
        source = r"""
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "acoustic_output.h"
struct vector3D_t { double x[3]; };
struct mysolver_t { double* pressure_tm1; double* pressure_tm2; };
struct plane_strip_element_t {
    int pressure_is_defined, acoustic_kinematics_ready;
    int nodestointerpolate[8]; vector3D_t localcoords; double h, rho;
    double acoustic_disp[3], acoustic_vel[3], acoustic_acc[3], acoustic_prev_acc[3];
    acoustic_output_cpml_t acoustic_cpml;
};
typedef plane_strip_element_t station_t;
static struct { mysolver_t* mySolver; } Global;
static struct {
    station_t* myStations; int myNumberOfStations, theAcousticVelocityPotential;
    double theDeltaT, theDeltaTSquared;
} Param;
static int g_acousticPressureStationDataReady = 0;
static int acoustic_pressure_enabled() { return 1; }
static const double* solver_acoustic_station_result(int) { return NULL; }
""" + planes[start:end] + station + r"""
static void close_value(double actual, double expected, const char* label) {
    if (!isfinite(actual) || fabs(actual - expected) > 2.e-10) {
        fprintf(stderr, "%s: actual=%.17g expected=%.17g\n", label, actual, expected);
        exit(1);
    }
}
int main(int argc, char** argv) {
    const int cpml = argc > 1 && strcmp(argv[1], "cpml") == 0;
    const int potential = cpml || (argc > 1 && strcmp(argv[1], "potential") == 0);
    const int missing = argc > 1 && strcmp(argv[1], "missing") == 0;
    double field0[8] = {}, field1[8] = {};
    mysolver_t solver = {field0, field1}; Global.mySolver = &solver;
    plane_strip_element_t plane = {};
    plane.h = 10; plane.rho = 1000; plane.pressure_is_defined = 1;
    plane.localcoords = {{0.31, -0.2, 0.13}};
    for (int k=0; k<8; ++k) plane.nodestointerpolate[k] = k;
    const double dt = .001, period = 2;
    if (cpml) {
        // Three different stretches: CFS, unshifted PML, and an identity
        // tangential axis. The integrated pulse gain is H(0)=alpha/(k*alpha+d).
        const double d[3] = {4, 3, 0}, k[3] = {2, 1, 1}, alpha[3] = {1, 0, 0};
        plane.acoustic_cpml.active = 1;
        for (int c=0; c<3; ++c) {
            const double lambda = d[c]/k[c] + alpha[c];
            plane.acoustic_cpml.inv_kappa[c] = 1/k[c];
            plane.acoustic_cpml.b[c] = exp(-lambda*dt);
            plane.acoustic_cpml.a[c] = lambda > 0
                ? d[c]/(k[c]*k[c]*lambda)*expm1(-lambda*dt) : 0;
        }
    }
    station_t station = plane; Param.myStations = &station; Param.myNumberOfStations = 1;
    const double amplitude[3] = {.2, -.1, .05};
    Param.theDeltaT = dt; Param.theDeltaTSquared = dt*dt;
    Param.theAcousticVelocityPotential = potential;
    if (missing) {
        solver.pressure_tm2 = NULL;
        plane_update_acoustic_point(&solver, &plane, 0, dt, dt*dt, 1/dt);
        if (plane.acoustic_kinematics_ready) return 1;
        return 0;
    }
    if (!potential) solver.pressure_tm2 = NULL; // pressure mode needs only p[n]
    const int steps = cpml ? 40000 : (potential ? 4000 : 1000);
    for (int n=0; n<=steps; ++n) {
        if (potential) {
            // Reproduce solver_swap_pressure_buffers + the current-field-only
            // resident GPU gather. The old patch is now the previous field.
            double* previous = solver.pressure_tm1;
            solver.pressure_tm1 = solver.pressure_tm2;
            solver.pressure_tm2 = previous;
        }
        const double t = n*dt;
        const double pulse = t < period ? pow(sin(M_PI*t/period), 2) : 0;
        for (int k=0; k<8; ++k) {
            const double x = (k&1) ? plane.h : 0;
            const double y = (k&2) ? plane.h : 0;
            const double z = (k&4) ? plane.h : 0;
            const double linear = plane.rho*(amplitude[0]*x + amplitude[1]*y + amplitude[2]*z);
            solver.pressure_tm1[k] = potential ? (linear + 1000)*pulse : -linear;
        }
        plane_update_acoustic_point(&solver, &plane, n, dt, dt*dt, potential ? 1/dt : 0);
        solver_update_acoustic_station_kinematics(n);
        for (int c=0; c<3; ++c) {
            close_value(plane.acoustic_disp[c], station.acoustic_disp[c], "station/plane displacement");
            close_value(plane.acoustic_vel[c], station.acoustic_vel[c], "station/plane velocity");
            close_value(plane.acoustic_acc[c], station.acoustic_acc[c], "station/plane acceleration");
            if (!potential) {
                close_value(plane.acoustic_disp[c], .5*amplitude[c]*t*t, "pressure displacement");
                close_value(plane.acoustic_vel[c], amplitude[c]*t, "pressure velocity");
            } else if (!cpml && t >= period + dt) {
                close_value(plane.acoustic_disp[c], amplitude[c]*period/2, "pulse final displacement");
                close_value(plane.acoustic_vel[c], 0, "pulse residual velocity");
            }
        }
    }
    if (cpml) {
        const double gain[3] = {1.0/6, 0, 1};
        for (int c=0; c<3; ++c) {
            close_value(plane.acoustic_disp[c], amplitude[c]*period/2*gain[c], "CPML integrated pulse gain");
            close_value(plane.acoustic_vel[c], 0, "CPML pulse residual velocity");
        }
    }
    printf("PASS %s: XYZ displacement %.12g %.12g %.12g m\n",
           potential ? "velocity_potential" : "pressure",
           plane.acoustic_disp[0], plane.acoustic_disp[1], plane.acoustic_disp[2]);
}
"""
        temp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(temp.cleanup)
        cpp = Path(temp.name) / "reconstruction.cc"
        cls.binary = Path(temp.name) / "reconstruction"
        cpp.write_text(source)
        subprocess.run(shlex.split(os.environ.get("CXX", "c++")) +
                       ["-std=c++11", "-O2", "-Wall", "-Wextra", "-Werror",
                        "-I", str(root / "quake/forward"),
                        str(cpp), "-o", str(cls.binary)], check=True)

    def run_case(self, mode):
        result = subprocess.run([str(self.binary), mode], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_velocity_potential_pulse_has_no_spurious_residual_velocity(self):
        self.run_case("potential")

    def test_pressure_constant_acceleration_and_null_previous_buffer(self):
        self.run_case("pressure")

    def test_missing_potential_history_does_not_fall_back_to_raw_field(self):
        self.run_case("missing")

    def test_cpml_pulse_has_analytic_integrated_gain_and_no_residual_velocity(self):
        self.run_case("cpml")


if __name__ == "__main__":
    unittest.main()
