"""Compile the production setup exporter with a tiny host-only mesh fixture."""

import csv
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest


class PlaneMetadataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(__file__).resolve().parents[2]
        source = (root / "quake/forward/io_planes.cc").read_text()
        start = source.index("static int32_t plane_element_is_acoustic(")
        end = source.index("static void plane_init_acoustic_metadata(", start)
        # The production material classifier and exporter run unchanged.
        harness = r"""
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#define ACOUSTIC_VS_EPS 1.0e-6
struct vector3D_t { double x[3]; };
struct edata_t { float Vp, Vs, rho, edgesize; };
struct octant_t { int id; };
struct plane_strip_element_t {
  int point_is_valid, topostrip, owner_is_acoustic, pressure_is_defined;
};
struct plane_t {
  vector3D_t origincoords;
  int numberofstepsdowndip, numberofstripsthisplane;
  double stepalongstrike, stepdowndip, dip, strike;
  int stripstart[1], stripend[1];
  plane_strip_element_t* strip[1];
};
static plane_t* thePlanes;
static int g_reconstructAcousticPlaneDisplacements = 1;
static void solver_abort(const char*, const char*, const char* format, ...) {
  va_list args; va_start(args, format); vfprintf(stderr, format, args); va_end(args);
  exit(86);
}
static FILE* hu_fopen(const char* path, const char* mode) {
  FILE* fp = fopen(path, mode); if (!fp) exit(87); return fp;
}
static vector3D_t compute_global_coords(vector3D_t origin, vector3D_t local,
                                        double, double, double) {
  for (int c=0; c<3; ++c) local.x[c] += origin.x[c]; return local;
}
static double get_thebase_topo() { return 100; }
static double point_elevation(double x, double y) { return 100 - 2*x - y; }
static int search_point(vector3D_t xyz, octant_t** octant) {
  static octant_t found;
  found.id = (int)(xyz.x[0]/10)*3 + (int)(xyz.x[1]/20); *octant = &found; return 1;
}
static int32_t solver_get_octant_material(octant_t* octant, edata_t* material) {
  material->Vp = octant->id == 5 ? -1 : 1500;
  material->Vs = octant->id < 3 ? 0 : 500;
  material->rho = 1000; material->edgesize = 25; return octant->id;
}
static int cpml_pml_element(int index) { return index % 3 == 0; }
""" + source[start:end] + r"""
int main(int, char** argv) {
  plane_t plane = {};
  plane_strip_element_t points[6] = {};
  plane.origincoords.x[2] = 5;
  plane.numberofstepsdowndip = 3; plane.numberofstripsthisplane = 1;
  plane.stepalongstrike = 10; plane.stepdowndip = 20;
  plane.stripstart[0] = 0; plane.stripend[0] = 5; plane.strip[0] = points;
  for (int k=0; k<6; ++k) {
    points[k].point_is_valid = k != 4;
    points[k].owner_is_acoustic = k < 3;
    points[k].pressure_is_defined = k == 1 || k == 2;
  }
  points[0].topostrip = 1;
  thePlanes = &plane;
  planes_write_point_metadata(0, 1, argv[1]);
}
"""
        temp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(temp.cleanup)
        cpp = Path(temp.name) / "metadata.cc"
        cls.binary = Path(temp.name) / "metadata"
        cpp.write_text(harness)
        subprocess.run(shlex.split(os.environ.get("CXX", "c++")) +
                       ["-std=c++11", "-Wall", "-Wextra", "-Werror", str(cpp),
                        "-o", str(cls.binary)], check=True)
        # CPML membership must be queried after operator setup, not plane setup.
        solver = (root / "quake/forward/psolve.cc").read_text()
        assert solver.index("  solver_init();") < solver.index("  planes_write_point_metadata(")

    def run_export(self, mode):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        path = Path(temp.name)
        env = os.environ.copy()
        env.pop("HERCULES_PLANE_METADATA", None)
        if mode is not None:
            env["HERCULES_PLANE_METADATA"] = mode
        result = subprocess.run([str(self.binary), str(path)], env=env,
                                capture_output=True, text=True, check=False)
        return result, path

    def test_disabled_by_default_and_rejects_invalid_values(self):
        for mode in (None, "", "0"):
            result, path = self.run_export(mode)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(list(path.iterdir()), [])
        result, _ = self.run_export("auto")
        self.assertEqual(result.returncode, 86)

    def test_containing_material_cpml_coordinates_and_reconstruction_flags(self):
        result, path = self.run_export("1")
        self.assertEqual(result.returncode, 0, result.stderr)
        with (path / "plane_point_metadata.0.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 6)
        self.assertEqual([r["medium"] for r in rows],
                         ["acoustic", "acoustic", "acoustic", "solid", "unknown", "inactive"])
        self.assertEqual([r["cpml"] for r in rows], ["1", "0", "0", "1", "-1", "0"])
        self.assertEqual((rows[5]["x_m"], rows[5]["y_m"], rows[5]["z_m"]), ("10", "40", "45"))
        self.assertEqual((rows[0]["acoustic_reconstruction"], rows[0]["pressure_defined"]), ("1", "0"))
        self.assertEqual((rows[1]["rho_kg_m3"], rows[1]["h_m"]), ("1000", "25"))


if __name__ == "__main__":
    unittest.main()
