"""Compile production surface selection and mass assembly without CUDA/MPI.

Translation must preserve the weak gravity operator. Buried water, the old
octree top in a shifted model, and a submerged face must not receive it.
"""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest

from test_acoustic_plane_reconstruction import function


class AcousticSurfaceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(__file__).resolve().parents[2]
        source = (root / 'quake/forward/psolve.cc').read_text()
        geometry = function(source, 'static double acoustic_free_surface_z() {')
        geometry += function(source, 'static int acoustic_element_has_free_surface(')
        start = source.index('      const double bulk = edata->rho * edata->Vp * edata->Vp;',
                             source.index('static void solver_init() {'))
        mass = source[start:source.index('      ep->c1 =', start)]
        start = source.index('      for (j = 0; j < 8; j++) {', start)
        assembly = source[start:source.index('      localPressureElements++;', start)]
        cpp = r'''
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define YES 1
#define ACOUSTIC_GRAVITY_ACCELERATION_M_S2 9.81
typedef double solver_float;
struct edata_t { double rho, Vp, edgesize; };
struct elem_t { int32_t lnid[8]; void* data; };
struct node_t { int x, y, z; };
struct mesh_t { node_t nodeTable[8]; };
struct tree_t { int nearendp[3]; double ticksize; };
struct solver_t { double pressure_mass[8], pressure_damp[8]; };
static struct { int includeTopography, includeBuildings, theAcousticFreeSurfaceGravity;
                double theDeltaT; } Param;
static struct { mesh_t* myMesh; tree_t* myOctree; solver_t* mySolver; } Global;
static double base, shift, terrain;
static double get_thebase_topo() { return base; }
static double get_surface_shift() { return shift; }
static double point_elevation(double, double) { return terrain; }
''' + geometry + r'''
static void require(int ok, const char* text) {
  if (!ok) { fprintf(stderr, "%s\n", text); exit(1); }
}
static void run(int topo, int buildings, double datum, double top,
                double terrainZ, int gravity, int expected) {
  mesh_t mesh = {}; tree_t tree = {{0,0,0}, 0.0625}; solver_t solver = {};
  Global.myMesh = &mesh; Global.myOctree = &tree; Global.mySolver = &solver;
  Param.includeTopography=topo; Param.includeBuildings=buildings;
  Param.theAcousticFreeSurfaceGravity=gravity; Param.theDeltaT=.0005;
  base = topo ? datum : 0; shift = buildings ? datum : 0; terrain=terrainZ;
  edata_t material = {1025, 1500, 10}; edata_t* edata=&material;
  elem_t element = {{0,1,2,3,4,5,6,7}, &material}; elem_t* elemp=&element;
  mesh.nodeTable[0].z = (int)llround(top/tree.ticksize);
  double pressureFreeSurfaceHost[8] = {};
  int localPressureSurfaceFaces = 0, j;
  const double acousticSurfaceZ = acoustic_free_surface_z();
  require(acousticSurfaceZ == datum, "wrong model surface datum");
''' + mass + r'''
  // Match the production tangential surface stretch, deliberately distinct
  // from the three-axis volume stretch to catch accidental use of the latter.
  double acoustic_kprod=24, acoustic_damp=9;
  double gravity_surface_kprod=6, gravity_surface_damp=5;
''' + assembly + r'''
  require(localPressureSurfaceFaces == expected, "wrong exposed surface-face count");
  for (int n=0; n<8; ++n) {
    const int surface = expected && n<4;
    const double added = surface && gravity ? 100/(4*1025*9.81) : 0;
    const double massExpected = 1000.0/(8.0*1025*1500*1500)*24 + added*6;
    const double dampExpected = 1000.0/(8.0*1025*1500*1500)*9*.0005 + added*5*.0005;
    require(pressureFreeSurfaceHost[n] == surface, "wrong pressure-release node flag");
    require(fabs(solver.pressure_mass[n] - massExpected) < 1e-15, "wrong gravity mass");
    require(fabs(solver.pressure_damp[n] - dampExpected) < 1e-15, "wrong gravity CPML damping");
  }
}
int main() {
  for (int gravity=0; gravity<=1; ++gravity) {
    run(0,0,0,0,0,gravity,1);
    run(1,0,2812.5,2812.5,2812.5,gravity,1);
    run(1,0,2812.5,0,2812.5,gravity,0);
    run(1,0,2812.5,2822.5,2812.5,gravity,0);
    run(1,0,2812.5,2812.5,2802.5,gravity,0);
    run(0,1,500,500,500,gravity,1);
  }
}
'''
        cls.temp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temp.cleanup)
        path = Path(cls.temp.name)
        (path / 'surface.cc').write_text(cpp)
        cls.binary = path / 'surface'
        subprocess.run(shlex.split(os.environ.get('CXX', 'c++')) +
                       ['-std=c++11', '-O2', '-Wall', '-Wextra', '-Werror',
                        str(path / 'surface.cc'), '-o', str(cls.binary)], check=True)

    def test_translated_pressure_release_and_gravity_surface_assembly(self):
        subprocess.run([str(self.binary)], check=True)


if __name__ == '__main__':
    unittest.main()
