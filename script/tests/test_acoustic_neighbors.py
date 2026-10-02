"""Check production remote-leaf lookup against spatial bounds and materials.

The Horizon shifted-surface fixture covers the collective exchange and compares
one/four ranks. This host test catches aggregate fallback and Morton-cache gaps.
"""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest

from test_acoustic_plane_reconstruction import function


class AcousticNeighborTests(unittest.TestCase):
    def test_remote_leaf_bounds_material_and_discovery(self):
        root = Path(__file__).resolve().parents[2]
        solver = (root / 'quake/forward/psolve.cc').read_text()
        octor = (root / 'octor/octor.c').read_text()
        code = r'''
#include <map>
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#define PIXELLEVEL 4
#define REMOTE 0
#define AGGREGATE_SEARCH 1
#define LOG2_32b(v) ((v) ? 31 - __builtin_clz((unsigned)(v)) : -1)
using tick_t = int32_t;
struct point_t { tick_t x,y,z; };
struct vector3D_t { double x[3]; };
struct edata_t { double Vp,Vs,rho; };
struct node_t { tick_t x,y,z; };
struct elem_t { int lnid[8]; int8_t level; int64_t geid; void* data; };
struct octant_t { tick_t lx,ly,lz; int8_t level; int where; void* appdata; };
struct mesh_t { double ticksize; int lenum; node_t* nodeTable; elem_t* elemTable; };
struct cvmpayload_t {};
static int missing_remote = 0;
static octant_t aggregate = {0,0,0,0,REMOTE,NULL};
static point_t intervals[2] = {{0,0,0},{0,0,0}};
static octant_t* octor_searchoctant(void*,tick_t,tick_t,tick_t,int,int) { return missing_remote ? NULL : &aggregate; }
static const point_t* octor_getintervaltable(void*) { return intervals; }
static int32_t mesh_handoff_owner(const point_t*,int32_t,const point_t*) { return 1; }
static int material_values_are_elastic(double vp,double vs,double rho) { return vp>0 && vs>0 && rho>0; }
static int poroelastic_element_candidate(mesh_t*,int,int,void*,void*) { return 0; }
static double poroelastic_biot_storage_modulus_for_element(const edata_t*,int) { return 0; }
static struct { double theDomainX,theDomainY,theDomainZ; } Param = {16,16,16};
static struct { mesh_t* myMesh; void* myOctree; int theGroupSize, myID; } Global;
static void require(bool condition, const char* text) {
  if (!condition) { fprintf(stderr,"%s\n",text); exit(1); }
}
''' + 'int32_t\n' + function(octor, 'octor_zcompare(const void *ptr1, const void *ptr2)')
        start = solver.index('struct acoustic_neighbor_t {')
        end = solver.index('static void mesh_exchange_acoustic_neighbors(', start)
        code += solver[start:end]
        code += function(solver, 'static int mesh_probe_face_elastic_info(')
        code += r'''
int main() {
  mesh_t mesh = {1,0,NULL,NULL};
  Global.myMesh=&mesh; Global.myOctree=&aggregate; Global.theGroupSize=2;
  acoustic_neighbor_t a = {{0,0,0},2,10,1500,0,1025,0,0};
  acoustic_neighbor_t b = {{4,4,4},2,20,2400,750,2200,1e9,1};
  acoustic_neighbor_t c = {{8,8,8},1,30,2400,750,2200,0,0};
  g_acousticNeighbors[a.lower]=a;
  g_acousticNeighbors[b.lower]=b;
  g_acousticNeighbors[c.lower]=c;
  for (int x=0;x<16;++x) for (int y=0;y<16;++y) for (int z=0;z<16;++z) {
    int64_t expected = x<4 && y<4 && z<4 ? 10 :
      x>=4 && x<8 && y>=4 && y<8 && z>=4 && z<8 ? 20 :
      x>=8 && y>=8 && z>=8 ? 30 : -1;
    const auto* actual=mesh_cached_acoustic_neighbor({x,y,z});
    require((actual ? actual->geid : -1)==expected,"Morton lookup crossed a leaf or gap");
  }
  vector3D_t p={{5.5,6.5,7.5}};
  tick_t lo[3],hi[3]; int8_t level=-1; int local=-1; int32_t index=-1;
  require(mesh_find_octant_at_domain_point(&p,lo,hi,&level,&local,&index),"leaf lookup failed");
  require(level==2 && !local && index==-1,"remote aggregate used as leaf");
  for(int d=0;d<3;++d) require(lo[d]==4 && hi[d]==8,"wrong final leaf bounds");
  int32_t neighbor_level=-1; int64_t geid=-1; int poro=-1; double storage=0;
  require(mesh_probe_face_elastic_info(&p,NULL,NULL,NULL,&index,&neighbor_level,&geid,&poro,&storage,4),"corrected solid material lost");
  require(neighbor_level==2 && geid==20 && poro==1 && storage==1e9,"owner properties lost");
  p={{1.5,2.5,3.5}};
  require(!mesh_probe_face_elastic_info(&p,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,4),"water classified as solid");
  g_acousticCollectNeighbors=1;
  p={{6.5,.5,.5}};
  require(!mesh_find_octant_at_domain_point(&p,NULL,NULL,NULL,NULL,NULL),"missing leaf accepted");
  require(g_acousticNeighborRequests.size()==1,"missing leaf not queued");
  require(g_acousticNeighborRequests.begin()->second==1,"wrong request owner");
  missing_remote=1;
  p={{5.5,6.5,7.5}};
  require(mesh_find_octant_at_domain_point(&p,lo,hi,&level,&local,&index),"cached leaf requires a local octree aggregate");
  p={{7.5,.5,.5}};
  require(!mesh_find_octant_at_domain_point(&p,NULL,NULL,NULL,NULL,NULL),"missing remote leaf accepted");
  require(g_acousticNeighborRequests.size()==2,"missing remote octant was not queued");
  p={{16,0,0}};
  require(!mesh_find_octant_at_domain_point(&p,NULL,NULL,NULL,NULL,NULL),"out-of-domain lookup accepted");
  require(g_acousticNeighborRequests.size()==2,"out-of-domain request queued");
}
'''
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / 'neighbors.cc'
            binary = Path(tmp) / 'neighbors'
            source.write_text(code)
            compiler = shlex.split(os.environ.get('CXX', 'c++'))
            subprocess.run(compiler + ['-std=c++17', '-O2', str(source), '-o', str(binary)], check=True)
            subprocess.run([str(binary)], check=True)


if __name__ == '__main__':
    unittest.main()
