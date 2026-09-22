/* Loading-branch research regression; no OpenSees or CUDA runtime required.
 * Usage: executable OUTPUT_DIR
 * The field/material-point entry point explicitly selects type 2 here.
 * Type 1 reproduces the diagnosed event loss on the mixed path.
 */
#define RIVA_STANDALONE
#include "../../quake/forward/nonlinear_riva_sand_intermediate_bias_research.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

static void require(bool ok, const char *message)
{ if (!ok) throw std::runtime_error(message); }

static bool same(const riva_ib_state_t &a, const riva_ib_state_t &b)
{
    double av[RIVA_IB_STATE_VALUE_COUNT], bv[RIVA_IB_STATE_VALUE_COUNT];
    riva_ib_state_values(&a,av); riva_ib_state_values(&b,bv);
    for (int i=0;i<RIVA_IB_STATE_VALUE_COUNT;++i)
        if (!std::isfinite(av[i]) || av[i]!=bv[i]) return false;
    return a.base.initialized==b.base.initialized &&
           a.base.geostatic_admitted==b.base.geostatic_admitted;
}

struct Snapshot {
    uint32_t revision;
    uint64_t state_bytes;
    char parameter_hash[65];
    double elastic_scale;
    riva_ib_state_t state;
};

static bool restore(const Snapshot &saved,int custom,
                    riva_material_parameters_t &m,riva_ib_state_t &s)
{
    if (!riva_ib_checkpoint_compatible(saved.revision,saved.state_bytes,
                                      saved.parameter_hash,custom,2) ||
        !std::isfinite(saved.elastic_scale) || saved.elastic_scale<=0.0)
        return false;
    m.E_ref=saved.elastic_scale;
    s=saved.state;
    return true;
}

static void checkpoint(riva_material_parameters_t &m,riva_ib_state_t &s,int custom)
{
    Snapshot saved={};
    saved.revision=riva_ib_checkpoint_revision(custom,2);
    saved.state_bytes=sizeof(s);
    std::memcpy(saved.parameter_hash,RIVA_IB_PARAMETER_SHA256,65);
    saved.elastic_scale=m.E_ref;
    saved.state=s;
    // A byte round trip retains the complete branch reference and cached G0.
    std::stringstream bytes(std::ios::in|std::ios::out|std::ios::binary);
    bytes.write(reinterpret_cast<const char*>(&saved),sizeof(saved));
    Snapshot read={};
    bytes.read(reinterpret_cast<char*>(&read),sizeof(read));
    auto restored=s; auto material=m;
    require(restore(read,custom,material,restored) && same(s,restored) &&
            material.E_ref==m.E_ref,"new checkpoint round trip failed");
    for (uint32_t old=0;old<=6;++old) {
        read.revision=old;
        require(!restore(read,custom,material,restored),"old checkpoint accepted");
        require(same(s,restored),"rejected checkpoint changed live state");
    }
    read=saved; read.revision=8;
    require(!restore(read,custom,material,restored),"future checkpoint accepted");
    read=saved; read.state_bytes--;
    require(!restore(read,custom,material,restored),"wrong state width accepted");
    read=saved; read.parameter_hash[0]='x';
    require(!restore(read,custom,material,restored),"wrong calibration accepted");
    require(restore(saved,custom,m,s),"restore failed");
}

struct Row { tensor_t stress; tensor_t reference; int64_t reversals; };
using History=std::vector<Row>;

static tensor_t rotated(tensor_t s,double angle)
{
    const double c=std::cos(angle),b=std::sin(angle);
    return {c*c*s.xx+b*b*s.yy-2*c*b*s.xy,
            b*b*s.xx+c*c*s.yy+2*c*b*s.xy,s.zz,
            c*b*(s.xx-s.yy)+(c*c-b*b)*s.xy,
            b*s.xz+c*s.yz,c*s.xz-b*s.yz};
}

static History run(int steps,int nsub,bool mixed,bool branch,
                   bool loose=false,double scale=1.0,double angle=0.0,
                   bool standard_row=false)
{
    auto p=riva_ib_reference_parameters(scale);
    auto m=riva_reference_material_parameters(&p.base);
    riva_material_set_G0(&p.base,&m,loose?450.0:900.0);
    m.h=loose?25.0:122.44207260468994;
    m.Q=loose?10.0:4.0; m.R=loose?1.5:4.0;
    if (standard_row) { m.Q=10.0; m.R=1.5; }
    const double e=riva_void_ratio_from_material_relative_density(&m,loose?0.36:0.90);
    const tensor_t initial={-100*scale,-100*scale,-100*scale,0,0,0};
    riva_ib_state_t s={},legacy={},restarted={};
    require(riva_ib_initialize_material(&p,&m,initial,e,&s) &&
            riva_ib_begin_dynamic_phase(&p,&m,&initial,&s),"initialization failed");
    legacy=restarted=s; auto restart_material=m;
    History history;
    history.push_back({s.base.stress,s.base.last_host_deviatoric_strain_direction,0});
    double previous=0.0;
    for (int i=1;i<=2*steps;++i) {
        const double t=double(i)/steps;
        const double gamma=0.003*std::sin(2*std::acos(-1.0)*t);
        const tensor_t deps=rotated({0,0,mixed?0.0001/steps:0,0,0,
                                     0.5*(gamma-previous)},angle);
        previous=gamma;
        const auto committed=s;
        riva_update_info_t info={};
        const int ok=branch ? riva_ib_update_material(&p,&m,deps,nsub,&s,nullptr,&info,2) :
            riva_ib_update_material_reference_ex(&p,&m,deps,nsub,&s,nullptr,&info,0);
        require(ok && info.accepted_substeps==nsub,"constitutive update failed");
        if (branch) {
            // Every accepted state must agree with a separately restored point.
            riva_update_info_t other={};
            require(riva_ib_update_material(&p,&restart_material,deps,nsub,
                     &restarted,nullptr,&other,2) && same(s,restarted),"restart trajectory changed");
            if (i==steps/4-1 || i==steps/4 || info.reversal_registered || i==steps)
                for (int custom : {0,1}) checkpoint(restart_material,restarted,custom);
            if (!mixed) {
                require(riva_ib_update_material_reference_ex(&p,&m,deps,nsub,
                          &legacy,nullptr,&other,0),"legacy pure-shear update failed");
                auto comparable=legacy;
                comparable.base.last_host_deviatoric_strain_direction=
                    s.base.last_host_deviatoric_strain_direction;
                require(same(comparable,s),"pure shear changed beyond reference slots");
            }
            if (i==steps/4) {
                auto trial=committed;
                require(riva_ib_update_material(&p,&m,riva_scale(deps,-1),nsub,
                          &trial,nullptr,&other,2),"speculative update failed");
                trial=committed;
                require(riva_ib_update_material(&p,&m,deps,nsub,&trial,nullptr,&other,2)
                        && same(trial,s),"discarded trial changed accepted state");
                auto failed=s;
                tensor_t bad=deps; bad.xx=std::numeric_limits<double>::infinity();
                require(!riva_ib_update_material(&p,&m,bad,nsub,&failed,nullptr,&other,2)
                        && same(failed,s),"failed trial changed committed state");
            }
        }
        double values[RIVA_IB_STATE_VALUE_COUNT]; riva_ib_state_values(&s,values);
        for (double value:values) require(std::isfinite(value),"nonfinite state");
        history.push_back({s.base.stress,s.base.last_host_deviatoric_strain_direction,
                           s.base.reversals});
    }
    return history;
}

static double stress_error(const History &a,const History &b,double scale=1.0,
                           double angle=0.0)
{
    require((b.size()-1)%(a.size()-1)==0,"incompatible time grids");
    const size_t stride=(b.size()-1)/(a.size()-1);
    double error=0.0;
    for (size_t i=0;i<a.size();++i) {
        const auto expected=riva_scale(rotated(a[i].stress,angle),scale);
        error=std::max(error,riva_norm(riva_sub(expected,b[i*stride].stress)));
    }
    return error;
}

int main(int argc,char **argv)
{
    try {
        require(argc==2,"usage: branch_reversal_test OUTPUT_DIR");
        const std::filesystem::path out(argv[1]); std::filesystem::create_directories(out);
        std::ofstream summary(out/"refinement.csv");
        summary << "steps_per_cycle,nsub,legacy_reversals,branch_reversals,error_to_8192_kPa\n";
        const auto fine=run(8192,4,true,true);
        // Same compact and explicit reference rows as the target CUDA fixtures.
        for (int steps : {128,512,2048})
            require(run(steps,4,true,true,false,1000.0,0.0,true).back().reversals==4,
                    "CUDA fixture row failed its host event-count check");
        for (int steps : {128,512,2048,8192}) {
            const auto old=run(steps,4,true,false);
            const auto corrected=run(steps,4,true,true);
            require(corrected.back().reversals==4,"corrected mixed path missed reversal");
            require(old.back().reversals==(steps==128?4:0),"legacy reproducer changed");
            require(run(steps,4,false,true).back().reversals==4,"pure shear count changed");
            const double error=stress_error(corrected,fine);
            summary << std::setprecision(17) << steps << ",4," << old.back().reversals
                    << ',' << corrected.back().reversals << ',' << error << '\n';
            if (steps==2048) require(error<0.5,"fine mixed-path stress refinement exceeds 0.5 kPa");
        }
        for (bool loose : {false,true}) {
            const auto shared_path=run(2048,4,true,true,loose);
            std::ofstream shared(out/(loose?"loose_reference.csv":"bottom_reference.csv"));
            shared << "time_s,sxx,syy,szz,sxy,syz,sxz,reversals\n" << std::setprecision(17);
            for (size_t i=1;i<shared_path.size();++i) {
                const auto &s=shared_path[i].stress;
                shared << double(i)/2048 << ',' << s.xx << ',' << s.yy << ',' << s.zz
                       << ',' << s.xy << ',' << s.yz << ',' << s.xz << ','
                       << shared_path[i].reversals << '\n';
            }
            const auto base=run(512,4,true,true,loose);
            const double rotation=stress_error(base,run(512,4,true,true,loose,1.0,0.61),1.0,0.61);
            const double units=stress_error(base,run(512,4,true,true,loose,1000.0),1000.0);
            require(rotation<1e-8,"rotation covariance failed");
            require(units<1e-5,"Pa/kPa covariance failed");
            std::cout << "layer=" << (loose?"loose":"bottom") << " rotation_kPa="
                      << rotation << " units_Pa=" << units << '\n';
            for (int nsub : {1,4,16})
                require(run(512,nsub,true,true,loose).back().reversals==4,
                        "substepping changed reversal count");
        }
        std::cout << "PASS: branch reversal refinement, pure-shear parity, rotation, units, "
                     "trial isolation, old-checkpoint rejection and restart continuation\n";
    } catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
