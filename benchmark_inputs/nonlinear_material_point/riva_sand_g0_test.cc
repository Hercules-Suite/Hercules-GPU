#define RIVA_STANDALONE
#include "../../quake/forward/nonlinear_riva_sand_intermediate_bias_research.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

static void require(bool ok, const char *what) {
    if (!ok) throw std::runtime_error(what);
}
static void close(double a, double b) {
    require(std::isfinite(a) && std::abs(a-b)<=1.e-12*(1+std::abs(b)), "stiffness/units mismatch");
}
int main() {
    try {
        for (double units : {1.0,1000.0}) {
            auto p=riva_ib_reference_parameters(units);
            const auto reference=riva_reference_material_parameters(&p.base);
            for (double g0 : {350.0,RIVA_REFERENCE_G0,700.0}) {
                auto material=reference;
                riva_material_set_G0(&p.base,&material,g0);
                require(riva_material_parameters_valid(&p.base,&material),"valid G0 rejected");
                if (g0==RIVA_REFERENCE_G0)
                    require(material.E_ref==reference.E_ref,"reference fingerprint drift");
                for (double pressure : {0.001,2.0,26.266666666666666,101.3}) {
                    double G,K;
                    riva_moduli(&p.base,&material,pressure*units,&G,&K);
                    close(G/units,g0*101.3*std::pow(pressure/101.3,0.65));
                    close(K/G,2*(1+p.base.nu)/(3*(1-2*p.base.nu)));
                }
                // Independent oracle: change the existing E_ref parameter,
                // then compare both update backbones with the per-material cache.
                auto oracle=p;
                oracle.base.E_ref=material.E_ref;
                const auto om=riva_reference_material_parameters(&oracle.base);
                tensor_t initial=riva_scale({-19.4,-19.4,-40,10,0,0},units);
                riva_state_t a={},b={};
                riva_ib_state_t c={},d={};
                require(riva_initialize_material(&p.base,&material,initial,.601,&a) &&
                        riva_initialize_material(&oracle.base,&om,initial,.601,&b),"original init");
                require(riva_ib_initialize_material(&p,&material,initial,.601,&c) &&
                        riva_ib_initialize_material(&oracle,&om,initial,.601,&d),"02 init");
                for (int step=0;step<128;++step) {
                    const double angle=2*3.14159265358979323846/32;
                    tensor_t de={0,0,0,0.0002*(std::sin((step+1)*angle)-std::sin(step*angle)),0,0};
                    require(riva_update_material(&p.base,&material,de,4,&a,nullptr,nullptr) &&
                            riva_update_material(&oracle.base,&om,de,4,&b,nullptr,nullptr),"original update");
                    require(riva_ib_update_material(&p,&material,de,4,&c,nullptr,nullptr) &&
                            riva_ib_update_material(&oracle,&om,de,4,&d,nullptr,nullptr),"02 update");
                    double av[RIVA_IB_STATE_VALUE_COUNT],bv[RIVA_IB_STATE_VALUE_COUNT];
                    riva_state_values(&a,av); riva_state_values(&b,bv);
                    for(int i=0;i<RIVA_STATE_VALUE_COUNT;++i)require(av[i]==bv[i],"original G0/E_ref mismatch");
                    riva_ib_state_values(&c,av); riva_ib_state_values(&d,bv);
                    for(int i=0;i<RIVA_IB_STATE_VALUE_COUNT;++i)require(av[i]==bv[i],"02 G0/E_ref mismatch");
                }
            }
            for(double bad : {0.0,-1.0,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
                auto m=reference; riva_material_set_G0(&p.base,&m,bad);
                require(!riva_material_parameters_valid(&p.base,&m),"invalid G0 accepted");
            }
        }
        std::cout<<"PASS: G0 reference identity, custom E_ref equivalence for both generations, pressure/Pa/kPa scaling and invalid input\n";
    } catch(const std::exception &e) { std::cerr<<e.what()<<'\n'; return 1; }
}
