/* Host input/default/restart contracts for the two existing reversal rules. */
#define RIVA_STANDALONE
#include "../../quake/forward/nonlinear_riva_sand_intermediate_bias_research.h"
#include "../../quake/forward/nonlinear_riva_sand_input.h"
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>

static void require(bool ok,const char *why)
{ if (!ok) throw std::runtime_error(why); }

static int parse(const std::string &text,int is_02,int32_t *value)
{
    FILE *f=tmpfile(); require(f!=nullptr,"tmpfile failed");
    require(fwrite(text.data(),1,text.size(),f)==text.size(),"fixture write failed");
    const int result=riva_read_reversal_type(f,"reversalType",is_02,value);
    fclose(f); return result;
}

static bool same(const riva_ib_state_t &a,const riva_ib_state_t &b)
{
    double av[RIVA_IB_STATE_VALUE_COUNT],bv[RIVA_IB_STATE_VALUE_COUNT];
    riva_ib_state_values(&a,av); riva_ib_state_values(&b,bv);
    for (int i=0;i<RIVA_IB_STATE_VALUE_COUNT;++i)
        if (!std::isfinite(av[i]) || av[i]!=bv[i]) return false;
    return a.base.initialized==b.base.initialized &&
        a.base.geostatic_admitted==b.base.geostatic_admitted;
}

int main()
{
    try {
        int32_t value=99;
        require(parse("# omitted\nmaterial_model=RIVASand02\n",1,&value)==0 && value==1,
                "omitted selector did not default to type 1");
        for (int type : {1,2}) {
            require(parse("  reversalType = "+std::to_string(type)+" # comment\n",1,&value)==0 && value==type,
                    "valid selector rejected");
            require(parse("reversalType="+std::to_string(type),1,&value)==0 && value==type,
                    "selector without final newline rejected");
        }
        for (const char *bad : {"0","3","-1","1.0","1e0","2junk","1 2","nan",""}) {
            value=99;
            require(parse(std::string("reversalType = ")+bad+"\n",1,&value)!=0 && value==99,
                    "invalid selector accepted or changed output");
        }
        for (const char *bad : {"reversalType\n","reversalType 2\n",
                               "reversalType=1\nreversalType=2\n"})
            require(parse(bad,1,&value)!=0,"malformed/duplicate selector accepted");
        require(parse("reversalType=1\n",0,&value)!=0,"selector accepted for another model");
        require(parse("material_model=RIVASand\n",0,&value)==0 && value==1,
                "original model without selector rejected");

        const auto p=riva_ib_reference_parameters(1000.0);
        const auto m=riva_reference_material_parameters(&p.base);
        const tensor_t stress={-100000,-100000,-100000,0,0,0};
        const double e=riva_void_ratio_from_material_relative_density(&m,.9);
        for (int nsub : {1,4,16}) {
            riva_ib_state_t initial={};
            require(riva_ib_initialize_material(&p,&m,stress,e,&initial) &&
                riva_ib_begin_dynamic_phase(&p,&m,&stress,&initial),"initialization failed");
            auto default_state=initial,explicit_one=initial,two=initial;
            double previous=0;
            for (int i=1;i<=1024;++i) {
                const double shear=.0015*std::sin(2*std::acos(-1.0)*i/512.0);
                const tensor_t de={0,0,.0001/512,0,0,shear-previous}; previous=shear;
                require(riva_ib_update_material(&p,&m,de,nsub,&default_state,nullptr,nullptr) &&
                        riva_ib_update_material(&p,&m,de,nsub,&explicit_one,nullptr,nullptr,1) &&
                        same(default_state,explicit_one),"default differs from explicit type 1");
                require(riva_ib_update_material(&p,&m,de,nsub,&two,nullptr,nullptr,2),"type 2 failed");
            }
            require(default_state.base.reversals==0 && two.base.reversals==4,
                    "selector did not distinguish the known mixed-path reversal counts");
            for (int invalid : {0,3,-1}) {
                auto copy=two;
                require(!riva_ib_update_material(&p,&m,{0,0,0,1e-5,0,0},nsub,
                        &copy,nullptr,nullptr,invalid) && same(copy,two),
                        "invalid direct selector mutated the committed state");
            }
            for (int custom : {0,1}) for (int type : {1,2}) {
                const auto state=type==1?default_state:two;
                const auto marker=riva_ib_checkpoint_revision(custom,type);
                require(marker==(type==1?8u:7u),"wrong checkpoint marker");
                require(riva_ib_checkpoint_compatible(marker,sizeof(state),RIVA_IB_PARAMETER_SHA256,custom,type),
                        "same-type checkpoint rejected");
                require(!riva_ib_checkpoint_compatible(marker,sizeof(state),RIVA_IB_PARAMETER_SHA256,custom,3-type),
                        "cross-type restart accepted");
                for (unsigned old=0;old<=6;++old)
                    require(!riva_ib_checkpoint_compatible(old,sizeof(state),RIVA_IB_PARAMETER_SHA256,custom,type),
                            "pre-selector checkpoint accepted");
                require(!riva_ib_checkpoint_compatible(marker,sizeof(state)-1,RIVA_IB_PARAMETER_SHA256,custom,type) &&
                        !riva_ib_checkpoint_compatible(marker,sizeof(state),"wrong hash",custom,type),
                        "corrupt checkpoint accepted");
                auto restored=initial;
                std::memcpy(&restored,&state,sizeof(state));
                auto continued=state;
                const tensor_t de={0,0,2e-7,0,0,-1e-5};
                require(riva_ib_update_material(&p,&m,de,nsub,&restored,nullptr,nullptr,type) &&
                        riva_ib_update_material(&p,&m,de,nsub,&continued,nullptr,nullptr,type) &&
                        same(restored,continued),"restart changed subsequent response");
            }
        }
        require(riva_ib_checkpoint_revision(0)==8 &&
                !riva_ib_checkpoint_compatible(0,sizeof(riva_ib_state_t),RIVA_IB_PARAMETER_SHA256,0,3),
                "default/invalid checkpoint contract failed");
        std::cout << "PASS: reversalType input, type-1 default, type-2 event counts, rejected selectors and restart isolation\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n'; return 1;
    }
    return 0;
}
