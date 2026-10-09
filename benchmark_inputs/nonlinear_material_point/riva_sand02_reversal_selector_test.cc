/* Host input/default/restart contracts for the three reversal rules. */
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

static int parse(const std::string &text,int is_02,int32_t *value,bool bias=false)
{
    FILE *f=tmpfile(); require(f!=nullptr,"tmpfile failed");
    require(fwrite(text.data(),1,text.size(),f)==text.size(),"fixture write failed");
    const int result=bias ? riva_read_bias_volume(f,"rivasand02_biasvolume",is_02,value) :
        riva_read_reversal_type(f,"rivasand02_reversaltype",is_02,value);
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
        require(parse("# omitted\nmaterial_model=RIVASand02\n",1,&value)==0 && value==3,
                "omitted selector did not default to type 3");
        for (int type : {1,2,3}) {
            require(parse("  rivasand02_reversaltype = "+std::to_string(type)+" # comment\n",1,&value)==0 && value==type,
                    "valid selector rejected");
            require(parse("rivasand02_reversaltype="+std::to_string(type),1,&value)==0 && value==type,
                    "selector without final newline rejected");
        }
        for (const char *bad : {"0","4","-1","1.0","1e0","2junk","1 2","nan",""}) {
            value=99;
            require(parse(std::string("rivasand02_reversaltype = ")+bad+"\n",1,&value)!=0 && value==99,
                    "invalid selector accepted or changed output");
        }
        for (const char *bad : {"rivasand02_reversaltype\n","rivasand02_reversaltype 2\n",
                               "rivasand02_reversaltype=1\nrivasand02_reversaltype=2\n"})
            require(parse(bad,1,&value)!=0,"malformed/duplicate selector accepted");
        require(parse("rivasand02_reversaltype=1\n",0,&value)!=0,"selector accepted for another model");
        require(parse("material_model=RIVASand\n",0,&value)==0 && value==3,
                "original model without selector rejected");
        for (const char *retired : {"reversalType=2\n", "  reversalType = 1 # old\n",
                                   "reversalType\n", "reversalType 2\n",
                                   "reversalType=2\nrivasand02_reversaltype=2\n",
                                   "rivasand02_reversaltype=2\nreversalType=2\n",
                                   "rivasand02_reversalType=2\n",
                                   "  rivasand02_reversalType = 1 # old\n",
                                   "rivasand02_reversalType\n", "rivasand02_reversalType 2\n",
                                   "rivasand02_reversalType=2\nrivasand02_reversaltype=2\n",
                                   "rivasand02_reversaltype=2\nrivasand02_reversalType=2\n"}) {
            value=99;
            require(parse(retired,1,&value)!=0 && value==99,
                    "retired key accepted or silently defaulted");
        }
        require(parse("# reversalType=2\n# rivasand02_reversalType=2\nrivasand02_reversaltype=2\n",1,&value)==0 && value==2,
                "comment mentioning retired key rejected");
        require(parse("reversalTypeExtra=2\n",1,&value)==0 && value==3,
                "unrelated key matched retired name");

        require(parse("# omitted\n",1,&value,true)==0 && value==1,
                "omitted BiasVolume did not default to mode 1");
        for (int mode : {0,1,2}) {
            require(parse("rivasand02_biasvolume="+std::to_string(mode),1,&value,true)==0 && value==mode,
                    "valid BiasVolume rejected");
        }
        for (const char *bad : {"-1","3","1.0","1e0","nan","1 2","","1junk"}) {
            value=99;
            require(parse(std::string("rivasand02_biasvolume=")+bad,1,&value,true)!=0 && value==99,
                    "invalid BiasVolume accepted or mutated output");
        }
        require(parse("rivasand02_biasvolume=1\nrivasand02_biasvolume=1\n",1,&value,true)!=0 &&
                parse("rivasand02_biasvolume=1\n",0,&value,true)!=0 &&
                parse("rivasand02_biasvolume 1\n",1,&value,true)!=0,
                "duplicate/model/malformed BiasVolume accepted");
        require(!riva_ib_bias_volume_valid(-1) && !riva_ib_bias_volume_valid(3),"invalid mode accepted");
        require(riva_ib_reference_parameters(1000.0).base.bias_reversible_volume_enabled==0 &&
                riva_ib_reference_parameters(1000.0).field_bias_mean_correction_enabled==0,
                "default parameters differ from mode 1");

        for (int mode : {0,1,2}) {
            const auto p=riva_ib_reference_parameters(1000.0,mode);
            require(p.base.bias_reversible_volume_enabled==(mode!=1) &&
                    p.field_bias_mean_correction_enabled==(mode==2),"incorrect mode flags");
            const auto m=riva_reference_material_parameters(&p.base);
            const tensor_t stress={-100000,-100000,-100000,0,0,0};
            const double e=riva_void_ratio_from_material_relative_density(&m,.9);
            for (int nsub : {1,4,16}) {
                riva_ib_state_t initial={};
                require(riva_ib_initialize_material(&p,&m,stress,e,&initial) &&
                    riva_ib_begin_dynamic_phase(&p,&m,&stress,&initial),"initialization failed");
                auto default_state=initial,one=initial,two=initial,three=initial;
                auto volume_probe=initial; volume_probe.phase_reversible_volume=0.001;
                if (mode==1) require(riva_ib_bias_volume_target(&p,&volume_probe,0)==0.0 &&
                    riva_ib_bias_volume_target(&p,&volume_probe,1)==0.001,
                    "mode 1 disabled separate phase volume");
                double previous=0;
                for (int i=1;i<=1024;++i) {
                    const double shear=.0015*std::sin(2*std::acos(-1.0)*i/512.0);
                    const tensor_t de={0,0,.0001/512,0,0,shear-previous}; previous=shear;
                    require(riva_ib_update_material(&p,&m,de,nsub,&default_state,nullptr,nullptr) &&
                            riva_ib_update_material(&p,&m,de,nsub,&one,nullptr,nullptr,1),"type 1/default failed");
                    require(riva_ib_update_material(&p,&m,de,nsub,&two,nullptr,nullptr,2),"type 2 failed");
                    riva_update_info_t info={};
                    const auto before=three.base.reversals;
                    require(riva_ib_update_material(&p,&m,de,nsub,&three,nullptr,&info,3),"type 3 failed");
                    require(same(default_state,three),"default differs from explicit type 3");
                    require(three.base.reversals-before>=0 && three.base.reversals-before<=1 &&
                            info.reversal_registered==three.base.reversals-before,"type 3 host event accounting");
                }
                require(one.base.reversals==0 && two.base.reversals==4 && three.base.reversals==4,
                        "selector did not distinguish the known mixed-path reversal counts");
                for (int invalid : {0,4,-1}) {
                    auto copy=two;
                    require(!riva_ib_update_material(&p,&m,{0,0,0,1e-5,0,0},nsub,
                            &copy,nullptr,nullptr,invalid) && same(copy,two),
                            "invalid direct selector mutated the committed state");
                }
                for (int custom : {0,1}) for (int type : {1,2,3}) {
                    const auto state=type==1?one:type==2?two:three;
                    const auto marker=riva_ib_checkpoint_revision(custom,type,mode);
                    require(marker==10u+3u*mode+type-1u,"wrong checkpoint marker");
                    require(riva_ib_checkpoint_compatible(marker,sizeof(state),RIVA_IB_PARAMETER_SHA256,custom,type,mode),
                            "same-type checkpoint rejected");
                    for (int other_mode : {0,1,2}) for (int other : {1,2,3})
                        if (other!=type || other_mode!=mode)
                            require(!riva_ib_checkpoint_compatible(marker,sizeof(state),RIVA_IB_PARAMETER_SHA256,custom,other,other_mode),
                                    "cross-control restart accepted");
                    for (unsigned old=0;old<=9;++old)
                        require(!riva_ib_checkpoint_compatible(old,sizeof(state),RIVA_IB_PARAMETER_SHA256,custom,type,mode),
                                "old-layout checkpoint accepted");
                    require(!riva_ib_checkpoint_compatible(marker,sizeof(state)-1,RIVA_IB_PARAMETER_SHA256,custom,type,mode) &&
                            !riva_ib_checkpoint_compatible(marker,sizeof(state),"wrong hash",custom,type,mode),
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
        }
        require(riva_ib_checkpoint_revision(0)==15 &&
                !riva_ib_checkpoint_compatible(0,sizeof(riva_ib_state_t),RIVA_IB_PARAMETER_SHA256,0,4),
                "default/invalid checkpoint contract failed");
        std::cout << "PASS: rivasand02_reversaltype input, type-3 default, BiasVolume modes and default, host events, rejected selectors and restart isolation\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n'; return 1;
    }
    return 0;
}
