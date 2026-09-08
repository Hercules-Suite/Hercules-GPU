#define RIVA_STANDALONE
#include "../../quake/forward/nonlinear_riva_sand_intermediate_bias_research.h"
#include "../../quake/forward/nonlinear_riva_sand_input.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <cstring>

static void require(bool ok,const char *why) { if(!ok) throw std::runtime_error(why); }
static void close(double a,double b) {
    require(std::isfinite(a) && std::abs(a-b)<2e-12*(1+std::abs(b)),"modulus mismatch");
}
static int parse(const std::string &body,int rows,int custom,int vs,double *out) {
    FILE *f=tmpfile(); require(f!=NULL,"tmpfile");
    const std::string input="# fixture\nmaterial_model = RIVASand\nmaterial_properties_list =\n"+body;
    fwrite(input.data(),1,input.size(),f);
    const auto p=riva_reference_parameters(1000);
    const int rc=riva_read_material_rows(f,"material_properties_list",rows,custom,vs,&p,out);
    fclose(f); return rc;
}
static void input_tests() {
    double v[26]={};
    require(parse("150 .5 # loose\n# comment\n450 .9\nnext_key = 1\n",2,0,1,v)==0,"Dr-only table");
    require(v[0]==150 && v[1]==.5 && v[2]==0 && v[3]==1.25 && v[14]==.9,"defaults/order");
    require(parse("150 .5 350\n450 .9 700\n",2,0,1,v)==0 && v[15]==700,"Dr G0 table");
    require(parse(".6\n",1,0,0,v)==0 && v[1]==0,"compact auto material point");
    require(parse(".6 700\n",1,0,0,v)==0 && v[1]==700,"compact explicit material point");
    const std::string full="150 .6 700 1.25 1.125 122.44207260468994 .945 .025 .819 .51 10 1.5 .65\n";
    require(parse(full,1,1,1,v)==0 && v[8]==.819 && v[12]==.65,"Custom order");
    require(parse("150 .6 0 1.25 1.125 122.4 .945 .025 .78 .51 10 1.5 .65\n",1,1,1,v)==0,"Custom auto");
    require(parse(full,1,0,1,v)!=0,"full row accepted by standard");
    require(parse("150 .6\n",1,1,1,v)!=0,"compact row accepted by Custom");
    require(parse("150 .6.700\n",1,0,1,v)!=0,"missing numeric delimiter accepted");
    for(const char *bad : {"150 .6 -1\n","150 .6 nan\n","150 .6 inf\n","150 2\n",
        "150 .6 350 extra\n","150 .6\n450 .9 700\n","150 .6 0\n450 .9 700\n",
        "150 .6\n450 .9\n600 .9\n"}) {
        require(parse(bad,2,0,1,v)!=0,"invalid table accepted");
    }
}

static void calibration_tests() {
    for(double units : {1.0,1000.0}) for(double pv : {1.6,40.0,100.0})
    for(double Dr : {.5,.663,.9}) for(double bias : {0.0,.0875,.3}) {
        auto p=riva_ib_reference_parameters(units);
        auto m=riva_reference_material_parameters(&p.base);
        m.e_max=.819; // exercise material-specific density in 02 as well
        const double e=riva_void_ratio_from_material_relative_density(&m,Dr);
        const auto original=m;
        const tensor_t stress=riva_scale({-.6*pv,-.6*pv,-pv,bias*pv,0,0},units);
        tensor_t reference=stress; reference.xy=0;
        riva_state_t a={}; riva_ib_state_t b={};
        require(riva_initialize_geostatic_material(&p.base,&m,stress,e,reference,&a,NULL),"original admission");
        require(riva_ib_initialize_material(&p,&m,stress,e,&b) &&
            riva_ib_admit_geostatic_state(&p,&m,&b,NULL) &&
            riva_ib_begin_dynamic_phase(&p,&m,&reference,&b),"02 admission");
        const auto before_a=a; const auto before_b=b;
        auto ma=m,mb=m;
        const double target=1800*200.0*200.0/1000*units;
        require(riva_fit_initial_shear(&p.base,&ma,&a,target),"original fit");
        require(riva_ib_fit_initial_shear(&p,&mb,&b,target),"02 fit");
        require(!memcmp(&a,&before_a,sizeof(a)) && !memcmp(&b,&before_b,sizeof(b)),"fit changed history");
        require(m.E_ref==original.E_ref,"fit changed shared material");
        double G,K;
        riva_moduli_for_state(&p.base,&ma,riva_cone_pressure(&p.base,a.stress),&a,&G,&K); close(G,target);
        riva_ib_moduli_for_state(&p,&mb,&b,riva_cone_pressure(&p.base,b.base.stress),&G,&K); close(G,target);
        // A checkpoint must carry E_ref, not reconstruct it from the now-changed p'.
        for(int nsub : {1,4,16}) {
            auto c=b,d=b; double cache=mb.E_ref;
            auto cm=mb,dm=mb;
            auto ac=a,ad=a; auto am=ma,an=ma;
            for(int step=0;step<24;++step) {
                tensor_t de={0,0,0,(step<6 || step>=18?1:-1)*1e-6,0,0};
                require(riva_ib_update_material(&p,&cm,de,nsub,&c,NULL,NULL) &&
                    riva_ib_update_material(&p,&dm,de,nsub,&d,NULL,NULL),"02 cached update");
                require(riva_update_material(&p.base,&am,de,nsub,&ac,NULL,NULL) &&
                    riva_update_material(&p.base,&an,de,nsub,&ad,NULL,NULL),"original cached update");
                if(step==11) {
                    FILE *f=tmpfile(); require(f!=NULL,"restart file");
                    fwrite(&d,sizeof(d),1,f); fwrite(&cache,sizeof(cache),1,f); rewind(f);
                    d={}; dm=m;
                    require(fread(&d,sizeof(d),1,f)==1 && fread(&dm.E_ref,sizeof(double),1,f)==1,"restart read");
                    fclose(f);
                    f=tmpfile(); require(f!=NULL,"original restart file");
                    fwrite(&ad,sizeof(ad),1,f); fwrite(&an.E_ref,sizeof(double),1,f); rewind(f);
                    ad={}; an=m;
                    require(fread(&ad,sizeof(ad),1,f)==1 && fread(&an.E_ref,sizeof(double),1,f)==1,"original restart read");
                    fclose(f);
                }
                double cv[RIVA_IB_STATE_VALUE_COUNT],dv[RIVA_IB_STATE_VALUE_COUNT];
                riva_ib_state_values(&c,cv); riva_ib_state_values(&d,dv);
                for(int j=0;j<RIVA_IB_STATE_VALUE_COUNT;++j) require(cv[j]==dv[j],"restart drift");
                riva_state_values(&ac,cv); riva_state_values(&ad,dv);
                for(int j=0;j<RIVA_STATE_VALUE_COUNT;++j) require(cv[j]==dv[j],"original restart drift");
                require(cm.E_ref==cache && dm.E_ref==cache,"cache changed during dynamics");
                require(am.E_ref==ma.E_ref && an.E_ref==ma.E_ref,"original cache changed");
            }
        }
        // Effective pressure can vary across the eight quadrature points.
        for(int qp=0;qp<8;++qp) {
            auto qm=m; auto qs=a; qs.stress=riva_scale(stress,1+.1*qp);
            require(riva_fit_initial_shear(&p.base,&qm,&qs,target),"QP fit");
            riva_moduli_for_state(&p.base,&qm,riva_cone_pressure(&p.base,qs.stress),&qs,&G,&K);
            close(G,target);
        }
        for(double bad : {0.,-1.,std::numeric_limits<double>::infinity()}) {
            auto badm=m;
            require(!riva_fit_initial_shear(&p.base,&badm,&a,bad) &&
                !riva_ib_fit_initial_shear(&p,&badm,&b,bad),"bad target accepted");
            require(badm.E_ref==m.E_ref,"bad fit mutated cache");
        }
        // Modulus regularization must be included in the calibration.
        p.base.p_min=2*units;
        require(riva_fit_initial_shear(&p.base,&ma,&a,target),"floor fit");
        riva_moduli_for_state(&p.base,&ma,riva_cone_pressure(&p.base,a.stress),&a,&G,&K); close(G,target);
    }
}
int main(int argc,char **argv) {
    try { input_tests(); calibration_tests();
        // Optional real-fixture parsing, using the same reader as nonlinear.cc.
        for(int i=1;i<argc;++i) {
            FILE *f=fopen(argv[i],"r"); require(f!=NULL,"fixture open");
            const char *name=strrchr(argv[i],'/'); name=name?name+1:argv[i];
            const int custom=strstr(name,"custom")!=NULL;
            const auto p=riva_reference_parameters(1000); double row[12];
            require(riva_read_material_rows(f,"material_point_constants",1,custom,0,&p,row)==0,"fixture parse");
            fclose(f);
        }
        std::cout<<"PASS: compact/Custom parsing, automatic/full initialized stiffness, Pa/kPa, low pressure, bias/density, per-QP isolation, fixed cache and 1/4/16 restart\n";
    } catch(const std::exception &e) { std::cerr<<e.what()<<'\n'; return 1; }
}
