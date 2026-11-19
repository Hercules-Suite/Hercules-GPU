#ifndef HERCULES_POROELASTIC_DRAIN_SETUP_H
#define HERCULES_POROELASTIC_DRAIN_SETUP_H

/* Weak zero-pressure flux: R=lump(P^T R_face P), R_face=integral(Y Ni Nj).
 * P preserves constants and has nonnegative weights, so row lumping retains
 * total boundary flux and positivity. Y is independent of the timestep.
 * The stored diagonal B=dt R/2 enters a trapezoidal pressure update. */
static std::vector<solver_float> g_biotSurfaceRate;
static std::vector<unsigned char> g_biotSurfaceFixed;
static double g_biotSurfaceArea=0;
static FILE* g_biotSurfaceBalanceFile=nullptr;
static long long g_biotSurfaceStep=0;

static bool biot_surface_weak() {
  return full_biot_enabled() && Param.thePoroelasticSurfacePressureEnforcement==1 &&
      Param.thePoroelasticSurfacePressureBC==POROELASTIC_SURFACE_PRESSURE_DRAINED &&
      Param.thePoroelasticSurfacePressurePenaltyFactor>0;
}

static bool biot_surface_strong() {
  return poroelastic_enabled() && Param.thePoroelasticSurfacePressureEnforcement==2 &&
      Param.thePoroelasticSurfacePressureBC==POROELASTIC_SURFACE_PRESSURE_DRAINED;
}

static bool biot_surface_enforced() {
  return biot_surface_weak() || biot_surface_strong();
}

static double biot_surface_admittance(int ei,const double normal[3]) {
  const edata_t& mat=*(edata_t*)Global.myMesh->elemTable[ei].data;
  const double k[3]={Param.thePoroelasticPermeabilityX,
      Param.thePoroelasticPermeabilityY,Param.thePoroelasticPermeabilityZ};
  double allowed=0;
  for (int d=0;d<3;++d) if (k[d]>0) allowed+=normal[d]*normal[d];
  if (!allowed) return 0; // a physically impermeable normal cannot drain
  // Strong enforcement needs only positive trace support. It is independent
  // of the weak penalty factor and does not introduce a numerical admittance.
  if (biot_surface_strong()) return 1;
  const double rf=Param.thePoroelasticFluidDensity;
  const double M=poroelastic_biot_storage_modulus_for_element(&mat,ei);
  double fast2,slow2;
  biot_compressional_speeds_squared(mat.rho,rf,
      Param.thePoroelasticTortuosity*rf/Param.thePoroelasticPorosity,
      mat.rho*mat.Vp*mat.Vp,Param.thePoroelasticBiotAlpha,M,&fast2,&slow2);
  return Param.thePoroelasticSurfacePressurePenaltyFactor*sqrt(fast2)/M*allowed;
}

static void biot_surface_add_face(int ei,int axis,int sign,
    double a,double b,double width,bool groundwater=false) {
  if (!biot_surface_enforced()) return;
  const elem_t& e=Global.myMesh->elemTable[ei];
  const double h=((edata_t*)e.data)->edgesize;
  const node_t& n0=Global.myMesh->nodeTable[e.lnid[0]];
  const double tick=Global.myMesh->ticksize;
  const int t0=(axis+1)%3,t1=(axis+2)%3;
  double normal[3]={};normal[axis]=sign;
  const double Y=biot_surface_admittance(ei,normal);
  if (!(Y>0)) return;
  for (int q=0;q<4;++q) {
    double x[3]={};x[axis]=sign;
    x[t0]=2*(a+0.5*width*(1+((q&1)?1:-1)*0.57735026918962576451))/h-1;
    x[t1]=2*(b+0.5*width*(1+((q&2)?1:-1)*0.57735026918962576451))/h-1;
    double bed;
    if (!groundwater && topo_seabed_enabled() && topo_seabed_sample(n0.x*tick+0.5*h*(1+x[0]),
                          n0.y*tick+0.5*h*(1+x[1]),&bed)) continue;
    g_biotSurfaceArea+=0.25*width*width;
    for (int j=0;j<8;++j) {
      double N=0.125;
      for (int d=0;d<3;++d) N*=1+(((j>>d)&1)?1:-1)*x[d];
      g_biotSurfaceRate[e.lnid[j]]+=Y*0.25*width*width*N;
    }
  }
}

/* A buried discrete groundwater boundary can cross the exterior face of a
 * native VT terrain cell. Integrate its retained P1 face triangles, clipped
 * to the neighbor fragment and the linear terrain contour. Q1 face weights
 * are not the pressure trace in these cells. No work enters the time loop. */
static void biot_surface_add_vt_groundwater_face(int ei,const biot_cut_t& cut,
    int axis,int sign,double a,double b,double width) {
  if (!biot_surface_enforced()) return;
  const elem_t& e=Global.myMesh->elemTable[ei];
  const double h=((edata_t*)e.data)->edgesize,tick=Global.myMesh->ticksize;
  const node_t& n0=Global.myMesh->nodeTable[e.lnid[0]];
  const double low[3]={n0.x*tick,n0.y*tick,n0.z*tick};
  const int t0=(axis+1)%3,t1=(axis+2)%3;
  double normal[3]={};normal[axis]=sign;
  const double Y=biot_surface_admittance(ei,normal);
  if (!(Y>0)) return;
  struct vertex { double x[3],N[4],terrain; };
  for (int t=0;t<5;++t) {
    if (!(cut.fraction[t]>0)) continue;
    int ids[4];double grad[4][3];
    biot_tet_geometry(cut.partition,t,ids,grad);
    std::vector<vertex> polygon;
    for (int j=0;j<4;++j) {
      if (((ids[j]>>axis)&1)!=(sign>0?1:0)) continue;
      vertex v={};v.N[j]=1;
      for (int d=0;d<3;++d) v.x[d]=(ids[j]>>d)&1;
      v.terrain=low[2]+h*v.x[2]-point_elevation(low[0]+h*v.x[0],low[1]+h*v.x[1]);
      polygon.push_back(v);
    }
    if (polygon.size()!=3) continue;
    auto clip=[&](int direction,double boundary,bool above) {
      std::vector<vertex> result;
      if (polygon.empty()) return;
      for (size_t j=0;j<polygon.size();++j) {
        const vertex& u=polygon[j];const vertex& v=polygon[(j+1)%polygon.size()];
        const double fu=(direction<0?u.terrain:u.x[direction])-boundary;
        const double fv=(direction<0?v.terrain:v.x[direction])-boundary;
        const bool iu=above?fu>=0:fu<=0,iv=above?fv>=0:fv<=0;
        if (iu) result.push_back(u);
        if (iu!=iv) {
          const double s=fu/(fu-fv);vertex q={};
          for (int d=0;d<3;++d) q.x[d]=(1-s)*u.x[d]+s*v.x[d];
          for (int k=0;k<4;++k) q.N[k]=(1-s)*u.N[k]+s*v.N[k];
          q.terrain=(1-s)*u.terrain+s*v.terrain;result.push_back(q);
        }
      }
      polygon.swap(result);
    };
    clip(t0,a/h,true);clip(t0,(a+width)/h,false);
    clip(t1,b/h,true);clip(t1,(b+width)/h,false);clip(-1,0,true);
    for (size_t j=1;j+1<polygon.size();++j) {
      const vertex* v[3]={&polygon[0],&polygon[j],&polygon[j+1]};
      const double area=.5*h*h*fabs((v[1]->x[t0]-v[0]->x[t0])*(v[2]->x[t1]-v[0]->x[t1])-
                                  (v[1]->x[t1]-v[0]->x[t1])*(v[2]->x[t0]-v[0]->x[t0]));
      if (!(area>0)) continue;
      g_biotSurfaceArea+=area;
      // P1 traces integrate exactly by the triangle centroid, with positive weights.
      for (int k=0;k<4;++k)
        g_biotSurfaceRate[e.lnid[ids[k]]]+=Y*area*(v[0]->N[k]+v[1]->N[k]+v[2]->N[k])/3;
    }
  }
}

/* Native VT does not retain a triangulated cut surface. Reconstruct a linear
 * terrain level set in each retained P1 tetrahedron for boundary quadrature.
 * Integrate its triangular/quadrilateral zero contour with P1 trace weights.
 * Interior curvature invisible to all four vertices requires mesh refinement;
 * reject that unresolved case instead of silently omitting a draining patch. */
static void biot_surface_add_cut(int ei,const biot_cut_t& cut) {
  if (!biot_surface_enforced()) return;
  const elem_t& e=Global.myMesh->elemTable[ei];
  const double h=((edata_t*)e.data)->edgesize,tick=Global.myMesh->ticksize;
  const node_t& n0=Global.myMesh->nodeTable[e.lnid[0]];
  const double low[3]={n0.x*tick,n0.y*tick,n0.z*tick};
  for (int t=0;t<5;++t) {
    if (!(cut.fraction[t]>0)) continue;
    int ids[4];double grad[4][3],f[4];
    biot_tet_geometry(cut.partition,t,ids,grad);
    for (int i=0;i<4;++i) {
      const double x=low[0]+h*(ids[i]&1),y=low[1]+h*((ids[i]>>1)&1);
      f[i]=low[2]+h*((ids[i]>>2)&1)-point_elevation(x,y);
      if (fabs(f[i])<1e-10*h) f[i]=0;
    }
    struct point { double x[3],N[4],angle; } v[6]={};int count=0;
    for (int i=0;i<4;++i) if (f[i]==0) {
      v[count].N[i]=1;
      for (int d=0;d<3;++d) v[count].x[d]=(ids[i]>>d)&1;
      ++count;
    }
    for (int i=0;i<4;++i) for (int j=i+1;j<4;++j) if (f[i]*f[j]<0) {
      const double s=f[i]/(f[i]-f[j]);
      v[count].N[i]=1-s;v[count].N[j]=s;
      for (int d=0;d<3;++d)
        v[count].x[d]=(1-s)*((ids[i]>>d)&1)+s*((ids[j]>>d)&1);
      ++count;
    }
    if (count<3) {
      if (cut.fraction[t]<1-1e-10)
        solver_abort(__FUNCTION_NAME,NULL,
            "VT drainage cannot resolve the terrain crossing from tetrahedron "
            "vertices (element %d, tet %d); refine the terrain mesh.\n",ei,t);
      continue;
    }
    double center[3]={},normal[3]={},u[3],w[3],norm=0,unorm=0;
    for (int d=0;d<3;++d) {
      for (int i=0;i<count;++i) center[d]+=v[i].x[d]/count;
      for (int i=0;i<4;++i) normal[d]+=f[i]*grad[i][d];
      norm+=normal[d]*normal[d];
      u[d]=v[0].x[d]-center[d];unorm+=u[d]*u[d];
    }
    if (!(norm>0 && unorm>0))
      solver_abort(__FUNCTION_NAME,NULL,"Degenerate VT drainage surface quadrature\n");
    for (int d=0;d<3;++d) {normal[d]/=-sqrt(norm);u[d]/=sqrt(unorm);}
    for (int d=0;d<3;++d)
      w[d]=normal[(d+1)%3]*u[(d+2)%3]-normal[(d+2)%3]*u[(d+1)%3];
    for (int i=0;i<count;++i) {
      double a=0,b=0;
      for (int d=0;d<3;++d) {a+=(v[i].x[d]-center[d])*u[d];b+=(v[i].x[d]-center[d])*w[d];}
      v[i].angle=atan2(b,a);
    }
    std::sort(v,v+count,[](const point& a,const point& b){return a.angle<b.angle;});
    const double Y=biot_surface_admittance(ei,normal);
    if (!(Y>0)) continue;
    for (int i=1;i<count-1;++i) {
      double cross[3],aa[3],bb[3],area=0;
      for (int d=0;d<3;++d) {aa[d]=v[i].x[d]-v[0].x[d];bb[d]=v[i+1].x[d]-v[0].x[d];}
      for (int d=0;d<3;++d) {cross[d]=aa[(d+1)%3]*bb[(d+2)%3]-aa[(d+2)%3]*bb[(d+1)%3];area+=cross[d]*cross[d];}
      area=0.5*h*h*sqrt(area);
      const point* triangle[3]={&v[0],&v[i],&v[i+1]};
      // Positive three-point quadrature. The explicit seabed mask excludes
      // submerged traces, including a shoreline crossing a native cut cell.
      // Shoreline integration, like the VT terrain contour, needs refinement.
      for (int q=0;q<3;++q) {
        double xy[3]={},N[4]={};
        for (int a=0;a<3;++a) {
          const double weight=a==q?2.0/3:1.0/6;
          for (int d=0;d<3;++d) xy[d]+=weight*triangle[a]->x[d];
          for (int j=0;j<4;++j) N[j]+=weight*triangle[a]->N[j];
        }
        double bed;
        if (topo_seabed_enabled() && topo_seabed_sample(low[0]+h*xy[0],low[1]+h*xy[1],&bed)) continue;
        g_biotSurfaceArea+=area/3;
        if (full_biot_enabled()) {
          for (int j=0;j<4;++j)
            g_biotSurfaceRate[e.lnid[ids[j]]]+=Y*area*N[j]/3;
        } else {
          // BrickUP retains its Q1 pressure field on a VT volume-fraction
          // cell; the native full-Biot P1 trace is not its pressure basis.
          // Positive Q1 support also guarantees zero station interpolation.
          for (int j=0;j<8;++j) {
            double shape=1;
            for (int d=0;d<3;++d) shape*=((j>>d)&1)?xy[d]:1-xy[d];
            g_biotSurfaceRate[e.lnid[j]]+=Y*area*shape/3;
          }
        }
      }
    }
  }
}

static void biot_surface_initialize() {
  if (!biot_surface_enforced()) return;
  mesh_t* mesh=Global.myMesh;mysolver_t* s=Global.mySolver;
  auto& rate=g_biotSurfaceRate;
  schedule_senddata(s->dn_sched,rate.data(),1,CONTRIBUTION,DN_PORE_PRESSURE_MASS_MSG);
  compute_adjust(rate.data(),1,DISTRIBUTION,ACOUSTIC_CONSTRAINT_PORE);
  schedule_senddata(s->an_sched,rate.data(),1,CONTRIBUTION,AN_PORE_PRESSURE_MASS_MSG);
  std::vector<unsigned char> owned(mesh->nharbored);
  for (int n=0;n<mesh->nharbored;++n) owned[n]=mesh->nodeTable[n].ismine && s->pore_pressure_mass[n]>0;
  const auto& map=g_acousticConstraintMaps[ACOUSTIC_CONSTRAINT_PORE];
  if (map.enabled) for (int base:map.node_base3) owned[base/3]=0;
  else for (int n=0;n<mesh->nharbored;++n) if (!mesh->nodeTable[n].isanchored) owned[n]=0;
  double local_sum=0,global_sum=0,area=0;
  for (int n=0;n<mesh->nharbored;++n) {
    if (!isfinite(rate[n]) || rate[n]<0)
      solver_abort(__FUNCTION_NAME,NULL,"Invalid weak surface drainage weight\n");
    if (owned[n]) local_sum+=rate[n];
    if (!biot_surface_strong()) rate[n]*=0.5*Param.theDeltaT;
  }
  MPI_Allreduce(&local_sum,&global_sum,1,MPI_DOUBLE,MPI_SUM,comm_solver);
  MPI_Reduce(&g_biotSurfaceArea,&area,1,MPI_DOUBLE,MPI_SUM,0,comm_solver);
  if (biot_surface_strong()) {
    // Nonnegative trace and hanging weights: fixing every independent parent
    // with positive P^T*r support guarantees a zero reconstructed trace.
    // This is deliberately a support restriction, not minimal cut constraints:
    // VT can lose the pressure gradient throughout a cut tetrahedron, and a
    // nonmatching face can also fix coarse trace support outside the fragment.
    g_biotSurfaceFixed.assign(mesh->nharbored,0);
    long long fixed=0,global_fixed=0;
    for (int n=0;n<mesh->nharbored;++n) if (owned[n] && rate[n]>0) {
      g_biotSurfaceFixed[n]=1;++fixed;
    }
    MPI_Allreduce(&fixed,&global_fixed,1,MPI_LONG_LONG,MPI_SUM,comm_solver);
    if (global_fixed) acoustic_constraint_upload(
        &Global.gpuData.biotSurfaceFixedDevice,g_biotSurfaceFixed);
    if (Global.myID==0 && full_biot_enabled()) monitor_print(
        "Full-Biot strong drained surface: %lld independent owned pressure DOFs fixed; "
        "permeable dry area=%.12g m^2. Pressure and stabilization auxiliary are zero "
        "on the exposed trace; normal relative fluid motion remains free. "
        "Positive trace support is closed through pore constraints. VT fixes cut-surface "
        "support nodes, suppressing pressure variation within affected tetrahedra; "
        "nonmatching faces/shorelines can extend the fixed support. Refine to assess "
        "this boundary approximation. Weak penalty factor and weak balance audit are inactive.\n",
        global_fixed,area);
    if (Global.myID==0 && !full_biot_enabled()) monitor_print(
        "BrickUP strong drained surface: %lld independent owned pressure DOFs fixed; "
        "permeable dry area=%.12g m^2. Zero excess pressure in the explicit lumped "
        "update; legacy compliance (including CPML drain histories) disabled. "
        "Positive Q1 trace support is closed through pore constraints. Interior "
        "VT/groundwater cuts can fix every active pressure node of a boundary "
        "brick and suppress its pressure gradient; nonmatching/shoreline support "
        "can extend beyond the dry patch. Experimental on cuts: refine and check "
        "drainage and interior pressure. The weak penalty factor is inactive.\n",
        global_fixed,area);
    std::vector<solver_float>().swap(rate);return;
  }
  if (!(global_sum>0)) {
    if (Global.myID==0) monitor_print("Full-Biot weak drainage: no permeable exposed surface; no additional flux.\n");
    std::vector<solver_float>().swap(rate);return;
  }
  acoustic_constraint_upload(&Global.gpuData.biotDrainCapacityDevice,rate);
  const char* env=getenv("HERCULES_BIOT_SURFACE_DIAGNOSTICS");
  if (env && strcmp(env,"0") && strcmp(env,"1"))
    solver_abort(__FUNCTION_NAME,NULL,"HERCULES_BIOT_SURFACE_DIAGNOSTICS must be 0 or 1\n");
  int diagnostic=env && !strcmp(env,"1"),lo,hi;
  MPI_Allreduce(&diagnostic,&lo,1,MPI_INT,MPI_MIN,comm_solver);
  MPI_Allreduce(&diagnostic,&hi,1,MPI_INT,MPI_MAX,comm_solver);
  if (lo!=hi) solver_abort(__FUNCTION_NAME,NULL,"Surface diagnostic settings differ between ranks\n");
  if (diagnostic) {
    acoustic_constraint_upload(&Global.gpuData.biotDrainOwnedDevice,owned);
    if (GPU_SETUP_MALLOC(&Global.gpuData.biotDrainDiagnosticsDevice,6*sizeof(double))!=cudaSuccess)
      solver_abort(__FUNCTION_NAME,NULL,"Cannot allocate surface diagnostics\n");
    if (Global.myID==0) {
      // Station setup is skipped when no receivers are requested. The audit
      // remains useful on its own and must still honor its output directory.
      if (!Param.theNumberOfStations) {
        FILE* input=hu_fopen(Param.parameters_input_file,"r");
        parsetext(input,"output_stations_directory",'s',Param.theStationsDirOut);
        hu_fclose(input);
        createDir(Param.theStationsDirOut,__FUNCTION_NAME);
      }
      char path[1024];
      snprintf(path,sizeof(path),"%s/biot_surface_balance.txt",Param.theStationsDirOut);
      g_biotSurfaceBalanceFile=hu_fopen(path,"w");
      fprintf(g_biotSurfaceBalanceFile,
          "# t_s pressure_balance_residual_m3 residual_scale_m3 "
          "penalty_outflow_m3 penalty_dissipation_J storage_change_m3 "
          "unpenalized_rhs_increment_m3\n"
          "# residual = storage_change + penalty_outflow - rhs_increment; "
          "rhs includes physical divergence, sources and stabilization. "
          "This is not a whole-system CPML energy diagnostic.\n");
      monitor_print("Full-Biot surface diagnostics: %s\n",path);
    }
  }
  if (Global.myID==0) monitor_print(
      "Full-Biot weak drained surface: beta=%.9g; area=%.12g m^2; "
      "integrated admittance=%.12g m^3/(Pa s). Y=beta*c_fast/M in permeable "
      "normal modes; row-lumped surface flux, time-centered C+dt*R/2 update. "
      "Normal relative motion is free; stabilization auxiliary retains C. "
      "Station q is the volume field; numerical boundary outflow adds Y*p. "
      "Balance diagnostics %s.\n",Param.thePoroelasticSurfacePressurePenaltyFactor,
      area,global_sum,diagnostic?"enabled":"disabled");
  std::vector<solver_float>().swap(rate);
}

static void biot_surface_finish() {
  cudaFree(Global.gpuData.biotDrainCapacityDevice);
  cudaFree(Global.gpuData.biotSurfaceFixedDevice);
  cudaFree(Global.gpuData.biotDrainOwnedDevice);
  cudaFree(Global.gpuData.biotDrainDiagnosticsDevice);
  if (g_biotSurfaceBalanceFile) hu_fclose(g_biotSurfaceBalanceFile);
  g_biotSurfaceBalanceFile=nullptr;g_biotSurfaceStep=0;g_biotSurfaceArea=0;
  std::vector<solver_float>().swap(g_biotSurfaceRate);
  std::vector<unsigned char>().swap(g_biotSurfaceFixed);
}
#endif
