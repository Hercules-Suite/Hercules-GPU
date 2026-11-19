#ifndef HERCULES_POROELASTIC_MASS_SETUP_H
#define HERCULES_POROELASTIC_MASS_SETUP_H

static biot_mass_t g_biotMass = {};
static std::vector<fvector_t> g_biotMassSolidDamping;
static double g_biotMassStiffnessDampingBound=HUGE_VAL;
static bool g_biotMassHasSolidDamping=false;
static std::vector<solver_float> g_biotMassHost;
static long long g_biotMassSolves[2]={}, g_biotMassIterations[2]={};
static int g_biotMassMaximum[2]={};
static double g_biotMassWorstResidual[2]={};
static int g_biotMassHistoryCount[4]={};
static int g_biotMassProductHistoryCount[4]={};
static double g_biotMassProductDragDt[4]={};
static long long g_biotMassCachedPredictions[2]={},g_biotMassCacheRefreshes[2]={};
static bool g_biotMassFusion=true, g_biotMassCompact=false, g_biotMassHaloBatch=true;
static bool g_biotMassHaloPack=false;
static bool g_biotMassBatchSolves=false;
static bool g_biotMassStabilizationReady=false;
static bool g_biotMassOneRoundRequested=false,g_biotMassOneRound=false;
// Product extrapolation subtracts nearby floating-point vectors. Keep the
// recomputed path in single precision until observable errors are qualified.
static const int g_biotMassPredictorCacheDefault=sizeof(solver_float)==sizeof(double);
static bool g_biotMassPredictorCache=g_biotMassPredictorCacheDefault;
static std::vector<int32_t> g_biotMassIndex;
static acoustic_constraint_map_t g_biotMassConstraints[2];
static int scalar_mpi_share_mask(schedule_t*,int,std::vector<uint8_t>&);
#include "poroelastic_mass_compact_setup.h"
#include "poroelastic_mass_one_round.h"

/* Construct partial assembled rows for globally uniform porous Hex8
 * material/spacing. CPML-adjacent rows have separate inertia/storage weights;
 * purely physical rows retain the compact integer-weight implementation.
 * A changed material, VT cut, or spacing selects the general element path. */
static void biot_mass_build_rows(mesh_t* mesh,const std::vector<int32_t>& elements) {
  if (!Param.thePoroelasticMassOperator) return;
  auto& m=g_biotMass;
  double lo[3]={HUGE_VAL,HUGE_VAL,HUGE_VAL},hi[3]={-HUGE_VAL,-HUGE_VAL,-HUGE_VAL};
  int cut=biot_vt_enabled()?1:0,any_cut;
  for (int pi:elements) {
    const int ei=g_poroelasticElemMapper[pi];
    const edata_t& mat=*(edata_t*)mesh->elemTable[ei].data;
    const double v[3]={mat.edgesize,mat.rho,poroelastic_biot_storage_modulus_for_element(&mat,ei)};
    for (int j=0;j<3;++j) {lo[j]=fmin(lo[j],v[j]);hi[j]=fmax(hi[j],v[j]);}
  }
  double global_lo[3],global_hi[3];
  MPI_Allreduce(lo,global_lo,3,MPI_DOUBLE,MPI_MIN,comm_solver);
  MPI_Allreduce(hi,global_hi,3,MPI_DOUBLE,MPI_MAX,comm_solver);
  MPI_Allreduce(&cut,&any_cut,1,MPI_INT,MPI_MAX,comm_solver);
  bool uniform=!any_cut;
  for (int j=0;j<3;++j) uniform &= isfinite(global_lo[j]) && global_lo[j]>0 && global_lo[j]==global_hi[j];
  if (!uniform) {
    if (Global.myID==0) monitor_print("Biot mass operator: general element fallback (VT, variable material or spacing).\n");
    return;
  }
  std::vector<int32_t> row_of(m.nodes,-1),nodes;
  std::vector<unsigned char> stretched(m.nodes,0);
  for (int pi:elements) {
    const int ei=g_poroelasticElemMapper[pi];
    if (!g_poroCpmlByElement.empty() && g_poroCpmlByElement[ei]>=0)
      for (int i=0;i<8;++i) stretched[mesh->elemTable[ei].lnid[i]]=1;
  }
  for (int pi:elements) {
    const int ei=g_poroelasticElemMapper[pi];
    for (int i=0;i<8;++i) {
      const int n=mesh->elemTable[ei].lnid[i];
      if (row_of[n]<0) {row_of[n]=nodes.size();nodes.push_back(n);}
    }
  }
  std::stable_partition(nodes.begin(),nodes.end(),[&](int n){return !stretched[n];});
  m.uniform_rows=0;
  for (size_t i=0;i<nodes.size();++i) {
    row_of[nodes[i]]=i;
    if (!stretched[nodes[i]]) ++m.uniform_rows;
  }
  const size_t weighted_rows=nodes.size()-m.uniform_rows;
  std::vector<double> mass_weight(27*weighted_rows),storage_weight(27*weighted_rows);
  std::vector<int32_t> neighbor(27*nodes.size(),-1);
  std::vector<unsigned char> weight(27*nodes.size());
  for (int pi:elements) {
    const int ei=g_poroelasticElemMapper[pi];
    const int ci=g_poroCpmlByElement.empty()?-1:g_poroCpmlByElement[ei];
    const double inertia=ci>=0?g_poroCpml[ci].inertia.mass:1;
    const double capacity=ci>=0?poro_cpml_capacity(g_poroCpml[ci]):1;
    for (int i=0;i<8;++i) for (int j=0;j<8;++j) {
      int slot=0,factor=1,stride=1;
      for (int d=0;d<3;++d) {
        const int delta=((j>>d)&1)-((i>>d)&1);
        slot+=stride*(delta+1);stride*=3;
        if (!delta) factor*=2;
      }
      const int n=mesh->elemTable[ei].lnid[i],other=mesh->elemTable[ei].lnid[j];
      const size_t index=(size_t)slot*nodes.size()+row_of[n];
      if ((neighbor[index]>=0 && neighbor[index]!=other) || weight[index]+factor>64)
        solver_abort(__FUNCTION_NAME,NULL,"Invalid regular Q1 mass row topology\n");
      neighbor[index]=other;weight[index]+=factor;
      if (stretched[n]) {
        const size_t wi=(size_t)slot*weighted_rows+row_of[n]-m.uniform_rows;
        mass_weight[wi]+=factor*inertia;
        storage_weight[wi]+=factor*capacity;
      }
    }
  }
  m.rows=nodes.size();m.row_volume=pow(global_lo[0],3);
  m.row_rho=global_lo[1];m.row_capacity=1/global_lo[2];
  int local_weighted=weighted_rows>0,any_weighted;
  MPI_Allreduce(&local_weighted,&any_weighted,1,MPI_INT,MPI_MAX,comm_solver);
  if(g_biotMassOneRoundRequested && !any_weighted)
    biot_mass_one_round_setup(mesh,nodes,neighbor,weight);
  m.rows=nodes.size();
  if (!any_weighted) m.uniform_rows=m.rows;
  acoustic_constraint_upload(&m.row_node,nodes);
  acoustic_constraint_upload(&m.neighbor,neighbor);
  acoustic_constraint_upload(&m.weight,weight);
  acoustic_constraint_upload(&m.row_mass_weight,mass_weight);
  acoustic_constraint_upload(&m.row_storage_weight,storage_weight);
  if (Global.myID==0) monitor_print("Biot mass operator: assembled uniform Q1 rows "
      "(27 neighbors, no element atomics); %s.\n",any_weighted?
      "weighted CPML inertia/storage rows, integer physical rows":"integer physical rows");
}

static void biot_mass_initialize() {
  if (!Param.thePoroelasticMassMatrix) return;
  const char* fusion=getenv("HERCULES_BIOT_MASS_FUSION");
  int choice=1,lo,hi;
  if (fusion && strcmp(fusion,"0") && strcmp(fusion,"1"))
    solver_abort(__FUNCTION_NAME,NULL,"HERCULES_BIOT_MASS_FUSION must be 0 or 1\n");
  if (fusion) choice=strcmp(fusion,"0")!=0;
  MPI_Allreduce(&choice,&lo,1,MPI_INT,MPI_MIN,comm_solver);
  MPI_Allreduce(&choice,&hi,1,MPI_INT,MPI_MAX,comm_solver);
  if (lo!=hi) solver_abort(__FUNCTION_NAME,NULL,"Biot mass fusion settings differ between ranks\n");
  g_biotMassFusion=choice;
  g_biotMassCompact=biot_mass_option("HERCULES_BIOT_MASS_COMPACT",0);
  g_biotMassHaloBatch=biot_mass_option("HERCULES_BIOT_MASS_HALO_BATCH");
  g_biotMassHaloPack=biot_mass_option("HERCULES_BIOT_MASS_HALO_PACK",0);
  g_biotMassBatchSolves=biot_mass_option("HERCULES_BIOT_MASS_BATCH_SOLVES",0);
  g_biotMassOneRoundRequested=biot_mass_option("HERCULES_BIOT_MASS_ONE_ROUND",0);
  if (g_biotMassBatchSolves && Param.thePoroelasticMassKrylov) {
    if(Global.myID==0) monitor_print("Biot solve batching: sequential fallback for single-reduction CG.\n");
    g_biotMassBatchSolves=false;
  }
  g_biotMassPredictorCache=biot_mass_option("HERCULES_BIOT_MASS_PREDICTOR_CACHE",g_biotMassPredictorCacheDefault) &&
      Param.thePoroelasticMassInitialGuess!=0;
  mesh_t* mesh=Global.myMesh;
  mysolver_t* s=Global.mySolver;
  const int N=mesh->nharbored;
  biot_mass_t& m=g_biotMass;
  m.nodes=m.full_nodes=N;
  m.dt=Param.theDeltaT;
  g_biotMassSolidDamping.resize(N);
  biot_assemble_nodes(g_biotMassSolidDamping);
  int local_damping=0,global_damping=0;
  for (const auto& d:g_biotMassSolidDamping) for (int axis=0;axis<3;++axis) {
    if (!isfinite(d.f[axis]) || d.f[axis]<0)
      solver_abort(__FUNCTION_NAME,NULL,"Consistent mass requires nonnegative physical solid damping\n");
    local_damping|=d.f[axis]>0;
  }
  MPI_Allreduce(&local_damping,&global_damping,1,MPI_INT,MPI_MAX,comm_solver);
  g_biotMassHasSolidDamping=global_damping;
  if (global_damping) acoustic_constraint_upload(&m.solid_damping,g_biotMassSolidDamping);
  std::vector<fvector_t>().swap(g_biotMassSolidDamping);
  m.rho_f=Param.thePoroelasticFluidDensity;
  m.fluid_mass=Param.thePoroelasticTortuosity*m.rho_f/Param.thePoroelasticPorosity;
  std::vector<fvector_t> removed(N), fluid(N);
  std::vector<solver_float> capacity(N);
  std::vector<double> rho;
  std::vector<int32_t> elements;
  double bound=g_biotMassStiffnessDampingBound;
  for (int pi=0;pi<g_poroelasticLocalElementCount;++pi) {
    const int ei=g_poroelasticElemMapper[pi];
    const edata_t& mat=*(edata_t*)mesh->elemTable[ei].data;
    const double h=mat.edgesize;
    double fast2,slow2;
    biot_compressional_speeds_squared(mat.rho,m.rho_f,m.fluid_mass,
        mat.rho*mat.Vp*mat.Vp,Param.thePoroelasticBiotAlpha,
        poroelastic_biot_storage_modulus_for_element(&mat,ei),&fast2,&slow2);
    /* Suggested consistent-mass wave estimate, including VT and gamma=0.
     * Exceeding a timestep estimate warns without changing the requested step.
     * CPML's geometry/interface restrictions remain independently in force. */
    bound=fmin(bound,(biot_vt_enabled()?0.02:0.2)*h/sqrt(fast2));
    const double k[3]={Param.thePoroelasticPermeabilityX,
        Param.thePoroelasticPermeabilityY,Param.thePoroelasticPermeabilityZ};
    for (int d=0;d<3;++d) if (k[d]>0) {
      const double nu=Param.thePoroelasticPressureStabilization*sqrt(slow2)*h*h*h/
          (1+Param.thePoroelasticFluidViscosity*h/(k[d]*m.fluid_mass*sqrt(slow2)));
      if (nu>0) bound=fmin(bound,pow(h,4)/((biot_vt_enabled()?14400:1296)*nu));
    }
    // CPML's complete inertia/storage history uses the same Q1 integration.
    const int cpml=g_poroCpmlByElement.empty()?-1:g_poroCpmlByElement[ei];
    const double inertia=cpml>=0?g_poroCpml[cpml].inertia.mass:1;
    const double storage=cpml>=0?poro_cpml_capacity(g_poroCpml[cpml]):1;
    elements.push_back(pi); rho.push_back(mat.rho);
    double volume[8],total;
    poroelastic_element_nodal_volumes(mesh,ei,volume,&total);
    const double M=poroelastic_biot_storage_modulus_for_element(&mat,ei);
    for (int j=0;j<8;++j) {
      const int n=mesh->elemTable[ei].lnid[j];
      removed[n].f[0]+=mat.rho*volume[j]*inertia;
      removed[n].f[1]+=m.rho_f*volume[j]*inertia;
      removed[n].f[2]+=m.fluid_mass*volume[j]*inertia;
      capacity[n]+=volume[j]/M*storage;
    }
  }
  double global_bound;
  MPI_Allreduce(&bound,&global_bound,1,MPI_DOUBLE,MPI_MIN,comm_solver);
  warn_timestep_estimate("consistent Biot wave/stabilization/stiffness-damping",global_bound);
  if (Global.myID==0)
    monitor_print("Consistent Biot mass estimated dt <= %.12g s (advisory only).\n",global_bound);
  biot_assemble_nodes(removed);
  schedule_senddata(s->dn_sched,capacity.data(),1,CONTRIBUTION,DN_PORE_PRESSURE_MASS_MSG);
  compute_adjust(capacity.data(),1,DISTRIBUTION,ACOUSTIC_CONSTRAINT_PORE);
  schedule_senddata(s->an_sched,capacity.data(),1,CONTRIBUTION,AN_PORE_PRESSURE_MASS_MSG);
  cudaMemcpy(fluid.data(),Global.gpuData.biotMassDevice,N*sizeof(fvector_t),cudaMemcpyDeviceToHost);
  std::vector<double> remainder(4*(size_t)N);
  std::vector<unsigned char> owned(N);
  std::vector<int32_t> boundary(N,-1);
  std::vector<biot_boundary_t> b(Global.gpuData.biotBoundaryCount);
  if (!b.empty()) cudaMemcpy(b.data(),Global.gpuData.biotBoundaryDevice,
      b.size()*sizeof(biot_boundary_t),cudaMemcpyDeviceToHost);
  for (size_t j=0;j<b.size();++j) boundary[b[j].node]=j;
  for (int n=0;n<N;++n) {
    if (mesh->nodeTable[n].ismine) owned[n]=3;
    const double base[4]={s->nTable[n].mass_simple,fluid[n].f[0],fluid[n].f[1],s->pore_pressure_mass[n]};
    const double sub[4]={removed[n].f[0],removed[n].f[1],removed[n].f[2],capacity[n]};
    for (int j=0;j<4;++j) {
      double v=base[j]-sub[j];
      const double roundoff=sizeof(solver_float)==8?1e-12:1e-5;
      if (fabs(v)<roundoff*fmax(fabs(base[j]),fabs(sub[j]))) v=0;
      if (mesh->nodeTable[n].ismine && v<0)
        solver_abort(__FUNCTION_NAME,NULL,"Negative consistent mass remainder at node %d\n",n);
      remainder[4*(size_t)n+j]=v;
    }
  }
  // Only independent owned unknowns enter Krylov inner products. Released
  // field-specific dangling nodes are independent; retained nodes are not.
  const int fields[2]={ACOUSTIC_CONSTRAINT_SOLID,ACOUSTIC_CONSTRAINT_PORE};
  for (int j=0;j<2;++j) {
    const auto& map=g_acousticConstraintMaps[fields[j]];
    if (map.enabled) {
      for (int base : map.node_base3) owned[base/3]&=~(1<<j);
    } else for (int n=0;n<N;++n)
      if (!mesh->nodeTable[n].isanchored) owned[n]&=~(1<<j);
  }
  m.elements=elements.size();
  if (!g_biotSurfaceFixed.empty())
    for (int n=0;n<N;++n) if (g_biotSurfaceFixed[n]) owned[n]&=~2;
  acoustic_constraint_upload(&m.element,elements);
  acoustic_constraint_upload(&m.rho,rho);
  acoustic_constraint_upload(&m.remainder,remainder);
  acoustic_constraint_upload(&m.owned,owned);
  acoustic_constraint_upload(&m.boundary,boundary);
  biot_mass_build_rows(mesh,elements);
  if(g_biotMassOneRoundRequested && !g_biotMassOneRound && Global.myID==0)
    monitor_print("Biot one-round mass: two-exchange fallback; requires uniform Cartesian rows, no porous CPML weights, full vectors and no geometric constraints.\n");
  biot_mass_compact_setup(mesh,elements,removed,capacity,owned);
  const size_t slots=std::max(1,m.nodes);
  solver_float** buffers[]={&m.x,&m.r,&m.z,&m.direction,&m.product,&m.rhs};
  for (auto ptr:buffers) {
    if (GPU_SETUP_MALLOC(ptr,7*slots*sizeof(solver_float))!=cudaSuccess)
      solver_abort(__FUNCTION_NAME,NULL,"Cannot allocate consistent Biot mass workspace\n");
    cudaMemset(*ptr,0,7*slots*sizeof(solver_float));
  }
  if (Param.thePoroelasticMassPreconditioner &&
      GPU_SETUP_MALLOC(&m.precondition,7*slots*sizeof(solver_float))!=cudaSuccess)
    solver_abort(__FUNCTION_NAME,NULL,"Cannot allocate Biot element-inverse workspace\n");
  if (Param.thePoroelasticMassKrylov &&
      GPU_SETUP_MALLOC(&m.applied,7*slots*sizeof(solver_float))!=cudaSuccess)
    solver_abort(__FUNCTION_NAME,NULL,"Cannot allocate single-reduction Biot workspace\n");
  size_t cache_bytes=0;
  if (Param.thePoroelasticMassInitialGuess) for (int j=0;j<4;++j) {
    if (j==1 && !g_biotMassHasSolidDamping &&
        !(Param.thePoroelasticFluidViscosity>0 && !poroelastic_zero_permeability())) continue;
    const size_t history_bytes=(j<2?6:1)*slots*sizeof(solver_float);
    for (int k=0;k<Param.thePoroelasticMassInitialGuess;++k) {
      if (GPU_SETUP_MALLOC(&m.history[j][k],history_bytes)!=cudaSuccess)
        solver_abort(__FUNCTION_NAME,NULL,"Cannot allocate Biot solution-history workspace\n");
      if (g_biotMassPredictorCache) {
        if (GPU_SETUP_MALLOC(&m.history_product[j][k],history_bytes)!=cudaSuccess ||
            GPU_SETUP_MALLOC(&m.history_preconditioned_product[j][k],history_bytes)!=cudaSuccess)
          solver_abort(__FUNCTION_NAME,NULL,"Cannot allocate Biot predictor-product workspace\n");
        cache_bytes+=2*history_bytes;
      }
    }
    g_biotMassHistoryCount[j]=0;
    g_biotMassProductHistoryCount[j]=0;
  }
  unsigned long long local_cache_bytes=cache_bytes,max_cache_bytes=0;
  MPI_Allreduce(&local_cache_bytes,&max_cache_bytes,1,MPI_UNSIGNED_LONG_LONG,MPI_MAX,comm_solver);
  if (GPU_SETUP_MALLOC(&m.dots,8*sizeof(double))!=cudaSuccess)
    solver_abort(__FUNCTION_NAME,NULL,"Cannot allocate mass reduction workspace\n");
  cudaMemset(m.dots,0,8*sizeof(double));
  if (!g_biotMassCompact) g_biotMassHost.resize(3*(size_t)N);
  long long local=m.elements,total;
  MPI_Allreduce(&local,&total,1,MPI_LONG_LONG,MPI_SUM,comm_solver);
  if (Global.myID==0) monitor_print(
      "Biot mass: consistent porous cells (%lld), including CPML, native Hex8/P1 integrals; "
      "lumped nonporous remainder. Projected MPI PCG, relative "
      "preconditioned residual %.3g, limit %d; conservative dt <= %.9g s. "
      "Darcy and physical porous-cell solid mass damping/dashpots remain lumped, "
      "integrated by two-stage L-stable SDIRK2. CPML and nonporous damping "
      "retain their existing integration; stiffness damping remains explicit.\n",
      total,Param.thePoroelasticMassTolerance,Param.thePoroelasticMassMaxIterations,global_bound);
  if (Global.myID==0) monitor_print("Biot mass solver: %s CG, %s preconditioner, %s initial guess.\n",
      Param.thePoroelasticMassKrylov?"single-reduction":"standard",
      Param.thePoroelasticMassPreconditioner?"element-inverse":"lumped",
      Param.thePoroelasticMassInitialGuess==3?"quadratic per-stage":
      Param.thePoroelasticMassInitialGuess==2?"linear per-stage":
      Param.thePoroelasticMassInitialGuess?"previous per-stage":"zero");
  if (Global.myID==0) monitor_print("Biot mass local kernels: %s; constraints and true-residual verification retained.\n",
      g_biotMassFusion?"fused":"separate (comparison path)");
  if (Global.myID==0) monitor_print("Biot predictor products: %s; maximum cache %.3f GiB/rank; "
      "fixed mass/preconditioner, separate stage histories, fresh residual before acceptance.\n",
      g_biotMassPredictorCache?"cached":"recomputed",
      (double)max_cache_bytes/(1024.*1024.*1024.));
  if(Global.myID==0) monitor_print("Biot solve batching: %s; independent inertia/stabilization tolerances, histories and fresh residual checks.\n",
      g_biotMassBatchSolves?"experimental shared iterations":"sequential");
}

#endif
