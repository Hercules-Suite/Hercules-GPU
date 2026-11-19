#ifndef HERCULES_POROELASTIC_MASS_PROFILE_H
#define HERCULES_POROELASTIC_MASS_PROFILE_H

/* Optional diagnostic intervals, not a hardware bandwidth measurement.
 * GPU events are consumed after the existing reduction stream wait. The fixed
 * pool is allocated at setup; no allocation or added synchronization in a step.
 * CPU exchange/ready intervals overlap GPU work and must not be added to it. */
enum { BM_PROFILE_CLEAR, BM_PROFILE_OPERATOR, BM_PROFILE_NODES,
       BM_PROFILE_LOCAL, BM_PROFILE_FINISH, BM_PROFILE_DOT,
       BM_PROFILE_DIRECT, BM_PROFILE_IO, BM_PROFILE_EXCHANGE, BM_PROFILE_READY, BM_PROFILE_ALLREDUCE,
       BM_PROFILE_CATEGORIES, BM_PROFILE_EVENTS=256 };
static struct {
  bool enabled=false;
  int used=0, field[BM_PROFILE_EVENTS]={}, kind[BM_PROFILE_EVENTS]={};
  cudaEvent_t begin[BM_PROFILE_EVENTS]={}, end[BM_PROFILE_EVENTS]={};
  double seconds[3][BM_PROFILE_CATEGORIES]={}, calls[3][BM_PROFILE_CATEGORIES]={};
} g_biotMassProfile;

static void biot_mass_profile_initialize(mysolver_t*) {
  if (!Param.thePoroelasticMassMatrix) return;
  const char* value=getenv("HERCULES_BIOT_MASS_PROFILE");
  if (value && strcmp(value,"0") && strcmp(value,"1"))
    solver_abort(__FUNCTION_NAME,NULL,"HERCULES_BIOT_MASS_PROFILE must be 0 or 1\n");
  int choice=value && !strcmp(value,"1"),lo,hi;
  MPI_Allreduce(&choice,&lo,1,MPI_INT,MPI_MIN,comm_solver);
  MPI_Allreduce(&choice,&hi,1,MPI_INT,MPI_MAX,comm_solver);
  if (lo!=hi) solver_abort(__FUNCTION_NAME,NULL,"Biot mass profile settings differ between ranks\n");
  auto& p=g_biotMassProfile;p.enabled=choice;
  if (!p.enabled) return;
  for (int i=0;i<BM_PROFILE_EVENTS;++i) {
    solver_pore_pressure_device_check(cudaEventCreate(&p.begin[i]),"Mass profile start event");
    solver_pore_pressure_device_check(cudaEventCreate(&p.end[i]),"Mass profile end event");
  }
  if (Global.myID==0) monitor_print("Biot mass profiling: GPU event intervals and CPU exchange/reduction waits; "
      "no extra step synchronization. CPU/GPU intervals overlap; do not sum them. "
      "HERCULES_GPU_MPI_PROFILE=1 additionally separates packet/helper counters.\n");
}

static int biot_mass_profile_begin(mysolver_t* s,int scalar,int kind) {
  auto& p=g_biotMassProfile;
  if (!p.enabled) return -1;
  if (p.used==BM_PROFILE_EVENTS)
    solver_abort(__FUNCTION_NAME,NULL,"Biot mass profile event pool exhausted before a reduction wait\n");
  const int id=p.used++;p.field[id]=scalar;p.kind[id]=kind;
  solver_pore_pressure_device_check(cudaEventRecord(p.begin[id],s->streams[CUDA_STREAM_MAIN]),"Mass profile begin");
  return id;
}
static void biot_mass_profile_end(mysolver_t* s,int id) {
  if (id>=0) solver_pore_pressure_device_check(
      cudaEventRecord(g_biotMassProfile.end[id],s->streams[CUDA_STREAM_MAIN]),"Mass profile end");
}
// Call only after the stream is already known complete.
static void biot_mass_profile_collect() {
  auto& p=g_biotMassProfile;
  for (int i=0;i<p.used;++i) {
    float ms=0;
    solver_pore_pressure_device_check(cudaEventElapsedTime(&ms,p.begin[i],p.end[i]),"Mass profile elapsed");
    p.seconds[p.field[i]][p.kind[i]]+=ms*1e-3;
    p.calls[p.field[i]][p.kind[i]]+=1;
  }
  p.used=0;
}
static void biot_mass_profile_cpu(int scalar,int kind,double start) {
  if (!g_biotMassProfile.enabled) return;
  g_biotMassProfile.seconds[scalar][kind]+=MPI_Wtime()-start;
  g_biotMassProfile.calls[scalar][kind]+=1;
}

static void biot_mass_profile_finish() {
  auto& p=g_biotMassProfile;
  if (!p.enabled) return;
  // Teardown only: collect any final history/RHS copies after the last solve.
  solver_pore_pressure_device_check(cudaStreamSynchronize(
      Global.mySolver->streams[CUDA_STREAM_MAIN]),"Mass profile final wait");
  biot_mass_profile_collect();
  double sum[3][BM_PROFILE_CATEGORIES],lo[3][BM_PROFILE_CATEGORIES],hi[3][BM_PROFILE_CATEGORIES];
  double counts[3][BM_PROFILE_CATEGORIES];
  const int n=3*BM_PROFILE_CATEGORIES;
  MPI_Reduce(p.seconds,sum,n,MPI_DOUBLE,MPI_SUM,0,comm_solver);
  MPI_Reduce(p.seconds,lo,n,MPI_DOUBLE,MPI_MIN,0,comm_solver);
  MPI_Reduce(p.seconds,hi,n,MPI_DOUBLE,MPI_MAX,0,comm_solver);
  MPI_Reduce(p.calls,counts,n,MPI_DOUBLE,MPI_SUM,0,comm_solver);
  const char* names[]={"clear/copy GPU","operator GPU","separate nodes GPU",
      "fused local GPU","fused finish GPU","dot GPU","direct blocks GPU","scalar gather/scatter GPU","exchange CPU",
      "reduction ready CPU","Allreduce CPU"};
  if (Global.myID==0) {
    monitor_print("Biot mass profile: cumulative seconds average/min/max and calls/rank; "
        "CPU intervals include pending GPU work/rank imbalance and overlap GPU intervals.\n");
    for (int f=0;f<3;++f) for (int k=0;k<BM_PROFILE_CATEGORIES;++k)
      if (counts[f][k]) monitor_print("Biot mass profile %s %s: %.9g / %.9g / %.9g; calls/rank %.9g\n",
          f==2?"joint":(f?"storage":"inertia"),names[k],sum[f][k]/Global.theGroupSize,
          lo[f][k],hi[f][k],counts[f][k]/Global.theGroupSize);
  }
  for (int i=0;i<BM_PROFILE_EVENTS;++i) {cudaEventDestroy(p.begin[i]);cudaEventDestroy(p.end[i]);}
  p.enabled=false;
}
#endif
