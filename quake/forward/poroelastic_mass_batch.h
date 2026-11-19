#ifndef HERCULES_POROELASTIC_MASS_BATCH_H
#define HERCULES_POROELASTIC_MASS_BATCH_H

/* Two independent standard-PCG recurrences. Only scheduling, packets and
 * reductions are shared; coefficients, stopping tests, restarts and histories
 * belong to each field. A converged field leaves the batch immediately. */
static void biot_mass_batch_reduce(mysolver_t* s,double value[2]) {
  double local[8]={},global[8]={};
  double started=g_biotMassProfile.enabled?MPI_Wtime():0;
  const auto stream=s->streams[CUDA_STREAM_MAIN];
  solver_pore_pressure_device_check(cudaMemcpyAsync(local,g_biotMass.dots,
      sizeof(local),cudaMemcpyDeviceToHost,stream),"Batched Biot residual download");
  solver_pore_pressure_device_check(cudaStreamSynchronize(stream),"Batched Biot residual wait");
  biot_mass_profile_cpu(2,BM_PROFILE_READY,started);biot_mass_profile_collect();
  started=g_biotMassProfile.enabled?MPI_Wtime():0;
  MPI_Allreduce(local,global,8,MPI_DOUBLE,MPI_SUM,comm_solver);
  biot_mass_profile_cpu(2,BM_PROFILE_ALLREDUCE,started);
  value[0]=global[0];value[1]=global[4];
}

static void biot_mass_batch_apply(mysolver_t* s,int mask,double drag_dt,
    solver_float* in,solver_float* out,int patch=0,
    solver_float* preconditioned=nullptr,const solver_float* left=nullptr,int dot=0) {
  if(mask!=3) {
    const int scalar=mask==2;
    biot_mass_apply(s,scalar,scalar?0:drag_dt,in,out,patch,preconditioned,left,dot);
    return;
  }
  auto& m=g_biotMass;const auto stream=s->streams[CUDA_STREAM_MAIN];
  double started=g_biotMassProfile.enabled?MPI_Wtime():0;
  if(g_biotMassOneRound)biot_mass_one_round_exchange(s,in,2);
  else biot_mass_comm_exchange_pair(s,in,false);
  biot_mass_profile_cpu(2,BM_PROFILE_EXCHANGE,started);
  biot_mass_zero(s,2,out,7*(size_t)m.nodes*sizeof(solver_float));
  const int id=biot_mass_profile_begin(s,2,BM_PROFILE_OPERATOR);
  psolve_biot_mass_operator(stream,s->gpuDataDevice,m,2,patch,in,out);
  biot_mass_profile_end(s,id);
  started=g_biotMassProfile.enabled?MPI_Wtime():0;
  if(g_biotMassOneRound)biot_mass_one_round_used(stream);
  else biot_mass_comm_exchange_pair(s,out,true);
  biot_mass_profile_cpu(2,BM_PROFILE_EXCHANGE,started);
  for(int scalar=0;scalar<2;++scalar) {
    if(g_biotMassFusion && (preconditioned || dot)) {
      const int fid=biot_mass_profile_begin(s,scalar,BM_PROFILE_FINISH);
      psolve_biot_mass_finish_product(stream,s->gpuDataDevice,biot_mass_view(scalar),
          scalar,scalar?0:drag_dt,in,out,preconditioned?preconditioned:out,
          preconditioned!=nullptr,left,dot);
      biot_mass_profile_end(s,fid);
    } else biot_mass_nodes(s,scalar,0,scalar?0:drag_dt,0,in,out);
  }
}

static void biot_mass_batch_local(mysolver_t* s,int mask,double drag_dt,
    int mode,const double alpha[2]) {
  auto& m=g_biotMass;
  for(int scalar=0;scalar<2;++scalar) if(mask&(1<<scalar)) {
    const double d=scalar?0:drag_dt;
    if(g_biotMassFusion) {
      const int id=biot_mass_profile_begin(s,scalar,BM_PROFILE_LOCAL);
      psolve_biot_mass_local(s->streams[CUDA_STREAM_MAIN],s->gpuDataDevice,
          biot_mass_view(scalar),scalar,mode,alpha[scalar],d,m.r,m.z,
          !Param.thePoroelasticMassPreconditioner);
      biot_mass_profile_end(s,id);
    } else {
      if(mode==1) {
        biot_mass_nodes(s,scalar,3,alpha[scalar],1,m.direction,m.x);
        biot_mass_nodes(s,scalar,3,-alpha[scalar],1,m.product,m.r);
      } else biot_mass_nodes(s,scalar,4,0,0,m.product,m.r);
      biot_mass_nodes(s,scalar,2,d,0,m.r,m.z);
    }
  }
  if(Param.thePoroelasticMassPreconditioner) {
    biot_mass_batch_apply(s,mask,drag_dt,m.z,m.precondition,1,m.z,m.r,g_biotMassFusion?1:0);
    if(!g_biotMassFusion) for(int scalar=0;scalar<2;++scalar) if(mask&(1<<scalar))
      biot_mass_nodes(s,scalar,2,scalar?0:drag_dt,0,m.precondition,m.z);
  }
  if(!g_biotMassFusion) for(int scalar=0;scalar<2;++scalar) if(mask&(1<<scalar)) {
    const int id=biot_mass_profile_begin(s,scalar,BM_PROFILE_DOT);
    psolve_biot_mass_dot(s->streams[CUDA_STREAM_MAIN],biot_mass_view(scalar),scalar,m.r,m.z);
    biot_mass_profile_end(s,id);
  }
}

static void biot_mass_solve_batch(mysolver_t* s,double drag_dt) {
  Timer_Start("Biot batched inertia/stabilization PCG");
  auto& m=g_biotMass;
  biot_mass_solve_state_t c[2]={
      biot_mass_start(s,0,drag_dt,0,g_biotMassCompact?0:-1,false),
      biot_mass_start(s,1,0,3,g_biotMassCompact?0:-1,false)};
  for(int scalar=0;scalar<2;++scalar)
    biot_mass_nodes(s,scalar,1,0,0,m.z,m.direction);
  for(;;) {
    int mask=0;
    for(int f=0;f<2;++f) if(!c[f].converged &&
        c[f].iterations<Param.thePoroelasticMassMaxIterations) mask|=1<<f;
    if(!mask)break;
    biot_mass_batch_apply(s,mask,drag_dt,m.direction,m.product,0,nullptr,m.direction,g_biotMassFusion?1:0);
    if(!g_biotMassFusion) for(int f=0;f<2;++f) if(mask&(1<<f))
      psolve_biot_mass_dot(s->streams[CUDA_STREAM_MAIN],biot_mass_view(f),f,m.direction,m.product);
    double pap[2],alpha[2]={},next[2];biot_mass_batch_reduce(s,pap);
    for(int f=0;f<2;++f) if(mask&(1<<f)) {
      if(!(pap[f]>0) || !isfinite(pap[f]))
        solver_abort(__FUNCTION_NAME,NULL,"Batched Biot mass operator lost positive definiteness\n");
      alpha[f]=c[f].rz/pap[f];
    }
    biot_mass_batch_local(s,mask,drag_dt,1,alpha);
    biot_mass_batch_reduce(s,next);
    int verify=0;
    for(int f=0;f<2;++f) if(mask&(1<<f)) {
      ++c[f].iterations;
      if(!isfinite(next[f]) || next[f]<0)
        solver_abort(__FUNCTION_NAME,NULL,"Invalid batched Biot residual\n");
      if(next[f]<=c[f].target || c[f].iterations==Param.thePoroelasticMassMaxIterations)
        verify|=1<<f;
    }
    if(verify) {
      biot_mass_batch_apply(s,verify,drag_dt,m.x,m.product);
      biot_mass_batch_local(s,verify,drag_dt,2,alpha);
      double fresh[2];biot_mass_batch_reduce(s,fresh);
      for(int f=0;f<2;++f) if(verify&(1<<f)) {
        next[f]=fresh[f];c[f].converged=isfinite(next[f]) && next[f]>=0 && next[f]<=c[f].target;
      }
    }
    for(int f=0;f<2;++f) if(mask&(1<<f)) {
      if(verify&(1<<f)) biot_mass_nodes(s,f,1,0,0,m.z,m.direction);
      else biot_mass_nodes(s,f,3,1,next[f]/c[f].rz,m.z,m.direction);
      c[f].rz=next[f];
    }
  }
  for(int f=0;f<2;++f) {
    if(!c[f].converged) {
      c[f].rz=biot_mass_true_residual(s,f,f?0:drag_dt);
      c[f].converged=isfinite(c[f].rz) && c[f].rz>=0 && c[f].rz<=c[f].target;
    }
    biot_mass_accept(s,c[f],false);
  }
  Timer_Stop("Biot batched inertia/stabilization PCG");
}
#endif
