#ifndef HERCULES_POROELASTIC_MASS_RUNTIME_H
#define HERCULES_POROELASTIC_MASS_RUNTIME_H

#include "poroelastic_mass_profile.h"
#include "poroelastic_mass_comm.h"

// Separate reduction slots let scalar and mechanics kernels be interleaved.
static biot_mass_t biot_mass_view(int scalar) {
  biot_mass_t m=g_biotMass;m.dots+=4*scalar;return m;
}
static void biot_stabilization_prepare(mysolver_t*,mesh_t*);

/* Reuse the established field-specific prolongation and its transpose, along
 * with both communication backends. Sequential exchanges use completed tags;
 * all ranks participate, including ranks with no local poroelastic cells. */
static void biot_mass_exchange(mysolver_t* s,solver_float* value,int scalar,bool rhs) {
  const double started=g_biotMassProfile.enabled?MPI_Wtime():0;
  if (g_biotMassCompact || g_biotMassHaloPack || (!scalar && g_biotMassHaloBatch)) {
    biot_mass_comm_exchange(s,value,scalar,rhs);
    biot_mass_profile_cpu(scalar,BM_PROFILE_EXCHANGE,started);return;
  }
  const int saved_scope=g_gpuMpiMassScope;
  g_gpuMpiMassScope=scalar?11:10;
  const int N=g_biotMass.nodes;
  schedule_t* first=rhs?s->dn_sched:s->an_sched;
  schedule_t* last=rhs?s->an_sched:s->dn_sched;
  const int direction=rhs?CONTRIBUTION:SHARING;
  const int adjustment=rhs?DISTRIBUTION:ASSIGNMENT;
  for (int field=scalar?2:0;field<(scalar?3:2);++field) {
    solver_float* device=value+3*(size_t)N*field;
    const int tag=scalar?(rhs?DN_PORE_PRESSURE_FORCE_MSG:AN_PORE_PRESSURE_MSG):(rhs?180:182);
    const int tag2=scalar?(rhs?AN_PORE_PRESSURE_FORCE_MSG:DN_PORE_PRESSURE_MSG):tag+1;
    if (gpu_device_sync_enabled()) {
      if (scalar) {
        solver_scalar_device_schedule_senddata(s,first,device,direction,tag);
        solver_scalar_device_adjust(s,device,adjustment,ACOUSTIC_CONSTRAINT_PORE);
        solver_scalar_device_schedule_senddata(s,last,device,direction,tag2);
      } else {
        fvector_t* vector=(fvector_t*)device;
        solver_vector_device_schedule_senddata(s,first,vector,direction,tag);
        solver_vector_device_adjust(s,vector,adjustment);
        solver_vector_device_schedule_senddata(s,last,vector,direction,tag2);
      }
    } else {
      const size_t bytes=(scalar?1:3)*(size_t)N*sizeof(solver_float);
      solver_float* host=g_biotMassHost.data();
      solver_pore_pressure_device_check(cudaMemcpyAsync(host,device,bytes,
          cudaMemcpyDeviceToHost,s->streams[CUDA_STREAM_MAIN]),"Mass vector download");
      solver_pore_pressure_device_check(cudaStreamSynchronize(s->streams[CUDA_STREAM_MAIN]),"Mass vector wait");
      schedule_senddata(first,host,scalar?1:3,direction,tag);
      compute_adjust(host,scalar?1:3,adjustment,scalar?ACOUSTIC_CONSTRAINT_PORE:ACOUSTIC_CONSTRAINT_SOLID);
      schedule_senddata(last,host,scalar?1:3,direction,tag2);
      solver_pore_pressure_device_check(cudaMemcpyAsync(device,host,bytes,
          cudaMemcpyHostToDevice,s->streams[CUDA_STREAM_MAIN]),"Mass vector upload");
    }
  }
  g_gpuMpiMassScope=saved_scope;
  biot_mass_profile_cpu(scalar,BM_PROFILE_EXCHANGE,started);
}

static void biot_mass_nodes(mysolver_t* s,int scalar,int op,double a,double b,
    const solver_float* in,solver_float* out) {
  const int id=biot_mass_profile_begin(s,scalar,BM_PROFILE_NODES);
  psolve_biot_mass_nodes(s->streams[CUDA_STREAM_MAIN],s->gpuDataDevice,biot_mass_view(scalar),
      scalar,op,a,b,in,out);
  biot_mass_profile_end(s,id);
}
static void biot_mass_zero(mysolver_t* s,int scalar,solver_float* out,size_t bytes) {
  const int id=biot_mass_profile_begin(s,scalar,BM_PROFILE_CLEAR);
  solver_pore_pressure_device_check(cudaMemsetAsync(out,0,bytes,s->streams[CUDA_STREAM_MAIN]),"Mass vector clear");
  biot_mass_profile_end(s,id);
}
static void biot_mass_copy(mysolver_t* s,int scalar,solver_float* out,const solver_float* in,size_t bytes) {
  const int id=biot_mass_profile_begin(s,scalar,BM_PROFILE_CLEAR);
  solver_pore_pressure_device_check(cudaMemcpyAsync(out,in,bytes,cudaMemcpyDeviceToDevice,
      s->streams[CUDA_STREAM_MAIN]),"Mass vector copy");
  biot_mass_profile_end(s,id);
}
static void biot_mass_apply(mysolver_t* s,int scalar,double drag_dt,
    solver_float* in,solver_float* out,int patch=0,
    solver_float* preconditioned=nullptr,const solver_float* left=nullptr,int dot_mode=0) {
  const biot_mass_t m=biot_mass_view(scalar);
  const cudaStream_t stream=s->streams[CUDA_STREAM_MAIN];
  if(g_biotMassOneRound) {
    const double started=g_biotMassProfile.enabled?MPI_Wtime():0;
    biot_mass_one_round_exchange(s,in,scalar);
    biot_mass_profile_cpu(scalar,BM_PROFILE_EXCHANGE,started);
  } else biot_mass_exchange(s,in,scalar,false);
  biot_mass_zero(s,scalar,out+(scalar?6*(size_t)m.nodes:0),
      (scalar?1:6)*(size_t)m.nodes*sizeof(solver_float));
  int id=biot_mass_profile_begin(s,scalar,BM_PROFILE_OPERATOR);
  psolve_biot_mass_operator(stream,s->gpuDataDevice,m,scalar,patch,in,out);
  biot_mass_profile_end(s,id);
  if(g_biotMassOneRound) biot_mass_one_round_used(stream);
  else biot_mass_exchange(s,out,scalar,true);
  // Dot products and nodal inverses must follow MPI accumulation/projection.
  if (g_biotMassFusion && (preconditioned || dot_mode)) {
    id=biot_mass_profile_begin(s,scalar,BM_PROFILE_FINISH);
    psolve_biot_mass_finish_product(stream,s->gpuDataDevice,m,scalar,drag_dt,
        in,out,preconditioned?preconditioned:out,preconditioned!=nullptr,left,dot_mode);
    biot_mass_profile_end(s,id);
  } else biot_mass_nodes(s,scalar,0,drag_dt,0,in,out);
}
static void biot_mass_reduce(mysolver_t* s,int scalar,int count,double* global) {
  double local[2]={};
  const cudaStream_t stream=s->streams[CUDA_STREAM_MAIN];
  double started=g_biotMassProfile.enabled?MPI_Wtime():0;
  solver_pore_pressure_device_check(cudaMemcpyAsync(local,g_biotMass.dots+4*scalar,count*sizeof(double),
      cudaMemcpyDeviceToHost,stream),"Mass residual reduction");
  solver_pore_pressure_device_check(cudaStreamSynchronize(stream),"Mass residual wait");
  biot_mass_profile_cpu(scalar,BM_PROFILE_READY,started);
  biot_mass_profile_collect();
  started=g_biotMassProfile.enabled?MPI_Wtime():0;
  MPI_Allreduce(local,global,count,MPI_DOUBLE,MPI_SUM,comm_solver);
  biot_mass_profile_cpu(scalar,BM_PROFILE_ALLREDUCE,started);
}
static double biot_mass_dot(mysolver_t* s,int scalar,
    const solver_float* a,const solver_float* b) {
  const int id=biot_mass_profile_begin(s,scalar,BM_PROFILE_DOT);
  psolve_biot_mass_dot(s->streams[CUDA_STREAM_MAIN],biot_mass_view(scalar),scalar,a,b);
  biot_mass_profile_end(s,id);
  double result;biot_mass_reduce(s,scalar,1,&result);return result;
}

/* Unfused comparison path: B = L^-1 H L^-1. Fused paths below keep these same
 * assembly stages and projections, combining only neighboring local passes. */
static void biot_mass_precondition(mysolver_t* s,int scalar,double drag_dt,
    solver_float* in,solver_float* out) {
  biot_mass_nodes(s,scalar,2,drag_dt,0,in,out);
  if (Param.thePoroelasticMassPreconditioner) {
    biot_mass_apply(s,scalar,drag_dt,out,g_biotMass.precondition,1);
    biot_mass_nodes(s,scalar,2,drag_dt,0,g_biotMass.precondition,out);
  }
}
static void biot_mass_local_stage(mysolver_t* s,int scalar,double drag_dt,
    int mode,double alpha,bool norm) {
  const auto m=biot_mass_view(scalar);
  if (g_biotMassFusion) {
    const int id=biot_mass_profile_begin(s,scalar,BM_PROFILE_LOCAL);
    psolve_biot_mass_local(s->streams[CUDA_STREAM_MAIN],s->gpuDataDevice,m,
        scalar,mode,alpha,drag_dt,m.r,m.z,norm && !Param.thePoroelasticMassPreconditioner);
    biot_mass_profile_end(s,id);
    if (Param.thePoroelasticMassPreconditioner)
      biot_mass_apply(s,scalar,drag_dt,m.z,m.precondition,1,m.z,norm?m.r:nullptr,norm?1:0);
  } else {
    if (mode==1) {
      biot_mass_nodes(s,scalar,3,alpha,1,m.direction,m.x);
      biot_mass_nodes(s,scalar,3,-alpha,1,m.product,m.r);
    } else if (mode==2) biot_mass_nodes(s,scalar,4,0,0,m.product,m.r);
    else if (mode==3) biot_mass_nodes(s,scalar,1,0,0,m.rhs,m.r);
    biot_mass_precondition(s,scalar,drag_dt,m.r,m.z);
  }
}
static double biot_mass_norm(mysolver_t* s,int scalar,int diagonal=-1,double* active=nullptr) {
  if (!g_biotMassFusion) {
    const int id=biot_mass_profile_begin(s,scalar,BM_PROFILE_DOT);
    psolve_biot_mass_dot(s->streams[CUDA_STREAM_MAIN],biot_mass_view(scalar),scalar,g_biotMass.r,g_biotMass.z);
    biot_mass_profile_end(s,id);
  }
  const bool combined=g_biotMassCompact && diagonal>=0;
  if (combined) solver_pore_pressure_device_check(cudaMemcpyAsync(g_biotMass.dots+4*scalar+1,
      g_biotMass.dots+4*scalar+2+diagonal,sizeof(double),cudaMemcpyDeviceToDevice,
      s->streams[CUDA_STREAM_MAIN]),"Biot direct-block normalization");
  double result[2]={};biot_mass_reduce(s,scalar,combined?2:1,result);
  if (active)*active=result[0];
  return result[0]+result[1];
}
static double biot_mass_true_residual(mysolver_t* s,int scalar,double drag_dt) {
  auto& m=g_biotMass;
  biot_mass_apply(s,scalar,drag_dt,m.x,m.product);
  biot_mass_local_stage(s,scalar,drag_dt,2,0,true);
  return biot_mass_norm(s,scalar);
}
static void biot_mass_pair(mysolver_t* s,int scalar,double pair[2]) {
  const auto m=biot_mass_view(scalar);
  if (!g_biotMassFusion) {
    const int id=biot_mass_profile_begin(s,scalar,BM_PROFILE_DOT);
    psolve_biot_mass_dot_pair(s->streams[CUDA_STREAM_MAIN],m,scalar,m.r,m.z,m.applied);
    biot_mass_profile_end(s,id);
  }
  biot_mass_reduce(s,scalar,2,pair);
}

/* Both CG recurrences require a recomputed true residual before acceptance.
 * Normalization always uses the projected RHS, even with a previous guess.
 * Separate histories prevent mixing inertia stages and the two storage RHSs. */
struct biot_mass_solve_state_t {
  int scalar,history,diagonal,iterations;
  double drag_dt,initial,target,rz;
  bool converged,zero_start;
};
static biot_mass_solve_state_t biot_mass_start(mysolver_t* s,int scalar,
    double drag_dt,int history,int diagonal,bool timed) {
  biot_mass_t& m=g_biotMass;
  const size_t offset=scalar?6*(size_t)m.nodes:0;
  const size_t bytes=(scalar?1:6)*(size_t)m.nodes*sizeof(solver_float);
  if(timed) Timer_Start(scalar?"Biot storage PCG":"Biot inertia PCG");
  biot_mass_local_stage(s,scalar,drag_dt,3,0,true);
  double active_initial;
  const double initial=biot_mass_norm(s,scalar,diagonal,&active_initial);
  if (!isfinite(initial) || initial<0)
    solver_abort(__FUNCTION_NAME,NULL,"Invalid Biot mass RHS norm\n");
  const double target=initial*Param.thePoroelasticMassTolerance*Param.thePoroelasticMassTolerance;
  const bool keep=history>=0 && Param.thePoroelasticMassInitialGuess;
  const bool cache=keep && g_biotMassPredictorCache;
  const int depth=Param.thePoroelasticMassInitialGuess;
  // Material, topology, projections and B are frozen by setup. drag_dt is
  // the only per-call operator coefficient; never reuse products across it.
  if(cache && g_biotMassProductDragDt[history]!=drag_dt) {
    g_biotMassProductHistoryCount[history]=0;
    g_biotMassProductDragDt[history]=drag_dt;
  }
  const int order=keep?std::min(depth,g_biotMassHistoryCount[history]):0;
  const bool cached_start=cache && active_initial>0 && order>0 &&
      g_biotMassProductHistoryCount[history]>=order;
  if(cache && !cached_start)
    biot_mass_copy(s,scalar,m.history_preconditioned_product[history][depth-1],m.z+offset,bytes);
  double rz=active_initial;
  bool zero_start=true;
  if (active_initial>0 && keep && g_biotMassHistoryCount[history]) {
    zero_start=false;
    if(cached_start) {
      const int id=biot_mass_profile_begin(s,scalar,BM_PROFILE_NODES);
      psolve_biot_mass_cached_predict(s->streams[CUDA_STREAM_MAIN],s->gpuDataDevice,
          biot_mass_view(scalar),scalar,history,order,depth);
      biot_mass_profile_end(s,id);
      biot_mass_reduce(s,scalar,1,&rz);
      ++g_biotMassCachedPredictions[scalar];
      // A cached residual is an initial approximation, never a certificate.
      // Refresh before zero-iteration acceptance or rejecting a suspect guess.
      if(!isfinite(rz) || rz<0 || rz<=target || rz>active_initial) {
        rz=biot_mass_true_residual(s,scalar,drag_dt);
        ++g_biotMassCacheRefreshes[scalar];
      }
    } else if (order==1) biot_mass_copy(s,scalar,m.x+offset,m.history[history][0],bytes);
    else {
      const int id=biot_mass_profile_begin(s,scalar,BM_PROFILE_NODES);
      psolve_biot_mass_predict(s->streams[CUDA_STREAM_MAIN],m,scalar,history,order);
      biot_mass_profile_end(s,id);
    }
    if(!cached_start) rz=biot_mass_true_residual(s,scalar,drag_dt);
    // Reject a poor extrapolation using the existing global norm. This may
    // occur at source onset; it must never weaken the acceptance criterion.
    if (order>1 && (!isfinite(rz) || rz<0 || rz>active_initial)) {
      biot_mass_zero(s,scalar,m.x+offset,bytes);
      biot_mass_local_stage(s,scalar,drag_dt,3,0,true);
      rz=biot_mass_norm(s,scalar);
      zero_start=true;
    }
  } else biot_mass_zero(s,scalar,m.x+offset,bytes);
  bool converged=isfinite(rz) && rz>=0 && rz<=target;
  return {scalar,history,diagonal,0,drag_dt,initial,target,rz,converged,zero_start};
}
static void biot_mass_accept(mysolver_t* s,biot_mass_solve_state_t& c,bool timed) {
  auto& m=g_biotMass;
  const int scalar=c.scalar,history=c.history;
  const double initial=c.initial;
  double& rz=c.rz;int& iterations=c.iterations;bool& converged=c.converged;
  const bool keep=history>=0 && Param.thePoroelasticMassInitialGuess;
  const bool cache=keep && g_biotMassPredictorCache,zero_start=c.zero_start;
  const int depth=Param.thePoroelasticMassInitialGuess;
  const size_t offset=scalar?6*(size_t)m.nodes:0;
  const size_t bytes=(scalar?1:6)*(size_t)m.nodes*sizeof(solver_float);
  const double relative=initial>0?sqrt(rz/initial):0;
  if (!converged) solver_abort(__FUNCTION_NAME,NULL,
      "Biot %s PCG did not converge in %d iterations (true relative residual %.9g, tolerance %.9g).\n",
      scalar?"storage":"inertia",iterations,relative,Param.thePoroelasticMassTolerance);
  if (keep) {
    auto& h=m.history[history];
    solver_float* oldest=h[Param.thePoroelasticMassInitialGuess-1];
    for (int k=Param.thePoroelasticMassInitialGuess-1;k>0;--k)h[k]=h[k-1];
    h[0]=oldest;
    if(cache) {
      auto& a=m.history_product[history];
      auto& ba=m.history_preconditioned_product[history];
      solver_float* oldest_a=a[depth-1];solver_float* oldest_ba=ba[depth-1];
      for(int k=depth-1;k>0;--k) {a[k]=a[k-1];ba[k]=ba[k-1];}
      a[0]=oldest_a;ba[0]=oldest_ba;
      const int id=biot_mass_profile_begin(s,scalar,BM_PROFILE_NODES);
      psolve_biot_mass_cache_store(s->streams[CUDA_STREAM_MAIN],m,scalar,history,
          zero_start && iterations==0);
      biot_mass_profile_end(s,id);
      g_biotMassProductHistoryCount[history]=std::min(3,g_biotMassProductHistoryCount[history]+1);
    } else biot_mass_copy(s,scalar,h[0],m.x+offset,bytes);
    g_biotMassHistoryCount[history]=std::min(3,g_biotMassHistoryCount[history]+1);
  }
  ++g_biotMassSolves[scalar];g_biotMassIterations[scalar]+=iterations;
  g_biotMassMaximum[scalar]=std::max(g_biotMassMaximum[scalar],iterations);
  g_biotMassWorstResidual[scalar]=fmax(g_biotMassWorstResidual[scalar],relative);
  if(timed) Timer_Stop(scalar?"Biot storage PCG":"Biot inertia PCG");
}
static void biot_mass_solve(mysolver_t* s,int scalar,double drag_dt,int history=-1,int diagonal=-1) {
  auto c=biot_mass_start(s,scalar,drag_dt,history,diagonal,true);
  auto& m=g_biotMass;
  const double target=c.target;
  double& rz=c.rz;int& iterations=c.iterations;bool& converged=c.converged;
  if (Param.thePoroelasticMassKrylov) {
    // Chronopoulos-Gear recurrence: one combined (r,z),(z,Az) reduction.
    // Residual replacement also resets both direction recurrences.
    double previous=0,previous_alpha=1;
    bool restart=true;
    while (!converged && iterations<Param.thePoroelasticMassMaxIterations) {
      biot_mass_apply(s,scalar,drag_dt,m.z,m.applied,0,nullptr,m.r,g_biotMassFusion?2:0);
      double pair[2];biot_mass_pair(s,scalar,pair);rz=pair[0];
      if (!isfinite(rz) || rz<0)
        solver_abort(__FUNCTION_NAME,NULL,"Invalid single-reduction mass residual\n");
      if (rz<=target) {
        rz=biot_mass_true_residual(s,scalar,drag_dt);
        converged=isfinite(rz) && rz>=0 && rz<=target;
        restart=true;
        if (converged) break;
        continue;
      }
      const double beta=restart?0:rz/previous;
      const double denominator=pair[1]-(restart?0:beta*rz/previous_alpha);
      if (!(denominator>0) || !isfinite(denominator)) {
        if (restart) solver_abort(__FUNCTION_NAME,NULL,"Biot mass operator lost positive definiteness\n");
        rz=biot_mass_true_residual(s,scalar,drag_dt);restart=true;continue;
      }
      if (restart) {
        biot_mass_nodes(s,scalar,1,0,0,m.z,m.direction);
        biot_mass_nodes(s,scalar,1,0,0,m.applied,m.product);
      } else {
        biot_mass_nodes(s,scalar,3,1,beta,m.z,m.direction);
        biot_mass_nodes(s,scalar,3,1,beta,m.applied,m.product);
      }
      const double alpha=rz/denominator;
      biot_mass_local_stage(s,scalar,drag_dt,1,alpha,false);
      previous=rz;previous_alpha=alpha;restart=false;++iterations;
    }
  } else {
    biot_mass_nodes(s,scalar,1,0,0,m.z,m.direction);
    while (!converged && iterations<Param.thePoroelasticMassMaxIterations) {
      biot_mass_apply(s,scalar,drag_dt,m.direction,m.product,0,nullptr,m.direction,g_biotMassFusion?1:0);
      double pap;
      if (g_biotMassFusion) biot_mass_reduce(s,scalar,1,&pap);
      else pap=biot_mass_dot(s,scalar,m.direction,m.product);
      if (!(pap>0) || !isfinite(pap))
        solver_abort(__FUNCTION_NAME,NULL,"Biot mass operator lost positive definiteness\n");
      const double alpha=rz/pap;
      biot_mass_local_stage(s,scalar,drag_dt,1,alpha,true);
      double next=biot_mass_norm(s,scalar);
      ++iterations;
      if (!isfinite(next) || next<0)
        solver_abort(__FUNCTION_NAME,NULL,"Invalid Biot mass residual\n");
      if (next<=target || iterations==Param.thePoroelasticMassMaxIterations) {
        next=biot_mass_true_residual(s,scalar,drag_dt);
        converged=isfinite(next) && next>=0 && next<=target;
        biot_mass_nodes(s,scalar,1,0,0,m.z,m.direction);
      } else biot_mass_nodes(s,scalar,3,1,next/rz,m.z,m.direction);
      rz=next;
    }
  }
  if (!converged) {
    rz=biot_mass_true_residual(s,scalar,drag_dt);
    converged=isfinite(rz) && rz>=0 && rz<=target;
  }
  biot_mass_accept(s,c,true);
}

#include "poroelastic_mass_batch.h"

static void biot_mass_displacement(mysolver_t* s,mesh_t* mesh) {
  biot_mass_t& m=g_biotMass;
  const cudaStream_t stream=s->streams[CUDA_STREAM_MAIN];
  const bool damping=g_biotMassHasSolidDamping ||
      (Param.thePoroelasticFluidViscosity>0 && !poroelastic_zero_permeability());
  const double gamma=1-1/sqrt(2.0),dt=Param.theDeltaT;
  if (g_biotMassCompact) {
    const int id=biot_mass_profile_begin(s,0,BM_PROFILE_DIRECT);
    psolve_biot_mass_diagonal(stream,s->gpuDataDevice,m,0,damping,nullptr);
    biot_mass_profile_end(s,id);
  }
  if (damping) {
    psolve_biot_mass_state(stream,s->gpuDataDevice,m,0,0);
    biot_mass_apply(s,0,0,m.x,m.rhs); // M times previous half-step increment
  } else biot_mass_zero(s,0,m.rhs,6*(size_t)m.nodes*sizeof(solver_float));
  psolve_biot_mass_state(stream,s->gpuDataDevice,m,1,damping?gamma:1);
  if(g_biotMassBatchSolves && biot_stabilization_active()) {
    if(g_biotMassStabilizationReady) solver_abort(__FUNCTION_NAME,NULL,"Unconsumed batched stabilization\n");
    biot_stabilization_prepare(s,mesh);
    solver_float* work=s->gpuData->biotStabilizationWorkDevice;
    if(g_biotMassCompact) {
      psolve_biot_mass_diagonal(stream,s->gpuDataDevice,biot_mass_view(1),1,0,work);
      psolve_biot_mass_scalar_io(stream,m,work,0);
    } else biot_mass_copy(s,1,m.rhs+6*(size_t)m.nodes,work,(size_t)m.nodes*sizeof(solver_float));
    biot_mass_solve_batch(s,damping?gamma*dt:0);
    if(g_biotMassCompact) psolve_biot_mass_scalar_io(stream,m,work,1);
    else biot_mass_copy(s,1,work,m.x+6*(size_t)m.nodes,(size_t)m.nodes*sizeof(solver_float));
    g_biotMassStabilizationReady=true;
  } else biot_mass_solve(s,0,damping?gamma*dt:0,0,g_biotMassCompact?0:-1);
  if (damping) {
    psolve_biot_mass_state(stream,s->gpuDataDevice,m,1,1-gamma);
    psolve_biot_mass_state(stream,s->gpuDataDevice,m,2,(1-gamma)*dt);
    biot_mass_solve(s,0,gamma*dt,1,g_biotMassCompact?1:-1);
  }
  psolve_biot_mass_state(stream,s->gpuDataDevice,m,damping?4:3,0);
}

static void biot_mass_pressure(mysolver_t* s,mesh_t*,solver_float* rhs,bool update) {
  biot_mass_t& m=g_biotMass;
  const size_t bytes=(size_t)m.nodes*sizeof(solver_float),offset=6*(size_t)m.nodes;
  const cudaStream_t stream=s->streams[CUDA_STREAM_MAIN];
  if (g_biotMassCompact) {
    int id=biot_mass_profile_begin(s,1,BM_PROFILE_DIRECT);
    psolve_biot_mass_diagonal(stream,s->gpuDataDevice,biot_mass_view(1),1,update,rhs);
    biot_mass_profile_end(s,id);
    id=biot_mass_profile_begin(s,1,BM_PROFILE_IO);
    psolve_biot_mass_scalar_io(stream,m,rhs,0);
    biot_mass_profile_end(s,id);
  } else biot_mass_copy(s,1,m.rhs+offset,rhs,bytes);
  // The physical pressure update uses C+dt*R/2. The stabilization auxiliary
  // retains C, preserving its original symmetric positive dissipation.
  const double surface=update && s->gpuData->biotDrainCapacityDevice?1:0;
  biot_mass_solve(s,1,surface,update?2:3,g_biotMassCompact?0:-1);
  if (update) psolve_biot_mass_state(stream,s->gpuDataDevice,m,5,0);
  else if (g_biotMassCompact) {
    const int id=biot_mass_profile_begin(s,1,BM_PROFILE_IO);
    psolve_biot_mass_scalar_io(stream,m,rhs,1);biot_mass_profile_end(s,id);
  } else biot_mass_copy(s,1,rhs,m.x+offset,bytes);
}

static void biot_mass_finish() {
  if (!Param.thePoroelasticMassMatrix) return;
  biot_mass_profile_finish();
  biot_mass_comm_finish();
  biot_mass_one_round_finish();
  if (Global.myID==0) for (int s=0;s<2;++s)
    monitor_print("Biot %s PCG: %lld solves, mean %.3f iterations, maximum %d, "
        "worst true relative residual %.9g.\n",s?"storage":"inertia",g_biotMassSolves[s],
        g_biotMassSolves[s]?(double)g_biotMassIterations[s]/g_biotMassSolves[s]:0,
        g_biotMassMaximum[s],g_biotMassWorstResidual[s]);
  if(Global.myID==0 && g_biotMassPredictorCache) for(int s=0;s<2;++s)
    monitor_print("Biot %s predictor cache: %lld starting residuals reused, "
        "%lld refreshed before initial acceptance/rejection; final true-residual checks retained.\n",
        s?"storage":"inertia",g_biotMassCachedPredictions[s],g_biotMassCacheRefreshes[s]);
  auto& m=g_biotMass;
  cudaFree(m.node);cudaFree(m.index);cudaFree(m.diagonal);
  cudaFree(m.element);cudaFree(m.boundary);cudaFree(m.owned);
  cudaFree(m.rho);cudaFree(m.remainder);cudaFree(m.solid_damping);
  cudaFree(m.x);cudaFree(m.r);cudaFree(m.z);cudaFree(m.direction);
  cudaFree(m.product);cudaFree(m.rhs);cudaFree(m.dots);
  cudaFree(m.precondition);cudaFree(m.applied);
  for (auto& h:m.history)for(auto ptr:h)cudaFree(ptr);
  for (auto& h:m.history_product)for(auto ptr:h)cudaFree(ptr);
  for (auto& h:m.history_preconditioned_product)for(auto ptr:h)cudaFree(ptr);
  cudaFree(m.row_node);cudaFree(m.neighbor);cudaFree(m.weight);
  cudaFree(m.row_mass_weight);cudaFree(m.row_storage_weight);
  m={};std::vector<solver_float>().swap(g_biotMassHost);
}
#endif
