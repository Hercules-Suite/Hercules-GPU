#ifndef HERCULES_POROELASTIC_MASS_COMM_H
#define HERCULES_POROELASTIC_MASS_COMM_H

/* Private consistent-mass packets reuse the established gather/scatter and
 * constraint kernels. Mechanics packets use [u vectors][w vectors]; joint
 * inertia/storage packets interleave seven components per node.
 * Host MPI stages only packets; CUDA-aware MPI never stages field data. */
struct biot_mass_packet_t {
  messenger_t* map=nullptr;
  solver_float* device=nullptr;
  std::vector<solver_float> host;
};
struct biot_mass_plan_t {
  std::vector<biot_mass_packet_t> send,recv;
  std::vector<MPI_Request> sr,rr;
  std::vector<MPI_Status> ss,rs;
  mpi_unpack_fence_t fence;
  biot_mass_packet_entry_t *gather=nullptr,*scatter=nullptr;
  int32_t* offsets=nullptr;
  int gather_count=0,scatter_nodes=0;
};
static schedule_t* g_biotMassSchedules[2]={};
static biot_mass_plan_t g_biotMassPlans[4];
static bool g_biotMassCommReady=false,g_biotMassCommDevice=false;

static void biot_mass_comm_initialize(mysolver_t* s) {
  if (!Param.thePoroelasticMassMatrix || (!g_biotMassCompact && !g_biotMassHaloBatch && !g_biotMassHaloPack && !g_biotMassBatchSolves)) return;
  std::vector<uint8_t> active(g_biotMass.full_nodes,1);
  if (g_biotMassCompact)
    for (int n=0;n<g_biotMass.full_nodes;++n)active[n]=g_biotMassIndex[n]>=0;
  // Verify IDs/order and owner/sharer agreement before removing any mapping.
  for (auto* sched:{s->dn_sched,s->an_sched}) {
    scalar_mpi_share_mask(sched,CONTRIBUTION,active);
    scalar_mpi_share_mask(sched,SHARING,active);
  }
  for (int i=0;i<2;++i) {
    auto* sched=g_biotMassSchedules[i]=scalar_mpi_filter(i?s->an_sched:s->dn_sched,active,0);
    for (auto* list:{sched->first_c,sched->first_s}) {
      for (auto* m=list;m;m=m->next) {
        if ((size_t)m->nodecount>INT_MAX/(7*sizeof(solver_float)))
          solver_abort(__FUNCTION_NAME,NULL,"Biot mass packet exceeds MPI byte count\n");
        if (g_biotMassCompact)for(int j=0;j<m->nodecount;++j) {
          int index=g_biotMassIndex[m->mapping[j]];
          if(index<0)solver_abort(__FUNCTION_NAME,NULL,"Inconsistent compact mass peer mask\n");
          m->mapping[j]=index;
        }
        std::vector<int32_t> sorted(m->mapping,m->mapping+m->nodecount);
        std::sort(sorted.begin(),sorted.end());
        m->mapping_has_duplicates=std::adjacent_find(sorted.begin(),sorted.end())!=sorted.end();
      }
      solver_init_device_messenger_mapping(list);
    }
  }
  std::vector<int32_t>().swap(g_biotMassIndex);
  g_biotMassCommDevice=gpu_device_sync_enabled();
  if (!g_biotMassCompact) solver_init_device_adjust_sync();
  for (int phase=0;phase<4;++phase) {
    auto& p=g_biotMassPlans[phase];auto* sched=g_biotMassSchedules[phase==1 || phase==2];
    for (int receiving=0;receiving<2;++receiving) {
      auto& packets=receiving?p.recv:p.send;
      const bool c=(phase<2)!=(receiving!=0);
      for(auto* map=c?sched->first_c:sched->first_s;map;map=map->next) {
        packets.emplace_back();auto& v=packets.back();v.map=map;
        const size_t count=(g_biotMassBatchSolves?7:6)*(size_t)map->nodecount;
        solver_pore_pressure_device_check(cudaMalloc((void**)&v.device,count*sizeof(solver_float)),"Biot mass packet");
        if(!g_biotMassCommDevice)v.host.resize(count);
      }
    }
    p.sr.resize(p.send.size());p.ss.resize(p.send.size());p.rr.resize(p.recv.size());p.rs.resize(p.recv.size());
    p.fence.init();
    if(g_biotMassHaloPack || g_biotMassBatchSolves) {
      for(int receiving=0;receiving<2;++receiving) {
        std::vector<biot_mass_packet_entry_t> entries;
        for(auto& v:receiving?p.recv:p.send) {
          if(entries.size()+(size_t)v.map->nodecount>INT_MAX)
            solver_abort(__FUNCTION_NAME,NULL,"Biot aggregate packet mapping exceeds index range\n");
          for(int j=0;j<v.map->nodecount;++j)
            entries.push_back({v.device,v.map->mapping[j],j,v.map->nodecount});
        }
        if(receiving) {
          std::stable_sort(entries.begin(),entries.end(),[](const biot_mass_packet_entry_t& a,
              const biot_mass_packet_entry_t& b){return a.node<b.node;});
          std::vector<int32_t> offsets;
          for(size_t j=0;j<entries.size();++j)
            if(j==0 || entries[j].node!=entries[j-1].node)offsets.push_back(j);
          p.scatter_nodes=offsets.size();offsets.push_back(entries.size());
          acoustic_constraint_upload(&p.offsets,offsets);
          acoustic_constraint_upload(&p.scatter,entries);
        } else {
          p.gather_count=entries.size();acoustic_constraint_upload(&p.gather,entries);
        }
      }
    }
  }
  g_biotMassCommReady=true;
  if(Global.myID==0)monitor_print("Biot mass halo: %s solid/fluid packets, %s transport, %s node mappings; "
      "fixed peer unpack order and receive-reuse events.\n",g_biotMassHaloBatch?"combined":"separate",
      g_biotMassCommDevice?"CUDA-aware":"host-staged",g_biotMassCompact?"compact":"full");
  if(Global.myID==0)monitor_print("Biot mass packet kernels: %s; stable per-node receive order.\n",
      g_biotMassHaloPack?"one gather/unpack launch per phase":"per peer/field (comparison path)");
}

static void biot_mass_comm_phase(mysolver_t* s,solver_float* values,int scalar,
    int phase,int field,int fields) {
  auto& p=g_biotMassPlans[phase];const int width=scalar==2?7:(scalar?1:3);
  const bool aggregate=g_biotMassHaloPack || scalar==2;
  const int tag=180+phase;
  const cudaStream_t stream=s->streams[CUDA_STREAM_MAIN];
  gpu_mpi_profile_t profile;
  if(g_gpuMpiProfile && (!p.send.empty() || !p.recv.empty())) {
    profile.metrics=g_gpuMpiMetrics[scalar?11:10];auto* v=profile.metrics;v[0]+=1;
    for(auto& t:p.send){v[1]+=(double)t.map->nodecount*width*fields*sizeof(solver_float);v[3]+=1;}
    for(auto& t:p.recv){v[2]+=(double)t.map->nodecount*width*fields*sizeof(solver_float);v[4]+=1;}
    profile.start=MPI_Wtime();
  }
  p.fence.wait();
  for(size_t i=0;i<p.recv.size();++i) {
    auto& t=p.recv[i];MPI_Irecv(g_biotMassCommDevice?t.device:t.host.data(),
        t.map->nodecount*width*fields*sizeof(solver_float),MPI_BYTE,t.map->procid,tag,comm_solver,&p.rr[i]);
  }
  profile.mark(5);
  if(aggregate)psolve_biot_mass_pack(stream,g_biotMass.nodes,scalar,field,fields,
      p.gather_count,p.gather,values);
  for(auto& t:p.send) {
    const int count=t.map->nodecount,grid=(count+255)/256;
    for(int f=0;!aggregate && f<fields;++f) {
      auto* in=values+(scalar?6:3*(field+f))*(size_t)g_biotMass.nodes;
      auto* out=t.device+(size_t)f*width*count;
      if(scalar)psolve_launch_scalar_boundary_gather(grid,256,0,stream,count,t.map->scalar_mapping_device,in,out);
      else psolve_launch_vector_boundary_gather(grid,256,0,stream,count,t.map->scalar_mapping_device,(fvector_t*)in,(fvector_t*)out);
    }
    if(!g_biotMassCommDevice)solver_pore_pressure_device_check(cudaMemcpyAsync(t.host.data(),t.device,
        (size_t)count*width*fields*sizeof(solver_float),cudaMemcpyDeviceToHost,stream),"Biot host packet download");
  }
  if(!p.send.empty())solver_pore_pressure_device_check(cudaStreamSynchronize(stream),"Biot mass packet ready");
  profile.mark(6);
  for(size_t i=0;i<p.send.size();++i) {
    auto& t=p.send[i];MPI_Isend(g_biotMassCommDevice?t.device:t.host.data(),
        t.map->nodecount*width*fields*sizeof(solver_float),MPI_BYTE,t.map->procid,tag,comm_solver,&p.sr[i]);
  }
  profile.mark(7);
  if(!p.recv.empty())MPI_Waitall(p.recv.size(),p.rr.data(),p.rs.data());
  profile.mark(8);
  for(auto& t:p.recv) {
    const int count=t.map->nodecount,grid=(count+255)/256;
    if(!g_biotMassCommDevice)solver_pore_pressure_device_check(cudaMemcpyAsync(t.device,t.host.data(),
        (size_t)count*width*fields*sizeof(solver_float),cudaMemcpyHostToDevice,stream),"Biot host packet upload");
    for(int f=0;!aggregate && f<fields;++f) {
      auto* out=values+(scalar?6:3*(field+f))*(size_t)g_biotMass.nodes;
      auto* in=t.device+(size_t)f*width*count;
      if(scalar) {
        if(phase<2)psolve_launch_scalar_boundary_accumulate(grid,256,0,stream,count,t.map->scalar_mapping_device,in,out);
        else if(t.map->mapping_has_duplicates)psolve_launch_scalar_boundary_scatter_serial(1,1,0,stream,count,t.map->scalar_mapping_device,in,out);
        else psolve_launch_scalar_boundary_scatter(grid,256,0,stream,count,t.map->scalar_mapping_device,in,out);
      } else {
        if(phase<2)psolve_launch_vector_boundary_accumulate(grid,256,0,stream,count,t.map->scalar_mapping_device,(fvector_t*)in,(fvector_t*)out);
        else if(t.map->mapping_has_duplicates)psolve_launch_vector_boundary_scatter_serial(1,1,0,stream,count,t.map->scalar_mapping_device,(fvector_t*)in,(fvector_t*)out);
        else psolve_launch_vector_boundary_scatter(grid,256,0,stream,count,t.map->scalar_mapping_device,(fvector_t*)in,(fvector_t*)out);
      }
    }
  }
  if(aggregate)psolve_biot_mass_unpack(stream,g_biotMass.nodes,scalar,field,fields,
      phase<2,p.scatter_nodes,p.offsets,p.scatter,values);
  solver_pore_pressure_device_check(cudaGetLastError(),"Biot mass packet kernels");
  if(!p.recv.empty()) {
    if(g_gpuMpiOverlap)p.fence.record(stream);
    else solver_pore_pressure_device_check(cudaStreamSynchronize(stream),"Biot mass unpack wait");
  }
  profile.mark(9);
  if(!p.send.empty())MPI_Waitall(p.send.size(),p.sr.data(),p.ss.data());
  profile.mark(10);
}

static void biot_mass_comm_adjust(mysolver_t* s,solver_float* v,int scalar,int how) {
  if(!g_biotMassCompact) {
    if(scalar)solver_scalar_device_adjust(s,v,how,ACOUSTIC_CONSTRAINT_PORE);
    else solver_vector_device_adjust(s,(fvector_t*)v,how);
    return;
  }
  const auto& m=g_biotMassConstraints[scalar];
  const int count=how==DISTRIBUTION?m.run_parent_base3.size():m.node_base3.size();
  if(!count)return;const int grid=(count+255)/256;const auto stream=s->streams[CUDA_STREAM_MAIN];
  if(scalar) {
    if(how==DISTRIBUTION)psolve_launch_scalar_adjust_distribution(grid,256,0,stream,count,m.d_runs,m.d_run_parent,m.d_edges,m.d_node,m.d_weights,v);
    else psolve_launch_scalar_adjust_assignment(grid,256,0,stream,count,m.d_order,m.d_offsets,m.d_parent,m.d_node,m.d_weights,v);
  } else {
    if(how==DISTRIBUTION)psolve_launch_vector_adjust_distribution(grid,256,0,stream,count,m.d_runs,m.d_run_parent,m.d_edges,m.d_node,m.d_weights,(fvector_t*)v);
    else psolve_launch_vector_adjust_assignment(grid,256,0,stream,count,m.d_order,m.d_offsets,m.d_parent,m.d_node,m.d_weights,(fvector_t*)v);
  }
}
static void biot_mass_comm_exchange(mysolver_t* s,solver_float* v,int scalar,bool rhs) {
  const int first=rhs?0:2,how=rhs?DISTRIBUTION:ASSIGNMENT;
  const int fields=scalar?1:(g_biotMassHaloBatch?2:1);
  for(int f=0;f<(scalar?1:2);f+=fields) {
    biot_mass_comm_phase(s,v,scalar,first,f,fields);
    for(int j=0;j<fields;++j)biot_mass_comm_adjust(s,v+(scalar?6:3*(f+j))*(size_t)g_biotMass.nodes,scalar,how);
    biot_mass_comm_phase(s,v,scalar,first+1,f,fields);
  }
}
static void biot_mass_comm_exchange_pair(mysolver_t* s,solver_float* v,bool rhs) {
  const int first=rhs?0:2,how=rhs?DISTRIBUTION:ASSIGNMENT;
  biot_mass_comm_phase(s,v,2,first,0,1);
  biot_mass_comm_adjust(s,v,0,how);
  biot_mass_comm_adjust(s,v+3*(size_t)g_biotMass.nodes,0,how);
  biot_mass_comm_adjust(s,v+6*(size_t)g_biotMass.nodes,1,how);
  biot_mass_comm_phase(s,v,2,first+1,0,1);
}
static void biot_mass_comm_finish() {
  if(!g_biotMassCommReady)return;
  for(auto& p:g_biotMassPlans) {
    p.fence.destroy();for(auto* list:{&p.send,&p.recv})for(auto& v:*list)cudaFree(v.device);
    cudaFree(p.gather);cudaFree(p.scatter);cudaFree(p.offsets);
    p={};
  }
  for(auto*& s:g_biotMassSchedules){schedule_delete(s);s=nullptr;}
  for(auto& m:g_biotMassConstraints) {
    for(auto* p:{m.d_node,m.d_offsets,m.d_parent,m.d_order,m.d_edges,m.d_runs,m.d_run_parent})cudaFree(p);
    cudaFree(m.d_weights);m={};
  }
  g_biotMassCommReady=false;
}
#endif
