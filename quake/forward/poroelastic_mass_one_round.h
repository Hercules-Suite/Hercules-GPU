#ifndef HERCULES_POROELASTIC_MASS_ONE_ROUND_H
#define HERCULES_POROELASTIC_MASS_ONE_ROUND_H
#include <unordered_map>

/* Experimental owner-complete Cartesian rows. Integer Q1 weights are summed
 * once at setup, including owners with no local porous cells. Runtime exchanges
 * only owner input values; no partial output sums remain. Hanging/field-specific
 * geometric constraints and compact vectors retain the general two-exchange
 * operator. No physical/interface constraint is approximated by this path. */
struct biot_mass_owner_row_t {
  int64_t gid,neighbor[27];
  int32_t owner[27];
  unsigned char weight[27];
};
struct biot_mass_input_peer_t {
  int rank,count,offset;
  int32_t* map=nullptr;
  solver_float* device=nullptr;
  std::vector<solver_float> host;
};
static struct {
  std::vector<biot_mass_input_peer_t> send,recv;
  std::vector<MPI_Request> sr,rr;
  cudaEvent_t used=nullptr;
  bool pending=false;
  bool device=false;
} g_biotMassInput;

/* Setup-only all-to-all. Byte counts and displacements are checked before MPI;
 * the data are node IDs/integer topology, never time-dependent field values. */
template<class T> static std::vector<T> biot_mass_setup_exchange(
    const std::vector<std::vector<T>>& lists,std::vector<int>& received) {
  const int ranks=Global.theGroupSize;
  std::vector<int> send(ranks),so(ranks),ro(ranks);received.resize(ranks);
  size_t ns=0,nr=0;
  for(int p=0;p<ranks;++p) {
    const size_t bytes=lists[p].size()*sizeof(T);
    if(bytes>INT_MAX || ns+bytes>INT_MAX)
      solver_abort(__FUNCTION_NAME,NULL,"Biot owner-row setup exceeds MPI byte count\n");
    send[p]=bytes;so[p]=ns;ns+=bytes;
  }
  MPI_Alltoall(send.data(),1,MPI_INT,received.data(),1,MPI_INT,comm_solver);
  for(int p=0;p<ranks;++p) {
    if(received[p]<0 || received[p]%sizeof(T) || nr+(size_t)received[p]>INT_MAX)
      solver_abort(__FUNCTION_NAME,NULL,"Invalid Biot owner-row setup receive count\n");
    ro[p]=nr;nr+=received[p];
  }
  std::vector<T> packed(ns/sizeof(T)),result(nr/sizeof(T));
  for(int p=0;p<ranks;++p)std::copy(lists[p].begin(),lists[p].end(),packed.begin()+so[p]/sizeof(T));
  MPI_Alltoallv(packed.data(),send.data(),so.data(),MPI_BYTE,
      result.data(),received.data(),ro.data(),MPI_BYTE,comm_solver);
  for(auto& n:received)n/=sizeof(T);
  return result;
}

static void biot_mass_one_round_setup(mesh_t* mesh,std::vector<int32_t>& nodes,
    std::vector<int32_t>& neighbor,std::vector<unsigned char>& weight) {
  int unsupported=g_biotMassCompact;
  for(int n=0;n<mesh->nharbored;++n)unsupported|=!mesh->nodeTable[n].isanchored;
  for(int f:{ACOUSTIC_CONSTRAINT_SOLID,ACOUSTIC_CONSTRAINT_PORE})
    unsupported|=!g_acousticConstraintMaps[f].node_base3.empty();
  int any;MPI_Allreduce(&unsupported,&any,1,MPI_INT,MPI_MAX,comm_solver);
  if(any)return;
  const int rank=Global.myID,ranks=Global.theGroupSize,N=mesh->nharbored;
  int64_t first=INT64_MAX,last=-1;size_t owned_count=0;
  for(int n=0;n<N;++n)if(mesh->nodeTable[n].ismine) {
    const int64_t id=mesh->nodeTable[n].gnid;
    if(id<0)solver_abort(__FUNCTION_NAME,NULL,"Negative owned Biot node ID\n");
    first=std::min(first,id);last=std::max(last,id);++owned_count;
  }
  const bool dense=owned_count && (uint64_t)(last-first)<2*owned_count && last-first<INT_MAX;
  std::vector<int32_t> lookup(dense?(size_t)(last-first+1):0,-1);
  std::unordered_map<int64_t,int32_t> sparse;
  if(!dense)sparse.reserve(owned_count);
  for(int n=0;n<N;++n)if(mesh->nodeTable[n].ismine) {
    if(dense)lookup[mesh->nodeTable[n].gnid-first]=n;
    else sparse.emplace(mesh->nodeTable[n].gnid,n);
  }
  const auto local_node=[&](int64_t gid) {
    int n=-1;
    if(dense && gid>=first && gid<=last)n=lookup[gid-first];
    if(!dense){auto i=sparse.find(gid);if(i!=sparse.end())n=i->second;}
    if(n<0)solver_abort(__FUNCTION_NAME,NULL,"Unknown Biot owner-row node %lld\n",(long long)gid);
    return n;
  };
  std::vector<std::vector<biot_mass_owner_row_t>> outgoing(ranks);
  std::vector<biot_mass_owner_row_t> rows;
  std::vector<int32_t> row_of(N,-1),owned_nodes;
  for(size_t i=0;i<nodes.size();++i) {
    const int n=nodes[i];const auto& node=mesh->nodeTable[n];
    const int owner=node.ismine?rank:node.proc.ownerid;
    if(owner<0 || owner>=ranks)solver_abort(__FUNCTION_NAME,NULL,"Invalid Biot row owner\n");
    biot_mass_owner_row_t row={};row.gid=node.gnid;
    for(int j=0;j<27;++j)if((row.weight[j]=weight[(size_t)j*nodes.size()+i])) {
      const auto& other=mesh->nodeTable[neighbor[(size_t)j*nodes.size()+i]];
      row.neighbor[j]=other.gnid;row.owner[j]=other.ismine?rank:other.proc.ownerid;
    }
    if(owner==rank) {row_of[n]=rows.size();owned_nodes.push_back(n);rows.push_back(row);}
    else outgoing[owner].push_back(row);
  }
  std::vector<int> counts;
  auto incoming=biot_mass_setup_exchange(outgoing,counts);
  for(const auto& row:incoming) {
    const int n=local_node(row.gid);
    if(row_of[n]<0) {row_of[n]=rows.size();owned_nodes.push_back(n);rows.push_back({});rows.back().gid=row.gid;}
    auto& target=rows[row_of[n]];
    for(int j=0;j<27;++j)if(row.weight[j]) {
      if((target.weight[j] && (target.neighbor[j]!=row.neighbor[j] || target.owner[j]!=row.owner[j])) ||
          target.weight[j]+row.weight[j]>64)
        solver_abort(__FUNCTION_NAME,NULL,"Inconsistent assembled Biot owner row\n");
      target.neighbor[j]=row.neighbor[j];target.owner[j]=row.owner[j];target.weight[j]+=row.weight[j];
    }
  }
  outgoing.clear();incoming.clear();
  std::map<std::pair<int,int64_t>,int> remote;
  for(const auto& row:rows)for(int j=0;j<27;++j)if(row.weight[j] && row.owner[j]!=rank) {
    if(row.owner[j]<0 || row.owner[j]>=ranks)solver_abort(__FUNCTION_NAME,NULL,"Invalid Biot column owner\n");
    remote.emplace(std::make_pair(row.owner[j],row.neighbor[j]),0);
  }
  if(remote.size()>INT_MAX/7 || rows.size()>INT_MAX/27)
    solver_abort(__FUNCTION_NAME,NULL,"Biot owner-row mapping exceeds integer limits\n");
  std::vector<std::vector<int64_t>> requests(ranks);
  int slot=0;for(auto& entry:remote) {entry.second=slot++;requests[entry.first.first].push_back(entry.first.second);}
  neighbor.assign(27*rows.size(),-1);weight.assign(27*rows.size(),0);
  for(size_t i=0;i<rows.size();++i)for(int j=0;j<27;++j)if(rows[i].weight[j]) {
    const auto& row=rows[i];const size_t k=(size_t)j*rows.size()+i;
    neighbor[k]=row.owner[j]==rank?local_node(row.neighbor[j]):
        -1-remote.at(std::make_pair(row.owner[j],row.neighbor[j]));
    weight[k]=row.weight[j];
  }
  nodes.swap(owned_nodes);
  auto wanted=biot_mass_setup_exchange(requests,counts);
  auto& plan=g_biotMassInput;plan.device=gpu_device_sync_enabled();solver_pore_pressure_device_check(cudaEventCreateWithFlags(&plan.used,cudaEventDisableTiming),"Biot input reuse event");
  size_t at=0;int offset=0;
  for(int p=0;p<ranks;++p) {
    if(counts[p]) {
      if(counts[p]>INT_MAX/(7*sizeof(solver_float)))solver_abort(__FUNCTION_NAME,NULL,"Biot input packet too large\n");
      plan.send.push_back({p,counts[p],0});auto& peer=plan.send.back();
      std::vector<int32_t> map(counts[p]);
      for(int i=0;i<counts[p];++i)map[i]=local_node(wanted[at+i]);
      acoustic_constraint_upload(&peer.map,map);
      solver_pore_pressure_device_check(cudaMalloc((void**)&peer.device,7*(size_t)peer.count*sizeof(solver_float)),"Biot input packet");
      if(!plan.device)peer.host.resize(7*(size_t)peer.count);
    }
    at+=counts[p];
    if(!requests[p].empty()) {
      if(requests[p].size()>INT_MAX/(7*sizeof(solver_float)))solver_abort(__FUNCTION_NAME,NULL,"Biot input receive too large\n");
      plan.recv.push_back({p,(int)requests[p].size(),offset});auto& peer=plan.recv.back();
      if(!plan.device)peer.host.resize(7*(size_t)peer.count);
      offset+=peer.count;
    }
  }
  auto& m=g_biotMass;
  solver_pore_pressure_device_check(cudaMalloc((void**)&m.remote_values,
      7*std::max((size_t)1,remote.size())*sizeof(solver_float)),"Biot remote row inputs");
  plan.sr.resize(plan.send.size());plan.rr.resize(plan.recv.size());
  m.owner_rows=1;g_biotMassOneRound=true;
  unsigned long long bytes=7*std::max((size_t)1,remote.size())*sizeof(solver_float),maximum=0;
  for(const auto& peer:plan.send)bytes+=7*(size_t)peer.count*sizeof(solver_float);
  MPI_Allreduce(&bytes,&maximum,1,MPI_UNSIGNED_LONG_LONG,MPI_MAX,comm_solver);
  if(rank==0)monitor_print("Biot one-round mass: experimental owner-complete Q1 rows; one input halo, no output assembly; %s MPI.\n",
      plan.device?"CUDA-aware":"host-staged");
  if(rank==0)monitor_print("Biot owner-row input packets: maximum %.3f GiB/rank on GPU; setup topology workspace is additional.\n",
      (double)maximum/(1024.*1024.*1024.));
}

static void biot_mass_one_round_exchange(mysolver_t* s,solver_float* in,int scalar) {
  auto& p=g_biotMassInput;auto& m=g_biotMass;
  const auto stream=s->streams[CUDA_STREAM_MAIN];const int width=scalar==1?1:(scalar==2?7:6);
  // The previous row kernel must release receive storage before MPI writes it.
  if(p.pending)solver_pore_pressure_device_check(cudaEventSynchronize(p.used),"Biot input reuse wait");
  for(size_t i=0;i<p.recv.size();++i) {
    auto& peer=p.recv[i];
    MPI_Irecv(p.device?m.remote_values+width*(size_t)peer.offset:peer.host.data(),
        width*peer.count*sizeof(solver_float),MPI_BYTE,peer.rank,184,comm_solver,&p.rr[i]);
  }
  for(auto& peer:p.send) {
    psolve_biot_mass_input_pack(stream,m.nodes,scalar,peer.count,peer.map,in,peer.device);
    if(!p.device)solver_pore_pressure_device_check(cudaMemcpyAsync(peer.host.data(),peer.device,
        width*(size_t)peer.count*sizeof(solver_float),cudaMemcpyDeviceToHost,stream),"Biot input host packet");
  }
  if(!p.send.empty())solver_pore_pressure_device_check(cudaStreamSynchronize(stream),"Biot input packet ready");
  for(size_t i=0;i<p.send.size();++i) {
    auto& peer=p.send[i];MPI_Isend(p.device?peer.device:peer.host.data(),
        width*peer.count*sizeof(solver_float),MPI_BYTE,peer.rank,184,comm_solver,&p.sr[i]);
  }
  if(!p.recv.empty())MPI_Waitall(p.recv.size(),p.rr.data(),MPI_STATUSES_IGNORE);
  if(!p.device)for(auto& peer:p.recv)solver_pore_pressure_device_check(cudaMemcpyAsync(
      m.remote_values+width*(size_t)peer.offset,peer.host.data(),width*(size_t)peer.count*sizeof(solver_float),
      cudaMemcpyHostToDevice,stream),"Biot input receive upload");
  if(!p.send.empty())MPI_Waitall(p.send.size(),p.sr.data(),MPI_STATUSES_IGNORE);
}
static void biot_mass_one_round_used(cudaStream_t stream) {
  solver_pore_pressure_device_check(cudaEventRecord(g_biotMassInput.used,stream),"Biot input reuse record");
  g_biotMassInput.pending=true;
}
static void biot_mass_one_round_finish() {
  if(!g_biotMassOneRound)return;
  auto& p=g_biotMassInput;
  if(p.pending)solver_pore_pressure_device_check(cudaEventSynchronize(p.used),"Biot input final wait");
  cudaEventDestroy(p.used);
  for(auto& peer:p.send){cudaFree(peer.map);cudaFree(peer.device);}
  cudaFree(g_biotMass.remote_values);g_biotMass.remote_values=nullptr;
  p={};g_biotMassOneRound=false;
}
#endif
