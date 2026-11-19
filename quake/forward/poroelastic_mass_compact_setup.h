#ifndef HERCULES_POROELASTIC_MASS_COMPACT_SETUP_H
#define HERCULES_POROELASTIC_MASS_COMPACT_SETUP_H

static bool biot_mass_option(const char* name,int default_value=1) {
  const char* value=getenv(name);
  if (value && strcmp(value,"0") && strcmp(value,"1"))
    solver_abort(__FUNCTION_NAME,NULL,"%s must be 0 or 1\n",name);
  int local=value?!strcmp(value,"1"):default_value,lo,hi;
  MPI_Allreduce(&local,&lo,1,MPI_INT,MPI_MIN,comm_solver);
  MPI_Allreduce(&local,&hi,1,MPI_INT,MPI_MAX,comm_solver);
  if (lo!=hi) solver_abort(__FUNCTION_NAME,NULL,"%s differs between ranks\n",name);
  return local;
}

/* Keep the union of mechanical and pressure dependencies. Only the directed
 * child->parent closure is needed: an unrelated child of an active parent has
 * no consistent-cell contribution. MPI OR/sharing retains remote owners even
 * on ranks with no physical porous cell. Never classify by local cell count. */
static void biot_mass_compact_setup(mesh_t* mesh,const std::vector<int32_t>& elements,
    const std::vector<fvector_t>& removed,const std::vector<solver_float>& capacity,
    const std::vector<unsigned char>& owned) {
  auto& m=g_biotMass;
  if (!g_biotMassCompact) return;
  std::vector<uint8_t> active(m.full_nodes,0);
  for (int pi:elements) for (int j=0;j<8;++j)
    active[mesh->elemTable[g_poroelasticElemMapper[pi]].lnid[j]]=1;
  compute_adjust_cache_build();
  const int fields[2]={ACOUSTIC_CONSTRAINT_SOLID,ACOUSTIC_CONSTRAINT_PORE};
  int changed=1;
  while (changed) {
    changed=0;
    for (int field:fields) {
      const auto& map=g_acousticConstraintMaps[field];
      const auto mark=[&](int n,int parent) {
        if (active[n] && !active[parent]) {active[parent]=1;changed=1;}
      };
      if (map.enabled) {
        for (size_t d=0;d<map.node_base3.size();++d)
          for (int j=map.offsets[d];j<map.offsets[d+1];++j)
            mark(map.node_base3[d]/3,map.parent_base3[j]/3);
      } else for (int d=0;d<g_adjustDnCount;++d)
        for (int j=g_adjustParentOffsets[d];j<g_adjustParentOffsets[d+1];++j)
          mark(g_adjustDanglingLnid[d],g_adjustParentIds[j]);
    }
    auto* s=Global.mySolver;
    changed|=scalar_mpi_share_mask(s->dn_sched,CONTRIBUTION,active);
    changed|=scalar_mpi_share_mask(s->an_sched,CONTRIBUTION,active);
    changed|=scalar_mpi_share_mask(s->an_sched,SHARING,active);
    changed|=scalar_mpi_share_mask(s->dn_sched,SHARING,active);
    int global;MPI_Allreduce(&changed,&global,1,MPI_INT,MPI_MAX,comm_solver);changed=global;
  }
  std::vector<int32_t> nodes,diagonal;
  g_biotMassIndex.assign(m.full_nodes,-1);
  for (int n=0;n<m.full_nodes;++n) {
    if (active[n]) {g_biotMassIndex[n]=nodes.size();nodes.push_back(n);}
    else {
      if (((owned[n]&1) && removed[n].f[0]!=0) || ((owned[n]&2) && capacity[n]!=0))
        solver_abort(__FUNCTION_NAME,NULL,"Biot compact closure lost an assembled consistent DOF\n");
      diagonal.push_back(n);
    }
  }
  m.nodes=nodes.size();m.diagonal_nodes=diagonal.size();
  acoustic_constraint_upload(&m.node,nodes);
  acoustic_constraint_upload(&m.index,g_biotMassIndex);
  acoustic_constraint_upload(&m.diagonal,diagonal);
  for (int f=0;f<2;++f) {
    auto& out=g_biotMassConstraints[f];const auto& in=g_acousticConstraintMaps[fields[f]];
    out.enabled=true;out.offsets.push_back(0);
    const int count=in.enabled?in.node_base3.size():g_adjustDnCount;
    for (int d=0;d<count;++d) {
      const int n=in.enabled?in.node_base3[d]/3:g_adjustDanglingLnid[d];
      if (!active[n]) continue;
      out.order.push_back(out.node_base3.size());out.node_base3.push_back(3*g_biotMassIndex[n]);
      out.weights.push_back(in.enabled?in.weights[d]:g_adjustInvDeps[d]);
      const int lo=in.enabled?in.offsets[d]:g_adjustParentOffsets[d];
      const int hi=in.enabled?in.offsets[d+1]:g_adjustParentOffsets[d+1];
      for (int j=lo;j<hi;++j) {
        const int parent=in.enabled?in.parent_base3[j]/3:g_adjustParentBase3[j]/3;
        if (g_biotMassIndex[parent]<0) solver_abort(__FUNCTION_NAME,NULL,"Missing compact constraint parent\n");
        out.parent_base3.push_back(3*g_biotMassIndex[parent]);
      }
      out.offsets.push_back(out.parent_base3.size());
    }
    acoustic_constraint_finish(out);
    acoustic_constraint_upload(&out.d_node,out.node_base3);
    acoustic_constraint_upload(&out.d_offsets,out.offsets);
    acoustic_constraint_upload(&out.d_parent,out.parent_base3);
    acoustic_constraint_upload(&out.d_order,out.order);
    acoustic_constraint_upload(&out.d_weights,out.weights);
    acoustic_constraint_upload(&out.d_edges,out.edges);
    acoustic_constraint_upload(&out.d_runs,out.runs);
    acoustic_constraint_upload(&out.d_run_parent,out.run_parent_base3);
  }
  long long counts[2]={m.nodes,m.full_nodes},sum[2],lo[2],hi[2];
  MPI_Reduce(counts,sum,2,MPI_LONG_LONG,MPI_SUM,0,comm_solver);
  MPI_Reduce(counts,lo,2,MPI_LONG_LONG,MPI_MIN,0,comm_solver);
  MPI_Reduce(counts,hi,2,MPI_LONG_LONG,MPI_MAX,0,comm_solver);
  if (Global.myID==0) monitor_print("Biot compact vectors: %lld / %lld harbored node copies retained; "
      "rank retained min/max %lld/%lld, full min/max %lld/%lld. "
      "Includes both field constraints and MPI owners; isolated nodal blocks solved directly.\n",
      sum[0],sum[1],lo[0],hi[0],lo[1],hi[1]);
}
#endif
