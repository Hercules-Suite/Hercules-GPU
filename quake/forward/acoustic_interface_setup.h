#ifndef HERCULES_ACOUSTIC_INTERFACE_SETUP_H
#define HERCULES_ACOUSTIC_INTERFACE_SETUP_H

/* Private psolve.cc implementation: include after material/neighbor discovery.
 * The original geometric constraints remain the default. Both alternatives
 * are opt-in experiments for full Hex8 acoustic/linear-solid meshes.
 */
struct acoustic_constraint_map_t {
  std::vector<int32_t> node_base3, offsets, parent_base3, order;
  std::vector<int32_t> edges, runs, run_parent_base3;
  std::vector<solver_float> weights;
  int32_t *d_node = NULL, *d_offsets = NULL, *d_parent = NULL, *d_order = NULL;
  int32_t *d_edges = NULL, *d_runs = NULL, *d_run_parent = NULL;
  solver_float* d_weights = NULL;
};
static acoustic_constraint_map_t g_acousticConstraintMaps[2];
static bool g_acousticConstraintsReady = false;

static void acoustic_constraint_append(acoustic_constraint_map_t& map,
                                        const dnode_t& node) {
  const int32_t id = (int32_t)map.node_base3.size();
  map.node_base3.push_back(3 * node.ldnid);
  map.weights.push_back((solver_float)1.0 / node.deps);
  map.order.push_back(id);
  if (map.offsets.empty()) map.offsets.push_back(0);
  for (int32link_t* parent = node.lanid; parent; parent = parent->next)
    map.parent_base3.push_back(3 * parent->id);
  if (map.parent_base3.size() > (size_t)INT_MAX)
    solver_abort(__FUNCTION_NAME, NULL, "Field constraint map exceeds INT_MAX\n");
  map.offsets.push_back((int32_t)map.parent_base3.size());
}

static void acoustic_constraint_finish(acoustic_constraint_map_t& map) {
  // One distribution thread per parent: deterministic summation, no atomics.
  std::vector<std::pair<int32_t, int32_t>> edges;
  for (int32_t i = 0; i < (int32_t)map.node_base3.size(); ++i)
    for (int32_t j = map.offsets[i]; j < map.offsets[i + 1]; ++j)
      edges.push_back({map.parent_base3[j], i});
  std::sort(edges.begin(), edges.end());
  for (const auto& edge : edges) {
    if (map.run_parent_base3.empty() || map.run_parent_base3.back() != edge.first) {
      map.run_parent_base3.push_back(edge.first);
      map.runs.push_back((int32_t)map.edges.size());
    }
    map.edges.push_back(edge.second);
  }
  map.runs.push_back((int32_t)map.edges.size());
}

// Probe immediately inside each of the eight octants incident on a node.
// Remote aggregate octants are resolved to actual leaves by their MPI owner.
static int acoustic_constraint_probe(const node_t& node, int octant,
                                      tick_t low[3], tick_t high[3], int* field) {
  const double tick = Global.myMesh->ticksize;
  vector3D_t p = {{(node.x + ((octant & 1) ? 0.25 : -0.25)) * tick,
                   (node.y + ((octant & 2) ? 0.25 : -0.25)) * tick,
                   (node.z + ((octant & 4) ? 0.25 : -0.25)) * tick}};
  if (p.x[0] < 0 || p.x[0] >= Param.theDomainX ||
      p.x[1] < 0 || p.x[1] >= Param.theDomainY ||
      p.x[2] < 0 || p.x[2] >= Param.theDomainZ) return 0;
  int local = 0;
  int32_t element = -1;
  if (!mesh_find_octant_at_domain_point(&p, low, high, NULL, &local, &element)) {
    if (!g_acousticCollectNeighbors)
      solver_abort(__FUNCTION_NAME, NULL, "Unresolved incident field element\n");
    return 0;
  }
  double vp, vs, rho;
  if (local) {
    if (element < 0)
      solver_abort(__FUNCTION_NAME, NULL, "Missing local field element\n");
    const edata_t* data = (edata_t*)Global.myMesh->elemTable[element].data;
    vp = data->Vp; vs = data->Vs; rho = data->rho;
  } else {
    const acoustic_neighbor_t* neighbor = mesh_cached_acoustic_neighbor_at(&p);
    if (!neighbor || neighbor->geid < 0)
      solver_abort(__FUNCTION_NAME, NULL, "Missing remote field element\n");
    vp = neighbor->vp; vs = neighbor->vs; rho = neighbor->rho;
  }
  if (material_values_are_acoustic(vp, vs, rho)) *field = ACOUSTIC_CONSTRAINT_SCALAR;
  else if (material_values_are_elastic(vp, vs, rho)) *field = ACOUSTIC_CONSTRAINT_SOLID;
  else solver_abort(__FUNCTION_NAME, NULL, "Unsupported material in field constraint\n");
  return 1;
}

template <class T>
static void acoustic_constraint_upload(T** device, const std::vector<T>& host) {
  if (host.empty()) return;
  cudaError_t error = cudaMalloc((void**)device, host.size() * sizeof(T));
  if (error == cudaSuccess)
    error = cudaMemcpy(*device, host.data(), host.size() * sizeof(T), cudaMemcpyHostToDevice);
  if (error != cudaSuccess)
    solver_abort(__FUNCTION_NAME, NULL, "Field constraint upload: %s\n", cudaGetErrorString(error));
}

static void acoustic_constraints_build() {
  if (Param.theAcousticInterfaceMethod != ACOUSTIC_INTERFACE_FIELDS) return;
  mesh_t* mesh = Global.myMesh;
  if (mesh->nharbored > INT_MAX / 3)
    solver_abort(__FUNCTION_NAME, NULL, "Field constraint node index overflow\n");
  for (int pass = 0; ; ++pass) {
    if (pass > PIXELLEVEL + 1)
      solver_abort(__FUNCTION_NAME, NULL, "Field neighbor discovery did not converge\n");
    g_acousticCollectNeighbors = 1;
    for (int32_t i = 0; i < mesh->ldnnum; ++i)
      for (int o = 0; o < 8; ++o) {
        tick_t low[3], high[3]; int field;
        acoustic_constraint_probe(mesh->nodeTable[mesh->dnodeTable[i].ldnid],
                                  o, low, high, &field);
      }
    g_acousticCollectNeighbors = 0;
    int pending = !g_acousticNeighborRequests.empty(), any = 0;
    MPI_Allreduce(&pending, &any, 1, MPI_INT, MPI_MAX, comm_solver);
    if (!any) break;
    mesh_exchange_acoustic_neighbors();
  }
  int64_t counts[6] = {0}, total[6]; // keep/release/absent for each field
  for (int32_t i = 0; i < mesh->ldnnum; ++i) {
    const dnode_t& dangling = mesh->dnodeTable[i];
    const node_t& node = mesh->nodeTable[dangling.ldnid];
    bool active[2] = {false, false}, keep[2] = {false, false};
    std::map<point_t, double, acoustic_point_less> given;
    if (dangling.deps == 0)
      solver_abort(__FUNCTION_NAME, NULL, "Field constraint has no parents\n");
    for (int32link_t* p = dangling.lanid; p; p = p->next) {
      const node_t& parent = mesh->nodeTable[p->id];
      if (!parent.isanchored)
        solver_abort(__FUNCTION_NAME, NULL, "Field constraint requires anchored parents\n");
      given[{parent.x, parent.y, parent.z}] += 1.0 / dangling.deps;
    }
    for (int o = 0; o < 8; ++o) {
      tick_t low[3], high[3]; int field;
      if (!acoustic_constraint_probe(node, o, low, high, &field)) continue;
      active[field] = true;
      const tick_t pos[3] = {node.x, node.y, node.z};
      bool vertex = true;
      for (int j = 0; j < 3; ++j)
        if (pos[j] != low[j] && pos[j] != high[j]) vertex = false;
      if (vertex) continue;
      keep[field] = true;
      // Retain a constraint only when it is the exact same-field Hex8 trace.
      // Mixed-corner meshes can require different stencils; do not guess them.
      std::map<point_t, double, acoustic_point_less> expected;
      for (int c = 0; c < 8; ++c) {
        double weight = 1.0; tick_t v[3];
        for (int j = 0; j < 3; ++j) {
          const bool upper = c & (1 << j);
          v[j] = upper ? high[j] : low[j];
          const double f = double(pos[j] - low[j]) / double(high[j] - low[j]);
          weight *= upper ? f : 1.0 - f;
        }
        if (weight > 0.0) expected[{v[0], v[1], v[2]}] += weight;
      }
      bool equal = given.size() == expected.size();
      for (const auto& p : expected) {
        const auto found = given.find(p.first);
        if (found == given.end() || fabs(found->second - p.second) > 1.e-12) equal = false;
      }
      if (!equal)
        solver_abort(__FUNCTION_NAME, NULL,
                     "Unsupported fields-constraint stencil at (%g,%g,%g) m, field %d. "
                     "Use original-constraint or remesh.\n",
                     node.x * mesh->ticksize, node.y * mesh->ticksize,
                     node.z * mesh->ticksize, field);
    }
    for (int field = 0; field < 2; ++field) {
      ++counts[3 * field + (keep[field] ? 0 : (active[field] ? 1 : 2))];
      if (keep[field]) acoustic_constraint_append(g_acousticConstraintMaps[field], dangling);
    }
  }
  for (auto& map : g_acousticConstraintMaps) {
    acoustic_constraint_finish(map);
    acoustic_constraint_upload(&map.d_node, map.node_base3);
    acoustic_constraint_upload(&map.d_offsets, map.offsets);
    acoustic_constraint_upload(&map.d_parent, map.parent_base3);
    acoustic_constraint_upload(&map.d_order, map.order);
    acoustic_constraint_upload(&map.d_edges, map.edges);
    acoustic_constraint_upload(&map.d_runs, map.runs);
    acoustic_constraint_upload(&map.d_run_parent, map.run_parent_base3);
    acoustic_constraint_upload(&map.d_weights, map.weights);
  }
  MPI_Allreduce(counts, total, 6, MPI_INT64_T, MPI_SUM, comm_solver);
  if (Global.myID == 0)
    monitor_print("Field constraints (keep/release/absent): acoustic %lld/%lld/%lld; "
                  "solid %lld/%lld/%lld. Same-field Hex8 traces checked.\n",
                  (long long)total[0], (long long)total[1], (long long)total[2],
                  (long long)total[3], (long long)total[4], (long long)total[5]);
  g_acousticConstraintsReady = true;
}

static void acoustic_constraints_free() {
  for (auto& map : g_acousticConstraintMaps) {
    cudaFree(map.d_node); cudaFree(map.d_offsets); cudaFree(map.d_parent);
    cudaFree(map.d_order); cudaFree(map.d_edges); cudaFree(map.d_runs);
    cudaFree(map.d_run_parent); cudaFree(map.d_weights);
    map = acoustic_constraint_map_t();
  }
  g_acousticConstraintsReady = false;
}

static int acoustic_constraints_host_adjust(void* values, int width, int how, int field) {
  if (!g_acousticConstraintsReady || field < 0) return 0;
  const auto& map = g_acousticConstraintMaps[field];
  solver_float* v = (solver_float*)values;
  for (int32_t i = 0; i < (int32_t)map.node_base3.size(); ++i) {
    solver_float* node = v + (map.node_base3[i] / 3) * (size_t)width;
    if (how != DISTRIBUTION) std::fill(node, node + width, (solver_float)0);
    for (int32_t j = map.offsets[i]; j < map.offsets[i + 1]; ++j) {
      solver_float* parent = v + (map.parent_base3[j] / 3) * (size_t)width;
      for (int c = 0; c < width; ++c) {
        if (how == DISTRIBUTION) parent[c] += node[c] * map.weights[i];
        else node[c] += parent[c] * map.weights[i];
      }
    }
  }
  return 1;
}

static int acoustic_constraints_device_adjust(mysolver_t* solver, void* values,
                                               int how, int field) {
  if (!g_acousticConstraintsReady || field < 0) return 0;
  const auto& map = g_acousticConstraintMaps[field];
  const int32_t count = (int32_t)(how == DISTRIBUTION ? map.run_parent_base3.size()
                                                         : map.node_base3.size());
  if (count == 0) return 1;
  const int block = 256, grid = (count + block - 1) / block;
  const cudaStream_t stream = solver->streams[CUDA_STREAM_MAIN];
  if (field == ACOUSTIC_CONSTRAINT_SCALAR) {
    if (how == DISTRIBUTION)
      psolve_launch_scalar_adjust_distribution(grid, block, 0, stream, count,
          map.d_runs, map.d_run_parent, map.d_edges, map.d_node, map.d_weights,
          (solver_float*)values);
    else
      psolve_launch_scalar_adjust_assignment(grid, block, 0, stream, count,
          map.d_order, map.d_offsets, map.d_parent, map.d_node, map.d_weights,
          (solver_float*)values);
  } else {
    if (how == DISTRIBUTION)
      psolve_launch_vector_adjust_distribution(grid, block, 0, stream, count,
          map.d_runs, map.d_run_parent, map.d_edges, map.d_node, map.d_weights,
          (fvector_t*)values);
    else
      psolve_launch_vector_adjust_assignment(grid, block, 0, stream, count,
          map.d_order, map.d_offsets, map.d_parent, map.d_node, map.d_weights,
          (fvector_t*)values);
  }
  const cudaError_t error = cudaGetLastError();
  if (error != cudaSuccess)
    solver_abort(__FUNCTION_NAME, NULL, "Field constraint kernel: %s\n", cudaGetErrorString(error));
  return 1;
}

// Collect coarser interface leaves while the existing interface builder visits
// every candidate face (including faces omitted by the selected force filter).
static void mesh_record_acoustic_mismatch(int32_t eindex, const vector3D_t* probe,
                                           int32_t adjacent_level) {
  const elem_t& element = Global.myMesh->elemTable[eindex];
  if (element.level == adjacent_level) return;
  point_t lower;
  int level;
  if (element.level < adjacent_level) {
    const node_t& node = Global.myMesh->nodeTable[element.lnid[0]];
    lower = {node.x, node.y, node.z}; level = element.level;
  } else {
    tick_t low[3];
    if (!mesh_find_octant_at_domain_point(probe, low, NULL, NULL, NULL, NULL))
      solver_abort(__FUNCTION_NAME, NULL, "Missing local-matching neighbor\n");
    lower = {low[0], low[1], low[2]}; level = adjacent_level;
  }
  g_acousticMatchingRequests[lower] = level;
}

struct acoustic_matching_request_t { point_t lower; int level; };
static std::map<point_t, int, acoustic_point_less> g_acousticMatchingOwned;

static int32_t acoustic_matching_refine(octant_t* leaf, double ticksize,
                                        const void* data) {
  const point_t lower = {leaf->lx, leaf->ly, leaf->lz};
  const auto found = g_acousticMatchingOwned.find(lower);
  // Only split the marked leaf once; a child can have the same lower corner.
  if (found != g_acousticMatchingOwned.end() && leaf->level == found->second)
    return 1;
  // New children of a mixed cell can expose slower material. Preserve the
  // normal wavespeed/interface refinement rules as well as face matching.
  return toexpand(leaf, ticksize, data);
}

static void mesh_exchange_matching_requests() {
  const int ranks = Global.theGroupSize;
  std::vector<int> send_count(ranks, 0), recv_count(ranks, 0);
  std::vector<int> send_offset(ranks, 0), recv_offset(ranks, 0);
  const point_t* intervals = octor_getintervaltable(Global.myOctree);
  if (g_acousticMatchingRequests.size() > (size_t)INT_MAX)
    solver_abort(__FUNCTION_NAME, NULL, "Too many matching requests\n");
  for (const auto& request : g_acousticMatchingRequests)
    ++send_count[mesh_handoff_owner(intervals, ranks, &request.first)];
  MPI_Alltoall(send_count.data(), 1, MPI_INT, recv_count.data(), 1, MPI_INT, comm_solver);
  int64_t nsend = 0, nrecv = 0;
  for (int r = 0; r < ranks; ++r) {
    send_offset[r] = (int)nsend; recv_offset[r] = (int)nrecv;
    nsend += send_count[r]; nrecv += recv_count[r];
    if (nsend > INT_MAX || nrecv > INT_MAX)
      solver_abort(__FUNCTION_NAME, NULL, "Matching MPI count overflow\n");
  }
  std::vector<acoustic_matching_request_t> outgoing(nsend), incoming(nrecv);
  std::vector<int> cursor = send_offset;
  for (const auto& request : g_acousticMatchingRequests) {
    const int owner = mesh_handoff_owner(intervals, ranks, &request.first);
    outgoing[cursor[owner]++] = {request.first, request.second};
  }
  MPI_Datatype type;
  MPI_Type_contiguous(sizeof(acoustic_matching_request_t), MPI_BYTE, &type);
  MPI_Type_commit(&type);
  MPI_Alltoallv(outgoing.data(), send_count.data(), send_offset.data(), type,
                incoming.data(), recv_count.data(), recv_offset.data(), type, comm_solver);
  MPI_Type_free(&type);
  g_acousticMatchingOwned.clear();
  for (const auto& request : incoming) {
    octant_t* leaf = octor_searchoctant(Global.myOctree, request.lower.x,
        request.lower.y, request.lower.z, request.level, EXACT_SEARCH);
    if (!leaf || leaf->where != LOCAL || leaf->type != LEAF ||
        leaf->level >= PIXELLEVEL - 1)
      solver_abort(__FUNCTION_NAME, NULL, "Invalid local-matching refinement request\n");
    g_acousticMatchingOwned[request.lower] = request.level;
  }
}

static void mesh_match_acoustic_interfaces() {
  if (Param.theAcousticInterfaceMethod != ACOUSTIC_INTERFACE_MATCHING) return;
  for (int pass = 0; ; ++pass) {
    mesh_exchange_matching_requests();
    int64_t count = g_acousticMatchingOwned.size(), total = 0;
    MPI_Allreduce(&count, &total, 1, MPI_INT64_T, MPI_SUM, comm_solver);
    if (Global.myID == 0)
      monitor_print("Local matching pass %d: %lld coarse interface leaves%s\n",
                    pass, (long long)total, total ? " to refine" : "; matching complete");
    if (total == 0) break;
    if (pass >= 128)
      solver_abort(__FUNCTION_NAME, NULL, "Local matching did not converge in 128 mesh passes\n");
    mesh_release_acoustic_interface_data();
    g_acousticNeighbors.clear();
    g_acousticNeighborRequests.clear();
    octor_deletemesh(Global.myMesh);
    Global.myMesh = NULL;
    if (octor_refinetree(Global.myOctree, acoustic_matching_refine, setrec) != 0 ||
        octor_balancetree(Global.myOctree, setrec, 0) != 0 ||
        mesh_partitiontree("acoustic_local_matching") != 0)
      solver_abort(__FUNCTION_NAME, NULL, "Local matching mesh refinement failed\n");
    Global.myMesh = octor_extractmesh(Global.myOctree, bldgs_nodesearch,
        pushdowns_nodesearch, bldgs_nodesearch_com, find_topoAirOct, topo_nodesearch);
    if (!Global.myMesh)
      solver_abort(__FUNCTION_NAME, NULL, "Local matching mesh extraction failed\n");
    // Re-query material: refining a mixed cell can expose a different interface.
    mesh_correct_properties(Global.theCVMEp);
  }
  g_acousticMatchingOwned.clear();
}
#endif
