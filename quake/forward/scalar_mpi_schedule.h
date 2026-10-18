#ifndef HERCULES_SCALAR_MPI_SCHEDULE_H
#define HERCULES_SCALAR_MPI_SCHEDULE_H

/* Private psolve.cc implementation. Build after final material/constraint and
 * receiver maps exist, before timestepping. Never infer participation from a
 * rank's local element count: a rank without cells may own another rank's DOF.
 * Keep the original schedules for initialization, vector fields and fallback.
 */
enum { SCALAR_MPI_ACOUSTIC, SCALAR_MPI_PORE, SCALAR_MPI_TARGET,
       SCALAR_MPI_FIELDS };
static schedule_t* g_scalarMpiSchedules[SCALAR_MPI_FIELDS][2] = {};
static int g_gpuMpiProfile = 0;
static double g_gpuMpiMetrics[7][11] = {};

static int scalar_mpi_field(int tag) {
  switch (tag) {
    case DN_PRESSURE_FORCE_MSG: case AN_PRESSURE_FORCE_MSG:
    case DN_PRESSURE_MSG: case AN_PRESSURE_MSG: return SCALAR_MPI_ACOUSTIC;
    case DN_PORE_PRESSURE_FORCE_MSG: case AN_PORE_PRESSURE_FORCE_MSG:
    case DN_PORE_PRESSURE_MSG: case AN_PORE_PRESSURE_MSG: return SCALAR_MPI_PORE;
    case DN_PORE_PRESSURE_OPEN_PORE_TARGET_MSG:
    case AN_PORE_PRESSURE_OPEN_PORE_TARGET_MSG: return SCALAR_MPI_TARGET;
    default: return -1;
  }
}

static schedule_t* scalar_mpi_schedule(schedule_t* original, int tag) {
  const int field = scalar_mpi_field(tag);
  if (field < 0) return original;
  const int which = original == Global.mySolver->dn_sched ? 0 :
                    original == Global.mySolver->an_sched ? 1 : -1;
  return which >= 0 && g_scalarMpiSchedules[field][which]
             ? g_scalarMpiSchedules[field][which] : original;
}

static void solver_free_scalar_schedules() {
  for (auto& pair : g_scalarMpiSchedules)
    for (auto& sched : pair) {
      if (sched) schedule_delete(sched);
      sched = NULL;
    }
}

/* Exchange bit masks on the ORIGINAL graph, OR at owners then broadcast to
 * sharers. Carry global IDs to verify that filtering preserves matched packing
 * order, including appended trace/receiver ghosts and duplicate entries.
 * All setup messages finish before their temporary storage is released.
 */
static int scalar_mpi_share_mask(schedule_t* sched, int direction,
                                  std::vector<uint8_t>& active) {
  struct entry_t { int64_t gid; uint8_t fields; };
  struct packet_t { messenger_t* messenger; std::vector<entry_t> entries; };
  std::vector<packet_t> sends, recvs;
  const auto packets = [&](messenger_t* first, std::vector<packet_t>& out) {
    for (auto* m = first; m; m = m->next) {
      if (m->nidx != m->nodecount || m->nodecount < 0 ||
          (size_t)m->nodecount > INT_MAX / sizeof(entry_t))
        solver_abort(__FUNCTION_NAME, NULL, "Invalid scalar MPI setup mapping\n");
      out.push_back({m, std::vector<entry_t>((size_t)m->nodecount)});
      for (int32_t i = 0; i < m->nodecount; ++i) {
        const int32_t n = m->mapping[i];
        out.back().entries[i] = {Global.myMesh->nodeTable[n].gnid, active[n]};
      }
    }
  };
  packets(direction == CONTRIBUTION ? sched->first_c : sched->first_s, sends);
  packets(direction == CONTRIBUTION ? sched->first_s : sched->first_c, recvs);
  std::vector<MPI_Request> requests(recvs.size() + sends.size());
  std::vector<MPI_Status> statuses(requests.size());
  size_t r = 0;
  for (auto& p : recvs)
    MPI_Irecv(p.entries.data(), (int)(p.entries.size() * sizeof(entry_t)),
              MPI_BYTE, p.messenger->procid, SCALAR_MPI_SETUP_MSG, comm_solver,
              &requests[r++]);
  for (auto& p : sends)
    MPI_Isend(p.entries.data(), (int)(p.entries.size() * sizeof(entry_t)),
              MPI_BYTE, p.messenger->procid, SCALAR_MPI_SETUP_MSG, comm_solver,
              &requests[r++]);
  if (!requests.empty())
    MPI_Waitall((int)requests.size(), requests.data(), statuses.data());
  int changed = 0;
  for (size_t k = 0; k < recvs.size(); ++k) {
    const auto& p = recvs[k];
    int bytes = 0;
    MPI_Get_count(&statuses[k], MPI_BYTE, &bytes);
    if (bytes != (int)(p.entries.size() * sizeof(entry_t)))
      solver_abort(__FUNCTION_NAME, NULL, "Scalar MPI peer count mismatch\n");
    for (size_t i = 0; i < p.entries.size(); ++i) {
      const int32_t n = p.messenger->mapping[i];
      if (p.entries[i].gid != Global.myMesh->nodeTable[n].gnid)
        solver_abort(__FUNCTION_NAME, NULL, "Scalar MPI peer node order mismatch\n");
      const uint8_t value = active[n] | p.entries[i].fields;
      changed |= value != active[n];
      active[n] = value;
    }
  }
  return changed;
}

/* Dependency closure is directed: an active dangling DOF needs its parents.
 * Activating every other child of those parents would unnecessarily flood the
 * graph. Inactive RHS entries are reset to zero by the unchanged solver.
 */
static int scalar_mpi_close_constraints(std::vector<uint8_t>& active) {
  int changed = 0;
  for (int f = 0; f < SCALAR_MPI_FIELDS; ++f) {
    const uint8_t bit = (uint8_t)(1u << f);
    const int field = f == SCALAR_MPI_ACOUSTIC ? ACOUSTIC_CONSTRAINT_SCALAR
                                              : ACOUSTIC_CONSTRAINT_PORE;
    const auto mark = [&](int32_t node, int32_t parent) {
      if ((active[node] & bit) && !(active[parent] & bit)) {
        active[parent] |= bit;
        changed = 1;
      }
    };
    const auto& map = g_acousticConstraintMaps[field];
    if (g_acousticConstraintsReady && map.enabled) {
      for (size_t d = 0; d < map.node_base3.size(); ++d)
        for (int32_t j = map.offsets[d]; j < map.offsets[d + 1]; ++j)
          mark(map.node_base3[d] / 3, map.parent_base3[j] / 3);
    } else {
      for (int32_t d = 0; d < Global.myMesh->ldnnum; ++d) {
        const auto& node = Global.myMesh->dnodeTable[d];
        for (auto* parent = node.lanid; parent; parent = parent->next)
          mark(node.ldnid, parent->id);
      }
    }
  }
  return changed;
}

static schedule_t* scalar_mpi_filter(schedule_t* original,
                                      const std::vector<uint8_t>& active,
                                      int field) {
  schedule_t* result = schedule_new();
  if (!result) solver_abort(__FUNCTION_NAME, NULL, "Scalar MPI schedule allocation\n");
  for (int side = 0; side < 2; ++side) {
    messenger_t** tail = side ? &result->first_s : &result->first_c;
    for (auto* src = side ? original->first_s : original->first_c;
         src; src = src->next) {
      std::vector<int32_t> nodes;
      for (int32_t i = 0; i < src->nidx; ++i)
        if (active[src->mapping[i]] & (1u << field)) nodes.push_back(src->mapping[i]);
      if (nodes.empty()) continue;
      auto* m = messenger_new(src->procid);
      if (!m) solver_abort(__FUNCTION_NAME, NULL, "Scalar MPI messenger allocation\n");
      m->nodecount = m->nidx = (int32_t)nodes.size();
      m->mapping = (int32_t*)malloc(nodes.size() * sizeof(int32_t));
      if (!m->mapping) solver_abort(__FUNCTION_NAME, NULL, "Scalar MPI node allocation\n");
      memcpy(m->mapping, nodes.data(), nodes.size() * sizeof(int32_t));
      *tail = m;
      tail = &m->next;
      if (side) { result->messenger_s[m->procid] = m; ++result->s_count; }
      else { result->messenger_c[m->procid] = m; ++result->c_count; }
    }
  }
  schedule_allocMPIctl(result);
  const int size = sizeof(solver_float) + (DO_DEBUG ? sizeof(int64_t) : 0);
  schedule_prepare(result, size, size, size, size);
  messenger_prepare_offsets(result->first_c, 1);
  messenger_prepare_offsets(result->first_s, 1);
  return result;
}

static void solver_init_scalar_schedules(mysolver_t* solver) {
  const char* mode = getenv("HERCULES_SCALAR_MPI_SCHEDULE");
  const char* profile = getenv("HERCULES_GPU_MPI_PROFILE");
  int options[2] = {1, 0}, lo[2], hi[2];
  if (mode && strcasecmp(mode, "full") == 0) options[0] = 0;
  else if (mode && strcasecmp(mode, "compact") != 0)
    solver_abort(__FUNCTION_NAME, NULL, "HERCULES_SCALAR_MPI_SCHEDULE must be compact or full\n");
  if (profile && strcmp(profile, "1") == 0) options[1] = 1;
  else if (profile && strcmp(profile, "0") != 0)
    solver_abort(__FUNCTION_NAME, NULL, "HERCULES_GPU_MPI_PROFILE must be 0 or 1\n");
  MPI_Allreduce(options, lo, 2, MPI_INT, MPI_MIN, comm_solver);
  MPI_Allreduce(options, hi, 2, MPI_INT, MPI_MAX, comm_solver);
  if (lo[0] != hi[0] || lo[1] != hi[1])
    solver_abort(__FUNCTION_NAME, NULL, "Scalar MPI schedule/profile settings differ between ranks\n");
  g_gpuMpiProfile = options[1];
  memset(g_gpuMpiMetrics, 0, sizeof(g_gpuMpiMetrics));
  if (Global.myID == 0)
    monitor_print("Scalar MPI schedules: %s; GPU MPI profile: %s\n",
                  options[0] ? "compact (field/constraint closure)" : "full",
                  g_gpuMpiProfile ? "on" : "off");
  if (!options[0] || (!acoustic_pressure_enabled() && !poroelastic_gpu_active())) return;

  mesh_t* mesh = Global.myMesh;
  std::vector<uint8_t> active((size_t)mesh->nharbored, 0);
  for (int32_t e = 0; e < mesh->lenum; ++e) {
    const int bits =
        ((acoustic_pressure_enabled() && element_is_acoustic((edata_t*)mesh->elemTable[e].data)) ? 1 : 0) |
        ((poroelastic_gpu_active() && g_poroelasticElementFlags[e]) ? 2 : 0);
    if (bits) for (int i = 0; i < 8; ++i) active[mesh->elemTable[e].lnid[i]] |= bits;
  }
  for (int32_t n = 0; n < mesh->nharbored; ++n) {
    if (solver->pressure_mass && solver->pressure_mass[n] > 0) active[n] |= 1;
    if (solver->pore_pressure_mass && solver->pore_pressure_mass[n] > 0) active[n] |= 2;
    if (solver->pore_pressure_open_pore_mass && solver->pore_pressure_open_pore_mass[n] != 0)
      active[n] |= 4;
  }
  if (acoustic_pressure_enabled()) {
    for (auto n : g_acousticReceiverNodes) active[n] |= 1;
    const int32_t* nodes = planes_acoustic_pressure_nodes();
    for (int32_t i = 0; i < planes_acoustic_pressure_node_count(); ++i)
      if (nodes[i] >= 0 && nodes[i] < mesh->nharbored) active[nodes[i]] |= 1;
  }
  for (int32_t s = 0; s < Param.myNumberOfStations; ++s) {
    const auto& station = Param.myStations[s];
    const int bits = (station.pressure_is_defined ? 1 : 0) |
                     (station.pore_pressure_is_defined ? 2 : 0);
    if (bits)
      for (int i = 0; i < 8; ++i) active[station.nodestointerpolate[i]] |= bits;
  }
  if (acoustic_poroelastic_open_pore_enabled())
    for (int64_t f = 0; f < Global.acousticInterfaceFaces; ++f) {
      if (!Global.acousticInterfaceAdjacentPoroelastic[f] ||
          Global.acousticInterfacePoroelasticBoundaryScales[f] <= 0) continue;
      const int32_t* nodes = acoustic_solid_trace_nodes(f);
      for (int i = 0; i < 8; ++i) active[nodes[i]] |= 4;
    }

  int any = 1;
  while (any) {
    int changed = scalar_mpi_close_constraints(active);
    for (auto* sched : {solver->dn_sched, solver->an_sched}) {
      changed |= scalar_mpi_share_mask(sched, CONTRIBUTION, active);
      changed |= scalar_mpi_share_mask(sched, SHARING, active);
    }
    MPI_Allreduce(&changed, &any, 1, MPI_INT, MPI_MAX, comm_solver);
  }
  const char* names[] = {"acoustic", "pore pressure", "open-pore target"};
  for (int f = 0; f < SCALAR_MPI_FIELDS; ++f) {
    if ((f == SCALAR_MPI_ACOUSTIC && !acoustic_pressure_enabled()) ||
        (f == SCALAR_MPI_PORE && !poroelastic_gpu_active()) ||
        (f == SCALAR_MPI_TARGET && !acoustic_poroelastic_open_pore_enabled())) continue;
    int64_t counts[4] = {}, total[4];
    for (int which = 0; which < 2; ++which) {
      auto* original = which ? solver->an_sched : solver->dn_sched;
      auto* compact = scalar_mpi_filter(original, active, f);
      g_scalarMpiSchedules[f][which] = compact;
      for (auto* m = original->first_c; m; m = m->next) counts[0] += m->nodecount;
      for (auto* m = compact->first_c; m; m = m->next) counts[1] += m->nodecount;
      counts[2] += original->c_count;
      counts[3] += compact->c_count;
    }
    MPI_Reduce(counts, total, 4, MPI_INT64_T, MPI_SUM, 0, comm_solver);
    if (Global.myID == 0)
      monitor_print("Scalar MPI %s: node entries %lld -> %lld; peer schedules %lld -> %lld "
                    "(global contributor totals, dangling+anchored; each shared edge counted once)\n",
                    names[f], (long long)total[0], (long long)total[1],
                    (long long)total[2], (long long)total[3]);
  }
}

/* Profiling adds CPU timestamps only when explicitly enabled. It adds no
 * synchronization: pack/unpack times include the existing stream waits and
 * pending work, so they are NOT isolated kernel or network timings.
 */
struct gpu_mpi_profile_t {
  double* metrics = NULL;
  double start = 0;
  gpu_mpi_profile_t(int tag, messenger_t* send, messenger_t* recv) {
    if (!g_gpuMpiProfile || (!send && !recv)) return;
    int field = -1;
    switch (tag) {
      case DN_FORCE_MSG: case AN_FORCE_MSG: field = 0; break;
      case DN_DISP_MSG: case AN_DISP_MSG: field = 1; break;
      case DN_PRESSURE_FORCE_MSG: case AN_PRESSURE_FORCE_MSG: field = 2; break;
      case DN_PRESSURE_MSG: case AN_PRESSURE_MSG: field = 3; break;
      case DN_PORE_PRESSURE_FORCE_MSG: case AN_PORE_PRESSURE_FORCE_MSG: field = 4; break;
      case DN_PORE_PRESSURE_MSG: case AN_PORE_PRESSURE_MSG: field = 5; break;
      case DN_PORE_PRESSURE_OPEN_PORE_TARGET_MSG:
      case AN_PORE_PRESSURE_OPEN_PORE_TARGET_MSG: field = 6; break;
    }
    if (field < 0) return;
    metrics = g_gpuMpiMetrics[field];
    metrics[0] += 1;
    const int width = field < 2 ? sizeof(fvector_t) : sizeof(solver_float);
    for (auto* m = send; m; m = m->next) { metrics[1] += (double)m->nodecount * width; metrics[3] += 1; }
    for (auto* m = recv; m; m = m->next) { metrics[2] += (double)m->nodecount * width; metrics[4] += 1; }
    start = MPI_Wtime();
  }
  void mark(int metric) {
    if (!metrics) return;
    const double now = MPI_Wtime();
    metrics[metric] += now - start;
    start = now;
  }
};

static void solver_report_gpu_mpi_profile() {
  if (!g_gpuMpiProfile) return;
  double sum[7][11], low[7][11], high[7][11];
  MPI_Reduce(g_gpuMpiMetrics, sum, 77, MPI_DOUBLE, MPI_SUM, 0, comm_solver);
  MPI_Reduce(g_gpuMpiMetrics, low, 77, MPI_DOUBLE, MPI_MIN, 0, comm_solver);
  MPI_Reduce(g_gpuMpiMetrics, high, 77, MPI_DOUBLE, MPI_MAX, 0, comm_solver);
  if (Global.myID != 0) return;
  const char* fields[] = {"solid force", "displacement", "acoustic RHS", "acoustic field",
                          "pore RHS", "pore field", "open-pore target"};
  const char* metrics[] = {"nonempty calls", "send bytes", "receive bytes", "send messages", "receive messages",
      "post receive seconds", "pack/ready seconds", "post send seconds", "wait receive seconds",
      "unpack/ready seconds", "wait send seconds"};
  monitor_print("GPU MPI profile: cumulative rank average/min/max, device helpers only; "
                "pack/unpack include stream waits and pending GPU work; scalar phases are nested in Compute Physics\n");
  for (int f = 0; f < 7; ++f) {
    if (sum[f][0] == 0) continue;
    for (int m = 0; m < 11; ++m)
      monitor_print("GPU MPI %s %s: %.9g / %.9g / %.9g\n", fields[f], metrics[m],
                    sum[f][m] / Global.theGroupSize, low[f][m], high[f][m]);
  }
}

#endif
