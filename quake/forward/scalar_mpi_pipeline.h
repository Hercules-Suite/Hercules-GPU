#ifndef HERCULES_SCALAR_MPI_PIPELINE_H
#define HERCULES_SCALAR_MPI_PIPELINE_H

/* Private psolve.cc transport. Plans and storage are built once. Each peer's
 * packet concatenates existing field mappings (no union-node padding), while
 * unpack visits the original field/peer order, including duplicate entries.
 * Separate phase buffers remain live until MPI and CUDA consumers finish. */
enum { SCALAR_BATCH_ACOUSTIC, SCALAR_BATCH_PORE, SCALAR_BATCH_SEALED,
       SCALAR_BATCH_OPEN, SCALAR_BATCH_TARGET, SCALAR_BATCH_GROUPS };
static int g_scalarMpiBatch = 1;
static int g_gpuMpiOverlap = 1;
static int g_scalarSealedBatch = 0;
static const int g_scalarBatchFields[SCALAR_BATCH_GROUPS][2] = {
    {0, -1}, {1, -1}, {0, 1}, {1, 2}, {2, -1}};

/* Waiting for this event only protects reuse of the receive storage. Consumers
 * of the unpacked field stay ordered on CUDA_STREAM_MAIN; no host-side wait is
 * required after unpack. Never use a stream event to infer MPI send completion. */
struct mpi_unpack_fence_t {
  cudaEvent_t event = NULL;
  bool pending = false;
  void init() {
    solver_pore_pressure_device_check(
        cudaEventCreateWithFlags(&event, cudaEventDisableTiming),
        "MPI receive reuse event");
  }
  void wait() {
    if (!pending) return;
    solver_pore_pressure_device_check(cudaEventSynchronize(event),
                                      "MPI receive storage reuse");
    pending = false;
  }
  void record(cudaStream_t stream) {
    solver_pore_pressure_device_check(cudaEventRecord(event, stream),
                                      "MPI unpack completion event");
    pending = true;
  }
  void destroy() {
    wait();
    if (event) cudaEventDestroy(event);
    event = NULL;
  }
};
static mpi_unpack_fence_t g_vectorUnpackFences[4];
static mpi_unpack_fence_t g_poreHostUploadFence;
static bool g_poreHostUploadTracked = false;

struct scalar_batch_packet_t {
  int peer = 0;
  int count = 0;
  solver_float* device = NULL;
  std::vector<solver_float> host;
};
struct scalar_batch_segment_t {
  messenger_t* mapping;
  int field, packet, offset;
};
struct scalar_batch_plan_t {
  bool enabled = false, device = false, active = false;
  int group = 0, phase = 0;
  std::vector<scalar_batch_packet_t> send, recv;
  std::vector<scalar_batch_segment_t> pack, unpack;
  std::vector<MPI_Request> send_requests, recv_requests;
  std::vector<MPI_Status> send_status, recv_status;
  mpi_unpack_fence_t unpack_fence;
  cudaEvent_t pack_ready = NULL;
  gpu_mpi_profile_t profile;
};
static scalar_batch_plan_t g_scalarBatchPlans[SCALAR_BATCH_GROUPS][4];
struct scalar_update_partition_t {
  int32_t count[2] = {};
  int32_t* nodes[2] = {};
};
static scalar_update_partition_t g_scalarUpdatePartitions[2];

static int scalar_batch_profile_field(int group, int phase) {
  if (group == SCALAR_BATCH_SEALED) return phase < 2 ? 7 : 8;
  if (group == SCALAR_BATCH_OPEN) return 9;
  if (group == SCALAR_BATCH_TARGET) return 6;
  return (group == SCALAR_BATCH_ACOUSTIC ? 2 : 4) + (phase >= 2);
}
static int scalar_batch_tag(int field, int phase) {
  static const int tags[3][4] = {
    {DN_PRESSURE_FORCE_MSG, AN_PRESSURE_FORCE_MSG, AN_PRESSURE_MSG, DN_PRESSURE_MSG},
    {DN_PORE_PRESSURE_FORCE_MSG, AN_PORE_PRESSURE_FORCE_MSG, AN_PORE_PRESSURE_MSG, DN_PORE_PRESSURE_MSG},
    {DN_PORE_PRESSURE_OPEN_PORE_TARGET_MSG, AN_PORE_PRESSURE_OPEN_PORE_TARGET_MSG, 0, 0}};
  return tags[field][phase];
}

/* Verify the final concatenated layout too, including full schedules and
 * debug builds. Runtime packets need no per-value IDs after this setup check. */
static void scalar_batch_verify_layout(const scalar_batch_plan_t& plan) {
  struct entry_t { int64_t gid; int field; };
  std::vector<std::vector<entry_t>> send(plan.send.size()), recv(plan.recv.size());
  for (size_t i = 0; i < send.size(); ++i) send[i].resize(plan.send[i].count);
  for (size_t i = 0; i < recv.size(); ++i) recv[i].resize(plan.recv[i].count);
  for (auto& s : plan.pack)
    for (int i = 0; i < s.mapping->nodecount; ++i)
      send[s.packet][s.offset + i] = {
          Global.myMesh->nodeTable[s.mapping->mapping[i]].gnid, s.field};
  std::vector<MPI_Request> requests(recv.size() + send.size());
  std::vector<MPI_Status> statuses(requests.size());
  size_t k = 0;
  for (size_t i = 0; i < recv.size(); ++i) {
    if (recv[i].size() > INT_MAX / sizeof(entry_t))
      solver_abort(__FUNCTION_NAME, NULL, "Scalar MPI setup packet exceeds MPI byte count\n");
    MPI_Irecv(recv[i].data(), (int)(recv[i].size() * sizeof(entry_t)), MPI_BYTE,
        plan.recv[i].peer, SCALAR_MPI_SETUP_MSG, comm_solver, &requests[k++]);
  }
  for (size_t i = 0; i < send.size(); ++i) {
    if (send[i].size() > INT_MAX / sizeof(entry_t))
      solver_abort(__FUNCTION_NAME, NULL, "Scalar MPI setup packet exceeds MPI byte count\n");
    MPI_Isend(send[i].data(), (int)(send[i].size() * sizeof(entry_t)), MPI_BYTE,
        plan.send[i].peer, SCALAR_MPI_SETUP_MSG, comm_solver, &requests[k++]);
  }
  if (k) MPI_Waitall((int)k, requests.data(), statuses.data());
  for (size_t i = 0; i < recv.size(); ++i) {
    int bytes = 0; MPI_Get_count(&statuses[i], MPI_BYTE, &bytes);
    if (bytes != (int)(recv[i].size() * sizeof(entry_t)))
      solver_abort(__FUNCTION_NAME, NULL, "Scalar MPI batch peer count mismatch\n");
  }
  for (auto& s : plan.unpack)
    for (int i = 0; i < s.mapping->nodecount; ++i) {
      const auto& e = recv[s.packet][s.offset + i];
      if (e.field != s.field || e.gid != Global.myMesh->nodeTable[s.mapping->mapping[i]].gnid)
        solver_abort(__FUNCTION_NAME, NULL, "Scalar MPI batch field/node order mismatch\n");
    }
}

static void scalar_batch_build(scalar_batch_plan_t& plan, int group, int phase,
                               bool device, mysolver_t* solver) {
  plan.enabled = true; plan.group = group; plan.phase = phase; plan.device = device;
  for (int f : g_scalarBatchFields[group]) {
    if (f < 0) continue;
    auto* original = (phase == 0 || phase == 3) ? solver->dn_sched : solver->an_sched;
    auto* sched = scalar_mpi_schedule(original, scalar_batch_tag(f, phase));
    for (int receiving = 0; receiving < 2; ++receiving) {
      auto& packets = receiving ? plan.recv : plan.send;
      auto& segments = receiving ? plan.unpack : plan.pack;
      const bool contributors = (phase < 2) != (receiving != 0);
      for (auto* m = contributors ? sched->first_c : sched->first_s; m; m = m->next) {
        if (m->nodecount <= 0) continue;
        size_t p = 0;
        while (p < packets.size() && packets[p].peer != m->procid) ++p;
        if (p == packets.size()) {
          packets.emplace_back(); packets.back().peer = m->procid;
        }
        auto& packet = packets[p];
        if (m->nodecount > INT_MAX / (int)sizeof(solver_float) - packet.count)
          solver_abort(__FUNCTION_NAME, NULL, "Scalar MPI batch exceeds MPI byte count\n");
        segments.push_back({m, f, (int)p, packet.count});
        packet.count += m->nodecount;
        if (device && !m->scalar_mapping_device)
          solver_pore_pressure_device_copy_int32(&m->scalar_mapping_device,
              m->mapping, m->nodecount, "Scalar MPI batch mapping");
      }
    }
  }
  scalar_batch_verify_layout(plan);
  for (auto* packets : {&plan.send, &plan.recv})
    for (auto& p : *packets) {
      if (device)
        solver_pore_pressure_device_check(cudaMalloc((void**)&p.device,
            (size_t)p.count * sizeof(solver_float)), "Scalar MPI batch buffer");
      else p.host.resize(p.count);
    }
  plan.send_requests.resize(plan.send.size()); plan.send_status.resize(plan.send.size());
  plan.recv_requests.resize(plan.recv.size()); plan.recv_status.resize(plan.recv.size());
  if (device && (!plan.send.empty() || !plan.recv.empty())) {
    plan.unpack_fence.init();
    solver_pore_pressure_device_check(cudaEventCreateWithFlags(&plan.pack_ready,
        cudaEventDisableTiming), "Scalar MPI pack event");
  }
}

static void scalar_batch_begin(scalar_batch_plan_t& plan,
                               solver_float* const values[3], mysolver_t* solver) {
  if (!plan.enabled || plan.active)
    solver_abort(__FUNCTION_NAME, NULL, "Invalid scalar MPI batch lifecycle\n");
  plan.active = true;
  plan.profile.metrics = NULL;
  if (plan.send.empty() && plan.recv.empty()) return;
  if (g_gpuMpiProfile) {
    plan.profile.metrics = g_gpuMpiMetrics[scalar_batch_profile_field(plan.group, plan.phase)];
    auto* m = plan.profile.metrics;
    m[0] += 1;
    for (auto& p : plan.send) { m[1] += (double)p.count * sizeof(solver_float); m[3] += 1; }
    for (auto& p : plan.recv) { m[2] += (double)p.count * sizeof(solver_float); m[4] += 1; }
    plan.profile.start = MPI_Wtime();
  }
  // The previous iteration's GPU unpack must release storage before MPI writes.
  if (plan.device) plan.unpack_fence.wait();
  for (size_t i = 0; i < plan.recv.size(); ++i) {
    auto& p = plan.recv[i];
    MPI_Irecv(plan.device ? p.device : p.host.data(), p.count * sizeof(solver_float),
        MPI_BYTE, p.peer, 150 + 4 * plan.group + plan.phase, comm_solver,
        &plan.recv_requests[i]);
  }
  plan.profile.mark(5);
  for (auto& s : plan.pack) {
    auto* m = s.mapping;
    auto& p = plan.send[s.packet];
    if (plan.device)
      psolve_launch_scalar_boundary_gather((m->nodecount + 255) / 256, 256, 0,
          solver->streams[CUDA_STREAM_MAIN], m->nodecount, m->scalar_mapping_device,
          values[s.field], p.device + s.offset);
    else
      for (int32_t i = 0; i < m->nodecount; ++i)
        p.host[s.offset + i] = values[s.field][m->mapping[i]];
  }
  if (plan.device && !plan.send.empty()) {
    solver_pore_pressure_device_check(cudaGetLastError(), "Scalar MPI batch pack");
    solver_pore_pressure_device_check(cudaEventRecord(plan.pack_ready,
        solver->streams[CUDA_STREAM_MAIN]), "Scalar MPI pack ready record");
    // CUDA-aware MPI is not assumed to understand CUDA stream ordering.
    solver_pore_pressure_device_check(cudaEventSynchronize(plan.pack_ready),
                                      "Scalar MPI pack ready wait");
  }
  plan.profile.mark(6);
  for (size_t i = 0; i < plan.send.size(); ++i) {
    auto& p = plan.send[i];
    MPI_Isend(plan.device ? p.device : p.host.data(), p.count * sizeof(solver_float),
        MPI_BYTE, p.peer, 150 + 4 * plan.group + plan.phase, comm_solver,
        &plan.send_requests[i]);
  }
  plan.profile.mark(7);
}

static void scalar_batch_finish(scalar_batch_plan_t& plan,
                                solver_float* const values[3], mysolver_t* solver) {
  if (!plan.active) solver_abort(__FUNCTION_NAME, NULL, "Scalar MPI finish without begin\n");
  if (!plan.recv.empty()) MPI_Waitall((int)plan.recv.size(),
      plan.recv_requests.data(), plan.recv_status.data());
  plan.profile.mark(8);
  // Preserve original per-field peer order, independent of arrival order.
  for (auto& s : plan.unpack) {
    auto* m = s.mapping;
    auto& p = plan.recv[s.packet];
    if (!plan.device) {
      for (int32_t i = 0; i < m->nodecount; ++i) {
        if (plan.phase < 2) values[s.field][m->mapping[i]] += p.host[s.offset + i];
        else values[s.field][m->mapping[i]] = p.host[s.offset + i];
      }
    } else if (plan.phase < 2)
      psolve_launch_scalar_boundary_accumulate((m->nodecount + 255) / 256, 256, 0,
          solver->streams[CUDA_STREAM_MAIN], m->nodecount, m->scalar_mapping_device,
          p.device + s.offset, values[s.field]);
    else if (m->mapping_has_duplicates)
      psolve_launch_scalar_boundary_scatter_serial(1, 1, 0,
          solver->streams[CUDA_STREAM_MAIN], m->nodecount, m->scalar_mapping_device,
          p.device + s.offset, values[s.field]);
    else
      psolve_launch_scalar_boundary_scatter((m->nodecount + 255) / 256, 256, 0,
          solver->streams[CUDA_STREAM_MAIN], m->nodecount, m->scalar_mapping_device,
          p.device + s.offset, values[s.field]);
  }
  if (plan.device && !plan.recv.empty()) {
    solver_pore_pressure_device_check(cudaGetLastError(), "Scalar MPI batch unpack");
    if (g_gpuMpiOverlap) plan.unpack_fence.record(solver->streams[CUDA_STREAM_MAIN]);
    else solver_pore_pressure_device_check(cudaStreamSynchronize(
        solver->streams[CUDA_STREAM_MAIN]), "Scalar MPI batch unpack wait");
  }
  plan.profile.mark(9);
  // Send storage cannot be repacked or freed until MPI releases it.
  if (!plan.send.empty()) MPI_Waitall((int)plan.send.size(),
      plan.send_requests.data(), plan.send_status.data());
  plan.profile.mark(10);
  plan.active = false;
}

static void scalar_partition_build(int field, int rhs_group, mysolver_t* solver) {
  auto& part = g_scalarUpdatePartitions[field];
  std::vector<uint8_t> boundary(Global.myMesh->nharbored, 0);
  // Dangling contribution + distribution are already complete. Only the
  // anchored receive destinations can still change during the interior update.
  // The pore update also reads the open-pore target, so include its destinations.
  for (auto& s : g_scalarBatchPlans[rhs_group][1].unpack) {
    if (s.field != field && !(field == SCALAR_MPI_PORE && s.field == SCALAR_MPI_TARGET)) continue;
    for (int i = 0; i < s.mapping->nodecount; ++i) boundary[s.mapping->mapping[i]] = 1;
  }
  std::vector<int32_t> nodes[2];
  for (int32_t n = 0; n < Global.myMesh->nharbored; ++n) nodes[boundary[n]].push_back(n);
  for (int k = 0; k < 2; ++k) {
    part.count[k] = (int32_t)nodes[k].size();
    solver_pore_pressure_device_copy_int32(&part.nodes[k], nodes[k].data(), part.count[k],
                                          "Scalar MPI update partition");
  }
  (void)solver;
}

static void solver_init_mpi_pipeline(mysolver_t* solver) {
  int opts[2] = {1, 1}, lo[2], hi[2];
  const char* names[] = {"HERCULES_SCALAR_MPI_BATCH", "HERCULES_GPU_MPI_OVERLAP"};
  for (int i = 0; i < 2; ++i) {
    const char* env = getenv(names[i]);
    if (env && strcmp(env, "0") == 0) opts[i] = 0;
    else if (env && strcmp(env, "1") != 0)
      solver_abort(__FUNCTION_NAME, NULL, "%s must be 0 or 1\n", names[i]);
  }
  MPI_Allreduce(opts, lo, 2, MPI_INT, MPI_MIN, comm_solver);
  MPI_Allreduce(opts, hi, 2, MPI_INT, MPI_MAX, comm_solver);
  if (lo[0] != hi[0] || lo[1] != hi[1])
    solver_abort(__FUNCTION_NAME, NULL, "Scalar MPI batch/overlap settings differ between ranks\n");
  g_scalarMpiBatch = opts[0]; g_gpuMpiOverlap = opts[1];
  const bool a = acoustic_pressure_enabled() && g_acousticPressureElementCount > 0;
  const bool p = poroelastic_gpu_active();
  const bool open = acoustic_poroelastic_open_pore_enabled();
  const bool ad = a && acoustic_pressure_device_sync_enabled();
  const bool pd = p && gpu_device_sync_enabled();
  g_scalarSealedBatch = g_scalarMpiBatch && a && p && !open && ad == pd;
  if (Global.myID == 0)
    monitor_print("Scalar MPI batching: %s; GPU MPI overlap: %s; sealed acoustic/pore batch: %s; "
                  "open-pore RHS/target batch: %s (host transport has no GPU update overlap)\n",
        g_scalarMpiBatch ? "on" : "off", g_gpuMpiOverlap ? "on" : "off",
        g_scalarSealedBatch ? "on" : "off", g_scalarMpiBatch && open ? "on" : "off");
  if (g_gpuMpiOverlap && gpu_device_sync_enabled())
    for (auto& fence : g_vectorUnpackFences) fence.init();
  if (!g_scalarMpiBatch && !g_gpuMpiOverlap) return;
  // Consecutive RHS/target uploads share pinned sparse staging storage.
  // CUDA stream ordering protects device storage, but the CPU must also wait
  // before rewriting a pinned source that an earlier H2D copy still reads.
  g_poreHostUploadTracked = p && open && !pd && pore_pressure_force_sync_mode() != 0;
  if (g_poreHostUploadTracked) g_poreHostUploadFence.init();
  for (int group = 0; group < SCALAR_BATCH_GROUPS; ++group)
    for (int phase = 0; phase < 4; ++phase) {
      const bool use = group == SCALAR_BATCH_ACOUSTIC ? a && (!g_scalarSealedBatch || g_poreDiagnosticDeferred) :
          group == SCALAR_BATCH_PORE ? p && !g_scalarSealedBatch && !(g_scalarMpiBatch && open && phase < 2) :
          group == SCALAR_BATCH_SEALED ? g_scalarSealedBatch :
          group == SCALAR_BATCH_OPEN ? g_scalarMpiBatch && open && phase < 2 :
          open && !g_scalarMpiBatch && phase < 2;
      if (use) scalar_batch_build(g_scalarBatchPlans[group][phase], group, phase,
          group == SCALAR_BATCH_ACOUSTIC ? ad : pd, solver);
    }
  if (g_gpuMpiOverlap && ad)
    scalar_partition_build(0, g_scalarSealedBatch ? SCALAR_BATCH_SEALED : SCALAR_BATCH_ACOUSTIC, solver);
  // With unbatched open-pore coupling, both target phases finish before the
  // pore RHS starts. Only pore RHS destinations remain in flight here.
  if (g_gpuMpiOverlap && pd)
    scalar_partition_build(1, g_scalarSealedBatch ? SCALAR_BATCH_SEALED :
        g_scalarMpiBatch && open ? SCALAR_BATCH_OPEN : SCALAR_BATCH_PORE, solver);
}

static void solver_free_mpi_pipeline() {
  for (auto& phases : g_scalarBatchPlans)
    for (auto& p : phases) {
      if (p.active) solver_abort(__FUNCTION_NAME, NULL, "Freeing active scalar MPI batch\n");
      p.unpack_fence.destroy();
      if (p.pack_ready) cudaEventDestroy(p.pack_ready);
      for (auto* packets : {&p.send, &p.recv})
        for (auto& packet : *packets) if (packet.device) cudaFree(packet.device);
      p = scalar_batch_plan_t{};
    }
  for (auto& p : g_scalarUpdatePartitions) {
    for (auto* nodes : p.nodes) if (nodes) cudaFree(nodes);
    p = scalar_update_partition_t{};
  }
  for (auto& f : g_vectorUnpackFences) f.destroy();
  g_poreHostUploadFence.destroy();
  g_poreHostUploadTracked = false;
}
#endif
