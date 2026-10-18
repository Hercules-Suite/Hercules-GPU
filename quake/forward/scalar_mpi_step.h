#ifndef HERCULES_SCALAR_MPI_STEP_H
#define HERCULES_SCALAR_MPI_STEP_H

/* Private psolve.cc step orchestration. Both sealed RHS assemblies read the
 * already-shared displacement and old scalar fields. Open-pore target assembly
 * instead consumes the NEW acoustic field, so its group must follow acoustic
 * update/sharing. No physical operator or integration formula changes here. */
static void scalar_group_values(mysolver_t* s, bool device, bool rhs,
                                solver_float* values[3]) {
  values[0] = device ? (rhs ? s->gpuData->pressureForceDevice : s->gpuData->pressureTm2Device)
                     : (rhs ? s->pressure_force : s->pressure_tm2);
  values[1] = device ? (rhs ? s->gpuData->porePressureForceDevice : s->gpuData->porePressureTm2Device)
                     : (rhs ? s->pore_pressure_force : s->pore_pressure_tm2);
  values[2] = device ? s->gpuData->porePressureOpenPoreTargetDevice : s->pore_pressure_open_pore_target;
}

static void scalar_group_adjust(int group, bool device, bool rhs,
                                solver_float* const values[3], mysolver_t* s) {
  for (int field : g_scalarBatchFields[group]) {
    if (field < 0) continue;
    const int constraint = field == 0 ? ACOUSTIC_CONSTRAINT_SCALAR : ACOUSTIC_CONSTRAINT_PORE;
    if (device) solver_scalar_device_adjust(s, values[field], rhs ? DISTRIBUTION : ASSIGNMENT, constraint);
    else compute_adjust(values[field], 1, rhs ? DISTRIBUTION : ASSIGNMENT, constraint);
  }
}

static void scalar_update_nodes(mysolver_t* s, int field, int part) {
  auto& p = g_scalarUpdatePartitions[field];
  const int32_t count = p.count[part];
  if (!count) return;
  if (field == 0)
    psolve_launch_acoustic_pressure_update_nodes((count + 255) / 256, 256, 0,
        s->streams[CUDA_STREAM_MAIN], count, s->gpuDataDevice, p.nodes[part]);
  else if (g_poreDiagnosticDeferred)
    psolve_launch_poroelastic_diagnostic_update((count + 255) / 256, 256, 0,
        s->streams[CUDA_STREAM_MAIN], count, s->gpuDataDevice, p.nodes[part]);
  else
    psolve_launch_poroelastic_pressure_update_nodes((count + 255) / 256, 256, 0,
        s->streams[CUDA_STREAM_MAIN], count, s->gpuDataDevice,
        Param.thePoroelasticSurfacePressureBC == POROELASTIC_SURFACE_PRESSURE_DRAINED &&
            Param.thePoroelasticSurfacePressurePenaltyFactor > 0.0,
        acoustic_poroelastic_open_pore_enabled() &&
            s->gpuData->porePressureOpenPoreMassDevice && s->gpuData->porePressureOpenPoreTargetDevice,
        p.nodes[part]);
  solver_pore_pressure_device_check(cudaGetLastError(), "Scalar MPI partitioned update");
}

static void scalar_run_group(mysolver_t* s, mesh_t* mesh, int group,
                             bool acoustic, bool pore, bool open, bool device) {
  Timer_Start("Scalar Pressure Pipeline");
  if (acoustic) scalar_acoustic_prepare(s, mesh, device);
  if (pore) scalar_pore_prepare(s, mesh, device);
  if (open) {
    scalar_target_prepare(s, mesh, device);
    if (!g_scalarMpiBatch) {
      // Complete the independent target first so the pore interior update can
      // overlap its RHS exchange even when packet batching is disabled.
      solver_float* values[3]; scalar_group_values(s, device, true, values);
      auto& dn = g_scalarBatchPlans[SCALAR_BATCH_TARGET][0];
      auto& an = g_scalarBatchPlans[SCALAR_BATCH_TARGET][1];
      scalar_batch_begin(dn, values, s); scalar_batch_finish(dn, values, s);
      scalar_group_adjust(SCALAR_BATCH_TARGET, device, true, values, s);
      scalar_batch_begin(an, values, s); scalar_batch_finish(an, values, s);
      scalar_target_upload(s, mesh, device);
    }
  }
  solver_float* values[3]; scalar_group_values(s, device, true, values);
  Timer_Start("Scalar MPI RHS exchange");
  auto& dn = g_scalarBatchPlans[group][0];
  auto& an = g_scalarBatchPlans[group][1];
  scalar_batch_begin(dn, values, s); scalar_batch_finish(dn, values, s);
  scalar_group_adjust(group, device, true, values, s);
  scalar_batch_begin(an, values, s);
  const bool overlap = device && g_gpuMpiOverlap;
  if (overlap) {
    Timer_Start("Scalar interior update enqueue");
    if (acoustic) scalar_update_nodes(s, 0, 0);
    if (pore) scalar_update_nodes(s, 1, 0);
    Timer_Stop("Scalar interior update enqueue");
  }
  scalar_batch_finish(an, values, s);
  Timer_Stop("Scalar MPI RHS exchange");
  if (acoustic) scalar_acoustic_rhs_upload(s, mesh, device);
  if (pore) scalar_pore_rhs_upload(s, mesh, device);
  if (open && g_scalarMpiBatch) scalar_target_upload(s, mesh, device);
  if (overlap) {
    Timer_Start("Scalar boundary update enqueue");
    if (acoustic) scalar_update_nodes(s, 0, 1);
    if (pore) scalar_update_nodes(s, 1, 1);
    Timer_Stop("Scalar boundary update enqueue");
  } else {
    if (acoustic) scalar_acoustic_update(s, mesh, device);
    if (pore) scalar_pore_update(s, mesh, device);
  }
  if (acoustic) scalar_acoustic_field_download(s, mesh, device);
  if (pore) scalar_pore_field_download(s, mesh, device);
  scalar_group_values(s, device, false, values);
  const int field_group = group == SCALAR_BATCH_OPEN ? SCALAR_BATCH_PORE : group;
  Timer_Start("Scalar MPI field exchange");
  auto& anchored = g_scalarBatchPlans[field_group][2];
  auto& dangling = g_scalarBatchPlans[field_group][3];
  scalar_batch_begin(anchored, values, s); scalar_batch_finish(anchored, values, s);
  scalar_group_adjust(field_group, device, false, values, s);
  scalar_batch_begin(dangling, values, s); scalar_batch_finish(dangling, values, s);
  Timer_Stop("Scalar MPI field exchange");
  if (acoustic) scalar_acoustic_field_upload(s, mesh, device);
  if (pore) scalar_pore_field_upload(s, mesh, device);
  Timer_Stop("Scalar Pressure Pipeline");
}

static void solver_run_acoustic_pressure_step(mysolver_t* s, mesh_t* mesh) {
  if (!acoustic_pressure_enabled() || g_acousticPressureElementCount <= 0) return;
  const int device = acoustic_pressure_device_sync_enabled();
  if (g_scalarMpiBatch || g_gpuMpiOverlap) {
    scalar_run_group(s, mesh, SCALAR_BATCH_ACOUSTIC, true, false, false, device);
    return;
  }
  Timer_Start("Acoustic Pressure Step");
  scalar_acoustic_prepare(s, mesh, device);
  scalar_acoustic_rhs_exchange(s, mesh, device);
  scalar_acoustic_rhs_upload(s, mesh, device);
  scalar_acoustic_update(s, mesh, device);
  scalar_acoustic_field_download(s, mesh, device);
  scalar_acoustic_field_exchange(s, mesh, device);
  scalar_acoustic_field_upload(s, mesh, device);
  Timer_Stop("Acoustic Pressure Step");
}

static void solver_run_poroelastic_pressure_step(mysolver_t* s, mesh_t* mesh) {
  if (!poroelastic_gpu_active()) return;
  const int device = gpu_device_sync_enabled();
  const bool open = acoustic_poroelastic_open_pore_enabled();
  if (g_scalarMpiBatch || g_gpuMpiOverlap) {
    scalar_run_group(s, mesh, open && g_scalarMpiBatch ? SCALAR_BATCH_OPEN : SCALAR_BATCH_PORE,
                     false, true, open, device);
    return;
  }
  Timer_Start("Poroelastic Pressure Step");
  scalar_pore_prepare(s, mesh, device);
  scalar_pore_rhs_exchange(s, mesh, device);
  scalar_pore_rhs_upload(s, mesh, device);
  if (open) {
    scalar_target_prepare(s, mesh, device);
    scalar_target_exchange(s, mesh, device);
    scalar_target_upload(s, mesh, device);
  }
  scalar_pore_update(s, mesh, device);
  scalar_pore_field_download(s, mesh, device);
  scalar_pore_field_exchange(s, mesh, device);
  scalar_pore_field_upload(s, mesh, device);
  Timer_Stop("Poroelastic Pressure Step");
}

static void solver_run_scalar_pressure_steps(mysolver_t* s, mesh_t* mesh) {
  if (g_poreDiagnosticDeferred && !g_poreDiagnosticDue) {
    solver_run_acoustic_pressure_step(s, mesh);
    return;
  }
  if (g_scalarSealedBatch)
    scalar_run_group(s, mesh, SCALAR_BATCH_SEALED, true, true, false,
                     acoustic_pressure_device_sync_enabled());
  else {
    solver_run_acoustic_pressure_step(s, mesh);
    solver_run_poroelastic_pressure_step(s, mesh);
  }
}
#endif
