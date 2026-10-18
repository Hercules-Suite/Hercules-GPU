#ifndef HERCULES_PORE_DIAGNOSTIC_SCHEDULE_H
#define HERCULES_PORE_DIAGNOSTIC_SCHEDULE_H

/* Private psolve.cc policy. Only the linear undrained fast path has no
 * pressure feedback/history dependency. Fresh runs start at u=p=0, so its
 * diagnostic recurrence telescopes to -M^-1 Q^T u, including the unchanged
 * linear contribution/distribution/assignment operators. Restart/checkpoint
 * runs retain incremental pressure and both historical time levels. */
static bool g_poreDiagnosticDeferred = false;
static bool g_poreDiagnosticStations = false;
static bool g_poreDiagnosticDue = true;
static int64_t g_poreDiagnosticUpdates = 0;
static int64_t g_poreDiagnosticSkipped = 0;

static bool pore_diagnostic_due(int step, int total, int rate, bool stations) {
  // Keep both final pressure buffers current even without a final output.
  return step >= total - 2 || (stations && (step + 1) % rate == 0);
}

static void solver_init_pore_diagnostics() {
  const char* env = getenv("HERCULES_PORE_PRESSURE_UPDATE");
  int requested = 1, low, high;
  if (env && strcmp(env, "step") == 0) requested = 0;
  else if (env && strcmp(env, "output") != 0)
    solver_abort(__FUNCTION_NAME, NULL,
        "HERCULES_PORE_PRESSURE_UPDATE must be output or step\n");
  MPI_Allreduce(&requested, &low, 1, MPI_INT, MPI_MIN, comm_solver);
  MPI_Allreduce(&requested, &high, 1, MPI_INT, MPI_MAX, comm_solver);
  if (low != high)
    solver_abort(__FUNCTION_NAME, NULL, "Pore pressure update mode differs between ranks\n");
  g_poreDiagnosticDeferred = requested && poroelastic_gpu_active() &&
      poroelastic_undrained_fast_path_active() && !checkpoint_host_state_required();
  if (!poroelastic_gpu_active()) return;
  int local = 0, global = 0;
  if (Param.stationBinary != YES && station_pressure_value_requested())
    for (int i = 0; i < Param.myNumberOfStations; ++i)
      if (Param.myStations[i].pore_pressure_is_defined) { local = 1; break; }
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, comm_solver);
  g_poreDiagnosticStations = global != 0;
  if (Global.myID == 0)
    monitor_print("Pore pressure update mode: %s (request=%s; pore-pressure/Ru stations=%s); "
                  "checkpoint/restart and non-diagnostic physics retain every-step updates\n",
        g_poreDiagnosticDeferred ? "output" : "step", requested ? "output" : "step",
        global ? "yes" : "no");
}
#endif
