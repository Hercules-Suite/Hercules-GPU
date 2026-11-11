/* Included by psolve.cc after its parameter, mesh, and monitor declarations.
 * All work is confined to completed-mesh setup and solver finalization. */
#ifndef HERCULES_CPML_PARTITION_AUTO_H
#define HERCULES_CPML_PARTITION_AUTO_H

static std::string cpml_cache_file_identity(const char* path) {
  if (!path || !*path) return "";
  char resolved[PATH_MAX];
  const char* name = realpath(path, resolved) ? resolved : path;
  std::ostringstream out;
  out << name;
  struct stat info;
  if (stat(path, &info) == 0) {
    out << ':' << info.st_size << ':' << info.st_mtime;
#ifdef __APPLE__
    out << ':' << info.st_mtimespec.tv_nsec;
#else
    out << ':' << info.st_mtim.tv_nsec;
#endif
  } else {
    out << ":unavailable";
  }
  return out.str();
}

/* Keep physics settings and referenced-file metadata in the key. Discard
 * comments, run-specific destinations and the requested partition weight.
 * Linear sources may share a calibration. Nonlinear learning is disabled:
 * yielding/substep work can depend on the evolving source-driven state. */
static bool cpml_cache_configuration(std::string& result) {
  std::ifstream file(Param.parameters_input_file);
  if (!file) return false;
  std::ostringstream out;
  std::string line;
  while (std::getline(file, line)) {
    line = line.substr(0, line.find('#'));
    if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;
    const size_t equal = line.find('=');
    if (equal != std::string::npos) {
      std::istringstream left(line.substr(0, equal));
      std::string key;
      left >> key;
      const bool destination = key == "monitor_file" ||
          key == "planes_input_file" || key == "mesh_etree_output_file" ||
          key == "mesh_coordinates_directory_for_matlab" ||
          (key.find("output") != std::string::npos &&
           (key.find("file") != std::string::npos ||
            key.find("directory") != std::string::npos));
      if (destination || key == "cpml_partition_weight" ||
          key == "cpml_partition_cache" || key == "mesh_only_mode" ||
          key == "output_mesh" || key == "source_directory_output" ||
          (key == "source_directory" && Param.includeNonlinearAnalysis != YES))
        continue;
      // Referenced file metadata detects replaced material/geometry inputs.
      if (key.find("file") != std::string::npos ||
          key.find("directory") != std::string::npos) {
        std::istringstream right(line.substr(equal + 1));
        std::string path;
        right >> path;
        out << "file:" << cpml_cache_file_identity(path.c_str()) << '\n';
      }
    }
    std::istringstream tokens(line);
    std::string token;
    while (tokens >> token) out << token << ' ';
    out << '\n';
  }
  if (!file.eof()) return false;
  if (Param.theMeshInputMode)
    out << "archive:" << cpml_cache_file_identity(Param.mesh_etree_output_file) << '\n';
  out << "build:" << hercules_build_commit() << ':' << sizeof(solver_float) << '\n';
#ifdef __linux__
  // Include the executable identity so rebuilds with site-specific compiler
  // options cannot reuse timings solely because the Git commit is unchanged.
  out << "executable:" << cpml_cache_file_identity("/proc/self/exe") << '\n';
#else
  out << "compiled:" << __DATE__ << ':' << __TIME__ << '\n';
#endif
  out << "solver-ranks:" << Global.theGroupSize << '\n';
  result = out.str();
  return true;
}

static void mesh_prepare_cpml_partition_auto(void) {
  if (g_cpmlPartitionPrepared || Param.theAbcType != ABC_CPML ||
      Param.theMeshOnlyMode ||
      (!g_cpmlPartitionAutomatic && !Param.theCpmlPartitionCache[0]))
    return;
  g_cpmlPartitionPrepared = true;
  if (Param.includeNonlinearAnalysis == YES) {
    if (Global.myID == 0)
      monitor_print("CPML cost cache disabled for nonlinear analysis: "
                    "source-dependent yielding/substep work needs manual tuning; "
                    "selected weight %.6g.\n", Param.theCpmlPartitionWeight);
    return;
  }

  // A shared mesh path makes the default cache reusable across source run dirs.
  if (!Param.theCpmlPartitionCache[0]) {
    const int n = snprintf(Param.theCpmlPartitionCache,
        sizeof(Param.theCpmlPartitionCache), "%s.cpml-cost",
        Param.mesh_etree_output_file[0] ? Param.mesh_etree_output_file : "hercules");
    if (n < 0 || n >= (int)sizeof(Param.theCpmlPartitionCache))
      solver_abort(__FUNCTION_NAME, NULL,
                   "Automatic CPML cache path is too long; set cpml_partition_cache.\n");
  }

  int devices = 0;
  const int device = gpu_device_index(Global.myID, &devices);
  cudaDeviceProp props;
  memset(&props, 0, sizeof(props));
  const int device_ok = cudaGetDeviceProperties(&props, device) == cudaSuccess;
  MPI_Comm local_comm;
  MPI_Comm_split_type(comm_solver, MPI_COMM_TYPE_SHARED, Global.myID,
                      MPI_INFO_NULL, &local_comm);
  int local_size = 0;
  MPI_Comm_size(local_comm, &local_size);
  MPI_Comm_free(&local_comm);
  std::ostringstream hardware;
  int driver = 0, runtime = 0;
  cudaDriverGetVersion(&driver);
  cudaRuntimeGetVersion(&runtime);
  char mpi_version[MPI_MAX_LIBRARY_VERSION_STRING] = {};
  int mpi_version_length = 0;
  MPI_Get_library_version(mpi_version, &mpi_version_length);
  hardware << props.name << ':' << props.major << ':' << props.minor << ':'
           << props.multiProcessorCount << ':' << props.totalGlobalMem << ':'
           << CUDART_VERSION << ':' << driver << ':' << runtime << ':' << local_size
           << ':' << gpu_device_sync_enabled() << ':' << mpi_version;
  // Include performance overrides, but never rank IDs, host names, or job IDs.
  const char* const controls[] = {
      "HERCULES_GPU_MPI_BACKEND", "HERCULES_PORE_PRESSURE_SYNC_BACKEND",
      "HERCULES_BKT2_FUSED", "HERCULES_FORCE_GENERIC_DAMPING",
      "HERCULES_DISP_SYNC_MODE", "HERCULES_PRESSURE_FORCE_SYNC_MODE",
      "HERCULES_PRESSURE_FIELD_SYNC_MODE", "HERCULES_PORE_PRESSURE_FORCE_SYNC_MODE",
      "HERCULES_PORE_PRESSURE_FIELD_SYNC_MODE",
      "HERCULES_TOPO_FORCE_MODE", "HERCULES_SPARSE_CLEAR_MODE",
      "HERCULES_ADJUST_DIST_OMP", "HERCULES_ADJUST_ASSIGN_OMP",
      "HERCULES_FULL_CLEAR_OMP", "HERCULES_CONTRIB_PACK_OMP",
      "HERCULES_CONTRIB_UNPACK_OMP", "HERCULES_CONTRIB_RECV_MODE",
      "HERCULES_SPARSE_UPLOAD_MODE", "HERCULES_REPRO_MODE",
      "OMP_NUM_THREADS", "OMP_PROC_BIND", "OMP_PLACES",
      "HERCULES_UNIFIED_PHASE_B_OMP_THREADS", "OMPI_MCA_pml", "UCX_TLS"};
  for (const char* name : controls) {
    const char* value = getenv(name);
    hardware << '\n' << name << '=' << (value ? value : "<default>");
  }
  uint64_t local_key = cpml_partition::hash(hardware.str());
  uint64_t min_key = 0, max_key = 0;
  int local_ok = device_ok && devices >= local_size;
  int all_ok = 0;
  MPI_Allreduce(&local_ok, &all_ok, 1, MPI_INT, MPI_MIN, comm_solver);
  MPI_Allreduce(&local_key, &min_key, 1, MPI_UINT64_T, MPI_MIN, comm_solver);
  MPI_Allreduce(&local_key, &max_key, 1, MPI_UINT64_T, MPI_MAX, comm_solver);

  int ready = 0;
  if (Global.myID == 0) {
    std::string config;
    ready = all_ok && min_key == max_key &&
        strcmp(hercules_build_commit(), "unknown") != 0 &&
        strcmp(hercules_build_tree_state(), "clean") == 0 &&
        cpml_cache_configuration(config);
    if (ready) {
      // No leaf traversal or mesh payload transfer is needed for the key.
      g_cpmlPartitionCacheKey = cpml_partition::hash(config, local_key);
      const bool found = g_cpmlPartitionAutomatic && cpml_partition::load(
          Param.theCpmlPartitionCache, g_cpmlPartitionCacheKey,
          Param.theCpmlPartitionWeight);
      monitor_print("CPML cost cache: %s; key %016llx; %s; selected weight %.6g\n",
          Param.theCpmlPartitionCache,
          (unsigned long long)g_cpmlPartitionCacheKey,
          !g_cpmlPartitionAutomatic ? "numeric override (collect only)" :
          found ? "matching calibration" : "cache miss/invalid; fallback 1",
          Param.theCpmlPartitionWeight);
    } else {
      monitor_print("CPML cost cache disabled: requires a clean identified build, "
                    "readable input and homogeneous non-oversubscribed GPU layout; "
                    "selected weight %.6g.\n", Param.theCpmlPartitionWeight);
    }
  }
  MPI_Bcast(&ready, 1, MPI_INT, 0, comm_solver);
  MPI_Bcast(&Param.theCpmlPartitionWeight, 1, MPI_DOUBLE, 0, comm_solver);
  g_cpmlPartitionCacheReady = ready != 0;
}

static void solver_learn_cpml_partition_weight(void) {
  if (!g_cpmlPartitionCacheReady) return;
  // Compute Physics already covers elastic/BKT, acoustic scalar and coupling,
  // pore-pressure work. Nonlinear learning is disabled at setup because its
  // evolving state changes the work. No new per-step timer is added.
  const double seconds = Timer_Exists("Compute Physics")
      ? Timer_Value("Compute Physics", (TimerKind)0) : 0.0;
  std::vector<double> times(Global.myID == 0 ? Global.theGroupSize : 0);
  MPI_Gather(&seconds, 1, MPI_DOUBLE, times.data(), 1, MPI_DOUBLE, 0, comm_solver);
  if (Global.myID != 0) return;
  if (Param.theTotalSteps < 100) {
    monitor_print("CPML cost calibration skipped: fewer than 100 timesteps.\n");
    return;
  }
  std::vector<cpml_partition::Sample> samples;
  for (int rank = 0; rank < g_cpmlRankDiagSize; ++rank) {
    const long long* counts = g_cpmlRankDiagCounts + rank * CPML_RANK_DIAG_COUNT_STRIDE;
    samples.push_back({(double)counts[CPML_RANK_DIAG_LOCAL_ELEMS],
                       (double)counts[CPML_RANK_DIAG_TOTAL_ELEMS], times[rank]});
  }
  const cpml_partition::Estimate fit = cpml_partition::estimate(samples);
  if (!fit.valid) {
    monitor_print("CPML cost calibration rejected: %s; samples=%zu, "
                  "relative RMS error=%.4g. Existing cache retained.\n",
                  fit.reason, fit.samples, fit.relative_error);
    return;
  }
  const bool saved = cpml_partition::save(Param.theCpmlPartitionCache,
                                         g_cpmlPartitionCacheKey, fit);
  monitor_print("CPML cost calibration: ordinary %.9g, PML %.9g "
                "seconds/element over this run; estimated weight %.6g "
                "(auto range 1..16); relative RMS error %.4g; %zu ranks; %s.\n",
                fit.ordinary_cost, fit.pml_cost, fit.weight, fit.relative_error,
                fit.samples, saved ? "cache saved for future runs" :
                "cache write failed; simulation results unaffected");
}

#endif
