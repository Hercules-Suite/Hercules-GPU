#ifndef HERCULES_BIOT_RUNTIME_H
#define HERCULES_BIOT_RUNTIME_H

static void biot_surface_finish_update(mysolver_t* s) {
  if (!s->gpuData->biotDrainDiagnosticsDevice) return;
  double local[6]={},global[6]={};
  solver_pore_pressure_device_check(cudaMemcpyAsync(local,
      s->gpuData->biotDrainDiagnosticsDevice,sizeof(local),cudaMemcpyDeviceToHost,
      s->streams[CUDA_STREAM_MAIN]),"Biot surface diagnostics download");
  solver_pore_pressure_device_check(cudaStreamSynchronize(s->streams[CUDA_STREAM_MAIN]),
      "Biot surface diagnostics wait");
  MPI_Reduce(local,global,6,MPI_DOUBLE,MPI_SUM,0,comm_solver);
  ++g_biotSurfaceStep;
  if (Global.myID==0) {
    fprintf(g_biotSurfaceBalanceFile,"%.12e",Param.theStartT+g_biotSurfaceStep*Param.theDeltaT);
    for (int j=0;j<6;++j) fprintf(g_biotSurfaceBalanceFile," %.12e",global[j]);
    fputc('\n',g_biotSurfaceBalanceFile);
  }
}

/* Private psolve.cc vector communication; mirrors the solid constraint map.
 * Exchanges are collective across the solver, including ranks with no pores.
 * Use a disjoint tag range from scalar batches (150..169). */
static void biot_exchange(mysolver_t* s, mesh_t* mesh, bool force) {
  if (!full_biot_enabled()) return;
  fvector_t* device = force ? s->gpuData->biotForceDevice : s->gpuData->biotW2Device;
  schedule_t* first = force ? s->dn_sched : s->an_sched;
  schedule_t* last = force ? s->an_sched : s->dn_sched;
  const int direction = force ? CONTRIBUTION : SHARING;
  const int adjustment = force ? DISTRIBUTION : ASSIGNMENT;
  const int tag = force ? 180 : 182;
  if (gpu_device_sync_enabled()) {
    solver_vector_device_schedule_senddata(s, first, device, direction, tag);
    solver_vector_device_adjust(s, device, adjustment);
    solver_vector_device_schedule_senddata(s, last, device, direction, tag+1);
  } else {
    const size_t bytes = mesh->nharbored * sizeof(fvector_t);
    fvector_t* host = g_biotHostScratch.data();
    solver_pore_pressure_device_check(cudaMemcpyAsync(host, device, bytes,
        cudaMemcpyDeviceToHost, s->streams[CUDA_STREAM_MAIN]), "Biot vector download");
    cudaStreamSynchronize(s->streams[CUDA_STREAM_MAIN]);
    schedule_senddata(first, host, 3, direction, tag);
    compute_adjust(host, 3, adjustment, ACOUSTIC_CONSTRAINT_SOLID);
    schedule_senddata(last, host, 3, direction, tag+1);
    solver_pore_pressure_device_check(cudaMemcpyAsync(device, host, bytes,
        cudaMemcpyHostToDevice, s->streams[CUDA_STREAM_MAIN]), "Biot vector upload");
  }
}

static void biot_stabilization_exchange(mysolver_t* s, mesh_t* mesh, bool rhs) {
  solver_float* device=s->gpuData->biotStabilizationWorkDevice;
  schedule_t* first=rhs?s->dn_sched:s->an_sched;
  schedule_t* last=rhs?s->an_sched:s->dn_sched;
  const int direction=rhs?CONTRIBUTION:SHARING, adjustment=rhs?DISTRIBUTION:ASSIGNMENT;
  const int first_tag=rhs?DN_PORE_PRESSURE_FORCE_MSG:AN_PORE_PRESSURE_MSG;
  const int last_tag=rhs?AN_PORE_PRESSURE_FORCE_MSG:DN_PORE_PRESSURE_MSG;
  /* Reuse the pore schedule tags only after completing the preceding stage.
   * This also selects the existing compact pore maps on both backends. */
  if (gpu_device_sync_enabled()) {
    solver_scalar_device_schedule_senddata(s,first,device,direction,first_tag);
    solver_scalar_device_adjust(s,device,adjustment,ACOUSTIC_CONSTRAINT_PORE);
    solver_scalar_device_schedule_senddata(s,last,device,direction,last_tag);
  } else {
    const size_t bytes=mesh->nharbored*sizeof(solver_float);
    solver_float* host=g_biotScalarScratch.data();
    cudaMemcpyAsync(host,device,bytes,cudaMemcpyDeviceToHost,s->streams[CUDA_STREAM_MAIN]);
    cudaStreamSynchronize(s->streams[CUDA_STREAM_MAIN]);
    schedule_senddata(first,host,1,direction,first_tag);
    compute_adjust(host,1,adjustment,ACOUSTIC_CONSTRAINT_PORE);
    schedule_senddata(last,host,1,direction,last_tag);
    cudaMemcpyAsync(device,host,bytes,cudaMemcpyHostToDevice,s->streams[CUDA_STREAM_MAIN]);
  }
}

static void biot_stabilization_prepare(mysolver_t* s,mesh_t* mesh) {
  Timer_Start("Biot stabilization RHS");
  const cudaStream_t stream=s->streams[CUDA_STREAM_MAIN];
  solver_float* work=s->gpuData->biotStabilizationWorkDevice;
  cudaMemsetAsync(work,0,mesh->nharbored*sizeof(solver_float),stream);
  psolve_launch_biot_stabilization(stream,g_poroelasticLocalElementCount,
      s->gpuDataDevice,s->gpuData->porePressureTm1Device,work,1);
  biot_stabilization_exchange(s,mesh,true);
  Timer_Stop("Biot stabilization RHS");
}
static void biot_stabilize_pressure(mysolver_t* s,mesh_t* mesh) {
  if(!biot_stabilization_active()) return;
  Timer_Start("Biot pressure stabilization");
  const auto stream=s->streams[CUDA_STREAM_MAIN];
  solver_float* work=s->gpuData->biotStabilizationWorkDevice;
  if(g_biotMassStabilizationReady) g_biotMassStabilizationReady=false;
  else {
    biot_stabilization_prepare(s,mesh);
    if(Param.thePoroelasticMassMatrix) biot_mass_pressure(s,mesh,work,false);
    else psolve_launch_biot_stabilization_normalize(stream,mesh->nharbored,s->gpuDataDevice);
  }
  biot_stabilization_exchange(s,mesh,false);
  psolve_launch_biot_stabilization(stream,g_poroelasticLocalElementCount,
      s->gpuDataDevice,work,s->gpuData->porePressureForceDevice,-Param.theDeltaT);
  solver_pore_pressure_device_check(cudaGetLastError(),"Biot pressure stabilization");
  Timer_Stop("Biot pressure stabilization");
}

static std::vector<FILE*> g_biotStationFiles;
static std::vector<double> g_biotStationValues;
static int32_t* g_biotStationNodesDevice = NULL;
static double* g_biotStationCoordsDevice = NULL;
static double* g_biotStationValuesDevice = NULL;
static bool g_biotStationsReady = false;

static void biot_output_stations(int step) {
  if (!full_biot_enabled()) return;
  mysolver_t* s = Global.mySolver;
  if (!g_biotStationsReady) {
    std::vector<int32_t> nodes;
    std::vector<double> coords;
    for (int i=0; i<Param.myNumberOfStations; ++i) {
      const station_t& st = Param.myStations[i];
      if (!st.pore_pressure_is_defined) continue;
      nodes.insert(nodes.end(), st.nodestointerpolate, st.nodestointerpolate+8);
      double weights[8]; biot_station_weights(i,weights);
      coords.insert(coords.end(),weights,weights+8);
      char path[PATH_MAX];
      snprintf(path, sizeof(path), "%s/station.fluid.%d", Param.theStationsDirOut, st.id);
      FILE* file = hu_fopen(path, "w");
      fprintf(file, "# Full Biot: w=porosity*(Uf-u), q=dw/dt, vr=q/porosity, vf=vs+vr.\n"
          "# Displacement at t; velocities at t-dt/2 (zero initial row). SI, X north, Y east, Z down.\n"
          "# t_s wx_m wy_m wz_m qx_m_s qy_m_s qz_m_s vrx_m_s vry_m_s vrz_m_s vfx_m_s vfy_m_s vfz_m_s\n");
      g_biotStationFiles.push_back(file);
    }
    g_biotStationValues.resize(12*g_biotStationFiles.size());
    if (!nodes.empty()) {
      acoustic_constraint_upload(&g_biotStationNodesDevice, nodes);
      acoustic_constraint_upload(&g_biotStationCoordsDevice, coords);
      solver_pore_pressure_device_check(cudaMalloc((void**)&g_biotStationValuesDevice,
          g_biotStationValues.size()*sizeof(double)), "Biot station allocation");
    }
    g_biotStationsReady = true;
  }
  const int count = g_biotStationFiles.size();
  if (!count) return;
  psolve_launch_biot_stations(s->streams[CUDA_STREAM_MAIN], count,
      g_biotStationNodesDevice, g_biotStationCoordsDevice, s->gpuDataDevice,
      Param.theDeltaT, Param.thePoroelasticPorosity, g_biotStationValuesDevice);
  solver_pore_pressure_device_check(cudaGetLastError(), "Biot station interpolation");
  cudaMemcpyAsync(g_biotStationValues.data(), g_biotStationValuesDevice,
      g_biotStationValues.size()*sizeof(double), cudaMemcpyDeviceToHost,
      s->streams[CUDA_STREAM_MAIN]);
  cudaStreamSynchronize(s->streams[CUDA_STREAM_MAIN]);
  for (int i=0; i<count; ++i) {
    FILE* file = g_biotStationFiles[i];
    fprintf(file, "%.12e", step*Param.theDeltaT);
    for (int d=0; d<12; ++d) fprintf(file, " %.12e", g_biotStationValues[12*i+d]);
    fputc('\n', file);
  }
}

static void biot_close_stations(void) {
  for (FILE* file : g_biotStationFiles) {
    if (fflush(file) || ferror(file))
      solver_abort(__FUNCTION_NAME, NULL, "Cannot write full-Biot fluid stations\n");
    hu_fclose(file);
  }
  g_biotStationFiles.clear();
  g_biotStationValues.clear();
  cudaFree(g_biotStationNodesDevice); g_biotStationNodesDevice = NULL;
  cudaFree(g_biotStationCoordsDevice); g_biotStationCoordsDevice = NULL;
  cudaFree(g_biotStationValuesDevice); g_biotStationValuesDevice = NULL;
  g_biotStationsReady = false;
}

#endif
