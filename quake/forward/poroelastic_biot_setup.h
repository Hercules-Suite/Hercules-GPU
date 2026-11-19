#ifndef HERCULES_BIOT_SETUP_H
#define HERCULES_BIOT_SETUP_H

/* Private psolve.cc initialization. No work here enters the production step. */
static std::vector<fvector_t> g_biotHostScratch;
static std::vector<solver_float> g_biotScalarScratch;
static double g_biotBoundaryAreas[2][3]; /* sealed/drained, by face-normal axis */
static double g_biotGroundwaterAreas[2][3]; /* buried dry-cap subset */
static std::vector<int32_t> g_biotCutByElement;
static std::vector<poro_cpml_t> g_poroCpml;
static std::vector<int32_t> g_poroCpmlByElement;

static void poro_cpml_register(int ei,const cpml_element_t& stretch,
    const double d[3],const double alpha[3],const double kappa[3]) {
  if (!poroelastic_enabled() ||
      !poroelastic_element_candidate(Global.myMesh,ei,0,NULL,NULL)) return;
  if (g_poroCpmlByElement.empty()) g_poroCpmlByElement.assign(Global.myMesh->lenum,-1);
  const edata_t& mat=*(edata_t*)Global.myMesh->elemTable[ei].data;
  poro_cpml_t c={}; c.stretch=stretch; c.h=mat.edgesize;
  c.rho=mat.rho;
  c.consistent=full_biot_enabled() && Param.thePoroelasticMassMatrix;
  c.volume=c.h*c.h*c.h/8; c.dt=Param.theDeltaT;
  c.storage=poroelastic_biot_storage_modulus_for_element(&mat,ei);
  c.rho_f=Param.thePoroelasticFluidDensity;
  c.fluid_mass=Param.thePoroelasticTortuosity*c.rho_f/Param.thePoroelasticPorosity;
  c.inertia=cpml_inertia(d,alpha,kappa,3,0,c.dt);
  c.darcy=poro_cpml_darcy(d,alpha,kappa,c.dt);
  const double k[3]={Param.thePoroelasticPermeabilityX,
      Param.thePoroelasticPermeabilityY,Param.thePoroelasticPermeabilityZ};
  for (int j=0;j<3;++j) c.drag[j]=k[j]>0?Param.thePoroelasticFluidViscosity/k[j]:0;
  g_poroCpmlByElement[ei]=g_poroCpml.size(); g_poroCpml.push_back(c);
}

static void poro_cpml_initialize(void) {
  if (!poroelastic_enabled()) return;
  long long local=g_poroCpml.size(),total=0;
  MPI_Allreduce(&local,&total,1,MPI_LONG_LONG,MPI_SUM,comm_solver);
  if (!total) return;
  if (Param.includeTopography==YES && get_topo_meth()==VT) {
    const int32_t* elements=topo_get_elements_mapping();
    const topoconstants_t* constants=topo_get_constants();
    for (int t=0;t<topo_get_elements_count();++t) {
      const int ei=elements[t];
      if (g_poroCpmlByElement.empty() || g_poroCpmlByElement[ei]<0) continue;
      for (int q=0;q<5;++q) if (fabs(constants[t].tetraVol[q]-1)>1e-12)
        solver_abort(__FUNCTION_NAME,NULL,
            "Poroelastic CPML requires a regular Cartesian buffer: VT cut cell %d "
            "intersects the layer. Keep terrain cuts inside the physical domain "
            "and extend a flat, mesh-aligned surface through the CPML.\n",ei);
    }
  }
  if (Param.theCpmlNativeStiffness!=1)
    solver_abort(__FUNCTION_NAME,NULL,"Poroelastic CPML requires cpml_stiffness_mode=native.\n");
  if (full_biot_enabled() && acoustic_poroelastic_open_pore_enabled() &&
      Param.theAcousticInterfaceMethod!=ACOUSTIC_INTERFACE_FIELDS)
    solver_abort(__FUNCTION_NAME,NULL,
        "Full-Biot open-pore coupling with porous CPML requires "
        "acoustic_interface_method=field-specific-constraints. Shared geometric "
        "constraints showed late-time growth; local-matching is not qualified "
        "for this combination. Field-specific constraints support flat, FEM "
        "and compatible native VT traces.\n");
  if (!full_biot_enabled() && (get_geostatic_total_time()>0 ||
      Param.theUseCheckPoint || Param.theCheckPointingRate))
    solver_abort(__FUNCTION_NAME,NULL,
        "Poroelastic CPML requires zero geostatic ramp and no checkpoint/restart.\n");
  std::vector<int32_t> map(g_poroelasticLocalElementCount,-1);
  for (int p=0;p<g_poroelasticLocalElementCount;++p)
    if (!g_poroCpmlByElement.empty()) map[p]=g_poroCpmlByElement[g_poroelasticElemMapper[p]];
  gpu_data_t& data=Global.gpuData;
  data.poroCpmlCount=g_poroCpml.size();
  acoustic_constraint_upload(&data.poroCpmlIndexDevice,map);
  acoustic_constraint_upload(&data.poroCpmlCoeffDevice,g_poroCpml);
  const size_t channels=full_biot_enabled()?PORO_CPML_BIOT_STATES:PORO_CPML_REDUCED_STATES;
  const size_t bytes=channels*g_poroCpml.size()*sizeof(solver_float);
  if (bytes) {
    if (GPU_SETUP_MALLOC(&data.poroCpmlStateDevice,bytes)!=cudaSuccess)
      solver_abort(__FUNCTION_NAME,NULL,"Cannot allocate poroelastic CPML histories\n");
    cudaMemset(data.poroCpmlStateDevice,0,bytes);
  }
  if (Global.myID==0)
    monitor_print("Poroelastic CPML: %lld active porous cells, %zu extra histories/cell "
        "(%zu bytes/cell); Cartesian J storage, cofactor pressure coupling, %s. "
        "Reduced undrained fast path is disabled with CPML.\n",total,channels,
        channels*sizeof(solver_float),full_biot_enabled()?"coupled J inertia and Darcy drag":"stretched Darcy diffusion");
  if (Global.myID==0 && full_biot_enabled() && Param.thePoroelasticMassMatrix)
    monitor_print("Full-Biot CPML mass: consistent Q1 inertia and pressure storage, "
        "including polynomial and memory terms; physical solid mass damping "
        "and Darcy drag retain paired lumped J operators.\n");
  if (Global.myID==0 && !full_biot_enabled() && !brickup_surface_strong())
    monitor_print("BrickUP CPML: fixed groundwater classification; dry cells use "
        "solid CPML, saturated cells use porous CPML. Drained-plane compliance "
        "and storage share the Cartesian J filter; no additional histories.\n");
  if (Global.myID==0 && brickup_surface_strong())
    monitor_print("BrickUP CPML: strong drained pressure uses the shared fixed-node "
        "mask; legacy drain compliance and its J-filter tail are disabled. "
        "Storage and diffusion histories retain their existing equations.\n");
}

static bool biot_vt_enabled(void) {
  return full_biot_enabled() && Param.includeTopography==YES && get_topo_meth()==VT;
}

static void biot_station_weights(int station, double weights[8]) {
  const station_t& s=Param.myStations[station];
  int tet,part;
  if (biot_vt_enabled() && !cpml_pml_flattened_element(s.owningElementIndex) &&
      get_topo_station_tet(station,&tet,&part)) {
    double x[3]; for (int d=0; d<3; ++d) x[d]=(1+s.localcoords.x[d])/2;
    const int containing=biot_tet_weights(part,x,weights);
    if (containing<0)
      solver_abort(__FUNCTION_NAME,NULL,"Full-Biot station outside its VT tetrahedron\n");
    double fractions[5]; int owner_part;
    get_tetraProps(s.owningElementIndex,fractions,&owner_part);
    if (owner_part!=part || !(fractions[containing]>0))
      solver_abort(__FUNCTION_NAME,NULL,
          "Full-Biot station %d lies in a discarded VT tetrahedron; refine the "
          "surface mesh or move the receiver into retained material.\n",s.id);
  } else {
    for (int j=0; j<8; ++j) {
      weights[j]=0.125;
      for (int d=0; d<3; ++d)
        weights[j]*=1+(((j>>d)&1)?1:-1)*s.localcoords.x[d];
    }
  }
}

static bool biot_stabilization_active(void) {
  return full_biot_enabled() && !poroelastic_zero_permeability() &&
         Param.thePoroelasticPressureStabilization > 0;
}

#include "poroelastic_drain_setup.h"

/* Resolve each face down to the adjacent leaf size, including remote leaves.
 * A sealed subface constrains its supporting Q1 trace nodes. */
static void biot_visit_face(int32_t ei, int axis, int sign, double a, double b,
                            double width, int level,
                            std::vector<fvector_t>* locks) {
  const elem_t& e = Global.myMesh->elemTable[ei];
  const edata_t& mat = *(edata_t*)e.data;
  const node_t& n0 = Global.myMesh->nodeTable[e.lnid[0]];
  const double tick = Global.myMesh->ticksize;
  const double low[3] = {n0.x*tick, n0.y*tick, n0.z*tick};
  const int t0 = (axis+1)%3, t1 = (axis+2)%3;
  vector3D_t probe = {{low[0], low[1], low[2]}};
  probe.x[axis] += (sign > 0 ? mat.edgesize : 0) +
                   sign * fmax(tick*0.25, width*1e-4);
  probe.x[t0] += a + 0.5*width;
  probe.x[t1] += b + 0.5*width;
  int local = 0;
  int32_t neighbor = -1;
  int8_t neighbor_level = 0;
  const int found = mesh_find_octant_at_domain_point(
      &probe, NULL, NULL, &neighbor_level, &local, &neighbor);
  if (found && neighbor_level > level) {
    for (int i=0; i<2; ++i) for (int j=0; j<2; ++j)
      biot_visit_face(ei, axis, sign, a+i*width/2, b+j*width/2,
                      width/2, level+1, locks);
    return;
  }
  if (!locks) return;
  bool porous = false, acoustic = false, exists = found, above_groundwater = false;
  if (found && local && neighbor >= 0) {
    const edata_t* other = (edata_t*)Global.myMesh->elemTable[neighbor].data;
    exists = other->Vp > 0 && other->rho > 0;
    porous = g_poroelasticElementFlags[neighbor] != 0;
    acoustic = element_is_acoustic(other);
    above_groundwater = !porous && Param.surfaceGroundwaterTable == YES &&
        poroelastic_element_map_status(Global.myMesh,neighbor,NULL,NULL) ==
            POROELASTIC_MAP_ABOVE_GROUNDWATER;
  } else if (found) {
    const point_t p = {(tick_t)(probe.x[0]/tick), (tick_t)(probe.x[1]/tick),
                       (tick_t)(probe.x[2]/tick)};
    const acoustic_neighbor_t* other = mesh_cached_acoustic_neighbor(p);
    if (!other)
      solver_abort(__FUNCTION_NAME, NULL, "Missing full-Biot neighbor after exchange\n");
    exists = other->geid >= 0 && other->vp > 0 && other->rho > 0;
    porous = exists && other->poroelastic_active > 0;
    acoustic = exists && material_values_are_acoustic(other->vp, other->vs, other->rho);
    above_groundwater = exists && other->poroelastic_above_groundwater;
  }
  if (porous) return;
  if (acoustic && acoustic_poroelastic_open_pore_enabled()) return;
  /* FEM terrain consists of whole Hex8 cells. Every exposed stair-step face
   * (including vertical risers) uses the selected hydraulic surface condition.
   * The domain's lateral/bottom truncations remain sealed. Resolve air using
   * the owner's material, just as for local air and carved voids. */
  const bool inside = probe.x[0] >= 0 && probe.x[0] < Param.theDomainX &&
                      probe.x[1] >= 0 && probe.x[1] < Param.theDomainY &&
                      probe.x[2] >= 0 && probe.x[2] < Param.theDomainZ;
  const bool exposed = !exists && ((axis == 2 && sign < 0) ||
      (Param.includeTopography == YES && inside));
  // A dry cap caused by the water table is a hydraulic surface, not the
  // impermeable material/depth-limit interface. Solid motion stays shared.
  const bool groundwater = exists && !acoustic && above_groundwater;
  if ((exposed || groundwater) &&
      Param.thePoroelasticSurfacePressureBC == POROELASTIC_SURFACE_PRESSURE_DRAINED) {
    g_biotBoundaryAreas[1][axis] += width*width;
    if (groundwater) g_biotGroundwaterAreas[1][axis] += width*width;
    if (groundwater && biot_vt_enabled() && g_biotCutByElement[ei]>=0) {
      biot_cut_t cut={};get_tetraProps(ei,cut.fraction,&cut.partition);
      biot_surface_add_vt_groundwater_face(ei,cut,axis,sign,a,b,width);
    } else if (!biot_vt_enabled() || g_biotCutByElement[ei]<0)
      biot_surface_add_face(ei,axis,sign,a,b,width,groundwater);
    return;
  }
  /* A VT cut cell exposes an internal, generally inclined surface, not the
   * fictitious cube faces. Its normal constraint is assembled separately. */
  if (exposed && biot_vt_enabled() && g_biotCutByElement[ei]>=0) return;
  g_biotBoundaryAreas[0][axis] += width*width;
  if (groundwater) g_biotGroundwaterAreas[0][axis] += width*width;
  for (int j=0; j<8; ++j)
    if (((j >> axis) & 1) == (sign > 0 ? 1 : 0))
      (*locks)[e.lnid[j]].f[axis] += 1;
}

static void biot_visit_boundaries(std::vector<fvector_t>* locks) {
  for (int32_t ei=0; ei<Global.myMesh->lenum; ++ei) {
    if (!g_poroelasticElementFlags[ei]) continue;
    const elem_t& e = Global.myMesh->elemTable[ei];
    const double h = ((edata_t*)e.data)->edgesize;
    for (int axis=0; axis<3; ++axis) for (int sign=-1; sign<=1; sign+=2)
      biot_visit_face(ei, axis, sign, 0, 0, h, e.level, locks);
  }
}

static void biot_assemble_nodes(std::vector<fvector_t>& value) {
  mysolver_t* s = Global.mySolver;
  schedule_senddata(s->dn_sched, value.data(), 3, CONTRIBUTION, 145);
  compute_adjust(value.data(), 3, DISTRIBUTION, ACOUSTIC_CONSTRAINT_SOLID);
  schedule_senddata(s->an_sched, value.data(), 3, CONTRIBUTION, 146);
}

static void biot_vt_interface_weights(void) {
  if (!biot_vt_enabled()) return;
  const int count=Global.acousticInterfaceFaces;
  std::vector<double> weights((size_t)count*32);
  for (int f=0; f<count; ++f) {
    if (!(Global.acousticInterfaceFaceAreas[f]>0)) continue;
    if (g_acousticSolidTraceNodes.size()<(size_t)(8*f+8))
      solver_abort(__FUNCTION_NAME,NULL,"Missing full-Biot VT solid interface trace\n");
    /* Remote trace storage omits off-face vertices (their IDs are dummies).
     * Reconstruct the owner-confirmed cube bounds from the face coordinates
     * and halfwidth ratio, never from opposite trace-node IDs. */
    const int water=Global.acousticInterfaceElementIds[f];
    const elem_t& element=Global.myMesh->elemTable[water];
    const node_t& n0=Global.myMesh->nodeTable[element.lnid[0]];
    const double tick=Global.myMesh->ticksize;
    double fractions[5]; int part,axis,sign,tangents[2];
    mesh_acoustic_face_axis(Global.acousticInterfaceFaceIds[f],&axis,&sign,tangents);
    const double wh=((edata_t*)element.data)->edgesize;
    const int tangent=tangents[0];
    const double h=wh*Global.acousticInterfaceFaceHalfWidths[3*f+tangent]/
                      Global.acousticInterfaceAdjacentHalfWidths[3*f+tangent];
    if (!(h>0) || !isfinite(h))
      solver_abort(__FUNCTION_NAME,NULL,"Invalid full-Biot VT solid trace size\n");
    double low[3]={n0.x*tick,n0.y*tick,n0.z*tick};
    for (int d=0; d<3; ++d)
      low[d]+=0.5*wh*(Global.acousticInterfaceFaceCenters[3*f+d]+1)-
              0.5*h*(Global.acousticInterfaceAdjacentCenters[3*f+d]+1);
    const bool pml=Global.acousticInterfacePmlSideMask &&
        (Global.acousticInterfacePmlSideMask[f]&ACOUSTIC_INTERFACE_PML_ELASTIC_SIDE);
    const bool cut=!pml && topo_cut_geometry(low[0],low[1],low[2],h,fractions,&part)!=0;
    for (int q=0; q<4; ++q) {
      double x[3],*N=weights.data()+32*f+8*q;
      for (int d=0; d<3; ++d) x[d]=Global.acousticInterfaceAdjacentCenters[3*f+d];
      for (int a=0; a<2; ++a)
        x[tangents[a]]+=Global.acousticInterfaceAdjacentHalfWidths[3*f+tangents[a]]*
                        (((q>>(1-a))&1)?1:-1)*0.57735026918962576451;
      if (cut) {
        for (int d=0; d<3; ++d) x[d]=(1+x[d])/2;
        const int t=biot_tet_weights(part,x,N);
        if (t<0) solver_abort(__FUNCTION_NAME,NULL,"Invalid VT interface quadrature point\n");
        /* Native VT is a volume-fraction approximation. Integrate only the
         * retained solid trace, using identical weights in both directions. */
        if (!(fractions[t]>0) || low[2]+h*x[2]<point_elevation(low[0]+h*x[0],low[1]+h*x[1]))
          std::fill(N,N+8,0.0);
      } else {
        for (int j=0; j<8; ++j) {
          N[j]=0.125;
          for (int d=0; d<3; ++d) N[j]*=1+(((j>>d)&1)?1:-1)*x[d];
        }
      }
    }
  }
  acoustic_constraint_upload(&Global.gpuData.biotInterfaceWeightsDevice,weights);
}

static void biot_initialize(void) {
  if (!full_biot_enabled()) return;
  if (g_poroelasticElementCount == 0)
    solver_abort(__FUNCTION_NAME, NULL, "full_biot selected but no active porous elements\n");
  if (biot_vt_enabled() && acoustic_pressure_enabled() &&
      Param.acousticInterfaceCoupling==YES &&
      Param.theAcousticInterfaceProjectionMode!=ACOUSTIC_INTERFACE_PROJECTION_FACE_MASS)
    solver_abort(__FUNCTION_NAME,NULL,
        "Full-Biot VT acoustic coupling requires acoustic_interface_projection=face_mass "
        "for reciprocal tetrahedral solid/fluid traces.\n");
  /* Both Hex8 projection kernels exchange the same pressure/relative-motion
   * weights in opposite directions. VT still needs its four-point tetrahedral
   * traces, enforced above. */
  if (acoustic_poroelastic_open_pore_enabled() &&
      Param.acousticInterfaceCoupling != YES)
    solver_abort(__FUNCTION_NAME, NULL,
        "Full-Biot open pores require acoustic_interface_coupling=yes "
        "for fluid traction/flux exchange.\n");
  if (Param.theUseCheckPoint || Param.theCheckPointingRate)
    solver_abort(__FUNCTION_NAME, NULL,
        "Full-Biot checkpoint/restart is unsupported by the legacy checkpoint "
        "format; use use_checkpoint=0 and checkpointing_rate=0.\n");
  if (get_geostatic_total_time() > 0 || get_fixedbase_flag() == YES ||
      Param.drmImplement == YES)
    solver_abort(__FUNCTION_NAME, NULL,
        "Full Biot does not yet support geostatic, fixed-base or DRM "
        "histories; these do not provide the independent fluid state.\n");
  const size_t count = Global.myMesh->nharbored;
  if (biot_surface_enforced()) g_biotSurfaceRate.assign(count,0);
  const size_t bytes = count * sizeof(fvector_t);
  std::vector<fvector_t> mass(count), drag(count), locks(count), surface;
  std::vector<biot_cut_t> cuts;
  std::vector<int32_t> cut_index;
  if (biot_vt_enabled()) {
    surface.resize(count);
    g_biotCutByElement.assign(Global.myMesh->lenum,-1);
    cut_index.assign(g_poroelasticLocalElementCount,-1);
    const topoconstants_t* constants=topo_get_constants();
    const int32_t* elements=topo_get_elements_mapping();
    for (int t=0; t<topo_get_elements_count(); ++t) {
      const int ei=elements[t];
      if (!g_poroelasticElementFlags[ei] || cpml_pml_flattened_element(ei)) continue;
      biot_cut_t cut={}; cut.partition=constants[t].cube_part;
      for (int j=0; j<5; ++j) cut.fraction[j]=constants[t].tetraVol[j];
      g_biotCutByElement[ei]=cuts.size(); cuts.push_back(cut);
    }
    for (int p=0; p<g_poroelasticLocalElementCount; ++p)
      cut_index[p]=g_biotCutByElement[g_poroelasticElemMapper[p]];
  }
  g_biotHostScratch.resize(count);
  const double rf = Param.thePoroelasticFluidDensity;
  const double inertia = Param.thePoroelasticTortuosity * rf / Param.thePoroelasticPorosity;
  const double k[3] = {Param.thePoroelasticPermeabilityX,
                       Param.thePoroelasticPermeabilityY, Param.thePoroelasticPermeabilityZ};
  const double viscosity = Param.thePoroelasticFluidViscosity;
  double vt_wave_dt=HUGE_VAL;
  if (!isfinite(rf) || !isfinite(inertia) || !isfinite(viscosity) ||
      !isfinite(Param.thePoroelasticBiotAlpha) ||
      !isfinite(Param.thePoroelasticFluidBulkModulus))
    solver_abort(__FUNCTION_NAME, NULL, "Full-Biot material values must be finite\n");
  for (int d=0; d<3; ++d)
    if (!isfinite(k[d]))
      solver_abort(__FUNCTION_NAME, NULL, "Full-Biot permeability must be finite\n");
  for (int32_t ei=0; ei<Global.myMesh->lenum; ++ei) {
    if (!g_poroelasticElementFlags[ei]) continue;
    const elem_t& e = Global.myMesh->elemTable[ei];
    const edata_t& material = *(edata_t*)e.data;
    if (!(material.rho*inertia > rf*rf))
      solver_abort(__FUNCTION_NAME, NULL, "Nonpositive full-Biot material mass matrix\n");
    if (biot_vt_enabled()) {
      double fast2,slow2;
      biot_compressional_speeds_squared(material.rho,rf,inertia,
          material.rho*material.Vp*material.Vp,Param.thePoroelasticBiotAlpha,
          poroelastic_biot_storage_modulus_for_element(&material,ei),&fast2,&slow2);
      vt_wave_dt=fmin(vt_wave_dt,0.1*material.edgesize/sqrt(fast2));
    }
    double volumes[8];
    if (biot_vt_enabled() && g_biotCutByElement[ei]>=0) {
      /* Match native VT solid mass exactly, including discarded tetrahedra. */
      std::fill(volumes,volumes+8,0.0);
      const biot_cut_t& cut=cuts[g_biotCutByElement[ei]];
      biot_surface_add_cut(ei,cut);
      for (int t=0; t<5; ++t) {
        int ids[4]; double grad[4][3];
        const double v=biot_tet_geometry(cut.partition,t,ids,grad)*
                       pow(material.edgesize,3)*cut.fraction[t]/4;
        for (int a=0; a<4; ++a) volumes[ids[a]]+=v;
      }
    } else std::fill(volumes,volumes+8,pow(material.edgesize,3)/8);
    for (int j=0; j<8; ++j) {
      const int32_t n=e.lnid[j];
      const int pml=g_poroCpmlByElement.empty()?-1:g_poroCpmlByElement[ei];
      const double volume=volumes[j]*(pml>=0?g_poroCpml[pml].inertia.mass:1);
      if (!(volume>0)) continue;
      if (biot_vt_enabled() && g_biotCutByElement[ei]>=0 &&
          Param.thePoroelasticSurfacePressureBC==POROELASTIC_SURFACE_PRESSURE_UNDRAINED)
        surface[n].f[0]+=volume;
      mass[n].f[0] += rf*volume;
      mass[n].f[1] += inertia*volume;
      for (int d=0; d<3; ++d) {
        if (k[d] == 0) locks[n].f[d] += 1;
        else drag[n].f[d] += viscosity/k[d]*volume;
      }
    }
  }
  if (biot_vt_enabled()) {
    double bound;
    MPI_Allreduce(&vt_wave_dt,&bound,1,MPI_DOUBLE,MPI_MIN,comm_solver);
    warn_timestep_estimate("full-Biot VT wave",bound);
    if (Global.myID==0) monitor_print("Full-Biot VT conservative wave dt <= %.9g s (advisory only).\n",bound);
  }
  /* Previous acoustic discovery can precede final CPML/porous activity. */
  g_acousticNeighbors.clear();
  for (;;) {
    g_acousticCollectNeighbors = 1;
    biot_visit_boundaries(NULL);
    g_acousticCollectNeighbors = 0;
    int pending = !g_acousticNeighborRequests.empty(), any = 0;
    MPI_Allreduce(&pending, &any, 1, MPI_INT, MPI_MAX, comm_solver);
    if (!any) break;
    mesh_exchange_acoustic_neighbors();
  }
  memset(g_biotBoundaryAreas, 0, sizeof(g_biotBoundaryAreas));
  memset(g_biotGroundwaterAreas, 0, sizeof(g_biotGroundwaterAreas));
  biot_visit_boundaries(&locks);
  biot_surface_initialize();
  double boundary_areas[2][3];
  MPI_Reduce(g_biotBoundaryAreas, boundary_areas, 6, MPI_DOUBLE, MPI_SUM,
             0, comm_solver);
  if (Global.myID == 0)
    monitor_print("Full-Biot hydraulic face areas (m^2), normal axes X/Y/Z: "
        "sealed [%.9g %.9g %.9g], drained [%.9g %.9g %.9g].\n",
        boundary_areas[0][0], boundary_areas[0][1], boundary_areas[0][2],
        boundary_areas[1][0], boundary_areas[1][1], boundary_areas[1][2]);
  if (Param.surfaceGroundwaterTable == YES) {
    double groundwater_areas[2][3];
    MPI_Reduce(g_biotGroundwaterAreas,groundwater_areas,6,MPI_DOUBLE,MPI_SUM,0,comm_solver);
    if (Global.myID==0)
      monitor_print("Full-Biot buried groundwater: fixed ground-relative depth grid; "
          "whole-element centroid saturation, mesh-face hydraulic boundary. "
          "Dry-cap face areas (m^2), X/Y/Z: sealed [%.9g %.9g %.9g], drained [%.9g %.9g %.9g]. "
          "Solid motion remains shared; material properties are not rescaled. "
          "Refine near the water table; subcell saturation and moving/unsaturated flow are not modeled.\n",
          groundwater_areas[0][0],groundwater_areas[0][1],groundwater_areas[0][2],
          groundwater_areas[1][0],groundwater_areas[1][1],groundwater_areas[1][2]);
  }
  biot_assemble_nodes(mass);
  biot_assemble_nodes(drag);
  biot_assemble_nodes(locks);
  if (biot_vt_enabled()) biot_assemble_nodes(surface);
  std::vector<biot_boundary_t> boundaries;
  for (size_t n=0; n<count; ++n) {
    if (Global.myMesh->nodeTable[n].ismine && mass[n].f[1] > 0) {
      const double ms = Global.mySolver->nTable[n].mass_simple;
      if (ms > 0 && !(ms*mass[n].f[1] > mass[n].f[0]*mass[n].f[0]))
        solver_abort(__FUNCTION_NAME, NULL, "Nonpositive assembled full-Biot mass matrix\n");
    }
    for (int d=0; d<3; ++d) if (locks[n].f[d] > 0) drag[n].f[d] = -1;
    if (biot_vt_enabled() && surface[n].f[0]>0 &&
        Global.myMesh->nodeTable[n].ismine && mass[n].f[1]>0) {
      double normal[3], D[3];
      const node_t& node=Global.myMesh->nodeTable[n];
      topo_surface_normal(node.x*Global.myMesh->ticksize,
                          node.y*Global.myMesh->ticksize,normal);
      for (int d=0; d<3; ++d) D[d]=drag[n].f[d];
      biot_boundary_t bc={}; bc.node=n;
      biot_boundary_modes(normal,D,&bc); boundaries.push_back(bc);
    }
  }
  gpu_data_t& data = Global.gpuData;
  acoustic_constraint_upload(&data.biotCutIndexDevice,cut_index);
  acoustic_constraint_upload(&data.biotCutsDevice,cuts);
  acoustic_constraint_upload(&data.biotBoundaryDevice,boundaries);
  data.biotBoundaryCount=boundaries.size();
  biot_vt_interface_weights();
  fvector_t** arrays[] = {&data.biotW1Device, &data.biotW2Device,
      &data.biotW3Device, &data.biotForceDevice, &data.biotMassDevice, &data.biotDragDevice};
  for (fvector_t** ptr : arrays) {
    if (GPU_SETUP_MALLOC(ptr, bytes) != cudaSuccess)
      solver_abort(__FUNCTION_NAME, NULL, "Cannot allocate full-Biot nodal state\n");
    cudaMemset(*ptr, 0, bytes);
  }
  cudaMemcpy(data.biotMassDevice, mass.data(), bytes, cudaMemcpyHostToDevice);
  cudaMemcpy(data.biotDragDevice, drag.data(), bytes, cudaMemcpyHostToDevice);
  if (biot_stabilization_active()) {
    std::vector<fvector_t> coefficients(g_poroelasticLocalElementCount);
    double local_dt = HUGE_VAL, global_dt = HUGE_VAL;
    for (int pi=0; pi<g_poroelasticLocalElementCount; ++pi) {
      const int ei = g_poroelasticElemMapper[pi];
      const edata_t& mat = *(edata_t*)Global.myMesh->elemTable[ei].data;
      const double M = poroelastic_biot_storage_modulus_for_element(&mat, ei);
      double fast2, slow2;
      biot_compressional_speeds_squared(mat.rho, rf, inertia,
          mat.rho*mat.Vp*mat.Vp, Param.thePoroelasticBiotAlpha, M, &fast2, &slow2);
      const double h=mat.edgesize, c=sqrt(slow2);
      double max_nu=0;
      for (int d=0; d<3; ++d) {
        const double nu = k[d] > 0 ? Param.thePoroelasticPressureStabilization*c*h*h*h /
            (1+viscosity*h/(k[d]*inertia*c)) : 0;
        coefficients[pi].f[d] = sqrt(nu);
        max_nu=fmax(max_nu,nu);
      }
      /* Conservative wave+dissipation bound, keeping each contribution below
       * unity in the scalar Fourier stability condition. */
      local_dt=fmin(local_dt,(biot_vt_enabled()?0.1:0.25)*h/sqrt(fast2));
      /* A native P1 tetrahedron has ||C^-1 K||_inf <=24 sqrt(nu)/h^2;
       * the Hex8 bound is 12. Positive volume fractions cancel against the
       * matching nodal capacity. Assembly/positive hanging weights preserve it. */
      if (max_nu>0) local_dt=fmin(local_dt,h*h*h*h/((biot_vt_enabled()?576:144)*max_nu));
    }
    MPI_Allreduce(&local_dt,&global_dt,1,MPI_DOUBLE,MPI_MIN,comm_solver);
    warn_timestep_estimate("full-Biot wave/stabilization",global_dt);
    if (!coefficients.empty()) {
      if (GPU_SETUP_MALLOC(&data.biotStabilizationCoeffDevice,
                          coefficients.size()*sizeof(fvector_t)) != cudaSuccess)
        solver_abort(__FUNCTION_NAME,NULL,"Cannot allocate Biot stabilization coefficients\n");
      cudaMemcpy(data.biotStabilizationCoeffDevice,coefficients.data(),
                 coefficients.size()*sizeof(fvector_t),cudaMemcpyHostToDevice);
    }
    if (GPU_SETUP_MALLOC(&data.biotStabilizationWorkDevice,
                        count*sizeof(solver_float)) != cudaSuccess)
      solver_abort(__FUNCTION_NAME,NULL,"Cannot allocate Biot stabilization work\n");
    g_biotScalarScratch.resize(count);
    if (Global.myID==0)
      monitor_print("Full-Biot pressure stabilization: gamma=%.6g; fourth-order "
          "numerical dissipation K C^-1 K, O(h^3) in the inviscid limit; "
          "drag-scaled, zero in impermeable directions. Conservative dt <= %.9g s (advisory only).\n",
          Param.thePoroelasticPressureStabilization,global_dt);
  }
  if (Global.myID == 0)
    monitor_print("Full Biot u-w-p: saturated mixture density; tortuosity %.6g; "
        "w=porosity*(fluid displacement-solid displacement); %s Darcy "
        "integration, wave CFL required. Lateral/bottom fluid truncations are sealed; "
        "mechanical dashpots are approximate; abc_type=cpml uses the mixed poroelastic layer.\n",
        Param.thePoroelasticTortuosity,Param.thePoroelasticMassMatrix?"SDIRK2":"exact local");
  if (Global.myID == 0 && Param.includeTopography == YES && get_topo_meth()==FEM &&
      Param.surfaceGroundwaterTable == NO)
    monitor_print("Full-Biot FEM topography: whole Hex8 stair-step cells; "
        "hydraulic surface conditions use each exposed mesh-face normal. "
        "Saturation follows the physical terrain, independently of full/sqd "
        "material sampling; retained surface cells are fully saturated.\n");
  if (biot_vt_enabled()) {
    long long local[2]={(long long)cuts.size(),(long long)boundaries.size()},global[2];
    MPI_Reduce(local,global,2,MPI_LONG_LONG,MPI_SUM,0,comm_solver);
    if (Global.myID==0)
      monitor_print("Full-Biot VT: %lld cut cells, %lld owned surface constraints. "
          "P1 tetrahedral coupling, mass and stabilization use native retained "
          "volume fractions. Undrained terrain uses nodal tangent modes in the "
          "cut-cell support band (O(h) boundary approximation); drained terrain "
          "has natural p=0. %s hanging constraints.\n",
          global[0],global[1],Param.theAcousticInterfaceMethod==ACOUSTIC_INTERFACE_FIELDS?
          "Field-specific native P1/Hex8":"Shared geometric");
  }
  if (Global.myID == 0 && acoustic_poroelastic_open_pore_enabled())
    monitor_print("Full-Biot open pores: water pressure supplies natural fluid "
        "traction; acoustic normal flux is solid velocity plus Darcy flux; "
        "%s interface integration; the reduced-u-p pressure penalty is inactive.\n",
        Param.theAcousticInterfaceProjectionMode == ACOUSTIC_INTERFACE_PROJECTION_FACE_MASS
            ? "four-point face_mass" : "one-point Hex8");
  if (Global.myID == 0 &&
      Param.thePoroelasticSurfacePressureBC == POROELASTIC_SURFACE_PRESSURE_DRAINED)
    monitor_print("Full-Biot drained hydraulic surface: natural excess fluid pressure "
        "p=0, normal relative flow permitted; reduced-u-p drain penalty inactive.\n");
}

#endif
