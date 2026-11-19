#ifndef HERCULES_POROELASTIC_BRICKUP_DRAIN_SETUP_H
#define HERCULES_POROELASTIC_BRICKUP_DRAIN_SETUP_H

/* Strong reduced-u-p drainage. All geometry/ownership work is setup-only;
 * the existing lumped pressure kernel and pressure exchanges apply the mask.
 * BrickUP's volume-scaled Q1 basis differs from full Biot's native VT P1.
 * Fix positive Q1 trace support, including the unscaled station trace. */
static std::vector<unsigned char> g_brickupSurfaceCut, g_brickupWaterCut;

static void brickup_surface_plane(int ei,double zeta) {
  const double normal[3]={0,0,-1};
  if (!(biot_surface_admittance(ei,normal)>0)) return;
  const elem_t& e=Global.myMesh->elemTable[ei];
  const double h=((edata_t*)e.data)->edgesize;
  g_biotSurfaceArea+=h*h;
  // The Q1 plane integral is exact. For a subcell plane, both node planes
  // have positive support: this conservative restriction fixes the brick.
  for (int j=0;j<8;++j)
    g_biotSurfaceRate[e.lnid[j]]+=h*h*(1+((j&4)?1:-1)*zeta)/8;
}

static void brickup_surface_visit_face(int ei,int axis,int sign,double a,
    double b,double width,int level,bool assemble) {
  const elem_t& e=Global.myMesh->elemTable[ei];
  const edata_t& m=*(edata_t*)e.data;
  const node_t& n0=Global.myMesh->nodeTable[e.lnid[0]];
  const double tick=Global.myMesh->ticksize;
  vector3D_t probe={{n0.x*tick,n0.y*tick,n0.z*tick}};
  probe.x[axis]+=(sign>0?m.edgesize:0)+sign*fmax(tick*.25,width*1e-4);
  probe.x[(axis+1)%3]+=a+.5*width;probe.x[(axis+2)%3]+=b+.5*width;
  int local=0;int32_t neighbor=-1;int8_t neighbor_level=0;
  const bool found=mesh_find_octant_at_domain_point(
      &probe,NULL,NULL,&neighbor_level,&local,&neighbor);
  if (found && neighbor_level>level) {
    for (int i=0;i<2;++i) for (int j=0;j<2;++j)
      brickup_surface_visit_face(ei,axis,sign,a+i*width/2,b+j*width/2,
          width/2,level+1,assemble);
    return;
  }
  if (!assemble) return;
  bool exists=found,porous=false,acoustic=false,dry=false;
  if (found && local && neighbor>=0) {
    const edata_t* other=(edata_t*)Global.myMesh->elemTable[neighbor].data;
    exists=other->Vp>0 && other->rho>0;
    porous=g_poroelasticElementFlags[neighbor]!=0;
    acoustic=element_is_acoustic(other);
    dry=exists && !porous && brickup_surface_dry_cap(neighbor);
  } else if (found) {
    const point_t p={(tick_t)(probe.x[0]/tick),(tick_t)(probe.x[1]/tick),
                     (tick_t)(probe.x[2]/tick)};
    const acoustic_neighbor_t* other=mesh_cached_acoustic_neighbor(p);
    if (!other) solver_abort(__FUNCTION_NAME,NULL,"Missing BrickUP surface neighbor\n");
    exists=other->geid>=0 && other->vp>0 && other->rho>0;
    porous=exists && other->poroelastic_active>0;
    acoustic=exists && material_values_are_acoustic(other->vp,other->vs,other->rho);
    dry=exists && !porous && other->poroelastic_above_groundwater;
  }
  // Internal acoustic contacts retain their existing open/no-flow law.
  if (porous || acoustic) return;
  const bool inside=probe.x[0]>=0 && probe.x[0]<Param.theDomainX &&
      probe.x[1]>=0 && probe.x[1]<Param.theDomainY &&
      probe.x[2]>=0 && probe.x[2]<Param.theDomainZ;
  const bool exposed=!exists && ((axis==2 && sign<0) ||
      (Param.includeTopography==YES && inside));
  // The cut contour/groundwater plane has already supplied the trace.
  // Otherwise a center-classified dry cap drains at its discrete mesh face.
  if ((dry && !g_brickupWaterCut[ei] && !g_brickupSurfaceCut[ei]) ||
      (exposed && !g_brickupSurfaceCut[ei] && !g_brickupWaterCut[ei]))
    biot_surface_add_face(ei,axis,sign,a,b,width,dry);
}

static void brickup_surface_visit(bool assemble) {
  for (int ei=0;ei<Global.myMesh->lenum;++ei) if (g_poroelasticElementFlags[ei]) {
    const elem_t& e=Global.myMesh->elemTable[ei];
    const double h=((edata_t*)e.data)->edgesize;
    for (int d=0;d<3;++d) for (int s=-1;s<=1;s+=2)
      brickup_surface_visit_face(ei,d,s,0,0,h,e.level,assemble);
  }
}

static void brickup_surface_initialize() {
  if (!brickup_surface_strong()) return;
  mesh_t* mesh=Global.myMesh;
  g_biotSurfaceRate.assign(mesh->nharbored,0);
  g_brickupSurfaceCut.assign(mesh->lenum,0);
  g_brickupWaterCut.assign(mesh->lenum,0);
  const bool vt=Param.includeTopography==YES && get_topo_meth()==VT;
  for (int ei=0;ei<mesh->lenum;++ei) if (g_poroelasticElementFlags[ei]) {
    const elem_t& e=mesh->elemTable[ei];
    const double h=((edata_t*)e.data)->edgesize;
    double ground,water;
    poroelastic_element_candidate(mesh,ei,0,&ground,&water);
    // Keep BrickUP's fixed center-sampled groundwater representation. A zero
    // table uses the terrain trace; buried planes can cut the upper half-cell.
    const bool buried=Param.surfaceGroundwaterTable==YES && ground-water>1e-10*h;
    if (buried && water<=.5*h+fmax(1e-6,.01*h)) {
      g_brickupWaterCut[ei]=1;
      brickup_surface_plane(ei,-2*fmin(fmax(water,0.0),.5*h)/h);
    }
    if (vt && !buried && !cpml_pml_flattened_element(ei)) {
      const node_t& n=mesh->nodeTable[e.lnid[0]];
      biot_cut_t cut={};
      if (topo_cut_geometry(n.x*mesh->ticksize,n.y*mesh->ticksize,
                           n.z*mesh->ticksize,h,cut.fraction,&cut.partition)) {
        g_brickupSurfaceCut[ei]=1;
        biot_surface_add_cut(ei,cut);
      }
    }
  }
  g_acousticNeighbors.clear();
  for (;;) {
    g_acousticCollectNeighbors=1;brickup_surface_visit(false);
    g_acousticCollectNeighbors=0;
    int pending=!g_acousticNeighborRequests.empty(),any=0;
    MPI_Allreduce(&pending,&any,1,MPI_INT,MPI_MAX,comm_solver);
    if (!any) break;
    mesh_exchange_acoustic_neighbors();
  }
  brickup_surface_visit(true);
  biot_surface_initialize();
  std::vector<unsigned char>().swap(g_brickupSurfaceCut);
  std::vector<unsigned char>().swap(g_brickupWaterCut);
}

static void brickup_surface_validate_restart() {
  if (!brickup_surface_strong()) return;
  // The checkpoint extension restores both pressure levels to host/device.
  // Reject a changed incompatible boundary instead of resetting saved fluid.
  int invalid=0,any=0;
  for (size_t n=0;n<g_biotSurfaceFixed.size();++n) if (g_biotSurfaceFixed[n])
    if (Global.mySolver->pore_pressure_tm1[n]!=0 ||
        Global.mySolver->pore_pressure_tm2[n]!=0) invalid=1;
  MPI_Allreduce(&invalid,&any,1,MPI_INT,MPI_MAX,comm_solver);
  if (any) solver_abort(__FUNCTION_NAME,NULL,
      "BrickUP strong surface is incompatible with checkpoint pressure; "
      "restart with the original hydraulic boundary and geometry.\n");
}

#endif
