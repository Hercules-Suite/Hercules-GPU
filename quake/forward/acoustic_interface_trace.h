#ifndef HERCULES_ACOUSTIC_INTERFACE_TRACE_H
#define HERCULES_ACOUSTIC_INTERFACE_TRACE_H

/* Field-specific constraints make the solid and acoustic traces independent.
 * A remote solid face must therefore use its own vertices, not the local
 * acoustic face. Add missing vertices as ordinary non-owned harbored nodes
 * before solver allocation. Existing force-contribution/displacement-sharing
 * schedules then carry the trace on both host and CUDA-aware MPI paths.
 * Every constraint method needs this trace for face_mass: quadrature must not
 * depend on which rank owns the solid. No new per-step exchange or host staging
 * is needed. Open-pore hydraulic quadrature also requires the solid trace for
 * mechanical point projection; a coarse acoustic trace cannot carry the
 * solid-side shape weights, including independent fine-solid pressure nodes. */
static std::vector<int32_t> g_acousticSolidTraceNodes;
static std::vector<int32link_t*> g_acousticTraceShares;

static bool acoustic_solid_trace_enabled() {
  return acoustic_pressure_enabled() && Param.acousticInterfaceCoupling == YES &&
      (Param.theAcousticInterfaceMethod != ACOUSTIC_INTERFACE_ORIGINAL ||
       Param.theAcousticInterfaceProjectionMode == ACOUSTIC_INTERFACE_PROJECTION_FACE_MASS ||
       acoustic_poroelastic_open_pore_enabled());
}

static const int32_t* acoustic_solid_trace_nodes(int64_t face) {
  return g_acousticSolidTraceNodes.empty() ? NULL
                                         : g_acousticSolidTraceNodes.data() + 8 * face;
}

template <class T>
static std::vector<T> acoustic_trace_exchange(
    const std::vector<std::vector<T>>& outgoing) {
  const int ranks = Global.theGroupSize;
  std::vector<int> sc(ranks), rc(ranks), so(ranks), ro(ranks);
  int64_t ns = 0, nr = 0;
  for (int r = 0; r < ranks; ++r) {
    if (outgoing[r].size() > (size_t)INT_MAX)
      solver_abort(__FUNCTION_NAME, NULL, "Solid trace MPI count overflow\n");
    sc[r] = (int)outgoing[r].size();
  }
  MPI_Alltoall(sc.data(), 1, MPI_INT, rc.data(), 1, MPI_INT, comm_solver);
  for (int r = 0; r < ranks; ++r) {
    so[r] = (int)ns; ro[r] = (int)nr;
    ns += sc[r]; nr += rc[r];
    if (ns > INT_MAX || nr > INT_MAX)
      solver_abort(__FUNCTION_NAME, NULL, "Solid trace MPI offset overflow\n");
  }
  std::vector<T> send;
  send.reserve((size_t)ns);
  for (const auto& rank : outgoing) send.insert(send.end(), rank.begin(), rank.end());
  std::vector<T> receive((size_t)nr);
  MPI_Datatype type;
  MPI_Type_contiguous(sizeof(T), MPI_BYTE, &type);
  MPI_Type_commit(&type);
  MPI_Alltoallv(send.data(), sc.data(), so.data(), type,
               receive.data(), rc.data(), ro.data(), type, comm_solver);
  MPI_Type_free(&type);
  return receive;
}

struct acoustic_trace_request_t {
  vector3D_t point;
  int64_t face, geid;
  int rank;
};
struct acoustic_trace_node_t {
  int64_t gnid;
  tick_t x, y, z;
  int owner, anchored;
};
struct acoustic_trace_reply_t {
  int64_t face;
  acoustic_trace_node_t nodes[8];
  double lower[3], h;
};
struct acoustic_trace_share_t { int64_t gnid; int rank; };

static void acoustic_solid_trace_build() {
  if (!acoustic_solid_trace_enabled()) return;
  mesh_t* mesh = Global.myMesh;
  const int ranks = Global.theGroupSize;
  const int64_t faces = Global.acousticInterfaceFaces;
  if (faces > INT_MAX / 8)
    solver_abort(__FUNCTION_NAME, NULL, "Solid trace index overflow\n");
  g_acousticSolidTraceNodes.resize((size_t)faces * 8);
  std::vector<std::vector<acoustic_trace_request_t>> requests(ranks);
  const point_t* intervals = octor_getintervaltable(Global.myOctree);
  for (int64_t f = 0; f < faces; ++f) {
    const int32_t adjacent = Global.acousticInterfaceAdjacentElementIndexes[f];
    if (adjacent >= 0) {
      for (int c = 0; c < 8; ++c)
        g_acousticSolidTraceNodes[8*f+c] = mesh->elemTable[adjacent].lnid[c];
      continue;
    }
    const elem_t& water = mesh->elemTable[Global.acousticInterfaceElementIds[f]];
    const node_t& lower = mesh->nodeTable[water.lnid[0]];
    const double h = ((edata_t*)water.data)->edgesize;
    const double xyz[3] = {lower.x*mesh->ticksize, lower.y*mesh->ticksize,
                           lower.z*mesh->ticksize};
    acoustic_trace_request_t request = {};
    request.face = f; request.rank = Global.myID;
    request.geid = Global.acousticInterfaceAdjacentElementIds[f];
    for (int j = 0; j < 3; ++j)
      request.point.x[j] = xyz[j] + .5*h*(Global.acousticInterfaceFaceCenters[3*f+j]+1);
    // Use the Cartesian face, not a reconstructed bathymetry normal.
    const int bit = Global.acousticInterfaceFaceIds[f];
    const int axis = bit <= ACOUSTIC_FACE_X_PLUS ? 0 :
                     (bit <= ACOUSTIC_FACE_Y_PLUS ? 1 : 2);
    const int sign = (bit & (ACOUSTIC_FACE_X_PLUS | ACOUSTIC_FACE_Y_PLUS |
                             ACOUSTIC_FACE_Z_PLUS)) ? 1 : -1;
    request.point.x[axis] += sign * .25 * mesh->ticksize;
    const point_t ticks = {(tick_t)(request.point.x[0]/mesh->ticksize),
                           (tick_t)(request.point.x[1]/mesh->ticksize),
                           (tick_t)(request.point.x[2]/mesh->ticksize)};
    requests[mesh_handoff_owner(intervals, ranks, &ticks)].push_back(request);
  }
  const auto incoming = acoustic_trace_exchange(requests);
  std::vector<std::vector<acoustic_trace_reply_t>> answers(ranks);
  for (const auto& request : incoming) {
    int local = 0; int32_t e = -1;
    if (!mesh_find_octant_at_domain_point(&request.point, NULL, NULL, NULL, &local, &e) ||
        !local || e < 0 || mesh->elemTable[e].geid != request.geid)
      solver_abort(__FUNCTION_NAME, NULL, "Unresolved solid trace element\n");
    const elem_t& element = mesh->elemTable[e];
    acoustic_trace_reply_t reply = {};
    reply.face = request.face; reply.h = ((edata_t*)element.data)->edgesize;
    for (int c = 0; c < 8; ++c) {
      const node_t& n = mesh->nodeTable[element.lnid[c]];
      reply.nodes[c] = {n.gnid, n.x, n.y, n.z,
                        n.ismine ? Global.myID : n.proc.ownerid, n.isanchored};
    }
    reply.lower[0] = reply.nodes[0].x*mesh->ticksize;
    reply.lower[1] = reply.nodes[0].y*mesh->ticksize;
    reply.lower[2] = reply.nodes[0].z*mesh->ticksize;
    answers[request.rank].push_back(reply);
  }
  const auto replies = acoustic_trace_exchange(answers);
  // Keep the lookup proportional to interface nodes, not the volume mesh.
  std::map<int64_t, int32_t> ids;
  for (const auto& reply : replies)
    for (const auto& n : reply.nodes) ids.emplace(n.gnid, -1);
  for (int32_t n = 0; n < mesh->nharbored; ++n) {
    const auto it = ids.find(mesh->nodeTable[n].gnid);
    if (it != ids.end()) it->second = n;
  }
  std::vector<node_t> added;
  std::vector<std::vector<acoustic_trace_share_t>> shares(ranks);
  for (const auto& reply : replies) {
    const int64_t f = reply.face;
    const elem_t& water = mesh->elemTable[Global.acousticInterfaceElementIds[f]];
    const node_t& low = mesh->nodeTable[water.lnid[0]];
    const double wh = ((edata_t*)water.data)->edgesize;
    const double xyz[3] = {low.x*mesh->ticksize, low.y*mesh->ticksize, low.z*mesh->ticksize};
    int axis = 0;
    const int bit = Global.acousticInterfaceFaceIds[f];
    if (bit > ACOUSTIC_FACE_Y_PLUS) axis = 2;
    else if (bit > ACOUSTIC_FACE_X_PLUS) axis = 1;
    for (int j = 0; j < 3; ++j) {
      const double p = xyz[j]+.5*wh*(Global.acousticInterfaceFaceCenters[3*f+j]+1);
      Global.acousticInterfaceAdjacentCenters[3*f+j] = 2*(p-reply.lower[j])/reply.h-1;
      Global.acousticInterfaceAdjacentHalfWidths[3*f+j] =
          Global.acousticInterfaceFaceHalfWidths[3*f+j]*wh/reply.h;
    }
    const double face_coordinate = Global.acousticInterfaceAdjacentCenters[3*f+axis];
    if (fabs(fabs(face_coordinate)-1) > 1.e-8)
      solver_abort(__FUNCTION_NAME, NULL, "Solid trace does not lie on a Hex8 face\n");
    Global.acousticInterfaceAdjacentCenters[3*f+axis] = face_coordinate > 0 ? 1 : -1;
    for (int c = 0; c < 8; ++c) {
      // Off-face shape functions vanish; use a valid dummy index there.
      if (bool(c & (1 << axis)) != (face_coordinate > 0)) {
        g_acousticSolidTraceNodes[8*f+c] = water.lnid[0];
        continue;
      }
      const auto& n = reply.nodes[c];
      int32_t& id = ids[n.gnid];
      if (id < 0) {
        if (n.owner == Global.myID || n.owner < 0 || n.owner >= ranks ||
            (int64_t)mesh->nharbored + (int64_t)added.size() >= INT_MAX/3)
          solver_abort(__FUNCTION_NAME, NULL, "Invalid solid trace node ownership/count\n");
        id = mesh->nharbored + (int32_t)added.size();
        node_t ghost = {};
        ghost.x = n.x; ghost.y = n.y; ghost.z = n.z;
        ghost.gnid = n.gnid; ghost.ismine = 0; ghost.isanchored = n.anchored;
        ghost.proc.ownerid = n.owner;
        added.push_back(ghost);
        shares[n.owner].push_back({n.gnid, Global.myID});
      }
      g_acousticSolidTraceNodes[8*f+c] = id;
    }
  }
  if (!added.empty()) {
    node_t* nodes = (node_t*)realloc(mesh->nodeTable,
        ((size_t)mesh->nharbored+added.size())*sizeof(node_t));
    if (!nodes) solver_abort(__FUNCTION_NAME, NULL, "Solid trace node allocation failed\n");
    mesh->nodeTable = nodes;
    std::copy(added.begin(), added.end(), nodes+mesh->nharbored);
    mesh->nharbored += (int32_t)added.size();
  }
  const auto subscriptions = acoustic_trace_exchange(shares);
  ids.clear();
  for (const auto& s : subscriptions) ids.emplace(s.gnid, -1);
  for (int32_t n = 0; n < mesh->nharbored; ++n) {
    const auto it = ids.find(mesh->nodeTable[n].gnid);
    if (it != ids.end() && mesh->nodeTable[n].ismine) it->second = n;
  }
  for (const auto& s : subscriptions) {
    const int32_t id = ids[s.gnid];
    if (id < 0) solver_abort(__FUNCTION_NAME, NULL, "Solid trace subscriber reached wrong node owner\n");
    node_t& node = mesh->nodeTable[id];
    int32link_t* link = (int32link_t*)malloc(sizeof(int32link_t));
    if (!link) solver_abort(__FUNCTION_NAME, NULL, "Solid trace sharing allocation failed\n");
    link->id = s.rank; link->next = node.proc.share; node.proc.share = link;
    g_acousticTraceShares.push_back(link);
  }
  int64_t local[2] = {(int64_t)replies.size(), (int64_t)added.size()}, total[2];
  MPI_Allreduce(local, total, 2, MPI_INT64_T, MPI_SUM, comm_solver);
  if (Global.myID == 0)
    monitor_print("Acoustic interface solid traces: %lld remote fragments, %lld additional harbored nodes; existing force/field MPI schedules.\n",
                  (long long)total[0], (long long)total[1]);
}
#endif
