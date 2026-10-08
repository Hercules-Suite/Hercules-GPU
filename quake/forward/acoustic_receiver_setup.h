#ifndef HERCULES_ACOUSTIC_RECEIVER_SETUP_H
#define HERCULES_ACOUSTIC_RECEIVER_SETUP_H
#include "acoustic_receiver_gradient.h"
#include <array>

/* Private psolve.cc implementation, after acoustic_interface_trace.h.
 * Resolve actual remote leaves once, before solver allocation, and register
 * missing vertices with the existing anchored/dangling MPI schedules. Both
 * host and CUDA-aware backends then share them after constraint assignment.
 * No mesh elements, nodal ownership or physical contributions are added.
 */
struct acoustic_receiver_request_t {
  vector3D_t point;
  int station, side, rank;
  double vp, rho;
};
struct acoustic_receiver_cell_t {
  int station, side, valid;
  acoustic_trace_node_t nodes[8];
  int32_t local_nodes[8];
  tick_t low[3], high[3];
  double h;
};
static std::vector<acoustic_receiver_cell_t> g_acousticReceiverCells;
static std::vector<int32_t> g_acousticReceiverCellSlots;
static std::vector<int32link_t*> g_acousticReceiverShares;
static std::vector<int32_t> g_acousticReceiverOffsets, g_acousticReceiverNodes;
static std::vector<double> g_acousticReceiverWeights;
static int32_t *g_acousticReceiverOffsetsDevice = NULL,
               *g_acousticReceiverNodesDevice = NULL;
static double* g_acousticReceiverWeightsDevice = NULL;

static void acoustic_receiver_prepare() {
  if (!acoustic_pressure_enabled() || !Param.acousticStationGradientRecovered ||
      Param.theNumberOfStations == 0) return;
  mesh_t* mesh = Global.myMesh;
  const int ranks = Global.theGroupSize;
  const double tick = mesh->ticksize;
  const point_t* intervals = octor_getintervaltable(Global.myOctree);
  const double extent[3] = {Param.theDomainX, Param.theDomainY, Param.theDomainZ};
  std::vector<std::vector<acoustic_receiver_request_t>> requests(ranks);
  g_acousticReceiverCellSlots.assign(Param.myNumberOfStations, -1);
  for (int s = 0; s < Param.myNumberOfStations; ++s) {
    const station_t& st = Param.myStations[s];
    if (!st.pressure_is_defined) continue;
    if (g_acousticReceiverCells.size() > (size_t)INT_MAX-6)
      solver_abort(__FUNCTION_NAME, NULL, "Acoustic receiver cell count overflow\n");
    g_acousticReceiverCellSlots[s] = (int32_t)g_acousticReceiverCells.size();
    g_acousticReceiverCells.resize(g_acousticReceiverCells.size()+6);
    const elem_t& element = mesh->elemTable[st.owningElementIndex];
    const edata_t& material = *(edata_t*)element.data;
    const node_t& n = mesh->nodeTable[element.lnid[0]];
    const double low[3] = {n.x*tick, n.y*tick, n.z*tick};
    for (int side = 0; side < 6; ++side) {
      acoustic_receiver_request_t q = {};
      q.station = s; q.side = side; q.rank = Global.myID;
      q.vp = material.Vp; q.rho = material.rho;
      const int axis = side/2;
      // Keep a face/edge receiver on the owning cell's transverse side.
      // The returned derivative is evaluated at the original coordinates.
      for (int d = 0; d < 3; ++d)
        q.point.x[d] = std::max(low[d]+.25*tick,
            std::min(low[d]+st.h-.25*tick, st.coords.x[d]));
      q.point.x[axis] = low[axis] + ((side & 1) ? st.h+.25*tick : -.25*tick);
      if (q.point.x[axis] < 0 || q.point.x[axis] >= extent[axis]) continue;
      const point_t p = {(tick_t)(q.point.x[0]/tick),
                         (tick_t)(q.point.x[1]/tick),
                         (tick_t)(q.point.x[2]/tick)};
      requests[mesh_handoff_owner(intervals, ranks, &p)].push_back(q);
    }
  }
  const auto incoming = acoustic_trace_exchange(requests);
  std::vector<std::vector<acoustic_receiver_cell_t>> answers(ranks);
  for (const auto& q : incoming) {
    acoustic_receiver_cell_t reply = {};
    reply.station = q.station; reply.side = q.side;
    int local = 0; int32_t e = -1;
    if (mesh_find_octant_at_domain_point(&q.point, reply.low, reply.high,
                                          NULL, &local, &e) && local && e >= 0) {
      const elem_t& element = mesh->elemTable[e];
      const edata_t& material = *(edata_t*)element.data;
      // Do not recover through a solid, void or discontinuous water material.
      if (element_is_acoustic(&material) &&
          fabs(material.Vp-q.vp) <= 1.e-6*std::max(1.0, fabs(q.vp)) &&
          fabs(material.rho-q.rho) <= 1.e-6*std::max(1.0, fabs(q.rho))) {
        reply.valid = 1; reply.h = material.edgesize;
        for (int c = 0; c < 8; ++c) {
          const node_t& n = mesh->nodeTable[element.lnid[c]];
          reply.nodes[c] = {n.gnid, n.x, n.y, n.z,
              n.ismine ? Global.myID : n.proc.ownerid, n.isanchored};
        }
      }
    }
    answers[q.rank].push_back(reply);
  }
  const auto replies = acoustic_trace_exchange(answers);
  std::map<int64_t, int32_t> ids;
  for (const auto& cell : replies)
    if (cell.valid) for (const auto& n : cell.nodes) ids.emplace(n.gnid, -1);
  for (int32_t n = 0; n < mesh->nharbored; ++n) {
    const auto it = ids.find(mesh->nodeTable[n].gnid);
    if (it != ids.end()) it->second = n;
  }
  std::vector<node_t> added;
  std::vector<std::vector<acoustic_trace_share_t>> shares(ranks);
  for (auto cell : replies) {
    if (!cell.valid) continue;
    for (int c = 0; c < 8; ++c) {
      const auto& n = cell.nodes[c];
      int32_t& id = ids[n.gnid];
      if (id < 0) {
        if (n.owner == Global.myID || n.owner < 0 || n.owner >= ranks ||
            (int64_t)mesh->nharbored+(int64_t)added.size() >= INT_MAX/3)
          solver_abort(__FUNCTION_NAME, NULL, "Invalid acoustic receiver halo ownership/count\n");
        id = mesh->nharbored+(int32_t)added.size();
        node_t ghost = {};
        ghost.gnid = n.gnid; ghost.x = n.x; ghost.y = n.y; ghost.z = n.z;
        ghost.isanchored = n.anchored; ghost.proc.ownerid = n.owner;
        added.push_back(ghost);
        shares[n.owner].push_back({n.gnid, Global.myID});
      }
      cell.local_nodes[c] = id;
    }
    g_acousticReceiverCells[g_acousticReceiverCellSlots[cell.station]+cell.side] = cell;
  }
  if (!added.empty()) {
    node_t* nodes = (node_t*)realloc(mesh->nodeTable,
        ((size_t)mesh->nharbored+added.size())*sizeof(node_t));
    if (!nodes) solver_abort(__FUNCTION_NAME, NULL, "Acoustic receiver halo allocation failed\n");
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
    if (id < 0) solver_abort(__FUNCTION_NAME, NULL, "Acoustic receiver halo owner missing\n");
    node_t& node = mesh->nodeTable[id];
    int32link_t* link = (int32link_t*)malloc(sizeof(int32link_t));
    if (!link) solver_abort(__FUNCTION_NAME, NULL, "Acoustic receiver halo allocation failed\n");
    link->id = s.rank; link->next = node.proc.share; node.proc.share = link;
    g_acousticReceiverShares.push_back(link);
  }
  int64_t local = added.size(), total = 0;
  MPI_Allreduce(&local, &total, 1, MPI_INT64_T, MPI_SUM, comm_solver);
  if (Global.myID == 0) monitor_print("Acoustic receiver recovery: %lld additional harbored nodes; existing host/device MPI schedules\n",
                (long long)total);
}

static bool acoustic_receiver_cell_is_pml(acoustic_receiver_cell_t& cell) {
#ifdef BOUNDARY
  return Param.theAbcType == ABC_CPML && compute_pml_face_mask(cell.low, cell.high,
      Global.myOctree->nearendp, Global.myOctree->farendp, Global.myMesh->ticksize,
      g_cpmlBuild.thickness_m) != 0;
#else
  (void)cell;
  return false;
#endif
}

static void acoustic_receiver_build_weights() {
  if (!acoustic_pressure_enabled() || Param.theNumberOfStations == 0) return;
  if (!Param.acousticStationGradientRecovered) {
    if (Global.myID == 0)
      monitor_print("Acoustic station gradient: element (historical Hex8 derivative)\n");
    return;
  }
  const double tick = Global.myMesh->ticksize;
  int64_t counts[3] = {0}, totals[3]; // recovered axes, fallback axes, PML stations
  g_acousticReceiverOffsets.push_back(0);
  for (int s = 0; s < Param.myNumberOfStations; ++s) {
    const station_t& st = Param.myStations[s];
    // Ordered by global ID so summation order is independent of partitioning.
    std::map<int64_t, std::pair<int32_t, std::array<double,3>>> stencil;
    if (st.pressure_is_defined) {
      const bool pml = cpml_pml_element(st.owningElementIndex);
      counts[2] += pml;
      for (int axis = 0; axis < 3 && !pml; ++axis) {
        const int preferred = 2*axis+(st.localcoords.x[axis] >= 0 ? 1 : 0);
        acoustic_receiver_cell_t* cell = NULL;
        if (!pml) for (int attempt = 0; attempt < 2; ++attempt) {
          auto& candidate = g_acousticReceiverCells[g_acousticReceiverCellSlots[s]+(preferred ^ attempt)];
          if (candidate.valid && !acoustic_receiver_cell_is_pml(candidate)) {
            cell = &candidate; break;
          }
        }
        double fraction = 0, other_local[3] = {}, own[8], neighbor[8] = {};
        if (cell) {
          const double center = st.coords.x[axis]-.5*st.h*st.localcoords.x[axis];
          const double other_center = .5*(double(cell->low[axis])+cell->high[axis])*tick;
          fraction = acoustic_receiver::neighbor_fraction(st.coords.x[axis], center, other_center);
          for (int d = 0; d < 3; ++d)
            other_local[d] = 2*(st.coords.x[d]-cell->low[d]*tick)/cell->h-1;
          // Normal coordinates may lie outside the adjacent element; its
          // normal derivative is independent of that coordinate.
          for (int d = 0; d < 3; ++d) if (d != axis && fabs(other_local[d]) > 1+1.e-8)
            solver_abort(__FUNCTION_NAME, NULL, "Acoustic receiver neighbor misses transverse point\n");
          acoustic_receiver::gradient_weights(axis, cell->h, other_local, neighbor);
          ++counts[0];
        } else ++counts[1];
        acoustic_receiver::gradient_weights(axis, st.h, st.localcoords.x, own);
        for (int c = 0; c < 8; ++c) {
          const int32_t node = st.nodestointerpolate[c];
          auto& entry = stencil[Global.myMesh->nodeTable[node].gnid];
          entry.first = node; entry.second[axis] += (1-fraction)*own[c];
          if (cell && fraction != 0) {
            const int32_t other = cell->local_nodes[c];
            auto& extra = stencil[Global.myMesh->nodeTable[other].gnid];
            extra.first = other; extra.second[axis] += fraction*neighbor[c];
          }
        }
      }
    }
    for (const auto& item : stencil) {
      g_acousticReceiverNodes.push_back(item.second.first);
      for (double w : item.second.second) g_acousticReceiverWeights.push_back(w);
    }
    if (g_acousticReceiverNodes.size() > (size_t)INT_MAX)
      solver_abort(__FUNCTION_NAME, NULL, "Acoustic receiver stencil count overflow\n");
    g_acousticReceiverOffsets.push_back((int32_t)g_acousticReceiverNodes.size());
  }
  MPI_Allreduce(counts, totals, 3, MPI_INT64_T, MPI_SUM, comm_solver);
  if (Global.myID == 0) monitor_print("Acoustic station gradient: linear; %lld reconstructed / %lld element-fallback components, %lld CPML stations retain element gradients; pressure interpolation unchanged\n",
      (long long)totals[0], (long long)totals[1], (long long)totals[2]);
  std::vector<acoustic_receiver_cell_t>().swap(g_acousticReceiverCells);
  std::vector<int32_t>().swap(g_acousticReceiverCellSlots);
}

static void acoustic_receiver_free() {
  cudaFree(g_acousticReceiverOffsetsDevice); g_acousticReceiverOffsetsDevice = NULL;
  cudaFree(g_acousticReceiverNodesDevice); g_acousticReceiverNodesDevice = NULL;
  cudaFree(g_acousticReceiverWeightsDevice); g_acousticReceiverWeightsDevice = NULL;
  g_acousticReceiverCells.clear(); g_acousticReceiverCellSlots.clear();
  g_acousticReceiverOffsets.clear();
  g_acousticReceiverNodes.clear(); g_acousticReceiverWeights.clear();
  for (auto* link : g_acousticReceiverShares) free(link);
  g_acousticReceiverShares.clear();
}
#endif
