/* -*- C -*- */

/* @copyright_notice_start
 *
 * This file is part of the CMU Hercules ground motion simulator developed
 * by the CMU Quake project.
 *
 * Copyright (C) Carnegie Mellon University. All rights reserved.
 *
 * This program is covered by the terms described in the 'LICENSE.txt' file
 * included with this software package.
 *
 * This program comes WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * 'LICENSE.txt' file for more details.
 *
 *  @copyright_notice_end
 */

/*
 * Generate an unstructured mesh and solve the linear system thereof
 * derived.
 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <float.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "util.h"
#include "commutil.h"
#include "timers.h"

#ifdef __cplusplus
extern "C" {
#endif
#include "etree.h"
#include "octor.h"
#include "cvm.h"
#ifdef __cplusplus
}
#endif
#include "psolve.h"
#include "nrutila.h"
#include "quakesource.h"
#include "output.h"
#include "nonlinear.h"
#include "io_planes.h"
#include "io_checkpoint.h"
#include "stiffness.h"
#include "damping.h"
#include "quake_util.h"
#include "buildings.h"
#include "drm.h"
#include "meshformatlab.h"

#include "kernel.h"

#include "topography.h"
#include "drm_planewaves.h"
#include "basin.h"

#include <cuda.h>
#include <cuda_runtime.h>


/* ONLY GLOBAL VARIABLES ALLOWED OUTSIDE OF PARAM. and GLOBAL. IN ALL OF
 * PSOLVE!! */
MPI_Comm comm_solver;
MPI_Comm comm_IO;

static int32_t g_dampingElemCountGPU = -1;
static int32_t* g_dampingElemMapperDevice = NULL;
static int32_t g_topoElemCountGPU = -1;
static int32_t* g_topoElemMapperDevice = NULL;
static topoconstants_t* g_topoConstantsDevice = NULL;
static int32_t g_topoForceMode = -1;
static int32_t g_topoForceModeReported = 0;
static int32_t g_topoForcePathReported = 0;
static int64_t g_acousticPressureElementCount = 0;
static int64_t g_acousticPressureNodeCount = 0;
static int64_t g_acousticPressureFreeSurfaceNodeCount = 0;
static int32_t g_dampingBkt2CoeffsReady = 0;
static double g_dampingBkt2ShearCoef1 = 0.0;
static double g_dampingBkt2ShearCoef2 = 0.0;
static double g_dampingBkt2ShearExpCoef0 = 0.0;
static double g_dampingBkt2ShearCoef3 = 0.0;
static double g_dampingBkt2ShearCoef4 = 0.0;
static double g_dampingBkt2ShearExpCoef1 = 0.0;
static double g_dampingBkt2KappaCoef1 = 0.0;
static double g_dampingBkt2KappaCoef2 = 0.0;
static double g_dampingBkt2KappaExpCoef0 = 0.0;
static double g_dampingBkt2KappaCoef3 = 0.0;
static double g_dampingBkt2KappaCoef4 = 0.0;
static double g_dampingBkt2KappaExpCoef1 = 0.0;
static int32_t* g_sourceNodeListDevice = NULL;
static vector3D_t* g_sourceForcesDevice = NULL;
static int32_t g_sourceNodeCapacity = 0;
static int32_t* g_sparseNodeListDevice = NULL;
static fvector_t* g_sparseForcesDevice = NULL;
static int32_t g_sparseNodeCapacity = 0;
static int32_t g_sparseActiveCount = 0;
static int32_t* g_sparseNodeListHost = NULL;
static fvector_t* g_sparseForcesHost = NULL;
static int32_t g_sparseHostPinned = 0;
static int32_t g_sparseNodeListReady = 0;
static int32_t g_sparseSourceCountCached = -1;
static int32_t g_sparseTopoCountCached = -1;
static int32_t g_sparseCopiedSourceCount = -1;
static int32_t g_sparseCopiedTopoCount = -1;
static int32_t* g_rangeUploadStartsHost = NULL;
static int32_t* g_rangeUploadLengthsHost = NULL;
static int32_t g_rangeUploadCapacity = 0;
static int32_t g_rangeUploadCount = 0;
static int64_t g_rangeUploadNodeCount = 0;
static int32_t g_rangeUploadReady = 0;
static int32_t g_rangeUploadSourceCountCached = -1;
static int32_t g_rangeUploadTopoCountCached = -1;
static int32_t g_sparseClearMode = -1;
static int32_t g_sparseClearModeReported = 0;
static int32_t g_sparseUploadMode = -1;
static int32_t g_sparseUploadModeReported = 0;
static int32_t g_sparseUploadPathReported = -1;
static vector3D_t g_istanbulVelModelOrigin;
static int32_t g_hasIstanbulVelModelOrigin = 0;
static int32_t* g_adjustParentOffsets = NULL;
static int32_t* g_adjustParentIds = NULL;
static int32_t* g_adjustDanglingLnid = NULL;
static int32_t* g_adjustParentBase3 = NULL;
static int32_t* g_adjustDanglingBase3 = NULL;
static solver_float* g_adjustInvDeps = NULL;
static int32_t* g_adjustDnOrder = NULL;
static int32_t* g_adjustDnSortKey = NULL;
static int32_t* g_adjustDistEdgeDn = NULL;
static int32_t* g_adjustDistEdgeParentBase3 = NULL;
static int32_t* g_adjustDistRunParentBase3 = NULL;
static int32_t* g_adjustDistRunOffsets = NULL;
static int32_t g_adjustDistRunCount = 0;
static int32_t g_adjustDnCount = 0;
static int32_t g_adjustParentCount = 0;
static int32_t g_adjustDistributionOmpMode = -1;
static int32_t g_adjustAssignmentOmpMode = -1;
static int32_t g_adjustDistributionOmpReported = 0;
static int32_t g_adjustAssignmentOmpReported = 0;
static int32_t g_stationOutputOmpMode = -1;
static int32_t g_stationOutputOmpReported = 0;
static int32_t g_stationIntervalFilterMode = -1;
static int32_t g_stationIntervalFilterReported = 0;
static int32_t g_fullClearOmpMode = -1;
static int32_t g_fullClearOmpReported = 0;
static int32_t g_contribPackOmpMode = -1;
static int32_t g_contribPackOmpReported = 0;
static int32_t g_contribUnpackOmpMode = -1;
static int32_t g_contribUnpackOmpReported = 0;
static int32_t g_contribRecvMode = -1;
static int32_t g_contribRecvReported = 0;
static int32_t g_reproMode = -1;
static int32_t g_meshImportBootstrapLevel = -1;
static int32_t g_forcePathReported = 0;
static cudaEvent_t g_forceUploadSourceSafeEvent;
static int32_t g_forceUploadSourceSafeEventInit = 0;
static int32_t g_forceUploadSourceSafePending = 0;
static int32_t g_dispSyncMode = -1;
static int32_t g_dispSyncReported = 0;
static int32_t g_dispBoundaryReady = 0;
static int32_t* g_dispBoundaryNodesHost = NULL;
static int32_t* g_dispBoundaryNodesDevice = NULL;
static fvector_t* g_dispBoundaryBufferHost = NULL;
static fvector_t* g_dispBoundaryBufferDevice = NULL;
static int32_t g_dispBoundaryNodeCount = 0;

#define DN_DIAG_METRIC_COUNT 9
enum {
  DN_DIAG_SEND_BYTES = 0,
  DN_DIAG_RECV_BYTES,
  DN_DIAG_SEND_NODES,
  DN_DIAG_RECV_NODES,
  DN_DIAG_SEND_MSGS,
  DN_DIAG_RECV_MSGS,
  DN_DIAG_ADJUST_NODES,
  DN_DIAG_ADJUST_DEPS,
  DN_DIAG_STEPS
};

static double g_dnDiagLocal[DN_DIAG_METRIC_COUNT] = {0.0};
static double g_dnDiagMin[DN_DIAG_METRIC_COUNT] = {0.0};
static double g_dnDiagMax[DN_DIAG_METRIC_COUNT] = {0.0};
static double g_dnDiagAvg[DN_DIAG_METRIC_COUNT] = {0.0};
static int32_t g_dnDiagReady = 0;

#define WAIT_RECV_STAGE_COUNT 4
enum {
  WAIT_RECV_STAGE_1_DN_FORCE = 0,
  WAIT_RECV_STAGE_2_AN_FORCE,
  WAIT_RECV_STAGE_3_AN_DISP,
  WAIT_RECV_STAGE_4_DN_DISP
};

static const char* g_waitRecvStageNames[WAIT_RECV_STAGE_COUNT] = {
    "1st send data (contribution)",
    "2nd send data (contribution)",
    "3rd send data (sharing)",
    "4th send data (sharing)"};

static double g_waitRecvStageLocal[WAIT_RECV_STAGE_COUNT] = {0.0};
static double g_waitRecvStageMin[WAIT_RECV_STAGE_COUNT] = {0.0};
static double g_waitRecvStageMax[WAIT_RECV_STAGE_COUNT] = {0.0};
static double g_waitRecvStageAvg[WAIT_RECV_STAGE_COUNT] = {0.0};
static int32_t g_waitRecvStageMaxRank[WAIT_RECV_STAGE_COUNT] = {-1, -1, -1,
                                                                 -1};

#define STAGE_COMM_METRIC_COUNT (WAIT_RECV_STAGE_COUNT * 6)
enum {
  STAGE_COMM_SEND_BYTES = 0,
  STAGE_COMM_RECV_BYTES,
  STAGE_COMM_SEND_MSGS,
  STAGE_COMM_RECV_MSGS,
  STAGE_COMM_SEND_NODES,
  STAGE_COMM_RECV_NODES
};

static double g_stageCommLocal[STAGE_COMM_METRIC_COUNT] = {0.0};
static double g_stageCommMin[STAGE_COMM_METRIC_COUNT] = {0.0};
static double g_stageCommMax[STAGE_COMM_METRIC_COUNT] = {0.0};
static double g_stageCommAvg[STAGE_COMM_METRIC_COUNT] = {0.0};

#define STAGE_TIME_METRIC_COUNT_PER_STAGE 6
#define STAGE_TIME_METRIC_COUNT \
  (WAIT_RECV_STAGE_COUNT * STAGE_TIME_METRIC_COUNT_PER_STAGE)
enum {
  STAGE_TIME_POST_RECV = 0,
  STAGE_TIME_PACK_OUTGOING,
  STAGE_TIME_POST_SEND,
  STAGE_TIME_WAIT_RECV,
  STAGE_TIME_WAIT_SEND,
  STAGE_TIME_UNPACK_INCOMING
};

static const char* g_stageTimeMetricNames[STAGE_TIME_METRIC_COUNT_PER_STAGE] = {
    "post recv",
    "pack outgoing",
    "post send",
    "wait recv",
    "wait send",
    "unpack incoming"};

static double g_stageTimeLocal[STAGE_TIME_METRIC_COUNT] = {0.0};
static double g_stageTimeMin[STAGE_TIME_METRIC_COUNT] = {0.0};
static double g_stageTimeMax[STAGE_TIME_METRIC_COUNT] = {0.0};
static double g_stageTimeAvg[STAGE_TIME_METRIC_COUNT] = {0.0};

#ifndef PROCPERNODE
#define PROCPERNODE 4
#endif

#define PI 3.14159265358979323846264338327

#define GOAHEAD_MSG 100
#define MESH_MSG 101
#define STAT_MSG 102
#define OUT4D_MSG 103
#define DN_MASS_MSG 104
#define AN_MASS_MSG 105
#define DN_FORCE_MSG 106
#define AN_FORCE_MSG 107
#define DN_DISP_MSG 108
#define AN_DISP_MSG 109
#define DN_PRESSURE_FORCE_MSG 120
#define AN_PRESSURE_FORCE_MSG 121
#define AN_PRESSURE_MSG 122
#define DN_PRESSURE_MSG 123
#define DN_PRESSURE_MASS_MSG 124
#define AN_PRESSURE_MASS_MSG 125
#define DN_PRESSURE_SURFACE_MSG 126
#define AN_PRESSURE_SURFACE_MSG 127
#define DN_PRESSURE_DAMP_MSG 128
#define AN_PRESSURE_DAMP_MSG 129
#define CVMRECORD_MSG 110

#define BATCH (1 << 20)
#define LINESIZE 512
#define FILEBUFSIZE (1 << 25)

#define CONTRIBUTION 901 /**< Harboring processors to owner processors */
#define SHARING 902      /**< Owner processors to harboring processors */

#define DISTRIBUTION 903 /**< Dangling nodes to anchored nodes */
#define ASSIGNMENT 904   /**< Anchored nodes to dangling nodes */

/*---------------Initialization and cleanup routines----------------------*/
static void read_parameters(int argc, char** argv);
static int32_t parse_parameters(const char* numericalin);
static void local_finalize(void);

/*---------------- Mesh generation data structures -----------------------*/
#ifdef USECVMDB

#ifndef CVMBUFSIZE
#define CVMBUFSIZE 100
#endif

static void replicateDB(const char* dbname);
static void open_cvmdb(void);

#else

static int32_t zsearch(void* base, int32_t count, int32_t recordsize,
                       const point_t* searchpt);
static cvmrecord_t* sliceCVM(const char* cvm_flatfile);

#endif
/** cvmrecord_t: cvm record.  Only used if USECVMDB not defined **/
typedef struct cvmrecord_t {
  char key[12];
  float Vp, Vs, density;
} cvmrecord_t;

/**
 * mrecord_t: Complete mesh database record
 *
 */
typedef struct mrecord_t {
  etree_addr_t addr;
  mdata_t mdata;
} mrecord_t;

/* Mesh generation related routines */
static inline int monitor_print(const char* format, ...);
static int32_t toexpand(octant_t* leaf, double ticksize, const void* data);
static void setrec(octant_t* leaf, double ticksize, void* data);
static void mesh_generate(void);
static int32_t toexpand_mesh_etree(octant_t* leaf, double ticksize,
                                   const void* data);
static void setrec_mesh_etree(octant_t* leaf, double ticksize, void* data);
static void mesh_import_apply_properties(void);
static void mesh_generate_from_etree(void);
static int32_t bulkload(etree_t* mep, mrecord_t* partTable, int32_t count);
static void mesh_output(void);
static void mesh_correct_properties(etree_t* cvm);
static void mesh_release_acoustic_interface_data(void);
static void mesh_build_acoustic_interface_flags(
    const vector3D_t* IstVelModel_origin, const vector3D_t* basinVelModel_origin,
    const vector3D_t* threeDVelModel_origin);
static void mesh_write_acoustic_interface_diagnostics(void);

static void gpu_init(int32_t myID);
static void solver_init_damping_gpu_mapper(int32_t myID, mesh_t* mesh);
static void solver_init(void);
static void solver_printstat(mysolver_t* solver);
static void solver_delete(void);
static void solver_run(void);
void solver_output_seq(void);
static int solver_print_schedules(mysolver_t* solver);

static schedule_t* schedule_new(void);
static void schedule_build(mesh_t* mesh, schedule_t* dnsched,
                           schedule_t* ansched);
static void schedule_allocMPIctl(schedule_t* sched);
static void schedule_allocmapping(schedule_t* sched);
static void schedule_delete(schedule_t* sched);
static void schedule_prepare(schedule_t* sched, int32_t c_outsize,
                             int32_t c_insize, int32_t s_outsize,
                             int32_t s_insize);
static void schedule_senddata(schedule_t* sched, void* valuetable,
                              int32_t itemsperentry, int32_t direction,
                              int32_t msgtag);
static void schedule_unpack_incoming_message(messenger_t* recv_messenger,
                                             void* valuetable,
                                             int32_t itemsperentry,
                                             int32_t direction);
static int schedule_print(schedule_t* sched, char type, FILE* out);
static int schedule_print_detail(schedule_t* sched, char type, FILE* out);
static int schedule_print_messenger_list(schedule_t* sched, messenger_t* msg,
                                         int count, char type, char cs,
                                         FILE* out);

static messenger_t* messenger_new(int32_t procid);
static void messenger_delete(messenger_t* messenger);
static void messenger_set(messenger_t* messenger, int32_t outsize,
                          int32_t insize);
static int32_t messenger_countnodes(messenger_t* first);
static void messenger_prepare_offsets(messenger_t* first,
                                      int32_t itemsperentry);
static int compare_int32_asc(const void* a, const void* b);
static int32_t messenger_detect_duplicate_offsets(messenger_t* messenger);

static void compute_K(void);
static void constract_Quality_Factor_Table(void);

#ifdef BOUNDARY
static char compute_setflag(tick_t ldb[3], tick_t ruf[3], tick_t nearendp[3],
                            tick_t farendp[3]);
static void compute_setboundary(float size, float Vp, float Vs, float rho,
                                int flag, double dashpot[8][3]);
static uint8_t compute_pml_face_mask(tick_t ldb[3], tick_t ruf[3],
                                     tick_t nearendp[3], tick_t farendp[3],
                                     double ticksize, double pml_thickness_m);
static void compute_pml_profile_coeffs(uint8_t face_mask, tick_t ldb[3],
                                       tick_t ruf[3], tick_t nearendp[3],
                                       tick_t farendp[3], double ticksize,
                                       double pml_thickness_m, double Vp,
                                       double dt, int poly_order, double R,
                                       double kappa_max, double alpha_max,
                                       double mpml_cross_fraction,
                                       double alpha_floor_fraction,
                                       double a_out[3], double b_out[3],
                                       double invKappa_out[3],
                                       double d_out[3],
                                       double alpha_out[3],
                                       double kappa_out[3]);
#endif /* BOUNDARY */

/* CPML build-time context. Populated by solver_init() while iterating
 * elements; the contents are uploaded into Global.gpuData and then freed
 * (except elem_list, which is moved into mySolver for checkpoint I/O).
 *
 * PML elements are partitioned by material type (elastic vs acoustic) at
 * setup time using element_is_acoustic(edata). Each side has its own
 * mapper, ADE coefficient arrays, and (later) memory variables. The build
 * arrays are sized to the upper bound (lenum) for simplicity; only the
 * first elastic_count / acoustic_count entries are populated and uploaded. */
typedef struct cpml_build_s {
  uint8_t* face_mask;            /* [lenum], NULL when abc_type != cpml */

  /* Elastic PML build state. The elastic PML runs as native effective
   * stiffness + a CPML correction in the native (c1,c2) normalization
   * (mu=9*c1, lambda=9*c2, dN=0.25*sign), so:
   *   a[j]        = raw a_j        (no 1/(4h) prescale)
   *   invKappa[j] = raw 1/kappa_j  (no 1/(4h) prescale)
   *   b[j]        = exp() decay, dimensionless (unchanged)
   *   muCoeff/lambdaCoeff = unused now (kept allocated; cleanup later)
   */
  int32_t  elastic_count;
  int32_t* elastic_elem_list;    /* sized [lenum] during build */
  double*  elastic_a[3];         /* [lenum] each, used [0..elastic_count) */
  double*  elastic_b[3];
  double*  elastic_invKappa[3];
  double*  elastic_muCoeff;
  double*  elastic_lambdaCoeff;
  /* CFS divergence-side coefficients (Matzen 2nd-order form), per axis. */
  double*  elastic_zetaDecay[3]; /* exp(-alpha_a*dt) */
  double*  elastic_zetaGain[3];  /* (1-exp(-alpha_a*dt))/alpha_a */
  double*  elastic_memForce[3];  /* dt^2 * (rho*h^3/8) * d_a * alpha_a^2 */
  double*  elastic_restoreCoeff; /* dt^2 * (rho*h^3/8) * sum_a d_a*alpha_a */

  /* Acoustic PML build state. Same pre-scaling for a / invKappa.
   *   forceCoeff = dt^2 * h^2 / (4 * rho)
   */
  int32_t  acoustic_count;
  int32_t* acoustic_elem_list;
  double*  acoustic_a[3];
  double*  acoustic_b[3];
  double*  acoustic_invKappa[3];
  double*  acoustic_forceCoeff;
  /* Step B divergence-side CFS coeffs (scalar mirror of elastic;
   * memForce/restoreCoeff use the acoustic nodal mass volume/(8*bulk)
   * in place of the elastic rho*h^3/8). */
  double*  acoustic_zetaDecay[3]; /* exp(-alpha_a*dt) */
  double*  acoustic_zetaGain[3];  /* (1-exp(-alpha_a*dt))/alpha_a */
  double*  acoustic_memForce[3];  /* dt^2 * nodalPressureMass * d_a*alpha_a^2 */
  double*  acoustic_restoreCoeff; /* dt^2 * nodalPressureMass * sum d_a*alpha_a */

  double   thickness_m;          /* realized PML thickness */
  double   alpha_max;            /* alpha_max used for profiles */
} cpml_build_t;
static cpml_build_t g_cpmlBuild = {0};

static void compute_setab(double freq, double* aBasePtr, double* bBasePtr);
static int32_t mesh_element_index_for_octant(octant_t* octant,
                                             const char* context);
static void compute_addforce_s(int32_t timestep);
static void solver_compute_force_source_gpu(int32_t timestep);
static void solver_compute_force_acoustic_pressure(mysolver_t* solver,
                                                   mesh_t* mesh);
static void solver_run_acoustic_pressure_step(mysolver_t* solver, mesh_t* mesh);
static void solver_compute_force_topography_gpu(mysolver_t* solver,
                                                mesh_t* mesh,
                                                double deltaT2);
static void solver_load_forces_sparse_gpu(mysolver_t* solver, mesh_t* mesh);
static void solver_load_forces_range_gpu(mysolver_t* solver, mesh_t* mesh,
                                         int32_t sourceCount);
static void solver_clear_sparse_host_forces(mysolver_t* solver);
static int32_t solver_sparse_clear_mode(void);
static int32_t use_gpu_topography_fastpath(void);
static int32_t topography_has_host_force_path(void);
static void solver_report_force_paths_once(void);
static int32_t need_host_force_upload(void);
static void solver_init_topography_gpu_mapper(int32_t myID, mesh_t* mesh);
static void compute_adjust(void* valuetable, int32_t itemsperentry,
                           int32_t how);
static void compute_adjust_cache_build(void);
static void compute_adjust_cache_free(void);
static void solver_diag_reset(void);
static void solver_diag_prepare_local(mesh_t* mesh);
static void solver_diag_reduce(void);
static int32_t wait_recv_stage_from_msgtag(int32_t msgtag);
static int32_t stage_comm_index(int32_t stage, int32_t metric);
static int32_t stage_time_index(int32_t stage, int32_t metric);
static int32_t adjust_distribution_use_omp(void);
static int32_t adjust_assignment_use_omp(void);
static int32_t station_output_use_omp(void);
static int32_t station_interval_filter_enabled(void);
static int32_t station_point_in_local_interval(vector3D_t point);
static int32_t full_clear_use_omp(void);
static int32_t contrib_pack_use_omp(void);
static int32_t contrib_unpack_use_omp(void);
static int32_t contrib_recv_use_waitsome(void);
static int32_t sparse_force_upload_mode(void);
static int32_t hercules_repro_mode(void);
static int32_t hercules_repro_strict(void);
static int32_t displacement_sync_mode(void);
static void solver_init_disp_boundary_sync(mysolver_t* solver);
static void solver_finalize_disp_boundary_download(mysolver_t* solver);
static void solver_upload_disp_boundary_gpu(mysolver_t* solver);
static void solver_record_force_upload_source_event(mysolver_t* solver);
static void solver_wait_force_upload_source_event(void);
static void solver_record_gpu_event(mysolver_t* solver, int eventID,
                                    const char* eventLabel);
static int compare_adjust_dn_order(const void* a, const void* b);
static int compare_adjust_dist_edge(const void* a, const void* b);
static void solver_snapshot_station_accel_history(mysolver_t* solver);
static void solver_update_acoustic_station_kinematics(int32_t step);
static int station_pressure_gradient(int32_t iStation, double grad[3]);
typedef struct station_medium_record_t station_medium_record_t;
static void station_write_mediums_file(
    const station_medium_record_t* localRecords);
static void solver_prepare_force_upload_ranges(mesh_t* mesh,
                                               int32_t sourceCount);

static int interpolate_station_displacements(int32_t step);
static double interpolate_station_pressure_value(int32_t iStation);

__global__ static void kernelSourceForceScatter(int32_t count,
                                                const int32_t* lnidList,
                                                const vector3D_t* srcForces,
                                                fvector_t* forceDevice,
                                                double dt2) {
  int32_t idx = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (idx >= count) {
    return;
  }

  int32_t lnid = lnidList[idx];
  const vector3D_t src = srcForces[idx];
  fvector_t out;
  out.f[0] = src.x[0] * dt2;
  out.f[1] = src.x[1] * dt2;
  out.f[2] = src.x[2] * dt2;
  forceDevice[lnid] = out;
}

__global__ static void kernelSparseForceScatter(int32_t count,
                                                const int32_t* lnidList,
                                                const fvector_t* srcForces,
                                                fvector_t* forceDevice) {
  int32_t idx = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (idx >= count) {
    return;
  }

  int32_t lnid = lnidList[idx];
  forceDevice[lnid] = srcForces[idx];
}

__global__ static void kernelDispBoundaryGather(int32_t count,
                                                const int32_t* lnidList,
                                                const fvector_t* dispDevice,
                                                fvector_t* bufferDevice) {
  int32_t idx = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (idx >= count) {
    return;
  }

  bufferDevice[idx] = dispDevice[lnidList[idx]];
}

__global__ static void kernelDispBoundaryScatter(int32_t count,
                                                 const int32_t* lnidList,
                                                 const fvector_t* bufferDevice,
                                                 fvector_t* dispDevice) {
  int32_t idx = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (idx >= count) {
    return;
  }

  dispDevice[lnidList[idx]] = bufferDevice[idx];
}

__device__ __constant__ int kTopoTetNodes[2][5][4] = {
    {{0, 2, 1, 4}, {3, 1, 2, 7}, {6, 4, 7, 2},
     {5, 7, 4, 1}, {2, 4, 7, 1}},
    {{0, 3, 1, 5}, {0, 2, 3, 6}, {4, 5, 6, 0},
     {6, 5, 7, 3}, {0, 6, 3, 5}}};

__device__ __constant__ double kTopoCubeCoord[8][3] = {
    {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
    {1.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0},
    {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}};

__device__ static double topo_det3(const double a[3][3]) {
  return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
         a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
         a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
}

__device__ static void topo_inv3(const double a[3][3], double inv[3][3]) {
  double d = topo_det3(a);

  inv[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / d;
  inv[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) / d;
  inv[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / d;
  inv[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) / d;
  inv[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) / d;
  inv[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) / d;
  inv[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) / d;
  inv[2][1] = (a[0][1] * a[2][0] - a[0][0] * a[2][1]) / d;
  inv[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / d;
}

__device__ static void topoTetraForcesEffective(const fvector_t* un,
                                                fvector_t* resVec,
                                                const topoconstants_t* topo_ec) {
  const int part = (topo_ec->cube_part == 1) ? 0 : 1;
  const double h = topo_ec->h;

  for (int k = 0; k < 5; k++) {
    if (topo_ec->tetraVol[k] == 0.0) {
      continue;
    }

    int nodes[4];
    double x[4][3];
    double a[3][3];
    double inv[3][3];
    double grad[4][3];
    double gradU[3][3];
    double sigma[3][3];

    for (int i = 0; i < 4; i++) {
      nodes[i] = kTopoTetNodes[part][k][i];
    }

    for (int i = 0; i < 4; i++) {
      int n = nodes[i];
      x[i][0] = kTopoCubeCoord[n][0] * h;
      x[i][1] = kTopoCubeCoord[n][1] * h;
      x[i][2] = kTopoCubeCoord[n][2] * h;
    }

    for (int col = 0; col < 3; col++) {
      for (int row = 0; row < 3; row++) {
        a[row][col] = x[col + 1][row] - x[0][row];
      }
    }

    topo_inv3(a, inv);

    for (int i = 1; i < 4; i++) {
      grad[i][0] = inv[i - 1][0];
      grad[i][1] = inv[i - 1][1];
      grad[i][2] = inv[i - 1][2];
    }
    for (int d = 0; d < 3; d++) {
      grad[0][d] = -(grad[1][d] + grad[2][d] + grad[3][d]);
    }

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        gradU[i][j] = 0.0;
      }
    }

    for (int aidx = 0; aidx < 4; aidx++) {
      int n = nodes[aidx];
      for (int comp = 0; comp < 3; comp++) {
        for (int d = 0; d < 3; d++) {
          gradU[comp][d] += un[n].f[comp] * grad[aidx][d];
        }
      }
    }

    double trace = gradU[0][0] + gradU[1][1] + gradU[2][2];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        sigma[i][j] = topo_ec->mu * (gradU[i][j] + gradU[j][i]);
        if (i == j) {
          sigma[i][j] += topo_ec->lambda * trace;
        }
      }
    }

    double volume = fabs(topo_det3(a)) * topo_ec->tetraVol[k] / 6.0;
    for (int aidx = 0; aidx < 4; aidx++) {
      int n = nodes[aidx];
      for (int comp = 0; comp < 3; comp++) {
        double value = 0.0;
        for (int d = 0; d < 3; d++) {
          value += sigma[comp][d] * grad[aidx][d];
        }
        resVec[n].f[comp] += (solver_float)(volume * value);
      }
    }
  }
}

__global__ static void kernelTopoForceEffective(int32_t lenum,
                                                gpu_data_t* gpuData,
                                                int32_t topoElemCount,
                                                const int32_t* topoElemMapper,
                                                const topoconstants_t* topoConstants,
                                                double dt2) {
  int32_t topoIdx = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (topoIdx >= topoElemCount) {
    return;
  }

  int32_t eindex = topoElemMapper[topoIdx];
  topoconstants_t topo_ec = topoConstants[topoIdx];
  fvector_t localForce[8];
  fvector_t curDisp[8];
  double c1 = gpuData->c1ArrayDevice[eindex];
  double b_over_dt;
  int do_contrib = 0;

  if (c1 == 0.0) {
    return;
  }

  b_over_dt = gpuData->c3ArrayDevice[eindex] / c1;

  for (int i = 0; i < 8; i++) {
    int32_t lnid = gpuData->lnidArrayDevice[i * lenum + eindex];
    const fvector_t* tm1Disp = gpuData->tm1Device + lnid;
    const fvector_t* tm2Disp = gpuData->tm2Device + lnid;

    localForce[i].f[0] = 0;
    localForce[i].f[1] = 0;
    localForce[i].f[2] = 0;

    curDisp[i].f[0] =
        tm1Disp->f[0] * (1.0 + b_over_dt) - b_over_dt * tm2Disp->f[0];
    curDisp[i].f[1] =
        tm1Disp->f[1] * (1.0 + b_over_dt) - b_over_dt * tm2Disp->f[1];
    curDisp[i].f[2] =
        tm1Disp->f[2] * (1.0 + b_over_dt) - b_over_dt * tm2Disp->f[2];

    do_contrib =
        (fabs(curDisp[i].f[0]) > 1e-20) ? 1 : do_contrib;
    do_contrib =
        (fabs(curDisp[i].f[1]) > 1e-20) ? 1 : do_contrib;
    do_contrib =
        (fabs(curDisp[i].f[2]) > 1e-20) ? 1 : do_contrib;
  }

  if (!do_contrib) {
    return;
  }

  topoTetraForcesEffective(curDisp, localForce, &topo_ec);

  for (int i = 0; i < 8; i++) {
    int32_t lnid = gpuData->lnidArrayDevice[i * lenum + eindex];
    fvector_t* nodalForce = gpuData->forceDevice + lnid;

    atomicAdd((solver_float*)&(nodalForce->f[0]),
              (solver_float)(-localForce[i].f[0] * dt2));
    atomicAdd((solver_float*)&(nodalForce->f[1]),
              (solver_float)(-localForce[i].f[1] * dt2));
    atomicAdd((solver_float*)&(nodalForce->f[2]),
              (solver_float)(-localForce[i].f[2] * dt2));
  }
}

double** malloc2dDouble(int n, int m) {
  /* NOTE: To MPI_Bcast a 2D array, the pre-allocated memory blocks have to be
  continuous. That's why this function should be used to do the pre-allocation.
*/
  static const char* fname = __FUNCTION_NAME;
  /* allocate the n*m contiguous items */
  double* p = (double*)malloc(n * m * sizeof(double));
  if (p == NULL) {
    solver_abort(fname, NULL,
                 "Error allocating memory for the n*m contiguous items\n");
  }
  /* allocate the row pointers into the memory */
  double** array = (double**)malloc(n * sizeof(double*));
  if (array == NULL) {
    free(p);
    solver_abort(fname, NULL, "Error allocating memory for the row pointers\n");
  }
  /* set up the pointers into the contiguous memory */
  for (int i = 0; i < n; i++) array[i] = &(p[i * m]);
  return array;
}

/* ---------- Reading parameters from files ---------- */
double** readBasin(char* filename, int* x_count, int* y_count, int* z_count,
                   double* increment) {
  FILE* fp;
  int i;
  // Open the file in read mode
  fp = fopen(filename, "r");
  if (fp == NULL) {
    printf("Error opening file %s.\n", filename);
    exit(1);
  }
  // Read the first line to get the counts
  fscanf(fp, "%d %d %d", x_count, y_count, z_count);
  // Calculate the total number of lines in the file
  int total_lines = (*x_count) * (*y_count) * (*z_count);
  // Create a 2D array to store the data
  double** data = malloc2dDouble(total_lines, 3);
  // Read the data and store it in the array
  for (i = 0; i < total_lines; i++) {
    fscanf(fp, "%lf %lf %lf", &data[i][0], &data[i][1], &data[i][2]);
  }
  // Assuming the increments for x and y values in the 'data' array are the same
  *increment = data[1][0] - data[0][0];
  // Close the file
  fclose(fp);
  return data;
}

/* ---------- Creating the directory if it does not exist ---------- */
void createDir(const char* path, const char* fname) {
  struct stat st = {0};  // initializing all of its fields to zero
  // Make a copy of the path because the ``strtok`` function modifies the string
  char* path_copy = strdup(path);
  char* dir = strtok(path_copy, "/");
  char accum_path[256] = "";
  /* NOTE: The ``strtok`` function uses a static buffer to keep track of what
  it's doing. The first time you call ``strtok``, you pass it the string you
  want to tokenize. strtok saves that string in its static buffer and returns
  a pointer to the first token.
  When you call ``strtok`` with ``NULL`` as the first argument, strtok doesn't
  actually try to split ``NULL``. Instead, it looks at its static buffer to
  see if there's any leftover string from the last time it was called. If
  there is, it continues tokenizing that string.
  Because it uses a static buffer, it's not thread-safe, and you can only
  tokenize one string at a time. If you try to use ``strtok`` to tokenize two
  strings at the same time, you'll get unexpected results because the calls to
  ``strtok`` will interfere with each other. */

  // Handle absolute paths
  if (path[0] == '/') {
    strcat(accum_path, "/");
  }
  // Create the directory recursively if it doesn't exist
  while (dir != NULL) {
    strcat(accum_path, dir);
    if (stat(accum_path, &st) == -1) {
      /* The full permission set for directories is 0777 (read, write,
      execute for owner, group, and others). It would be subtracted by
      the umask value to get the actual permission set. */
      if (mkdir(accum_path, 0777) == -1) {
        solver_abort(fname, NULL, "Error creating directory %s\n", accum_path);
      }
    }
    strcat(accum_path, "/");
    dir = strtok(NULL, "/");
  }
  free(path_copy);
}

static void create_parent_dir_for_file(const char* path, const char* fname) {
  char dir[256];
  char* slash;

  if ((path == NULL) || (path[0] == '\0')) {
    return;
  }

  strncpy(dir, path, sizeof(dir));
  dir[sizeof(dir) - 1] = '\0';
  slash = strrchr(dir, '/');
  if (slash == NULL) {
    return;
  }

  if (slash == dir) {
    dir[1] = '\0';
  } else {
    *slash = '\0';
  }

  if (dir[0] != '\0') {
    createDir(dir, fname);
  }
}

/* ---------- Static global variables ------------------------------------ */

/* These are all of the input parameters - add new ones here */
static struct Param_t {
  char FourDOutFile[256];
  FILE* FourDOutFp;
  FILE* theMonitorFileFp;
  char* theMonitorFileName;
  char parameters_input_file[256];
  char cvmdb_input_file[256];
  char mesh_etree_output_file[256];
  char planes_input_file[256];
  char basin_input_file[256];
  char thePlaneDirOut[256];
  char theSourceOutputDir[256];
  char theMeshDirMatlab[256];
  char the3DVelModelDir[256];
  char theTopoDir[256];
  char theStationMediumsFile[256];
  char theAcousticInterfaceDiagnosticsFile[256];
  double theVsCut;
  double theFactor;
  double theAcousticMeshRefinementFactor;
  double theFreq;
  double theFreq_Vel;
  double theDeltaT;
  double theDeltaTSquared;
  double theEndT;
  double theStartT;
  double theDomainAzimuth;
  int monitor_stats_rate;
  double theSofteningFactor;
  int theStepMeshingFactor;
  int theAcousticInterfaceMaxLevel;
  int32_t theTotalSteps;
  int32_t theRate;
  damping_type_t theTypeOfDamping;
  double theThresholdDamping;
  double theThresholdVpVs;
  int theDampingStatisticsFlag;
  int theSchedulePrintErrorCheckFlag;
  int theSchedulePrintToStdout;
  int theSchedulePrintToFile;
  char* theSchedulePrintFilename;
  char* theScheduleStatFilename;
  char* theMeshStatFilename;
  char* theTopoStatFilename;
  noyesflag_t printStationVelocities;
  noyesflag_t printK;
  noyesflag_t printStationAccelerations;
  noyesflag_t printStationPressures;
  noyesflag_t includeBuildings;
  noyesflag_t includeTopography;
  noyesflag_t useParametricQ;
  noyesflag_t includeNonlinearAnalysis;
  noyesflag_t useInfQk;
  noyesflag_t includeIncidentPlaneWaves;
  noyesflag_t includeHomogeneousHalfSpace;
  noyesflag_t IstanbulVelocityModel;
  noyesflag_t basinVelocityModel;
  noyesflag_t threeDVelocityModel;
  noyesflag_t includeAcoustic;
  noyesflag_t acousticInterfaceCoupling;
  noyesflag_t acousticPressureFormulation;
  noyesflag_t reconstructAcousticPlaneDisplacements;

  int theTimingBarriersFlag;
  stiffness_type_t theStiffness;
  int theStationsPrintRate;
  double* theStationX;
  double* theStationY;
  double* theStationZ;
  int32_t theNumberOfStations;
  int32_t myNumberOfStations;
  int IO_pool_pe_count;
  int32_t thePlanePrintRate;
  int theNumberOfPlanes;
  char theStationsDirOut[256];
  station_t* myStations;
  int theCheckPointingRate;
  int theUseCheckPoint;
  char theCheckPointingDirOut[256];
  noyesflag_t storeMeshCoordinatesForMatlab;
  int theMeshOnlyMode;
  int theMeshInputMode;
  double the4DOutSize;
  int theMeshOutFlag;
  char theCVMFlatFile[128];
  output_parameters_t theOutputParameters;
  double theRegionLong;
  double theRegionLat;
  double theRegionDepth;
  double region_depth_deep_m;
  double theSurfaceCornersLong[4];
  double theSurfaceCornersLat[4];
  double theDomainX;
  double theDomainY;
  double theDomainZ;
  noyesflag_t drmImplement;
  drm_part_t theDrmPart;
  noyesflag_t useProfile;
  int32_t theNumberOfLayers;
  double* theProfileZ;
  double* theProfileVp;
  double* theProfileVs;
  double* theProfileRho;
  double* theProfileQp;
  double* theProfileQs;
  double theQConstant;
  double theQAlpha;
  double theQBeta;
  double theBasinLong;
  double theBasinLat;
  int basinXCount, basinYCount, basinZCount;
  int basinMaterialType;
  double basinIncrement;
  double basinWaterVp;
  double basinWaterRho;
  double** basinData;
  double the3DVelocityModelLong;
  double the3DVelocityModelLat;
  UTMZone_t utmZone;

  /* Absorbing boundary configuration. Placed at the end of the struct so the
   * designated initializer below can keep these fields in declaration order
   * (NVCC compiles psolve.cu as C++, which requires that). */
  abc_type_t theAbcType;
  double theCpmlThicknessM;
  int theCpmlLayers;
  double theCpmlReflectionTarget;
  int theCpmlPolyOrder;
  double theCpmlKappaMax;
  double theCpmlAlphaMaxHz;
  double theCpmlCflSafety;
  double theCpmlMpmlCrossFraction;
  double theCpmlAlphaFloorFraction;
  int theCpmlPmlMaterialDamping;        /* 0 none(default) 1 bkt 2 rayleigh */
  noyesflag_t theCpmlClampOuter;
  /* YES = native effective stiffness + CPML correction (default);
   * NO  = standalone 1-point reduced CPML stiffness. Declared/initialised
   * AFTER theCpmlClampOuter to preserve C++ designated-init order. */
  noyesflag_t theCpmlNativeStiffness;
  /* Acoustic PML stiffness mode. Defaults to theCpmlNativeStiffness
   * (so existing configs are bit-unchanged); overridable independently
   * via cpml_stiffness_mode_acoustic. Declared AFTER theCpmlNativeStiffness
   * to preserve C++ designated-init order. */
  noyesflag_t theCpmlNativeStiffnessAcoustic;
  /* Acoustic-elastic coupling formulation + interface source filter.
   * DEFAULT is velocity_potential (1) — the cross-validated stable
   * coupling; pressure (0) is the legacy explicit opt-in. Only affects
   * acoustic / acoustic-elastic runs; pure-elastic is unaffected. The
   * interface source filter stays gated off (found destabilizing). */
  int theAcousticVelocityPotential;          /* 0 pressure, 1 vel-pot (DEFAULT) */
  noyesflag_t theAcousticInterfaceFilter;    /* EMA low-pass on/off */
  double theAcousticInterfaceFilterAlpha;    /* EMA coeff (0,1]      */
  double theAcousticInterfaceFilterStrength; /* blend [0,1]          */
  double theAcousticQualityFactor;           /* water Q; 0=lossless  */
} Param = {.FourDOutFp = NULL,
           .theMonitorFileFp = NULL,
           .theMonitorFileName = NULL,
           .thePlaneDirOut = "outputfiles/planes",
           .theSourceOutputDir = "outputfiles/srctmp",
           .theMeshDirMatlab = "outputfiles/For_Matlab",
           .the3DVelModelDir = "inputfiles/materialfiles",
           .theTopoDir = "inputfiles/topography",
           .theStationMediumsFile = "",
           .theAcousticInterfaceDiagnosticsFile = "",
           .theAcousticMeshRefinementFactor = 1.0,
           .theFreq_Vel = 0,
           .monitor_stats_rate = 50,
           .theSofteningFactor = 0,
           .theStepMeshingFactor = 0,
           .theAcousticInterfaceMaxLevel = -1,
           .theSchedulePrintErrorCheckFlag = 0,
           .theSchedulePrintToStdout = 0,
           .theSchedulePrintToFile = 0,
           .theSchedulePrintFilename = "schedule_info.txt",
           .theScheduleStatFilename = NULL,
           .theMeshStatFilename = NULL,
           .theTopoStatFilename = NULL,
           .printStationPressures = NO,
           .includeAcoustic = NO,
           .acousticInterfaceCoupling = YES,
           .acousticPressureFormulation = NO,
           .reconstructAcousticPlaneDisplacements = NO,
           .theTimingBarriersFlag = 0,
           .myNumberOfStations = 0,
           .theStationsDirOut = "outputfiles/stations",
           .theUseCheckPoint = 0,
           .theMeshOnlyMode = 0,
           .theMeshInputMode = 0,
           .the4DOutSize = 0,
           .theMeshOutFlag = DO_OUTPUT,
           .useProfile = NO,
           .theNumberOfLayers = 0,
           .basinMaterialType = BASIN_MATERIAL_ELASTIC,
           .basinWaterVp = 1500.0,
           .basinWaterRho = 1000.0,
           .theAbcType = ABC_DASHPOT,
           .theCpmlThicknessM = 0.0,
           .theCpmlLayers = 10,
           .theCpmlReflectionTarget = 1.0e-3,
           .theCpmlPolyOrder = 2,
           .theCpmlKappaMax = 1.0,
           .theCpmlAlphaMaxHz = -1.0,
           .theCpmlCflSafety = 0.9,
           .theCpmlMpmlCrossFraction = 0.1,
           .theCpmlAlphaFloorFraction = 0.1,
           .theCpmlPmlMaterialDamping = 0,
           .theCpmlClampOuter = NO,
           .theCpmlNativeStiffness = YES,
           .theCpmlNativeStiffnessAcoustic = YES,
           .theAcousticVelocityPotential = 1,
           .theAcousticInterfaceFilter = NO,
           .theAcousticInterfaceFilterAlpha = 0.15,
           .theAcousticInterfaceFilterStrength = 1.0,
           .theAcousticQualityFactor = 0.0};

/* These are all of the remaining global variables - this list should not grow
 */
static struct Global_t {
  int32_t myID;
  int32_t theGroupSize;
  octree_t* myOctree;
  mesh_t* myMesh;
  int64_t theETotal;
  int64_t theNTotal;
  mysolver_t* mySolver;
  fvector_t* myVelocityTable;
  fmatrix_t theK1[8][8];
  fmatrix_t theK2[8][8];
  fmatrix_t theK3[8][8];
  double theQTABLE[26][6];
  double theABase;
  double theBBase;
  double theCriticalT;
  double fastestTimeSteps;
  double slowestTimeSteps;
  numerics_info_t theNumericsInformation;
  mpi_info_t theMPIInformation;
  int32_t theNodesLoaded;
  int32_t* theNodesLoadedList;
  vector3D_t* myForces;
  FILE* fpsource;
  etree_t* theCVMEp;
  etree_t* theMeshInputEp;
  int32_t theCVMQueryStage;
  double theXForMeshOrigin;
  double theYForMeshOrigin;
  double theZForMeshOrigin;
  cvmrecord_t* theCVMRecord;
  int theCVMRecordSize;
  int theCVMRecordCount;
  gpu_spec_t gpu_spec;
  gpu_kernel_t gpu_kernel[MAX_CUDA_KERNEL];
  gpu_data_t gpuData;
  uint8_t* acousticFaceFlags;
  int32_t* acousticInterfaceElementIds;
  uint8_t* acousticInterfaceFaceIds;
  double* acousticInterfaceFaceAreas;
  double* acousticInterfaceFaceCenters;
  double* acousticInterfaceFaceNormals;
  int32_t* acousticInterfaceAdjacentElementIndexes;
  double* acousticInterfaceAdjacentCenters;
  int64_t* acousticInterfaceAdjacentElementIds;
  int64_t acousticInterfaceElements;
  int64_t acousticInterfaceFaces;

  /* Host scratch for the CPML-aware interface coupling correction. Filled
   * in solver_init while g_cpmlBuild.face_mask is still alive; uploaded
   * alongside the other interface arrays. NULL when not applicable. */
  uint8_t* acousticInterfacePmlFlag;
  double* acousticInterfacePml_a;
  double* acousticInterfacePml_b;
  double* acousticInterfacePml_invKappa;
  int64_t acousticInterfacePmlFaceCount;

} Global = {.myID = -1,
            .theGroupSize = -1,
            .myOctree = NULL,
            .myMesh = NULL,
            .theETotal = 0,
            .theNTotal = 0,
            .mySolver = NULL,
            .myVelocityTable = NULL,
            .theCriticalT = 0,
            .fastestTimeSteps = 10000000,
            .slowestTimeSteps = 0,
            .theNodesLoaded = -1,
            .theNodesLoadedList = NULL,
            .myForces = NULL,
            .fpsource = NULL,
            .theXForMeshOrigin = 0.0,
            .theYForMeshOrigin = 0.0,
            .theZForMeshOrigin = 0.0,
            .theCVMRecordSize = sizeof(cvmrecord_t)};

struct station_medium_record_t {
  int32_t stationId;
  int32_t ownerRank;
  int32_t ownerIsAcoustic;
  int32_t pressureDefined;
  int32_t elementLevel;
  int32_t reserved;
  int64_t elementId;
  double x;
  double y;
  double z;
  double csi;
  double eta;
  double dzeta;
  double h;
  double Vp;
  double Vs;
  double rho;
};

/* ------------------------------End of
 * declarations------------------------------ */

#define ACOUSTIC_INTERFACE_FRAGMENT_MAX_LEVEL_DELTA 5

static int material_values_are_acoustic(double Vp, double Vs, double rho) {
  return (Param.includeAcoustic == YES) && (Vp > 0.0) && (rho > 0.0) &&
         (Vs <= ACOUSTIC_VS_EPS);
}

static int element_is_acoustic(const edata_t* edata) {
  return (edata != NULL) &&
         material_values_are_acoustic(edata->Vp, edata->Vs, edata->rho);
}

static int acoustic_pressure_enabled(void) {
  return (Param.includeAcoustic == YES) &&
         (Param.acousticPressureFormulation == YES);
}

static void acoustic_zero_damping(edata_t* edata) {
  if (edata == NULL) {
    return;
  }

  edata->a0_shear = 0.0f;
  edata->a1_shear = 0.0f;
  edata->a2_shear = 0.0f;
  edata->a0_kappa = 0.0f;
  edata->a1_kappa = 0.0f;
  edata->a2_kappa = 0.0f;
  edata->g0_shear = 0.0f;
  edata->g1_shear = 0.0f;
  edata->g2_shear = 0.0f;
  edata->g0_kappa = 0.0f;
  edata->g1_kappa = 0.0f;
  edata->g2_kappa = 0.0f;
  edata->b_shear = 0.0f;
  edata->b_kappa = 0.0f;
}

static void acoustic_sanitize_element(edata_t* edata) {
  if (!element_is_acoustic(edata)) {
    return;
  }

  edata->Vs = 0.0f;
  acoustic_zero_damping(edata);
}

static int acoustic_has_damping(const edata_t* edata) {
  return (edata != NULL) &&
         ((edata->a0_shear != 0.0f) || (edata->a1_shear != 0.0f) ||
          (edata->a2_shear != 0.0f) || (edata->a0_kappa != 0.0f) ||
          (edata->a1_kappa != 0.0f) || (edata->a2_kappa != 0.0f) ||
          (edata->g0_shear != 0.0f) || (edata->g1_shear != 0.0f) ||
          (edata->g2_shear != 0.0f) || (edata->g0_kappa != 0.0f) ||
          (edata->g1_kappa != 0.0f) || (edata->g2_kappa != 0.0f) ||
          (edata->b_shear != 0.0f) || (edata->b_kappa != 0.0f));
}

static void mesh_sanitize_and_report_acoustic(const char* context) {
  int32_t eindex;
  int64_t local_counts[2] = {0, 0};
  int64_t global_counts[2] = {0, 0};

  if ((Param.includeAcoustic != YES) || (Global.myMesh == NULL)) {
    return;
  }

  for (eindex = 0; eindex < Global.myMesh->lenum; eindex++) {
    edata_t* edata = (edata_t*)Global.myMesh->elemTable[eindex].data;

    if (element_is_acoustic(edata)) {
      local_counts[0]++;
      if (acoustic_has_damping(edata)) {
        local_counts[1]++;
      }
      acoustic_sanitize_element(edata);
    }
  }

  MPI_Reduce(local_counts, global_counts, 2, MPI_INT64_T, MPI_SUM, 0,
             comm_solver);
  if (Global.myID == 0) {
    monitor_print("%s acoustic elements: %.1fM, sanitized damping: %.1fM\n",
                  context, (double)global_counts[0] / 1.0e6,
                  (double)global_counts[1] / 1.0e6);
    fflush(stdout);
  }
}

static double element_meshing_velocity(const edata_t* edata) {
  if (element_is_acoustic(edata)) {
    return edata->Vp;
  }

  return (edata != NULL) ? edata->Vs : 0.0;
}

static int32_t wavespeed_rule(const edata_t* edata, double theFactor) {
  double speed = element_meshing_velocity(edata);
  double effective_factor = theFactor;

  if ((edata == NULL) || (edata->edgesize <= 0.0) || (speed <= 0.0)) {
    return 1;
  }

  if (element_is_acoustic(edata)) {
    effective_factor *= Param.theAcousticMeshRefinementFactor;
  }

  return (edata->edgesize <= speed / effective_factor) ? 0 : 1;
}

void getGPUHardware(int32_t myID, int device, gpu_spec_t* gpuSpecs) {
  int i;
  const int kb = 1024;
  const int mb = kb * kb;
  const long long gb = (long long)mb * mb;

  cudaDeviceProp props;
  memset(&props, 0, sizeof(cudaDeviceProp));
  if (cudaGetDeviceProperties(&props, device) != cudaSuccess) {
    fprintf(stderr, "Failed to get GPU device properties\n");
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  gpuSpecs->device = device;
  gpuSpecs->max_threads = props.maxThreadsPerBlock;
  for (i = 0; i < 3; i++) {
    gpuSpecs->max_block_dim[i] = props.maxThreadsDim[i];
    gpuSpecs->max_grid_dim[i] = props.maxGridSize[i];
  }
  gpuSpecs->regs_per_block = props.regsPerBlock;
  gpuSpecs->shared_per_block = props.sharedMemPerBlock;
  gpuSpecs->warp_size = props.warpSize;

  /* TODO: Customize these parameters for compute capability */
  gpuSpecs->warp_allocation_size = 4; // original value: 4
  gpuSpecs->register_allocation_size = 256; // original value: 256

  if (myID == 0) {
    printf("\nGPU Device %d: %s: %d.%d\n", device, props.name, props.major, props.minor);
    printf("  Global memory (HBM 3):   %lu mb (%lu gb)\n",
           props.totalGlobalMem / mb, (long long)props.totalGlobalMem / gb);
    printf("  Shared memory:           %lu kb\n", props.sharedMemPerBlock / kb);
    printf("  Constant memory:         %lu kb\n", props.totalConstMem / kb);
    printf("  Block registers: %d\n\n", props.regsPerBlock);

    printf("  Warp size:         %d\n", props.warpSize);
    printf("  Threads per block: %d\n", props.maxThreadsPerBlock);
    printf("  Max block dimensions: [ %d, %d, %d ]\n",
           props.maxThreadsDim[0], props.maxThreadsDim[1], props.maxThreadsDim[2]);
    printf("  Max grid dimensions:  [ %d, %d, %d ]\n",
           props.maxGridSize[0], props.maxGridSize[1], props.maxGridSize[2]);
    printf("  Registers per block: %d\n\n", props.regsPerBlock);
  }
}

static inline int monitor_print(const char* format, ...) {
  int ret = 0;

  if (format != NULL) {
    va_list args;

    va_start(args, format);

    if (Param.theMonitorFileFp == NULL) {
      ret = vfprintf(stderr, format, args);
    } else {
      ret = hu_print_tee_va(Param.theMonitorFileFp, format, args);
    }

    va_end(args);
  }

  return ret;
}

/*-----------Parameter input routines---------------------------------*/

static void read_parameters(int argc, char** argv) {
#define LOCAL_INIT_DOUBLE_MESSAGE_LENGTH \
  39 /* Must adjust this if adding double params */
#define LOCAL_INIT_INT_MESSAGE_LENGTH \
  49 /* Must adjust this if adding int params */

  double double_message[LOCAL_INIT_DOUBLE_MESSAGE_LENGTH];
  int int_message[LOCAL_INIT_INT_MESSAGE_LENGTH];

  strcpy(Param.parameters_input_file, argv[1]);

  /* PE 0 reads all params from disk */
  if (Global.myID == 0) {
    if (parse_parameters(Param.parameters_input_file) != 0) {
      fprintf(stderr, "Thread 0: Problem reading parameters!\n");
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
  }

  /*Broadcast all double params*/
  double_message[0] = Param.theVsCut;
  double_message[1] = Param.theFactor;
  double_message[2] = Param.theFreq;
  double_message[3] = Param.theDeltaT;
  double_message[4] = Param.theDeltaTSquared;
  double_message[5] = Param.theEndT;
  double_message[6] = Param.theStartT;
  double_message[7] = Param.theDomainX;
  double_message[8] = Param.theDomainY;
  double_message[9] = Param.theDomainZ;
  double_message[10] = Param.theDomainAzimuth;
  double_message[11] = Param.theThresholdDamping;
  double_message[12] = Param.theThresholdVpVs;
  double_message[13] = Param.theSofteningFactor;
  double_message[14] = Param.theFreq_Vel;
  double_message[15] = Param.theRegionLat;
  double_message[16] = Param.theRegionLong;
  double_message[17] = Param.theRegionDepth;
  double_message[18] = Param.theQConstant;
  double_message[19] = Param.theQAlpha;
  double_message[20] = Param.theQBeta;
  double_message[21] = Param.theBasinLong;
  double_message[22] = Param.theBasinLat;
  double_message[23] = Param.basinIncrement;
  double_message[24] = Param.the3DVelocityModelLong;
  double_message[25] = Param.the3DVelocityModelLat;
  double_message[26] = Param.basinWaterVp;
  double_message[27] = Param.basinWaterRho;
  double_message[28] = Param.theAcousticMeshRefinementFactor;
  double_message[29] = Param.theCpmlThicknessM;
  double_message[30] = Param.theCpmlReflectionTarget;
  double_message[31] = Param.theCpmlKappaMax;
  double_message[32] = Param.theCpmlAlphaMaxHz;
  double_message[33] = Param.theCpmlCflSafety;
  double_message[34] = Param.theCpmlMpmlCrossFraction;
  double_message[35] = Param.theCpmlAlphaFloorFraction;
  double_message[36] = Param.theAcousticInterfaceFilterAlpha;
  double_message[37] = Param.theAcousticInterfaceFilterStrength;
  double_message[38] = Param.theAcousticQualityFactor;

  MPI_Bcast(double_message, LOCAL_INIT_DOUBLE_MESSAGE_LENGTH, MPI_DOUBLE, 0,
            comm_solver);

  Param.theVsCut = double_message[0];
  Param.theFactor = double_message[1];
  Param.theFreq = double_message[2];
  Param.theDeltaT = double_message[3];
  Param.theDeltaTSquared = double_message[4];
  Param.theEndT = double_message[5];
  Param.theStartT = double_message[6];
  Param.theDomainX = double_message[7];
  Param.theDomainY = double_message[8];
  Param.theDomainZ = double_message[9];
  Param.theDomainAzimuth = double_message[10];
  Param.theThresholdDamping = double_message[11];
  Param.theThresholdVpVs = double_message[12];
  Param.theSofteningFactor = double_message[13];
  Param.theFreq_Vel = double_message[14];
  Param.theRegionLat = double_message[15];
  Param.theRegionLong = double_message[16];
  Param.theRegionDepth = double_message[17];
  Param.theQConstant = double_message[18];
  Param.theQAlpha = double_message[19];
  Param.theQBeta = double_message[20];
  Param.theBasinLong = double_message[21];
  Param.theBasinLat = double_message[22];
  Param.basinIncrement = double_message[23];
  Param.the3DVelocityModelLong = double_message[24];
  Param.the3DVelocityModelLat = double_message[25];
  Param.basinWaterVp = double_message[26];
  Param.basinWaterRho = double_message[27];
  Param.theAcousticMeshRefinementFactor = double_message[28];
  Param.theCpmlThicknessM = double_message[29];
  Param.theCpmlReflectionTarget = double_message[30];
  Param.theCpmlKappaMax = double_message[31];
  Param.theCpmlAlphaMaxHz = double_message[32];
  Param.theCpmlCflSafety = double_message[33];
  Param.theCpmlMpmlCrossFraction = double_message[34];
  Param.theCpmlAlphaFloorFraction = double_message[35];
  Param.theAcousticInterfaceFilterAlpha = double_message[36];
  Param.theAcousticInterfaceFilterStrength = double_message[37];
  Param.theAcousticQualityFactor = double_message[38];

  /*Broadcast all integer params*/
  int_message[0] = Param.theTotalSteps;
  int_message[1] = Param.theRate;
  int_message[2] = Param.theNumberOfPlanes;
  int_message[3] = Param.theNumberOfStations;
  int_message[4] = (int)Param.theTypeOfDamping;
  int_message[5] = Param.theDampingStatisticsFlag;
  int_message[6] = Param.theMeshOutFlag;
  int_message[7] = Param.theCheckPointingRate;
  int_message[8] = Param.theUseCheckPoint;
  int_message[9] = (int)Param.includeNonlinearAnalysis;
  int_message[10] = (int)Param.theStiffness;
  int_message[11] = (int)Param.printK;
  int_message[12] = (int)Param.printStationVelocities;
  int_message[13] = (int)Param.printStationAccelerations;
  int_message[14] = Param.theTimingBarriersFlag;
  int_message[15] = (int)Param.includeBuildings;
  int_message[16] = (int)Param.storeMeshCoordinatesForMatlab;
  int_message[17] = (int)Param.drmImplement;
  int_message[18] = (int)Param.useInfQk;
  int_message[19] = Param.theStepMeshingFactor;
  int_message[20] = (int)Param.includeTopography;
  int_message[21] = (int)Param.includeIncidentPlaneWaves;
  int_message[22] = (int)Param.includeHomogeneousHalfSpace;
  int_message[23] = (int)Param.IstanbulVelocityModel;
  int_message[24] = (int)Param.useProfile;
  int_message[25] = (int)Param.useParametricQ;
  int_message[26] = (int)Param.basinVelocityModel;
  int_message[27] = (int)Param.basinXCount;
  int_message[28] = (int)Param.basinYCount;
  int_message[29] = (int)Param.basinZCount;
  int_message[30] = (int)Param.threeDVelocityModel;
  int_message[31] = Param.theMeshOnlyMode;
  int_message[32] = Param.theMeshInputMode;
  int_message[33] = (int)Param.includeAcoustic;
  int_message[34] = Param.theAcousticInterfaceMaxLevel;
  int_message[35] = (int)Param.acousticInterfaceCoupling;
  int_message[36] = (int)Param.acousticPressureFormulation;
  int_message[37] = (int)Param.printStationPressures;
  int_message[38] = (int)Param.reconstructAcousticPlaneDisplacements;
  int_message[39] = Param.basinMaterialType;
  int_message[40] = (int)Param.theAbcType;
  int_message[41] = Param.theCpmlLayers;
  int_message[42] = Param.theCpmlPolyOrder;
  int_message[43] = (int)Param.theCpmlClampOuter;
  int_message[44] = (int)Param.theCpmlNativeStiffness;
  int_message[45] = (int)Param.theCpmlNativeStiffnessAcoustic;
  int_message[46] = Param.theAcousticVelocityPotential;
  int_message[47] = (int)Param.theAcousticInterfaceFilter;
  int_message[48] = Param.theCpmlPmlMaterialDamping;

  MPI_Bcast(int_message, LOCAL_INIT_INT_MESSAGE_LENGTH, MPI_INT, 0,
            comm_solver);

  Param.theTotalSteps = int_message[0];
  Param.theRate = int_message[1];
  Param.theNumberOfPlanes = int_message[2];
  Param.theNumberOfStations = int_message[3];
  Param.theTypeOfDamping = (damping_type_t)int_message[4];
  Param.theDampingStatisticsFlag = int_message[5];
  Param.theMeshOutFlag = int_message[6];
  Param.theCheckPointingRate = int_message[7];
  Param.theUseCheckPoint = int_message[8];
  Param.includeNonlinearAnalysis = (noyesflag_t)int_message[9];
  Param.theStiffness = (stiffness_type_t)int_message[10];
  Param.printK = (noyesflag_t)int_message[11];
  Param.printStationVelocities = (noyesflag_t)int_message[12];
  Param.printStationAccelerations = (noyesflag_t)int_message[13];
  Param.theTimingBarriersFlag = int_message[14];
  Param.includeBuildings = (noyesflag_t)int_message[15];
  Param.storeMeshCoordinatesForMatlab = (noyesflag_t)int_message[16];
  Param.drmImplement = (noyesflag_t)int_message[17];
  Param.useInfQk = (noyesflag_t)int_message[18];
  Param.theStepMeshingFactor = int_message[19];
  Param.includeTopography = (noyesflag_t)int_message[20];
  Param.includeIncidentPlaneWaves = (noyesflag_t)int_message[21];
  Param.includeHomogeneousHalfSpace = (noyesflag_t)int_message[22];
  Param.IstanbulVelocityModel = (noyesflag_t)int_message[23];
  Param.useProfile = (noyesflag_t)int_message[24];
  Param.useParametricQ = (noyesflag_t)int_message[25];
  Param.basinVelocityModel = (noyesflag_t)int_message[26];
  Param.basinXCount = int_message[27];
  Param.basinYCount = int_message[28];
  Param.basinZCount = int_message[29];
  Param.threeDVelocityModel = (noyesflag_t)int_message[30];
  Param.theMeshOnlyMode = int_message[31];
  Param.theMeshInputMode = int_message[32];
  Param.includeAcoustic = (noyesflag_t)int_message[33];
  Param.theAcousticInterfaceMaxLevel = int_message[34];
  Param.acousticInterfaceCoupling = (noyesflag_t)int_message[35];
  Param.acousticPressureFormulation = (noyesflag_t)int_message[36];
  Param.printStationPressures = (noyesflag_t)int_message[37];
  Param.reconstructAcousticPlaneDisplacements =
      (noyesflag_t)int_message[38];
  Param.basinMaterialType = int_message[39];
  Param.theAbcType = (abc_type_t)int_message[40];
  Param.theCpmlLayers = int_message[41];
  Param.theCpmlPolyOrder = int_message[42];
  Param.theCpmlClampOuter = (noyesflag_t)int_message[43];
  Param.theCpmlNativeStiffness = (noyesflag_t)int_message[44];
  Param.theCpmlNativeStiffnessAcoustic = (noyesflag_t)int_message[45];
  Param.theAcousticVelocityPotential = int_message[46];
  Param.theAcousticInterfaceFilter = (noyesflag_t)int_message[47];
  Param.theCpmlPmlMaterialDamping = int_message[48];

  /*Broadcast all string params*/
  MPI_Bcast(Param.parameters_input_file, 256, MPI_CHAR, 0, comm_solver);
  MPI_Bcast(Param.theCheckPointingDirOut, 256, MPI_CHAR, 0, comm_solver);
  MPI_Bcast(Param.FourDOutFile, 256, MPI_CHAR, 0, comm_solver);
  MPI_Bcast(Param.cvmdb_input_file, 256, MPI_CHAR, 0, comm_solver);
  MPI_Bcast(Param.mesh_etree_output_file, 256, MPI_CHAR, 0, comm_solver);
  MPI_Bcast(Param.planes_input_file, 256, MPI_CHAR, 0, comm_solver);
  MPI_Bcast(Param.basin_input_file, 256, MPI_CHAR, 0, comm_solver);
  MPI_Bcast(Param.thePlaneDirOut, 256, MPI_CHAR, 0, comm_solver);
  MPI_Bcast(Param.theSourceOutputDir, 256, MPI_CHAR, 0, comm_solver);
  MPI_Bcast(Param.theMeshDirMatlab, 256, MPI_CHAR, 0, comm_solver);
  MPI_Bcast(Param.the3DVelModelDir, 256, MPI_CHAR, 0, comm_solver);
  MPI_Bcast(Param.theTopoDir, 256, MPI_CHAR, 0, comm_solver);
  MPI_Bcast(Param.theAcousticInterfaceDiagnosticsFile, 256, MPI_CHAR, 0,
            comm_solver);

  /*Broadcast domain's coords */
  MPI_Bcast(Param.theSurfaceCornersLong, 4, MPI_DOUBLE, 0, comm_solver);
  MPI_Bcast(Param.theSurfaceCornersLat, 4, MPI_DOUBLE, 0, comm_solver);

  /* Broadcast basin data */
  if (Param.basinVelocityModel == YES) {
    int total_lines = Param.basinXCount * Param.basinYCount * Param.basinZCount;
    if (Global.myID != 0) {
      Param.basinData = malloc2dDouble(total_lines, 3);
    }
    MPI_Bcast(&Param.basinData[0][0], 3 * total_lines, MPI_DOUBLE, 0,
              comm_solver);
  }

#ifdef PROJ
  /* Compute UTM Zone */
  Param.utmZone = getUTMZone(Param.theRegionLat, Param.theRegionLong);
  Param.utmZone.refPoint = convertLatLonToUTM(
      Param.theRegionLat, Param.theRegionLong, Param.utmZone.P);
  /* Compute the corners of the surface */
  /* NOTE: This is not necessary since the latitude and longitude of the
  corners (except the origin) are not used when the PROJ library is used.
  But it can be printed for post-processing purposes. */
  UTMCoordinates_t cornersUTM[4];
  cornersUTM[0].n = Param.utmZone.refPoint.n;
  cornersUTM[0].e = Param.utmZone.refPoint.e;
  cornersUTM[1].n = cornersUTM[0].n + Param.theDomainX;
  cornersUTM[1].e = cornersUTM[0].e;
  cornersUTM[2].n = cornersUTM[0].n + Param.theDomainX;
  cornersUTM[2].e = cornersUTM[0].e + Param.theDomainY;
  cornersUTM[3].n = cornersUTM[0].n;
  cornersUTM[3].e = cornersUTM[0].e + Param.theDomainY;
  double lat, lon;
  for (int iCorner = 0; iCorner < 4; iCorner++) {
    convertUTMToLatLon(cornersUTM[iCorner], Param.utmZone.P, &lat, &lon);
    Param.theSurfaceCornersLong[iCorner] = lon;
    Param.theSurfaceCornersLat[iCorner] = lat;
  }
  /* Print the UTM Zone and the corners of the surface */
  if (Global.myID == 0) {
    printf("UTM Zone: %d%c\n", Param.utmZone.zone, Param.utmZone.hemisphere);
    printf("UTM Reference Point (X-Northing, Y-Easting): %f, %f\n",
           Param.utmZone.refPoint.n, Param.utmZone.refPoint.e);
    printf("UTM Corners:\n");
    for (int iCorner = 0; iCorner < 4; iCorner++) {
      printf(" %f, %f\n", cornersUTM[iCorner].n, cornersUTM[iCorner].e);
    }
    printf("Longitude and Latitude of the Corners:\n");
    for (int iCorner = 0; iCorner < 4; iCorner++) {
      printf(" %f, %f\n", Param.theSurfaceCornersLong[iCorner],
             Param.theSurfaceCornersLat[iCorner]);
    }
    printf("\n");
  }
#endif

  return;
}

static void local_finalize() {
  /* Free memory associated with the octree and mesh */
  octor_deletetree(Global.myOctree);
  octor_deletemesh(Global.myMesh);
  mesh_release_acoustic_interface_data();

  /* Free the memory used by the source */
  free(Global.myForces);

  if (Global.myVelocityTable != NULL) free(Global.myVelocityTable);

  /* Free memory associated with the stiffness module */
  stiffness_delete(Global.myID);

  /* Free memory associated with the solver */
  solver_delete();

  return;
}

/**
 * Parse a text file and return the value of a match string.
 *
 * \return 0 if OK, -1 on error.
 */
int parsetext(FILE* fp, const char* querystring, const char type,
              void* result) {
  const static char delimiters[] = " =\n\t";

  int32_t res = 0, found = 0;

  /* Start from the beginning */
  rewind(fp);

  /* Look for the string until found */
  while (!found) {
    char line[LINESIZE];
    char *name, *value;

    /* Read in one line */
    if (fgets(line, LINESIZE, fp) == NULL) break;

    name = strtok(line, delimiters);
    if ((name != NULL) && (strcmp(name, querystring) == 0)) {
      found = 1;
      value = strtok(NULL, delimiters);

      switch (type) {
        case 'i':
          res = sscanf(value, "%d", (int*)result);
          break;
        case 'f':
          res = sscanf(value, "%f", (float*)result);
          break;
        case 'd':
          res = sscanf(value, "%lf", (double*)result);
          break;
        case 's':
          res = 1;
          strcpy((char*)result, value);
          break;
        case 'u':
          res = sscanf(value, "%u", (uint32_t*)result);
          break;
        default:
          fprintf(stderr, "parsetext: unknown type %c\n", type);
          return -1;
      }
    }
  }

  return (res == 1) ? 0 : -1;
}

/**
 * This is like \c parsetext with the following differences:
 * - works only for strings;
 * - avoids writting past the end of the string;
 * - return convention is different, it distinguishes between "key not found"
 *   and other type of errors.
 *
 * \return
 *  1 if the key name is found and the value is stored in \c value;
 *  0 if the key name was not found; -1 on error.
 */
static int read_config_string(FILE* fp, const char* key, char* value_ptr,
                              size_t size) {
  static const char delimiters[] = " =\n\t";

  int ret;
  char line[LINESIZE];
  char state[LINESIZE];
  char *name, *value, *state_ptr;

  HU_ASSERT_PTR_ALWAYS(fp);
  HU_ASSERT_PTR_ALWAYS(value_ptr);
  HU_ASSERT_PTR_ALWAYS(key);

  rewind(fp);
  ret = 0;
  *value_ptr = '\0';

  while (0 == ret && !ferror(fp)) {
    if (fgets(line, LINESIZE, fp) == NULL) {
      if (!feof(fp)) {
        ret = -1; /* input error */
      }
      break;
    }

    state_ptr = state;
    name = strtok_r(line, delimiters, &state_ptr);

    if ((name != NULL) && (strcmp(name, key) == 0)) {
      size_t value_len;

      value = strtok_r(NULL, delimiters, &state_ptr);

      if (NULL != value) {
        value_len = strlen(value);

        if (value_len >= size) {
          ret = -2; /* return buffer is too short */
        } else {
          strncpy(value_ptr, value, size);
          ret = 1;
        }
      }

      break;
    }
  }

  return ret;
}

/**
 * Open material database and initialize various static global variables.
 *
 * \return 0 if OK, -1 on error.
 */
static int32_t parse_parameters(const char* numericalin) {
  FILE* fp;

  /*
   * This used to be a separate file, now aliased to numericalin.
   * This fakes a separate file as per old code
   */
  const char* physicsin = numericalin;

  int32_t samples;

  double freq, vscut, region_origin_latitude_deg, region_origin_longitude_deg,
      region_azimuth_leftface_deg, region_depth_shallow_m, region_length_east_m,
      region_length_north_m, region_depth_deep_m, endT, deltaT,
      threshold_damping, threshold_VpVs, qconstant, qalpha, qbeta;
  char type_of_damping[64], use_parametricq[64], use_infinite_qk[64],
      which_drm_part[64], drm_directory_full_path[256];
  /* Optional parameters */
  int32_t rate = 1000000;
  int number_output_planes = 0, number_output_stations = 0,
      damping_statistics = 0, use_checkpoint = 0, checkpointing_rate = 0,
      step_meshing = 1, acoustic_interface_max_level = -1;
  double freq_vel = 0., softening_factor = 0, startT = 0.,
         acoustic_mesh_refinement_factor = 1.0;
  char print_matrix_k[64] = "no", print_station_velocities[64] = "no",
       print_station_accelerations[64] = "no",
       print_station_pressures[64] = "no",
       mesh_coordinates_for_matlab[64] = "no", include_buildings[64] = "no",
       implement_drm[64] = "no", include_topography[64] = "no",
       include_acoustic[64] = "no", acoustic_interface_coupling[64] = "yes",
       acoustic_pressure_formulation[64] = "no",
       reconstruct_acoustic_plane_displacements[64] = "no",
       basin_material_type[64] = "elastic",
       IstanbulModel[64] = "no", basinModel[64] = "no", threeDModel[64] = "no",
       include_nonlinear_analysis[64] = "no",
       include_incident_planewaves[64] = "no", include_hmgHalfSpace[64] = "no",
       stiffness_calculation_method[64] = "effective",
       checkpoint_path[256] = "outputfiles/checkpoints",
       drm_directory[256] = "outputfiles/DRM";

  damping_type_t typeOfDamping = UNSET_DAMP;
  stiffness_type_t stiffness_method = UNSET_STIFF;
  noyesflag_t have_buildings = UNSET;
  noyesflag_t have_parametricq = UNSET;
  noyesflag_t includeNonlinear = UNSET;
  noyesflag_t printMatrixK = UNSET;
  noyesflag_t printStationVels = UNSET;
  noyesflag_t printStationAccs = UNSET;
  noyesflag_t printStationPressures = NO;
  noyesflag_t useInfQk = UNSET;

  noyesflag_t meshCoordinatesForMatlab = UNSET;
  noyesflag_t implementdrm = UNSET;
  noyesflag_t haveTopography = UNSET;
  noyesflag_t includeAcoustic = UNSET;
  noyesflag_t acousticInterfaceCoupling = YES;
  noyesflag_t acousticPressureFormulation = NO;
  noyesflag_t reconstructAcousticPlaneDisplacements = NO;
  noyesflag_t includePlaneWaves = UNSET;
  noyesflag_t includeHmgHalfSpace = UNSET;
  noyesflag_t includeIstanbulVelModel = UNSET;
  noyesflag_t includeBasinVelModel = UNSET;
  noyesflag_t include3DVelModel = UNSET;
  basin_material_type_t basinMaterialType = BASIN_MATERIAL_ELASTIC;

  /* Obtain the specification of the simulation */
  if ((fp = fopen(physicsin, "r")) == NULL) {
    fprintf(stderr, "Error opening %s\n", physicsin);
    return -1;
  }

  /* Julio, I know this violates the 80 chars printing limit, but we never
   * print this code and is a heck of a lot easier to read it like this
   * --Ricardo
   */
  if ((parsetext(fp, "region_origin_latitude_deg", 'd',
                 &region_origin_latitude_deg) != 0) ||
      (parsetext(fp, "region_origin_longitude_deg", 'd',
                 &region_origin_longitude_deg) != 0) ||
      (parsetext(fp, "region_depth_shallow_m", 'd', &region_depth_shallow_m) !=
       0) ||
      (parsetext(fp, "region_length_east_m", 'd', &region_length_east_m) !=
       0) ||
      (parsetext(fp, "region_length_north_m", 'd', &region_length_north_m) !=
       0) ||
      (parsetext(fp, "region_depth_deep_m", 'd', &region_depth_deep_m) != 0) ||
      (parsetext(fp, "region_azimuth_leftface_deg", 'd',
                 &region_azimuth_leftface_deg) != 0) ||
      (parsetext(fp, "type_of_damping", 's', &type_of_damping) != 0)) {
    fprintf(stderr, "Error reading region origin from %s\n", physicsin);
    return -1;
  }

  if (strcasecmp(type_of_damping, "rayleigh") == 0) {
    typeOfDamping = RAYLEIGH;
  } else if (strcasecmp(type_of_damping, "mass") == 0) {
    typeOfDamping = MASS;
  } else if (strcasecmp(type_of_damping, "none") == 0) {
    typeOfDamping = NONE;
  } else if (strcasecmp(type_of_damping, "bkt") == 0) {
    typeOfDamping = BKT;
  } else if (strcasecmp(type_of_damping, "bkt2") == 0) {
    typeOfDamping = BKT2;
  } else if (strcasecmp(type_of_damping, "bkt3") == 0) {
    typeOfDamping = BKT3;
  } else if (strcasecmp(type_of_damping, "bkt3f") == 0) {
    typeOfDamping = BKT3F;
  } else {
    solver_abort(__FUNCTION_NAME, NULL, "Unknown damping type: %s\n",
                 type_of_damping);
  }

  if (typeOfDamping >= BKT) {
    /* Conditional Parameters */
    if ((parsetext(fp, "use_parametric_q_factor", 's', &use_parametricq) !=
         0) ||
        (parsetext(fp, "use_infinite_qk", 's', &use_infinite_qk) != 0)) {
      fprintf(stderr, "Error reading q factor related parameters from %s\n",
              numericalin);
      return -1;
    }

    if (strcasecmp(use_parametricq, "yes") == 0) {
      have_parametricq = YES;
    } else if (strcasecmp(use_parametricq, "no") == 0) {
      have_parametricq = NO;
    } else {
      solver_abort(__FUNCTION_NAME, NULL,
                   "Unknown response for use new q factor (yes or no): %s\n",
                   use_parametricq);
    }

    if (strcasecmp(use_infinite_qk, "yes") == 0) {
      useInfQk = YES;
    } else if (strcasecmp(use_infinite_qk, "no") == 0) {
      useInfQk = NO;
    } else {
      solver_abort(__FUNCTION_NAME, NULL,
                   "Unknown response using infinite Qk (yes or no): %s\n",
                   use_infinite_qk);
    }

    if (have_parametricq == YES) {
      if ((parsetext(fp, "parametric_q_factor_constant", 'd', &qconstant) !=
           0) ||
          (parsetext(fp, "parametric_q_factor_alpha", 'd', &qalpha) != 0) ||
          (parsetext(fp, "parametric_q_factor_beta", 'd', &qbeta) != 0)) {
        fprintf(stderr, "Error reading parametric q factor from %s\n",
                numericalin);
        return -1;
      }
    }
  }

  /* Absorbing boundary condition selection (optional; defaults to dashpot). */
  {
    char type_of_abc[64] = "dashpot";
    abc_type_t abcType = ABC_DASHPOT;
    if (parsetext(fp, "abc_type", 's', &type_of_abc) == 0) {
      if (strcasecmp(type_of_abc, "dashpot") == 0) {
        abcType = ABC_DASHPOT;
      } else if (strcasecmp(type_of_abc, "cpml") == 0) {
        abcType = ABC_CPML;
      } else {
        solver_abort(__FUNCTION_NAME, NULL,
                     "Unknown abc_type: %s (expected dashpot or cpml)\n",
                     type_of_abc);
      }
    }
    Param.theAbcType = abcType;

    if (abcType == ABC_CPML) {
      /* All CPML parameters are optional; defaults set in Param initializer. */
      double cpml_thickness_m = Param.theCpmlThicknessM;
      int cpml_layers = Param.theCpmlLayers;
      double cpml_reflection_target = Param.theCpmlReflectionTarget;
      int cpml_poly_order = Param.theCpmlPolyOrder;
      double cpml_kappa_max = Param.theCpmlKappaMax;
      double cpml_alpha_max_hz = Param.theCpmlAlphaMaxHz;
      double cpml_cfl_safety = Param.theCpmlCflSafety;
      double cpml_mpml_cross_fraction = Param.theCpmlMpmlCrossFraction;
      double cpml_alpha_floor_fraction = Param.theCpmlAlphaFloorFraction;

      (void)parsetext(fp, "cpml_thickness_m", 'd', &cpml_thickness_m);
      (void)parsetext(fp, "cpml_layers", 'i', &cpml_layers);
      (void)parsetext(fp, "cpml_reflection_target", 'd',
                      &cpml_reflection_target);
      (void)parsetext(fp, "cpml_poly_order", 'i', &cpml_poly_order);
      (void)parsetext(fp, "cpml_kappa_max", 'd', &cpml_kappa_max);
      (void)parsetext(fp, "cpml_alpha_max_hz", 'd', &cpml_alpha_max_hz);
      (void)parsetext(fp, "cpml_cfl_safety", 'd', &cpml_cfl_safety);
      (void)parsetext(fp, "cpml_mpml_cross_fraction", 'd',
                      &cpml_mpml_cross_fraction);
      (void)parsetext(fp, "cpml_alpha_floor_fraction", 'd',
                      &cpml_alpha_floor_fraction);

      char cpml_clamp_outer[64] = "";
      if (parsetext(fp, "cpml_clamp_outer", 's', &cpml_clamp_outer) == 0) {
        if (strcasecmp(cpml_clamp_outer, "yes") == 0) {
          Param.theCpmlClampOuter = YES;
        } else if (strcasecmp(cpml_clamp_outer, "no") == 0) {
          Param.theCpmlClampOuter = NO;
        } else {
          solver_abort(__FUNCTION_NAME, NULL,
                       "cpml_clamp_outer must be 'yes' or 'no': got '%s'\n",
                       cpml_clamp_outer);
        }
      }

      /* Elastic PML stiffness mode (default native). 'native' = full
       * native effective stiffness + CPML correction (exact interior
       * operator at the inner edge, native hourglass control). 'reduced'
       * = standalone 1-point reduced CPML stiffness (no higher modes;
       * the pre-a1c0fb8 behaviour, kept for comparison). */
      char cpml_stiffness_mode[64] = "";
      if (parsetext(fp, "cpml_stiffness_mode", 's',
                    &cpml_stiffness_mode) == 0) {
        if (strcasecmp(cpml_stiffness_mode, "native") == 0) {
          Param.theCpmlNativeStiffness = YES;
        } else if (strcasecmp(cpml_stiffness_mode, "reduced") == 0) {
          Param.theCpmlNativeStiffness = NO;
        } else {
          solver_abort(__FUNCTION_NAME, NULL,
                       "cpml_stiffness_mode must be 'native' or "
                       "'reduced': got '%s'\n", cpml_stiffness_mode);
        }
      }

      /* Acoustic PML stiffness mode. Defaults to the (resolved) elastic
       * value so existing configs are bit-unchanged, but can be set
       * independently via cpml_stiffness_mode_acoustic -- e.g. run the
       * acoustic PML 'native' (consistent base + CPML correction, fixes
       * the lossless-fluid reduced-operator hourglass) while keeping the
       * elastic PML 'reduced'. */
      Param.theCpmlNativeStiffnessAcoustic = Param.theCpmlNativeStiffness;
      char cpml_stiffness_mode_ac[64] = "";
      if (parsetext(fp, "cpml_stiffness_mode_acoustic", 's',
                    &cpml_stiffness_mode_ac) == 0) {
        if (strcasecmp(cpml_stiffness_mode_ac, "native") == 0) {
          Param.theCpmlNativeStiffnessAcoustic = YES;
        } else if (strcasecmp(cpml_stiffness_mode_ac, "reduced") == 0) {
          Param.theCpmlNativeStiffnessAcoustic = NO;
        } else {
          solver_abort(__FUNCTION_NAME, NULL,
                       "cpml_stiffness_mode_acoustic must be 'native' or "
                       "'reduced': got '%s'\n", cpml_stiffness_mode_ac);
        }
      }

      if (cpml_thickness_m <= 0.0 && cpml_layers <= 0) {
        solver_abort(__FUNCTION_NAME, NULL,
                     "abc_type=cpml requires cpml_thickness_m > 0 or "
                     "cpml_layers > 0.\n");
      }
      if (cpml_reflection_target <= 0.0 || cpml_reflection_target >= 1.0) {
        solver_abort(__FUNCTION_NAME, NULL,
                     "cpml_reflection_target=%g must be in (0,1).\n",
                     cpml_reflection_target);
      }
      if (cpml_poly_order < 1 || cpml_poly_order > 4) {
        solver_abort(__FUNCTION_NAME, NULL,
                     "cpml_poly_order=%d must be in [1,4].\n",
                     cpml_poly_order);
      }
      if (cpml_kappa_max < 1.0) {
        solver_abort(__FUNCTION_NAME, NULL,
                     "cpml_kappa_max=%g must be >= 1.0.\n", cpml_kappa_max);
      }
      if (cpml_cfl_safety <= 0.0 || cpml_cfl_safety > 1.0) {
        solver_abort(__FUNCTION_NAME, NULL,
                     "cpml_cfl_safety=%g must be in (0,1].\n",
                     cpml_cfl_safety);
      }
      if (cpml_mpml_cross_fraction < 0.0 ||
          cpml_mpml_cross_fraction > 0.5) {
        solver_abort(__FUNCTION_NAME, NULL,
                     "cpml_mpml_cross_fraction=%g must be in [0, 0.5].\n",
                     cpml_mpml_cross_fraction);
      }
      if (cpml_alpha_floor_fraction < 0.0 ||
          cpml_alpha_floor_fraction > 1.0) {
        solver_abort(__FUNCTION_NAME, NULL,
                     "cpml_alpha_floor_fraction=%g must be in [0, 1].\n",
                     cpml_alpha_floor_fraction);
      }

      Param.theCpmlThicknessM = cpml_thickness_m;
      Param.theCpmlLayers = cpml_layers;
      Param.theCpmlReflectionTarget = cpml_reflection_target;
      Param.theCpmlPolyOrder = cpml_poly_order;
      Param.theCpmlKappaMax = cpml_kappa_max;
      Param.theCpmlAlphaMaxHz = cpml_alpha_max_hz;
      Param.theCpmlCflSafety = cpml_cfl_safety;
      Param.theCpmlMpmlCrossFraction = cpml_mpml_cross_fraction;
      Param.theCpmlAlphaFloorFraction = cpml_alpha_floor_fraction;

      /* Opt-in material damping inside the elastic PML zone (viscoelastic
       * CFS-CPML). Default "none" => all material damping skipped on PML
       * elements (the 5d91b09 baseline, byte-identical). "bkt"/"rayleigh"
       * feed the BKT / Rayleigh anelastic-displacement vectors through
       * the unsplit-displacement CPML stretch (stretch-consistent,
       * centroid-projected). Staged: parsed/plumbed here; only "none" is
       * wired in this commit (kernels added in later stages). */
      {
        char cpml_pml_mat[32] = "none";
        (void)parsetext(fp, "cpml_pml_material_damping", 's', &cpml_pml_mat);
        if (strcasecmp(cpml_pml_mat, "none") == 0) {
          Param.theCpmlPmlMaterialDamping = 0;
        } else if (strcasecmp(cpml_pml_mat, "bkt") == 0) {
          Param.theCpmlPmlMaterialDamping = 1;
        } else if (strcasecmp(cpml_pml_mat, "rayleigh") == 0) {
          Param.theCpmlPmlMaterialDamping = 2;
        } else {
          solver_abort(__FUNCTION_NAME, NULL,
                       "cpml_pml_material_damping must be none|bkt|"
                       "rayleigh: got '%s'\n", cpml_pml_mat);
        }
      }

      /* The CPML kernels and time-loop wiring are not yet committed; the
       * solver aborts just before the time-stepping loop in phase B so that
       * phase A (mesh_only_mode) and the CPML element-classification log
       * still complete normally. See the abort site near the main
       * "march forward in time" loop. */
    }
  }

  fclose(fp); /* physics.in */

  if ((fp = fopen(numericalin, "r")) == NULL) {
    fprintf(stderr, "Error opening %s\n", numericalin);
    return -1;
  }

  size_t monitor_name_len = 0;
  hu_config_get_string_def(fp, "monitor_file", &Param.theMonitorFileName,
                           &monitor_name_len, "monitor.txt");

  /* open the monitor file of the simulation in pe 0 */
  Param.theMonitorFileFp = fopen(Param.theMonitorFileName, "w");
  if (Param.theMonitorFileFp == NULL) {
    fprintf(stderr, "\n Err opening the monitor file");
  } else {
    setlinebuf(Param.theMonitorFileFp);
  }

  xfree_char(&Param.theMonitorFileName);

  /* numerical.in parse texts */
  if ((parsetext(fp, "simulation_wave_max_freq_hz", 'd', &freq) != 0) ||
      (parsetext(fp, "simulation_node_per_wavelength", 'i', &samples) != 0) ||
      (parsetext(fp, "simulation_shear_velocity_min", 'd', &vscut) != 0) ||
      (parsetext(fp, "simulation_end_time_sec", 'd', &endT) != 0) ||
      (parsetext(fp, "simulation_delta_time_sec", 'd', &deltaT) != 0) ||
      (parsetext(fp, "the_threshold_damping", 'd', &threshold_damping) != 0) ||
      (parsetext(fp, "the_threshold_Vp_over_Vs", 'd', &threshold_VpVs) != 0) ||
      /*  This parameter is used in solver_output_seq() function in psolve.c,
          but this function is not called in the code anywhere. */
      // (parsetext(fp, "4D_output_file", 's', &Param.FourDOutFile) != 0) ||
      (parsetext(fp, "cvmdb_input_file", 's', &Param.cvmdb_input_file) != 0)) {
    fprintf(stderr, "Error parsing simulation parameters from %s\n",
            numericalin);
    return -1;
  }
  /* Optional parameters */
  parsetext(fp, "mesh_etree_output_file", 's', &Param.mesh_etree_output_file);
  parsetext(fp, "simulation_output_rate", 'i', &rate);
  parsetext(fp, "checkpoint_path", 's', &checkpoint_path);
  /* If no plane input file specified, the plane input file is numercalin itself
   */
  if (parsetext(fp, "planes_input_file", 's', &Param.planes_input_file) != 0) {
    strcpy(Param.planes_input_file, numericalin);
  };
  /* These parameters have been initialized as 0. If they are
  not set in the input file, their values will just be 0. */
  parsetext(fp, "number_output_planes", 'i', &number_output_planes);
  parsetext(fp, "number_output_stations", 'i', &number_output_stations);
  parsetext(fp, "simulation_start_time_sec", 'd', &startT);
  parsetext(fp, "softening_factor", 'd', &softening_factor);
  parsetext(fp, "do_damping_statistics", 'i', &damping_statistics);
  parsetext(fp, "simulation_velocity_profile_freq_hz", 'd', &freq_vel);
  parsetext(fp, "use_checkpoint", 'i', &use_checkpoint);
  parsetext(fp, "checkpointing_rate", 'i', &checkpointing_rate);
  /* These parameters have been initialized as "no". If they are
  not set in the input file, their values will just be "no". */
  parsetext(fp, "print_matrix_k", 's', &print_matrix_k);
  parsetext(fp, "print_station_velocities", 's', &print_station_velocities);
  parsetext(fp, "print_station_accelerations", 's',
            &print_station_accelerations);
  parsetext(fp, "print_station_pressures", 's', &print_station_pressures);
  parsetext(fp, "mesh_coordinates_for_matlab", 's',
            &mesh_coordinates_for_matlab);
  parsetext(fp, "include_buildings", 's', &include_buildings);
  parsetext(fp, "implement_drm", 's', &implement_drm);
  parsetext(fp, "include_topography", 's', &include_topography);
  parsetext(fp, "include_acoustic", 's', &include_acoustic);
  parsetext(fp, "acoustic_interface_coupling", 's',
            &acoustic_interface_coupling);
  parsetext(fp, "acoustic_pressure_formulation", 's',
            &acoustic_pressure_formulation);
  parsetext(fp, "reconstruct_acoustic_plane_displacements", 's',
            &reconstruct_acoustic_plane_displacements);
  parsetext(fp, "acoustic_interface_diagnostics_file", 's',
            Param.theAcousticInterfaceDiagnosticsFile);
  parsetext(fp, "Istanbul_velocity_model", 's', &IstanbulModel);
  parsetext(fp, "basin_velocity_model", 's', &basinModel);
  parsetext(fp, "basin_material_type", 's', &basin_material_type);
  parsetext(fp, "3D_velocity_model", 's', &threeDModel);
  parsetext(fp, "include_nonlinear_analysis", 's', &include_nonlinear_analysis);
  parsetext(fp, "include_incident_planewaves", 's',
            &include_incident_planewaves);
  parsetext(fp, "include_hmg_halfspace", 's', &include_hmgHalfSpace);
  /* These parameters have been initialized with different default values
  and do not necessary to be set in the parameter input file. */
  parsetext(fp, "use_progressive_meshing", 'i', &step_meshing);
  parsetext(fp, "stiffness_calculation_method", 's',
            &stiffness_calculation_method);

  hu_config_get_int_opt(fp, "output_mesh", &Param.theMeshOutFlag);
  hu_config_get_int_opt(fp, "mesh_only_mode", &Param.theMeshOnlyMode);
  hu_config_get_int_opt(fp, "mesh_input_mode", &Param.theMeshInputMode);

  /* Acoustic-elastic coupling formulation. DEFAULT = velocity_potential:
   * it is the cross-validated, stable coupling (Z-disp corr 0.99-0.998 vs
   * CPS, eliminates the 0 km secular growth the pressure form exhibited).
   * The legacy pressure formulation remains available as an explicit
   * opt-in via `acoustic_formulation = pressure`. This flag only affects
   * simulations with acoustic elements / acoustic-elastic interfaces;
   * pure-elastic runs are unaffected and bit-identical. The interface
   * source filter stays gated off (found destabilizing). */
  {
    char acform[64] = "velocity_potential";
    char afilt[64] = "no";
    (void)parsetext(fp, "acoustic_formulation", 's', &acform);
    if (strcasecmp(acform, "pressure") == 0) {
      Param.theAcousticVelocityPotential = 0;
    } else if ((strcasecmp(acform, "velocity_potential") == 0) ||
               (strcasecmp(acform, "potential") == 0)) {
      Param.theAcousticVelocityPotential = 1;
    } else {
      solver_abort(__FUNCTION_NAME, NULL,
                   "acoustic_formulation must be 'pressure' or "
                   "'velocity_potential': got '%s'\n", acform);
    }
    (void)parsetext(fp, "acoustic_interface_filter", 's', &afilt);
    Param.theAcousticInterfaceFilter =
        (strcasecmp(afilt, "yes") == 0) ? YES : NO;
    hu_config_get_double_opt(fp, "acoustic_interface_filter_alpha",
                             &Param.theAcousticInterfaceFilterAlpha);
    hu_config_get_double_opt(fp, "acoustic_interface_filter_strength",
                             &Param.theAcousticInterfaceFilterStrength);
    /* Opt-in acoustic (water) attenuation. acoustic_quality_factor = a
     * quality factor Q for acoustic elements; > 0 enables a small
     * stiffness-proportional (Kelvin-Voigt) damping in the acoustic
     * force kernel (NOT the elastic zeta/a*M route). Unset / <= 0 keeps
     * acoustic exactly lossless and bit-identical to the validated
     * baseline. A numerical regularizer for the lossless-cavity grid
     * ring; not physical ocean attenuation (real seawater Q at the
     * seismic band is ~1e5). See README. */
    hu_config_get_double_opt(fp, "acoustic_quality_factor",
                             &Param.theAcousticQualityFactor);
  }
  hu_config_get_int_opt(fp, "maximum_octant_level",
                        &acoustic_interface_max_level);
  hu_config_get_int_opt(fp, "enable_timing_barriers",
                        &Param.theTimingBarriersFlag);
  hu_config_get_int_opt(fp, "forces_buffer_size", &theForcesBufferSize);
  hu_config_get_int_opt(fp, "schedule_print_file",
                        &Param.theSchedulePrintToFile);
  hu_config_get_double_opt(fp, "basin_water_vp", &Param.basinWaterVp);
  hu_config_get_double_opt(fp, "basin_water_rho", &Param.basinWaterRho);
  hu_config_get_double_opt(fp, "acoustic_mesh_refinement_factor",
                           &acoustic_mesh_refinement_factor);

  hu_config_get_int_opt(fp, "schedule_print_error_check",
                        &Param.theSchedulePrintErrorCheckFlag);
  hu_config_get_int_opt(fp, "schedule_print_stdout",
                        &Param.theSchedulePrintToStdout);

  size_t schedule_stat_len = 0;
  hu_config_get_string_def(fp, "stat_schedule_filename",
                           &Param.theScheduleStatFilename, &schedule_stat_len,
                           "stat-sched.txt");
  size_t mesh_stat_len = 0;
  hu_config_get_string_def(fp, "stat_mesh_filename", &Param.theMeshStatFilename,
                           &mesh_stat_len, "stat-mesh.txt");

  size_t topo_stat_len = 0;
  hu_config_get_string_def(fp, "stat_topo_filename", &Param.theTopoStatFilename,
                           &topo_stat_len, "stat-topo.txt");

  /* sanity check */
  if (freq <= 0) {
    fprintf(stderr, "Illegal frequency value %f\n", freq);
    return -1;
  }

  if (freq_vel < 0 || freq_vel > freq) {
    fprintf(stderr,
            "Illegal frequency value, velocity profile frequency can not be "
            "smaller than zero or bigger than max freq %f\n",
            freq_vel);
    return -1;
  }

  if (samples <= 0) {
    fprintf(stderr, "Illegal samples value %d\n", samples);
    return -1;
  }

  if (vscut <= 0) {
    fprintf(stderr, "Illegal vscut value %f\n", vscut);
    return -1;
  }

  if ((startT < 0) || (endT < 0) || (startT > endT)) {
    fprintf(stderr, "Illegal startT %f or endT %f\n", startT, endT);
    return -1;
  }

  if (deltaT <= 0) {
    fprintf(stderr, "Illegal deltaT %f\n", deltaT);
    return -1;
  }

  if ((softening_factor <= 1) && (softening_factor != 0)) {
    fprintf(stderr, "Illegal softening factor %f\n", softening_factor);
    return -1;
  }

  if (step_meshing < 0) {
    fprintf(stderr, "Illegal progressive meshing factor %d\n", step_meshing);
    return -1;
  }

  if ((acoustic_interface_max_level < -1) ||
      (acoustic_interface_max_level > PIXELLEVEL)) {
    fprintf(stderr, "Illegal acoustic interface maximum octant level %d\n",
            acoustic_interface_max_level);
    return -1;
  }

  if (acoustic_mesh_refinement_factor < 1.0) {
    fprintf(stderr, "Illegal acoustic mesh refinement factor %f\n",
            acoustic_mesh_refinement_factor);
    return -1;
  }

  if (rate <= 0) {
    fprintf(stderr, "Illegal output rate %d\n", rate);
    return -1;
  }

  if (number_output_planes < 0) {
    fprintf(stderr, "Illegal number of output planes %d\n",
            number_output_planes);
    return -1;
  } else if (number_output_planes > 0) {
    parsetext(fp, "output_planes_directory", 's', &Param.thePlaneDirOut);
    createDir(Param.thePlaneDirOut, __FUNCTION_NAME);
  }

  if (number_output_stations < 0) {
    fprintf(stderr, "Illegal number of output stations %d\n",
            number_output_stations);
    return -1;
  }

  if (threshold_damping < 0) {
    fprintf(stderr, "Illegal threshold damping %f\n", threshold_damping);
    return -1;
  }

  if (threshold_VpVs < 0) {
    fprintf(stderr, "Illegal threshold Vp over Vs %f\n", threshold_VpVs);
    return -1;
  }

  if ((damping_statistics < 0) || (damping_statistics > 1)) {
    fprintf(stderr, "Illegal do damping statistics flag %d\n",
            damping_statistics);
    return -1;
  }

  if ((use_checkpoint < 0) || (use_checkpoint > 1)) {
    fprintf(stderr, "Illegal use checkpoint flag %d\n", use_checkpoint);
    return -1;
  }

  if (checkpointing_rate < 0) {
    fprintf(stderr, "Illegal checkpointing rate %d\n", use_checkpoint);
    return -1;
  }

  if (strcasecmp(include_nonlinear_analysis, "yes") == 0) {
    includeNonlinear = YES;
  } else if (strcasecmp(include_nonlinear_analysis, "no") == 0) {
    includeNonlinear = NO;
  } else {
    solver_abort(__FUNCTION_NAME, NULL,
                 "Unknown response for including"
                 "nonlinear analysis (yes or no): %s\n",
                 include_nonlinear_analysis);
  }

  if (strcasecmp(stiffness_calculation_method, "effective") == 0) {
    stiffness_method = EFFECTIVE;
  } else if (strcasecmp(stiffness_calculation_method, "conventional") == 0) {
    stiffness_method = CONVENTIONAL;
  } else {
    solver_abort(__FUNCTION_NAME, NULL,
                 "Unknown response for stiffness"
                 "calculation method (effective or conventional): %s\n",
                 stiffness_calculation_method);
  }

  if (strcasecmp(print_matrix_k, "yes") == 0) {
    printMatrixK = YES;
  } else if (strcasecmp(print_matrix_k, "no") == 0) {
    printMatrixK = NO;
  } else {
    solver_abort(__FUNCTION_NAME, NULL,
                 "Unknown response for printing K matrix (yes or no): %s\n",
                 print_matrix_k);
  }

  if (strcasecmp(print_station_velocities, "yes") == 0) {
    printStationVels = YES;
  } else if (strcasecmp(print_station_velocities, "no") == 0) {
    printStationVels = NO;
  } else {
    solver_abort(__FUNCTION_NAME, NULL,
                 "Unknown response for printing"
                 "station velocities (yes or no): %s\n",
                 print_station_velocities);
  }

  if (strcasecmp(print_station_accelerations, "yes") == 0) {
    printStationAccs = YES;
  } else if (strcasecmp(print_station_accelerations, "no") == 0) {
    printStationAccs = NO;
  } else {
    solver_abort(__FUNCTION_NAME, NULL,
                 "Unknown response for printing"
                 "station accelerations (yes or no): %s\n",
                 print_station_accelerations);
  }

  if (strcasecmp(print_station_pressures, "yes") == 0) {
    printStationPressures = YES;
  } else if (strcasecmp(print_station_pressures, "no") == 0) {
    printStationPressures = NO;
  } else {
    solver_abort(__FUNCTION_NAME, NULL,
                 "Unknown response for printing"
                 "station pressures (yes or no): %s\n",
                 print_station_pressures);
  }

  if (strcasecmp(mesh_coordinates_for_matlab, "yes") == 0) {
    meshCoordinatesForMatlab = YES;
    parsetext(fp, "mesh_coordinates_directory_for_matlab", 's',
              &Param.theMeshDirMatlab);
    createDir(Param.theMeshDirMatlab, __FUNCTION_NAME);
  } else if (strcasecmp(mesh_coordinates_for_matlab, "no") == 0) {
    meshCoordinatesForMatlab = NO;
  } else {
    solver_abort(__FUNCTION_NAME, NULL,
                 "Unknown response for mesh coordinates"
                 "for matlab (yes or no): %s\n",
                 mesh_coordinates_for_matlab);
  }

  if (strcasecmp(include_buildings, "yes") == 0) {
    have_buildings = YES;
  } else if (strcasecmp(include_buildings, "no") == 0) {
    have_buildings = NO;
  } else {
    solver_abort(__FUNCTION_NAME, NULL,
                 "Unknown response for including buildings (yes or no): %s\n",
                 include_buildings);
  }

  if (strcasecmp(implement_drm, "yes") == 0) {
    implementdrm = YES;
    parsetext(fp, "drm_directory", 's', &drm_directory);
    if (parsetext(fp, "which_drm_part", 's', &which_drm_part) != 0) {
      fprintf(stderr, "Error reading drm part from %s\n", numericalin);
      return -1;
    }
    // Make sure ``which_drm_part`` is either ``part0``, ``part1``, or ``part2``
    if (strcasecmp(which_drm_part, "part0") != 0 &&
        strcasecmp(which_drm_part, "part1") != 0 &&
        strcasecmp(which_drm_part, "part2") != 0) {
      solver_abort(
          __FUNCTION_NAME, NULL,
          "Unknown response for which_drm_part (part0, part1, or part2): %s\n",
          which_drm_part);
    }
    // Create the directory if it does not exist
    snprintf(drm_directory_full_path, sizeof(drm_directory_full_path), "%s/%s",
             drm_directory, which_drm_part);
    createDir(drm_directory_full_path, __FUNCTION_NAME);
  } else if (strcasecmp(implement_drm, "no") == 0) {
    implementdrm = NO;
  } else {
    solver_abort(__FUNCTION_NAME, NULL,
                 "Unknown response for impelement_drm (yes or no): %s\n",
                 implement_drm);
  }

  if (strcasecmp(include_topography, "yes") == 0) {
    haveTopography = YES;
    parsetext(fp, "topography_directory", 's', &Param.theTopoDir);
  } else if (strcasecmp(include_topography, "no") == 0) {
    haveTopography = NO;
  } else {
    solver_abort(__FUNCTION_NAME, NULL,
                 "Unknown response for include_topography (yes or no): %s\n",
                 include_topography);
  }

  if (strcasecmp(include_acoustic, "yes") == 0) {
    includeAcoustic = YES;
  } else if (strcasecmp(include_acoustic, "no") == 0) {
    includeAcoustic = NO;
  } else {
    solver_abort(__FUNCTION_NAME, NULL,
                 "Unknown response for include_acoustic (yes or no): %s\n",
                 include_acoustic);
  }

  if (strcasecmp(acoustic_interface_coupling, "yes") == 0) {
    acousticInterfaceCoupling = YES;
  } else if (strcasecmp(acoustic_interface_coupling, "no") == 0) {
    acousticInterfaceCoupling = NO;
  } else {
    solver_abort(__FUNCTION_NAME, NULL,
                 "Unknown response for acoustic_interface_coupling (yes or no): "
                 "%s\n",
                 acoustic_interface_coupling);
  }

  if (strcasecmp(acoustic_pressure_formulation, "yes") == 0) {
    acousticPressureFormulation = YES;
  } else if (strcasecmp(acoustic_pressure_formulation, "no") == 0) {
    acousticPressureFormulation = NO;
  } else {
    solver_abort(__FUNCTION_NAME, NULL,
                 "Unknown response for acoustic_pressure_formulation "
                 "(yes or no): %s\n",
                 acoustic_pressure_formulation);
  }

  if (strcasecmp(reconstruct_acoustic_plane_displacements, "yes") == 0) {
    reconstructAcousticPlaneDisplacements = YES;
  } else if (strcasecmp(reconstruct_acoustic_plane_displacements, "no") == 0) {
    reconstructAcousticPlaneDisplacements = NO;
  } else {
    solver_abort(__FUNCTION_NAME, NULL,
                 "Unknown response for reconstruct_acoustic_plane_displacements "
                 "(yes or no): %s\n",
                 reconstruct_acoustic_plane_displacements);
  }

  if ((strcasecmp(basin_material_type, "elastic") == 0) ||
      (strcasecmp(basin_material_type, "default") == 0)) {
    basinMaterialType = BASIN_MATERIAL_ELASTIC;
  } else if ((strcasecmp(basin_material_type, "acoustic_water") == 0) ||
             (strcasecmp(basin_material_type, "water") == 0)) {
    basinMaterialType = BASIN_MATERIAL_ACOUSTIC_WATER;
  } else {
    solver_abort(__FUNCTION_NAME, NULL,
                 "Unknown response for basin_material_type "
                 "(elastic or acoustic_water): %s\n",
                 basin_material_type);
  }

  if ((acousticPressureFormulation == YES) && (includeAcoustic != YES)) {
    solver_abort(__FUNCTION_NAME, NULL,
                 "acoustic_pressure_formulation=yes requires "
                 "include_acoustic=yes\n");
  }
  if ((printStationPressures == YES) &&
      ((includeAcoustic != YES) || (acousticPressureFormulation != YES))) {
    solver_abort(__FUNCTION_NAME, NULL,
                 "print_station_pressures=yes requires include_acoustic=yes "
                 "and acoustic_pressure_formulation=yes\n");
  }
  if ((reconstructAcousticPlaneDisplacements == YES) &&
      ((includeAcoustic != YES) || (acousticPressureFormulation != YES))) {
    solver_abort(__FUNCTION_NAME, NULL,
                 "reconstruct_acoustic_plane_displacements=yes requires "
                 "include_acoustic=yes and acoustic_pressure_formulation=yes\n");
  }
  if ((basinMaterialType == BASIN_MATERIAL_ACOUSTIC_WATER) &&
      (includeAcoustic != YES)) {
    solver_abort(__FUNCTION_NAME, NULL,
                 "basin_material_type=acoustic_water requires "
                 "include_acoustic=yes\n");
  }
  if ((basinMaterialType == BASIN_MATERIAL_ACOUSTIC_WATER) &&
      ((Param.basinWaterVp <= 0.0) || (Param.basinWaterRho <= 0.0))) {
    solver_abort(__FUNCTION_NAME, NULL,
                 "basin_water_vp and basin_water_rho must be positive "
                 "for basin_material_type=acoustic_water\n");
  }

  if (strcasecmp(include_incident_planewaves, "yes") == 0) {
    includePlaneWaves = YES;
  } else if (strcasecmp(include_incident_planewaves, "no") == 0) {
    includePlaneWaves = NO;
  } else {
    solver_abort(
        __FUNCTION_NAME, NULL,
        "Unknown response for include_incident_planewaves (yes or no): %s\n",
        include_incident_planewaves);
  }

  if (strcasecmp(include_hmgHalfSpace, "yes") == 0) {
    includeHmgHalfSpace = YES;
  } else if (strcasecmp(include_hmgHalfSpace, "no") == 0) {
    includeHmgHalfSpace = NO;
  } else {
    solver_abort(__FUNCTION_NAME, NULL,
                 "Unknown response for include_hmg_halfspace (yes or no): %s\n",
                 include_hmgHalfSpace);
  }

  if (strcasecmp(IstanbulModel, "yes") == 0) {
    includeIstanbulVelModel = YES;
  } else if (strcasecmp(IstanbulModel, "no") == 0) {
    includeIstanbulVelModel = NO;
  } else {
    solver_abort(
        __FUNCTION_NAME, NULL,
        "Unknown response for Istanbul_velocity_model (yes or no): %s\n",
        includeIstanbulVelModel);
  }

  if (strcasecmp(basinModel, "yes") == 0) {
    includeBasinVelModel = YES;
    if ((parsetext(fp, "basin_origin_latitude_deg", 'd', &Param.theBasinLat) !=
         0) ||
        (parsetext(fp, "basin_origin_longitude_deg", 'd',
                   &Param.theBasinLong) != 0)) {
      fprintf(stderr, "Error reading basin origin from %s\n", numericalin);
      return -1;
    }
    if (parsetext(fp, "basin_input_file", 's', &Param.basin_input_file) == 0) {
      Param.basinData = readBasin(Param.basin_input_file, &Param.basinXCount,
                                  &Param.basinYCount, &Param.basinZCount,
                                  &Param.basinIncrement);
    } else {
      fprintf(stderr, "Error reading basin input file from %s\n", numericalin);
      return -1;
    }
  } else if (strcasecmp(basinModel, "no") == 0) {
    includeBasinVelModel = NO;
  } else {
    solver_abort(__FUNCTION_NAME, NULL,
                 "Unknown response for basin_velocity_model (yes or no): %s\n",
                 includeBasinVelModel);
  }

  /* 3D Velocity Model */
  if (strcasecmp(threeDModel, "yes") == 0) {
    include3DVelModel = YES;
    parsetext(fp, "3D_velocity_model_directory", 's', &Param.the3DVelModelDir);
    if ((parsetext(fp, "3D_velocity_model_origin_latitude_deg", 'd',
                   &Param.the3DVelocityModelLat) != 0) ||
        (parsetext(fp, "3D_velocity_model_origin_longitude_deg", 'd',
                   &Param.the3DVelocityModelLong) != 0)) {
      fprintf(stderr, "Error reading 3D velocity model origin from %s\n",
              numericalin);
      return -1;
    }
  } else if (strcasecmp(threeDModel, "no") == 0) {
    include3DVelModel = NO;
  } else {
    solver_abort(__FUNCTION_NAME, NULL,
                 "Unknown response for 3D_velocity_model (yes or no): %s\n",
                 include3DVelModel);
  }

  parsetext(fp, "source_directory_output", 's', &Param.theSourceOutputDir);
  createDir(Param.theSourceOutputDir, __FUNCTION_NAME);
  if (Param.theAcousticInterfaceDiagnosticsFile[0] != '\0') {
    create_parent_dir_for_file(Param.theAcousticInterfaceDiagnosticsFile,
                               __FUNCTION_NAME);
  }

/* read domain corners */
#ifndef PROJ
  int iCorner;
  double* auxiliar;
  auxiliar = (double*)malloc(sizeof(double) * 8);

  if (parsedarray(fp, "domain_surface_corners", 8, auxiliar) != 0) {
    solver_abort(__FUNCTION_NAME, NULL,
                 "Error parsing domain_surface_corners field from %s\n",
                 numericalin);
  }

  for (iCorner = 0; iCorner < 4; iCorner++) {
    Param.theSurfaceCornersLong[iCorner] = auxiliar[iCorner * 2];
    Param.theSurfaceCornersLat[iCorner] = auxiliar[iCorner * 2 + 1];
  }
  free(auxiliar);
#endif

  fclose(fp);

  /* Check whether a profile will be used and init global flag to YES */
  if (strcasecmp(Param.cvmdb_input_file, "profile") == 0) {
    Param.useProfile = YES;
  }

  /* Init the static global variables */

  Param.theRegionLat = region_origin_latitude_deg;
  Param.theRegionLong = region_origin_longitude_deg;
  Param.theRegionDepth = region_depth_shallow_m;

  Param.theVsCut = vscut;
  Param.theQConstant = qconstant;
  Param.theQAlpha = qalpha;
  Param.theQBeta = qbeta;
  Param.theFactor = freq * samples;
  Param.theFreq = freq;
  Param.theFreq_Vel = freq_vel;
  Param.theDeltaT = deltaT;
  Param.theDeltaTSquared = deltaT * deltaT;
  Param.theStartT = startT;
  Param.theEndT = endT;
  Param.theTotalSteps = (int)(((endT - startT) / deltaT));

  Param.theDomainX = region_length_north_m;
  Param.theDomainY = region_length_east_m;
  Param.region_depth_deep_m = region_depth_deep_m;
  Param.theDomainZ = region_depth_deep_m - region_depth_shallow_m;
  Param.theDomainAzimuth = region_azimuth_leftface_deg;
  Param.theTypeOfDamping = typeOfDamping;
  Param.useInfQk = useInfQk;

  /* Auto-enable viscoelastic (BKT) PML for the BKT*+CPML combination.
   * type_of_damping=bkt(2/3) + abc_type=cpml is catastrophically
   * UNSTABLE on the material-damping-skipped-on-PML baseline: the field
   * diverges to NaN by t~40s (confirmed: v5 @23f58bf == v3, byte-
   * identical blow-up). The viscoelastic CFS-CPML
   * (cpml_pml_material_damping=bkt) is the validated fix (Stage 2:
   * stable to 240s, corr 0.994 vs the matching-attenuation CPS). There
   * is no valid reason to run bkt*+cpml without it, so when it is left
   * unset/none(0) for that combo we auto-promote to bkt(1) and log
   * loudly. Explicit =bkt is already 1 (no-op); explicit =rayleigh(2)
   * with BKT damping stays a real misconfig caught by the abort below.
   * Only the previously-NaN combo changes behavior; type_of_damping<BKT
   * (incl. all prior =none CPML validation) is untouched. */
  if (Param.theAbcType == ABC_CPML && Param.theTypeOfDamping >= BKT &&
      Param.theCpmlPmlMaterialDamping == 0) {
    Param.theCpmlPmlMaterialDamping = 1;
    if (Global.myID == 0) {
      monitor_print(
          "CPML: auto-enabled viscoelastic PML "
          "(cpml_pml_material_damping=bkt) -- type_of_damping is "
          "BKT-family with abc_type=cpml, which is numerically unstable "
          "(field -> NaN ~t40s) without it; using the validated "
          "viscoelastic CFS-CPML.\n");
    }
  }

  /* Stage-2 viscoelastic CPML cross-validation: cpml_pml_material_damping
   * only does something if the interior actually carries the matching
   * material damping (so the PML elements have nonzero anelastic vectors
   * to stretch). Fail fast on a misconfiguration rather than waste a run
   * (the kernel NULL-guard would otherwise silently no-op it). Both
   * params are set by here (CPML parse ~line 2140; typeOfDamping above). */
  if (Param.theCpmlPmlMaterialDamping == 1 && Param.theTypeOfDamping < BKT) {
    solver_abort(__FUNCTION_NAME, NULL,
                 "cpml_pml_material_damping=bkt requires a BKT "
                 "type_of_damping (bkt/bkt2/bkt3); got type_of_damping "
                 "enum %d\n", (int)Param.theTypeOfDamping);
  }
  if (Param.theCpmlPmlMaterialDamping == 2 &&
      Param.theTypeOfDamping != RAYLEIGH) {
    solver_abort(__FUNCTION_NAME, NULL,
                 "cpml_pml_material_damping=rayleigh requires "
                 "type_of_damping=rayleigh; got enum %d\n",
                 (int)Param.theTypeOfDamping);
  }

  Param.theRate = rate;

  Param.theNumberOfPlanes = number_output_planes;
  Param.theNumberOfStations = number_output_stations;

  Param.theSofteningFactor = softening_factor;
  Param.theStepMeshingFactor = step_meshing;
  Param.theAcousticInterfaceMaxLevel = acoustic_interface_max_level;
  Param.theAcousticMeshRefinementFactor = acoustic_mesh_refinement_factor;
  Param.theThresholdDamping = threshold_damping;
  Param.theThresholdVpVs = threshold_VpVs;
  Param.theDampingStatisticsFlag = damping_statistics;

  Param.theCheckPointingRate = checkpointing_rate;
  Param.theUseCheckPoint = use_checkpoint;

  Param.includeNonlinearAnalysis = includeNonlinear;
  Param.theStiffness = stiffness_method;

  Param.printK = printMatrixK;
  Param.printStationVelocities = printStationVels;
  Param.printStationAccelerations = printStationAccs;
  Param.printStationPressures = printStationPressures;

  Param.includeBuildings = have_buildings;
  Param.useParametricQ = have_parametricq;

  Param.storeMeshCoordinatesForMatlab = meshCoordinatesForMatlab;

  Param.drmImplement = implementdrm;

  Param.includeTopography = haveTopography;
  Param.includeAcoustic = includeAcoustic;
  Param.acousticInterfaceCoupling = acousticInterfaceCoupling;
  Param.acousticPressureFormulation = acousticPressureFormulation;
  Param.reconstructAcousticPlaneDisplacements =
      reconstructAcousticPlaneDisplacements;
  Param.includeIncidentPlaneWaves = includePlaneWaves;
  Param.includeHomogeneousHalfSpace = includeHmgHalfSpace;

  Param.IstanbulVelocityModel = includeIstanbulVelModel;
  Param.basinVelocityModel = includeBasinVelModel;
  Param.basinMaterialType = basinMaterialType;
  Param.threeDVelocityModel = include3DVelModel;

  strcpy(Param.theCheckPointingDirOut, checkpoint_path);
  if (Param.theCheckPointingRate != 0) {
    // Create the directory if it does not exist
    createDir(Param.theCheckPointingDirOut, __FUNCTION_NAME);
  }

  monitor_print("\n\n---------------- Some Input Data ----------------\n\n");
  monitor_print("Vs cut:                             %f\n", Param.theVsCut);
  monitor_print("Softening factor:                   %f\n",
                Param.theSofteningFactor);
  monitor_print("Number of stations:                 %d\n",
                Param.theNumberOfStations);
  monitor_print("Number of planes:                   %d\n",
                Param.theNumberOfPlanes);
  monitor_print("Stiffness calculation method:       %s\n",
                stiffness_calculation_method);
  monitor_print("Include buildings:                  %s\n", include_buildings);
  monitor_print("Include nonlinear analysis:         %s\n",
                include_nonlinear_analysis);
  monitor_print("Printing velocities on stations:    %s\n",
                print_station_velocities);
  monitor_print("Printing accelerations on stations: %s\n",
                print_station_accelerations);
  monitor_print("Printing pressures on stations:     %s\n",
                print_station_pressures);
  monitor_print("Mesh Coordinates For Matlab:        %s\n",
                mesh_coordinates_for_matlab);
  monitor_print("Mesh-Mode Only:                     %s\n",
                (Param.theMeshOnlyMode != 0) ? "yes" : "no");
  monitor_print("cvmdb_input_file:                   %s\n",
                Param.cvmdb_input_file);
  monitor_print("Implement drm:                      %s\n", implement_drm);
  monitor_print("Include Topography:                 %s\n", include_topography);
  monitor_print("Include Acoustic:                   %s\n", include_acoustic);
  if (Param.includeAcoustic == YES) {
    monitor_print("Acoustic interface coupling:        %s\n",
                  acoustic_interface_coupling);
    monitor_print("Acoustic pressure formulation:      %s\n",
                  acoustic_pressure_formulation);
    monitor_print("Acoustic plane displacement recon.: %s\n",
                  reconstruct_acoustic_plane_displacements);
    monitor_print("Acoustic mesh refinement factor:    %f\n",
                  Param.theAcousticMeshRefinementFactor);
    monitor_print("Acoustic formulation:               %s\n",
                  (Param.theAcousticVelocityPotential == 1)
                      ? "velocity_potential (p = -dphi/dt; vel continuity) [DEFAULT]"
                      : "pressure (legacy; explicit opt-in)");
    if (Param.theAcousticQualityFactor > 0.0) {
      monitor_print("Acoustic quality factor Q:          %.3f "
                    "(stiffness-proportional; numerical regularizer)\n",
                    Param.theAcousticQualityFactor);
    } else {
      monitor_print("Acoustic quality factor Q:          %s\n",
                    "lossless (default; bit-identical baseline)");
    }
    monitor_print("Acoustic interface source filter:   %s",
                  (Param.theAcousticInterfaceFilter == YES) ? "yes" : "no");
    if (Param.theAcousticInterfaceFilter == YES) {
      monitor_print(" (alpha=%.4f strength=%.4f)",
                    Param.theAcousticInterfaceFilterAlpha,
                    Param.theAcousticInterfaceFilterStrength);
    }
    monitor_print("\n");
    if (Param.theAcousticInterfaceDiagnosticsFile[0] != '\0') {
      monitor_print("Acoustic interface diagnostics:     %s\n",
                    Param.theAcousticInterfaceDiagnosticsFile);
    }
    if (Param.theAcousticInterfaceMaxLevel >= 0) {
      monitor_print("Acoustic interface max level:       %d\n",
                    Param.theAcousticInterfaceMaxLevel);
    } else {
      monitor_print("Acoustic interface max level:       disabled\n");
    }
  }
  monitor_print("Include PML:                        %s\n",
                (Param.theAbcType == ABC_CPML) ? "yes" : "no");
  if (Param.theAbcType == ABC_CPML) {
    monitor_print("PML layers:                         %d\n",
                  Param.theCpmlLayers);
    monitor_print("PML material damping:               %s\n",
                  (Param.theCpmlPmlMaterialDamping == 1)   ? "bkt"
                  : (Param.theCpmlPmlMaterialDamping == 2) ? "rayleigh"
                  : "none (skipped; baseline)");
    monitor_print("PML thickness override (m):         %f\n",
                  Param.theCpmlThicknessM);
    monitor_print("PML reflection coefficient:         %e\n",
                  Param.theCpmlReflectionTarget);
    monitor_print("PML polynomial order:               %d\n",
                  Param.theCpmlPolyOrder);
    monitor_print("PML kappa max:                      %f\n",
                  Param.theCpmlKappaMax);
    if (Param.theCpmlAlphaMaxHz > 0.0) {
      monitor_print("PML alpha max (Hz):                 %f\n",
                    Param.theCpmlAlphaMaxHz);
    } else {
      monitor_print("PML alpha max (Hz):                 auto "
                    "(pi * simulation_wave_max_freq_hz)\n");
    }
    monitor_print("PML multiaxial cross fraction:      %f\n",
                  Param.theCpmlMpmlCrossFraction);
    monitor_print("PML alpha floor fraction:           %f\n",
                  Param.theCpmlAlphaFloorFraction);
    monitor_print("PML CFL safety factor:              %f\n",
                  Param.theCpmlCflSafety);
    monitor_print("PML outer Dirichlet clamp:          %s\n",
                  (Param.theCpmlClampOuter == YES) ? "yes" : "no");
    monitor_print("PML stiffness mode (elastic):       %s\n",
                  (Param.theCpmlNativeStiffness == YES)
                      ? "native+correction"
                      : "reduced (standalone 1-pt)");
    monitor_print("PML stiffness mode (acoustic):      %s\n",
                  (Param.theCpmlNativeStiffnessAcoustic == YES)
                      ? "native+correction"
                      : "reduced (standalone 1-pt)");
    monitor_print("PML formulation:                    "
                  "two-sided CFS-CPML (Matzen 2011)\n");
  }
  monitor_print("Include Incident Plane Waves:       %s\n",
                include_incident_planewaves);
  monitor_print("Include Homogeneous Halfspace:      %s\n",
                include_hmgHalfSpace);
  monitor_print("Include Istanbul Velocity model:    %s\n", IstanbulModel);
  monitor_print("Include Basin Velocity model:       %s\n", basinModel);
  if (Param.basinVelocityModel == YES) {
    monitor_print("Basin material type:                %s\n",
                  basin_material_type);
    if (Param.basinMaterialType == BASIN_MATERIAL_ACOUSTIC_WATER) {
      monitor_print("Basin water Vp:                     %f\n",
                    Param.basinWaterVp);
      monitor_print("Basin water rho:                    %f\n",
                    Param.basinWaterRho);
    }
  }
  monitor_print("Include 3D Velocity model:          %s\n", threeDModel);
#ifdef PROJ
  monitor_print("Use PROJ library:                   yes\n");
#endif
  if (typeOfDamping >= BKT) {
    monitor_print("Use Parametric Q factor:            %s\n", use_parametricq);
    if (have_parametricq == YES) {
      monitor_print("Constant value of Q factor:         %f\n",
                    Param.theQConstant);
      monitor_print("Alpha value of Q factor:            %f\n",
                    Param.theQAlpha);
      monitor_print("Beta value of Q factor:             %f\n", Param.theQBeta);
    }
  }
  monitor_print("\n-------------------------------------------------\n\n");

  fflush(Param.theMonitorFileFp);

  return 0;
}

/*-----------Mesh generation related routines------------------------------*/

/* In case a profile is used... */

typedef enum {
  PROFILE_VP = 0,
  PROFILE_VS,
  PROFILE_RHO,
  PROFILE_QP,
  PROFILE_QS
} profile_flag_t;

/**
 * Interpolates the values of the profile property given the corresponding flag
 * and the indices in the tables.
 *
 * @return float of property
 */
float interpolate_profile_property(double depth, int lowerIndex, int upperIndex,
                                   int flag) {
  static const char* fname = __FUNCTION_NAME;
  double* theProfileProperty = NULL;
  double theProperty;

  switch (flag) {
    case PROFILE_VP:
      theProfileProperty = Param.theProfileVp;
      break;
    case PROFILE_VS:
      theProfileProperty = Param.theProfileVs;
      break;
    case PROFILE_RHO:
      theProfileProperty = Param.theProfileRho;
      break;
    case PROFILE_QP:
      theProfileProperty = Param.theProfileQp;
      break;
    case PROFILE_QS:
      theProfileProperty = Param.theProfileQs;
      break;
    default:
      solver_abort(fname, NULL,
                   "Error on interpolation of profile properties\n");
  }
  /* uncomment if linear interpolation is needed
   theProperty = theProfileProperty[lowerIndex]
               + ( depth - Param.theProfileZ[lowerIndex] )
               * ( theProfileProperty[upperIndex] -
   theProfileProperty[lowerIndex] ) / ( Param.theProfileZ[upperIndex] -
   Param.theProfileZ[lowerIndex] ); */

  theProperty = theProfileProperty[lowerIndex];

  return (float)theProperty;
}

/**
 * Finds the indices between which the depth queried is located and calls
 * interpolate_profile to return the properties
 *
 * @return 0 on success
 * @return 1 on failure
 */
int profile_query(double depth, cvmpayload_t* props) {
  int i, last;

  last = Param.theNumberOfLayers - 1;

  if (depth < 0) return 1;

  if (depth >= Param.theProfileZ[last]) {
    props->Vp = Param.theProfileVp[last];
    props->Vs = Param.theProfileVs[last];
    props->rho = Param.theProfileRho[last];
    props->Qp = Param.theProfileQp[last];
    props->Qs = Param.theProfileQs[last];
    return 0;
  }

  for (i = 0; i < last; i++) {
    if ((depth >= Param.theProfileZ[i]) && (depth < Param.theProfileZ[i + 1])) {
      // props->Vp = interpolate_profile_property(depth, i, i + 1, PROFILE_VP);
      // props->Vs = interpolate_profile_property(depth, i, i + 1, PROFILE_VS);
      // props->rho = interpolate_profile_property(depth, i, i + 1,
      // PROFILE_RHO); props->Qp = interpolate_profile_property(depth, i, i + 1,
      // PROFILE_QP); props->Qs = interpolate_profile_property(depth, i, i + 1,
      // PROFILE_QS);

      // NOTE: Since we are not interpolating now, we can just use the value at
      // the lower index and make the overall performance a bit better.
      props->Vp = Param.theProfileVp[i];
      props->Vs = Param.theProfileVs[i];
      props->rho = Param.theProfileRho[i];
      props->Qp = Param.theProfileQp[i];
      props->Qs = Param.theProfileQs[i];

      return 0;
    }
  }

  return 1;
}

/* In case an etree is used...*/

static void open_cvmdb(void) {
#ifdef USECVMDB

  MPI_Barrier(comm_solver);
  replicateDB(Param.cvmdb_input_file);

  MPI_Barrier(comm_solver);
  Global.theCVMEp =
      etree_open(Param.cvmdb_input_file, O_RDONLY, CVMBUFSIZE, 0, 0);

  if (Global.theCVMEp == NULL) {
    fprintf(stderr, "Thread %d: open_cvmdb: error opening CVM etree %s\n",
            Global.myID, Param.cvmdb_input_file);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  dbctl_t* myctl;
  /* Obtain the material database application control/meta data */
  if ((myctl = cvm_getdbctl(Global.theCVMEp)) == NULL) {
    fprintf(stderr, "Error reading CVM etree control data\n");
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  /* Check the ranges of the mesh and the scope of the CVM etree */
  if ((Param.theRegionLat < myctl->region_origin_latitude_deg) ||
      (Param.theRegionLong < myctl->region_origin_longitude_deg) ||
      (Param.theRegionDepth < myctl->region_depth_shallow_m) ||
      (Param.region_depth_deep_m > myctl->region_depth_deep_m) ||
      (Param.theRegionLat + Param.theDomainX / DIST1LAT >
       myctl->region_origin_latitude_deg +
           myctl->region_length_north_m / DIST1LAT) ||
      (Param.theRegionLong + Param.theDomainY / DIST1LON >
       myctl->region_origin_longitude_deg +
           myctl->region_length_east_m / DIST1LON)) {
    fprintf(stderr, "Mesh area out of the CVM etree\n");
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  /* Compute the coordinates of the origin of the mesh coordinate
     system in the CVM etree domain coordinate system */
  Global.theXForMeshOrigin =
      (Param.theRegionLat - myctl->region_origin_latitude_deg) * DIST1LAT;
  Global.theYForMeshOrigin =
      (Param.theRegionLong - myctl->region_origin_longitude_deg) * DIST1LON;
  Global.theZForMeshOrigin =
      Param.theRegionDepth - myctl->region_depth_shallow_m;

  /* Free memory used by myctl */
  cvm_freedbctl(myctl);

  double double_message_extra[3];

  double_message_extra[0] = Global.theXForMeshOrigin;
  double_message_extra[1] = Global.theYForMeshOrigin;
  double_message_extra[2] = Global.theZForMeshOrigin;

  MPI_Bcast(double_message_extra, 3, MPI_DOUBLE, 0, comm_solver);

  Global.theXForMeshOrigin = double_message_extra[0];
  Global.theYForMeshOrigin = double_message_extra[1];
  Global.theZForMeshOrigin = double_message_extra[2];

#else
  strcpy(Param.theCVMFlatFile, Param.cvmdb_input_file);
#endif
}

#ifdef USECVMDB

/**
 * replicateDB: Copy the material database to local disks.
 *
 */
static void replicateDB(const char* dbname) {
  char* destdir;

#ifndef SCEC
  char* srcpath;
  MPI_Comm replica_comm;
#endif /* SCEC */

  /* Change the working directory to $LOCAL */
#ifndef CVM_DESTDIR
  char curdir[256];
  destdir = getenv("CVM_DESTDIR");
  if (destdir == NULL) { /* then use current directory */
    destdir = getcwd(curdir, 256);
  }
#else
  destdir = CVM_DESTDIR;
#endif

  /* Clean the disks:
   * NOTE: Guys! cleanup the disk in your job script before launching
   * psolve, using rm -rf.
   * E.g., on Lemieux add the following line to your PBS job script,
   * before launching psolve.
   *   prun -m cyclic -n ${RMS_NODES} -N ${RMS_NODES} rm -rf <dirname>/
   * where dirname is the directory you want to wipe out.
   * This will take care of this issue.
   *
   * On BigBen it is even easier, it can be done from the front end
   * since everything is a shared parallel file system.
   *
   * Unfortunately the 'system' function is not supported on all platforms,
   * e.g., Cray's XT3 catamount platform (BigBen) does not have it.
   */
#ifndef SCEC
  if (chdir(destdir) != 0) {
    fprintf(stderr, "Thread %d: replicateDB: cannot chdir to %s\n", Global.myID,
            destdir);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  MPI_Barrier(comm_solver);

  /* Replicate the material database among the processors */
  if (Global.myID % PROCPERNODE != 0) {
    MPI_Comm_split(comm_solver, MPI_UNDEFINED, Global.myID, &replica_comm);
  } else {
    int replica_id;
    off_t filesize, remains, batchsize;
    void* filebuf;
    int src_fd = -1, dest_fd;

    MPI_Comm_split(comm_solver, 0, Global.myID, &replica_comm);
    MPI_Comm_rank(replica_comm, &replica_id);

    if (replica_id == 0) {
      struct stat statbuf;

#ifndef CVM_SRCPATH
      srcpath = getenv("CVM_SRCPATH");
#else
      srcpath = CVM_SRCPATH;
#endif

      if (stat(srcpath, &statbuf) != 0) {
        fprintf(stderr, "Thread 0: replicateDB: Cannot get stat of %s\n",
                srcpath);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }

      filesize = statbuf.st_size;
      src_fd = open(srcpath, O_RDONLY);
      if (src_fd == -1) {
        fprintf(stderr, "Thread 0: replicateDB: Cannot open cvm source db\n");
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }
    }

    MPI_Bcast(&filesize, sizeof(off_t), MPI_CHAR, 0, replica_comm);
    // theDBSize = filesize; // this is not used anywhere and will make building
    // process fail

    if ((filebuf = malloc(FILEBUFSIZE)) == NULL) {
      fprintf(stderr, "Thread %d: replicateDB: ", Global.myID);
      fprintf(stderr, "run out of memory while ");
      fprintf(stderr, "preparing to receive material database\n");
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

    /* Everyone opens a replicate db */
    dest_fd = open(dbname, O_CREAT | O_TRUNC | O_WRONLY, S_IRUSR | S_IWUSR);
    if (dest_fd == -1) {
      fprintf(stderr, "Thread %d: replicateDB: ", Global.myID);
      fprintf(stderr, "cannot create replica database\n");
      perror("open");
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

    remains = filesize;
    while (remains > 0) {
      batchsize = (remains > FILEBUFSIZE) ? FILEBUFSIZE : remains;

      if (replica_id == 0) {
        if (read(src_fd, filebuf, batchsize) != batchsize) {
          fprintf(stderr, "Thread 0: replicateDB: ");
          fprintf(stderr, "Cannot read database\n");
          perror("read");
          MPI_Abort(MPI_COMM_WORLD, ERROR);
          exit(1);
        }
      }

      MPI_Bcast(filebuf, batchsize, MPI_CHAR, 0, replica_comm);

      if (write(dest_fd, filebuf, batchsize) != batchsize) {
        fprintf(stderr, "Thread %d: replicateDB: ", Global.myID);
        fprintf(stderr, "Cannot write replica database\n");
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }

      remains -= batchsize;
    }

    free(filebuf);

    if (close(dest_fd) != 0) {
      fprintf(stderr, "Thread %d: replicateDB: ", Global.myID);
      fprintf(stderr, "cannot close replica database\n");
      perror("close");
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

    if (replica_id == 0) {
      close(src_fd);
    }

    MPI_Comm_free(&replica_comm);

  } /* processors participating in the replication */

#endif /* SCEC */

  return;
}

/*
 * Originally developed by Wenyang and Doriam, refactored into this function by
 * Clifford
 */
int getIstanbulMaterial(cvmpayload_t* g_props, vector3D_t IstVelModel_origin,
                        double x_m, double y_m, double z_m) {
  double output[4] = {0.0}, x_rel, y_rel;

  x_rel = 4536400.00 + x_m - IstVelModel_origin.x[0];
  y_rel = 388500.00 + y_m - IstVelModel_origin.x[1];

  int res = material_property_relative_V10_local(y_rel, x_rel, -z_m, output,
                                                 IstVelModel_origin.x[1],
                                                 IstVelModel_origin.x[0]);
  g_props->Qs = output[0] * 0.1;  // Same simple expression as in la Habra runs
  g_props->Qp = 2.0 * g_props->Qs;

  if (res != 0) {
    if (Param.useProfile == NO) {
      res = cvm_query(Global.theCVMEp, y_m, x_m, z_m, g_props);
    } else {
      res = profile_query(z_m, g_props);
    }
  } else {
    g_props->Vs = output[0];
    g_props->Vp = output[1];
    g_props->rho = output[2];
  }

  return res;
}

// Function to compute the interpolated z value based on the nearest four points
double interpolateZ(double** data, int x_count, double increment, double x_t,
                    double y_t) {
  // Find the index of the nearest point in the 'data' array based on x_t and
  // y_t
  int x1_index = (int)(x_t / increment);
  int x2_index = x1_index + 1;
  int y1_index = (int)(y_t / increment);
  int y2_index = y1_index + 1;
  // The four nearest points are labelled from the bottom left corner in a
  // counter-clockwise direction
  int nearest_point_index[4];
  nearest_point_index[0] = (y1_index * x_count) + x1_index;
  nearest_point_index[1] = (y1_index * x_count) + x2_index;
  nearest_point_index[2] = (y2_index * x_count) + x2_index;
  nearest_point_index[3] = (y2_index * x_count) + x1_index;
  // Corresponding x, y, and z
  double x1 = data[nearest_point_index[0]][0];
  double y1 = data[nearest_point_index[0]][1];
  double z[4];
  for (int i = 0; i < 4; i++) {
    z[i] = data[nearest_point_index[i]][2];
  }
  // Linear interpolation
  double xi = (x_t - x1) / increment * 2.0 - 1.0;
  double eta = (y_t - y1) / increment * 2.0 - 1.0;
  double N1 = 0.25 * (1.0 - xi) * (1.0 - eta);
  double N2 = 0.25 * (1.0 + xi) * (1.0 - eta);
  double N3 = 0.25 * (1.0 + xi) * (1.0 + eta);
  double N4 = 0.25 * (1.0 - xi) * (1.0 + eta);
  return N1 * z[0] + N2 * z[1] + N3 * z[2] + N4 * z[3];
}

// Function to determine if a point is in the basin
int isInBasin(double** data, int x_count, int y_count, double increment,
              double x_t, double y_t, double z_t) {
  if (x_t < 0 || x_t > (x_count - 1) * increment || y_t < 0 ||
      y_t > (y_count - 1) * increment) {
    return -1;  // Outside the basin
  }
  double interpolatedZ = interpolateZ(data, x_count, increment, x_t, y_t);
  // Compare z_t with the z value of the nearest point to determine if the point
  // is above the surface
  if (z_t < interpolatedZ) {
    return 0;  // Above the surface (in the basin)
  } else {
    return -1;  // Below or on the surface
  }
}

/*
 * Implemented by Clifford
 */
int getBasinMaterial(cvmpayload_t* g_props, vector3D_t basinVelModel_origin,
                     double x_m, double y_m, double z_m) {
  int res = 0;
  double x_rel = x_m - basinVelModel_origin.x[0];
  double y_rel = y_m - basinVelModel_origin.x[1];
  if (isInBasin(Param.basinData, Param.basinXCount, Param.basinYCount,
                Param.basinIncrement, x_rel, y_rel, z_m) == 0) {
    getBasinMaterialProperties(g_props, z_m,
                               (basin_material_type_t)Param.basinMaterialType,
                               Param.basinWaterVp, Param.basinWaterRho);
  } else {
    if (Param.useProfile == NO) {
      res = cvm_query(Global.theCVMEp, y_m, x_m, z_m, g_props);
    } else {
      res = profile_query(z_m, g_props);
    }
  }
  return res;
}

/*
 * Implemented by Clifford
 */
int get3DMaterial(cvmpayload_t* g_props, vector3D_t threeDVelModel_origin,
                  double x_m, double y_m, double z_m) {
  double output[3] = {0.0};
  double x_rel = x_m - threeDVelModel_origin.x[0];
  double y_rel = y_m - threeDVelModel_origin.x[1];
  int res = getMaterialFrom3DVelocityModel(x_rel, y_rel, z_m, output,
                                           threeDVelModel_origin.x[0],
                                           threeDVelModel_origin.x[1]);
  if (res != 0) {
    if (Param.useProfile == NO) {
      res = cvm_query(Global.theCVMEp, y_m, x_m, z_m, g_props);
    } else {
      res = profile_query(z_m, g_props);
    }
  } else {
    g_props->Vs = output[0];
    g_props->Vp = output[1];
    g_props->rho = output[2];

    if ((Param.useProfile == YES) &&
        !material_values_are_acoustic(g_props->Vp, g_props->Vs, g_props->rho) &&
        (g_props->Vs <= Param.theVsCut)) {
      g_props->Vs = Param.theVsCut;
    }

    g_props->Qs =
        g_props->Vs * 0.1;  // Same simple expression as in la Habra runs
    g_props->Qp = 2.0 * g_props->Qs;
  }
  return res;
}

static int material_values_are_elastic(double Vp, double Vs, double rho) {
  return (Vp > 0.0) && (rho > 0.0) && (Vs > ACOUSTIC_VS_EPS);
}

static int mesh_acoustic_face_axis(uint8_t face, int* axis, int* sign,
                                   int tangent_axes[2]) {
  if ((axis == NULL) || (sign == NULL) || (tangent_axes == NULL)) {
    return 0;
  }

  switch (face) {
    case ACOUSTIC_FACE_X_MINUS:
      *axis = 0;
      *sign = -1;
      tangent_axes[0] = 1;
      tangent_axes[1] = 2;
      return 1;
    case ACOUSTIC_FACE_X_PLUS:
      *axis = 0;
      *sign = 1;
      tangent_axes[0] = 1;
      tangent_axes[1] = 2;
      return 1;
    case ACOUSTIC_FACE_Y_MINUS:
      *axis = 1;
      *sign = -1;
      tangent_axes[0] = 0;
      tangent_axes[1] = 2;
      return 1;
    case ACOUSTIC_FACE_Y_PLUS:
      *axis = 1;
      *sign = 1;
      tangent_axes[0] = 0;
      tangent_axes[1] = 2;
      return 1;
    case ACOUSTIC_FACE_Z_MINUS:
      *axis = 2;
      *sign = -1;
      tangent_axes[0] = 0;
      tangent_axes[1] = 1;
      return 1;
    case ACOUSTIC_FACE_Z_PLUS:
      *axis = 2;
      *sign = 1;
      tangent_axes[0] = 0;
      tangent_axes[1] = 1;
      return 1;
    default:
      return 0;
  }
}

static void mesh_acoustic_face_normal(uint8_t face, double normal[3]) {
  normal[0] = 0.0;
  normal[1] = 0.0;
  normal[2] = 0.0;

  switch (face) {
    case ACOUSTIC_FACE_X_MINUS:
      normal[0] = -1.0;
      break;
    case ACOUSTIC_FACE_X_PLUS:
      normal[0] = 1.0;
      break;
    case ACOUSTIC_FACE_Y_MINUS:
      normal[1] = -1.0;
      break;
    case ACOUSTIC_FACE_Y_PLUS:
      normal[1] = 1.0;
      break;
    case ACOUSTIC_FACE_Z_MINUS:
      normal[2] = -1.0;
      break;
    case ACOUSTIC_FACE_Z_PLUS:
      normal[2] = 1.0;
      break;
  }
}

static void mesh_release_acoustic_interface_data(void) {
  free(Global.acousticFaceFlags);
  free(Global.acousticInterfaceElementIds);
  free(Global.acousticInterfaceFaceIds);
  free(Global.acousticInterfaceFaceAreas);
  free(Global.acousticInterfaceFaceCenters);
  free(Global.acousticInterfaceFaceNormals);
  free(Global.acousticInterfaceAdjacentElementIndexes);
  free(Global.acousticInterfaceAdjacentCenters);
  free(Global.acousticInterfaceAdjacentElementIds);
  Global.acousticFaceFlags = NULL;
  Global.acousticInterfaceElementIds = NULL;
  Global.acousticInterfaceFaceIds = NULL;
  Global.acousticInterfaceFaceAreas = NULL;
  Global.acousticInterfaceFaceCenters = NULL;
  Global.acousticInterfaceFaceNormals = NULL;
  Global.acousticInterfaceAdjacentElementIndexes = NULL;
  Global.acousticInterfaceAdjacentCenters = NULL;
  Global.acousticInterfaceAdjacentElementIds = NULL;
  Global.acousticInterfaceElements = 0;
  Global.acousticInterfaceFaces = 0;
}

static int mesh_find_local_element_at_domain_point(const vector3D_t* point,
                                                   int32_t* eindex_out) {
  octant_t* octant;
  tick_t xTick, yTick, zTick;
  int32_t eindex = -1;

  if ((Global.myMesh == NULL) || (Global.myOctree == NULL) || (point == NULL) ||
      (eindex_out == NULL)) {
    return 0;
  }

  if ((point->x[0] < 0.0) || (point->x[0] >= Param.theDomainX) ||
      (point->x[1] < 0.0) || (point->x[1] >= Param.theDomainY) ||
      (point->x[2] < 0.0) || (point->x[2] >= Param.theDomainZ)) {
    return 0;
  }

  xTick = (tick_t)(point->x[0] / Global.myMesh->ticksize);
  yTick = (tick_t)(point->x[1] / Global.myMesh->ticksize);
  zTick = (tick_t)(point->x[2] / Global.myMesh->ticksize);

  octant = octor_searchoctant(Global.myOctree, xTick, yTick, zTick, PIXELLEVEL,
                              AGGREGATE_SEARCH);
  if ((octant == NULL) || (octant->where == REMOTE)) {
    return 0;
  }

  if (octant->appdata != NULL) {
    int32_t cached = (int32_t)((intptr_t)octant->appdata - 1);
    if ((cached >= 0) && (cached < Global.myMesh->lenum)) {
      eindex = cached;
    }
  }

  if (eindex < 0) {
    for (eindex = 0; eindex < Global.myMesh->lenum; eindex++) {
      int32_t lnid0 = Global.myMesh->elemTable[eindex].lnid[0];

      if ((Global.myMesh->nodeTable[lnid0].x == octant->lx) &&
          (Global.myMesh->nodeTable[lnid0].y == octant->ly) &&
          (Global.myMesh->nodeTable[lnid0].z == octant->lz) &&
          (Global.myMesh->elemTable[eindex].level == octant->level)) {
        break;
      }
    }
  }

  if ((eindex < 0) || (eindex >= Global.myMesh->lenum)) {
    return 0;
  }

  *eindex_out = eindex;
  return 1;
}

static int mesh_query_material_at_model_point(
    double north_m, double east_m, double depth_m,
    const vector3D_t* IstVelModel_origin, const vector3D_t* basinVelModel_origin,
    const vector3D_t* threeDVelModel_origin, cvmpayload_t* g_props);

static int mesh_probe_face_elastic_info(
    const vector3D_t* domain_point, const vector3D_t* IstVelModel_origin,
    const vector3D_t* basinVelModel_origin,
    const vector3D_t* threeDVelModel_origin, int32_t* neighbor_eindex_out,
    int32_t* neighbor_level_out, int64_t* neighbor_geid_out) {
  int32_t neighbor_eindex = -1;
  int found_local_neighbor = 0;
  int local_neighbor_is_elastic = 0;
  cvmpayload_t g_props;
  double north_m;
  double east_m;
  double depth_m;

  if (neighbor_eindex_out != NULL) {
    *neighbor_eindex_out = -1;
  }
  if (neighbor_level_out != NULL) {
    *neighbor_level_out = -1;
  }
  if (neighbor_geid_out != NULL) {
    *neighbor_geid_out = -1;
  }

  if (domain_point == NULL) {
    return 0;
  }

  if ((domain_point->x[0] < 0.0) || (domain_point->x[0] >= Param.theDomainX) ||
      (domain_point->x[1] < 0.0) || (domain_point->x[1] >= Param.theDomainY) ||
      (domain_point->x[2] < 0.0) || (domain_point->x[2] >= Param.theDomainZ)) {
    return 0;
  }

  found_local_neighbor =
      mesh_find_local_element_at_domain_point(domain_point, &neighbor_eindex);
  if (found_local_neighbor) {
    elem_t* neighbor_elemp = &Global.myMesh->elemTable[neighbor_eindex];
    edata_t* neighbor_edata = (edata_t*)neighbor_elemp->data;

    local_neighbor_is_elastic = material_values_are_elastic(
        neighbor_edata->Vp, neighbor_edata->Vs, neighbor_edata->rho);

    if (!local_neighbor_is_elastic) {
      return 0;
    }

    if (neighbor_eindex_out != NULL) {
      *neighbor_eindex_out = neighbor_eindex;
    }
    if (neighbor_level_out != NULL) {
      *neighbor_level_out = neighbor_elemp->level;
    }
    if (neighbor_geid_out != NULL) {
      *neighbor_geid_out = neighbor_elemp->geid;
    }
    return 1;
  }

  north_m = Global.theXForMeshOrigin + domain_point->x[0];
  east_m = Global.theYForMeshOrigin + domain_point->x[1];
  depth_m = Global.theZForMeshOrigin + domain_point->x[2];

  if (mesh_query_material_at_model_point(north_m, east_m, depth_m,
                                         IstVelModel_origin, basinVelModel_origin,
                                         threeDVelModel_origin, &g_props) == 0) {
    return material_values_are_elastic(g_props.Vp, g_props.Vs, g_props.rho);
  }

  return 0;
}

static int mesh_query_material_at_model_point(
    double north_m, double east_m, double depth_m,
    const vector3D_t* IstVelModel_origin, const vector3D_t* basinVelModel_origin,
    const vector3D_t* threeDVelModel_origin, cvmpayload_t* g_props) {
  if (g_props == NULL) {
    return 1;
  }

  if (Param.includeBuildings == YES) {
    depth_m -= get_surface_shift();
  }

  if ((Param.includeTopography == YES) && (get_theetree_type() == SQD)) {
    depth_m -= point_elevation(north_m, east_m);
    if (depth_m < 0.0) {
      return 1;
    }
    if (depth_m > Param.theDomainZ) {
      depth_m = Param.theDomainZ;
    }
  }

  if (belongs2hmgHalfspace(east_m, north_m, depth_m)) {
    return get_halfspaceproperties(g_props);
  }

  if ((Param.IstanbulVelocityModel == YES) && (IstVelModel_origin != NULL)) {
    return getIstanbulMaterial(g_props, *IstVelModel_origin, north_m, east_m,
                               depth_m);
  }
  if ((Param.basinVelocityModel == YES) && (basinVelModel_origin != NULL)) {
    return getBasinMaterial(g_props, *basinVelModel_origin, north_m, east_m,
                            depth_m);
  }
  if ((Param.threeDVelocityModel == YES) && (threeDVelModel_origin != NULL)) {
    return get3DMaterial(g_props, *threeDVelModel_origin, north_m, east_m,
                         depth_m);
  }

  if (Param.useProfile == YES) {
    return profile_query(depth_m, g_props);
  }

  if (Global.theCVMEp != NULL) {
    return cvm_query(Global.theCVMEp, east_m, north_m, depth_m, g_props);
  }

  return 1;
}

static void mesh_reserve_acoustic_interface_fragments(int64_t needed,
                                                      int64_t current_count,
                                                      int64_t* capacity) {
  int64_t new_capacity;
  int32_t* new_elements;
  uint8_t* new_faces;
  double* new_areas;
  double* new_centers;
  double* new_normals;
  int32_t* new_adjacent_indexes;
  double* new_adjacent_centers;
  int64_t* new_adjacent;

  if ((capacity == NULL) || (needed <= *capacity)) {
    return;
  }
  if (needed > INT32_MAX) {
    fprintf(stderr,
            "Thread %d: too many local acoustic interface fragments for GPU table\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  new_capacity = (*capacity == 0) ? 1024 : *capacity;
  while (new_capacity < needed) {
    new_capacity *= 2;
  }

  new_elements = (int32_t*)malloc((size_t)new_capacity * sizeof(int32_t));
  new_faces = (uint8_t*)malloc((size_t)new_capacity * sizeof(uint8_t));
  new_areas = (double*)malloc((size_t)new_capacity * sizeof(double));
  new_centers = (double*)malloc(3 * (size_t)new_capacity * sizeof(double));
  new_normals = (double*)malloc(3 * (size_t)new_capacity * sizeof(double));
  new_adjacent_indexes =
      (int32_t*)malloc((size_t)new_capacity * sizeof(int32_t));
  new_adjacent_centers =
      (double*)malloc(3 * (size_t)new_capacity * sizeof(double));
  new_adjacent = (int64_t*)malloc((size_t)new_capacity * sizeof(int64_t));
  if ((new_elements == NULL) || (new_faces == NULL) || (new_areas == NULL) ||
      (new_centers == NULL) || (new_normals == NULL) ||
      (new_adjacent_indexes == NULL) || (new_adjacent_centers == NULL) ||
      (new_adjacent == NULL)) {
    free(new_elements);
    free(new_faces);
    free(new_areas);
    free(new_centers);
    free(new_normals);
    free(new_adjacent_indexes);
    free(new_adjacent_centers);
    free(new_adjacent);
    fprintf(stderr,
            "Thread %d: failed to allocate acoustic interface fragment table\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  if (current_count > 0) {
    size_t count = (size_t)current_count;
    memcpy(new_elements, Global.acousticInterfaceElementIds,
           count * sizeof(int32_t));
    memcpy(new_faces, Global.acousticInterfaceFaceIds, count * sizeof(uint8_t));
    memcpy(new_areas, Global.acousticInterfaceFaceAreas, count * sizeof(double));
    memcpy(new_centers, Global.acousticInterfaceFaceCenters,
           3 * count * sizeof(double));
    memcpy(new_normals, Global.acousticInterfaceFaceNormals,
           3 * count * sizeof(double));
    memcpy(new_adjacent_indexes, Global.acousticInterfaceAdjacentElementIndexes,
           count * sizeof(int32_t));
    memcpy(new_adjacent_centers, Global.acousticInterfaceAdjacentCenters,
           3 * count * sizeof(double));
    memcpy(new_adjacent, Global.acousticInterfaceAdjacentElementIds,
           count * sizeof(int64_t));
  }

  free(Global.acousticInterfaceElementIds);
  free(Global.acousticInterfaceFaceIds);
  free(Global.acousticInterfaceFaceAreas);
  free(Global.acousticInterfaceFaceCenters);
  free(Global.acousticInterfaceFaceNormals);
  free(Global.acousticInterfaceAdjacentElementIndexes);
  free(Global.acousticInterfaceAdjacentCenters);
  free(Global.acousticInterfaceAdjacentElementIds);

  Global.acousticInterfaceElementIds = new_elements;
  Global.acousticInterfaceFaceIds = new_faces;
  Global.acousticInterfaceFaceAreas = new_areas;
  Global.acousticInterfaceFaceCenters = new_centers;
  Global.acousticInterfaceFaceNormals = new_normals;
  Global.acousticInterfaceAdjacentElementIndexes = new_adjacent_indexes;
  Global.acousticInterfaceAdjacentCenters = new_adjacent_centers;
  Global.acousticInterfaceAdjacentElementIds = new_adjacent;
  *capacity = new_capacity;
}

static void mesh_append_acoustic_interface_fragment(
    int64_t* count, int64_t* capacity, int32_t eindex, uint8_t face,
    double area, const double center[3], const double normal[3],
    int32_t adjacent_eindex, int64_t adjacent_geid,
    const double adjacent_center[3]) {
  int64_t index;

  mesh_reserve_acoustic_interface_fragments(*count + 1, *count, capacity);
  index = *count;

  Global.acousticInterfaceElementIds[index] = eindex;
  Global.acousticInterfaceFaceIds[index] = face;
  Global.acousticInterfaceFaceAreas[index] = area;
  Global.acousticInterfaceFaceCenters[3 * index + 0] = center[0];
  Global.acousticInterfaceFaceCenters[3 * index + 1] = center[1];
  Global.acousticInterfaceFaceCenters[3 * index + 2] = center[2];
  Global.acousticInterfaceFaceNormals[3 * index + 0] = normal[0];
  Global.acousticInterfaceFaceNormals[3 * index + 1] = normal[1];
  Global.acousticInterfaceFaceNormals[3 * index + 2] = normal[2];
  Global.acousticInterfaceAdjacentElementIndexes[index] = adjacent_eindex;
  Global.acousticInterfaceAdjacentElementIds[index] = adjacent_geid;
  Global.acousticInterfaceAdjacentCenters[3 * index + 0] = adjacent_center[0];
  Global.acousticInterfaceAdjacentCenters[3 * index + 1] = adjacent_center[1];
  Global.acousticInterfaceAdjacentCenters[3 * index + 2] = adjacent_center[2];
  *count = index + 1;
}

static const char* mesh_acoustic_face_name(uint8_t face) {
  switch (face) {
    case ACOUSTIC_FACE_X_MINUS:
      return "x-";
    case ACOUSTIC_FACE_X_PLUS:
      return "x+";
    case ACOUSTIC_FACE_Y_MINUS:
      return "y-";
    case ACOUSTIC_FACE_Y_PLUS:
      return "y+";
    case ACOUSTIC_FACE_Z_MINUS:
      return "z-";
    case ACOUSTIC_FACE_Z_PLUS:
      return "z+";
    default:
      return "unknown";
  }
}

static void mesh_write_acoustic_interface_diagnostics(void) {
  const char* path = Param.theAcousticInterfaceDiagnosticsFile;
  int64_t local_count;
  int64_t global_count = 0;
  int64_t fragment_offset = 0;

  if ((path[0] == '\0') || !acoustic_pressure_enabled() ||
      (Global.myMesh == NULL)) {
    return;
  }

  local_count = Global.acousticInterfaceFaces;
  MPI_Allreduce(&local_count, &global_count, 1, MPI_INT64_T, MPI_SUM,
                comm_solver);

  if (Global.theGroupSize > 1) {
    MPI_Exscan(&local_count, &fragment_offset, 1, MPI_INT64_T, MPI_SUM,
               comm_solver);
    if (Global.myID == 0) {
      fragment_offset = 0;
    }
  }

  if (Global.myID == 0) {
    FILE* fp = hu_fopen(path, "w");
    fprintf(fp,
            "# fragment_id rank local_fragment acoustic_element_index "
            "acoustic_element_id acoustic_level face face_name "
            "adjacent_elastic_element_index adjacent_elastic_element_id "
            "area_m2 center_csi center_eta center_dzeta adjacent_csi "
            "adjacent_eta adjacent_dzeta center_domain_x_m center_domain_y_m "
            "center_domain_z_m center_model_north_m center_model_east_m "
            "center_model_depth_m normal_x normal_y normal_z h_m\n");
    fclose(fp);
  }

  MPI_Barrier(comm_solver);
  for (int rank = 0; rank < Global.theGroupSize; rank++) {
    if ((Global.myID == rank) && (local_count > 0)) {
      FILE* fp = hu_fopen(path, "a");

      for (int64_t i = 0; i < local_count; i++) {
        int32_t eindex = Global.acousticInterfaceElementIds[i];
        uint8_t face = Global.acousticInterfaceFaceIds[i];
        elem_t* elemp;
        edata_t* edata;
        int32_t lnid0;
        double h;
        double center_csi = Global.acousticInterfaceFaceCenters[3 * i + 0];
        double center_eta = Global.acousticInterfaceFaceCenters[3 * i + 1];
        double center_dzeta = Global.acousticInterfaceFaceCenters[3 * i + 2];
        double adjacent_csi =
            Global.acousticInterfaceAdjacentCenters[3 * i + 0];
        double adjacent_eta =
            Global.acousticInterfaceAdjacentCenters[3 * i + 1];
        double adjacent_dzeta =
            Global.acousticInterfaceAdjacentCenters[3 * i + 2];
        double normal_x = Global.acousticInterfaceFaceNormals[3 * i + 0];
        double normal_y = Global.acousticInterfaceFaceNormals[3 * i + 1];
        double normal_z = Global.acousticInterfaceFaceNormals[3 * i + 2];
        double lower_x, lower_y, lower_z;
        double center_x, center_y, center_z;

        if ((eindex < 0) || (eindex >= Global.myMesh->lenum)) {
          continue;
        }

        elemp = &Global.myMesh->elemTable[eindex];
        edata = (edata_t*)elemp->data;
        lnid0 = elemp->lnid[0];
        h = edata->edgesize;
        lower_x = Global.myMesh->ticksize * Global.myMesh->nodeTable[lnid0].x;
        lower_y = Global.myMesh->ticksize * Global.myMesh->nodeTable[lnid0].y;
        lower_z = Global.myMesh->ticksize * Global.myMesh->nodeTable[lnid0].z;
        center_x = lower_x + 0.5 * h * (center_csi + 1.0);
        center_y = lower_y + 0.5 * h * (center_eta + 1.0);
        center_z = lower_z + 0.5 * h * (center_dzeta + 1.0);

        fprintf(fp,
                "%" PRId64 " %d %" PRId64 " %d %" PRId64 " %d %u %s "
                "%d %" PRId64 " %.9e %.9e %.9e %.9e %.9e %.9e %.9e "
                "%.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e\n",
                fragment_offset + i, Global.myID, i, eindex, elemp->geid,
                elemp->level, (unsigned int)face,
                mesh_acoustic_face_name(face),
                Global.acousticInterfaceAdjacentElementIndexes[i],
                Global.acousticInterfaceAdjacentElementIds[i],
                Global.acousticInterfaceFaceAreas[i], center_csi, center_eta,
                center_dzeta, adjacent_csi, adjacent_eta, adjacent_dzeta,
                center_x, center_y, center_z,
                Global.theXForMeshOrigin + center_x,
                Global.theYForMeshOrigin + center_y,
                Global.theZForMeshOrigin + center_z, normal_x, normal_y,
                normal_z, h);
      }

      fclose(fp);
    }
    MPI_Barrier(comm_solver);
  }

  if (Global.myID == 0) {
    monitor_print(
        "Acoustic interface diagnostics: wrote %.1fM pressure fragments to %s\n",
        (double)global_count / 1.0e6, path);
  }
}

static int32_t mesh_interface_fragment_split(int32_t element_level,
                                             int32_t target_level) {
  int32_t level_delta;
  int32_t split = 1;

  if (target_level <= element_level) {
    return 1;
  }

  level_delta = target_level - element_level;
  if (level_delta > ACOUSTIC_INTERFACE_FRAGMENT_MAX_LEVEL_DELTA) {
    level_delta = ACOUSTIC_INTERFACE_FRAGMENT_MAX_LEVEL_DELTA;
  }
  split <<= level_delta;
  return split;
}

static int64_t mesh_build_acoustic_interface_face_fragments(
    int32_t eindex, uint8_t face, int append_fragments,
    int64_t* table_count, int64_t* table_capacity,
    const vector3D_t* IstVelModel_origin,
    const vector3D_t* basinVelModel_origin,
    const vector3D_t* threeDVelModel_origin) {
  elem_t* elemp = &Global.myMesh->elemTable[eindex];
  edata_t* edata = (edata_t*)elemp->data;
  int32_t lnid0 = elemp->lnid[0];
  int axis, sign, tangent_axes[2];
  int32_t target_level = elemp->level;
  int32_t split;
  int changed;
  int64_t fragment_count = 0;
  double h = edata->edgesize;
  double ldb[3] = {Global.myMesh->ticksize * Global.myMesh->nodeTable[lnid0].x,
                   Global.myMesh->ticksize * Global.myMesh->nodeTable[lnid0].y,
                   Global.myMesh->ticksize * Global.myMesh->nodeTable[lnid0].z};
  double eps = fmax(Global.myMesh->ticksize * 0.25, h * 1.0e-3);
  double normal[3];

  if (!mesh_acoustic_face_axis(face, &axis, &sign, tangent_axes)) {
    return 0;
  }
  mesh_acoustic_face_normal(face, normal);

  if ((Param.theAcousticInterfaceMaxLevel >= elemp->level) &&
      (Param.theAcousticInterfaceMaxLevel > target_level)) {
    target_level = Param.theAcousticInterfaceMaxLevel;
  }

  do {
    changed = 0;
    split = mesh_interface_fragment_split(elemp->level, target_level);

    for (int iu = 0; iu < split; iu++) {
      for (int iv = 0; iv < split; iv++) {
        vector3D_t probe;
        int32_t neighbor_level = -1;
        double u = ((double)iu + 0.5) / (double)split;
        double v = ((double)iv + 0.5) / (double)split;

        probe.x[0] = ldb[0] + 0.5 * h;
        probe.x[1] = ldb[1] + 0.5 * h;
        probe.x[2] = ldb[2] + 0.5 * h;
        probe.x[axis] = ldb[axis] + ((sign > 0) ? h : 0.0) + sign * eps;
        probe.x[tangent_axes[0]] = ldb[tangent_axes[0]] + u * h;
        probe.x[tangent_axes[1]] = ldb[tangent_axes[1]] + v * h;

        if (mesh_probe_face_elastic_info(
                &probe, IstVelModel_origin, basinVelModel_origin,
                threeDVelModel_origin, NULL, &neighbor_level, NULL) &&
            (neighbor_level > target_level)) {
          target_level = neighbor_level;
          changed = 1;
        }
      }
    }
  } while (changed);

  split = mesh_interface_fragment_split(elemp->level, target_level);
  for (int iu = 0; iu < split; iu++) {
    for (int iv = 0; iv < split; iv++) {
      vector3D_t fragment_point;
      vector3D_t probe;
      int32_t adjacent_eindex = -1;
      int64_t adjacent_geid = -1;
      double center[3] = {0.0, 0.0, 0.0};
      double adjacent_center[3] = {0.0, 0.0, 0.0};
      double u = ((double)iu + 0.5) / (double)split;
      double v = ((double)iv + 0.5) / (double)split;

      fragment_point.x[0] = ldb[0] + 0.5 * h;
      fragment_point.x[1] = ldb[1] + 0.5 * h;
      fragment_point.x[2] = ldb[2] + 0.5 * h;
      fragment_point.x[axis] = ldb[axis] + ((sign > 0) ? h : 0.0);
      fragment_point.x[tangent_axes[0]] = ldb[tangent_axes[0]] + u * h;
      fragment_point.x[tangent_axes[1]] = ldb[tangent_axes[1]] + v * h;

      probe = fragment_point;
      probe.x[axis] += sign * eps;

      if (!mesh_probe_face_elastic_info(
              &probe, IstVelModel_origin, basinVelModel_origin,
              threeDVelModel_origin, &adjacent_eindex, NULL, &adjacent_geid)) {
        continue;
      }

      center[axis] = (double)sign;
      center[tangent_axes[0]] = -1.0 + 2.0 * u;
      center[tangent_axes[1]] = -1.0 + 2.0 * v;
      fragment_count++;

      if (append_fragments) {
        const double fragment_h = h / (double)split;

        if ((adjacent_eindex >= 0) && (adjacent_eindex < Global.myMesh->lenum)) {
          elem_t* adjacent_elemp = &Global.myMesh->elemTable[adjacent_eindex];
          edata_t* adjacent_edata = (edata_t*)adjacent_elemp->data;
          int32_t adjacent_lnid0 = adjacent_elemp->lnid[0];
          double adjacent_h = adjacent_edata->edgesize;
          double adjacent_ldb[3] = {
              Global.myMesh->ticksize *
                  Global.myMesh->nodeTable[adjacent_lnid0].x,
              Global.myMesh->ticksize *
                  Global.myMesh->nodeTable[adjacent_lnid0].y,
              Global.myMesh->ticksize *
                  Global.myMesh->nodeTable[adjacent_lnid0].z};

          for (int dim = 0; dim < 3; dim++) {
            adjacent_center[dim] =
                2.0 * (fragment_point.x[dim] - adjacent_ldb[dim]) /
                    adjacent_h -
                1.0;
            if ((adjacent_center[dim] > 1.0) &&
                (adjacent_center[dim] < 1.0 + 1.0e-8)) {
              adjacent_center[dim] = 1.0;
            }
            if ((adjacent_center[dim] < -1.0) &&
                (adjacent_center[dim] > -1.0 - 1.0e-8)) {
              adjacent_center[dim] = -1.0;
            }
          }
        }

        mesh_append_acoustic_interface_fragment(
            table_count, table_capacity, eindex, face, fragment_h * fragment_h,
            center, normal, adjacent_eindex, adjacent_geid, adjacent_center);
      }
    }
  }

  return fragment_count;
}

static void mesh_build_acoustic_interface_flags(
    const vector3D_t* IstVelModel_origin, const vector3D_t* basinVelModel_origin,
    const vector3D_t* threeDVelModel_origin) {
  static const uint8_t face_bits[6] = {
      ACOUSTIC_FACE_X_MINUS, ACOUSTIC_FACE_X_PLUS, ACOUSTIC_FACE_Y_MINUS,
      ACOUSTIC_FACE_Y_PLUS,  ACOUSTIC_FACE_Z_MINUS, ACOUSTIC_FACE_Z_PLUS};
  int64_t local_counts[3] = {0, 0, 0};
  int64_t global_counts[3] = {0, 0, 0};
  int64_t table_count = 0;
  int64_t table_capacity = 0;

  mesh_release_acoustic_interface_data();

  if ((Param.includeAcoustic != YES) ||
      (Param.acousticInterfaceCoupling != YES) || (Global.myMesh == NULL)) {
    return;
  }

  Global.acousticFaceFlags =
      (uint8_t*)calloc((size_t)Global.myMesh->lenum, sizeof(uint8_t));
  if (Global.acousticFaceFlags == NULL) {
    fprintf(stderr,
            "Thread %d: failed to allocate acoustic interface face flags\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  for (int32_t eindex = 0; eindex < Global.myMesh->lenum; eindex++) {
    edata_t* edata = (edata_t*)Global.myMesh->elemTable[eindex].data;
    uint8_t flags = 0;

    if (!element_is_acoustic(edata)) {
      continue;
    }

    for (int iface = 0; iface < 6; iface++) {
      uint8_t face = face_bits[iface];
      int64_t fragments = mesh_build_acoustic_interface_face_fragments(
          eindex, face, acoustic_pressure_enabled(), &table_count,
          &table_capacity, IstVelModel_origin, basinVelModel_origin,
          threeDVelModel_origin);

      if (fragments == 0) {
        continue;
      }

      flags |= face;
      local_counts[1]++;
      local_counts[2] += fragments;
    }

    Global.acousticFaceFlags[eindex] = flags;

    if (flags != 0) {
      local_counts[0]++;
    }
  }

  Global.acousticInterfaceElements = local_counts[0];
  if (acoustic_pressure_enabled()) {
    assert(table_count == local_counts[2]);
  }
  Global.acousticInterfaceFaces =
      acoustic_pressure_enabled() ? table_count : local_counts[1];

  if (Global.acousticInterfaceFaces > INT32_MAX) {
    fprintf(stderr,
            "Thread %d: too many local acoustic interface fragments for GPU table\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  MPI_Reduce(local_counts, global_counts, 3, MPI_INT64_T, MPI_SUM, 0,
             comm_solver);
  if (Global.myID == 0) {
    if (acoustic_pressure_enabled()) {
      monitor_print(
          "Acoustic interface coupling: interface faces %.1fM, "
          "face fragments %.1fM, acoustic interface elements %.1fM\n",
          (double)global_counts[1] / 1.0e6,
          (double)global_counts[2] / 1.0e6,
          (double)global_counts[0] / 1.0e6);
    } else {
      monitor_print(
          "Acoustic interface coupling: normal-projected acoustic faces %.1fM, "
          "acoustic interface elements %.1fM\n",
          (double)global_counts[1] / 1.0e6,
          (double)global_counts[0] / 1.0e6);
    }
    fflush(stdout);
  }

  mesh_write_acoustic_interface_diagnostics();
}

/**
 * Assign values (material properties) to a leaf octant specified by
 * octleaf.  In order to refine the mesh, select the minimum Vs of 27
 * sample points: 8 near the corners and 19 midpoints.
 */
static void setrec(octant_t* leaf, double ticksize, void* data) {
  double x_m, y_m, z_m; /* x:south-north, y:east-west, z:depth */
  tick_t halfticks;
  cvmpayload_t g_props;     /* cvm record with ground properties */
  cvmpayload_t g_props_min; /* cvm record with the min Vs found */
  vector3D_t basinVelModel_origin;
  vector3D_t threeDVelModel_origin;

  int i_x, i_y, i_z, n_points = 3;
  int acoustic_samples = 0, elastic_samples = 0;
  double points[3];

  int res = 0;
  edata_t* edata = (edata_t*)data;

  points[0] = 0.01;
  points[1] = 1;
  points[2] = 1.99;

  halfticks = (tick_t)1 << (PIXELLEVEL - leaf->level - 1);
  edata->edgesize = ticksize * halfticks * 2;

  if (Param.IstanbulVelocityModel == YES && g_hasIstanbulVelModelOrigin == 0) {
    solver_abort(__FUNCTION_NAME, NULL,
                 "Istanbul origin not initialized before meshing\n");
  }
  if (Param.basinVelocityModel == YES) {
    basinVelModel_origin = compute_domain_coords_linearinterp(
        Param.theBasinLong, Param.theBasinLat, Param.theSurfaceCornersLong,
        Param.theSurfaceCornersLat, Param.theDomainY, Param.theDomainX,
        &Param.utmZone);
  }
  if (Param.threeDVelocityModel == YES) {
    threeDVelModel_origin = compute_domain_coords_linearinterp(
        Param.the3DVelocityModelLong, Param.the3DVelocityModelLat,
        Param.theSurfaceCornersLong, Param.theSurfaceCornersLat,
        Param.theDomainY, Param.theDomainX, &Param.utmZone);
  }

  /* Check for buildings and proceed according to the buildings setrec */
  if ((Param.includeBuildings == YES) && (Param.useProfile == NO)) {
    if (bldgs_setrec(leaf, ticksize, edata, Global.theCVMEp,
                     Global.theXForMeshOrigin, Global.theYForMeshOrigin,
                     Global.theZForMeshOrigin)) {
      return;
    }
  }

  if (Param.includeTopography == YES) {
    if (topo_setrec(leaf, ticksize, edata, Global.theCVMEp)) {
      return;
    }
  }

  g_props_min.Vs = FLT_MAX;
  g_props_min.Vp = NAN;
  g_props_min.rho = NAN;

  for (i_x = 0; i_x < n_points; i_x++) {
    x_m = (Global.theXForMeshOrigin +
           (leaf->lx + points[i_x] * halfticks) * ticksize);

    for (i_y = 0; i_y < n_points; i_y++) {
      y_m = Global.theYForMeshOrigin +
            (leaf->ly + points[i_y] * halfticks) * ticksize;

      for (i_z = 0; i_z < n_points; i_z++) {
        z_m = Global.theZForMeshOrigin +
              (leaf->lz + points[i_z] * halfticks) * ticksize;

        /* Aug 15, 2024 Update: This should be addressed by implementing PROJ.
         */
        // Clifford's NOTE: Try debugging. I realized that it's possible the
        // point is outside the domain. But it seems like it's normal.
        // oct_expand() will take care of the points outside the domain.
        // ===== Print the coordinates of the point if it's outside the domain
        // ===== if (x_m < 0 || x_m > Param.theDomainX || y_m < 0 || y_m >
        // Param.theDomainY || z_m < 0 || z_m > Param.theDomainZ) {
        //     printf("Point outside the domain: (%f, %f, %f)\n", x_m, y_m,
        //     z_m);
        // }

        /* Shift the domain if topography with squashed etree is considered */
        if ((Param.includeTopography == YES) && (get_theetree_type() == SQD)) {
          z_m -= point_elevation(x_m, y_m);
          if (z_m < 0) /* get air properties */
          {
            edata_t* airedata = (edata_t*)data;
            get_airprops_topo(airedata);
            g_props.Vp = airedata->Vp;
            g_props.Vs = airedata->Vs;
            g_props.rho = airedata->rho;

            if (g_props.Vs < g_props_min.Vs) {
              /* assign minimum value of vs to produce elements
               * that are small enough to rightly represent the model */
              g_props_min = g_props;
            }

            continue;
          }
        }

        /* Shift the domain if buildings are considered */
        if (Param.includeBuildings == YES) {
          z_m -= get_surface_shift();
        }

        if (belongs2hmgHalfspace(y_m, x_m, z_m))
          res = get_halfspaceproperties(&g_props);
        else if (Param.IstanbulVelocityModel == YES) {
          res = getIstanbulMaterial(&g_props, g_istanbulVelModelOrigin, x_m,
                                    y_m, z_m);
        } else if (Param.basinVelocityModel == YES) {
          res = getBasinMaterial(&g_props, basinVelModel_origin, x_m, y_m, z_m);
        } else if (Param.threeDVelocityModel == YES) {
          res = get3DMaterial(&g_props, threeDVelModel_origin, x_m, y_m, z_m);
        } else if (Param.useProfile == NO) {
          res = cvm_query(Global.theCVMEp, y_m, x_m, z_m, &g_props);
        } else {
          res = profile_query(z_m, &g_props);
        }

        if (res != 0) {
          continue;
        }

        if (material_values_are_acoustic(g_props.Vp, g_props.Vs,
                                         g_props.rho)) {
          acoustic_samples++;
        } else if (g_props.Vs > 0.0) {
          elastic_samples++;
        }

        if (g_props.Vs < g_props_min.Vs) {
          /* assign minimum value of vs to produce elements
           * that are small enough to rightly represent the model */
          g_props_min = g_props;
        }

        if ((Param.includeAcoustic != YES) && (g_props.Vs <= Param.theVsCut)) {
          /* stop early if needed, completely break out of all
           * the loops, the label is just outside the loop */
          goto outer_loop_label;
        }
      }
    }
  }
outer_loop_label: /* in order to completely break out from the inner loop */

  if ((Param.includeAcoustic == YES) && (acoustic_samples > 0) &&
      (elastic_samples > 0) && (Param.theAcousticInterfaceMaxLevel >= 0) &&
      (leaf->level < Param.theAcousticInterfaceMaxLevel)) {
    edata->Vp = isnan(g_props_min.Vp) ? Param.theVsCut : g_props_min.Vp;
    edata->rho = isnan(g_props_min.rho) ? 1.0f : g_props_min.rho;
    edata->Vs = Param.theFactor * edata->edgesize / 2.0;
    return;
  }

  edata->Vp = g_props_min.Vp;
  edata->Vs = g_props_min.Vs;
  edata->rho = g_props_min.rho;

  if (res != 0 && g_props_min.Vs == FLT_MAX) {
    /* all the queries failed, then center out of bound point. Set Vs
     * to force split */
    edata->Vs = Param.theFactor * edata->edgesize / 2;
  } else if (element_is_acoustic(edata)) {
    edata->Vs = 0.0f;
  } else if (edata->Vs <= 0.0f) {
    solver_abort(__FUNCTION_NAME, NULL,
                 "Encountered Vs <= 0 with include_acoustic = no. "
                 "Set include_acoustic = yes for water/acoustic regions.\n");
  } else if (edata->Vs <= Param.theVsCut) { /* adjust Vs, Vp, and rho */
    double VpVsRatio = edata->Vp / edata->Vs;
    double RhoVpRatio = edata->rho / edata->Vp;

    edata->Vs = Param.theVsCut;
    edata->Vp = Param.theVsCut * VpVsRatio;
    edata->rho = edata->Vp * RhoVpRatio;
  }

  return;
}

#else  /* USECVMDB */

static int32_t zsearch(void* base, int32_t count, int32_t recordsize,
                       const point_t* searchpt) {
  int32_t start, end, offset, found;

  start = 0;
  end = count - 1;
  offset = (start + end) / 2;

  found = 0;
  do {
    if (end < start) {
      /* the two pointer crossed each other */
      offset = end;
      found = 1;
    } else {
      const void* pivot = (char*)base + offset * recordsize;

      switch (octor_zcompare(searchpt, pivot)) {
        case (0): /* equal */
          found = 1;
          break;
        case (1): /* searchpoint larger than the pivot */
          start = offset + 1;
          offset = (start + end) / 2;
          break;
        case (-1): /* searchpoint smaller than the pivot */
          end = offset - 1;
          offset = (start + end) / 2;
          break;
      }
    }
  } while (!found);

  return offset;
}

static cvmrecord_t* sliceCVM(const char* cvm_flatfile) {
  cvmrecord_t* cvmrecord;
  int32_t bufferedbytes, bytecount, recordcount;
  if (Global.myID == Global.theGroupSize - 1) {
    /* the last processor reads data and
       distribute to other processors*/

    struct timeval starttime, endtime;
    float iotime = 0, memmovetime = 0;
    MPI_Request* isendreqs;
    MPI_Status* isendstats;
    FILE* fp;
    int fd, procid;
    struct stat statbuf;
    void* maxbuf;
    const point_t* intervaltable;
    off_t bytesent;
    int32_t offset;
    const int maxcount = (1 << 29) / sizeof(cvmrecord_t);
    const int maxbufsize = maxcount * sizeof(cvmrecord_t);

    fp = fopen(cvm_flatfile, "r");
    if (fp == NULL) {
      fprintf(stderr, "Thread %d: Cannot open flat CVM file\n", Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

    fd = fileno(fp);
    if (fstat(fd, &statbuf) != 0) {
      fprintf(stderr, "Thread %d: Cannot get the status of CVM file\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

    intervaltable = octor_getintervaltable(Global.myOctree);

    /*
    for (procid = 0; procid <= Global.myID; procid++) {
        fprintf(stderr, "interval[%d] = {%d, %d, %d}\n", procid,
                intervaltable[procid].x << 1, intervaltable[procid].y << 1,
                intervaltable[procid].z << 1);
    }
    */

    bytesent = 0;
    maxbuf = malloc(maxbufsize);
    if (maxbuf == NULL) {
      fprintf(stderr, "Thread %d: Cannot allocate send buffer\n", Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

    isendreqs = (MPI_Request*)malloc(sizeof(MPI_Request) * Global.theGroupSize);
    isendstats = (MPI_Status*)malloc(sizeof(MPI_Status) * Global.theGroupSize);
    if ((isendreqs == NULL) || (isendstats == NULL)) {
      fprintf(stderr, "Thread %d: Cannot allocate isend controls\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

    /* Try to read max number of CVM records as allowed */
    // gettimeofday(&starttime, NULL);
    recordcount = fread(maxbuf, sizeof(cvmrecord_t),
                        maxbufsize / sizeof(cvmrecord_t), fp);
    // gettimeofday(&endtime, NULL);

    iotime += (endtime.tv_sec - starttime.tv_sec) * 1000.0 +
              (endtime.tv_usec - starttime.tv_usec) / 1000.0;

    if (recordcount != maxbufsize / sizeof(cvmrecord_t)) {
      fprintf(stderr, "Thread %d: Cannot read-init buffer\n", Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

    /* start with proc 0 */
    procid = 0;

    while (procid < Global.myID) { /* repeatedly fill the buffer */
      point_t searchpoint, *point;
      int newreads;
      int isendcount = 0;

      /* we have recordcount to work with */
      cvmrecord = (cvmrecord_t*)maxbuf;

      while (procid < Global.myID) { /* repeatedly send out data */

        searchpoint.x = intervaltable[procid + 1].x << 1;
        searchpoint.y = intervaltable[procid + 1].y << 1;
        searchpoint.z = intervaltable[procid + 1].z << 1;

        offset = zsearch(cvmrecord, recordcount, Global.theCVMRecordSize,
                         &searchpoint);

        point = (point_t*)(cvmrecord + offset);

        if ((point->x != searchpoint.x) || (point->y != searchpoint.y) ||
            (point->z != searchpoint.z)) {
          break;
        } else {
          bytecount = offset * sizeof(cvmrecord_t);
          MPI_Isend(cvmrecord, bytecount, MPI_CHAR, procid, CVMRECORD_MSG,
                    comm_solver, &isendreqs[isendcount]);
          isendcount++;

          /*
            fprintf(stderr,
            "Procid = %d offset = %qd bytecount = %d\n",
            procid, (int64_t)bytesent, bytecount);
          */

          bytesent += bytecount;

          /* prepare for the next processor */
          recordcount -= offset;
          cvmrecord = (cvmrecord_t*)point;
          procid++;
        }
      }

      /* Wait till the data in the buffer has been sent */
      MPI_Waitall(isendcount, isendreqs, isendstats);

      /* Move residual data to the beginning of the buffer
         and try to fill the newly free space */
      bufferedbytes = sizeof(cvmrecord_t) * recordcount;

      // gettimeofday(&starttime, NULL);
      memmove(maxbuf, cvmrecord, bufferedbytes);
      // gettimeofday(&endtime, NULL);
      memmovetime += (endtime.tv_sec - starttime.tv_sec) * 1000.0 +
                     (endtime.tv_usec - starttime.tv_usec) / 1000.0;

      // gettimeofday(&starttime, NULL);
      newreads = fread((char*)maxbuf + bufferedbytes, sizeof(cvmrecord_t),
                       maxcount - recordcount, fp);
      // gettimeofday(&endtime, NULL);
      iotime += (endtime.tv_sec - starttime.tv_sec) * 1000.0 +
                (endtime.tv_usec - starttime.tv_usec) / 1000.0;

      recordcount += newreads;

      if (newreads == 0) break;
    }

    free(maxbuf);
    free(isendreqs);
    free(isendstats);

    /* I am supposed to accomodate the remaining octants */
    bytecount = statbuf.st_size - bytesent;

    cvmrecord = (cvmrecord_t*)malloc(bytecount);
    if (cvmrecord == NULL) {
      fprintf(stderr, "Thread %d: out of memory\n", Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

    /* fseek exiting the for loop has file cursor propertly */
    if (fseeko(fp, bytesent, SEEK_SET) != 0) {
      fprintf(stderr, "Thread %d: fseeko failed\n", Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

    //  gettimeofday(&starttime, NULL);
    if (fread(cvmrecord, 1, bytecount, fp) != (size_t)bytecount) {
      fprintf(stderr, "Thread %d: fail to read the last chunk\n", Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    //  gettimeofday(&endtime, NULL);
    iotime += (endtime.tv_sec - starttime.tv_sec) * 1000.0 +
              (endtime.tv_usec - starttime.tv_usec) / 1000.0;

    /*
    fprintf(stderr, "Procid = %d offset = %qd bytecount = %d\n",
            Global.myID, (int64_t)bytesent, bytecount);
    */

    fclose(fp);

    fprintf(stdout, "Read %s (%.2fMB) in %.2f seconds (%.2fMB/sec)\n",
            cvm_flatfile, (float)statbuf.st_size / (1 << 20), iotime / 1000,
            (float)statbuf.st_size / (1 << 20) / (iotime / 1000));

    monitor_print("Memmove takes %.2f seconds\n", (float)memmovetime / 1000);
  } else {
    /* wait for my turn till PE(n - 1) tells me to go ahead */

    MPI_Status status;

    MPI_Probe(Global.theGroupSize - 1, CVMRECORD_MSG, comm_solver, &status);
    MPI_Get_count(&status, MPI_CHAR, &bytecount);

    cvmrecord = (cvmrecord_t*)malloc(bytecount);
    if (cvmrecord == NULL) {
      fprintf(stderr, "Thread %d: out of memory\n", Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

    MPI_Recv(cvmrecord, bytecount, MPI_CHAR, Global.theGroupSize - 1,
             CVMRECORD_MSG, comm_solver, &status);
  }

  /* Every processor should set these parameters correctly */
  Global.theCVMRecordCount = bytecount / sizeof(cvmrecord_t);
  if (Global.theCVMRecordCount * sizeof(cvmrecord_t) != (size_t)bytecount) {
    fprintf(stderr, "Thread %d: received corrupted CVM data\n", Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  return cvmrecord;
}

/**
 * setrec: Search the CVM record array to obtain the material property of
 *     a leaf octant.
 *
 */
void setrec(octant_t* leaf, double ticksize, void* data) {
  cvmrecord_t* agghit;
  edata_t* edata;
  etree_tick_t x, y, z;
  etree_tick_t halfticks;
  point_t searchpoint;

  edata = (edata_t*)data;

  halfticks = (tick_t)1 << (PIXELLEVEL - leaf->level - 1);

  edata->edgesize = ticksize * halfticks * 2;

  searchpoint.x = x = leaf->lx + halfticks;
  searchpoint.y = y = leaf->ly + halfticks;
  searchpoint.z = z = leaf->lz + halfticks;

  if ((x * ticksize >= Param.theDomainX) ||
      (y * ticksize >= Param.theDomainY) ||
      (z * ticksize >= Param.theDomainZ)) {
    /* Center point out the bound. Set Vs to force split */
    edata->Vs = Param.theFactor * edata->edgesize / 2;
  } else {
    int offset;

    /* map the coordinate from the octor address space to the
       etree address space */
    searchpoint.x = x << 1;
    searchpoint.y = y << 1;
    searchpoint.z = z << 1;

    /* Inbound */
    offset = zsearch(Global.theCVMRecord, Global.theCVMRecordCount,
                     Global.theCVMRecordSize, &searchpoint);
    if (offset < 0) {
      fprintf(stderr, "setrec: fatal error\n");
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

    agghit = Global.theCVMRecord + offset;
    edata->Vs = agghit->Vs;
    edata->Vp = agghit->Vp;
    edata->rho = agghit->density;

    /* Adjust the Vs */
    if (element_is_acoustic(edata)) {
      edata->Vs = 0.0f;
    } else if (edata->Vs <= 0.0f) {
      solver_abort(__FUNCTION_NAME, NULL,
                   "Encountered Vs <= 0 with include_acoustic = no. "
                   "Set include_acoustic = yes for water/acoustic regions.\n");
    } else {
      edata->Vs = (edata->Vs < Param.theVsCut) ? Param.theVsCut : edata->Vs;
    }
  }

  return;
}
#endif /* USECVMDB */

/**
 * mesh_generate: Generate and partition an unstructured octree mesh.
 *
 */
static void mesh_generate() {
  int mstep, step = 1;
  double originalFactor = Param.theFactor;
  double ppwl = Param.theFactor / Param.theFreq;
  double prevtref = 0, prevtbal = 0, prevtpar = 0;
  int64_t tote, mine, maxe;

  if (Global.myID == 0) {
    monitor_print("Meshing: ");
    if (Param.theStepMeshingFactor == 0) {
      monitor_print("Conventional\n\n");
    } else {
      monitor_print("Progressive\n\n");
    }
    monitor_print("Stage %14s Min %7s Max %5s Total    Time(s)", "", "", "");
    if (Param.theStepMeshingFactor == 0) {
      monitor_print("\n\n");
    } else {
      monitor_print("   Step  f(Hz)\n\n");
    }
  }

  /*----  Generate and partition an unstructured octree mesh ----*/
  MPI_Barrier(comm_solver);
  Timer_Start("Octor Newtree");
  if (Global.myID == 0) {
    monitor_print("New tree %41s", "");
  }
  Global.myOctree = octor_newtree(
      Param.theDomainX, Param.theDomainY, Param.theDomainZ, sizeof(edata_t),
      Global.myID, Global.theGroupSize, comm_solver, get_surface_shift());

  /* NOTE:
   * If you want to see the carving process, replace by:
   *     Global.myOctree = octor_newtree(
   *             Param.theDomainX, Param.theDomainY, Param.theDomainZ,
   *             sizeof(edata_t), Global.myID, Global.theGroupSize,
   *             comm_solver, 0);
   */

  if (Global.myOctree == NULL) {
    monitor_print("Thread %d: mesh_generate: fail to create octree\n",
                  Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  MPI_Barrier(comm_solver);
  Timer_Stop("Octor Newtree");
  if (Global.myID == 0) {
    monitor_print("%9.2f\n\n", Timer_Value("Octor Newtree", (TimerKind)0));
  }

  /* Essential for DRM implementation */
  if (Param.drmImplement == YES) {
    drm_fix_coordinates(Global.myOctree->ticksize);
  }

#ifdef USECVMDB
  Global.theCVMQueryStage = 0; /* Query CVM database to refine the mesh */
#else
  /* Use flat data record file and distibute the data in memories */
  if (Global.myID == 0) {
    monitor_print("slicing CVMDB ...");
  }
  Timer_Start("Slice CVM");
  Global.theCVMRecord = sliceCVM(Param.theCVMFlatFile);
  MPI_Barrier(comm_solver);
  Timer_Stop("Slice CVM");
  if (Global.theCVMRecord == NULL) {
    fprintf(stderr, "Thread %d: Error obtaining the CVM records from %s\n",
            Global.myID, Param.theCVMFlatFile);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  };
  if (Global.myID == 0) {
    monitor_print("done : %9.2f seconds\n",
                  Timer_Value("Slice CVM", (TimerKind)0));
  }
#endif

  for (mstep = Param.theStepMeshingFactor; mstep >= 0; mstep--) {
    double myFactor = (double)(1 << mstep);  // 2^mstep
    Param.theFactor = originalFactor / myFactor;

    /* Refinement */
    Timer_Start("Octor Refinetree");
    if (Global.myID == 0) {
      monitor_print("Refining     ");
      fflush(stdout);
    }
    if (octor_refinetree(Global.myOctree, toexpand, setrec) != 0) {
      fprintf(stderr, "Thread %d: mesh_generate: fail to refine octree\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    MPI_Barrier(comm_solver);
    tote = octor_getleavescount(Global.myOctree, GLOBAL);
    mine = octor_getminleavescount(Global.myOctree, GLOBAL);
    maxe = octor_getmaxleavescount(Global.myOctree, GLOBAL);
    if (Global.myID == 0) {
      monitor_print("%11" INT64_FMT " %11" INT64_FMT " %11" INT64_FMT, mine,
                    maxe, tote);
      fflush(stdout);
    }
    Timer_Stop("Octor Refinetree");
    if (Global.myID == 0) {
      monitor_print("%11.2f",
                    Timer_Value("Octor Refinetree", (TimerKind)0) - prevtref);
      if (Param.theStepMeshingFactor == 0) {
        monitor_print("\n");
      } else {
        monitor_print("   %4d %6.2f\n", step, Param.theFactor / ppwl);
      }
      prevtref = Timer_Value("Octor Refinetree", (TimerKind)0);
      fflush(stdout);
    }

    /* Balancing */
    Timer_Start("Octor Balancetree");
    if (Global.myID == 0) {
      monitor_print("Balancing    ");
      fflush(stdout);
    }
    if (octor_balancetree(Global.myOctree, setrec,
                          Param.theStepMeshingFactor) != 0) {
      fprintf(stderr, "Thread %d: mesh_generate: fail to balance octree\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    MPI_Barrier(comm_solver);
    tote = octor_getleavescount(Global.myOctree, GLOBAL);
    mine = octor_getminleavescount(Global.myOctree, GLOBAL);
    maxe = octor_getmaxleavescount(Global.myOctree, GLOBAL);
    if (Global.myID == 0) {
      monitor_print("%11" INT64_FMT " %11" INT64_FMT " %11" INT64_FMT, mine,
                    maxe, tote);
      fflush(stdout);
    }
    Timer_Stop("Octor Balancetree");
    if (Global.myID == 0) {
      monitor_print("%11.2f\n",
                    Timer_Value("Octor Balancetree", (TimerKind)0) - prevtbal);
      prevtbal = Timer_Value("Octor Balancetree", (TimerKind)0);
      fflush(stdout);
    }

    /* Partitioning */
    Timer_Start("Octor Partitiontree");
    if (Global.myID == 0) {
      monitor_print("Partitioning ");
      fflush(stdout);
    }
    if (octor_partitiontree(Global.myOctree, bldgs_nodesearch_com,
                            pushdowns_nodesearch) != 0) {
      fprintf(stderr, "Thread %d: mesh_generate: fail to balance load\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    MPI_Barrier(comm_solver);
    tote = octor_getleavescount(Global.myOctree, GLOBAL);
    mine = octor_getminleavescount(Global.myOctree, GLOBAL);
    maxe = octor_getmaxleavescount(Global.myOctree, GLOBAL);
    if (Global.myID == 0) {
      monitor_print("%11" INT64_FMT " %11" INT64_FMT " %11" INT64_FMT, mine,
                    maxe, tote);
      fflush(stdout);
    }
    Timer_Stop("Octor Partitiontree");
    if (Global.myID == 0) {
      monitor_print(
          "%11.2f\n\n",
          Timer_Value("Octor Partitiontree", (TimerKind)0) - prevtpar);
      prevtpar = Timer_Value("Octor Partitiontree", (TimerKind)0);
      fflush(stdout);
    }

    step++;
    fflush(stdout);
    MPI_Barrier(comm_solver);
  }

  /* gggg */

  /* Buildings/Topography Carving */
  //    if ( Param.includeBuildings == YES ||  Param.includeTopography == YES )
  //    {
  if (Param.includeBuildings == YES) {
    int flag = 0;
    Timer_Start("Carve Buildings");
    if (Global.myID == 0) {
      monitor_print("Carving Buildings/Topography");
      fflush(stdout);
    }

    if (Param.includeBuildings == YES) flag = 1;

    /* NOTE: If you want to see the carving process, comment next line */

    octor_carvebuildings(Global.myOctree, flag, bldgs_nodesearch_com,
                         pushdowns_search);
    MPI_Barrier(comm_solver);
    Timer_Stop("Carve Buildings");
    if (Global.myID == 0) {
      monitor_print("%9.2f\n", Timer_Value("Carve Buildings", (TimerKind)0));
      fflush(stdout);
    }

    Timer_Start("Octor Partitiontree");
    if (Global.myID == 0) {
      monitor_print("Repartitioning");
      fflush(stdout);
    }
    if (octor_partitiontree(Global.myOctree, bldgs_nodesearch_com,
                            pushdowns_nodesearch) != 0) {
      fprintf(stderr, "Thread %d: mesh_generate: fail to balance load\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    MPI_Barrier(comm_solver);
    Timer_Stop("Octor Partitiontree");
    if (Global.myID == 0) {
      monitor_print("%9.2f\n",
                    Timer_Value("Octor Partitiontree", (TimerKind)0));
      fflush(stdout);
    }
  }

  if (Global.myID == 0 && Param.theStepMeshingFactor != 0) {
    monitor_print("Total refine    %33s %9.2f\n", "",
                  Timer_Value("Octor Refinetree", (TimerKind)0));
    monitor_print("Total balance   %33s %9.2f\n", "",
                  Timer_Value("Octor Balancetree", (TimerKind)0));
    monitor_print("Total partition %33s %9.2f\n\n", "",
                  Timer_Value("Octor Partitiontree", (TimerKind)0));
    fflush(stdout);
  }

  Timer_Start("Octor Extractmesh");
  if (Global.myID == 0) {
    monitor_print("Extracting the mesh %30s", "");
    fflush(stdout);
  }
  Global.myMesh =
      octor_extractmesh(Global.myOctree, bldgs_nodesearch, pushdowns_nodesearch,
                        bldgs_nodesearch_com, find_topoAirOct, topo_nodesearch);

  if (Global.myMesh == NULL) {
    fprintf(stderr, "Thread %d: mesh_generate: fail to extract mesh\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  MPI_Barrier(comm_solver);
  Timer_Stop("Octor Extractmesh");
  if (Global.myID == 0) {
    monitor_print("%9.2f\n", Timer_Value("Octor Partitiontree", (TimerKind)0));
  }

  Timer_Start("Mesh correct properties");
  /* Re-populates the mesh with actual values from the CVM-etree */
  if (Global.myID == 0) {
    monitor_print("Correcting mesh properties %23s", "");
    fflush(stdout);
  }

  mesh_correct_properties(Global.theCVMEp);

  MPI_Barrier(comm_solver);
  Timer_Stop("Mesh correct properties");
  if (Global.myID == 0) {
    monitor_print("%9.2f\n\n",
                  Timer_Value("Mesh correct properties", (TimerKind)0));
    fflush(stdout);
  }

#ifdef USECVMDB
  if (Param.useProfile == NO) {
    /* Close the material database */
    etree_close(Global.theCVMEp);
  }
#else
  free(Global.theCVMRecord);
#endif /* USECVMDB */
}

/**
 * toexpand: Instruct the Octor library whether a leaf octant needs to
 *       be expanded or not. Return 1 if true, 0 otherwise.
 *
 */
static int32_t toexpand(octant_t* leaf, double ticksize, const void* data) {
  if (data == NULL) {
    return 1;
  }

  int res;
  edata_t* edata = (edata_t*)data;

  if (Param.includeBuildings == YES) {
    res = bldgs_toexpand(leaf, ticksize, edata, Param.theFactor);
    if (res != -1) {
      return res;
    }
  }

  if (Param.drmImplement == YES) {
    // if( Param.drmImplement == YES && Param.theDrmPart != PART1 ) {
    res = drm_toexpand(leaf, ticksize, edata);
    if (res != -1) {
      return res;
    }
  }

  if (Param.includeTopography == YES) {
    res = topo_toexpand(leaf, ticksize, edata, Param.theFactor);
    if (res != -1) {
      return res;
    }
  }

  return wavespeed_rule(edata, Param.theFactor);
}

static int64_t g_meshImportLookupCalls = 0;
static int64_t g_meshImportSearchAttempts = 0;
static int64_t g_meshImportFallbackLookups = 0;
static int64_t g_meshImportPayloadReadySkips = 0;
static int64_t g_meshImportPendingTargetSkips = 0;

static int32_t mesh_etree_lookup_try(const octant_t* leaf, tick_t edge_ticks,
                                     int32_t fx_num, int32_t fy_num,
                                     int32_t fz_num, int32_t f_den,
                                     etree_addr_t* hit_addr,
                                     mdata_t* hit_data,
                                     etree_addr_t* search_addr_out) {
  etree_addr_t search_addr;
  int64_t dx = ((int64_t)edge_ticks * (int64_t)fx_num) / (int64_t)f_den;
  int64_t dy = ((int64_t)edge_ticks * (int64_t)fy_num) / (int64_t)f_den;
  int64_t dz = ((int64_t)edge_ticks * (int64_t)fz_num) / (int64_t)f_den;

  if (dx <= 0) {
    dx = 1;
  }
  if (dy <= 0) {
    dy = 1;
  }
  if (dz <= 0) {
    dz = 1;
  }
  if (dx >= edge_ticks) {
    dx = edge_ticks - 1;
  }
  if (dy >= edge_ticks) {
    dy = edge_ticks - 1;
  }
  if (dz >= edge_ticks) {
    dz = edge_ticks - 1;
  }

  search_addr.x = leaf->lx + (tick_t)dx;
  search_addr.y = leaf->ly + (tick_t)dy;
  search_addr.z = leaf->lz + (tick_t)dz;
  search_addr.level = ETREE_MAXLEVEL;
  search_addr.type = ETREE_LEAF;

  if (search_addr_out != NULL) {
    *search_addr_out = search_addr;
  }

  g_meshImportSearchAttempts++;
  if (etree_search(Global.theMeshInputEp, search_addr, hit_addr, NULL,
                   hit_data) != 0) {
    return -1;
  }

  return 0;
}

static int32_t mesh_etree_lookup(const octant_t* leaf, etree_addr_t* hit_addr,
                                 mdata_t* hit_data) {
  tick_t edge_ticks = (tick_t)1 << (PIXELLEVEL - leaf->level);
  etree_addr_t search_addr;
  int32_t i;

  g_meshImportLookupCalls++;

  /*
   * Robust lookup for topography-trimmed meshes: center can be in air even
   * when this cell intersects valid mesh below. Try multiple in-cell samples.
   */
  const int32_t sample_fracs[][4] = {
      {2, 2, 2, 4},  /* center */
      {2, 2, 1, 4},  /* lower-center */
      {1, 1, 1, 4},  /* lower corners */
      {3, 1, 1, 4},
      {1, 3, 1, 4},
      {3, 3, 1, 4},
      {2, 2, 3, 4},  /* upper-center */
  };

  for (i = 0; i < (int32_t)(sizeof(sample_fracs) / sizeof(sample_fracs[0]));
       i++) {
    if (mesh_etree_lookup_try(leaf, edge_ticks, sample_fracs[i][0],
                              sample_fracs[i][1], sample_fracs[i][2],
                              sample_fracs[i][3], hit_addr, hit_data,
                              &search_addr) == 0) {
      if (i > 0) {
        g_meshImportFallbackLookups++;
      }
      return 0;
    }
  }

  fprintf(stderr,
          "Thread %d: mesh_etree_lookup failed: file=%s leaf=(%d,%d,%d,l=%d) "
          "last-search=(%d,%d,%d,l=%d,t=%d) edge_ticks=%d\n",
          Global.myID, Param.mesh_etree_output_file, leaf->lx, leaf->ly,
          leaf->lz, leaf->level, (int32_t)search_addr.x,
          (int32_t)search_addr.y, (int32_t)search_addr.z,
          (int32_t)search_addr.level, (int32_t)search_addr.type,
          (int32_t)edge_ticks);
  return -1;
}

static void mesh_import_store_pending_target(void* data,
                                             int32_t target_level) {
  edata_t* edata = (edata_t*)data;

  memset(edata, 0, sizeof(edata_t));
  edata->edgesize = -1.0f;
  edata->topo_eindex = target_level;
}

static int32_t mesh_import_pending_target(const void* data) {
  const edata_t* edata = (const edata_t*)data;

  if (edata == NULL) {
    return -1;
  }

  if ((edata->edgesize < 0.0f) && (edata->topo_eindex >= 0) &&
      (edata->topo_eindex <= PIXELLEVEL)) {
    return edata->topo_eindex;
  }

  return -1;
}

static int mesh_import_payload_ready(const edata_t* edata, int32_t level,
                                     double ticksize) {
  double expected_edgesize;
  double diff;
  double tol;

  if (edata == NULL || edata->edgesize <= 0.0f) {
    return 0;
  }

  expected_edgesize = ticksize * (double)((tick_t)1 << (PIXELLEVEL - level));
  diff = fabs((double)edata->edgesize - expected_edgesize);
  tol = fmax(1.0e-3, expected_edgesize * 1.0e-5);

  return diff <= tol;
}

static void mesh_import_store_payload(void* data, const mdata_t* hit_data) {
  memcpy(data, &hit_data->edgesize, sizeof(edata_t));
  acoustic_sanitize_element((edata_t*)data);
}

static int mesh_import_use_cursor_tree(void) {
  const char* env = getenv("HERCULES_MESH_IMPORT_TREE_MODE");

  if ((env == NULL) || (env[0] == '\0')) {
    return 1;
  }

  if ((strcasecmp(env, "refine") == 0) || (strcasecmp(env, "search") == 0) ||
      (strcasecmp(env, "legacy") == 0)) {
    return 0;
  }

  return 1;
}

static void mesh_import_reset_stats(void) {
  g_meshImportLookupCalls = 0;
  g_meshImportSearchAttempts = 0;
  g_meshImportFallbackLookups = 0;
  g_meshImportPayloadReadySkips = 0;
  g_meshImportPendingTargetSkips = 0;
}

static void mesh_import_print_stats(void) {
  int64_t local_stats[5];
  int64_t max_stats[5];
  int64_t sum_stats[5];

  local_stats[0] = g_meshImportLookupCalls;
  local_stats[1] = g_meshImportSearchAttempts;
  local_stats[2] = g_meshImportFallbackLookups;
  local_stats[3] = g_meshImportPayloadReadySkips;
  local_stats[4] = g_meshImportPendingTargetSkips;

  MPI_Reduce(local_stats, max_stats, 5, MPI_INT64_T, MPI_MAX, 0, comm_solver);
  MPI_Reduce(local_stats, sum_stats, 5, MPI_INT64_T, MPI_SUM, 0, comm_solver);

  if (Global.myID == 0) {
    double ranks = (Global.theGroupSize > 0) ? (double)Global.theGroupSize : 1.0;
    double avg_attempts_per_lookup =
        (sum_stats[0] > 0)
            ? ((double)sum_stats[1] / (double)sum_stats[0])
            : 0.0;

    monitor_print(
        "Mesh import lookup stats: lookups avg %.1fM max %.1fM, "
        "search attempts avg %.1fM max %.1fM, attempts/lookup %.2f, "
        "fallback lookups max %.1fM, ready-skip max %.1fM, pending-skip max %.1fM\n",
        (double)sum_stats[0] / ranks / 1.0e6,
        (double)max_stats[0] / 1.0e6,
        (double)sum_stats[1] / ranks / 1.0e6,
        (double)max_stats[1] / 1.0e6,
        avg_attempts_per_lookup,
        (double)max_stats[2] / 1.0e6,
        (double)max_stats[3] / 1.0e6,
        (double)max_stats[4] / 1.0e6);
    fflush(stdout);
  }
}

static void mesh_import_cursor_addr_point(const etree_addr_t* addr,
                                          point_t* point) {
  point->x = (tick_t)addr->x;
  point->y = (tick_t)addr->y;
  point->z = (tick_t)addr->z;
}

static int mesh_import_cursor_before_start(const etree_addr_t* addr,
                                           const point_t* start_point) {
  point_t point;

  mesh_import_cursor_addr_point(addr, &point);
  return octor_zcompare(&point, start_point) < 0;
}

static int mesh_import_cursor_past_end(const etree_addr_t* addr,
                                       const point_t* end_point) {
  point_t point;

  if (end_point == NULL) {
    return 0;
  }

  mesh_import_cursor_addr_point(addr, &point);
  return octor_zcompare(&point, end_point) >= 0;
}

static void mesh_import_cursor_start_addr(etree_addr_t* addr,
                                          const point_t** end_point_out) {
  const point_t* intervals = octor_getintervaltable(Global.myOctree);
  const point_t* end_point = NULL;

  memset(addr, 0, sizeof(*addr));

  if ((intervals != NULL) && (Global.theGroupSize > 1)) {
    addr->x = (etree_tick_t)intervals[Global.myID].x;
    addr->y = (etree_tick_t)intervals[Global.myID].y;
    addr->z = (etree_tick_t)intervals[Global.myID].z;
    if (Global.myID + 1 < Global.theGroupSize) {
      end_point = &intervals[Global.myID + 1];
    }
  }

  addr->level = 0;
  addr->type = ETREE_INTERIOR;

  *end_point_out = end_point;
}

static int64_t mesh_import_cursor_load_tree(void) {
  etree_addr_t start_addr;
  etree_addr_t addr;
  mdata_t hit_data;
  edata_t edata;
  point_t start_point;
  const point_t* end_point = NULL;
  int64_t imported = 0;
  int cursor_active = 0;

  mesh_import_cursor_start_addr(&start_addr, &end_point);
  start_point.x = (tick_t)start_addr.x;
  start_point.y = (tick_t)start_addr.y;
  start_point.z = (tick_t)start_addr.z;

  if (octor_beginleafimport(Global.myOctree) != 0) {
    fprintf(stderr, "Thread %d: mesh_import_cursor_load_tree: begin failed\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  if (etree_initcursor(Global.theMeshInputEp, start_addr) != 0) {
    fprintf(stderr,
            "Thread %d: mesh_import_cursor_load_tree: cannot initialize etree cursor\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  cursor_active = 1;

  while (etree_getcursor(Global.theMeshInputEp, &addr, NULL, &hit_data) == 0) {
    if (mesh_import_cursor_before_start(&addr, &start_point)) {
      if (etree_advcursor(Global.theMeshInputEp) != 0) {
        if (etree_errno(Global.theMeshInputEp) != ET_END_OF_TREE) {
          fprintf(stderr,
                  "Thread %d: mesh_import_cursor_load_tree: etree cursor advance failed\n",
                  Global.myID);
          MPI_Abort(MPI_COMM_WORLD, ERROR);
          exit(1);
        }
        cursor_active = 0;
        break;
      }
      continue;
    }

    if (mesh_import_cursor_past_end(&addr, end_point)) {
      break;
    }

    if (addr.type == ETREE_LEAF) {
      mesh_import_store_payload(&edata, &hit_data);
      if (octor_importleaf(Global.myOctree, (tick_t)addr.x, (tick_t)addr.y,
                           (tick_t)addr.z, (int8_t)addr.level, &edata) != 0) {
        fprintf(stderr,
                "Thread %d: mesh_import_cursor_load_tree: cannot import leaf "
                "(%d,%d,%d,l=%d)\n",
                Global.myID, (int32_t)addr.x, (int32_t)addr.y,
                (int32_t)addr.z, (int32_t)addr.level);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }
      imported++;
    }

    if (etree_advcursor(Global.theMeshInputEp) != 0) {
      if (etree_errno(Global.theMeshInputEp) != ET_END_OF_TREE) {
        fprintf(stderr,
                "Thread %d: mesh_import_cursor_load_tree: etree cursor advance failed\n",
                Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }
      cursor_active = 0;
      break;
    }
  }

  if (cursor_active && (etree_stopcursor(Global.theMeshInputEp) != 0)) {
    fprintf(stderr,
            "Thread %d: mesh_import_cursor_load_tree: etree cursor stop failed\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  if (octor_finishleafimport(Global.myOctree) != 0) {
    fprintf(stderr, "Thread %d: mesh_import_cursor_load_tree: finish failed\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  return imported;
}

static void mesh_import_print_cursor_stats(int64_t imported) {
  int64_t min_imported = 0;
  int64_t max_imported = 0;
  int64_t total_imported = 0;

  MPI_Reduce(&imported, &min_imported, 1, MPI_INT64_T, MPI_MIN, 0,
             comm_solver);
  MPI_Reduce(&imported, &max_imported, 1, MPI_INT64_T, MPI_MAX, 0,
             comm_solver);
  MPI_Reduce(&imported, &total_imported, 1, MPI_INT64_T, MPI_SUM, 0,
             comm_solver);

  if (Global.myID == 0) {
    monitor_print("Mesh import cursor leaves: min %.1fM max %.1fM total %.1fM\n",
                  (double)min_imported / 1.0e6,
                  (double)max_imported / 1.0e6,
                  (double)total_imported / 1.0e6);
    fflush(stdout);
  }
}

static int32_t toexpand_mesh_etree(octant_t* leaf, double ticksize,
                                   const void* data) {
  etree_addr_t hit_addr;
  mdata_t hit_data;
  int32_t target_level;

  if ((data != NULL) &&
      mesh_import_payload_ready((const edata_t*)data, leaf->level, ticksize)) {
    g_meshImportPayloadReadySkips++;
    return 0;
  }

  target_level = mesh_import_pending_target(data);
  if (target_level >= 0) {
    if ((target_level >= leaf->level) && (target_level <= PIXELLEVEL)) {
      g_meshImportPendingTargetSkips++;
      return (target_level > leaf->level) ? 1 : 0;
    }
  }

  if (mesh_etree_lookup(leaf, &hit_addr, &hit_data) != 0) {
    fprintf(stderr, "Thread %d: toexpand_mesh_etree: cannot query mesh etree\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  if (hit_addr.level < leaf->level) {
    fprintf(stderr,
            "Thread %d: toexpand_mesh_etree: inconsistent mesh level (%d < %d)\n",
            Global.myID, hit_addr.level, leaf->level);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  if ((hit_addr.level == leaf->level) && (data != NULL)) {
    mesh_import_store_payload((void*)data, &hit_data);
  } else if (data != NULL) {
    mesh_import_store_pending_target((void*)data, hit_addr.level);
  }

  return (hit_addr.level > leaf->level) ? 1 : 0;
}

static void setrec_mesh_etree(octant_t* leaf, double ticksize, void* data) {
  etree_addr_t hit_addr;
  mdata_t hit_data;

  (void)ticksize;

  if (mesh_etree_lookup(leaf, &hit_addr, &hit_data) != 0) {
    fprintf(stderr, "Thread %d: setrec_mesh_etree: cannot query mesh etree\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  if (hit_addr.level < leaf->level) {
    fprintf(stderr,
            "Thread %d: setrec_mesh_etree: invalid containing level (%d < %d)\n",
            Global.myID, hit_addr.level, leaf->level);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  if (hit_addr.level == leaf->level) {
    mesh_import_store_payload(data, &hit_data);
  } else {
    mesh_import_store_pending_target(data, hit_addr.level);
  }
}

static int mesh_import_all_properties_ready(void) {
  int32_t eindex;

  for (eindex = 0; eindex < Global.myMesh->lenum; eindex++) {
    elem_t* elemp = &Global.myMesh->elemTable[eindex];
    edata_t* edata = (edata_t*)elemp->data;

    if (!mesh_import_payload_ready(edata, elemp->level,
                                   Global.myMesh->ticksize)) {
      return 0;
    }
  }

  return 1;
}

static void mesh_import_apply_properties(void) {
  int32_t eindex;

  for (eindex = 0; eindex < Global.myMesh->lenum; eindex++) {
    elem_t* elemp = &Global.myMesh->elemTable[eindex];
    edata_t* edata = (edata_t*)elemp->data;
    int32_t lnid0 = elemp->lnid[0];
    tick_t edge_ticks = (tick_t)1 << (PIXELLEVEL - elemp->level);
    etree_addr_t search_addr;
    etree_addr_t hit_addr;
    mdata_t hit_data;

    search_addr.x = Global.myMesh->nodeTable[lnid0].x + edge_ticks / 2;
    search_addr.y = Global.myMesh->nodeTable[lnid0].y + edge_ticks / 2;
    search_addr.z = Global.myMesh->nodeTable[lnid0].z + edge_ticks / 2;
    search_addr.level = ETREE_MAXLEVEL;
    search_addr.type = ETREE_LEAF;

    if (etree_search(Global.theMeshInputEp, search_addr, &hit_addr, NULL,
                     &hit_data) != 0) {
      fprintf(stderr,
              "Thread %d: mesh_import_apply_properties: etree_search failed\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

    if (hit_addr.level != elemp->level) {
      fprintf(stderr,
              "Thread %d: mesh_import_apply_properties: level mismatch (%d != %d)\n",
              Global.myID, hit_addr.level, elemp->level);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

    mesh_import_store_payload(edata, &hit_data);
  }
}

static int32_t toexpand_mesh_import_bootstrap(octant_t* leaf, double ticksize,
                                              const void* data) {
  (void)ticksize;
  (void)data;
  return (g_meshImportBootstrapLevel >= 0) &&
         (leaf->level < g_meshImportBootstrapLevel);
}

static void setrec_mesh_import_bootstrap(octant_t* leaf, double ticksize,
                                         void* data) {
  (void)leaf;
  (void)ticksize;

  if (data != NULL) {
    mesh_import_store_pending_target(data, -1);
  }
}

static void mesh_generate_from_etree(void) {
  int64_t tote, mine, maxe;
  vector3D_t IstVelModel_origin;
  vector3D_t basinVelModel_origin;
  vector3D_t threeDVelModel_origin;
  const char* verify_env = getenv("HERCULES_MESH_IMPORT_VERIFY");
  int use_cursor_tree = mesh_import_use_cursor_tree();
  int verify_import =
      (verify_env != NULL) && (strcmp(verify_env, "0") != 0) &&
      (strcasecmp(verify_env, "no") != 0) &&
      (strcasecmp(verify_env, "false") != 0) &&
      (strcasecmp(verify_env, "off") != 0);

  setenv("OCTOR_SETCOM_ALL2ALL", "1", 1);
  setenv("OCTOR_EXTRACT_PROFILE", "1", 0);
  mesh_import_reset_stats();

  if (Param.IstanbulVelocityModel == YES) {
    IstVelModel_origin = compute_domain_coords_linearinterp(
        28.675, 40.9565, Param.theSurfaceCornersLong,
        Param.theSurfaceCornersLat, Param.theDomainY, Param.theDomainX,
        &Param.utmZone);
  }
  if (Param.basinVelocityModel == YES) {
    basinVelModel_origin = compute_domain_coords_linearinterp(
        Param.theBasinLong, Param.theBasinLat, Param.theSurfaceCornersLong,
        Param.theSurfaceCornersLat, Param.theDomainY, Param.theDomainX,
        &Param.utmZone);
  }
  if (Param.threeDVelocityModel == YES) {
    threeDVelModel_origin = compute_domain_coords_linearinterp(
        Param.the3DVelocityModelLong, Param.the3DVelocityModelLat,
        Param.theSurfaceCornersLong, Param.theSurfaceCornersLat,
        Param.theDomainY, Param.theDomainX, &Param.utmZone);
  }

  if (Global.myID == 0) {
    monitor_print("Meshing: Import from existing mesh etree\n\n");
    monitor_print("Mesh import tree mode: %s\n",
                  use_cursor_tree ? "cursor" : "refine");
    monitor_print("Stage %14s Min %7s Max %5s Total    Time(s)\n\n", "", "",
                  "");
  }

  MPI_Barrier(comm_solver);
  Timer_Start("Octor Newtree");
  if (Global.myID == 0) {
    monitor_print("New tree %41s", "");
  }
  Global.myOctree = octor_newtree(
      Param.theDomainX, Param.theDomainY, Param.theDomainZ, sizeof(edata_t),
      Global.myID, Global.theGroupSize, comm_solver, get_surface_shift());
  if (Global.myOctree == NULL) {
    fprintf(stderr,
            "Thread %d: mesh_generate_from_etree: fail to create octree\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  MPI_Barrier(comm_solver);
  Timer_Stop("Octor Newtree");
  if (Global.myID == 0) {
    monitor_print("%9.2f\n", Timer_Value("Octor Newtree", (TimerKind)0));
  }

  if (Param.drmImplement == YES) {
    drm_fix_coordinates(Global.myOctree->ticksize);
  }

  if ((Global.theGroupSize == 1) && !use_cursor_tree) {
    /*
     * Multi-rank imports start from a pre-refined tree inside octor_newtree().
     * Mirror that on single-rank runs so mesh-etree refinement starts from
     * comparable leaf sizes instead of probing very coarse root cells.
     */
    g_meshImportBootstrapLevel = 4;
    if (octor_refinetree(Global.myOctree, toexpand_mesh_import_bootstrap,
                         setrec_mesh_import_bootstrap) != 0) {
      fprintf(stderr,
              "Thread %d: mesh_generate_from_etree: fail to bootstrap single-rank import tree\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    g_meshImportBootstrapLevel = -1;
  }

  Global.theMeshInputEp =
      etree_open(Param.mesh_etree_output_file, O_RDONLY, 0, 0, 0);
  if (Global.theMeshInputEp == NULL) {
    fprintf(stderr,
            "Thread %d: mesh_generate_from_etree: cannot open mesh etree %s\n",
            Global.myID, Param.mesh_etree_output_file);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  Timer_Start("Octor Refinetree");
  if (use_cursor_tree) {
    int64_t imported = mesh_import_cursor_load_tree();
    mesh_import_print_cursor_stats(imported);
  } else {
    if (octor_refinetree(Global.myOctree, toexpand_mesh_etree,
                         setrec_mesh_etree) != 0) {
      fprintf(stderr,
              "Thread %d: mesh_generate_from_etree: fail to refine octree\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
  }

  MPI_Barrier(comm_solver);
  tote = octor_getleavescount(Global.myOctree, GLOBAL);
  mine = octor_getminleavescount(Global.myOctree, GLOBAL);
  maxe = octor_getmaxleavescount(Global.myOctree, GLOBAL);

  if (Global.myID == 0) {
    monitor_print("Refining     %11" INT64_FMT " %11" INT64_FMT
                  " %11" INT64_FMT " %11.2f\n",
                  mine, maxe, tote,
                  Timer_Value("Octor Refinetree", (TimerKind)0));
    fflush(stdout);
  }

  if (verify_import && !use_cursor_tree) {
    int32_t pass;
    int converged = 0;
    int64_t prev_total = tote;

    for (pass = 0; pass < PIXELLEVEL + 1; pass++) {
      if (octor_refinetree(Global.myOctree, toexpand_mesh_etree,
                           setrec_mesh_etree) != 0) {
        fprintf(stderr,
                "Thread %d: mesh_generate_from_etree: fail to verify refined octree\n",
                Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }

      MPI_Barrier(comm_solver);
      tote = octor_getleavescount(Global.myOctree, GLOBAL);
      mine = octor_getminleavescount(Global.myOctree, GLOBAL);
      maxe = octor_getmaxleavescount(Global.myOctree, GLOBAL);

      if (Global.myID == 0) {
        monitor_print("Refine check %11" INT64_FMT " %11" INT64_FMT
                      " %11" INT64_FMT " %11.2f\n",
                      mine, maxe, tote,
                      Timer_Value("Octor Refinetree", (TimerKind)0));
        fflush(stdout);
      }

      if (tote == prev_total) {
        converged = 1;
        break;
      }

      prev_total = tote;
    }

    if (!converged) {
      fprintf(stderr,
              "Thread %d: mesh_generate_from_etree: refinement did not converge\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
  }
  Timer_Stop("Octor Refinetree");
  if (!use_cursor_tree) {
    mesh_import_print_stats();
  }

  Timer_Start("Octor Partitiontree");
  if (Global.myID == 0) {
    monitor_print("Partitioning ");
    fflush(stdout);
  }
  if (octor_partitiontree(Global.myOctree, bldgs_nodesearch_com,
                          pushdowns_nodesearch) != 0) {
    fprintf(stderr,
            "Thread %d: mesh_generate_from_etree: fail to partition tree\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  MPI_Barrier(comm_solver);
  Timer_Stop("Octor Partitiontree");
  tote = octor_getleavescount(Global.myOctree, GLOBAL);
  mine = octor_getminleavescount(Global.myOctree, GLOBAL);
  maxe = octor_getmaxleavescount(Global.myOctree, GLOBAL);
  if (Global.myID == 0) {
    monitor_print("%11" INT64_FMT " %11" INT64_FMT " %11" INT64_FMT
                  " %11.2f\n\n",
                  mine, maxe, tote,
                  Timer_Value("Octor Partitiontree", (TimerKind)0));
    fflush(stdout);
  }

  Timer_Start("Octor Extractmesh");
  if (Global.myID == 0) {
    monitor_print("Extracting the mesh %30s", "");
    fflush(stdout);
  }
  Global.myMesh =
      octor_extractmesh(Global.myOctree, bldgs_nodesearch, pushdowns_nodesearch,
                        bldgs_nodesearch_com, find_topoAirOct, topo_nodesearch);
  if (Global.myMesh == NULL) {
    fprintf(stderr,
            "Thread %d: mesh_generate_from_etree: fail to extract mesh\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  MPI_Barrier(comm_solver);
  Timer_Stop("Octor Extractmesh");
  if (Global.myID == 0) {
    monitor_print("%9.2f\n\n",
                  Timer_Value("Octor Extractmesh", (TimerKind)0));
    fflush(stdout);
  }

  Timer_Start("Mesh correct properties");
  int import_properties_ready = mesh_import_all_properties_ready();
  if (import_properties_ready) {
    if (Global.myID == 0) {
      monitor_print("Mesh properties reused from etree import\n\n");
      fflush(stdout);
    }
  } else {
    if (Global.myID == 0) {
      monitor_print("Mesh import property backfill %20s", "");
      fflush(stdout);
    }
    mesh_import_apply_properties();
  }
  MPI_Barrier(comm_solver);
  Timer_Stop("Mesh correct properties");
  if (!import_properties_ready && Global.myID == 0) {
    monitor_print("%9.2f\n\n",
                  Timer_Value("Mesh correct properties", (TimerKind)0));
    fflush(stdout);
  }

  mesh_sanitize_and_report_acoustic("Imported mesh");
  mesh_build_acoustic_interface_flags(
      (Param.IstanbulVelocityModel == YES) ? &IstVelModel_origin : NULL,
      (Param.basinVelocityModel == YES) ? &basinVelModel_origin : NULL,
      (Param.threeDVelocityModel == YES) ? &threeDVelModel_origin : NULL);

  etree_close(Global.theMeshInputEp);
  Global.theMeshInputEp = NULL;
}

/**
 * bulkload: Append the data to the end of the mesh database. Return 0 if OK,
 *       -1 on error.
 *
 */
static int32_t bulkload(etree_t* mep, mrecord_t* partTable, int32_t count) {
  int index;

  for (index = 0; index < count; index++) {
    void* payload = &partTable[index].mdata;

    if (etree_append(mep, partTable[index].addr, payload) != 0) {
      /* Append error */
      return -1;
    }
  }

  return 0;
}

/** Enumeration of the counters used in the mesh statistics */
enum mesh_stat_count_t {
  ELEMENT_COUNT,
  NODE_COUNT,
  DANGLING_COUNT,
  HARBORED_COUNT,
  MESH_COUNT_LENGTH
};

static void mesh_print_stat_imp(int32_t* st, int group_size, FILE* out) {
  int pid;
  global_id_t total_count[MESH_COUNT_LENGTH] = {0, 0, 0, 0};

  fputs(
      "\n"
      "# ------------------------------------------------------------\n"
      "# Mesh statistics:\n"
      "# ------------------------------------------------------------\n"
      "# Rank    Elements       Nodes     D-nodes     H-nodes\n",
      out);

  for (pid = 0; pid < group_size; pid++) {
    fprintf(out, "%06d %11d %11d %11d %11d\n", pid, st[ELEMENT_COUNT],
            st[NODE_COUNT], st[DANGLING_COUNT], st[HARBORED_COUNT]);

    /* add to total count */
    total_count[ELEMENT_COUNT] += st[ELEMENT_COUNT];
    total_count[NODE_COUNT] += st[NODE_COUNT];
    total_count[DANGLING_COUNT] += st[DANGLING_COUNT];
    total_count[HARBORED_COUNT] += st[HARBORED_COUNT];

    st += MESH_COUNT_LENGTH; /* move to next row */
  }

  fputs(
      "\n\n"
      "# ------------------------------------------------------------\n"
      "# Total\n"
      "# ------------------------------------------------------------\n",
      out);

  fprintf(out,
          "       %11" INT64_FMT " %11" INT64_FMT " %11" INT64_FMT
          " %11" INT64_FMT "\n\n",
          total_count[ELEMENT_COUNT], total_count[NODE_COUNT],
          total_count[DANGLING_COUNT], total_count[HARBORED_COUNT]);

  /* copy totals to static globals */
  /* TODO this should be computed through different means */
  Global.theETotal = total_count[ELEMENT_COUNT];
  Global.theNTotal = total_count[NODE_COUNT];

  /* output aggregate information to the monitor file / stdout */
  monitor_print("Total elements:                      %11" INT64_FMT
                "\n"
                "Total nodes:                         %11" INT64_FMT
                "\n"
                "Total dangling nodes:                %11" INT64_FMT "\n\n",
                total_count[ELEMENT_COUNT], total_count[NODE_COUNT],
                total_count[DANGLING_COUNT]);
}

static int mesh_collect_print_stats(local_id_t mesh_stat[MESH_COUNT_LENGTH],
                                    int my_id, int group_size,
                                    const char* fname) {
  local_id_t* st = NULL;

  if (0 == my_id) { /* only the root node allocates memory */
    XMALLOC_VAR_N(st, local_id_t, (group_size * MESH_COUNT_LENGTH));
  }

  MPI_Gather(mesh_stat, MESH_COUNT_LENGTH, MPI_INT, st, MESH_COUNT_LENGTH,
             MPI_INT, 0, comm_solver);

  if (0 == my_id) {                  /* the root node prints the stats */
    const size_t bufsize = 1048576;  // 1MB
    FILE* out = hu_fopen(Param.theMeshStatFilename, "w");

    setvbuf(out, NULL, _IOFBF, bufsize);
    mesh_print_stat_imp(st, group_size, out);
    hu_fclosep(&out);
    xfree_int32_t(&st);
  }

  return 0;
}

static void mesh_printstat_imp(const mesh_t* mesh, int my_id, int group_size,
                               const char* fname) {
  local_id_t mesh_stat[MESH_COUNT_LENGTH];

  mesh_stat[ELEMENT_COUNT] = mesh->lenum;
  mesh_stat[NODE_COUNT] = mesh->lnnum;
  mesh_stat[DANGLING_COUNT] = mesh->ldnnum;
  mesh_stat[HARBORED_COUNT] = mesh->nharbored;

  mesh_collect_print_stats(mesh_stat, my_id, group_size, fname);
}

/**
 * Gather and print mesh statistics to a file with the given name.
 *
 * \param fname Name of the file where the statistics should be stored.
 */
static void mesh_print_stat(const octree_t* oct, const mesh_t* mesh, int my_id,
                            int group_size, const char* fname) {
  /* collective function calls */
  int32_t gmin = octor_getminleaflevel(oct, GLOBAL);
  int32_t gmax = octor_getmaxleaflevel(oct, GLOBAL);

  mesh_printstat_imp(mesh, my_id, group_size, fname);

  if (Global.myID == 0) {
    monitor_print("Maximum leaf level = %d\n", gmax);
    monitor_print("Minimum leaf level = %d\n", gmin);
  }
}

/**
 * Join elements and nodes, and send to Thread 0 for output.
 */
static void mesh_output() {
  int32_t eindex;
  int32_t remains, batch, batchlimit, idx;
  mrecord_t* partTable;

  Timer_Start("Mesh Out");

  batchlimit = BATCH;

  /* Allocate a fixed size buffer space to store the join results */
  partTable = (mrecord_t*)calloc(batchlimit, sizeof(mrecord_t));
  if (partTable == NULL) {
    fprintf(stderr, "Thread %d: mesh_output: out of memory\n", Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  if (Global.myID == 0) {
    etree_t* mep;
    int32_t procid;

    printf("mesh_output ... ");

    mep = etree_open(Param.mesh_etree_output_file, O_CREAT | O_RDWR | O_TRUNC,
                     0, sizeof(mdata_t), 3);
    if (mep == NULL) {
      fprintf(stderr, "Thread 0: mesh_output: ");
      fprintf(stderr, "cannot create mesh etree\n");
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

    /* Begin an appending operation */
    if (etree_beginappend(mep, 1) != 0) {
      fprintf(stderr, "Thread 0: mesh_output: \n");
      fprintf(stderr, "cannot begin an append operation\n");
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

    eindex = 0;
    while (eindex < Global.myMesh->lenum) {
      remains = Global.myMesh->lenum - eindex;
      batch = (remains < batchlimit) ? remains : batchlimit;

      for (idx = 0; idx < batch; idx++) {
        mrecord_t* mrecord;
        int32_t whichnode;
        int32_t localnid0;

        mrecord = &partTable[idx];

        /* Fill the address field */
        localnid0 = Global.myMesh->elemTable[eindex].lnid[0];

        mrecord->addr.x = Global.myMesh->nodeTable[localnid0].x;
        mrecord->addr.y = Global.myMesh->nodeTable[localnid0].y;
        mrecord->addr.z = Global.myMesh->nodeTable[localnid0].z;
        mrecord->addr.level = Global.myMesh->elemTable[eindex].level;
        mrecord->addr.type = ETREE_LEAF;

        /* Find the global node ids for the vertices */
        for (whichnode = 0; whichnode < 8; whichnode++) {
          int32_t localnid;
          int64_t globalnid;

          localnid = Global.myMesh->elemTable[eindex].lnid[whichnode];
          globalnid = Global.myMesh->nodeTable[localnid].gnid;

          mrecord->mdata.nid[whichnode] = globalnid;
        }

        /* data points to mdata_t type */
        memcpy(&mrecord->mdata.edgesize, Global.myMesh->elemTable[eindex].data,
               sizeof(edata_t));

        eindex++;
      } /* for a batch */

      if (bulkload(mep, partTable, batch) != 0) {
        fprintf(stderr, "Thread 0: mesh_output: ");
        fprintf(stderr, "error bulk-loading data\n");
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }
    } /* for all the elements Thread 0 has */

    /* Receive data from other processors */
    for (procid = 1; procid < Global.theGroupSize; procid++) {
      MPI_Status status;
      int32_t rcvbytecount;

      /* Signal the next processor to go ahead */
      MPI_Send(NULL, 0, MPI_CHAR, procid, GOAHEAD_MSG, comm_solver);

      while (1) {
        MPI_Probe(procid, MESH_MSG, comm_solver, &status);
        MPI_Get_count(&status, MPI_CHAR, &rcvbytecount);

        batch = rcvbytecount / sizeof(mrecord_t);

        MPI_Recv(partTable, rcvbytecount, MPI_CHAR, procid, MESH_MSG,
                 comm_solver, &status);

        if (batch == 0) {
          /* Done */
          break;
        }

        if (bulkload(mep, partTable, batch) != 0) {
          fprintf(stderr, "Thread 0: mesh_output: ");
          fprintf(stderr, "cannot bulkloading data from ");
          fprintf(stderr, "Thread %d\n", procid);
          MPI_Abort(MPI_COMM_WORLD, ERROR);
          exit(1);
        }
      } /* while there is more data to be received from procid */
    } /* for all the processors */

    /* End the appending operation */
    etree_endappend(mep);

    /* Close the mep to ensure the data is on disk */
    if (etree_close(mep) != 0) {
      fprintf(stderr, "Thread 0: mesh_output ");
      fprintf(stderr, "error closing the etree database\n");
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
  } else {
    /* Processors other than 0 needs to send data to 0 */
    int32_t sndbytecount;
    MPI_Status status;

    /* Wait for me turn */
    MPI_Recv(NULL, 0, MPI_CHAR, 0, GOAHEAD_MSG, comm_solver, &status);

    eindex = 0;
    while (eindex < Global.myMesh->lenum) {
      remains = Global.myMesh->lenum - eindex;
      batch = (remains < batchlimit) ? remains : batchlimit;

      for (idx = 0; idx < batch; idx++) {
        mrecord_t* mrecord;
        int32_t whichnode;
        int32_t localnid0;

        mrecord = &partTable[idx];

        /* Fill the address field */
        localnid0 = Global.myMesh->elemTable[eindex].lnid[0];

        mrecord->addr.x = Global.myMesh->nodeTable[localnid0].x;
        mrecord->addr.y = Global.myMesh->nodeTable[localnid0].y;
        mrecord->addr.z = Global.myMesh->nodeTable[localnid0].z;
        mrecord->addr.level = Global.myMesh->elemTable[eindex].level;
        mrecord->addr.type = ETREE_LEAF;

        /* Find the global node ids for the vertices */
        for (whichnode = 0; whichnode < 8; whichnode++) {
          int32_t localnid;
          int64_t globalnid;

          localnid = Global.myMesh->elemTable[eindex].lnid[whichnode];
          globalnid = Global.myMesh->nodeTable[localnid].gnid;

          mrecord->mdata.nid[whichnode] = globalnid;
        }

        memcpy(&mrecord->mdata.edgesize, Global.myMesh->elemTable[eindex].data,
               sizeof(edata_t));

        eindex++;
      } /* for a batch */

      /* Send data to proc 0 */
      sndbytecount = batch * sizeof(mrecord_t);
      MPI_Send(partTable, sndbytecount, MPI_CHAR, 0, MESH_MSG, comm_solver);
    } /* While there is data left to be sent */

    /* Send an empty message to indicate the end of my transfer */
    MPI_Send(NULL, 0, MPI_CHAR, 0, MESH_MSG, comm_solver);
  }

  /* Free the memory for the partial join results */
  free(partTable);

  Timer_Stop("Mesh Out");

  if (Global.myID == 0) {
    printf("done : %9.2f seconds\n", Timer_Value("Mesh Out", (TimerKind)0));
  }

  return;
}

/*-----------Computation routines ------------------------------*/

/*
   Macros to facitilate computation

   INTEGRAL_1: Integral_Delfixl_Delfjxl()
   INTEGRAL_2: Integral_Delfixl_Delfjxm()
 */

#define INTEGRAL_1(xki, xkj, xli, xlj, xmi, xmj) \
  (4.5 * xki * xkj * (1 + xli * xlj / 3) * (1 + xmi * xmj / 3) / 8)

#define INTEGRAL_2(xki, xlj, xmi, xmj) \
  (4.5 * xki * xlj * (1 + xmi * xmj / 3) / 8)

#define DS_TOTAL_INTERVALS 40
#define DS_TOTAL_PARAMETERS 6

/**
 * Compute histograms for xi, zeta and associated values to understand
 * what is happenning with those parameters involved the damping and
 * delta T.
 */
static void damping_statistics(double min_xi, double max_xi, double min_zeta,
                               double max_zeta, double min_VsVp,
                               double max_VsVp, double min_VpVsZ,
                               double max_VpVsZ, double min_Vs, double max_Vs) {
  static const int totalintervals = DS_TOTAL_INTERVALS;
  static const int totalparameters = DS_TOTAL_PARAMETERS;

  int interval, parameter, row, col, matrixelements;

  double min_VpVs, max_VpVs;

  double themins[DS_TOTAL_PARAMETERS], themaxs[DS_TOTAL_PARAMETERS],
      spacing[DS_TOTAL_PARAMETERS];

  int32_t counters[DS_TOTAL_PARAMETERS][DS_TOTAL_INTERVALS],
      global_counter[DS_TOTAL_PARAMETERS][DS_TOTAL_INTERVALS],
      global_total[DS_TOTAL_PARAMETERS], eindex;

  /* Initializing clue values and variables */

  min_VpVs = (max_VsVp > 0.0) ? (1 / max_VsVp) : 0.0;
  max_VpVs = (min_VsVp > 0.0) ? (1 / min_VsVp) : 0.0;

  themins[0] = min_zeta;
  themins[1] = min_xi;
  themins[2] = min_VsVp;
  themins[3] = min_VpVs;
  themins[4] = min_VpVsZ;
  themins[5] = min_Vs;

  themaxs[0] = max_zeta;
  themaxs[1] = max_xi;
  themaxs[2] = max_VsVp;
  themaxs[3] = max_VpVs;
  themaxs[4] = max_VpVsZ;
  themaxs[5] = max_Vs;

  for (row = 0; row < totalparameters; row++) {
    for (col = 0; col < totalintervals; col++) {
      counters[row][col] = 0;
    }
    global_total[row] = 0;
  }

  for (row = 0; row < totalparameters; row++) {
    spacing[row] = (themaxs[row] - themins[row]) / totalintervals;
  }

  /* loop over the elements */
  for (eindex = 0; eindex < Global.myMesh->lenum; eindex++) {
    /* loop variables */
    elem_t* elemp;
    edata_t* edata;
    double a, b, omega, elementvalues[6];

    /* capturing the elements */
    elemp = &Global.myMesh->elemTable[eindex];
    edata = (edata_t*)elemp->data;

    /* the parameters */
    elementvalues[0] = element_is_acoustic(edata) ? 0.0 : (25.0 / edata->Vs);
    /* (edata->Vs < 1500) ? 25 / edata->Vs : 5 / edata->Vs; */
    /* zeta */
    omega = 3.46410161514 / (edata->edgesize / edata->Vp);  /* freq in rad */
    a = elementvalues[0] * Global.theABase;                 /* a     */
    b = elementvalues[0] * Global.theBBase;                 /* b     */
    elementvalues[1] = (a / (2 * omega)) + (b * omega / 2); /* xi    */
    elementvalues[2] = (edata->Vs / edata->Vp); /* Vs/Vp */
    elementvalues[3] =
        element_is_acoustic(edata) ? 0.0 : (edata->Vp / edata->Vs); /* Vp/Vs */
    elementvalues[4] = element_is_acoustic(edata)
                           ? 0.0
                           : elementvalues[0] * (edata->Vp / edata->Vs);
    /* Vp/Vs*zeta  */
    elementvalues[5] = edata->Vs; /* Vs    */

    /* loop over the parameters */
    for (parameter = 0; parameter < totalparameters; parameter++) {
      /* loop over each interval */
      for (interval = 0; interval < totalintervals; interval++) {
        /* loop variables */
        double liminf, limsup;

        /* histogram limits */
        liminf = themins[parameter] + (interval * spacing[parameter]);
        limsup = liminf + spacing[parameter];

        /* for the last interval adjust to the max value */
        if (interval == totalintervals - 1) {
          limsup = themaxs[parameter];
        }

        /* counting elements within the interval */
        if ((elementvalues[parameter] > liminf) &&
            (elementvalues[parameter] <= limsup)) {
          counters[parameter][interval]++;
        }
      } /* ends loop on intervals */
    } /* ends loop on parameters */
  } /*ends loop on elements */

  /* add all counting results from each processor */
  matrixelements = totalparameters * totalintervals;
  MPI_Reduce(&counters[0][0], &global_counter[0][0], matrixelements, MPI_INT,
             MPI_SUM, 0, comm_solver);
  MPI_Bcast(&global_counter[0][0], matrixelements, MPI_INT, 0, comm_solver);

  /* sums the total of elements for each histogram */
  if (Global.myID == 0) {
    for (parameter = 0; parameter < totalparameters; parameter++) {
      global_counter[parameter][0]++;
      for (interval = 0; interval < totalintervals; interval++) {
        global_total[parameter] =
            global_total[parameter] + global_counter[parameter][interval];
      }
    }
  }

  /* MPI Barrier */
  MPI_Barrier(comm_solver);

  /* prints to the terminal the histograms */
  if (Global.myID == 0) {
    /* header to identify each column */
    printf("\n\n\tThe histograms of the following parameters: \n\n");
    printf("\t 1. Zeta\n");
    printf("\t 2. Xi\n");
    printf("\t 3. Vs/Vp\n");
    printf("\t 4. Vp/Vs\n");
    printf("\t 5. Vp/Vs*zeta\n");
    printf("\t 6. Vs\n\n");
    printf("\tAre given in the following table\n");
    printf("\t(each column is one of the parameters)\n\n\t");

    /* printing the histograms */
    for (interval = 0; interval < totalintervals; interval++) {
      for (parameter = 0; parameter < totalparameters; parameter++) {
        printf("%12d", global_counter[parameter][interval]);
      }
      printf("\n\t");
    }

    /* prints the total of elements for each column */
    printf("\n\tTotals:\n\n\t");
    for (parameter = 0; parameter < totalparameters; parameter++) {
      printf("%12d", global_total[parameter]);
    }

    /* prints the interval witdth */
    printf("\n\n\tAnd the intervals width is:");
    for (parameter = 0; parameter < totalparameters; parameter++) {
      printf("\n\t %2d. %.6f ", parameter + 1, spacing[parameter]);
    }
    printf("\n\n");
    fflush(stdout);
  }

  return;
} /* end damping_statistics */

/**
 * Determine the limit values associated with the damping and delta_t problem.
 */
static void solver_set_critical_T() {
  int32_t eindex; /* element index */

  double min_h_over_Vp = 1e32; /* the min h/Vp group */
  double min_h_over_Vp_global;
  int32_t min_h_over_Vp_elem_index = -1;

  double min_dt_factor_X = 1e32; /* the min delta_t group */
  double min_dt_factor_Z = 1e32, min_dt_factor_X_global, min_dt_factor_Z_global;
  int32_t min_dt_factor_X_elem_index = -1, min_dt_factor_Z_elem_index = -1;

  double min_zeta = 1e32; /* the zeta group */
  double max_zeta = 0, min_zeta_global, max_zeta_global;
  int32_t min_zeta_elem_index = -1, max_zeta_elem_index = -1;

  double min_xi = 1e32; /* the xi group */
  double max_xi = 0, min_xi_global, max_xi_global;
  int32_t min_xi_elem_index = -1, max_xi_elem_index = -1;

  double min_VsVp = 1e32; /* the Vs/Vp group */
  double min_VsVp_global, max_VsVp = 0, max_VsVp_global;
  int32_t min_VsVp_elem_index = -1, max_VsVp_elem_index = -1;

  double min_VpVsZ = 1e32; /* the Vp/Vs group */
  double min_VpVsZ_global, max_VpVsZ = 0, max_VpVsZ_global;
  int32_t min_VpVsZ_elem_index = -1, max_VpVsZ_elem_index = -1;

  double min_Vs = 1e32; /* the Vs group */
  double min_Vs_global, max_Vs = 0, max_Vs_global;
  int32_t min_Vs_elem_index = -1, max_Vs_elem_index = -1;

  /* Find the minima and maxima for all needed coefficients */
  /* Start loop over the mesh elements */
  for (eindex = 0; eindex < Global.myMesh->lenum; eindex++) {
    /* Loop local variables */

    elem_t* elemp;  /* pointer to the mesh database */
    edata_t* edata; /* pointer to the element data  */

    double ratio;       /* the h/Vp ratio                   */
    double zeta;        /* the time domain zeta-damping     */
    double xi;          /* the freq domain xi-damping       */
    double omega;       /* the element associated freq from w=3.46*Vp/h */
    double a, b;        /* the same constants we use for C = aM + bK    */
    double dt_factor_X; /* the factor of 0.577(1-xi)*h/Vp           */
    double dt_factor_Z; /* the factor of 0.577(1-zeta)*h/Vp         */
    double VsVp;        /* the quotient Vs/Vp                       */
    double VpVsZ;       /* the result of Vp / Vs * zeta             */
    double Vs;          /* the Vs                                   */

    /* yigit says -- Buildings assume 5% damping */
    int64_t n_0;
    double z_m;

    /* Captures the element */

    elemp = &Global.myMesh->elemTable[eindex];
    edata = (edata_t*)elemp->data;

    /* Calculate the clue quantities */

    if (edata->Vp != -1)
      ratio = edata->edgesize / edata->Vp;
    else
      ratio = 1e32;

    /* Old formula for damping */
    /* zeta        = (edata->Vs < 1500) ? 25 / edata->Vs : 5 / edata->Vs; */
    /* New formula acording to Graves */

    n_0 = Global.myMesh->elemTable[eindex].lnid[0];
    z_m = Global.theZForMeshOrigin +
          (Global.myOctree->ticksize) * Global.myMesh->nodeTable[n_0].z;

    /* Shift the domain if buildings are considered */
    if (Param.includeBuildings == YES) {
      z_m -= get_surface_shift();
    }

    zeta = element_is_acoustic(edata) ? 0.0 : (25.0 / edata->Vs);

    /*If element is in the building, use 5% damping.*/
    if (z_m < 0) {
      zeta = 0.05;
    }

    // Clifford's NOTE: 3.46410161514 = sqrt(12) = 2*sqrt(3)
    omega = 3.46410161514 / ratio;
    a = zeta * Global.theABase;
    b = zeta * Global.theBBase;
    xi = (a / (2 * omega)) + (b * omega / 2);
    // Clifford's NOTE: 0.57735026919 = 1/sqrt(3)
    dt_factor_X = 0.57735026919 * (1 - xi) * ratio;
    dt_factor_Z = 0.57735026919 * (1 - zeta) * ratio;
    if (element_is_acoustic(edata)) {
      VsVp = 0.0;
      VpVsZ = 0.0;
    } else {
      VsVp = edata->Vs / edata->Vp;
      VpVsZ = zeta * (edata->Vp / edata->Vs);
    }
    Vs = edata->Vs;

    /* Updating for extreme values */

    /* ratio */
    if (ratio < min_h_over_Vp) {
      min_h_over_Vp = ratio;
      min_h_over_Vp_elem_index = eindex;
    }

    /* dt_factors */
    if (dt_factor_X < min_dt_factor_X) {
      min_dt_factor_X = dt_factor_X;
      min_dt_factor_X_elem_index = eindex;
    }
    if (dt_factor_Z < min_dt_factor_Z) {
      min_dt_factor_Z = dt_factor_Z;
      min_dt_factor_Z_elem_index = eindex;
    }

    /* min_zeta and max_zeta */
    if (zeta < min_zeta) {
      min_zeta = zeta;
      min_zeta_elem_index = eindex;
    }
    if (zeta > max_zeta) {
      max_zeta = zeta;
      max_zeta_elem_index = eindex;
    }

    /* min_xi and max_xi */
    if (xi < min_xi) {
      min_xi = xi;
      min_xi_elem_index = eindex;
    }
    if (xi > max_xi) {
      max_xi = xi;
      max_xi_elem_index = eindex;
    }

    /* min Vs/Vp */
    if (VsVp < min_VsVp) {
      min_VsVp = VsVp;
      min_VsVp_elem_index = eindex;
    }
    if (VsVp > max_VsVp) {
      max_VsVp = VsVp;
      max_VsVp_elem_index = eindex;
    }

    /* min and max VpVsZ */
    if (VpVsZ < min_VpVsZ) {
      min_VpVsZ = VpVsZ;
      min_VpVsZ_elem_index = eindex;
    }
    if (VpVsZ > max_VpVsZ) {
      max_VpVsZ = VpVsZ;
      max_VpVsZ_elem_index = eindex;
    }

    /* min Vs */
    if (Vs < min_Vs) {
      min_Vs = Vs;
      min_Vs_elem_index = eindex;
    }
    if (Vs > max_Vs) {
      max_Vs = Vs;
      max_Vs_elem_index = eindex;
    }

  } /* End of the loop over the mesh elements */

  /* Reducing to global values */
  MPI_Reduce(&min_h_over_Vp, &min_h_over_Vp_global, 1, MPI_DOUBLE, MPI_MIN, 0,
             comm_solver);
  if (Param.theDampingStatisticsFlag == 1) {
    MPI_Reduce(&min_dt_factor_X, &min_dt_factor_X_global, 1, MPI_DOUBLE,
               MPI_MIN, 0, comm_solver);
    MPI_Reduce(&min_dt_factor_Z, &min_dt_factor_Z_global, 1, MPI_DOUBLE,
               MPI_MIN, 0, comm_solver);
    MPI_Reduce(&min_zeta, &min_zeta_global, 1, MPI_DOUBLE, MPI_MIN, 0,
               comm_solver);
    MPI_Reduce(&max_zeta, &max_zeta_global, 1, MPI_DOUBLE, MPI_MAX, 0,
               comm_solver);
    MPI_Reduce(&min_xi, &min_xi_global, 1, MPI_DOUBLE, MPI_MIN, 0, comm_solver);
    MPI_Reduce(&max_xi, &max_xi_global, 1, MPI_DOUBLE, MPI_MAX, 0, comm_solver);
    MPI_Reduce(&min_VsVp, &min_VsVp_global, 1, MPI_DOUBLE, MPI_MIN, 0,
               comm_solver);
    MPI_Reduce(&max_VsVp, &max_VsVp_global, 1, MPI_DOUBLE, MPI_MAX, 0,
               comm_solver);
    MPI_Reduce(&min_VpVsZ, &min_VpVsZ_global, 1, MPI_DOUBLE, MPI_MIN, 0,
               comm_solver);
    MPI_Reduce(&max_VpVsZ, &max_VpVsZ_global, 1, MPI_DOUBLE, MPI_MAX, 0,
               comm_solver);
    MPI_Reduce(&min_Vs, &min_Vs_global, 1, MPI_DOUBLE, MPI_MIN, 0, comm_solver);
    MPI_Reduce(&max_Vs, &max_Vs_global, 1, MPI_DOUBLE, MPI_MAX, 0, comm_solver);
  }

  /* Inform everyone about the global values */
  MPI_Bcast(&min_h_over_Vp_global, 1, MPI_DOUBLE, 0, comm_solver);
  if (Param.theDampingStatisticsFlag == 1) {
    MPI_Bcast(&min_dt_factor_X_global, 1, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(&min_dt_factor_Z_global, 1, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(&min_zeta_global, 1, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(&max_zeta_global, 1, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(&min_xi_global, 1, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(&max_xi_global, 1, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(&min_VsVp_global, 1, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(&max_VsVp_global, 1, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(&min_VpVsZ_global, 1, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(&max_VpVsZ_global, 1, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(&min_Vs_global, 1, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(&max_Vs_global, 1, MPI_DOUBLE, 0, comm_solver);
  }

  /* go for damping statistics */
  if (Param.theDampingStatisticsFlag == 1) {
    damping_statistics(min_xi_global, max_xi_global, min_zeta_global,
                       max_zeta_global, min_VsVp_global, max_VsVp_global,
                       min_VpVsZ_global, max_VpVsZ_global, min_Vs_global,
                       max_Vs_global);
  }

  /* Static global variable for the critical delta t */
  Global.theCriticalT = min_h_over_Vp_global;

  /* CPML imposes a stricter CFL constraint to keep the corner-element
   * multi-axis ADE recurrence stable. Default safety factor is 0.9. */
  if (Param.theAbcType == ABC_CPML) {
    Global.theCriticalT *= Param.theCpmlCflSafety;
    if (Global.myID == 0) {
      monitor_print("\n CPML CFL safety factor %.3f applied: critical dt = %.6e s\n",
                    Param.theCpmlCflSafety, Global.theCriticalT);
    }
  }

  /* Printing of information */
  MPI_Barrier(comm_solver);
  if (Global.myID == 0) {
    if (Param.theDampingStatisticsFlag == 1) {
      monitor_print("\n\n Critical delta t related information: \n\n");
      monitor_print("\t 1. The minimum h/Vp     = %.6f \n",
                    min_h_over_Vp_global);
      monitor_print("\t 2. The minimum dt X     = %.6f \n",
                    min_dt_factor_X_global);
      monitor_print("\t 3. The minimum dt Z     = %.6f \n",
                    min_dt_factor_Z_global);
      monitor_print("\t 4. The minimum zeta     = %.6f \n", min_zeta_global);
      monitor_print("\t 5. The maximum zeta     = %.6f \n", max_zeta_global);
      monitor_print("\t 6. The minimum xi       = %.6f \n", min_xi_global);
      monitor_print("\t 7. The maximum xi       = %.6f \n", max_xi_global);
      monitor_print("\t 8. The minimum Vs/Vp    = %.6f \n", min_VsVp_global);
      monitor_print("\t 9. The maximum Vs/Vp    = %.6f \n", max_VsVp_global);
      monitor_print("\t10. The minimum (Vp/Vs)*zeta = %.6f \n",
                    min_VpVsZ_global);
      monitor_print("\t11. The maximum (Vp/Vs)*zeta = %.6f \n",
                    max_VpVsZ_global);
      monitor_print("\t12. The minimum Vs       = %.6f \n", min_Vs_global);
      monitor_print("\t13. The maximum Vs       = %.6f \n", max_Vs_global);
    } else {
      monitor_print("\n\n Critical delta t related information: \n\n");
      monitor_print("\t The minimum h/Vp = %.6f \n\n", min_h_over_Vp_global);
    }
  }

#ifdef AUTO_DELTA_T
  /* Set the critical delta T */
  Param.theDeltaT = Global.theCriticalT;
  Param.theDeltaTSquared = Param.theDeltaT * Param.theDeltaT;

  /* Set the total steps */
  Param.theTotalSteps =
      (int)(((Param.theEndT - Param.theStartT) / Param.theDeltaT));
#endif /* AUTO_DELTA_T */

  /* Printing location and element properties of the maximum values */
  if (Param.theDampingStatisticsFlag == 1) {
    /* Local variables */

    double local_extremes[13], global_extremes[13];
    int32_t element_indices[13];
    int32_t extreme_index;

    local_extremes[0] = min_h_over_Vp;
    local_extremes[1] = min_dt_factor_X;
    local_extremes[2] = min_dt_factor_Z;
    local_extremes[3] = min_zeta;
    local_extremes[4] = max_zeta;
    local_extremes[5] = min_xi;
    local_extremes[6] = max_xi;
    local_extremes[7] = min_VsVp;
    local_extremes[8] = max_VsVp;
    local_extremes[9] = min_VpVsZ;
    local_extremes[10] = max_VpVsZ;
    local_extremes[11] = min_Vs;
    local_extremes[12] = max_Vs;

    global_extremes[0] = min_h_over_Vp_global;
    global_extremes[1] = min_dt_factor_X_global;
    global_extremes[2] = min_dt_factor_Z_global;
    global_extremes[3] = min_zeta_global;
    global_extremes[4] = max_zeta_global;
    global_extremes[5] = min_xi_global;
    global_extremes[6] = max_xi_global;
    global_extremes[7] = min_VsVp_global;
    global_extremes[8] = max_VsVp_global;
    global_extremes[9] = min_VpVsZ_global;
    global_extremes[10] = max_VpVsZ_global;
    global_extremes[11] = min_Vs_global;
    global_extremes[12] = max_Vs_global;

    element_indices[0] = min_h_over_Vp_elem_index;
    element_indices[1] = min_dt_factor_X_elem_index;
    element_indices[2] = min_dt_factor_Z_elem_index;
    element_indices[3] = min_zeta_elem_index;
    element_indices[4] = max_zeta_elem_index;
    element_indices[5] = min_xi_elem_index;
    element_indices[6] = max_xi_elem_index;
    element_indices[7] = min_VsVp_elem_index;
    element_indices[8] = max_VsVp_elem_index;
    element_indices[9] = min_VpVsZ_elem_index;
    element_indices[10] = max_VpVsZ_elem_index;
    element_indices[11] = min_Vs_elem_index;
    element_indices[12] = max_Vs_elem_index;

    /* Printing section title */
    MPI_Barrier(comm_solver);
    if (Global.myID == 0) {
      printf(
          "\n\t Their corresponding element properties and coordinates are: "
          "\n\n");
    }

    /* Loop over the six extreme values */
    MPI_Barrier(comm_solver);
    for (extreme_index = 0; extreme_index < 6; extreme_index++) {
      MPI_Barrier(comm_solver);
      if (local_extremes[extreme_index] == global_extremes[extreme_index]) {
        tick_t ldb[3];
        elem_t* elemp;
        edata_t* edata;
        int lnid0 =
            Global.myMesh->elemTable[element_indices[extreme_index]].lnid[0];

        ldb[0] = Global.myMesh->nodeTable[lnid0].x;
        ldb[1] = Global.myMesh->nodeTable[lnid0].y;
        ldb[2] = Global.myMesh->nodeTable[lnid0].z;

        elemp = &Global.myMesh->elemTable[element_indices[extreme_index]];
        edata = (edata_t*)elemp->data;

        printf("\t For extreme value No. %d:", extreme_index + 1);
        printf("\n\t\t h = %.6f, Vp = %.6f Vs = %.6f rho = %.6f",
               edata->edgesize, edata->Vp, edata->Vs, edata->rho);
        printf("\n\t\t by thread %d, element_coord = (%.6f, %.6f, %.6f)\n\n",
               Global.myID, ldb[0] * Global.myMesh->ticksize,
               ldb[1] * Global.myMesh->ticksize,
               ldb[2] * Global.myMesh->ticksize);
      }
      MPI_Barrier(comm_solver);
    } /* End of loop over the extreme values */

    if (Global.myID == 0) {
      fflush(stdout);
    }
  } /* end if damping statistics */

  return;
} /* End solver_set_critical_T */

/**
 * Iterates through all processor to obtain the minimum * edgesize of the
 * mesh.
 */
static void get_minimum_edge() {
  int32_t eindex;
  double min_h = 1e32, min_h_global;
  // int32_t min_h_elem_index;

  /* Find the minimal h/Vp in the domain */
  for (eindex = 0; eindex < Global.myMesh->lenum; eindex++) {
    elem_t* elemp; /* pointer to the mesh database */
    edata_t* edata;
    double h;

    elemp = &Global.myMesh->elemTable[eindex];
    edata = (edata_t*)elemp->data;

    /* Update the min_h value. (h = edgesize)  */
    h = edata->edgesize;
    if (h < min_h) {
      min_h = h;
      // min_h_elem_index = eindex;
    }
  }

  MPI_Reduce(&min_h, &min_h_global, 1, MPI_DOUBLE, MPI_MIN, 0, comm_solver);
  /* Inform everyone of this value */
  MPI_Bcast(&min_h_global, 1, MPI_DOUBLE, 0, comm_solver);

  if (Global.myID == 0) {
    monitor_print("\nThe minimal h  = %.6f\n\n\n", min_h_global);
  }

  Global.theNumericsInformation.minimumh = min_h_global;

  return;
}

/**
 * Print the stiffness matrix K1, K2 and K3 to the given output stream.
 */
static void print_K_stdoutput() {
  int i, j, iloc, jloc;

  monitor_print("\n\nStiffness Matrix K1 \n\n");
  for (i = 0; i < 8; i++) {
    for (iloc = 0; iloc < 3; iloc++) {
      for (j = 0; j < 8; j++) {
        for (jloc = 0; jloc < 3; jloc++) {
          monitor_print("%10.2e", Global.theK1[i][j].f[iloc][jloc]);
        }
      }
      monitor_print("\n");
    }
  }

  monitor_print("\n\nStiffness Matrix K2 \n\n");
  for (i = 0; i < 8; i++) {
    for (iloc = 0; iloc < 3; iloc++) {
      for (j = 0; j < 8; j++) {
        for (jloc = 0; jloc < 3; jloc++) {
          monitor_print("%10.2e", Global.theK2[i][j].f[iloc][jloc]);
        }
      }
      monitor_print("\n");
    }
  }

  monitor_print("\n\nStiffness Matrix K3 \n\n");
  for (i = 0; i < 8; i++) {
    for (iloc = 0; iloc < 3; iloc++) {
      for (j = 0; j < 8; j++) {
        for (jloc = 0; jloc < 3; jloc++) {
          monitor_print("%10.2e", Global.theK3[i][j].f[iloc][jloc]);
        }
      }
      monitor_print("\n");
    }
  }

  monitor_print("\n\n");
}

/**
 * mu_and_lambda: Calculates mu and lambda according to the element values
 *                of Vp, Vs, and Rho and verifies/applies some rules defined
 *                by Jacobo, Leonardo and Ricardo.  It was originally within
 *                solver_init but was moved out because it needs to be used
 *                in other places as well (nonlinear)
 */
void mu_and_lambda(double* theMu, double* theLambda, edata_t* edata,
                   int32_t eindex) {
  double mu, lambda;

  // This is added for the treatment of air elements
  if ((edata->Vp == -1) && (edata->rho == 0)) {
    *theMu = 0;
    *theLambda = 0;
    return;
  }

  if (element_is_acoustic(edata)) {
    *theMu = 0.0;
    *theLambda = edata->rho * edata->Vp * edata->Vp;
    return;
  }

  mu = edata->rho * edata->Vs * edata->Vs;

  if (edata->Vp > (edata->Vs * Param.theThresholdVpVs)) {
    lambda = edata->rho * edata->Vs * edata->Vs * Param.theThresholdVpVs *
                 Param.theThresholdVpVs -
             2 * mu;
  } else {
    lambda = edata->rho * edata->Vp * edata->Vp - 2 * mu;
  }

  /* Adjust Vs, Vp to fix Poisson ratio problem, formula provided by Jacobo */
  if (lambda < 0) {
    if (edata->Vs < 500)
      edata->Vp = 2.45 * edata->Vs;
    else if (edata->Vs < 1200)
      edata->Vp = 2 * edata->Vs;
    else
      edata->Vp = 1.87 * edata->Vs;

    lambda = edata->rho * edata->Vp * edata->Vp;
  }

  if (lambda < 0) {
    fprintf(stderr,
            "\nThread %d: %d element produces negative lambda = %.6f; Vp = %f; "
            "Vs = %f; Rho = %f",
            Global.myID, eindex, lambda, edata->Vp, edata->Vs, edata->rho);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
  }

  /* assign results to returning values */
  *theMu = mu;
  *theLambda = lambda;
}

/**
 * gpu_init: Select one of the attached GPU devices and retrieve its
 *           specifications.
 *
 */
static void gpu_init(int32_t myID) {
  int myDevID, devCount;
  int localRank = 0;
  const char* rankEnv;
  const char* allowOversubEnv;
  int allowOversub = 0;

  /* Select a GPU device for this process */
  if (cudaGetDeviceCount(&devCount) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to get GPU device count\n", myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  if (devCount == 0) {
    fprintf(stderr, "Thread %d: NO GPU devices found\n", myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  rankEnv = getenv("OMPI_COMM_WORLD_LOCAL_RANK");
  if (rankEnv == NULL) {
    rankEnv = getenv("MV2_COMM_WORLD_LOCAL_RANK");
  }
  if (rankEnv == NULL) {
    rankEnv = getenv("SLURM_LOCALID");
  }
  if (rankEnv == NULL) {
    rankEnv = getenv("MPI_LOCALRANKID");
  }
  if (rankEnv != NULL) {
    localRank = atoi(rankEnv);
  }

  allowOversubEnv = getenv("HERCULES_ALLOW_MULTI_RANK_PER_GPU");
  if ((allowOversubEnv != NULL) &&
      ((strcmp(allowOversubEnv, "1") == 0) ||
       (strcasecmp(allowOversubEnv, "yes") == 0) ||
       (strcasecmp(allowOversubEnv, "true") == 0))) {
    allowOversub = 1;
  }

  if ((localRank >= devCount) && (allowOversub == 0)) {
    fprintf(stderr,
            "Thread %d: local MPI rank %d oversubscribes %d GPU(s). "
            "Run one MPI rank per GPU or set HERCULES_ALLOW_MULTI_RANK_PER_GPU=1 to override.\n",
            myID, localRank, devCount);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  myDevID = (localRank % devCount);
  if (cudaSetDevice(myDevID) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to set device\n", myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  /* Reset the device */
  if (cudaDeviceReset() != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to reset device\n", myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  /* Set cache preference for larger L1 */
  // if (cudaDeviceSetCacheConfig(cudaFuncCachePreferL1) != cudaSuccess) {
  //   fprintf(stderr, "Thread %d: Failed to reset device cache config\n",
  //   myID); MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
  // }

  /* Set the shared memory bank size to 8 */
  cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);

  /* Retrieve device hardware details */
  getGPUHardware(myID, myDevID, &(Global.gpu_spec));

  /* Initialize GPU metrics */
  Global.gpu_spec.numdevices = devCount;
  Global.gpu_spec.numbytespci = 0;
  Global.gpu_spec.numflops = 0;
  Global.gpu_spec.numbytes = 0;

  return;
}

/**
 * Init matrices and constants, build comm schedule, allocate/init space
 * for the solver.
 */
static void solver_init() {
  /* local variables */
  int32_t eindex;
  int32_t c_outsize, c_insize, s_outsize, s_insize;
  int64_t localPressureElements = 0;
  int64_t globalPressureElements = 0;
  solver_float* pressureFreeSurfaceHost = NULL;

  /* Initialize the GPU */
  gpu_init(Global.myID);

  /* compute the damping parameters a/zeta and b/zeta */
  compute_setab(Param.theFreq, &Global.theABase, &Global.theBBase);
  if (Global.myID == 0) {
    printf("ABase = %f, BBase = %f\n", Global.theABase, Global.theBBase);
  }

  /* find out the critical delta T of the current simulation */
  /* and goes for the damping statistics if falg is == 1     */
  MPI_Barrier(comm_solver);
  solver_set_critical_T();

  /* find the minimum edge size */

  get_minimum_edge();

  /* Init stiffness matrices and other constants */
  compute_K();

  /* For debugging */
  if ((Param.printK == YES) && (Global.myID == 0)) {
    print_K_stdoutput();
  }

  /* allocation of memory */
  Global.mySolver = (mysolver_t*)malloc(sizeof(mysolver_t));
  if (Global.mySolver == NULL) {
    fprintf(stderr, "Thread %d: solver_init: out of memory\n", Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  memset(Global.mySolver, 0, sizeof(mysolver_t));

  /* Allocate memory */
  Global.mySolver->eTable = (e_t*)calloc(Global.myMesh->lenum, sizeof(e_t));
  Global.mySolver->nTable = (n_t*)calloc(Global.myMesh->nharbored, sizeof(n_t));
  /* Global.mySolver->tm1 = */
  /*     (fvector_t*)calloc(Global.myMesh->nharbored, sizeof(fvector_t)); */
  /* Global.mySolver->tm2 = */
  /*     (fvector_t*)calloc(Global.myMesh->nharbored, sizeof(fvector_t)); */
  /* Global.mySolver->force = */
  /*     (fvector_t*)calloc(Global.myMesh->nharbored, sizeof(fvector_t)); */
  if (cudaMallocHost(&(Global.mySolver->tm1),
                     Global.myMesh->nharbored * sizeof(fvector_t)) !=
      cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to allocate tm1 pinned memory on host\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  memset(Global.mySolver->tm1, 0, Global.myMesh->nharbored * sizeof(fvector_t));
  if (cudaMallocHost(&(Global.mySolver->tm2),
                     Global.myMesh->nharbored * sizeof(fvector_t)) !=
      cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to allocate tm2 pinned memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  memset(Global.mySolver->tm2, 0, Global.myMesh->nharbored * sizeof(fvector_t));
  if (cudaMallocHost(&(Global.mySolver->force),
                     Global.myMesh->nharbored * sizeof(fvector_t)) !=
      cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to allocate force pinned memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  memset(Global.mySolver->force, 0,
         Global.myMesh->nharbored * sizeof(fvector_t));

  if (acoustic_pressure_enabled()) {
    if (cudaMallocHost(&(Global.mySolver->pressure_tm1),
                       Global.myMesh->nharbored * sizeof(solver_float)) !=
        cudaSuccess) {
      fprintf(stderr,
              "Thread %d: Failed to allocate acoustic pressure tm1 pinned memory\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    if (cudaMallocHost(&(Global.mySolver->pressure_tm2),
                       Global.myMesh->nharbored * sizeof(solver_float)) !=
        cudaSuccess) {
      fprintf(stderr,
              "Thread %d: Failed to allocate acoustic pressure tm2 pinned memory\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    if (cudaMallocHost(&(Global.mySolver->pressure_force),
                       Global.myMesh->nharbored * sizeof(solver_float)) !=
        cudaSuccess) {
      fprintf(stderr,
              "Thread %d: Failed to allocate acoustic pressure force pinned memory\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    if (cudaMallocHost(&(Global.mySolver->pressure_mass),
                       Global.myMesh->nharbored * sizeof(solver_float)) !=
        cudaSuccess) {
      fprintf(stderr,
              "Thread %d: Failed to allocate acoustic pressure mass pinned memory\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    /* Step B: acoustic CPML nodal damping C. Allocated unconditionally
     * (acoustic path exists regardless of BOUNDARY/CPML); stays 0
     * everywhere unless an acoustic PML element contributes, so the
     * pressure update is bit-identical when there is no acoustic PML. */
    if (cudaMallocHost(&(Global.mySolver->pressure_damp),
                       Global.myMesh->nharbored * sizeof(solver_float)) !=
        cudaSuccess) {
      fprintf(stderr,
              "Thread %d: Failed to allocate acoustic pressure damp pinned memory\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    memset(Global.mySolver->pressure_tm1, 0,
           Global.myMesh->nharbored * sizeof(solver_float));
    memset(Global.mySolver->pressure_tm2, 0,
           Global.myMesh->nharbored * sizeof(solver_float));
    memset(Global.mySolver->pressure_force, 0,
           Global.myMesh->nharbored * sizeof(solver_float));
    memset(Global.mySolver->pressure_mass, 0,
           Global.myMesh->nharbored * sizeof(solver_float));
    memset(Global.mySolver->pressure_damp, 0,
           Global.myMesh->nharbored * sizeof(solver_float));

    pressureFreeSurfaceHost =
        (solver_float*)calloc(Global.myMesh->nharbored, sizeof(solver_float));
    if (pressureFreeSurfaceHost == NULL) {
      fprintf(stderr,
              "Thread %d: Failed to allocate acoustic pressure free-surface flags\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
  }

  if (Param.theTypeOfDamping >= BKT) {
    Global.mySolver->conv_shear_1 =
        (fvector_t*)calloc(8 * Global.myMesh->lenum, sizeof(fvector_t));
    Global.mySolver->conv_shear_2 =
        (fvector_t*)calloc(8 * Global.myMesh->lenum, sizeof(fvector_t));
    Global.mySolver->conv_kappa_1 =
        (fvector_t*)calloc(8 * Global.myMesh->lenum, sizeof(fvector_t));
    Global.mySolver->conv_kappa_2 =
        (fvector_t*)calloc(8 * Global.myMesh->lenum, sizeof(fvector_t));

    if ((Global.mySolver->conv_shear_1 == NULL) ||
        (Global.mySolver->conv_shear_2 == NULL) ||
        (Global.mySolver->conv_kappa_1 == NULL) ||
        (Global.mySolver->conv_kappa_2 == NULL)) {
      fprintf(stderr, "Thread %d: solver_init: out of memory\n", Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
  }

  Global.mySolver->dn_sched = schedule_new();
  Global.mySolver->an_sched = schedule_new();

  if (Param.theTypeOfDamping >= BKT3) {
    Global.mySolver->conv_shear_3 =
        (fvector_t*)calloc(8 * Global.myMesh->lenum, sizeof(fvector_t));
    Global.mySolver->conv_kappa_3 =
        (fvector_t*)calloc(8 * Global.myMesh->lenum, sizeof(fvector_t));

    if ((Global.mySolver->conv_shear_3 == NULL) ||
        (Global.mySolver->conv_kappa_3 == NULL)) {
      fprintf(stderr, "Thread %d: solver_init: out of memory\n", Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
  }

  if ((Global.mySolver->eTable == NULL) || (Global.mySolver->nTable == NULL) ||
      (Global.mySolver->tm1 == NULL) || (Global.mySolver->tm2 == NULL) ||
      (Global.mySolver->force == NULL)) {
    fprintf(stderr, "Thread %d: solver_init: out of memory\n", Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  if (Param.printStationAccelerations == YES) {
    Global.mySolver->tm3 =
        (fvector_t*)calloc(Global.myMesh->nharbored, sizeof(fvector_t));

    if (Global.mySolver->tm3 == NULL) {
      fprintf(stderr, "Thread %d: solver_init: out of memory for accs\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
  }

#ifdef BOUNDARY
  /* CPML build-time scratch. Sized to the upper bound (lenum) per type;
   * only the first elastic_count / acoustic_count entries are uploaded. */
  if (Param.theAbcType == ABC_CPML) {
    g_cpmlBuild.elastic_count = 0;
    g_cpmlBuild.acoustic_count = 0;

    /* Publish the elastic PML stiffness mode to the device (read by
     * kernelStiffnessCalcLocal's PML skip and kernelPMLElasticStep). */
    Global.gpuData.pmlNativeStiffness =
        (Param.theCpmlNativeStiffness == YES) ? 1 : 0;
    Global.gpuData.pmlNativeStiffnessAcoustic =
        (Param.theCpmlNativeStiffnessAcoustic == YES) ? 1 : 0;
    /* Stage-1 plumbing: 0 none (default; PML material-damping skip
     * unchanged, byte-identical baseline). Consumed from Stage 2 on. */
    Global.gpuData.cpmlPmlMaterialDamping = Param.theCpmlPmlMaterialDamping;

    /* Resolve PML thickness. Priority: explicit cpml_thickness_m wins; else
     * cpml_layers × (max edgesize among elements ON an absorbing face).
     *
     * Octree meshes refine by local Vs, so the global max element size is
     * usually a coarse *interior* element (deep fast rock) that is NOT
     * representative of the boundary. Sizing the PML from the global max
     * over-thickens it and pushes the inner edge deep into the model.
     * The PML thickness should instead span cpml_layers element layers
     * AT THE BOUNDARY, so the max is taken only over elements touching
     * one of the five absorbing faces (x_min, x_max, y_min, y_max,
     * z_max); z_min is the free surface and is excluded. The MPI max-
     * reduce keeps the realized thickness identical on every rank, which
     * face_mask (computed independently per element) requires for shared
     * corners across partitions. (Mirrors the v3-branch behaviour.) */
    if (Param.theCpmlThicknessM > 0.0) {
      g_cpmlBuild.thickness_m = Param.theCpmlThicknessM;
    } else {
      const double ticksize = Global.myOctree->ticksize;
      const double eps = ticksize * 0.25;
      const double lo[3] = {
          (double)Global.myOctree->nearendp[0] * ticksize,
          (double)Global.myOctree->nearendp[1] * ticksize,
          (double)Global.myOctree->nearendp[2] * ticksize};
      const double hi[3] = {
          (double)Global.myOctree->farendp[0] * ticksize,
          (double)Global.myOctree->farendp[1] * ticksize,
          (double)Global.myOctree->farendp[2] * ticksize};
      double local_bnd_h = 0.0;   /* max edgesize on an absorbing face */
      double local_any_h = 0.0;   /* global max (defensive fallback)   */
      for (int32_t ei = 0; ei < Global.myMesh->lenum; ei++) {
        elem_t* ep = &Global.myMesh->elemTable[ei];
        edata_t* ed = (edata_t*)ep->data;
        double h = (double)ed->edgesize;
        if (h <= 0.0) continue;
        if (h > local_any_h) local_any_h = h;
        int32_t lnid0 = ep->lnid[0];
        const double c0[3] = {
            (double)Global.myMesh->nodeTable[lnid0].x * ticksize,
            (double)Global.myMesh->nodeTable[lnid0].y * ticksize,
            (double)Global.myMesh->nodeTable[lnid0].z * ticksize};
        int on_face = 0;
        /* x, y: both min and max faces are absorbing. */
        for (int ax = 0; ax < 2; ax++) {
          if ((c0[ax] - lo[ax]) <= eps) on_face = 1;
          if ((hi[ax] - (c0[ax] + h)) <= eps) on_face = 1;
        }
        /* z: only z_max (deep/far) absorbs; z_min is the free surface. */
        if ((hi[2] - (c0[2] + h)) <= eps) on_face = 1;
        if (on_face && h > local_bnd_h) local_bnd_h = h;
      }
      double global_bnd_h = 0.0, global_any_h = 0.0;
      MPI_Allreduce(&local_bnd_h, &global_bnd_h, 1, MPI_DOUBLE, MPI_MAX,
                    comm_solver);
      MPI_Allreduce(&local_any_h, &global_any_h, 1, MPI_DOUBLE, MPI_MAX,
                    comm_solver);
      /* Defensive: a non-empty domain always has boundary elements, but
       * if the boundary scan somehow found none, fall back to the global
       * element max rather than abort with a zero-thickness PML. */
      int used_fallback = (global_bnd_h <= 0.0);
      double basis_h = used_fallback ? global_any_h : global_bnd_h;
      g_cpmlBuild.thickness_m = (double)Param.theCpmlLayers * basis_h;
      /* Cache the realized thickness so the log line and the MPI broadcast
       * (if a checkpoint or replay needs it later) reflect what was used. */
      Param.theCpmlThicknessM = g_cpmlBuild.thickness_m;
      if (Global.myID == 0) {
        fprintf(stdout,
                "CPML: derived thickness %.3f m from %d layers x %.3f m "
                "(%s max element edgesize)\n",
                g_cpmlBuild.thickness_m, Param.theCpmlLayers, basis_h,
                used_fallback ? "global-fallback" : "boundary");
      }
    }

    /* Default alpha_max = pi * f_max if user didn't specify */
    g_cpmlBuild.alpha_max = (Param.theCpmlAlphaMaxHz > 0.0)
                                ? (2.0 * M_PI * Param.theCpmlAlphaMaxHz)
                                : (M_PI * Param.theFreq);
    g_cpmlBuild.face_mask =
        (uint8_t*)calloc(Global.myMesh->lenum, sizeof(uint8_t));
    g_cpmlBuild.elastic_elem_list =
        (int32_t*)malloc(Global.myMesh->lenum * sizeof(int32_t));
    g_cpmlBuild.acoustic_elem_list =
        (int32_t*)malloc(Global.myMesh->lenum * sizeof(int32_t));
    g_cpmlBuild.elastic_muCoeff =
        (double*)malloc(Global.myMesh->lenum * sizeof(double));
    g_cpmlBuild.elastic_lambdaCoeff =
        (double*)malloc(Global.myMesh->lenum * sizeof(double));
    g_cpmlBuild.elastic_restoreCoeff =
        (double*)malloc(Global.myMesh->lenum * sizeof(double));
    g_cpmlBuild.acoustic_forceCoeff =
        (double*)malloc(Global.myMesh->lenum * sizeof(double));
    g_cpmlBuild.acoustic_restoreCoeff =
        (double*)malloc(Global.myMesh->lenum * sizeof(double));
    for (int axis = 0; axis < 3; axis++) {
      g_cpmlBuild.elastic_a[axis] =
          (double*)malloc(Global.myMesh->lenum * sizeof(double));
      g_cpmlBuild.elastic_b[axis] =
          (double*)malloc(Global.myMesh->lenum * sizeof(double));
      g_cpmlBuild.elastic_invKappa[axis] =
          (double*)malloc(Global.myMesh->lenum * sizeof(double));
      g_cpmlBuild.elastic_zetaDecay[axis] =
          (double*)malloc(Global.myMesh->lenum * sizeof(double));
      g_cpmlBuild.elastic_zetaGain[axis] =
          (double*)malloc(Global.myMesh->lenum * sizeof(double));
      g_cpmlBuild.elastic_memForce[axis] =
          (double*)malloc(Global.myMesh->lenum * sizeof(double));
      g_cpmlBuild.acoustic_a[axis] =
          (double*)malloc(Global.myMesh->lenum * sizeof(double));
      g_cpmlBuild.acoustic_b[axis] =
          (double*)malloc(Global.myMesh->lenum * sizeof(double));
      g_cpmlBuild.acoustic_invKappa[axis] =
          (double*)malloc(Global.myMesh->lenum * sizeof(double));
      g_cpmlBuild.acoustic_zetaDecay[axis] =
          (double*)malloc(Global.myMesh->lenum * sizeof(double));
      g_cpmlBuild.acoustic_zetaGain[axis] =
          (double*)malloc(Global.myMesh->lenum * sizeof(double));
      g_cpmlBuild.acoustic_memForce[axis] =
          (double*)malloc(Global.myMesh->lenum * sizeof(double));
    }
    if ((g_cpmlBuild.face_mask == NULL) ||
        (g_cpmlBuild.elastic_elem_list == NULL) ||
        (g_cpmlBuild.acoustic_elem_list == NULL) ||
        (g_cpmlBuild.elastic_a[0] == NULL) ||
        (g_cpmlBuild.elastic_b[0] == NULL) ||
        (g_cpmlBuild.elastic_invKappa[0] == NULL) ||
        (g_cpmlBuild.elastic_muCoeff == NULL) ||
        (g_cpmlBuild.elastic_lambdaCoeff == NULL) ||
        (g_cpmlBuild.elastic_zetaDecay[0] == NULL) ||
        (g_cpmlBuild.elastic_zetaGain[0] == NULL) ||
        (g_cpmlBuild.elastic_memForce[0] == NULL) ||
        (g_cpmlBuild.elastic_restoreCoeff == NULL) ||
        (g_cpmlBuild.acoustic_a[0] == NULL) ||
        (g_cpmlBuild.acoustic_b[0] == NULL) ||
        (g_cpmlBuild.acoustic_invKappa[0] == NULL) ||
        (g_cpmlBuild.acoustic_forceCoeff == NULL) ||
        (g_cpmlBuild.acoustic_zetaDecay[0] == NULL) ||
        (g_cpmlBuild.acoustic_zetaGain[0] == NULL) ||
        (g_cpmlBuild.acoustic_memForce[0] == NULL) ||
        (g_cpmlBuild.acoustic_restoreCoeff == NULL)) {
      fprintf(stderr,
              "Thread %d: solver_init: out of memory for CPML build state\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
  }
#endif /* BOUNDARY */

  /* For each element:
   * Initialize the data structures. tm1, tm2 and force have
   * already been initialized to 0 by calloc(). */
  for (eindex = 0; eindex < Global.myMesh->lenum; eindex++) {
    elem_t* elemp; /* pointer to the mesh database */
    edata_t* edata;
    e_t* ep; /* pointer to the element constant table */
    double mass, M, mu, lambda;
    int j;

    /* Product of kappa over the element's active stretch axes
     * (kappa_x*kappa_y*kappa_z; the two transverse factors are 1 for a
     * face element). The matched kappa-scaled inertia: the gradient side
     * stretches strain by 1/kappa, so the second-order displacement
     * CFS-CPML requires rho*kappa_prod*d2u/dt2 (leading term of
     * multiplying the momentum eq by s_x s_y s_z). 1.0 for non-PML and
     * dashpot (exact no-op; kappa_max=1 -> kappa_prod=1 everywhere, so
     * runs that never raised kappa are bit-identical). Elastic CPML
     * only. Scaling inertia up also raises the local CFL margin.
     * Declared OUTSIDE #ifdef BOUNDARY because it scales the
     * unconditional nodal inertia accumulation below. */
    double cpml_kappa_prod = 1.0;

#ifdef BOUNDARY
    tick_t edgeticks;
    tick_t ldb[3], ruf[3];
    int32_t lnid0;
    char flag;
    double dashpot[8][3];
    uint8_t pml_face_mask = 0;
    /* Summed CPML damping rate sum_j d_j (1/s) for this element, used to
     * build the Matzen second-order divergence-side mass-proportional
     * damping rho*d*du/dt. Non-zero only for elastic CPML elements. */
    double cpml_damp_rate = 0.0;
#endif

    double zeta, a, b;

    /* yigit says -- Buildings assume 5% damping */

    int32_t n_0;
    double z_m;

    /* Note the difference between the two tables */
    elemp = &Global.myMesh->elemTable[eindex];
    edata = (edata_t*)elemp->data;
    ep = &Global.mySolver->eTable[eindex];
    acoustic_sanitize_element(edata);

    if (acoustic_pressure_enabled() && element_is_acoustic(edata)) {
      const double bulk = edata->rho * edata->Vp * edata->Vp;
      const double volume =
          edata->edgesize * edata->edgesize * edata->edgesize;
      const double nodalPressureMass = volume / (8.0 * bulk);

      ep->c1 = 0.0;
      ep->c2 = 0.0;
      ep->c3 = 0.0;
      ep->c4 = 0.0;

      /* Step B divergence-side scalars. Defaults are a STRICT no-op:
       * acoustic_kprod = 1 (pressure mass unscaled, bit-identical) and
       * acoustic_damp = 0 (no damping). They are set only for an
       * ABC_CPML acoustic element inside the PML zone. Declared here
       * (outside #ifdef BOUNDARY) because they scale the unconditional
       * nodal pressure-mass / pressure-damp accumulation below. */
      double acoustic_kprod = 1.0;
      double acoustic_damp = 0.0;

#ifdef BOUNDARY
      /* Step A activation + Step B divergence-side coefficient build.
       * Acoustic elements 'continue' below before the elastic-path CPML
       * block, so this is where acoustic PML is enrolled (Step A) and
       * its divergence-side coeffs are built (Step B). Relocated ABOVE
       * the pressure-mass accumulation so the matched kappa-inertia can
       * scale it. Strict no-op for dashpot / non-PML / non-BOUNDARY ->
       * every previously-validated configuration is bit-identical. */
      if (Param.theAbcType == ABC_CPML) {
        lnid0 = elemp->lnid[0];
        ldb[0] = Global.myMesh->nodeTable[lnid0].x;
        ldb[1] = Global.myMesh->nodeTable[lnid0].y;
        ldb[2] = Global.myMesh->nodeTable[lnid0].z;
        edgeticks = (tick_t)1 << (PIXELLEVEL - elemp->level);
        ruf[0] = ldb[0] + edgeticks;
        ruf[1] = ldb[1] + edgeticks;
        ruf[2] = ldb[2] + edgeticks;
        pml_face_mask = compute_pml_face_mask(
            ldb, ruf, Global.myOctree->nearendp, Global.myOctree->farendp,
            Global.myOctree->ticksize, g_cpmlBuild.thickness_m);
        g_cpmlBuild.face_mask[eindex] = pml_face_mask;
        if (pml_face_mask != 0) {
          double a_axis[3], b_axis[3], invKappa_axis[3], d_axis[3];
          double alpha_axis[3], kappa_axis[3];
          compute_pml_profile_coeffs(
              pml_face_mask, ldb, ruf, Global.myOctree->nearendp,
              Global.myOctree->farendp, Global.myOctree->ticksize,
              g_cpmlBuild.thickness_m, (double)edata->Vp, Param.theDeltaT,
              Param.theCpmlPolyOrder, Param.theCpmlReflectionTarget,
              Param.theCpmlKappaMax, g_cpmlBuild.alpha_max,
              Param.theCpmlMpmlCrossFraction,
              Param.theCpmlAlphaFloorFraction,
              a_axis, b_axis, invKappa_axis, d_axis, alpha_axis, kappa_axis);
          const double h = (double)edata->edgesize;
          const double dt = Param.theDeltaT;
          const double dt2 = Param.theDeltaTSquared;
          int32_t pml_idx = g_cpmlBuild.acoustic_count++;
          g_cpmlBuild.acoustic_elem_list[pml_idx] = eindex;
          /* Gradient side -- store RAW profile coeffs (natural signed-sum
           * normalization), mirroring the validated elastic native+
           * correction (elastic_a/b/invKappa = a/b/invKappa raw, the
           * 1/(4h) folded into the single stiffness coeff below = the
           * corrected consistent interior coeff dt^2*h/(16 rho), exactly
           * as elastic folds it into mu=9*c1). kernelPMLAcousticStep then
           * forms the native delta -(P_pml - P_raw), zero at d=0 -> exact
           * match to the hourglass-free consistent interior operator. */
          for (int axis = 0; axis < 3; axis++) {
            g_cpmlBuild.acoustic_a[axis][pml_idx] = a_axis[axis];
            g_cpmlBuild.acoustic_b[axis][pml_idx] = b_axis[axis];
            g_cpmlBuild.acoustic_invKappa[axis][pml_idx] =
                invKappa_axis[axis];
          }
          const double rho = (double)edata->rho;
          /* Cstiff = dt^2 * h / (16 rho): the SAME corrected consistent
           * interior acoustic stiffness coefficient (commit 20180b6), so
           * P_raw exactly equals the constant-gradient content of the
           * interior operator and the native delta vanishes at the PML
           * inner edge. rho > 0 guaranteed by mesh validation. (Array
           * name kept; semantics are now the stiffness coeff, not the
           * old dt^2 h^2/4 forceCoeff.) */
          g_cpmlBuild.acoustic_forceCoeff[pml_idx] = dt2 * h / (16.0 * rho);

          /* Step B divergence side -- scalar mirror of the validated
           * elastic two-sided. memForce/restoreCoeff use the acoustic
           * nodal mass (nodalPressureMass) in place of the elastic
           * rho*h^3/8. */
          acoustic_kprod =
              kappa_axis[0] * kappa_axis[1] * kappa_axis[2];
          double restore = 0.0;
          for (int axis = 0; axis < 3; axis++) {
            const double da = d_axis[axis];
            const double aa = alpha_axis[axis];
            const double decay = exp(-aa * dt);
            g_cpmlBuild.acoustic_zetaDecay[axis][pml_idx] = decay;
            g_cpmlBuild.acoustic_zetaGain[axis][pml_idx] =
                (aa > 0.0) ? ((1.0 - decay) / aa) : dt;
            g_cpmlBuild.acoustic_memForce[axis][pml_idx] =
                dt2 * nodalPressureMass * da * aa * aa;
            restore += da * aa;
          }
          g_cpmlBuild.acoustic_restoreCoeff[pml_idx] =
              dt2 * nodalPressureMass * restore;
          acoustic_damp = d_axis[0] + d_axis[1] + d_axis[2];
        }
      }
#endif /* BOUNDARY */

      for (j = 0; j < 8; j++) {
        const int32_t lnid = elemp->lnid[j];
        Global.mySolver->pressure_mass[lnid] +=
            (solver_float)(nodalPressureMass * acoustic_kprod);
        /* Bake dt into pressure_damp so kernelAcousticPressureUpdate
         * needs no dt: stored value is the assembled dt*C term, exactly
         * as the elastic path bakes dt*C into mass_minusaM/mass2. */
        Global.mySolver->pressure_damp[lnid] +=
            (solver_float)(nodalPressureMass * acoustic_damp *
                           Param.theDeltaT);
      }
      if (Global.myMesh->nodeTable[elemp->lnid[0]].z ==
          Global.myOctree->nearendp[2]) {
        for (j = 0; j < 4; j++) {
          pressureFreeSurfaceHost[elemp->lnid[j]] = 1.0;
        }
      }

      localPressureElements++;
      continue;
    }

    /* Calculate the Lame constants */
    mu_and_lambda(&mu, &lambda, edata, eindex);

    /* coefficients for term (deltaT_squared * Ke * Ut) */
    ep->c1 = Param.theDeltaTSquared * edata->edgesize * mu / 9;
    ep->c2 = Param.theDeltaTSquared * edata->edgesize * lambda / 9;

    /* coefficients for term (b * deltaT * Ke_off * (Ut-1 - Ut)) */
    /* Anelastic attenuation (material damping) */

    /* Old formula for damping */
    /* zeta = (edata->Vs < 1500) ? 25 / edata->Vs : 5 / edata->Vs; */

    n_0 = Global.myMesh->elemTable[eindex].lnid[0];
    z_m = Global.theZForMeshOrigin +
          (Global.myOctree->ticksize) * Global.myMesh->nodeTable[n_0].z;

    /* Shift the domain if buildings are considered */
    if (Param.includeBuildings == YES) {
      z_m -= get_surface_shift();
    }

    /* New formula for damping according to Graves. Acoustic elements do not
     * use the shear-velocity damping law. */
    zeta = element_is_acoustic(edata) ? 0.0 : (25.0 / edata->Vs);

    /* If element is in a building, use 5% damping */
    if ((z_m < 0) &&
        (Param.includeBuildings ==
         YES)) { /* Dorian says: This is only valid for buildings  */
      zeta = 0.05;
    }

    if (zeta > Param.theThresholdDamping) {
      zeta = Param.theThresholdDamping;
    }

    /* the a,b coefficients */
    a = zeta * Global.theABase;
    b = zeta * Global.theBBase;

    /* coefficients for term (b * deltaT * Ke_off * (Ut-1 - Ut)) */
    ep->c3 = b * Param.theDeltaT * edata->edgesize * mu / 9;
    ep->c4 = b * Param.theDeltaT * edata->edgesize * lambda / 9;

#ifdef BOUNDARY

    /* Set the flag for the element */
    lnid0 = elemp->lnid[0];

    ldb[0] = Global.myMesh->nodeTable[lnid0].x;
    ldb[1] = Global.myMesh->nodeTable[lnid0].y;
    ldb[2] = Global.myMesh->nodeTable[lnid0].z;

    edgeticks = (tick_t)1 << (PIXELLEVEL - elemp->level);
    ruf[0] = ldb[0] + edgeticks;
    ruf[1] = ldb[1] + edgeticks;
    ruf[2] = ldb[2] + edgeticks;

    flag = compute_setflag(ldb, ruf, Global.myOctree->nearendp,
                           Global.myOctree->farendp);

    if (Param.theAbcType == ABC_CPML) {
      pml_face_mask = compute_pml_face_mask(
          ldb, ruf, Global.myOctree->nearendp, Global.myOctree->farendp,
          Global.myOctree->ticksize, g_cpmlBuild.thickness_m);
      g_cpmlBuild.face_mask[eindex] = pml_face_mask;
      if (pml_face_mask != 0) {
        double a_axis[3], b_axis[3], invKappa_axis[3], d_axis[3];
        double alpha_axis[3], kappa_axis[3];
        compute_pml_profile_coeffs(
            pml_face_mask, ldb, ruf, Global.myOctree->nearendp,
            Global.myOctree->farendp, Global.myOctree->ticksize,
            g_cpmlBuild.thickness_m, (double)edata->Vp, Param.theDeltaT,
            Param.theCpmlPolyOrder, Param.theCpmlReflectionTarget,
            Param.theCpmlKappaMax, g_cpmlBuild.alpha_max,
            Param.theCpmlMpmlCrossFraction,
            Param.theCpmlAlphaFloorFraction,
            a_axis, b_axis, invKappa_axis, d_axis, alpha_axis, kappa_axis);

        /* Native (c1,c2-style) normalization: store RAW profile a/b/
         * invKappa; the physical scaling is folded into a single
         * stiffness coeff per type (elastic: mu=9*c1, lambda=9*c2;
         * acoustic: Cstiff = dt^2*h/(16 rho)). h2_over_4 is still used
         * by the elastic muCoeff/lambdaCoeff below. */
        const double h = (double)edata->edgesize;
        const double dt2 = Param.theDeltaTSquared;
        const double h2_over_4 = h * h * 0.25;

        if (element_is_acoustic(edata)) {
          int32_t pml_idx = g_cpmlBuild.acoustic_count++;
          g_cpmlBuild.acoustic_elem_list[pml_idx] = eindex;
          /* RAW profile coeffs + Cstiff = dt^2*h/(16 rho): native+
           * correction mirror of the elastic path (see the matching
           * acoustic block above and kernelPMLAcousticStep). At d=0 the
           * PML delta -(P_pml - P_raw) is exactly 0 -> the element is
           * purely the validated consistent interior operator. */
          for (int axis = 0; axis < 3; axis++) {
            g_cpmlBuild.acoustic_a[axis][pml_idx] = a_axis[axis];
            g_cpmlBuild.acoustic_b[axis][pml_idx] = b_axis[axis];
            g_cpmlBuild.acoustic_invKappa[axis][pml_idx] =
                invKappa_axis[axis];
          }
          const double rho = (double)edata->rho;
          /* rho > 0 guaranteed by mesh validation. (Array name kept;
           * semantics now the stiffness coeff, not the old dt^2 h^2/4. */
          g_cpmlBuild.acoustic_forceCoeff[pml_idx] = dt2 * h / (16.0 * rho);
        } else {
          int32_t pml_idx = g_cpmlBuild.elastic_count++;
          g_cpmlBuild.elastic_elem_list[pml_idx] = eindex;
          for (int axis = 0; axis < 3; axis++) {
            /* RAW a / 1/kappa (no 1/(4h) prescale). The elastic PML now
             * runs as native effective stiffness + a CPML correction
             * evaluated in the native (c1,c2) normalization with
             * dN = 0.25*sign centre derivatives, so the psi recurrence
             * convolves the natural-coord gradient with the raw a, and
             * the modified gradient uses the raw 1/kappa -- matching the
             * codex v3 reduced-correction (which is proven self-
             * consistent with aTransposeU/firstVector/au). b is the
             * dimensionless exp() decay, unchanged. (Acoustic PML still
             * uses the 1-pt /(4h) form above.) */
            g_cpmlBuild.elastic_a[axis][pml_idx] = a_axis[axis];
            g_cpmlBuild.elastic_b[axis][pml_idx] = b_axis[axis];
            g_cpmlBuild.elastic_invKappa[axis][pml_idx] =
                invKappa_axis[axis];
          }
          /* muCoeff/lambdaCoeff are no longer read by kernelPMLElasticStep
           * (it uses mu = 9*c1, lambda = 9*c2 like the native operator);
           * kept allocated for now, slated for cleanup. */
          g_cpmlBuild.elastic_muCoeff[pml_idx] = dt2 * h2_over_4 * mu;
          g_cpmlBuild.elastic_lambdaCoeff[pml_idx] = dt2 * h2_over_4 * lambda;
          /* Matzen second-order divergence-side term: total volumetric
           * CPML damping rate sum_j d_j (1/s). Injected below as the
           * mass-proportional damping rho*d*du/dt at this element's
           * nodes, exactly mirroring the Lysmer dashpot discretization
           * (mass_minusaM / mass2_minusaM -= dt * C). Elastic only;
           * the acoustic-pressure update path is left gradient-side
           * only for now (documented asymmetry). */
          cpml_damp_rate = d_axis[0] + d_axis[1] + d_axis[2];
          /* Matched kappa-scaled inertia (rho*kappa_prod*d2u/dt2). For a
           * face element the two transverse kappa are 1; at the inner
           * edge (xi=0) all kappa=1 so kappa_prod=1 -> reduces to the
           * native operator (continuous interface). kappa_max=1 ->
           * kappa_prod=1 everywhere (exact no-op, back-compat). */
          cpml_kappa_prod =
              kappa_axis[0] * kappa_axis[1] * kappa_axis[2];

          /* CFS sub-terms (Matzen 2011; Xie et al. 2014). Per damping
           * axis a the divergence-side stretch also contributes
           *   +rho*d_a*alpha_a u_i           (restoring)
           *   -rho*d_a*alpha_a^2 zeta_i^(a)  (displacement memory)
           * with zeta the exact exponential convolution of the element
           * displacement, d zeta/dt = -alpha_a zeta + u:
           *   zetaDecay_a = exp(-alpha_a*dt)
           *   zetaGain_a  = (1 - zetaDecay_a)/alpha_a   (alpha_a>0: floor)
           *   memForce_a  = dt^2 * M_node * d_a * alpha_a^2
           *   restoreCoeff= dt^2 * M_node * sum_a d_a*alpha_a
           * M_node = rho*h^3/8 (per-node share). BOTH the restoring and
           * the memory are applied in kernelPMLElasticStep from the
           * element-centroid displacement so the CFS pair cancels
           * exactly at DC discretely. (v10 applied the restoring per
           * node via mass2_minusaM while zeta was centroid-based; the
           * mismatch left a net anti-restoring forcing that grew
           * exponentially at PML corners. Using ubar for both halves
           * removes that net injection.) */
          const double dt = Param.theDeltaT;
          const double rho_e = (double)edata->rho;
          const double M_node = rho_e * h * h * h / 8.0;
          double restore = 0.0;
          for (int axis = 0; axis < 3; axis++) {
            const double da = d_axis[axis];
            const double aa = alpha_axis[axis];
            const double decay = exp(-aa * dt);
            g_cpmlBuild.elastic_zetaDecay[axis][pml_idx] = decay;
            g_cpmlBuild.elastic_zetaGain[axis][pml_idx] =
                (aa > 0.0) ? ((1.0 - decay) / aa) : dt;
            g_cpmlBuild.elastic_memForce[axis][pml_idx] =
                dt2 * M_node * da * aa * aa;
            restore += da * aa;
          }
          g_cpmlBuild.elastic_restoreCoeff[pml_idx] = dt2 * M_node * restore;
        }
      }
    }

    /* Compute the Lysmer dashpot contribution unless CPML claims this
     * element (i.e., pml_face_mask != 0 in CPML mode). Interior elements
     * (flag == 13) always skip; that is unchanged. */
    if ((flag != 13) &&
        ((Param.theAbcType == ABC_DASHPOT) || (pml_face_mask == 0))) {
      compute_setboundary(edata->edgesize, edata->Vp, edata->Vs, edata->rho,
                          flag, dashpot);
    } else {
      memset(dashpot, 0, sizeof(dashpot));
    }
#endif /* BOUNDARY */

    /* Assign the element mass to its vertices */
    /* mass is the total mass of the element   */
    /* and M is the mass assigned to each node */
    mass = edata->rho * edata->edgesize * edata->edgesize * edata->edgesize;
    M = mass / 8;

    /* Essential for topography module:   */
    /* Here I correct the nodal mass in the case that topography exists  */
    double nodesMass[8] = {0};

    /* get element coordinates */
    double xo = Global.myMesh->ticksize * Global.myMesh->nodeTable[lnid0].x;
    double yo = Global.myMesh->ticksize * Global.myMesh->nodeTable[lnid0].y;
    double zo = Global.myMesh->ticksize * Global.myMesh->nodeTable[lnid0].z;

    if (Param.includeTopography == YES) {
      toponodes_mass(eindex, nodesMass, M, xo, yo, zo);
    } else {
      for (j = 0; j < 8; j++) {
        nodesMass[j] = M;
      }
    }

    /* For each node */
    for (j = 0; j < 8; j++) {
      int32_t lnid;
      int axis;
      n_t* np;

      lnid = elemp->lnid[j];
      np = &Global.mySolver->nTable[lnid];

      M = nodesMass[j];

      /* Matched kappa-scaled inertia (rho*kappa_prod*d2u/dt2). Only the
       * INERTIA contributions (mass_simple and the +M / +2M below) are
       * scaled; the Rayleigh/dashpot/CPML-damping terms stay on plain M
       * (they come from the d-part of s, not the kappa-part). For
       * non-PML / dashpot / kappa_max=1, cpml_kappa_prod == 1.0 -> exact
       * no-op (bit-identical). */
      const double Mk = M * cpml_kappa_prod;

      np->mass_simple += Mk;

      /* Material mass-proportional Rayleigh damping a*M. In CPML mode
       * this is now skipped on PML elements, CONSISTENTLY with the BKT/Q,
       * stiffness-proportional (b*K), and Lysmer-dashpot material-damping
       * skips already applied there: the CPML ADE recurrence + the
       * divergence-side mass-proportional damping/CFS provide the PML
       * elements' attenuation, and layering physical Rayleigh damping
       * into the same mass_minusaM/mass2_minusaM arrays double-counts
       * decay and perturbs the delicate divergence-side/CFS DC balance.
       * Predicate mirrors the adjacent dashpot skip exactly (minus the
       * dashpot-only flag!=13). Dashpot mode, non-PML elements, and
       * BOUNDARY-undefined builds are bit-identical; only attenuating-
       * CPML PML elements change (and become consistent/correct).
       * Acoustic elements have zeta=0 (a=0) so this is an elastic-PML
       * concern; no-attenuation runs (a~0) are unaffected. */
      int apply_material_mass_damping = 1;
#ifdef BOUNDARY
      if ((Param.theAbcType != ABC_DASHPOT) && (pml_face_mask != 0)) {
        apply_material_mass_damping = 0;
      }
#endif

      /* loop for each axis */
      for (axis = 0; axis < 3; axis++) {
        if (apply_material_mass_damping) {
          np->mass_minusaM[axis] -= (Param.theDeltaT * a * M);
          np->mass2_minusaM[axis] -= (Param.theDeltaT * a * M);
        }

#ifdef BOUNDARY
        if ((flag != 13) &&
            ((Param.theAbcType == ABC_DASHPOT) || (pml_face_mask == 0))) {
          /* boundary impact (Lysmer dashpot) — skipped for CPML elements */
          np->mass_minusaM[axis] -= (Param.theDeltaT * dashpot[j][axis]);
          np->mass2_minusaM[axis] -= (Param.theDeltaT * dashpot[j][axis]);
        }

        /* Matzen 2011 / Xie et al. 2014 second-order displacement CPML:
         * the divergence-side coordinate stretch, obtained by multiplying
         * the momentum equation by s_j, contributes (to leading order,
         * kappa_max = 1 default) a mass-proportional volumetric damping
         *     rho * d_j * du_i/dt
         * Per node the damping coefficient is C = (rho * V_node) * sum_j d_j
         * = M * cpml_damp_rate. It is discretized exactly like the Lysmer
         * dashpot above (same backward-difference, mass_minusaM and
         * mass2_minusaM each -= dt*C), so it inherits the dashpot path's
         * stability. cpml_damp_rate > 0 only for elastic CPML elements;
         * interior and acoustic elements are untouched. The gradient/
         * strain side is handled separately in kernelPMLElasticStep via
         * the psi memory (modified strain); the CFS restoring/memory
         * sub-terms are added just below + in the kernel. */
        if ((Param.theAbcType == ABC_CPML) && (cpml_damp_rate > 0.0)) {
          double cpml_C = M * cpml_damp_rate;
          np->mass_minusaM[axis]  -= (Param.theDeltaT * cpml_C);
          np->mass2_minusaM[axis] -= (Param.theDeltaT * cpml_C);
        }

        /* (CFS restoring is NOT applied here. v10 added it as a per-node
         * mass2_minusaM u^n coefficient while the paired zeta memory was
         * element-centroid based; the discrete DC cancellation failed
         * and corners grew exponentially. Both CFS halves are now
         * evaluated together from the centroid displacement in
         * kernelPMLElasticStep.) */
#endif /* BOUNDARY */

        np->mass_minusaM[axis] += Mk;
        np->mass2_minusaM[axis] += (Mk * 2);

      } /* end loop for each axis */

    } /* end loop for each node */

  } /* eindex for elements */

#ifdef BOUNDARY
  /* CPML post-loop: publish per-type counts, copy elem_lists to mySolver
   * (for checkpoint I/O), allocate GPU coefficient/psi arrays per type,
   * and upload. The GPU psi arrays are zero-initialized. Build-time host
   * scratch is freed after upload. */
  if (Param.theAbcType == ABC_CPML) {
    int64_t local_counts[2] = {(int64_t)g_cpmlBuild.elastic_count,
                               (int64_t)g_cpmlBuild.acoustic_count};
    int64_t global_counts[2] = {0, 0};
    MPI_Reduce(local_counts, global_counts, 2, MPI_LONG_LONG, MPI_SUM, 0,
               comm_solver);
    if (Global.myID == 0) {
      fprintf(stdout,
              "CPML: %lld elastic + %lld acoustic PML elements globally "
              "(thickness %.3f m, alpha_max = %.3g rad/s)\n",
              (long long)global_counts[0], (long long)global_counts[1],
              g_cpmlBuild.thickness_m, g_cpmlBuild.alpha_max);
      if (Param.theMonitorFileFp != NULL) {
        fprintf(Param.theMonitorFileFp,
                "CPML: %lld elastic + %lld acoustic PML elements globally, "
                "thickness %.3f m, alpha_max %.3g rad/s\n",
                (long long)global_counts[0], (long long)global_counts[1],
                g_cpmlBuild.thickness_m, g_cpmlBuild.alpha_max);
      }
    }

    /* --- Elastic PML upload --- */
    {
      int32_t pmlN = g_cpmlBuild.elastic_count;
      Global.gpuData.pmlElasticElemCount = pmlN;
      Global.mySolver->pmlElasticElemCount = pmlN;
      if (pmlN > 0) {
        Global.mySolver->pmlElasticElemList =
            (int32_t*)malloc((size_t)pmlN * sizeof(int32_t));
        if (Global.mySolver->pmlElasticElemList == NULL) {
          fprintf(stderr,
                  "Thread %d: solver_init: oom pmlElasticElemList\n",
                  Global.myID);
          MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
        }
        memcpy(Global.mySolver->pmlElasticElemList,
               g_cpmlBuild.elastic_elem_list,
               (size_t)pmlN * sizeof(int32_t));

        for (int i = 0; i < 3; i++) {
          for (int jx = 0; jx < 3; jx++) {
            Global.mySolver->pml_elastic_psi[i][jx] = (solver_float*)calloc(
                (size_t)pmlN, sizeof(solver_float));
            if (Global.mySolver->pml_elastic_psi[i][jx] == NULL) {
              fprintf(stderr,
                      "Thread %d: solver_init: oom pml_elastic_psi[%d][%d]\n",
                      Global.myID, i, jx);
              MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
            }
          }
        }

        if (cudaMalloc((void**)&(Global.gpuData.pmlElasticElemMapperDevice),
                       (size_t)pmlN * sizeof(int32_t)) != cudaSuccess) {
          fprintf(stderr, "Thread %d: cudaMalloc pmlElasticElemMapper failed\n",
                  Global.myID);
          MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
        }
        double** dev_a[3]   = {&Global.gpuData.pml_elastic_a_xArrayDevice,
                               &Global.gpuData.pml_elastic_a_yArrayDevice,
                               &Global.gpuData.pml_elastic_a_zArrayDevice};
        double** dev_b[3]   = {&Global.gpuData.pml_elastic_b_xArrayDevice,
                               &Global.gpuData.pml_elastic_b_yArrayDevice,
                               &Global.gpuData.pml_elastic_b_zArrayDevice};
        double** dev_invK[3]= {&Global.gpuData.pml_elastic_invKappa_xArrayDevice,
                               &Global.gpuData.pml_elastic_invKappa_yArrayDevice,
                               &Global.gpuData.pml_elastic_invKappa_zArrayDevice};
        for (int axis = 0; axis < 3; axis++) {
          if (cudaMalloc((void**)dev_a[axis],
                         (size_t)pmlN * sizeof(double)) != cudaSuccess ||
              cudaMalloc((void**)dev_b[axis],
                         (size_t)pmlN * sizeof(double)) != cudaSuccess ||
              cudaMalloc((void**)dev_invK[axis],
                         (size_t)pmlN * sizeof(double)) != cudaSuccess) {
            fprintf(stderr,
                    "Thread %d: cudaMalloc elastic CPML coeff arrays failed\n",
                    Global.myID);
            MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
          }
        }
        for (int i = 0; i < 3; i++) {
          for (int jx = 0; jx < 3; jx++) {
            if (cudaMalloc(
                    (void**)&(Global.gpuData.pml_elastic_psiArrayDevice[i][jx]),
                    (size_t)pmlN * sizeof(solver_float)) != cudaSuccess) {
              fprintf(stderr,
                      "Thread %d: cudaMalloc pml_elastic_psi[%d][%d] failed\n",
                      Global.myID, i, jx);
              MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
            }
            cudaMemset(Global.gpuData.pml_elastic_psiArrayDevice[i][jx], 0,
                       (size_t)pmlN * sizeof(solver_float));
          }
        }
        /* Stage-2 viscoelastic CPML: relaxed-strain psi memory for the
         * BKT shear/kappa anelastic vectors. Allocated unconditionally
         * (so cudaFree stays symmetric); zeroed; used by
         * kernelPMLElasticStep only when cpmlPmlMaterialDamping==1. */
        for (int i = 0; i < 3; i++) {
          for (int jx = 0; jx < 3; jx++) {
            if (cudaMalloc(
                    (void**)&(Global.gpuData
                                  .pml_elastic_psiShearArrayDevice[i][jx]),
                    (size_t)pmlN * sizeof(solver_float)) != cudaSuccess ||
                cudaMalloc(
                    (void**)&(Global.gpuData
                                  .pml_elastic_psiKappaArrayDevice[i][jx]),
                    (size_t)pmlN * sizeof(solver_float)) != cudaSuccess) {
              fprintf(stderr,
                      "Thread %d: cudaMalloc pml_elastic_psiShear/Kappa"
                      "[%d][%d] failed\n",
                      Global.myID, i, jx);
              MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
            }
            cudaMemset(Global.gpuData.pml_elastic_psiShearArrayDevice[i][jx],
                       0, (size_t)pmlN * sizeof(solver_float));
            cudaMemset(Global.gpuData.pml_elastic_psiKappaArrayDevice[i][jx],
                       0, (size_t)pmlN * sizeof(solver_float));
          }
        }

        if (cudaMemcpy(Global.gpuData.pmlElasticElemMapperDevice,
                       g_cpmlBuild.elastic_elem_list,
                       (size_t)pmlN * sizeof(int32_t),
                       cudaMemcpyHostToDevice) != cudaSuccess) {
          fprintf(stderr,
                  "Thread %d: cudaMemcpy pmlElasticElemMapper failed\n",
                  Global.myID);
          MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
        }
        for (int axis = 0; axis < 3; axis++) {
          if (cudaMemcpy(*dev_a[axis], g_cpmlBuild.elastic_a[axis],
                         (size_t)pmlN * sizeof(double),
                         cudaMemcpyHostToDevice) != cudaSuccess ||
              cudaMemcpy(*dev_b[axis], g_cpmlBuild.elastic_b[axis],
                         (size_t)pmlN * sizeof(double),
                         cudaMemcpyHostToDevice) != cudaSuccess ||
              cudaMemcpy(*dev_invK[axis], g_cpmlBuild.elastic_invKappa[axis],
                         (size_t)pmlN * sizeof(double),
                         cudaMemcpyHostToDevice) != cudaSuccess) {
            fprintf(stderr,
                    "Thread %d: cudaMemcpy elastic CPML coeffs failed\n",
                    Global.myID);
            MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
          }
        }

        /* Elastic material coefficients per PML element. */
        if (cudaMalloc((void**)&(Global.gpuData.pml_elastic_muCoeffArrayDevice),
                       (size_t)pmlN * sizeof(double)) != cudaSuccess ||
            cudaMalloc((void**)&(Global.gpuData.pml_elastic_lambdaCoeffArrayDevice),
                       (size_t)pmlN * sizeof(double)) != cudaSuccess) {
          fprintf(stderr,
                  "Thread %d: cudaMalloc elastic CPML mu/lambda failed\n",
                  Global.myID);
          MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
        }
        if (cudaMemcpy(Global.gpuData.pml_elastic_muCoeffArrayDevice,
                       g_cpmlBuild.elastic_muCoeff,
                       (size_t)pmlN * sizeof(double),
                       cudaMemcpyHostToDevice) != cudaSuccess ||
            cudaMemcpy(Global.gpuData.pml_elastic_lambdaCoeffArrayDevice,
                       g_cpmlBuild.elastic_lambdaCoeff,
                       (size_t)pmlN * sizeof(double),
                       cudaMemcpyHostToDevice) != cudaSuccess) {
          fprintf(stderr,
                  "Thread %d: cudaMemcpy elastic CPML mu/lambda failed\n",
                  Global.myID);
          MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
        }

        /* CFS divergence-side zeta state: 9 per-element arrays
         * (zero-initialised; advanced in kernelPMLElasticStep). */
        for (int i = 0; i < 3; i++) {
          for (int jx = 0; jx < 3; jx++) {
            if (cudaMalloc(
                    (void**)&(Global.gpuData.pml_elastic_divPsiArrayDevice[i][jx]),
                    (size_t)pmlN * sizeof(solver_float)) != cudaSuccess) {
              fprintf(stderr,
                      "Thread %d: cudaMalloc pml_elastic_divPsi[%d][%d] failed\n",
                      Global.myID, i, jx);
              MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
            }
            cudaMemset(Global.gpuData.pml_elastic_divPsiArrayDevice[i][jx], 0,
                       (size_t)pmlN * sizeof(solver_float));
          }
        }

        /* CFS divergence-side coefficients (zetaDecay, zetaGain,
         * memForce) per axis. divPsi above is reused as the zeta state
         * (zero-initialised by the cudaMemset). */
        double** cfs_dev[3] = {
            Global.gpuData.pml_elastic_zetaDecayArrayDevice,
            Global.gpuData.pml_elastic_zetaGainArrayDevice,
            Global.gpuData.pml_elastic_memForceArrayDevice};
        double* cfs_host[3][3] = {
            {g_cpmlBuild.elastic_zetaDecay[0],
             g_cpmlBuild.elastic_zetaDecay[1],
             g_cpmlBuild.elastic_zetaDecay[2]},
            {g_cpmlBuild.elastic_zetaGain[0],
             g_cpmlBuild.elastic_zetaGain[1],
             g_cpmlBuild.elastic_zetaGain[2]},
            {g_cpmlBuild.elastic_memForce[0],
             g_cpmlBuild.elastic_memForce[1],
             g_cpmlBuild.elastic_memForce[2]}};
        for (int k = 0; k < 3; k++) {
          for (int axis = 0; axis < 3; axis++) {
            if (cudaMalloc((void**)&(cfs_dev[k][axis]),
                           (size_t)pmlN * sizeof(double)) != cudaSuccess ||
                cudaMemcpy(cfs_dev[k][axis], cfs_host[k][axis],
                           (size_t)pmlN * sizeof(double),
                           cudaMemcpyHostToDevice) != cudaSuccess) {
              fprintf(stderr,
                      "Thread %d: cudaMalloc/Memcpy elastic CPML CFS "
                      "coeff [%d][%d] failed\n", Global.myID, k, axis);
              MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
            }
          }
        }
        if (cudaMalloc(
                (void**)&(Global.gpuData.pml_elastic_restoreCoeffArrayDevice),
                (size_t)pmlN * sizeof(double)) != cudaSuccess ||
            cudaMemcpy(Global.gpuData.pml_elastic_restoreCoeffArrayDevice,
                       g_cpmlBuild.elastic_restoreCoeff,
                       (size_t)pmlN * sizeof(double),
                       cudaMemcpyHostToDevice) != cudaSuccess) {
          fprintf(stderr,
                  "Thread %d: cudaMalloc/Memcpy elastic CPML "
                  "restoreCoeff failed\n", Global.myID);
          MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
        }
      }
    }

    /* --- Acoustic PML upload --- */
    {
      int32_t pmlN = g_cpmlBuild.acoustic_count;
      Global.gpuData.pmlAcousticElemCount = pmlN;
      Global.mySolver->pmlAcousticElemCount = pmlN;
      if (pmlN > 0) {
        Global.mySolver->pmlAcousticElemList =
            (int32_t*)malloc((size_t)pmlN * sizeof(int32_t));
        if (Global.mySolver->pmlAcousticElemList == NULL) {
          fprintf(stderr,
                  "Thread %d: solver_init: oom pmlAcousticElemList\n",
                  Global.myID);
          MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
        }
        memcpy(Global.mySolver->pmlAcousticElemList,
               g_cpmlBuild.acoustic_elem_list,
               (size_t)pmlN * sizeof(int32_t));

        for (int axis = 0; axis < 3; axis++) {
          Global.mySolver->pml_acoustic_psi[axis] = (solver_float*)calloc(
              (size_t)pmlN, sizeof(solver_float));
          if (Global.mySolver->pml_acoustic_psi[axis] == NULL) {
            fprintf(stderr,
                    "Thread %d: solver_init: oom pml_acoustic_psi[%d]\n",
                    Global.myID, axis);
            MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
          }
        }

        if (cudaMalloc((void**)&(Global.gpuData.pmlAcousticElemMapperDevice),
                       (size_t)pmlN * sizeof(int32_t)) != cudaSuccess) {
          fprintf(stderr,
                  "Thread %d: cudaMalloc pmlAcousticElemMapper failed\n",
                  Global.myID);
          MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
        }
        double** dev_a[3] = {&Global.gpuData.pml_acoustic_a_xArrayDevice,
                             &Global.gpuData.pml_acoustic_a_yArrayDevice,
                             &Global.gpuData.pml_acoustic_a_zArrayDevice};
        double** dev_b[3] = {&Global.gpuData.pml_acoustic_b_xArrayDevice,
                             &Global.gpuData.pml_acoustic_b_yArrayDevice,
                             &Global.gpuData.pml_acoustic_b_zArrayDevice};
        double** dev_invK[3] = {
            &Global.gpuData.pml_acoustic_invKappa_xArrayDevice,
            &Global.gpuData.pml_acoustic_invKappa_yArrayDevice,
            &Global.gpuData.pml_acoustic_invKappa_zArrayDevice};
        for (int axis = 0; axis < 3; axis++) {
          if (cudaMalloc((void**)dev_a[axis],
                         (size_t)pmlN * sizeof(double)) != cudaSuccess ||
              cudaMalloc((void**)dev_b[axis],
                         (size_t)pmlN * sizeof(double)) != cudaSuccess ||
              cudaMalloc((void**)dev_invK[axis],
                         (size_t)pmlN * sizeof(double)) != cudaSuccess) {
            fprintf(stderr,
                    "Thread %d: cudaMalloc acoustic CPML coeff arrays failed\n",
                    Global.myID);
            MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
          }
          if (cudaMalloc(
                  (void**)&(Global.gpuData.pml_acoustic_psiArrayDevice[axis]),
                  (size_t)pmlN * sizeof(solver_float)) != cudaSuccess) {
            fprintf(stderr,
                    "Thread %d: cudaMalloc pml_acoustic_psi[%d] failed\n",
                    Global.myID, axis);
            MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
          }
          cudaMemset(Global.gpuData.pml_acoustic_psiArrayDevice[axis], 0,
                     (size_t)pmlN * sizeof(solver_float));
        }

        if (cudaMemcpy(Global.gpuData.pmlAcousticElemMapperDevice,
                       g_cpmlBuild.acoustic_elem_list,
                       (size_t)pmlN * sizeof(int32_t),
                       cudaMemcpyHostToDevice) != cudaSuccess) {
          fprintf(stderr,
                  "Thread %d: cudaMemcpy pmlAcousticElemMapper failed\n",
                  Global.myID);
          MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
        }
        for (int axis = 0; axis < 3; axis++) {
          if (cudaMemcpy(*dev_a[axis], g_cpmlBuild.acoustic_a[axis],
                         (size_t)pmlN * sizeof(double),
                         cudaMemcpyHostToDevice) != cudaSuccess ||
              cudaMemcpy(*dev_b[axis], g_cpmlBuild.acoustic_b[axis],
                         (size_t)pmlN * sizeof(double),
                         cudaMemcpyHostToDevice) != cudaSuccess ||
              cudaMemcpy(*dev_invK[axis], g_cpmlBuild.acoustic_invKappa[axis],
                         (size_t)pmlN * sizeof(double),
                         cudaMemcpyHostToDevice) != cudaSuccess) {
            fprintf(stderr,
                    "Thread %d: cudaMemcpy acoustic CPML coeffs failed\n",
                    Global.myID);
            MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
          }
        }

        /* Acoustic per-element force coefficient (dt^2 * h^2 / (4 rho)). */
        if (cudaMalloc(
                (void**)&(Global.gpuData.pml_acoustic_forceCoeffArrayDevice),
                (size_t)pmlN * sizeof(double)) != cudaSuccess) {
          fprintf(stderr,
                  "Thread %d: cudaMalloc acoustic CPML forceCoeff failed\n",
                  Global.myID);
          MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
        }
        if (cudaMemcpy(Global.gpuData.pml_acoustic_forceCoeffArrayDevice,
                       g_cpmlBuild.acoustic_forceCoeff,
                       (size_t)pmlN * sizeof(double),
                       cudaMemcpyHostToDevice) != cudaSuccess) {
          fprintf(stderr,
                  "Thread %d: cudaMemcpy acoustic CPML forceCoeff failed\n",
                  Global.myID);
          MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
        }

        /* Step B: divergence-side CFS coeffs (zetaDecay/Gain/memForce
         * per axis + scalar restoreCoeff) and the zeta_p memory state
         * (3 per element, zero-initialised; advanced in
         * kernelPMLAcousticStep). Scalar mirror of the elastic CFS
         * upload. */
        double** cfs_dev[3] = {
            Global.gpuData.pml_acoustic_zetaDecayArrayDevice,
            Global.gpuData.pml_acoustic_zetaGainArrayDevice,
            Global.gpuData.pml_acoustic_memForceArrayDevice};
        double* cfs_host[3][3] = {
            {g_cpmlBuild.acoustic_zetaDecay[0],
             g_cpmlBuild.acoustic_zetaDecay[1],
             g_cpmlBuild.acoustic_zetaDecay[2]},
            {g_cpmlBuild.acoustic_zetaGain[0],
             g_cpmlBuild.acoustic_zetaGain[1],
             g_cpmlBuild.acoustic_zetaGain[2]},
            {g_cpmlBuild.acoustic_memForce[0],
             g_cpmlBuild.acoustic_memForce[1],
             g_cpmlBuild.acoustic_memForce[2]}};
        for (int k = 0; k < 3; k++) {
          for (int axis = 0; axis < 3; axis++) {
            if (cudaMalloc((void**)&(cfs_dev[k][axis]),
                           (size_t)pmlN * sizeof(double)) != cudaSuccess ||
                cudaMemcpy(cfs_dev[k][axis], cfs_host[k][axis],
                           (size_t)pmlN * sizeof(double),
                           cudaMemcpyHostToDevice) != cudaSuccess) {
              fprintf(stderr,
                      "Thread %d: cudaMalloc/Memcpy acoustic CPML CFS "
                      "coeff [%d][%d] failed\n", Global.myID, k, axis);
              MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
            }
          }
          if (cudaMalloc(
                  (void**)&(Global.gpuData.pml_acoustic_zetaPArrayDevice[k]),
                  (size_t)pmlN * sizeof(solver_float)) != cudaSuccess) {
            fprintf(stderr,
                    "Thread %d: cudaMalloc pml_acoustic_zetaP[%d] failed\n",
                    Global.myID, k);
            MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
          }
          cudaMemset(Global.gpuData.pml_acoustic_zetaPArrayDevice[k], 0,
                     (size_t)pmlN * sizeof(solver_float));
        }
        if (cudaMalloc(
                (void**)&(Global.gpuData.pml_acoustic_restoreCoeffArrayDevice),
                (size_t)pmlN * sizeof(double)) != cudaSuccess ||
            cudaMemcpy(Global.gpuData.pml_acoustic_restoreCoeffArrayDevice,
                       g_cpmlBuild.acoustic_restoreCoeff,
                       (size_t)pmlN * sizeof(double),
                       cudaMemcpyHostToDevice) != cudaSuccess) {
          fprintf(stderr,
                  "Thread %d: cudaMalloc/Memcpy acoustic CPML "
                  "restoreCoeff failed\n", Global.myID);
          MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
        }
      }
    }

    /* --- Per-element PML flag (lenum-sized) --- */
    /* Lets the non-PML kernels (kernelStiffnessCalcLocal,
     * kernelAcousticPressureForce, BKT damping) early-exit on PML elements
     * so each element is processed by exactly one stiffness/force kernel.
     * Allocated even if both PML counts are zero so the kernels can
     * unconditionally null-check the pointer. */
    if (cudaMalloc((void**)&(Global.gpuData.pmlFlagArrayDevice),
                   (size_t)Global.myMesh->lenum * sizeof(uint8_t)) !=
        cudaSuccess) {
      fprintf(stderr, "Thread %d: cudaMalloc pmlFlag failed\n", Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
    }
    /* g_cpmlBuild.face_mask already has 1 bit per touched face; cast to a
     * boolean "is any PML face" array for the kernel-side fast check. */
    {
      uint8_t* flag_host =
          (uint8_t*)malloc((size_t)Global.myMesh->lenum * sizeof(uint8_t));
      if (flag_host == NULL) {
        fprintf(stderr, "Thread %d: oom pmlFlag host buffer\n", Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
      }
      for (int32_t ei = 0; ei < Global.myMesh->lenum; ei++) {
        flag_host[ei] = (g_cpmlBuild.face_mask[ei] != 0) ? 1u : 0u;
      }
      if (cudaMemcpy(Global.gpuData.pmlFlagArrayDevice, flag_host,
                     (size_t)Global.myMesh->lenum * sizeof(uint8_t),
                     cudaMemcpyHostToDevice) != cudaSuccess) {
        fprintf(stderr, "Thread %d: cudaMemcpy pmlFlag failed\n",
                Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
      }
      free(flag_host);
    }

    /* --- Optional Dirichlet clamp at outer PML face --- */
    /* When cpml_clamp_outer=yes, mark every node that sits on a domain-
     * extent face touched by a PML (x_min, x_max, y_min, y_max, z_max).
     * z_min is always the free surface so it is never clamped. The clamp
     * is applied per-step in kernelDispCalc by zeroing tm2Disp at flagged
     * nodes.
     *
     * The check is geometric (node coordinate equals nearendp / farendp on
     * the relevant axis) rather than element-based, so it correctly
     * captures every shared node on the outer face regardless of which
     * PML element owns it. */
    if (Param.theCpmlClampOuter == YES) {
      if (cudaMalloc((void**)&(Global.gpuData.pmlClampNodeMaskDevice),
                     (size_t)Global.myMesh->nharbored * sizeof(uint8_t)) !=
          cudaSuccess) {
        fprintf(stderr, "Thread %d: cudaMalloc pmlClampNodeMask failed\n",
                Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
      }
      uint8_t* clamp_host = (uint8_t*)calloc(
          (size_t)Global.myMesh->nharbored, sizeof(uint8_t));
      if (clamp_host == NULL) {
        fprintf(stderr,
                "Thread %d: oom pmlClampNodeMask host buffer\n",
                Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
      }
      const tick_t xn = Global.myOctree->nearendp[0];
      const tick_t yn = Global.myOctree->nearendp[1];
      const tick_t xf = Global.myOctree->farendp[0];
      const tick_t yf = Global.myOctree->farendp[1];
      const tick_t zf = Global.myOctree->farendp[2];
      int64_t local_clamped = 0;
      for (int32_t ni = 0; ni < Global.myMesh->nharbored; ni++) {
        tick_t nx = Global.myMesh->nodeTable[ni].x;
        tick_t ny = Global.myMesh->nodeTable[ni].y;
        tick_t nz = Global.myMesh->nodeTable[ni].z;
        if (nx == xn || nx == xf || ny == yn || ny == yf || nz == zf) {
          clamp_host[ni] = 1u;
          local_clamped++;
        }
      }
      int64_t global_clamped = 0;
      MPI_Reduce(&local_clamped, &global_clamped, 1, MPI_LONG_LONG, MPI_SUM,
                 0, comm_solver);
      if (Global.myID == 0) {
        fprintf(stdout,
                "CPML: cpml_clamp_outer=yes; %lld outer-face nodes clamped "
                "to u=0 globally\n",
                (long long)global_clamped);
      }
      if (cudaMemcpy(Global.gpuData.pmlClampNodeMaskDevice, clamp_host,
                     (size_t)Global.myMesh->nharbored * sizeof(uint8_t),
                     cudaMemcpyHostToDevice) != cudaSuccess) {
        fprintf(stderr,
                "Thread %d: cudaMemcpy pmlClampNodeMask failed\n",
                Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
      }
      free(clamp_host);
    }

    /* --- CPML-aware acoustic-elastic interface coupling correction --- */
    /* For each interface face whose acoustic or elastic side is a PML
     * element, precompute the face-normal-axis ADE coefficients (a, b,
     * 1/kappa) from that PML element so the coupling kernels can apply the
     * same unsplit-CPML modification to the interface flux as the volume
     * kernels apply to the strain/gradient. No-op when there are no
     * interface faces or none of them touch the PML. */
    Global.acousticInterfacePmlFaceCount = 0;
    int64_t pml_face_count = 0;
    if (Global.acousticInterfaceFaces > 0) {
      int64_t nface = Global.acousticInterfaceFaces;
      uint8_t* iface_flag = (uint8_t*)calloc((size_t)nface, sizeof(uint8_t));
      double*  iface_a    = (double*)calloc((size_t)nface, sizeof(double));
      double*  iface_b    = (double*)malloc((size_t)nface * sizeof(double));
      double*  iface_iK   = (double*)malloc((size_t)nface * sizeof(double));
      if (iface_flag == NULL || iface_a == NULL || iface_b == NULL ||
          iface_iK == NULL) {
        fprintf(stderr,
                "Thread %d: oom CPML interface-correction host arrays\n",
                Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR); exit(1);
      }
      /* Default = identity (no-op): b=1, a=0, invKappa=1. */
      for (int64_t i = 0; i < nface; i++) {
        iface_b[i] = 1.0;
        iface_iK[i] = 1.0;
      }
      for (int64_t i = 0; i < nface; i++) {
        int32_t ea = Global.acousticInterfaceElementIds[i];
        int32_t ee = (Global.acousticInterfaceAdjacentElementIndexes != NULL)
                         ? Global.acousticInterfaceAdjacentElementIndexes[i]
                         : -1;
        uint8_t ma = (ea >= 0 && ea < Global.myMesh->lenum)
                         ? g_cpmlBuild.face_mask[ea]
                         : 0;
        uint8_t me = (ee >= 0 && ee < Global.myMesh->lenum)
                         ? g_cpmlBuild.face_mask[ee]
                         : 0;
        if (ma == 0 && me == 0) continue; /* face not in PML: identity */

        /* Use whichever side is the PML element to source the profile.
         * Both sides being PML is rare (interface deep in PML); the
         * acoustic side is preferred since the coupling pressure sample
         * is taken there. */
        int32_t pml_eindex = (ma != 0) ? ea : ee;
        elem_t* pe = &Global.myMesh->elemTable[pml_eindex];
        edata_t* ped = (edata_t*)pe->data;
        int32_t pln0 = pe->lnid[0];
        tick_t ldb[3] = {Global.myMesh->nodeTable[pln0].x,
                         Global.myMesh->nodeTable[pln0].y,
                         Global.myMesh->nodeTable[pln0].z};
        tick_t edgeticks = (tick_t)1 << (PIXELLEVEL - pe->level);
        tick_t ruf[3] = {ldb[0] + edgeticks, ldb[1] + edgeticks,
                         ldb[2] + edgeticks};
        double a3[3], b3[3], iK3[3], d3[3]; /* d3 unused at interface */
        double al3[3], ka3[3]; /* alpha/kappa unused at interface */
        compute_pml_profile_coeffs(
            (ma != 0) ? ma : me, ldb, ruf, Global.myOctree->nearendp,
            Global.myOctree->farendp, Global.myOctree->ticksize,
            g_cpmlBuild.thickness_m, (double)ped->Vp, Param.theDeltaT,
            Param.theCpmlPolyOrder, Param.theCpmlReflectionTarget,
            Param.theCpmlKappaMax, g_cpmlBuild.alpha_max,
            Param.theCpmlMpmlCrossFraction,
            Param.theCpmlAlphaFloorFraction, a3, b3, iK3, d3, al3, ka3);
        (void)d3;
        (void)al3;
        (void)ka3;

        /* Pick the axis most aligned with the interface normal. */
        double nx = fabs(Global.acousticInterfaceFaceNormals[3 * i + 0]);
        double ny = fabs(Global.acousticInterfaceFaceNormals[3 * i + 1]);
        double nz = fabs(Global.acousticInterfaceFaceNormals[3 * i + 2]);
        int axis = (nx >= ny && nx >= nz) ? 0 : ((ny >= nz) ? 1 : 2);

        iface_flag[i] = 1u;
        iface_a[i]    = a3[axis];
        iface_b[i]    = b3[axis];
        iface_iK[i]   = iK3[axis];
        pml_face_count++;
      }

      if (pml_face_count > 0) {
        Global.acousticInterfacePmlFlag = iface_flag;
        Global.acousticInterfacePml_a = iface_a;
        Global.acousticInterfacePml_b = iface_b;
        Global.acousticInterfacePml_invKappa = iface_iK;
        Global.acousticInterfacePmlFaceCount = pml_face_count;
      } else {
        /* No interface face is in the PML: keep everything NULL so the
         * coupling kernels skip the correction entirely. */
        free(iface_flag);
        free(iface_a);
        free(iface_b);
        free(iface_iK);
      }
    }

    /* All solver ranks must enter this collective. Many 3D decompositions have
     * ranks with no local water-solid interface faces, especially for a flat
     * water layer, and letting those ranks skip the reduction deadlocks the
     * interface-owning ranks during Phase B startup. */
    int64_t global_pml_iface = 0;
    MPI_Reduce(&pml_face_count, &global_pml_iface, 1, MPI_LONG_LONG,
               MPI_SUM, 0, comm_solver);
    if (Global.myID == 0 && global_pml_iface > 0) {
      fprintf(stdout,
              "CPML: %lld acoustic-elastic interface faces lie inside the "
              "PML zone; applying unsplit-CPML interface correction\n",
              (long long)global_pml_iface);
    }

    /* Free the build-time host scratch. The compact per-type elem_lists
     * survive in mySolver; coefficients live only on the GPU from here on. */
    free(g_cpmlBuild.face_mask);
    g_cpmlBuild.face_mask = NULL;
    free(g_cpmlBuild.elastic_elem_list);
    g_cpmlBuild.elastic_elem_list = NULL;
    free(g_cpmlBuild.acoustic_elem_list);
    g_cpmlBuild.acoustic_elem_list = NULL;
    free(g_cpmlBuild.elastic_muCoeff);
    g_cpmlBuild.elastic_muCoeff = NULL;
    free(g_cpmlBuild.elastic_lambdaCoeff);
    g_cpmlBuild.elastic_lambdaCoeff = NULL;
    free(g_cpmlBuild.elastic_restoreCoeff);
    g_cpmlBuild.elastic_restoreCoeff = NULL;
    free(g_cpmlBuild.acoustic_forceCoeff);
    g_cpmlBuild.acoustic_forceCoeff = NULL;
    free(g_cpmlBuild.acoustic_restoreCoeff);
    g_cpmlBuild.acoustic_restoreCoeff = NULL;
    for (int axis = 0; axis < 3; axis++) {
      free(g_cpmlBuild.elastic_a[axis]);
      g_cpmlBuild.elastic_a[axis] = NULL;
      free(g_cpmlBuild.elastic_b[axis]);
      g_cpmlBuild.elastic_b[axis] = NULL;
      free(g_cpmlBuild.elastic_invKappa[axis]);
      g_cpmlBuild.elastic_invKappa[axis] = NULL;
      free(g_cpmlBuild.elastic_zetaDecay[axis]);
      g_cpmlBuild.elastic_zetaDecay[axis] = NULL;
      free(g_cpmlBuild.elastic_zetaGain[axis]);
      g_cpmlBuild.elastic_zetaGain[axis] = NULL;
      free(g_cpmlBuild.elastic_memForce[axis]);
      g_cpmlBuild.elastic_memForce[axis] = NULL;
      free(g_cpmlBuild.acoustic_a[axis]);
      g_cpmlBuild.acoustic_a[axis] = NULL;
      free(g_cpmlBuild.acoustic_b[axis]);
      g_cpmlBuild.acoustic_b[axis] = NULL;
      free(g_cpmlBuild.acoustic_invKappa[axis]);
      g_cpmlBuild.acoustic_invKappa[axis] = NULL;
      free(g_cpmlBuild.acoustic_zetaDecay[axis]);
      g_cpmlBuild.acoustic_zetaDecay[axis] = NULL;
      free(g_cpmlBuild.acoustic_zetaGain[axis]);
      g_cpmlBuild.acoustic_zetaGain[axis] = NULL;
      free(g_cpmlBuild.acoustic_memForce[axis]);
      g_cpmlBuild.acoustic_memForce[axis] = NULL;
    }
  }
#endif /* BOUNDARY */

  /* Build the communication schedules */
  schedule_build(Global.myMesh, Global.mySolver->dn_sched,
                 Global.mySolver->an_sched);

#ifdef DEBUG
  /* For debug purpose, add gnid into the data field. */
  c_outsize = sizeof(n_t) + sizeof(int64_t);
  c_insize = sizeof(n_t) + sizeof(int64_t);
  s_outsize = sizeof(n_t) + sizeof(int64_t);
  s_insize = sizeof(n_t) + sizeof(int64_t);
#else
  c_outsize = sizeof(n_t);
  c_insize = sizeof(n_t);
  s_outsize = sizeof(n_t);
  s_insize = sizeof(n_t);
#endif /* DEBUG */

  schedule_prepare(Global.mySolver->dn_sched, c_outsize, c_insize, s_outsize,
                   s_insize);

  schedule_prepare(Global.mySolver->an_sched, c_outsize, c_insize, s_outsize,
                   s_insize);

  /* Send mass information of dangling nodes to their owners */
  schedule_senddata(Global.mySolver->dn_sched, Global.mySolver->nTable,
                    sizeof(n_t) / sizeof(solver_float), CONTRIBUTION,
                    DN_MASS_MSG);

  /* Distribute the mass from dangling nodes to anchored nodes. (local) */
  compute_adjust(Global.mySolver->nTable, sizeof(n_t) / sizeof(solver_float),
                 DISTRIBUTION);

  /* Send mass information of anchored nodes to their owner processors*/
  schedule_senddata(Global.mySolver->an_sched, Global.mySolver->nTable,
                    sizeof(n_t) / sizeof(solver_float), CONTRIBUTION,
                    AN_MASS_MSG);

  if (acoustic_pressure_enabled()) {
    int64_t localPressureNodes = 0;
    int64_t globalPressureNodes = 0;
    int64_t localPressureFreeSurfaceNodes = 0;
    int64_t globalPressureFreeSurfaceNodes = 0;

    schedule_senddata(Global.mySolver->dn_sched, Global.mySolver->pressure_mass,
                      1, CONTRIBUTION, DN_PRESSURE_MASS_MSG);
    compute_adjust(Global.mySolver->pressure_mass, 1, DISTRIBUTION);
    schedule_senddata(Global.mySolver->an_sched, Global.mySolver->pressure_mass,
                      1, CONTRIBUTION, AN_PRESSURE_MASS_MSG);
    /* Step B: assemble the acoustic CPML nodal damping across ranks
     * exactly like pressure_mass -- it is a node-accumulated quantity,
     * so partition-boundary nodes must get the full cross-rank sum,
     * otherwise damping is under-counted (asymmetric) there. All ranks
     * call these collectives unconditionally (pressure_damp is allocated
     * + zeroed for every rank; it is simply all-zero when there is no
     * acoustic PML, so this is a strict no-op then). */
    schedule_senddata(Global.mySolver->dn_sched, Global.mySolver->pressure_damp,
                      1, CONTRIBUTION, DN_PRESSURE_DAMP_MSG);
    compute_adjust(Global.mySolver->pressure_damp, 1, DISTRIBUTION);
    schedule_senddata(Global.mySolver->an_sched, Global.mySolver->pressure_damp,
                      1, CONTRIBUTION, AN_PRESSURE_DAMP_MSG);
    schedule_senddata(Global.mySolver->dn_sched, pressureFreeSurfaceHost, 1,
                      CONTRIBUTION, DN_PRESSURE_SURFACE_MSG);
    compute_adjust(pressureFreeSurfaceHost, 1, DISTRIBUTION);
    schedule_senddata(Global.mySolver->an_sched, pressureFreeSurfaceHost, 1,
                      CONTRIBUTION, AN_PRESSURE_SURFACE_MSG);

    for (int32_t nindex = 0; nindex < Global.myMesh->nharbored; nindex++) {
      if (Global.mySolver->pressure_mass[nindex] > 0.0) {
        localPressureNodes++;
      }
      if (pressureFreeSurfaceHost[nindex] > 0.0) {
        localPressureFreeSurfaceNodes++;
      }
    }

    MPI_Allreduce(&localPressureElements, &globalPressureElements, 1,
                  MPI_INT64_T, MPI_SUM, comm_solver);
    MPI_Allreduce(&localPressureNodes, &globalPressureNodes, 1, MPI_INT64_T,
                  MPI_SUM, comm_solver);
    MPI_Allreduce(&localPressureFreeSurfaceNodes,
                  &globalPressureFreeSurfaceNodes, 1, MPI_INT64_T, MPI_SUM,
                  comm_solver);
    g_acousticPressureElementCount = globalPressureElements;
    g_acousticPressureNodeCount = globalPressureNodes;
    g_acousticPressureFreeSurfaceNodeCount = globalPressureFreeSurfaceNodes;
    if (Global.myID == 0) {
      monitor_print(
          "Acoustic pressure formulation: pressure elements %.1fM, pressure nodes %.1fM, free-surface nodes %.1fM\n",
          (double)globalPressureElements / 1.0e6,
          (double)globalPressureNodes / 1.0e6,
          (double)globalPressureFreeSurfaceNodes / 1.0e6);
    }
  } else {
    g_acousticPressureElementCount = 0;
    g_acousticPressureNodeCount = 0;
    g_acousticPressureFreeSurfaceNodeCount = 0;
  }

  if (Global.myID == 0) {
    fprintf(stdout, "GPU device memory allocation: start\n");
  }

  /* Save device configuration */
  Global.mySolver->gpu_spec = &(Global.gpu_spec);

  /* Allocate memory for GPU */
  if (cudaMalloc((void**)&(Global.gpuData.tm1Device),
                 Global.myMesh->nharbored * sizeof(fvector_t)) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to allocate tm1 device memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMalloc((void**)&(Global.gpuData.tm2Device),
                 Global.myMesh->nharbored * sizeof(fvector_t)) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to allocate tm2 device memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if ((Param.printStationAccelerations == YES) || acoustic_pressure_enabled()) {
    if (cudaMalloc((void**)&(Global.gpuData.tm3Device),
                   Global.myMesh->nharbored * sizeof(fvector_t)) !=
        cudaSuccess) {
      fprintf(stderr, "Thread %d: Failed to allocate tm3 device memory\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
  }

  if (cudaMalloc((void**)&(Global.gpuData.forceDevice),
                 Global.myMesh->nharbored * sizeof(fvector_t)) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to allocate force memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (acoustic_pressure_enabled()) {
    if (cudaMalloc((void**)&(Global.gpuData.pressureTm1Device),
                   Global.myMesh->nharbored * sizeof(solver_float)) !=
        cudaSuccess) {
      fprintf(stderr, "Thread %d: Failed to allocate pressure tm1 memory\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    if (cudaMalloc((void**)&(Global.gpuData.pressureTm2Device),
                   Global.myMesh->nharbored * sizeof(solver_float)) !=
        cudaSuccess) {
      fprintf(stderr, "Thread %d: Failed to allocate pressure tm2 memory\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    if (cudaMalloc((void**)&(Global.gpuData.pressureForceDevice),
                   Global.myMesh->nharbored * sizeof(solver_float)) !=
        cudaSuccess) {
      fprintf(stderr, "Thread %d: Failed to allocate pressure force memory\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    if (cudaMalloc((void**)&(Global.gpuData.pressureMassDevice),
                   Global.myMesh->nharbored * sizeof(solver_float)) !=
        cudaSuccess) {
      fprintf(stderr, "Thread %d: Failed to allocate pressure mass memory\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    if (cudaMalloc((void**)&(Global.gpuData.pressureDampDevice),
                   Global.myMesh->nharbored * sizeof(solver_float)) !=
        cudaSuccess) {
      fprintf(stderr, "Thread %d: Failed to allocate pressure damp memory\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    if (cudaMalloc(
            (void**)&Global.gpuData.pressureFreeSurfaceNodeFlagsDevice,
            Global.myMesh->nharbored * sizeof(uint8_t)) != cudaSuccess) {
      fprintf(stderr,
              "Thread %d: Failed to allocate pressure free-surface flag memory\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
  }
  if (cudaMalloc((void**)&(Global.gpuData.conv_shear_1Device),
                 8 * Global.myMesh->lenum * sizeof(fvector_t)) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to allocate conv_shear_1 memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMalloc((void**)&(Global.gpuData.conv_shear_2Device),
                 8 * Global.myMesh->lenum * sizeof(fvector_t)) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to allocate conv_shear_2 memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
    if (cudaMalloc((void**)&(Global.gpuData.conv_shear_3Device),
       8 * Global.myMesh->lenum * sizeof(fvector_t)) != cudaSuccess) {
      fprintf(stderr, "Thread %d: Failed to allocate conv_shear_3 memory\n",
        Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
  if (cudaMalloc((void**)&(Global.gpuData.conv_kappa_1Device),
                 8 * Global.myMesh->lenum * sizeof(fvector_t)) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to allocate conv_kappa_1 memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMalloc((void**)&(Global.gpuData.conv_kappa_2Device),
                 8 * Global.myMesh->lenum * sizeof(fvector_t)) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to allocate conv_kappa_2 memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
    if (cudaMalloc((void**)&(Global.gpuData.conv_kappa_3Device),
       8 * Global.myMesh->lenum * sizeof(fvector_t)) != cudaSuccess) {
      fprintf(stderr, "Thread %d: Failed to allocate conv_kappa_3 memory\n",
        Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

  if (cudaMalloc((void**)&(Global.gpuData.c1ArrayDevice),
                 Global.myMesh->lenum * sizeof(double)) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to allocate c1 array memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMalloc((void**)&(Global.gpuData.c2ArrayDevice),
                 Global.myMesh->lenum * sizeof(double)) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to allocate c1 array memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMalloc((void**)&(Global.gpuData.c3ArrayDevice),
                 Global.myMesh->lenum * sizeof(double)) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to allocate c3 array memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  if (cudaMalloc((void**)&(Global.gpuData.lnidArrayDevice),
                 Global.myMesh->lenum * 8 * sizeof(int32_t)) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to allocate lnid array memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  if (cudaMalloc((void**)&(Global.gpuData.acousticFaceFlagsArrayDevice),
                 Global.myMesh->lenum * sizeof(uint8_t)) != cudaSuccess) {
    fprintf(stderr,
            "Thread %d: Failed to allocate acoustic interface flag memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  Global.gpuData.acousticInterfaceFaceCount = 0;
  /* Independent velocity-potential coupling + interface source filter. */
  Global.gpuData.acousticVelocityPotential =
      Param.theAcousticVelocityPotential;
  Global.gpuData.acousticInterfaceFilterEnabled =
      (Param.theAcousticInterfaceFilter == YES) ? 1 : 0;
  Global.gpuData.acousticInterfaceFilterAlpha =
      Param.theAcousticInterfaceFilterAlpha;
  Global.gpuData.acousticInterfaceFilterStrength =
      Param.theAcousticInterfaceFilterStrength;
  /* Stiffness-proportional (Kelvin-Voigt) acoustic damping coefficient,
   * derived self-containedly from the target quality factor Q -- it must
   * NOT depend on theBBase: compute_setab() zeroes theBBase for every
   * non-RAYLEIGH/MASS run (NONE and all BKT*, psolve.cu:15312), which is
   * exactly the no-attenuation lossless benchmark this knob targets, so
   * a theBBase-scaled beta silently no-ops there (observed: Q=500 v22 ==
   * v20 within GPU atomic float noise).
   *
   * The acoustic force adds  -coeff*beta*K*(p^n - p^{n-1}). With the
   * leapfrog M(p^{n+1}-2p^n+p^{n-1})/dt^2 = F and modal cK = w^2 M, the
   * added term gives modal EOM  p'' + (beta*w^2*dt) p' + w^2 p = 0, i.e.
   * damping ratio  zeta(w) = beta*w*dt/2  and  Q(w) = 1/(beta*w*dt).
   * Pinning Q at the simulation reference angular frequency
   * w_ref = 2*pi*Param.theFreq (the same freq compute_setab uses, always
   * set regardless of damping type) gives
   *     beta = 1 / (Q * 2*pi*theFreq * theDeltaT).
   * zeta ~ w, so the high-frequency lossless-cavity ring is damped far
   * more strongly than the sub-reference-band physical signal. Q unset
   * / <= 0 (or degenerate freq/dt) => beta = 0 => acoustic stays exactly
   * lossless and bit-identical to the validated baseline. */
  {
    double q_ac = Param.theAcousticQualityFactor;
    double w_ref = 2.0 * PI * Param.theFreq;
    Global.gpuData.acousticStiffnessDampingBeta =
        (q_ac > 0.0 && Param.theFreq > 0.0 && Param.theDeltaT > 0.0)
            ? (1.0 / (q_ac * w_ref * Param.theDeltaT))
            : 0.0;
    if (Global.myID == 0) {
      double bA = Global.gpuData.acousticStiffnessDampingBeta;
      if (bA != 0.0) {
        monitor_print(
            "Acoustic stiffness damping:         Q=%.3f -> beta=%.6e "
            "(zeta@%.3fHz=%.3e; ACTIVE)\n",
            q_ac, bA, Param.theFreq,
            0.5 * bA * w_ref * Param.theDeltaT);
      } else {
        monitor_print("Acoustic stiffness damping:         beta=0 "
                      "(lossless; baseline bit-identical)\n");
      }
    }
  }
  Global.gpuData.acousticInterfaceFilterStateDevice = NULL;
  Global.gpuData.acousticInterfaceElementIdsDevice = NULL;
  Global.gpuData.acousticInterfaceFaceIdsDevice = NULL;
  Global.gpuData.acousticInterfaceFaceAreasDevice = NULL;
  Global.gpuData.acousticInterfaceFaceCentersDevice = NULL;
  Global.gpuData.acousticInterfaceFaceNormalsDevice = NULL;
  Global.gpuData.acousticInterfaceAdjacentElementIndexesDevice = NULL;
  Global.gpuData.acousticInterfaceAdjacentCentersDevice = NULL;
  if (acoustic_pressure_enabled()) {
    if (cudaMalloc((void**)&(Global.gpuData.acousticElementFlagsArrayDevice),
                   Global.myMesh->lenum * sizeof(uint8_t)) != cudaSuccess) {
      fprintf(stderr,
              "Thread %d: Failed to allocate acoustic pressure element flag memory\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    if (cudaMalloc((void**)&(Global.gpuData.acousticPressureCoeffArrayDevice),
                   Global.myMesh->lenum * sizeof(double)) != cudaSuccess) {
      fprintf(stderr,
              "Thread %d: Failed to allocate acoustic pressure coefficient memory\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    if (cudaMalloc(
            (void**)&(Global.gpuData.acousticPressureBoundaryCoeffArrayDevice),
            Global.myMesh->lenum * sizeof(double)) != cudaSuccess) {
      fprintf(stderr,
              "Thread %d: Failed to allocate acoustic pressure boundary coefficient memory\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    Global.gpuData.acousticInterfaceFaceCount =
        (int32_t)Global.acousticInterfaceFaces;
    if (Global.gpuData.acousticInterfaceFaceCount > 0) {
      size_t face_count =
          (size_t)Global.gpuData.acousticInterfaceFaceCount;

      /* Per-face EMA state for the interface source filter (gated; only
       * read/written when acousticInterfaceFilterEnabled). */
      if (cudaMalloc(
              (void**)&(Global.gpuData.acousticInterfaceFilterStateDevice),
              face_count * sizeof(solver_float)) != cudaSuccess) {
        fprintf(stderr,
                "Thread %d: Failed to allocate acoustic interface filter "
                "state memory\n", Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }
      cudaMemset(Global.gpuData.acousticInterfaceFilterStateDevice, 0,
                 face_count * sizeof(solver_float));

      if (cudaMalloc(
              (void**)&(Global.gpuData.acousticInterfaceElementIdsDevice),
              face_count * sizeof(int32_t)) != cudaSuccess) {
        fprintf(stderr,
                "Thread %d: Failed to allocate acoustic interface element table\n",
                Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }
      if (cudaMalloc((void**)&(Global.gpuData.acousticInterfaceFaceIdsDevice),
                     face_count * sizeof(uint8_t)) != cudaSuccess) {
        fprintf(stderr,
                "Thread %d: Failed to allocate acoustic interface face table\n",
                Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }
      if (cudaMalloc((void**)&(Global.gpuData.acousticInterfaceFaceAreasDevice),
                     face_count * sizeof(double)) != cudaSuccess) {
        fprintf(stderr,
                "Thread %d: Failed to allocate acoustic interface face areas\n",
                Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }
      if (cudaMalloc(
              (void**)&(Global.gpuData.acousticInterfaceFaceCentersDevice),
              3 * face_count * sizeof(double)) != cudaSuccess) {
        fprintf(stderr,
                "Thread %d: Failed to allocate acoustic interface fragment centers\n",
                Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }
      if (cudaMalloc(
              (void**)&(Global.gpuData.acousticInterfaceFaceNormalsDevice),
              3 * face_count * sizeof(double)) != cudaSuccess) {
        fprintf(stderr,
                "Thread %d: Failed to allocate acoustic interface face normals\n",
                Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }
      if (cudaMalloc(
              (void**)&(
                  Global.gpuData.acousticInterfaceAdjacentElementIndexesDevice),
              face_count * sizeof(int32_t)) != cudaSuccess) {
        fprintf(stderr,
                "Thread %d: Failed to allocate acoustic interface adjacent element indexes\n",
                Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }
      if (cudaMalloc(
              (void**)&(Global.gpuData.acousticInterfaceAdjacentCentersDevice),
              3 * face_count * sizeof(double)) != cudaSuccess) {
        fprintf(stderr,
                "Thread %d: Failed to allocate acoustic interface adjacent centers\n",
                Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }
    }
  }

  if (cudaMalloc((void**)&(Global.gpuData.g0_shearArrayDevice),
                 Global.myMesh->lenum * sizeof(double)) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to allocate g0 shear array memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMalloc((void**)&(Global.gpuData.g1_shearArrayDevice),
                 Global.myMesh->lenum * sizeof(double)) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to allocate g1 shear array memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
    if (cudaMalloc((void**)&(Global.gpuData.g2_shearArrayDevice),
       Global.myMesh->lenum * sizeof(double)) != cudaSuccess) {
      fprintf(stderr, "Thread %d: Failed to allocate g2 shear array memory\n",
        Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

  if (cudaMalloc((void**)&(Global.gpuData.g0_kappaArrayDevice),
                 Global.myMesh->lenum * sizeof(double)) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to allocate g0 kappa array memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMalloc((void**)&(Global.gpuData.g1_kappaArrayDevice),
                 Global.myMesh->lenum * sizeof(double)) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to allocate g1 kappa array memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
    if (cudaMalloc((void**)&(Global.gpuData.g2_kappaArrayDevice),
       Global.myMesh->lenum * sizeof(double)) != cudaSuccess) {
      fprintf(stderr, "Thread %d: Failed to allocate g2 kappa array memory\n",
        Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

  if (cudaMalloc((void**)&(Global.gpuData.b_shearArrayDevice),
                 Global.myMesh->lenum * sizeof(double)) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to allocate b shear array memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMalloc((void**)&(Global.gpuData.a0_shearArrayDevice),
                 Global.myMesh->lenum * sizeof(double)) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to allocate a0 shear array memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMalloc((void**)&(Global.gpuData.a1_shearArrayDevice),
                 Global.myMesh->lenum * sizeof(double)) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to allocate a1 shear array memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
    if (cudaMalloc((void**)&(Global.gpuData.a2_shearArrayDevice),
       Global.myMesh->lenum * sizeof(double)) != cudaSuccess) {
      fprintf(stderr, "Thread %d: Failed to allocate a2 shear array memory\n",
        Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
  if (cudaMalloc((void**)&(Global.gpuData.b_kappaArrayDevice),
                 Global.myMesh->lenum * sizeof(double)) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to allocate b kappa array memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMalloc((void**)&(Global.gpuData.a0_kappaArrayDevice),
                 Global.myMesh->lenum * sizeof(double)) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to allocate a0 kappa array memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMalloc((void**)&(Global.gpuData.a1_kappaArrayDevice),
                 Global.myMesh->lenum * sizeof(double)) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to allocate a1 kappa array memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
    if (cudaMalloc((void**)&(Global.gpuData.a2_kappaArrayDevice),
       Global.myMesh->lenum * sizeof(double)) != cudaSuccess) {
      fprintf(stderr, "Thread %d: Failed to allocate a2 kappa array memory\n",
        Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

  if (cudaMalloc((void**)&(Global.gpuData.shearVectorDevice),
                 8 * Global.myMesh->lenum * sizeof(fvector_t)) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to allocate shear vector memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMalloc((void**)&(Global.gpuData.kappaVectorDevice),
                 8 * Global.myMesh->lenum * sizeof(fvector_t)) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to allocate kappa vector memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  if (cudaMalloc((void**)&(Global.gpuData.mass_simpleArrayDevice),
                 Global.myMesh->nharbored * sizeof(solver_float)) !=
      cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to allocate mass_simple array memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

    for (int i = 0; i < 3; i++) {
    if (cudaMalloc((void**)&(Global.gpuData.mass2_minusaMArrayDevice[i]),
                   Global.myMesh->nharbored * sizeof(solver_float)) !=
        cudaSuccess) {
      fprintf(stderr,
              "Thread %d: Failed to allocate mass2_minusaM array memory\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    if (cudaMalloc((void**)&(Global.gpuData.mass_minusaMArrayDevice[i]),
                   Global.myMesh->nharbored * sizeof(solver_float)) !=
        cudaSuccess) {
      fprintf(stderr,
              "Thread %d: Failed to allocate mass_minusaM array memory\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
  }

  Global.gpuData.nharbored = Global.myMesh->nharbored;
  Global.gpuData.lenum = Global.myMesh->lenum;

  /* Temporary structure of arrays */
  double* c1Array = (double*)malloc(Global.myMesh->lenum * sizeof(double));
  double* c2Array = (double*)malloc(Global.myMesh->lenum * sizeof(double));
  double* c3Array = (double*)malloc(Global.myMesh->lenum * sizeof(double));
  int32_t* lnidArray =
      (int32_t*)malloc(Global.myMesh->lenum * 8 * sizeof(int32_t));
  uint8_t* acousticFaceFlagsArray =
      (uint8_t*)calloc((size_t)Global.myMesh->lenum, sizeof(uint8_t));
  if (acousticFaceFlagsArray == NULL) {
    fprintf(stderr,
            "Thread %d: Failed to allocate acoustic interface flag host array\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  uint8_t* acousticElementFlagsArray = NULL;
  uint8_t* pressureFreeSurfaceNodeFlagsArray = NULL;
  double* acousticPressureCoeffArray = NULL;
  double* acousticPressureBoundaryCoeffArray = NULL;
  if (acoustic_pressure_enabled()) {
    acousticElementFlagsArray =
        (uint8_t*)calloc((size_t)Global.myMesh->lenum, sizeof(uint8_t));
    pressureFreeSurfaceNodeFlagsArray =
        (uint8_t*)calloc((size_t)Global.myMesh->nharbored, sizeof(uint8_t));
    acousticPressureCoeffArray =
        (double*)calloc((size_t)Global.myMesh->lenum, sizeof(double));
    acousticPressureBoundaryCoeffArray =
        (double*)calloc((size_t)Global.myMesh->lenum, sizeof(double));
    if ((acousticElementFlagsArray == NULL) ||
        (pressureFreeSurfaceNodeFlagsArray == NULL) ||
        (acousticPressureCoeffArray == NULL) ||
        (acousticPressureBoundaryCoeffArray == NULL)) {
      fprintf(stderr,
              "Thread %d: Failed to allocate acoustic pressure host arrays\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
  }

  double* g0_shearArray =
      (double*)malloc(Global.myMesh->lenum * sizeof(double));
  double* g1_shearArray =
      (double*)malloc(Global.myMesh->lenum * sizeof(double));
    double* g2_shearArray =
      (double*)malloc(Global.myMesh->lenum * sizeof(double));
  double* g0_kappaArray =
      (double*)malloc(Global.myMesh->lenum * sizeof(double));
  double* g1_kappaArray =
      (double*)malloc(Global.myMesh->lenum * sizeof(double));
    double* g2_kappaArray =
      (double*)malloc(Global.myMesh->lenum * sizeof(double));

  double* b_shearArray = (double*)malloc(Global.myMesh->lenum * sizeof(double));
  double* a0_shearArray =
      (double*)malloc(Global.myMesh->lenum * sizeof(double));
  double* a1_shearArray =
      (double*)malloc(Global.myMesh->lenum * sizeof(double));
    double* a2_shearArray =
      (double*)malloc(Global.myMesh->lenum * sizeof(double));
  double* b_kappaArray = (double*)malloc(Global.myMesh->lenum * sizeof(double));
  double* a0_kappaArray =
      (double*)malloc(Global.myMesh->lenum * sizeof(double));
  double* a1_kappaArray =
      (double*)malloc(Global.myMesh->lenum * sizeof(double));
    double* a2_kappaArray =
      (double*)malloc(Global.myMesh->lenum * sizeof(double));

  solver_float* mass_simpleArray =
      (solver_float*)malloc(Global.myMesh->nharbored * sizeof(solver_float));
  solver_float *mass2_minusaMArray[3], *mass_minusaMArray[3];
  for (int i = 0; i < 3; i++) {
    mass2_minusaMArray[i] =
        (solver_float*)malloc(Global.myMesh->nharbored * sizeof(solver_float));
    mass_minusaMArray[i] =
        (solver_float*)malloc(Global.myMesh->nharbored * sizeof(solver_float));
  }

  double rmax;
  if (Param.theTypeOfDamping == BKT) {
    rmax = 2. * M_PI * Param.theFreq * Param.theDeltaT;
  } else {
    rmax = Param.theDeltaT;
  }

  for (int i = 0; i < Global.myMesh->lenum; i++) {
    // if (cudaMemcpy(&Global.gpuData.matPropsDevice[i],
    //	     Global.myMesh->elemTable[i].data,
    //  		     sizeof(edata_t), cudaMemcpyHostToDevice) !=
    //  cudaSuccess) { 	fprintf(stderr,
    //	"Thread %d: Failed to copy mat props for element %d\n",
    //		Global.myID, i);
    //	MPI_Abort(MPI_COMM_WORLD, ERROR);
    //	exit(1);
    // }

    c1Array[i] = Global.mySolver->eTable[i].c1;
    c2Array[i] = Global.mySolver->eTable[i].c2;
    c3Array[i] = Global.mySolver->eTable[i].c3;
    if (Global.acousticFaceFlags != NULL) {
      acousticFaceFlagsArray[i] = Global.acousticFaceFlags[i];
    }

    for (int j = 0; j < 8; j++) {
      lnidArray[j * Global.myMesh->lenum + i] =
          Global.myMesh->elemTable[i].lnid[j];
    }

    /* loop variables */
    elem_t* elemp;
    edata_t* edata;

    /* capturing the elements */
    elemp = &Global.myMesh->elemTable[i];
    edata = (edata_t*)elemp->data;
    acoustic_sanitize_element(edata);

    if (acoustic_pressure_enabled() && element_is_acoustic(edata)) {
      acousticElementFlagsArray[i] = 1;
      /* Interior acoustic stiffness coefficient for the consistent
       * (fully integrated) trilinear-hex Laplacian assembled in
       * kernelAcousticPressureForce as -coeff * sum_j shape(i,j)*p[j],
       * where acoustic_scalar_stiffness_shape(i,j) = 8 * \int_ref
       * grad(N_i).grad(N_j).  The physical consistent stiffness for a
       * cube of edge h with density rho is
       *   K_ij = (1/rho) \int_phys grad(N_i).grad(N_j) dV
       *        = (1/rho) * (h/2) * (shape/8) = shape * h / (16 rho),
       * cross-checked by the tensor product Kx (x) My (x) Mz which gives
       * the exact K_00 = h/(3 rho) = [h/(16 rho)] * shape(0,0)=16/3.
       * The lumped (1/kappa) mass volume/(8*kappa) is exact, so the
       * stiffness prefactor alone sets the discrete wave speed.  The
       * original /(32 rho) was a factor of 2 too small, giving a
       * discrete acoustic speed of Vp/sqrt(2) (~0.707 Vp); /(16 rho)
       * restores the correct leading-order speed Vp. */
      acousticPressureCoeffArray[i] =
          Param.theDeltaTSquared * edata->edgesize / (16.0 * edata->rho);
      acousticPressureBoundaryCoeffArray[i] =
          0.25 * edata->edgesize * edata->edgesize;
    }

    g0_shearArray[i] = edata->g0_shear * rmax;
    g1_shearArray[i] = edata->g1_shear * rmax;
    g2_shearArray[i] = edata->g2_shear * rmax;
    g0_kappaArray[i] = edata->g0_kappa * rmax;
    g1_kappaArray[i] = edata->g1_kappa * rmax;
    g2_kappaArray[i] = edata->g2_kappa * rmax;

    b_shearArray[i] = edata->b_shear;
    a0_shearArray[i] = edata->a0_shear;
    a1_shearArray[i] = edata->a1_shear;
    a2_shearArray[i] = edata->a2_shear;
    b_kappaArray[i] = edata->b_kappa;
    a0_kappaArray[i] = edata->a0_kappa;
    a1_kappaArray[i] = edata->a1_kappa;
    a2_kappaArray[i] = edata->a2_kappa;
  }

  for (int i = 0; i < Global.myMesh->nharbored; i++) {
    mass_simpleArray[i] = Global.mySolver->nTable[i].mass_simple;
    mass2_minusaMArray[0][i] = Global.mySolver->nTable[i].mass2_minusaM[0];
    mass2_minusaMArray[1][i] = Global.mySolver->nTable[i].mass2_minusaM[1];
    mass2_minusaMArray[2][i] = Global.mySolver->nTable[i].mass2_minusaM[2];
    mass_minusaMArray[0][i] = Global.mySolver->nTable[i].mass_minusaM[0];
    mass_minusaMArray[1][i] = Global.mySolver->nTable[i].mass_minusaM[1];
    mass_minusaMArray[2][i] = Global.mySolver->nTable[i].mass_minusaM[2];
    if (acoustic_pressure_enabled()) {
      pressureFreeSurfaceNodeFlagsArray[i] =
          (pressureFreeSurfaceHost[i] > 0.0) ? 1u : 0u;
    }
  }

  if (cudaMemcpy(Global.gpuData.c1ArrayDevice, c1Array,
                 Global.myMesh->lenum * sizeof(double),
                 cudaMemcpyHostToDevice) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to copy c1 array\n", Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMemcpy(Global.gpuData.c2ArrayDevice, c2Array,
                 Global.myMesh->lenum * sizeof(double),
                 cudaMemcpyHostToDevice) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to copy c2 array\n", Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMemcpy(Global.gpuData.c3ArrayDevice, c3Array,
                 Global.myMesh->lenum * sizeof(double),
                 cudaMemcpyHostToDevice) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to copy c3 array\n", Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMemcpy(Global.gpuData.lnidArrayDevice, lnidArray,
                 Global.myMesh->lenum * 8 * sizeof(int32_t),
                 cudaMemcpyHostToDevice) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to copy lnid array\n", Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMemcpy(Global.gpuData.acousticFaceFlagsArrayDevice,
                 acousticFaceFlagsArray, Global.myMesh->lenum * sizeof(uint8_t),
                 cudaMemcpyHostToDevice) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to copy acoustic interface flags\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (acoustic_pressure_enabled()) {
    if (cudaMemcpy(Global.gpuData.acousticElementFlagsArrayDevice,
                   acousticElementFlagsArray,
                   Global.myMesh->lenum * sizeof(uint8_t),
                   cudaMemcpyHostToDevice) != cudaSuccess) {
      fprintf(stderr,
              "Thread %d: Failed to copy acoustic pressure element flags\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    if (cudaMemcpy(Global.gpuData.acousticPressureCoeffArrayDevice,
                   acousticPressureCoeffArray,
                   Global.myMesh->lenum * sizeof(double),
                   cudaMemcpyHostToDevice) != cudaSuccess) {
      fprintf(stderr,
              "Thread %d: Failed to copy acoustic pressure coefficients\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    if (cudaMemcpy(Global.gpuData.acousticPressureBoundaryCoeffArrayDevice,
                   acousticPressureBoundaryCoeffArray,
                   Global.myMesh->lenum * sizeof(double),
                   cudaMemcpyHostToDevice) != cudaSuccess) {
      fprintf(stderr,
              "Thread %d: Failed to copy acoustic pressure boundary coefficients\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
  }

  if (cudaMemcpy(Global.gpuData.g0_shearArrayDevice, g0_shearArray,
                 Global.myMesh->lenum * sizeof(double),
                 cudaMemcpyHostToDevice) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to copy g0 shear array\n", Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMemcpy(Global.gpuData.g1_shearArrayDevice, g1_shearArray,
                 Global.myMesh->lenum * sizeof(double),
                 cudaMemcpyHostToDevice) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to copy g1 shear array\n", Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMemcpy(Global.gpuData.g2_shearArrayDevice, g2_shearArray,
                 Global.myMesh->lenum * sizeof(double),
                 cudaMemcpyHostToDevice) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to copy g2 shear array\n", Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMemcpy(Global.gpuData.g0_kappaArrayDevice, g0_kappaArray,
                 Global.myMesh->lenum * sizeof(double),
                 cudaMemcpyHostToDevice) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to copy g0 kappa array\n", Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMemcpy(Global.gpuData.g1_kappaArrayDevice, g1_kappaArray,
                 Global.myMesh->lenum * sizeof(double),
                 cudaMemcpyHostToDevice) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to copy g1 kappa array\n", Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMemcpy(Global.gpuData.g2_kappaArrayDevice, g2_kappaArray,
                 Global.myMesh->lenum * sizeof(double),
                 cudaMemcpyHostToDevice) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to copy g2 kappa array\n", Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  if (cudaMemcpy(Global.gpuData.b_shearArrayDevice, b_shearArray,
                 Global.myMesh->lenum * sizeof(double),
                 cudaMemcpyHostToDevice) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to copy b shear array\n", Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMemcpy(Global.gpuData.a0_shearArrayDevice, a0_shearArray,
                 Global.myMesh->lenum * sizeof(double),
                 cudaMemcpyHostToDevice) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to copy a0 shear array\n", Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMemcpy(Global.gpuData.a1_shearArrayDevice, a1_shearArray,
                 Global.myMesh->lenum * sizeof(double),
                 cudaMemcpyHostToDevice) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to copy a1 shear array\n", Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMemcpy(Global.gpuData.a2_shearArrayDevice, a2_shearArray,
                 Global.myMesh->lenum * sizeof(double),
                 cudaMemcpyHostToDevice) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to copy a2 shear array\n", Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMemcpy(Global.gpuData.b_kappaArrayDevice, b_kappaArray,
                 Global.myMesh->lenum * sizeof(double),
                 cudaMemcpyHostToDevice) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to copy b kappa array\n", Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMemcpy(Global.gpuData.a0_kappaArrayDevice, a0_kappaArray,
                 Global.myMesh->lenum * sizeof(double),
                 cudaMemcpyHostToDevice) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to copy a0 kappa array\n", Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMemcpy(Global.gpuData.a1_kappaArrayDevice, a1_kappaArray,
                 Global.myMesh->lenum * sizeof(double),
                 cudaMemcpyHostToDevice) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to copy a1 kappa array\n", Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMemcpy(Global.gpuData.a2_kappaArrayDevice, a2_kappaArray,
                 Global.myMesh->lenum * sizeof(double),
                 cudaMemcpyHostToDevice) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to copy a2 kappa array\n", Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  for (int i = 0; i < 3; i++) {
    if (cudaMemcpy(Global.gpuData.mass2_minusaMArrayDevice[i],
                   mass2_minusaMArray[i],
                   Global.myMesh->nharbored * sizeof(solver_float),
                   cudaMemcpyHostToDevice) != cudaSuccess) {
      fprintf(stderr, "Thread %d: Failed to copy mass2 array\n", Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    if (cudaMemcpy(Global.gpuData.mass_minusaMArrayDevice[i],
                   mass_minusaMArray[i],
                   Global.myMesh->nharbored * sizeof(solver_float),
                   cudaMemcpyHostToDevice) != cudaSuccess) {
      fprintf(stderr, "Thread %d: Failed to copy mass array\n", Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
  }
  if (cudaMemcpy(Global.gpuData.mass_simpleArrayDevice, mass_simpleArray,
                 Global.myMesh->nharbored * sizeof(solver_float),
                 cudaMemcpyHostToDevice) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to copy mass simple array\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (acoustic_pressure_enabled()) {
    if (cudaMemcpy(Global.gpuData.pressureMassDevice,
                   Global.mySolver->pressure_mass,
                   Global.myMesh->nharbored * sizeof(solver_float),
                   cudaMemcpyHostToDevice) != cudaSuccess) {
      fprintf(stderr, "Thread %d: Failed to copy pressure mass array\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    if (cudaMemcpy(Global.gpuData.pressureDampDevice,
                   Global.mySolver->pressure_damp,
                   Global.myMesh->nharbored * sizeof(solver_float),
                   cudaMemcpyHostToDevice) != cudaSuccess) {
      fprintf(stderr, "Thread %d: Failed to copy pressure damp array\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    if (cudaMemcpy(Global.gpuData.pressureTm1Device,
                   Global.mySolver->pressure_tm1,
                   Global.myMesh->nharbored * sizeof(solver_float),
                   cudaMemcpyHostToDevice) != cudaSuccess) {
      fprintf(stderr, "Thread %d: Failed to copy pressure tm1 array\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    if (cudaMemcpy(Global.gpuData.pressureTm2Device,
                   Global.mySolver->pressure_tm2,
                   Global.myMesh->nharbored * sizeof(solver_float),
                   cudaMemcpyHostToDevice) != cudaSuccess) {
      fprintf(stderr, "Thread %d: Failed to copy pressure tm2 array\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    if (cudaMemcpy(Global.gpuData.pressureForceDevice,
                   Global.mySolver->pressure_force,
                   Global.myMesh->nharbored * sizeof(solver_float),
                   cudaMemcpyHostToDevice) != cudaSuccess) {
      fprintf(stderr, "Thread %d: Failed to copy pressure force array\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    if (cudaMemcpy(Global.gpuData.pressureFreeSurfaceNodeFlagsDevice,
                   pressureFreeSurfaceNodeFlagsArray,
                   Global.myMesh->nharbored * sizeof(uint8_t),
                   cudaMemcpyHostToDevice) != cudaSuccess) {
      fprintf(stderr,
              "Thread %d: Failed to copy pressure free-surface flags\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    if (Global.gpuData.acousticInterfaceFaceCount > 0) {
      size_t face_count =
          (size_t)Global.gpuData.acousticInterfaceFaceCount;

      if (cudaMemcpy(Global.gpuData.acousticInterfaceElementIdsDevice,
                     Global.acousticInterfaceElementIds,
                     face_count * sizeof(int32_t),
                     cudaMemcpyHostToDevice) != cudaSuccess) {
        fprintf(stderr,
                "Thread %d: Failed to copy acoustic interface element table\n",
                Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }
      if (cudaMemcpy(Global.gpuData.acousticInterfaceFaceIdsDevice,
                     Global.acousticInterfaceFaceIds,
                     face_count * sizeof(uint8_t),
                     cudaMemcpyHostToDevice) != cudaSuccess) {
        fprintf(stderr,
                "Thread %d: Failed to copy acoustic interface face table\n",
                Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }
      if (cudaMemcpy(Global.gpuData.acousticInterfaceFaceAreasDevice,
                     Global.acousticInterfaceFaceAreas,
                     face_count * sizeof(double),
                     cudaMemcpyHostToDevice) != cudaSuccess) {
        fprintf(stderr,
                "Thread %d: Failed to copy acoustic interface face areas\n",
                Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }
      if (cudaMemcpy(Global.gpuData.acousticInterfaceFaceCentersDevice,
                     Global.acousticInterfaceFaceCenters,
                     3 * face_count * sizeof(double),
                     cudaMemcpyHostToDevice) != cudaSuccess) {
        fprintf(stderr,
                "Thread %d: Failed to copy acoustic interface fragment centers\n",
                Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }
      if (cudaMemcpy(Global.gpuData.acousticInterfaceFaceNormalsDevice,
                     Global.acousticInterfaceFaceNormals,
                     3 * face_count * sizeof(double),
                     cudaMemcpyHostToDevice) != cudaSuccess) {
        fprintf(stderr,
                "Thread %d: Failed to copy acoustic interface face normals\n",
                Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }
      if (cudaMemcpy(
              Global.gpuData.acousticInterfaceAdjacentElementIndexesDevice,
              Global.acousticInterfaceAdjacentElementIndexes,
              face_count * sizeof(int32_t), cudaMemcpyHostToDevice) !=
          cudaSuccess) {
        fprintf(stderr,
                "Thread %d: Failed to copy acoustic interface adjacent element indexes\n",
                Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }
      if (cudaMemcpy(Global.gpuData.acousticInterfaceAdjacentCentersDevice,
                     Global.acousticInterfaceAdjacentCenters,
                     3 * face_count * sizeof(double),
                     cudaMemcpyHostToDevice) != cudaSuccess) {
        fprintf(stderr,
                "Thread %d: Failed to copy acoustic interface adjacent centers\n",
                Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }

      /* CPML-aware interface coupling correction arrays. Only allocated
       * when at least one interface face is inside the PML zone (host
       * arrays non-NULL). Two psi memory variables per face, zero-init. */
      if (Global.acousticInterfacePmlFlag != NULL) {
        if (cudaMalloc(
                (void**)&(Global.gpuData.acousticInterfacePmlFlagDevice),
                face_count * sizeof(uint8_t)) != cudaSuccess ||
            cudaMalloc(
                (void**)&(Global.gpuData.acousticInterfacePml_aDevice),
                face_count * sizeof(double)) != cudaSuccess ||
            cudaMalloc(
                (void**)&(Global.gpuData.acousticInterfacePml_bDevice),
                face_count * sizeof(double)) != cudaSuccess ||
            cudaMalloc(
                (void**)&(Global.gpuData.acousticInterfacePml_invKappaDevice),
                face_count * sizeof(double)) != cudaSuccess ||
            cudaMalloc(
                (void**)&(Global.gpuData.acousticInterfacePsiPressureDevice),
                face_count * sizeof(solver_float)) != cudaSuccess ||
            cudaMalloc(
                (void**)&(Global.gpuData.acousticInterfacePsiAccelDevice),
                face_count * sizeof(solver_float)) != cudaSuccess) {
          fprintf(stderr,
                  "Thread %d: Failed to allocate CPML interface-correction "
                  "arrays\n",
                  Global.myID);
          MPI_Abort(MPI_COMM_WORLD, ERROR);
          exit(1);
        }
        cudaMemset(Global.gpuData.acousticInterfacePsiPressureDevice, 0,
                   face_count * sizeof(solver_float));
        cudaMemset(Global.gpuData.acousticInterfacePsiAccelDevice, 0,
                   face_count * sizeof(solver_float));
        if (cudaMemcpy(Global.gpuData.acousticInterfacePmlFlagDevice,
                       Global.acousticInterfacePmlFlag,
                       face_count * sizeof(uint8_t),
                       cudaMemcpyHostToDevice) != cudaSuccess ||
            cudaMemcpy(Global.gpuData.acousticInterfacePml_aDevice,
                       Global.acousticInterfacePml_a,
                       face_count * sizeof(double),
                       cudaMemcpyHostToDevice) != cudaSuccess ||
            cudaMemcpy(Global.gpuData.acousticInterfacePml_bDevice,
                       Global.acousticInterfacePml_b,
                       face_count * sizeof(double),
                       cudaMemcpyHostToDevice) != cudaSuccess ||
            cudaMemcpy(Global.gpuData.acousticInterfacePml_invKappaDevice,
                       Global.acousticInterfacePml_invKappa,
                       face_count * sizeof(double),
                       cudaMemcpyHostToDevice) != cudaSuccess) {
          fprintf(stderr,
                  "Thread %d: Failed to copy CPML interface-correction "
                  "arrays\n",
                  Global.myID);
          MPI_Abort(MPI_COMM_WORLD, ERROR);
          exit(1);
        }
      }
    }
  }

  free(c1Array);
  free(c2Array);
  free(c3Array);
  free(lnidArray);
  free(acousticFaceFlagsArray);
  free(acousticElementFlagsArray);
  free(pressureFreeSurfaceNodeFlagsArray);
  free(acousticPressureCoeffArray);
  free(acousticPressureBoundaryCoeffArray);
  free(pressureFreeSurfaceHost);
  free(g0_shearArray);
  free(g1_shearArray);
  free(g2_shearArray);
  free(g0_kappaArray);
  free(g1_kappaArray);
  free(g2_kappaArray);
  free(b_shearArray);
  free(a0_shearArray);
  free(a1_shearArray);
  free(a2_shearArray);
  free(b_kappaArray);
  free(a0_kappaArray);
  free(a1_kappaArray);
  free(a2_kappaArray);
  free(mass_simpleArray);
  for (int i = 0; i < 3; i++) {
    free(mass2_minusaMArray[i]);
    free(mass_minusaMArray[i]);
  }

  Global.gpu_spec.numbytespci += 19 * (Global.myMesh->lenum * sizeof(double));
  Global.gpu_spec.numbytespci += 1 * (Global.myMesh->lenum * sizeof(int32_t));

    if (Global.myID == 0) {
    fprintf(stdout, "GPU Device memory allocation: initialize convolution\n");
  }

  /* Initialize device displacements */
  cudaMemset(Global.gpuData.tm1Device, 0,
             Global.myMesh->nharbored * sizeof(fvector_t));
  cudaMemset(Global.gpuData.tm2Device, 0,
             Global.myMesh->nharbored * sizeof(fvector_t));
  if ((Param.printStationAccelerations == YES) || acoustic_pressure_enabled()) {
    cudaMemset(Global.gpuData.tm3Device, 0,
               Global.myMesh->nharbored * sizeof(fvector_t));
  }

  /* Initialize device forces */
  cudaMemset(Global.gpuData.forceDevice, 0,
             Global.myMesh->nharbored * sizeof(fvector_t));
  if (acoustic_pressure_enabled()) {
    cudaMemset(Global.gpuData.pressureForceDevice, 0,
               Global.myMesh->nharbored * sizeof(solver_float));
  }

  /* Initialize device convolution */
  cudaMemset(Global.gpuData.conv_shear_1Device, 0,
             8 * Global.myMesh->lenum * sizeof(fvector_t));
  cudaMemset(Global.gpuData.conv_shear_2Device, 0,
             8 * Global.myMesh->lenum * sizeof(fvector_t));
  cudaMemset(Global.gpuData.conv_shear_3Device, 0,
             8 * Global.myMesh->lenum * sizeof(fvector_t));
  cudaMemset(Global.gpuData.conv_kappa_1Device, 0,
             8 * Global.myMesh->lenum * sizeof(fvector_t));
  cudaMemset(Global.gpuData.conv_kappa_2Device, 0,
             8 * Global.myMesh->lenum * sizeof(fvector_t));
  cudaMemset(Global.gpuData.conv_kappa_3Device, 0,
             8 * Global.myMesh->lenum * sizeof(fvector_t));

  /* Initialize damping buffers */
  cudaMemset(Global.gpuData.shearVectorDevice, 0,
             8 * Global.myMesh->lenum * sizeof(fvector_t));
  cudaMemset(Global.gpuData.kappaVectorDevice, 0,
             8 * Global.myMesh->lenum * sizeof(fvector_t));

  if (Global.myID == 0) {
    fprintf(stdout, "GPU Device memory allocation: load constant memory\n");
  }

  if (gpu_copy_constant_symbols(Global.myID, Global.theK1, Global.theK2) != 0) {
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  /* Copy GPU data to device */
  if (cudaMalloc((void**)&(Global.mySolver->gpuDataDevice),
                 sizeof(gpu_data_t)) != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to allocate GPU data memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  cudaMemcpy(Global.mySolver->gpuDataDevice, &Global.gpuData,
             sizeof(gpu_data_t), cudaMemcpyHostToDevice);
  Global.gpu_spec.numbytespci += sizeof(gpu_data_t);

  /* Create streams */
  if (Global.myID == 0) {
    fprintf(stdout, "GPU Device memory allocation: create streams\n");
  }
  for (int i = 0; i < MAX_CUDA_STREAM; i++) {
    cudaStreamCreate(&(Global.mySolver->streams[i]));
  }

  /* Create events */
  for (int i = 0; i < MAX_CUDA_EVENT; i++) {
    cudaEventCreate(&(Global.mySolver->events[i]));
  }

  /* Check for CUDA errors */
  cudaError_t cerror = cudaGetLastError();
  if (cerror != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to setup data in GPU memory\n",
            Global.myID, cudaGetErrorString(cerror));
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  /* Save reference to device data in solver */
  Global.mySolver->gpuData = &Global.gpuData;

  if (Global.myID == 0) {
    fprintf(stdout, "GPU Device memory allocation: end\n");
  }

  /* Calculate kernel block sizes */
  gpu_kernel_t* k;

  k = &Global.gpu_kernel[CUDA_KERNEL_STIFFNESS_FORCE];
  strcpy(k->name, CUDA_KERNEL_NAME_STIFFNESS_FORCE);
  k->handle = (char*)kernelStiffnessCalcLocal;
  k->blocksize = gpu_get_blocksize(Global.mySolver->gpu_spec,
                                   (char*)kernelStiffnessCalcLocal, 0);
  k->gridsize = (Global.myMesh->lenum + k->blocksize - 1) / k->blocksize;
  k->sharedmem = 0;

  k = &Global.gpu_kernel[CUDA_KERNEL_DAMPING_CONV];
  if (Param.theTypeOfDamping == BKT2) {
    strcpy(k->name, CUDA_KERNEL_NAME_DAMPING_CONV_BKT2);
    k->handle = (char*)kernelDampingCalcConvBkt2;
    k->blocksize = gpu_get_blocksize(Global.mySolver->gpu_spec,
                                     (char*)kernelDampingCalcConvBkt2, 0);
  } else {
    strcpy(k->name, CUDA_KERNEL_NAME_DAMPING_CONV);
    k->handle = (char*)kernelDampingCalcConv;
    k->blocksize = gpu_get_blocksize(Global.mySolver->gpu_spec,
                                     (char*)kernelDampingCalcConv, 0);
  }
  k->gridsize = (Global.myMesh->lenum + k->blocksize - 1) / k->blocksize;
  k->sharedmem = 0;

  k = &Global.gpu_kernel[CUDA_KERNEL_DAMPING_FORCE];
  strcpy(k->name, CUDA_KERNEL_NAME_DAMPING_FORCE);
  k->handle = (char*)kernelDampingCalcLocal;
  k->blocksize = gpu_get_blocksize(Global.mySolver->gpu_spec,
                                   (char*)kernelDampingCalcLocal, 0);
  k->gridsize = (Global.myMesh->lenum + k->blocksize - 1) / k->blocksize;
  k->sharedmem = 0;

  k = &Global.gpu_kernel[CUDA_KERNEL_DISPLACEMENT];
  strcpy(k->name, CUDA_KERNEL_NAME_DISPLACEMENT);
  k->handle = (char*)kernelDispCalc;
  k->blocksize =
      gpu_get_blocksize(Global.mySolver->gpu_spec, (char*)kernelDispCalc, 0);
  k->gridsize = (Global.myMesh->nharbored + k->blocksize - 1) / k->blocksize;
    k->sharedmem = 0;

  /* Save reference to kernel data in solver */
  Global.mySolver->gpu_kernel = &(Global.gpu_kernel[0]);

  if (Global.myID == 0) {
    fprintf(stdout, "GPU Device kernel specs: calculate block sizes\n");
    for (int k = 0; k < MAX_CUDA_KERNEL; k++) {
      printf("\t%s  : %d\n", Global.mySolver->gpu_kernel[k].name,
             Global.mySolver->gpu_kernel[k].blocksize);
    }
  }

  /* Calculate kernel register counts */
  if (Global.myID == 0) {
    fprintf(stdout, "GPU Device kernel specs: calculate register counts\n");

    int count;
    for (int k = 0; k < MAX_CUDA_KERNEL; k++) {
      count = gpu_get_reg_count(Global.mySolver->gpu_kernel[k].handle);
      printf("\t%s  : %d\n", Global.mySolver->gpu_kernel[k].name, count);
    }
  }

  return;
}

/**
 * Print the communication schedule of each processor to the given output
 * stream.
 */
static void solver_printstat_to_stream(mysolver_t* solver, FILE* out) {
  /**
   * Enumeration for the indices in the counts array
   * the first 4 values are for the dangling nodes:
   * the last 4 values are for the anchored nodes.
   */
  enum solver_msg_idx_t {
    I_DN_C_P,    /**< 0: Dangling contribution processor count */
    I_DN_C_N,    /**< 1: Dangling contribution node count */
    I_DN_S_P,    /**< 2: Dangling shared processor count */
    I_DN_S_N,    /**< 3: Dangling shared node count */
    I_AN_C_P,    /**< 4: Shared contribution processor count */
    I_AN_C_N,    /**< 5: Shared contribution node count */
    I_AN_S_P,    /**< 6: Shared shared processor count */
    I_AN_S_N,    /**< 7: Shared shared node count */
    I_SCHED_LAST /**< 8: Number of counters */
  };

  int32_t* recv_counts = NULL;
  int32_t send_counts[I_SCHED_LAST];

  /* number of processors this PE communicates with */
  send_counts[I_DN_C_P] = solver->dn_sched->c_count;
  send_counts[I_DN_S_P] = solver->dn_sched->s_count;
  send_counts[I_AN_C_P] = solver->an_sched->c_count;
  send_counts[I_AN_S_P] = solver->an_sched->s_count;

  /* number of nodes this PE communicates */
  send_counts[I_DN_C_N] = messenger_countnodes(solver->dn_sched->first_c);
  send_counts[I_DN_S_N] = messenger_countnodes(solver->dn_sched->first_s);
  send_counts[I_AN_C_N] = messenger_countnodes(solver->an_sched->first_c);
  send_counts[I_AN_S_N] = messenger_countnodes(solver->an_sched->first_s);

  if (Global.myID == 0) {
    /* allocate a buffer in PE 0 to receive stats from all other PEs */
    XMALLOC_VAR_N(recv_counts, int32_t, I_SCHED_LAST* Global.theGroupSize);
  }

  MPI_Gather(send_counts, I_SCHED_LAST, MPI_INT, recv_counts, I_SCHED_LAST,
             MPI_INT, 0, comm_solver);

  if (Global.myID == 0) {
    int pe_id;
    int32_t* pe_counts = recv_counts;

    fprintf(out,
            "# Solver communication schedule summary:\n"
            "#   PE dc_p     dc_n ds_p     ds_n ac_p     ac_n as_p"
            "     as_n total_ncnt\n");

    for (pe_id = 0; pe_id < Global.theGroupSize; pe_id++) {
      /* print a row for a PE, each row has I_SCHED_LAST (8) items
       * besides the pe_id and the total node sum at the end */
      int pe_comm_count = pe_counts[I_DN_C_N] + pe_counts[I_DN_S_N] +
                          pe_counts[I_AN_C_N] + pe_counts[I_AN_S_N];

      fprintf(out, "%6d %4d %8d %4d %8d %4d %8d %4d %8d %10d\n", pe_id,
              pe_counts[I_DN_C_P], pe_counts[I_DN_C_N], pe_counts[I_DN_S_P],
              pe_counts[I_DN_S_N], pe_counts[I_AN_C_P], pe_counts[I_AN_C_N],
              pe_counts[I_AN_S_P], pe_counts[I_AN_S_N], pe_comm_count);

      pe_counts += I_SCHED_LAST; /* advance a row of size I_SCHED_LAST */
    }

    fprintf(out, "\n\n");
    fflush(out);

    free(recv_counts);
  }

  return;
}

/**
 * Print the communication schedule of each processor to the given output
 * stream.
 */
static void solver_printstat(mysolver_t* solver) {
  FILE* stat_out = NULL;

  if (Global.myID == 0) {
    stat_out = hu_fopen(Param.theScheduleStatFilename, "w");
  }

  solver_printstat_to_stream(solver, stat_out);

  if (Global.myID == 0) {
    hu_fclose(stat_out);
    xfree_char(&Param.theScheduleStatFilename);
  }
}

/**
 * solver_delete: Release all the memory associate with Global.mySolver.
 *
 */
static void solver_delete() {
  if (Global.mySolver == NULL) {
    return;
  }

  /* Free device memory */
  // cudaFree(Global.gpuData.elemTableDevice);
  // cudaFree(Global.gpuData.eTableDevice);
  // cudaFree(Global.gpuData.nTableDevice);
  cudaFree(Global.gpuData.tm1Device);
  cudaFree(Global.gpuData.tm2Device);
  if ((Param.printStationAccelerations == YES) || acoustic_pressure_enabled()) {
    cudaFree(Global.gpuData.tm3Device);
  }
  cudaFree(Global.gpuData.forceDevice);
  cudaFree(Global.gpuData.pressureTm1Device);
  cudaFree(Global.gpuData.pressureTm2Device);
  cudaFree(Global.gpuData.pressureForceDevice);
  cudaFree(Global.gpuData.pressureMassDevice);
  cudaFree(Global.gpuData.pressureFreeSurfaceNodeFlagsDevice);
  cudaFree(Global.gpuData.conv_shear_1Device);
  cudaFree(Global.gpuData.conv_shear_2Device);
  cudaFree(Global.gpuData.conv_shear_3Device);
  cudaFree(Global.gpuData.conv_kappa_1Device);
  cudaFree(Global.gpuData.conv_kappa_2Device);
  cudaFree(Global.gpuData.conv_kappa_3Device);
  cudaFree(Global.gpuData.c1ArrayDevice);
  cudaFree(Global.gpuData.c2ArrayDevice);
  cudaFree(Global.gpuData.c3ArrayDevice);
  cudaFree(Global.gpuData.lnidArrayDevice);
  cudaFree(Global.gpuData.acousticElementFlagsArrayDevice);
  cudaFree(Global.gpuData.acousticFaceFlagsArrayDevice);
  cudaFree(Global.gpuData.acousticInterfaceElementIdsDevice);
  cudaFree(Global.gpuData.acousticInterfaceFaceIdsDevice);
  cudaFree(Global.gpuData.acousticInterfaceFilterStateDevice);
  cudaFree(Global.gpuData.acousticInterfaceFaceAreasDevice);
  cudaFree(Global.gpuData.acousticInterfaceFaceCentersDevice);
  cudaFree(Global.gpuData.acousticInterfaceFaceNormalsDevice);
  cudaFree(Global.gpuData.acousticInterfaceAdjacentElementIndexesDevice);
  cudaFree(Global.gpuData.acousticInterfaceAdjacentCentersDevice);
  cudaFree(Global.gpuData.acousticInterfacePmlFlagDevice);
  cudaFree(Global.gpuData.acousticInterfacePml_aDevice);
  cudaFree(Global.gpuData.acousticInterfacePml_bDevice);
  cudaFree(Global.gpuData.acousticInterfacePml_invKappaDevice);
  cudaFree(Global.gpuData.acousticInterfacePsiPressureDevice);
  cudaFree(Global.gpuData.acousticInterfacePsiAccelDevice);
  free(Global.acousticInterfacePmlFlag);
  Global.acousticInterfacePmlFlag = NULL;
  free(Global.acousticInterfacePml_a);
  Global.acousticInterfacePml_a = NULL;
  free(Global.acousticInterfacePml_b);
  Global.acousticInterfacePml_b = NULL;
  free(Global.acousticInterfacePml_invKappa);
  Global.acousticInterfacePml_invKappa = NULL;
  cudaFree(Global.gpuData.acousticPressureCoeffArrayDevice);
  cudaFree(Global.gpuData.acousticPressureBoundaryCoeffArrayDevice);
  cudaFree(Global.gpuData.g0_shearArrayDevice);
  cudaFree(Global.gpuData.g1_shearArrayDevice);
  cudaFree(Global.gpuData.g2_shearArrayDevice);
  cudaFree(Global.gpuData.g0_kappaArrayDevice);
  cudaFree(Global.gpuData.g1_kappaArrayDevice);
  cudaFree(Global.gpuData.g2_kappaArrayDevice);
  cudaFree(Global.gpuData.b_shearArrayDevice);
  cudaFree(Global.gpuData.a0_shearArrayDevice);
  cudaFree(Global.gpuData.a1_shearArrayDevice);
  cudaFree(Global.gpuData.a2_shearArrayDevice);
  cudaFree(Global.gpuData.b_kappaArrayDevice);
  cudaFree(Global.gpuData.a0_kappaArrayDevice);
  cudaFree(Global.gpuData.a1_kappaArrayDevice);
  cudaFree(Global.gpuData.a2_kappaArrayDevice);
  cudaFree(Global.gpuData.shearVectorDevice);
  cudaFree(Global.gpuData.kappaVectorDevice);
  cudaFree(Global.gpuData.mass_simpleArrayDevice);
  // cudaFree(Global.gpuData.matPropsDevice);
  for (int i = 0; i < 3; i++) {
    cudaFree(Global.gpuData.mass2_minusaMArrayDevice[i]);
    cudaFree(Global.gpuData.mass_minusaMArrayDevice[i]);
  }

  /* CPML device arrays — both elastic and acoustic. cudaFree(NULL) and
   * free(NULL) are no-ops, so this is safe in dashpot mode too. */
  cudaFree(Global.gpuData.pmlFlagArrayDevice);
  cudaFree(Global.gpuData.pmlClampNodeMaskDevice);
  cudaFree(Global.gpuData.pmlElasticElemMapperDevice);
  cudaFree(Global.gpuData.pml_elastic_a_xArrayDevice);
  cudaFree(Global.gpuData.pml_elastic_a_yArrayDevice);
  cudaFree(Global.gpuData.pml_elastic_a_zArrayDevice);
  cudaFree(Global.gpuData.pml_elastic_b_xArrayDevice);
  cudaFree(Global.gpuData.pml_elastic_b_yArrayDevice);
  cudaFree(Global.gpuData.pml_elastic_b_zArrayDevice);
  cudaFree(Global.gpuData.pml_elastic_invKappa_xArrayDevice);
  cudaFree(Global.gpuData.pml_elastic_invKappa_yArrayDevice);
  cudaFree(Global.gpuData.pml_elastic_invKappa_zArrayDevice);
  cudaFree(Global.gpuData.pml_elastic_muCoeffArrayDevice);
  cudaFree(Global.gpuData.pml_elastic_lambdaCoeffArrayDevice);
  for (int i = 0; i < 3; i++) {
    for (int jx = 0; jx < 3; jx++) {
      cudaFree(Global.gpuData.pml_elastic_psiArrayDevice[i][jx]);
      cudaFree(Global.gpuData.pml_elastic_divPsiArrayDevice[i][jx]);
      cudaFree(Global.gpuData.pml_elastic_psiShearArrayDevice[i][jx]);
      cudaFree(Global.gpuData.pml_elastic_psiKappaArrayDevice[i][jx]);
    }
  }
  for (int axis = 0; axis < 3; axis++) {
    cudaFree(Global.gpuData.pml_elastic_zetaDecayArrayDevice[axis]);
    cudaFree(Global.gpuData.pml_elastic_zetaGainArrayDevice[axis]);
    cudaFree(Global.gpuData.pml_elastic_memForceArrayDevice[axis]);
  }
  cudaFree(Global.gpuData.pml_elastic_restoreCoeffArrayDevice);
  cudaFree(Global.gpuData.pmlAcousticElemMapperDevice);
  cudaFree(Global.gpuData.pml_acoustic_a_xArrayDevice);
  cudaFree(Global.gpuData.pml_acoustic_a_yArrayDevice);
  cudaFree(Global.gpuData.pml_acoustic_a_zArrayDevice);
  cudaFree(Global.gpuData.pml_acoustic_b_xArrayDevice);
  cudaFree(Global.gpuData.pml_acoustic_b_yArrayDevice);
  cudaFree(Global.gpuData.pml_acoustic_b_zArrayDevice);
  cudaFree(Global.gpuData.pml_acoustic_invKappa_xArrayDevice);
  cudaFree(Global.gpuData.pml_acoustic_invKappa_yArrayDevice);
  cudaFree(Global.gpuData.pml_acoustic_invKappa_zArrayDevice);
  cudaFree(Global.gpuData.pml_acoustic_forceCoeffArrayDevice);
  cudaFree(Global.gpuData.pml_acoustic_restoreCoeffArrayDevice);
  for (int axis = 0; axis < 3; axis++) {
    cudaFree(Global.gpuData.pml_acoustic_psiArrayDevice[axis]);
    cudaFree(Global.gpuData.pml_acoustic_zetaDecayArrayDevice[axis]);
    cudaFree(Global.gpuData.pml_acoustic_zetaGainArrayDevice[axis]);
    cudaFree(Global.gpuData.pml_acoustic_memForceArrayDevice[axis]);
    cudaFree(Global.gpuData.pml_acoustic_zetaPArrayDevice[axis]);
  }
  free(Global.mySolver->pmlElasticElemList);
  Global.mySolver->pmlElasticElemList = NULL;
  for (int i = 0; i < 3; i++) {
    for (int jx = 0; jx < 3; jx++) {
      free(Global.mySolver->pml_elastic_psi[i][jx]);
      Global.mySolver->pml_elastic_psi[i][jx] = NULL;
    }
  }
  free(Global.mySolver->pmlAcousticElemList);
  Global.mySolver->pmlAcousticElemList = NULL;
  for (int axis = 0; axis < 3; axis++) {
    free(Global.mySolver->pml_acoustic_psi[axis]);
    Global.mySolver->pml_acoustic_psi[axis] = NULL;
  }

  cudaFree(Global.mySolver->gpuDataDevice);

  if (g_dampingElemMapperDevice != NULL) {
    cudaFree(g_dampingElemMapperDevice);
    g_dampingElemMapperDevice = NULL;
  }
  g_dampingElemCountGPU = -1;
  g_dampingBkt2CoeffsReady = 0;

  if (g_topoElemMapperDevice != NULL) {
    cudaFree(g_topoElemMapperDevice);
    g_topoElemMapperDevice = NULL;
  }
  if (g_topoConstantsDevice != NULL) {
    cudaFree(g_topoConstantsDevice);
    g_topoConstantsDevice = NULL;
  }
  g_topoElemCountGPU = -1;
  g_topoForceMode = -1;
  g_topoForceModeReported = 0;
  g_topoForcePathReported = 0;

  if (g_sourceNodeListDevice != NULL) {
    cudaFree(g_sourceNodeListDevice);
    g_sourceNodeListDevice = NULL;
  }
  if (g_sourceForcesDevice != NULL) {
    cudaFree(g_sourceForcesDevice);
    g_sourceForcesDevice = NULL;
  }
  g_sourceNodeCapacity = 0;

  if (g_sparseNodeListDevice != NULL) {
    cudaFree(g_sparseNodeListDevice);
    g_sparseNodeListDevice = NULL;
  }
  if (g_sparseForcesDevice != NULL) {
    cudaFree(g_sparseForcesDevice);
    g_sparseForcesDevice = NULL;
  }
  g_sparseActiveCount = 0;
  if (g_sparseNodeListHost != NULL) {
    if (g_sparseHostPinned) {
      cudaFreeHost(g_sparseNodeListHost);
    } else {
      free(g_sparseNodeListHost);
    }
    g_sparseNodeListHost = NULL;
  }
  if (g_sparseForcesHost != NULL) {
    if (g_sparseHostPinned) {
      cudaFreeHost(g_sparseForcesHost);
    } else {
      free(g_sparseForcesHost);
    }
    g_sparseForcesHost = NULL;
  }
  g_sparseHostPinned = 0;
  g_sparseNodeCapacity = 0;
  g_sparseNodeListReady = 0;
  g_sparseSourceCountCached = -1;
  g_sparseTopoCountCached = -1;
  g_sparseCopiedSourceCount = -1;
  g_sparseCopiedTopoCount = -1;
  free(g_rangeUploadStartsHost);
  free(g_rangeUploadLengthsHost);
  g_rangeUploadStartsHost = NULL;
  g_rangeUploadLengthsHost = NULL;
  g_rangeUploadCapacity = 0;
  g_rangeUploadCount = 0;
  g_rangeUploadNodeCount = 0;
  g_rangeUploadReady = 0;
  g_rangeUploadSourceCountCached = -1;
  g_rangeUploadTopoCountCached = -1;

  compute_adjust_cache_free();

  /* Destroy streams */
  for (int i = 0; i < MAX_CUDA_STREAM; i++) {
    cudaStreamDestroy(Global.mySolver->streams[i]);
  }

  /* Destroy events */
  for (int i = 0; i < MAX_CUDA_EVENT; i++) {
    cudaEventDestroy(Global.mySolver->events[i]);
  }

  free(Global.mySolver->eTable);
  free(Global.mySolver->nTable);

  cudaFreeHost(Global.mySolver->tm1);
  cudaFreeHost(Global.mySolver->tm2);
  cudaFreeHost(Global.mySolver->force);
  if (Global.mySolver->tm3 != NULL) {
    free(Global.mySolver->tm3);
  }
  if (Global.mySolver->pressure_tm1 != NULL) {
    cudaFreeHost(Global.mySolver->pressure_tm1);
  }
  if (Global.mySolver->pressure_tm2 != NULL) {
    cudaFreeHost(Global.mySolver->pressure_tm2);
  }
  if (Global.mySolver->pressure_force != NULL) {
    cudaFreeHost(Global.mySolver->pressure_force);
  }
  if (Global.mySolver->pressure_mass != NULL) {
    cudaFreeHost(Global.mySolver->pressure_mass);
  }
  if (Global.mySolver->pressure_damp != NULL) {
    cudaFreeHost(Global.mySolver->pressure_damp);
  }

  if (Param.theTypeOfDamping >= BKT) {
    free(Global.mySolver->conv_shear_1);
    free(Global.mySolver->conv_shear_2);
    free(Global.mySolver->conv_kappa_1);
    free(Global.mySolver->conv_kappa_2);

    if (Param.theTypeOfDamping >= BKT3) {
      free(Global.mySolver->conv_shear_3);
      free(Global.mySolver->conv_kappa_3);
    }
  }

  schedule_delete(Global.mySolver->dn_sched);
  schedule_delete(Global.mySolver->an_sched);

  free(Global.mySolver);
}

static int read_myForces(int32_t timestep, double dt) {
  off_t whereToRead;
  size_t to_read;
  // size_t read_count;
  int interval, i;

  if (get_sourceType() == SRFH) {
    double source_dt = get_srfhdt();
    double T = dt * timestep;

    vector3D_t* aux1 =
        (vector3D_t*)calloc(Global.theNodesLoaded, sizeof(vector3D_t));
    vector3D_t* aux2 =
        (vector3D_t*)calloc(Global.theNodesLoaded, sizeof(vector3D_t));

    if (aux1 == NULL || aux2 == NULL) {
      solver_abort(
          "read_myforces", "memory allocation failed",
          "Cannot allocate memory for Global.myForces interpolation \n");
    }

    interval = floor(T / source_dt);

    /* read first interval */
    whereToRead = ((off_t)sizeof(int32_t)) +
                  Global.theNodesLoaded * sizeof(int32_t) +
                  Global.theNodesLoaded * interval * sizeof(double) * 3;

    hu_fseeko(Global.fpsource, whereToRead, SEEK_SET);

    to_read = Global.theNodesLoaded * 3;
    // read_count = hu_fread(aux1, sizeof(double), to_read, Global.fpsource);
    // read_count = hu_fread(aux2, sizeof(double), to_read, Global.fpsource);
    hu_fread(aux1, sizeof(double), to_read, Global.fpsource);
    hu_fread(aux2, sizeof(double), to_read, Global.fpsource);

    for (i = 0; i < Global.theNodesLoaded; i++) {
      Global.myForces[i].x[0] = aux1[i].x[0] + (aux2[i].x[0] - aux1[i].x[0]) /
                                                   source_dt *
                                                   (T - interval * source_dt);
      Global.myForces[i].x[1] = aux1[i].x[1] + (aux2[i].x[1] - aux1[i].x[1]) /
                                                   source_dt *
                                                   (T - interval * source_dt);
      Global.myForces[i].x[2] = aux1[i].x[2] + (aux2[i].x[2] - aux1[i].x[2]) /
                                                   source_dt *
                                                   (T - interval * source_dt);
    }
    free(aux1);
    free(aux2);
  } else {  // Todo: Doriam says: Optimize other source types (POINT, PLANE,
            // PLANEWITHKINKS) as we did for SRFH if we plan to keep them
    whereToRead = ((off_t)sizeof(int32_t)) +
                  Global.theNodesLoaded * sizeof(int32_t) +
                  Global.theNodesLoaded * timestep * sizeof(double) * 3;

    hu_fseeko(Global.fpsource, whereToRead, SEEK_SET);

    to_read = Global.theNodesLoaded * 3;
    // read_count = hu_fread(Global.myForces, sizeof(double), to_read,
    // Global.fpsource);
    hu_fread(Global.myForces, sizeof(double), to_read, Global.fpsource);
  }

  return 0; /* if we got here everything went OK */
}

/**
 * check the max and min value of displacement.
 */
static void solver_debug_overflow(mysolver_t* solver, mesh_t* mesh, int step) {
  int nindex;
  double max_disp, min_disp, global_max_disp, global_min_disp;

  max_disp = DBL_MIN;
  min_disp = DBL_MAX;

  /* find the min and max X displacement components */
  for (nindex = 0; nindex < mesh->nharbored; nindex++) {
    fvector_t* tm2Disp;
    // n_t *np;

    // np = &solver->nTable[nindex];
    tm2Disp = solver->tm2 + nindex;

    max_disp = (max_disp > tm2Disp->f[0]) ? max_disp : tm2Disp->f[0];
    min_disp = (min_disp < tm2Disp->f[0]) ? min_disp : tm2Disp->f[0];
  }

  /* get global min and max values */
  MPI_Reduce(&min_disp, &global_min_disp, 1, MPI_DOUBLE, MPI_MIN, 0,
             comm_solver);

  MPI_Reduce(&max_disp, &global_max_disp, 1, MPI_DOUBLE, MPI_MAX, 0,
             comm_solver);

  if (Global.myID == 0) {
    printf("Timestep %d: max_dx = %.6f min_dx = %.6f\n", step, global_max_disp,
           global_min_disp);
  }
}

int darray_has_nan_nd(const double* v, int dim, const size_t len[]) {
  HU_ASSERT_PTR_ALWAYS(v);
  HU_ASSERT_PTR_ALWAYS(len);
  HU_ASSERT_ALWAYS(dim > 0);

  size_t* idx = XMALLOC_N(size_t, dim);
  int ret = hu_darray_has_nan_nd(v, dim, len, idx);

  if (ret != 0) {
    int i;

    fputs("WARNING!: Found NAN value at index", stderr);

    for (i = 0; i < dim; i++) {
      fprintf(stderr, " %zu", idx[i]);
    }

    fputc('\n', stderr);
  }

  free(idx);
  idx = NULL;

  return ret;
}

/**
 * Check a fvector_t array for NAN values.
 */
static int fvector_array_has_nan(const fvector_t* f, size_t len,
                                 const char* varname) {
  size_t lengths[2];
  const double* v = (double*)f;

  lengths[0] = len;
  lengths[1] = 3;

  if (darray_has_nan_nd(v, 2, lengths) != 0) {
    if (NULL == varname) {
      varname = "(unknown)";
    }
    fputs("fvector_t varname=", stderr);
    fputs(varname, stderr);
    fputc('\n', stderr);

    return -1;
  }

  return 0;
}

/**
 * Debug function to check the solver structures for NAN.
 * The checked structures are the displacement arrays (tm1 & tm2) and
 * the forces array.
 */
void solver_check_nan(mysolver_t* solver, int node_count, int time_step) {
  HU_ASSERT_PTR_ALWAYS(solver);

  int ret1 = fvector_array_has_nan(solver->tm1, node_count, "tm1");
  int ret2 = fvector_array_has_nan(solver->tm2, node_count, "tm2");
  int ret3 = fvector_array_has_nan(solver->force, node_count, "force");

  if ((ret1 | ret2 | ret3) != 0) {
    hu_solver_abort(__FUNCTION_NAME, NULL, "Found a NAN value at timestep %d",
                    time_step);
  }
}

static void solver_run_init_comm(mysolver_t* solver) {
  /* The int64_t (global node id) is for debug purpose */
  static const int debug_size = (DO_DEBUG) ? sizeof(int64_t) : 0;

  /* properly setup the communication schedule */
  int c_outsize = sizeof(fvector_t) + debug_size; /* force */
  int c_insize = sizeof(fvector_t) + debug_size;  /* displacement */
  int s_outsize = sizeof(fvector_t) + debug_size; /* displacement */
  int s_insize = sizeof(fvector_t) + debug_size;  /* force */

  schedule_prepare(solver->dn_sched, c_outsize, c_insize, s_outsize, s_insize);
  schedule_prepare(solver->an_sched, c_outsize, c_insize, s_outsize, s_insize);

  solver_print_schedules(solver);
}

/**
 * \note Globals used:
 * - Global.myID (read)
 * - Param.theDeltaT (read)
 * - startingStep (read)
 */
/* show a progress bar to make the process less anxious! */
static void solver_update_status(int step, const int start_step) {
  static double lastCheckedTime = 0;
  double interval = 0;
  double CurrTime;

  if (Global.myID == 0) {
    CurrTime = Timer_Value("Total Wall Clock", (TimerKind)0);

    if (lastCheckedTime == 0) {
      lastCheckedTime = CurrTime;
      return;
    }

    monitor_print("*");

    if (step % Param.monitor_stats_rate == 0) {
      interval = CurrTime - lastCheckedTime;
      if (interval > Global.slowestTimeSteps)
        Global.slowestTimeSteps = interval;
      if (interval < Global.fastestTimeSteps)
        Global.fastestTimeSteps = interval;
      monitor_print(
          "     Sim=% 12.6f     ETA=% 6.1f    WC=% 6.1f\n",
          step * Param.theDeltaT,
          ((Param.theTotalSteps - step) / Param.monitor_stats_rate) * interval,
          CurrTime);
      lastCheckedTime = CurrTime;
    }
  }
}

static void solver_write_checkpoint(int step, int start_step) {
  if ((Param.theCheckPointingRate != 0) && (step != start_step) &&
      ((step % Param.theCheckPointingRate) == 0)) {
    checkpoint_write(step, Global.myID, Global.myMesh,
                     Param.theCheckPointingDirOut, Global.theGroupSize,
                     Global.mySolver, comm_solver);
  }
}

/**
 * \note Globals used:
 * - Param.theRate (read)
 */
// static void solver_output_wavefield(int step)
// {
//     if (DO_OUTPUT && (step % Param.theRate == 0))
//     {
//         /* output the current timestep */
//         do_solver_output();
//     }
// }

/**
 * \note Globals used:
 * - thePlanePrintRate (read)
 */
static void solver_output_planes(mysolver_t* solver, int my_id, int step) {
  if (Param.theNumberOfPlanes != 0) {
    if (step % Param.thePlanePrintRate == 0) {
      Timer_Start("Print Planes");
      planes_print(my_id, Param.IO_pool_pe_count, Param.theNumberOfPlanes,
                   solver);
      Timer_Stop("Print Planes");
    }
  }
}

static void solver_output_stations(int step) {
  if (Param.theNumberOfStations != 0) {
    if (step % Param.theStationsPrintRate == 0) {
      Timer_Start("Print Stations");
      interpolate_station_displacements(step);
      Timer_Stop("Print Stations");
    }
  }
}

static void solver_reset_gpu(mysolver_t* solver, mesh_t* mesh) {
  fvector_t* tmp;

  /* Swap displacements to prepare for new time step */
  tmp = solver->gpuData->tm2Device;
  solver->gpuData->tm2Device = solver->gpuData->tm1Device;
  solver->gpuData->tm1Device = tmp;

  /* Update data structure on device */
  cudaMemcpyAsync(solver->gpuDataDevice, solver->gpuData, sizeof(gpu_data_t),
                  cudaMemcpyHostToDevice, solver->streams[CUDA_STREAM_MAIN]);
  solver->gpu_spec->numbytespci += sizeof(gpu_data_t);

  /* Reset forces */
  cudaMemsetAsync(Global.gpuData.forceDevice, 0,
                  Global.myMesh->nharbored * sizeof(fvector_t),
                  solver->streams[CUDA_STREAM_MAIN]);

  /* Pre-calculate data for next time step */
  // solver_precalc_gpu(solver, mesh);

  /* Uncomment for timing tests */
  // cudaStreamSynchronize(solver->streams[CUDA_STREAM_MAIN]);

  return;
}

/**
 * Calculate the nonlinear entities necessary for the next step computation
 * of force correction.
 *
 * \note Globals used:
 * - Param.theNumberOfStations
 * - Param.myNumberOfStations
 * - Param.myStations
 * - Param.theDeltaT
 */
static void solver_nonlinear_state(mysolver_t* solver, mesh_t* mesh,
                                   fmatrix_t k1[8][8], fmatrix_t k2[8][8],
                                   int step) {
  if (Param.includeNonlinearAnalysis == YES) {
    Timer_Start("Compute Non-linear Entities");
    compute_nonlinear_state(mesh, solver, Param.theNumberOfStations,
                            Param.myNumberOfStations, Param.myStations,
                            Param.theDeltaT, step);
    if (get_geostatic_total_time() > 0) {
      compute_bottom_reactions(mesh, solver, k1, k2, step, Param.theDeltaT);
    }
    Timer_Stop("Compute Non-linear Entities");
    if (Param.theNumberOfStations != 0) {
      Timer_Start("Print Stations");
      print_nonlinear_stations(mesh, solver, Param.myStations,
                               Param.myNumberOfStations, Param.theDeltaT, step,
                               Param.theStationsPrintRate);
      Timer_Stop("Print Stations");
    }
  }
}

static void solver_read_source_forces(int step, double dt) {
  Timer_Start("Read My Forces");
  if (Global.theNodesLoaded > 0) {
    read_myForces(step, dt);
  }
  Timer_Stop("Read My Forces");
}

/**
 * TODO: This method uses global variable Param.theDeltaT
 */
static void solver_load_fixedbase_displacements(mysolver_t* solver, int step) {
  Timer_Start("Load Fixedbase Disps");
  if (get_fixedbase_flag() == YES) {
    bldgs_load_fixedbase_disps(solver, Param.theDeltaT, step);
  }
  Timer_Stop("Load Fixedbase Disps");
}

static void solver_init_damping_gpu_mapper(int32_t myID, mesh_t* mesh) {
  int32_t eindex;
  int32_t count;
  int32_t* mapperHost;
  const int need_bkt2_coeffs =
      ((Param.theTypeOfDamping == BKT2) && (g_dampingBkt2CoeffsReady == 0));
  const double bkt2_rmax = Param.theDeltaT;
  int have_bkt2_shear = 0;
  int have_bkt2_kappa = 0;

  if (g_dampingElemCountGPU >= 0) {
    return;
  }

  if (need_bkt2_coeffs) {
    g_dampingBkt2ShearCoef1 = 0.0;
    g_dampingBkt2ShearCoef2 = 0.0;
    g_dampingBkt2ShearExpCoef0 = 0.0;
    g_dampingBkt2ShearCoef3 = 0.0;
    g_dampingBkt2ShearCoef4 = 0.0;
    g_dampingBkt2ShearExpCoef1 = 0.0;
    g_dampingBkt2KappaCoef1 = 0.0;
    g_dampingBkt2KappaCoef2 = 0.0;
    g_dampingBkt2KappaExpCoef0 = 0.0;
    g_dampingBkt2KappaCoef3 = 0.0;
    g_dampingBkt2KappaCoef4 = 0.0;
    g_dampingBkt2KappaExpCoef1 = 0.0;
  }

  count = 0;
  for (eindex = 0; eindex < mesh->lenum; eindex++) {
    edata_t* edata = (edata_t*)mesh->elemTable[eindex].data;

    acoustic_sanitize_element(edata);

    if ((isThisElementNonLinear(mesh, eindex) == NO) &&
        (IsDampingElement(mesh, eindex) == YES)) {
      if (need_bkt2_coeffs) {
        if ((!have_bkt2_shear) &&
            ((edata->a0_shear != 0.0f) || (edata->a1_shear != 0.0f))) {
          const double g0 = edata->g0_shear * bkt2_rmax;
          const double g1 = edata->g1_shear * bkt2_rmax;
          g_dampingBkt2ShearCoef1 = g0 / 2.0;
          g_dampingBkt2ShearCoef2 =
              g_dampingBkt2ShearCoef1 * (1.0 - g0);
          g_dampingBkt2ShearExpCoef0 = exp(-g0);
          g_dampingBkt2ShearCoef3 = g1 / 2.0;
          g_dampingBkt2ShearCoef4 =
              g_dampingBkt2ShearCoef3 * (1.0 - g1);
          g_dampingBkt2ShearExpCoef1 = exp(-g1);
          have_bkt2_shear = 1;
        }
        if ((!have_bkt2_kappa) &&
            ((edata->a0_kappa != 0.0f) || (edata->a1_kappa != 0.0f))) {
          const double g0 = edata->g0_kappa * bkt2_rmax;
          const double g1 = edata->g1_kappa * bkt2_rmax;
          g_dampingBkt2KappaCoef1 = g0 / 2.0;
          g_dampingBkt2KappaCoef2 =
              g_dampingBkt2KappaCoef1 * (1.0 - g0);
          g_dampingBkt2KappaExpCoef0 = exp(-g0);
          g_dampingBkt2KappaCoef3 = g1 / 2.0;
          g_dampingBkt2KappaCoef4 =
              g_dampingBkt2KappaCoef3 * (1.0 - g1);
          g_dampingBkt2KappaExpCoef1 = exp(-g1);
          have_bkt2_kappa = 1;
        }
      }
      count++;
    }
  }

  if (need_bkt2_coeffs) {
    g_dampingBkt2CoeffsReady = 1;
  }

  g_dampingElemCountGPU = count;
  if (count == 0) {
    g_dampingElemMapperDevice = NULL;
    return;
  }

  mapperHost = (int32_t*)malloc(sizeof(int32_t) * count);
  if (mapperHost == NULL) {
    fprintf(stderr, "Thread %d: damping mapper host allocation failed\n", myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  count = 0;
  for (eindex = 0; eindex < mesh->lenum; eindex++) {
    edata_t* edata = (edata_t*)mesh->elemTable[eindex].data;

    acoustic_sanitize_element(edata);

    if ((isThisElementNonLinear(mesh, eindex) == NO) &&
        (IsDampingElement(mesh, eindex) == YES)) {
      mapperHost[count++] = eindex;
    }
  }

  if (cudaMalloc((void**)&g_dampingElemMapperDevice,
                 sizeof(int32_t) * g_dampingElemCountGPU) != cudaSuccess) {
    free(mapperHost);
    fprintf(stderr,
            "Thread %d: damping mapper device allocation failed\n", myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  if (cudaMemcpy(g_dampingElemMapperDevice, mapperHost,
                 sizeof(int32_t) * g_dampingElemCountGPU,
                 cudaMemcpyHostToDevice) != cudaSuccess) {
    free(mapperHost);
    fprintf(stderr, "Thread %d: damping mapper copy to device failed\n", myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  free(mapperHost);
}

/** Compute the force due to the earthquake source. */
static void solver_compute_force_source(int step) {
  Timer_Start("Compute addforces s");
  compute_addforce_s(step);
  Timer_Stop("Compute addforces s");
}

static int32_t topo_force_mode(void) {
  if (g_topoForceMode < 0) {
    const char* env = getenv("HERCULES_TOPO_FORCE_MODE");
    int32_t mode = 1;

    if (env != NULL) {
      if (strcasecmp(env, "gpu") == 0 || strcasecmp(env, "device") == 0 ||
          strcasecmp(env, "on") == 0 || strcasecmp(env, "yes") == 0 ||
          strcasecmp(env, "true") == 0 || strcmp(env, "1") == 0) {
        mode = 1;
      } else if (strcasecmp(env, "cpu") == 0 ||
                 strcasecmp(env, "host") == 0 ||
                 strcasecmp(env, "off") == 0 ||
                 strcasecmp(env, "no") == 0 ||
                 strcasecmp(env, "false") == 0 || strcmp(env, "0") == 0) {
        mode = 0;
      }
    }

    g_topoForceMode = mode;
    if ((Global.myID == 0) && (g_topoForceModeReported == 0)) {
      monitor_print("Topography force mode request: %s\n",
                    g_topoForceMode ? "gpu" : "cpu");
      g_topoForceModeReported = 1;
    }
  }

  return g_topoForceMode;
}

static int32_t use_gpu_topography_fastpath(void) {
  if (Param.includeTopography != YES) {
    return 0;
  }
  if (topo_force_mode() == 0) {
    return 0;
  }
  if (get_topo_meth() == FEM) {
    return 0;
  }
  return 1;
}

static int32_t topography_has_host_force_path(void) {
  if (Param.includeTopography != YES) {
    return 0;
  }
  if (get_topo_meth() == FEM) {
    return 0;
  }
  if (use_gpu_topography_fastpath()) {
    return 0;
  }
  return 1;
}

static int use_gpu_source_fastpath(void) {
  if (Param.drmImplement == YES) {
    return 0;
  }
  if (topography_has_host_force_path()) {
    return 0;
  }
  if ((Param.theTypeOfDamping != BKT) && (Param.theStiffness != EFFECTIVE)) {
    return 0;
  }
  return 1;
}

static int use_sparse_force_upload_topo(void) {
  if (!topography_has_host_force_path()) {
    return 0;
  }
  if (Param.drmImplement == YES) {
    return 0;
  }
  if ((Param.theTypeOfDamping != BKT) && (Param.theStiffness != EFFECTIVE)) {
    return 0;
  }
  return 1;
}

static int32_t need_host_force_upload(void) {
  if (Param.drmImplement == YES) {
    return 1;
  }
  if (!use_gpu_source_fastpath()) {
    return 1;
  }
  if (topography_has_host_force_path()) {
    return 1;
  }
  return 0;
}

static void solver_report_force_paths_once(void) {
  const char* topoMethodName = "none";
  const char* topoEffectivePath = "none";
  int32_t topoRequest = 0;

  if ((Global.myID != 0) || (g_forcePathReported != 0)) {
    return;
  }

  if (Param.includeTopography == YES) {
    topometh_t topoMethod = get_topo_meth();
    topoRequest = topo_force_mode();
    topoMethodName = (topoMethod == FEM) ? "fem" : "vt";
    if (topoMethod == FEM) {
      topoEffectivePath = "none-fem";
    } else if (use_gpu_topography_fastpath()) {
      topoEffectivePath = "gpu";
    } else {
      topoEffectivePath = "cpu";
    }
  }

  monitor_print(
      "Force path: source=%s, host-force-upload=%s, topography-effective=%s "
      "(topography=%s, method=%s, topo-request=%s)\n",
      use_gpu_source_fastpath() ? "gpu" : "cpu",
      need_host_force_upload() ? "yes" : "no", topoEffectivePath,
      (Param.includeTopography == YES) ? "yes" : "no", topoMethodName,
      topoRequest ? "gpu" : "cpu");
  g_forcePathReported = 1;
}

static void solver_compute_force_source_gpu(int32_t step) {
  int32_t count = Global.theNodesLoaded;
  mysolver_t* solver = Global.mySolver;
  mesh_t* mesh = Global.myMesh;
  int blocksize;
  int gridsize;

  (void)step;

  Timer_Start("Compute addforces s");

  cudaMemsetAsync(solver->gpuData->forceDevice, 0,
                  mesh->nharbored * sizeof(fvector_t),
                  solver->streams[CUDA_STREAM_MAIN]);

  if (count <= 0) {
    if (!use_gpu_topography_fastpath()) {
      solver_record_gpu_event(solver, CUDA_EVENT_FORCE_UPLOAD_DONE,
                              "force upload");
    }
    Timer_Stop("Compute addforces s");
    return;
  }

  if (count > g_sourceNodeCapacity) {
    if (g_sourceNodeListDevice != NULL) {
      cudaFree(g_sourceNodeListDevice);
      g_sourceNodeListDevice = NULL;
    }
    if (g_sourceForcesDevice != NULL) {
      cudaFree(g_sourceForcesDevice);
      g_sourceForcesDevice = NULL;
    }

    if (cudaMalloc((void**)&g_sourceNodeListDevice,
                   count * sizeof(int32_t)) != cudaSuccess) {
      fprintf(stderr, "Thread %d: source node list device alloc failed\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    if (cudaMalloc((void**)&g_sourceForcesDevice,
                   count * sizeof(vector3D_t)) != cudaSuccess) {
      fprintf(stderr, "Thread %d: source force device alloc failed\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

    g_sourceNodeCapacity = count;

    cudaMemcpyAsync(g_sourceNodeListDevice, Global.theNodesLoadedList,
                    count * sizeof(int32_t), cudaMemcpyHostToDevice,
                    solver->streams[CUDA_STREAM_MAIN]);
    solver->gpu_spec->numbytespci += count * sizeof(int32_t);
  }

  cudaMemcpyAsync(g_sourceForcesDevice, Global.myForces,
                  count * sizeof(vector3D_t), cudaMemcpyHostToDevice,
                  solver->streams[CUDA_STREAM_MAIN]);
  solver->gpu_spec->numbytespci += count * sizeof(vector3D_t);

  blocksize = 256;
  gridsize = (count + blocksize - 1) / blocksize;
  kernelSourceForceScatter<<<gridsize, blocksize, 0,
                             solver->streams[CUDA_STREAM_MAIN]>>>(
      count, g_sourceNodeListDevice, g_sourceForcesDevice,
      solver->gpuData->forceDevice, Param.theDeltaTSquared);

  cudaError_t cerror = cudaGetLastError();
  if (cerror != cudaSuccess) {
    fprintf(stderr, "Thread %d: source force scatter kernel - %s\n",
            Global.myID, cudaGetErrorString(cerror));
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  if (!use_gpu_topography_fastpath()) {
    solver_record_gpu_event(solver, CUDA_EVENT_FORCE_UPLOAD_DONE,
                            "force upload");
  }

  Timer_Stop("Compute addforces s");
}

/** Compute the force due to element stiffness matrices. */
static void solver_compute_force_stiffness(mysolver_t* solver, mesh_t* mesh,
                                           fmatrix_t k1[8][8],
                                           fmatrix_t k2[8][8]) {
  Timer_Start("Compute addforces e");
  if ((Param.theTypeOfDamping < BKT)) {
    if (Param.theStiffness == EFFECTIVE) {
      /* compute_addforce_effective(mesh, solver); */
      compute_addforce_effective_gpu(Global.myID, mesh, solver);
    } else if (Param.theStiffness == CONVENTIONAL) {
      compute_addforce_conventional(mesh, solver, k1, k2);
    }
  }

  /* CPML elastic-side force (ADE memory-variable update fused with the
   * modified-strain Galerkin assembly). Operates on its own mapper, so
   * the non-PML kernels above have already early-exited on these
   * elements via the pmlFlagArrayDevice check. */
  if (Param.theAbcType == ABC_CPML &&
      solver->gpuData->pmlElasticElemCount > 0) {
    const int32_t pmlN = solver->gpuData->pmlElasticElemCount;
    const int blocksize = 128;
    const int gridsize = (pmlN + blocksize - 1) / blocksize;
    cudaGetLastError();
    kernelPMLElasticStep<<<gridsize, blocksize, 0,
                           solver->streams[CUDA_STREAM_MAIN]>>>(
        mesh->lenum, solver->gpuDataDevice, pmlN);
    cudaError_t cerror = cudaGetLastError();
    if (cerror != cudaSuccess) {
      fprintf(stderr, "Thread %d: PML elastic step kernel - %s\n",
              Global.myID, cudaGetErrorString(cerror));
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
  }

  solver_record_gpu_event(solver, CUDA_EVENT_FORCE_STIFFNESS_DONE,
                          "force stiffness");
  Timer_Stop("Compute addforces e");
}

static void solver_compute_force_damping_gpu(int32_t myID, mysolver_t* solver,
                                             mesh_t* mesh, double theFreq,
                                             double theDeltaT,
                                             double theDeltaTSquared) {
  (void)theDeltaTSquared;

  solver_init_damping_gpu_mapper(myID, mesh);
  if (g_dampingElemCountGPU <= 0) {
    solver_record_gpu_event(solver, CUDA_EVENT_FORCE_DAMPING_CONV_DONE,
                            "force damping conv");
    solver_record_gpu_event(solver, CUDA_EVENT_FORCE_DAMPING_LOCAL_DONE,
                            "force damping local");
    return;
  }

  double rmax;
  if (Param.theTypeOfDamping == BKT) {
    rmax = 2. * M_PI * theFreq * theDeltaT;
  } else {
    rmax = theDeltaT;
  }

  int blocksize = solver->gpu_kernel[CUDA_KERNEL_DAMPING_CONV].blocksize;
  int gridsize = solver->gpu_kernel[CUDA_KERNEL_DAMPING_CONV].gridsize;
  int sharedmem = solver->gpu_kernel[CUDA_KERNEL_DAMPING_CONV].sharedmem;

  cudaGetLastError();
  if (Param.theTypeOfDamping == BKT2) {
    kernelDampingCalcConvBkt2<<<gridsize, blocksize, sharedmem,
                                solver->streams[CUDA_STREAM_MAIN]>>>(
        mesh->lenum, mesh->lenum * 8, solver->gpuDataDevice, rmax,
        g_dampingBkt2ShearCoef1, g_dampingBkt2ShearCoef2,
        g_dampingBkt2ShearExpCoef0, g_dampingBkt2ShearCoef3,
        g_dampingBkt2ShearCoef4, g_dampingBkt2ShearExpCoef1,
        g_dampingBkt2KappaCoef1, g_dampingBkt2KappaCoef2,
        g_dampingBkt2KappaExpCoef0, g_dampingBkt2KappaCoef3,
        g_dampingBkt2KappaCoef4, g_dampingBkt2KappaExpCoef1,
        g_dampingElemCountGPU, g_dampingElemMapperDevice);
  } else {
    kernelDampingCalcConv<<<gridsize, blocksize, sharedmem,
                            solver->streams[CUDA_STREAM_MAIN]>>>(
        mesh->lenum, mesh->lenum * 8, solver->gpuDataDevice, rmax,
        (int32_t)Param.theTypeOfDamping, g_dampingElemCountGPU,
        g_dampingElemMapperDevice);
  }

  cudaError_t cerror = cudaGetLastError();
  if (cerror != cudaSuccess) {
    fprintf(stderr, "Thread %d: Calc damping conv kernel - %s\n", myID,
            cudaGetErrorString(cerror));
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  solver_record_gpu_event(solver, CUDA_EVENT_FORCE_DAMPING_CONV_DONE,
                          "force damping conv");

  solver->gpu_spec->numflops +=
      kernel_flops_per_thread(FLOP_CALCCONV_KERNEL) * g_dampingElemCountGPU;
  solver->gpu_spec->numbytes +=
      kernel_mem_per_thread(FLOP_CALCCONV_KERNEL) * g_dampingElemCountGPU;

  blocksize = solver->gpu_kernel[CUDA_KERNEL_DAMPING_FORCE].blocksize;
  gridsize = solver->gpu_kernel[CUDA_KERNEL_DAMPING_FORCE].gridsize;
  sharedmem = solver->gpu_kernel[CUDA_KERNEL_DAMPING_FORCE].sharedmem;

  cudaGetLastError();
  kernelDampingCalcLocal<<<gridsize, blocksize, sharedmem,
                           solver->streams[CUDA_STREAM_MAIN]>>>(
      mesh->lenum, solver->gpuDataDevice, g_dampingElemCountGPU,
      g_dampingElemMapperDevice);

  cerror = cudaGetLastError();
  if (cerror != cudaSuccess) {
    fprintf(stderr, "Thread %d: Calc damping local kernel - %s\n", myID,
            cudaGetErrorString(cerror));
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  solver_record_gpu_event(solver, CUDA_EVENT_FORCE_DAMPING_LOCAL_DONE,
                          "force damping local");

  solver->gpu_spec->numflops +=
      kernel_flops_per_thread(FLOP_DAMPING_KERNEL) * g_dampingElemCountGPU;
  solver->gpu_spec->numbytes +=
      kernel_mem_per_thread(FLOP_DAMPING_KERNEL) * g_dampingElemCountGPU;
}

/** Compute contribution of damping to the force vector */
static void solver_compute_force_damping(mysolver_t* solver, mesh_t* mesh,
                                         fmatrix_t k1[8][8],
                                         fmatrix_t k2[8][8]) {
  Timer_Start("Damping addforce");

  if (Param.theTypeOfDamping == RAYLEIGH || Param.theTypeOfDamping == MASS) {
    /* Dorian says: The Mass and Rayleigh damping are implicitly considered in
     * the nodal mass and the stiffness computation module  */

    // damping_addforce(Global.myMesh, Global.mySolver, Global.theK1,
    // Global.theK2);
    solver_record_gpu_event(solver, CUDA_EVENT_FORCE_DAMPING_CONV_DONE,
                            "force damping conv");
    solver_record_gpu_event(solver, CUDA_EVENT_FORCE_DAMPING_LOCAL_DONE,
                            "force damping local");
  } else if (Param.theTypeOfDamping >= BKT) {
    solver_compute_force_damping_gpu(Global.myID, solver, mesh, Param.theFreq,
                                     Param.theDeltaT,
                                     Param.theDeltaTSquared);

    // conv_and_bktForceCombined(Global.myMesh, Global.mySolver, Param.theFreq,
    // Param.theDeltaT, Param.theDeltaTSquared, Param.theTypeOfDamping);
  } else {
    //      /* Should never reach this point */
    //        solver_abort( __FUNCTION_NAME, NULL,
    //                "Unknown damping type: %d\n",
    //                Param.theTypeOfDamping);
    solver_record_gpu_event(solver, CUDA_EVENT_FORCE_DAMPING_CONV_DONE,
                            "force damping conv");
    solver_record_gpu_event(solver, CUDA_EVENT_FORCE_DAMPING_LOCAL_DONE,
                            "force damping local");
  }

  Timer_Stop("Damping addforce");
}

/**
 * Compute the nonlinear contribution to the force.
 * \param deltaT2 Delta t^2 (i.e., squared).
 */
static void solver_compute_force_nonlinear(mysolver_t* solver, mesh_t* mesh,
                                           double deltaT2) {
  if (Param.includeNonlinearAnalysis == YES) {
    Timer_Start("Compute addforces Non-linear");
    compute_addforce_nl(mesh, solver, deltaT2);
    Timer_Stop("Compute addforces Non-linear");
  }
}

/**
 * Compute the plane waves contribution to the force.
 */
static void solver_compute_force_planewaves(mesh_t* myMesh,
                                            mysolver_t* mySolver,
                                            double theDeltaT, int step,
                                            fmatrix_t (*theK1)[8],
                                            fmatrix_t (*theK2)[8]) {
  if (Param.includeIncidentPlaneWaves == YES) {
    Timer_Start("Compute addforces Incident Plane Waves");
    compute_addforce_PlaneWaves(myMesh, mySolver, theDeltaT, step, theK1,
                                theK2);
    Timer_Stop("Compute addforces Incident Plane Waves");
  }
}

/**
 * Compute the topography contribution to the force.
 * \param deltaT2 Delta t^2 (i.e., squared).
 */
static void solver_init_topography_gpu_mapper(int32_t myID, mesh_t* mesh) {
  int32_t topoCount;
  const int32_t* topoMapperHost;
  const topoconstants_t* topoConstantsHost;
  int32_t* filteredMapperHost;
  topoconstants_t* filteredConstantsHost;
  int32_t writeIdx = 0;
  int32_t skippedNonlinear = 0;
  int32_t skippedTopoBkt = 0;
  noyesflag_t topoNonlinFlag;

  if (g_topoElemCountGPU >= 0) {
    return;
  }

  topoCount = topo_get_elements_count();
  topoMapperHost = topo_get_elements_mapping();
  topoConstantsHost = topo_get_constants();
  topoNonlinFlag = get_topo_nonlin_flag();

  if ((topoCount <= 0) || (topoMapperHost == NULL) ||
      (topoConstantsHost == NULL)) {
    g_topoElemCountGPU = 0;
    return;
  }

  filteredMapperHost = (int32_t*)malloc((size_t)topoCount * sizeof(int32_t));
  filteredConstantsHost =
      (topoconstants_t*)malloc((size_t)topoCount * sizeof(topoconstants_t));
  if ((filteredMapperHost == NULL) || (filteredConstantsHost == NULL)) {
    free(filteredMapperHost);
    free(filteredConstantsHost);
    fprintf(stderr,
            "Thread %d: topography GPU mapper host allocation failed\n",
            myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  for (int32_t topoIdx = 0; topoIdx < topoCount; topoIdx++) {
    int32_t eindex = topoMapperHost[topoIdx];
    edata_t* edata = (edata_t*)mesh->elemTable[eindex].data;
    const topoconstants_t* topo_ec = topoConstantsHost + topoIdx;

    if ((topo_ec->isTopoNonlin != 0) && (topoNonlinFlag == YES)) {
      skippedNonlinear++;
      continue;
    }
    if (edata->topoBkt != 0) {
      skippedTopoBkt++;
      continue;
    }

    filteredMapperHost[writeIdx] = eindex;
    filteredConstantsHost[writeIdx] = *topo_ec;
    writeIdx++;
  }

  g_topoElemCountGPU = writeIdx;
  if (g_topoElemCountGPU == 0) {
    free(filteredMapperHost);
    free(filteredConstantsHost);
    if ((Global.myID == 0) && (g_topoForcePathReported == 0)) {
      monitor_print(
          "Topography force path: gpu (eligible=0, total=%d, nonlinear-skipped=%d, topo-bkt-skipped=%d)\n",
          topoCount, skippedNonlinear, skippedTopoBkt);
      g_topoForcePathReported = 1;
    }
    return;
  }

  if (cudaMalloc((void**)&g_topoElemMapperDevice,
                 (size_t)g_topoElemCountGPU * sizeof(int32_t)) !=
      cudaSuccess) {
    free(filteredMapperHost);
    free(filteredConstantsHost);
    fprintf(stderr,
            "Thread %d: topography GPU mapper device allocation failed\n",
            myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMalloc((void**)&g_topoConstantsDevice,
                 (size_t)g_topoElemCountGPU * sizeof(topoconstants_t)) !=
      cudaSuccess) {
    free(filteredMapperHost);
    free(filteredConstantsHost);
    fprintf(stderr,
            "Thread %d: topography GPU constants device allocation failed\n",
            myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  if (cudaMemcpy(g_topoElemMapperDevice, filteredMapperHost,
                 (size_t)g_topoElemCountGPU * sizeof(int32_t),
                 cudaMemcpyHostToDevice) != cudaSuccess) {
    free(filteredMapperHost);
    free(filteredConstantsHost);
    fprintf(stderr, "Thread %d: topography GPU mapper copy failed\n", myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  if (cudaMemcpy(g_topoConstantsDevice, filteredConstantsHost,
                 (size_t)g_topoElemCountGPU * sizeof(topoconstants_t),
                 cudaMemcpyHostToDevice) != cudaSuccess) {
    free(filteredMapperHost);
    free(filteredConstantsHost);
    fprintf(stderr, "Thread %d: topography GPU constants copy failed\n", myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  free(filteredMapperHost);
  free(filteredConstantsHost);

  if ((Global.myID == 0) && (g_topoForcePathReported == 0)) {
    monitor_print(
        "Topography force path: gpu (eligible=%d, total=%d, nonlinear-skipped=%d, topo-bkt-skipped=%d)\n",
        g_topoElemCountGPU, topoCount, skippedNonlinear, skippedTopoBkt);
    g_topoForcePathReported = 1;
  }
}

static void solver_compute_force_topography_gpu(mysolver_t* solver, mesh_t* mesh,
                                                double deltaT2) {
  int blocksize;
  int gridsize;
  cudaError_t cerror;

  solver_init_topography_gpu_mapper(Global.myID, mesh);

  if (g_topoElemCountGPU <= 0) {
    solver_record_gpu_event(solver, CUDA_EVENT_FORCE_UPLOAD_DONE,
                            "force upload");
    return;
  }

  blocksize = 128;
  gridsize = (g_topoElemCountGPU + blocksize - 1) / blocksize;

  cudaGetLastError();
  kernelTopoForceEffective<<<gridsize, blocksize, 0,
                             solver->streams[CUDA_STREAM_MAIN]>>>(
      mesh->lenum, solver->gpuDataDevice, g_topoElemCountGPU,
      g_topoElemMapperDevice, g_topoConstantsDevice, deltaT2);

  cerror = cudaGetLastError();
  if (cerror != cudaSuccess) {
    fprintf(stderr, "Thread %d: topography force kernel - %s\n", Global.myID,
            cudaGetErrorString(cerror));
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  solver_record_gpu_event(solver, CUDA_EVENT_FORCE_UPLOAD_DONE,
                          "force upload");
}

static void solver_compute_force_topography(mysolver_t* solver, mesh_t* mesh,
                                            double deltaT2) {
  if (Param.includeTopography == YES) {
    Timer_Start("Compute addforces Topography");
    if (use_gpu_topography_fastpath()) {
      solver_compute_force_topography_gpu(solver, mesh, deltaT2);
    } else {
      compute_addforce_topoEffective(mesh, solver, deltaT2);
    }
    Timer_Stop("Compute addforces Topography");
  }
}

static void solver_compute_force_gravity(mysolver_t* solver, mesh_t* mesh,
                                         int step) {
  if (Param.includeNonlinearAnalysis == YES) {
    Timer_Start("Compute addforces gravity");
    if (get_geostatic_total_time() > 0) {
      compute_addforce_gravity(mesh, solver, step, Param.theDeltaT);
    }
    Timer_Stop("Compute addforces gravity");
  }
}

/** Send the forces on dangling nodes to their owner processors */
static void solver_send_force_dangling(mysolver_t* solver) {
  //    HU_COND_GLOBAL_BARRIER( Param.theTimingBarriersFlag );

  Timer_Start("1st schedule send data (contribution)");
  schedule_senddata(solver->dn_sched, solver->force,
                    sizeof(fvector_t) / sizeof(solver_float), CONTRIBUTION,
                    DN_FORCE_MSG);
  Timer_Stop("1st schedule send data (contribution)");
}

static void solver_adjust_forces(mysolver_t* solver) {
  Timer_Start("1st compute adjust (distribution)");
  /* Distribute the forces to LOCAL anchored nodes */
  compute_adjust(solver->force, sizeof(fvector_t) / sizeof(solver_float),
                 DISTRIBUTION);
  Timer_Stop("1st compute adjust (distribution)");
}

static void solver_send_force_anchored(mysolver_t* solver) {
  //    HU_COND_GLOBAL_BARRIER( Param.theTimingBarriersFlag );

  Timer_Start("2nd schedule send data (contribution)");
  /* Send the forces on anchored nodes to their owner processors */
  schedule_senddata(solver->an_sched, solver->force,
                    sizeof(fvector_t) / sizeof(solver_float), CONTRIBUTION,
                    AN_FORCE_MSG);
  Timer_Stop("2nd schedule send data (contribution)");
}

static void solver_load_pressure_gpu(mysolver_t* solver, mesh_t* mesh) {
  if (!acoustic_pressure_enabled()) {
    return;
  }

  cudaMemcpyAsync(solver->gpuData->pressureTm1Device, solver->pressure_tm1,
                  mesh->nharbored * sizeof(solver_float),
                  cudaMemcpyHostToDevice, solver->streams[CUDA_STREAM_MAIN]);
  cudaMemcpyAsync(solver->gpuData->pressureTm2Device, solver->pressure_tm2,
                  mesh->nharbored * sizeof(solver_float),
                  cudaMemcpyHostToDevice, solver->streams[CUDA_STREAM_MAIN]);
  cudaMemsetAsync(solver->gpuData->pressureForceDevice, 0,
                  mesh->nharbored * sizeof(solver_float),
                  solver->streams[CUDA_STREAM_MAIN]);
  solver->gpu_spec->numbytespci +=
      2 * mesh->nharbored * (int64_t)sizeof(solver_float);
}

static void solver_swap_pressure_buffers(mysolver_t* solver) {
  solver_float* tmp;

  if (!acoustic_pressure_enabled()) {
    return;
  }

  tmp = solver->pressure_tm2;
  solver->pressure_tm2 = solver->pressure_tm1;
  solver->pressure_tm1 = tmp;
}

static void solver_compute_force_acoustic_pressure(mysolver_t* solver,
                                                   mesh_t* mesh) {
  int blocksize;
  int32_t face_count;
  int gridsize;
  cudaError_t cerror;
  (void)mesh;

  if ((!acoustic_pressure_enabled()) || (g_acousticPressureElementCount <= 0)) {
    return;
  }

  Timer_Start("Acoustic pressure traction");
  face_count = solver->gpuData->acousticInterfaceFaceCount;
  if (face_count <= 0) {
    Timer_Stop("Acoustic pressure traction");
    return;
  }

  blocksize = 256;
  gridsize = (face_count + blocksize - 1) / blocksize;
  kernelAcousticPressureDisplacementCoupling<<<
      gridsize, blocksize, 0, solver->streams[CUDA_STREAM_MAIN]>>>(
      face_count, solver->gpuDataDevice, Param.theDeltaTSquared,
      Param.theDeltaT);

  cerror = cudaGetLastError();
  if (cerror != cudaSuccess) {
    fprintf(stderr, "Thread %d: acoustic pressure traction kernel - %s\n",
            Global.myID, cudaGetErrorString(cerror));
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  Timer_Stop("Acoustic pressure traction");
}

static void solver_run_acoustic_pressure_step(mysolver_t* solver,
                                              mesh_t* mesh) {
  int blocksize;
  int32_t face_count;
  int gridsize;
  cudaError_t cerror;

  if ((!acoustic_pressure_enabled()) || (g_acousticPressureElementCount <= 0)) {
    return;
  }

  Timer_Start("Acoustic Pressure Step");

  blocksize = 256;
  gridsize = (mesh->lenum + blocksize - 1) / blocksize;
  cudaMemsetAsync(solver->gpuData->pressureForceDevice, 0,
                  mesh->nharbored * sizeof(solver_float),
                  solver->streams[CUDA_STREAM_MAIN]);
  kernelAcousticPressureForce<<<gridsize, blocksize, 0,
                                solver->streams[CUDA_STREAM_MAIN]>>>(
      mesh->lenum, solver->gpuDataDevice, Param.theDeltaTSquared);

  cerror = cudaGetLastError();
  if (cerror != cudaSuccess) {
    fprintf(stderr, "Thread %d: acoustic pressure force kernel - %s\n",
            Global.myID, cudaGetErrorString(cerror));
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  /* CPML acoustic-side force (ADE update + modified Laplacian assembly). */
  if (Param.theAbcType == ABC_CPML &&
      solver->gpuData->pmlAcousticElemCount > 0) {
    const int32_t pmlN = solver->gpuData->pmlAcousticElemCount;
    const int pmlBlock = 128;
    const int pmlGrid = (pmlN + pmlBlock - 1) / pmlBlock;
    kernelPMLAcousticStep<<<pmlGrid, pmlBlock, 0,
                            solver->streams[CUDA_STREAM_MAIN]>>>(
        mesh->lenum, solver->gpuDataDevice, pmlN);
    cerror = cudaGetLastError();
    if (cerror != cudaSuccess) {
      fprintf(stderr, "Thread %d: PML acoustic step kernel - %s\n",
              Global.myID, cudaGetErrorString(cerror));
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
  }

  face_count = solver->gpuData->acousticInterfaceFaceCount;
  if (face_count > 0) {
    gridsize = (face_count + blocksize - 1) / blocksize;
    kernelAcousticPressureInterfaceForce<<<
        gridsize, blocksize, 0, solver->streams[CUDA_STREAM_MAIN]>>>(
        face_count, solver->gpuDataDevice, Param.theDeltaTSquared,
        Param.theDeltaT);

    cerror = cudaGetLastError();
    if (cerror != cudaSuccess) {
      fprintf(stderr,
              "Thread %d: acoustic pressure interface force kernel - %s\n",
              Global.myID, cudaGetErrorString(cerror));
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
  }

  cudaMemcpyAsync(solver->pressure_force, solver->gpuData->pressureForceDevice,
                  mesh->nharbored * sizeof(solver_float),
                  cudaMemcpyDeviceToHost, solver->streams[CUDA_STREAM_MAIN]);
  solver->gpu_spec->numbytespci +=
      mesh->nharbored * (int64_t)sizeof(solver_float);
  cudaStreamSynchronize(solver->streams[CUDA_STREAM_MAIN]);

  schedule_senddata(solver->dn_sched, solver->pressure_force, 1, CONTRIBUTION,
                    DN_PRESSURE_FORCE_MSG);
  compute_adjust(solver->pressure_force, 1, DISTRIBUTION);
  schedule_senddata(solver->an_sched, solver->pressure_force, 1, CONTRIBUTION,
                    AN_PRESSURE_FORCE_MSG);

  cudaMemcpyAsync(solver->gpuData->pressureForceDevice, solver->pressure_force,
                  mesh->nharbored * sizeof(solver_float),
                  cudaMemcpyHostToDevice, solver->streams[CUDA_STREAM_MAIN]);
  solver->gpu_spec->numbytespci +=
      mesh->nharbored * (int64_t)sizeof(solver_float);

  blocksize = 256;
  gridsize = (mesh->nharbored + blocksize - 1) / blocksize;
  kernelAcousticPressureUpdate<<<gridsize, blocksize, 0,
                                 solver->streams[CUDA_STREAM_MAIN]>>>(
      mesh->nharbored, solver->gpuDataDevice);

  cerror = cudaGetLastError();
  if (cerror != cudaSuccess) {
    fprintf(stderr, "Thread %d: acoustic pressure update kernel - %s\n",
            Global.myID, cudaGetErrorString(cerror));
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  cudaMemcpyAsync(solver->pressure_tm2, solver->gpuData->pressureTm2Device,
                  mesh->nharbored * sizeof(solver_float),
                  cudaMemcpyDeviceToHost, solver->streams[CUDA_STREAM_MAIN]);
  solver->gpu_spec->numbytespci +=
      mesh->nharbored * (int64_t)sizeof(solver_float);
  cudaStreamSynchronize(solver->streams[CUDA_STREAM_MAIN]);

  memset(solver->pressure_force, 0,
         mesh->nharbored * sizeof(solver_float));
  schedule_senddata(solver->an_sched, solver->pressure_tm2, 1, SHARING,
                    AN_PRESSURE_MSG);
  compute_adjust(solver->pressure_tm2, 1, ASSIGNMENT);
  schedule_senddata(solver->dn_sched, solver->pressure_tm2, 1, SHARING,
                    DN_PRESSURE_MSG);

  Timer_Stop("Acoustic Pressure Step");
}

/** Compute new displacements of my harbored nodes */
static void solver_compute_displacement(mysolver_t* solver, mesh_t* mesh) {
  lnid_t nindex;

  Timer_Start("Compute new displacement");
  solver_snapshot_station_accel_history(solver);
  for (nindex = 0; nindex < mesh->nharbored; nindex++) {
    const n_t* np = &solver->nTable[nindex];
    fvector_t nodalForce = solver->force[nindex];
    const fvector_t* tm1Disp = solver->tm1 + nindex;
    fvector_t* tm2Disp = solver->tm2 + nindex;

    /* total nodal forces */
    nodalForce.f[0] += np->mass2_minusaM[0] * tm1Disp->f[0] -
                       np->mass_minusaM[0] * tm2Disp->f[0];
    nodalForce.f[1] += np->mass2_minusaM[1] * tm1Disp->f[1] -
                       np->mass_minusaM[1] * tm2Disp->f[1];
    nodalForce.f[2] += np->mass2_minusaM[2] * tm1Disp->f[2] -
                       np->mass_minusaM[2] * tm2Disp->f[2];

    /* overwrite tm2 */
    /* mass sanity check */
    if (np->mass_simple < 0) {
      fprintf(stderr, "Negative lumped mass m=%f. Node ID=%d \n",
              np->mass_simple, nindex);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

    /* Dorian. correct Displacement */
    // TODO: Think of a better place for this
    if (np->mass_simple != 0) {
      tm2Disp->f[0] = nodalForce.f[0] / np->mass_simple;
      tm2Disp->f[1] = nodalForce.f[1] / np->mass_simple;
      tm2Disp->f[2] = nodalForce.f[2] / np->mass_simple;
    } else {
      tm2Disp->f[0] = 0.0;
      tm2Disp->f[1] = 0.0;
      tm2Disp->f[2] = 0.0;
    }

  } /* for (nindex ...): all my harbored nodes */

  /* zero out the force vector for all nodes */
  memset(solver->force, 0, sizeof(fvector_t) * mesh->nharbored);

  Timer_Stop("Compute new displacement");
}

static int32_t displacement_boundary_sync_safe(void) {
  if (Param.theUseCheckPoint == 1 || Param.theCheckPointingRate != 0) {
    return 0;
  }
  if (Param.theNumberOfStations != 0 || Param.theNumberOfPlanes != 0) {
    return 0;
  }
  if (Param.drmImplement == YES || Param.includeNonlinearAnalysis == YES ||
      Param.includeIncidentPlaneWaves == YES) {
    return 0;
  }
  if ((get_fixedbase_flag() == YES) || (get_geostatic_total_time() > 0)) {
    return 0;
  }
  return 1;
}

static int32_t displacement_sync_mode(void) {
  if (g_dispSyncMode < 0) {
    const char* env = getenv("HERCULES_DISP_SYNC_MODE");
    int32_t requestedBoundary = -1;

    if (env != NULL) {
      if (strcasecmp(env, "boundary") == 0 || strcasecmp(env, "sparse") == 0 ||
          strcasecmp(env, "gpu") == 0 || strcasecmp(env, "on") == 0 ||
          strcasecmp(env, "yes") == 0 || strcasecmp(env, "true") == 0 ||
          strcmp(env, "1") == 0) {
        requestedBoundary = 1;
      } else if (strcasecmp(env, "full") == 0 || strcasecmp(env, "cpu") == 0 ||
                 strcasecmp(env, "off") == 0 || strcasecmp(env, "no") == 0 ||
                 strcasecmp(env, "false") == 0 || strcmp(env, "0") == 0) {
        requestedBoundary = 0;
      }
    }

    if (requestedBoundary < 0) {
      requestedBoundary = displacement_boundary_sync_safe();
    }

    g_dispSyncMode =
        (requestedBoundary != 0 && displacement_boundary_sync_safe()) ? 1 : 0;

    if ((Global.myID == 0) && (g_dispSyncReported == 0)) {
      monitor_print(
          "Displacement sync path: %s%s\n",
          g_dispSyncMode ? "boundary" : "full",
          (requestedBoundary != 0 && g_dispSyncMode == 0)
              ? " (boundary requested but unsafe for this run)"
              : "");
      g_dispSyncReported = 1;
    }
  }

  return g_dispSyncMode;
}

static void disp_boundary_append_messenger_nodes(int32_t* nodes,
                                                 int32_t* count,
                                                 messenger_t* messenger) {
  while (messenger != NULL) {
    for (int32_t i = 0; i < messenger->nidx; i++) {
      nodes[(*count)++] = messenger->mapping[i];
    }
    messenger = messenger->next;
  }
}

static void solver_init_disp_boundary_sync(mysolver_t* solver) {
  int32_t capacity;
  int32_t count = 0;
  int32_t uniqueCount = 0;
  int32_t* nodes = NULL;
  cudaError_t cerror;

  if (g_dispBoundaryReady != 0 || displacement_sync_mode() == 0) {
    return;
  }

  compute_adjust_cache_build();

  capacity =
      messenger_countnodes(solver->an_sched->first_s) +
      messenger_countnodes(solver->an_sched->first_c) +
      messenger_countnodes(solver->dn_sched->first_s) +
      messenger_countnodes(solver->dn_sched->first_c) + g_adjustDnCount +
      g_adjustParentCount;

  if (capacity <= 0) {
    g_dispBoundaryReady = 1;
    g_dispBoundaryNodeCount = 0;
    return;
  }

  nodes = (int32_t*)malloc(sizeof(int32_t) * capacity);
  if (nodes == NULL) {
    fprintf(stderr, "Thread %d: displacement boundary list alloc failed\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  disp_boundary_append_messenger_nodes(nodes, &count, solver->an_sched->first_s);
  disp_boundary_append_messenger_nodes(nodes, &count, solver->an_sched->first_c);
  disp_boundary_append_messenger_nodes(nodes, &count, solver->dn_sched->first_s);
  disp_boundary_append_messenger_nodes(nodes, &count, solver->dn_sched->first_c);

  for (int32_t i = 0; i < g_adjustDnCount; i++) {
    nodes[count++] = g_adjustDanglingLnid[i];
  }
  for (int32_t i = 0; i < g_adjustParentCount; i++) {
    nodes[count++] = g_adjustParentIds[i];
  }

  qsort(nodes, (size_t)count, sizeof(int32_t), compare_int32_asc);
  for (int32_t i = 0; i < count; i++) {
    if ((i == 0) || (nodes[i] != nodes[i - 1])) {
      nodes[uniqueCount++] = nodes[i];
    }
  }

  g_dispBoundaryNodeCount = uniqueCount;
  if (g_dispBoundaryNodeCount <= 0) {
    free(nodes);
    g_dispBoundaryReady = 1;
    return;
  }

  g_dispBoundaryNodesHost =
      (int32_t*)malloc(sizeof(int32_t) * g_dispBoundaryNodeCount);
  if (g_dispBoundaryNodesHost == NULL) {
    free(nodes);
    fprintf(stderr, "Thread %d: displacement boundary unique alloc failed\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  memcpy(g_dispBoundaryNodesHost, nodes,
         sizeof(int32_t) * g_dispBoundaryNodeCount);
  free(nodes);

  cerror = cudaMalloc((void**)&g_dispBoundaryNodesDevice,
                      sizeof(int32_t) * g_dispBoundaryNodeCount);
  if (cerror != cudaSuccess) {
    fprintf(stderr, "Thread %d: displacement boundary node device alloc - %s\n",
            Global.myID, cudaGetErrorString(cerror));
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  cerror = cudaMalloc((void**)&g_dispBoundaryBufferDevice,
                      sizeof(fvector_t) * g_dispBoundaryNodeCount);
  if (cerror != cudaSuccess) {
    fprintf(stderr,
            "Thread %d: displacement boundary buffer device alloc - %s\n",
            Global.myID, cudaGetErrorString(cerror));
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  cerror = cudaHostAlloc((void**)&g_dispBoundaryBufferHost,
                         sizeof(fvector_t) * g_dispBoundaryNodeCount,
                         cudaHostAllocDefault);
  if (cerror != cudaSuccess) {
    g_dispBoundaryBufferHost =
        (fvector_t*)malloc(sizeof(fvector_t) * g_dispBoundaryNodeCount);
    if (g_dispBoundaryBufferHost == NULL) {
      fprintf(stderr,
              "Thread %d: displacement boundary buffer host alloc failed\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
  }

  cerror = cudaMemcpy(g_dispBoundaryNodesDevice, g_dispBoundaryNodesHost,
                      sizeof(int32_t) * g_dispBoundaryNodeCount,
                      cudaMemcpyHostToDevice);
  if (cerror != cudaSuccess) {
    fprintf(stderr, "Thread %d: displacement boundary node copy - %s\n",
            Global.myID, cudaGetErrorString(cerror));
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  if (Global.myID == 0) {
    monitor_print("Displacement boundary sync nodes: %d of %d harbored\n",
                  g_dispBoundaryNodeCount, Global.myMesh->nharbored);
  }

  g_dispBoundaryReady = 1;
}

static void solver_load_disp_gpu(mysolver_t* solver, mesh_t* mesh) {
  if (displacement_sync_mode() != 0) {
    return;
  }

  /* Synchronize for beginning of timestep */
  cudaStreamSynchronize(solver->streams[CUDA_STREAM_MAIN]);

  /* Copy current tm1 displacement only.
   * tm2 is maintained on device by the end-of-step device pointer swap in
   * solver_reset_gpu(), matching legacy GPU behavior and reducing PCIe/NVLink
   * traffic by one full displacement vector per step. */
  cudaMemcpyAsync(solver->gpuData->tm1Device, solver->tm1,
                  solver->gpuData->nharbored * sizeof(fvector_t),
                  cudaMemcpyHostToDevice, solver->streams[CUDA_STREAM_MAIN]);
  solver->gpu_spec->numbytespci +=
      solver->gpuData->nharbored * sizeof(fvector_t);

  /* Uncomment for timing tests */
  // cudaStreamSynchronize(solver->streams[CUDA_STREAM_MAIN]);

  return;
}

static void solver_load_forces_gpu(mysolver_t* solver, mesh_t* mesh) {
  /* Copy forces to device */
  cudaMemcpyAsync(solver->gpuData->forceDevice, solver->force,
                  mesh->nharbored * sizeof(fvector_t), cudaMemcpyHostToDevice,
                  solver->streams[CUDA_STREAM_MAIN]);
  solver->gpu_spec->numbytespci += mesh->nharbored * sizeof(fvector_t);
  if (!use_gpu_topography_fastpath()) {
    solver_record_gpu_event(solver, CUDA_EVENT_FORCE_UPLOAD_DONE,
                            "force upload");
  }
  solver_record_force_upload_source_event(solver);

  /* Uncomment for timing tests */
  // cudaStreamSynchronize(solver->streams[CUDA_STREAM_MAIN]);

  return;
}

static void solver_prepare_force_upload_ranges(mesh_t* mesh,
                                               int32_t sourceCount) {
  const int32_t topoRangeCount = topo_get_node_range_count();
  const int32_t* topoRangeStarts = topo_get_node_range_starts();
  const int32_t* topoRangeLengths = topo_get_node_range_lengths();
  const int32_t totalIntervals = topoRangeCount + sourceCount;
  const int32_t mergeGapNodes = 256;
  int32_t* sortedSourceNodes = NULL;
  int32_t topoIdx = 0;
  int32_t sourceIdx = 0;
  int32_t outCount = 0;
  int64_t copiedNodes = 0;
  int32_t currentStart = 0;
  int32_t currentEnd = 0;
  int haveCurrent = 0;

  (void)mesh;

  if ((g_rangeUploadReady == 1) &&
      (g_rangeUploadSourceCountCached == sourceCount) &&
      (g_rangeUploadTopoCountCached == topoRangeCount)) {
    return;
  }

  if (totalIntervals > g_rangeUploadCapacity) {
    free(g_rangeUploadStartsHost);
    free(g_rangeUploadLengthsHost);
    g_rangeUploadStartsHost = NULL;
    g_rangeUploadLengthsHost = NULL;

    g_rangeUploadStartsHost =
        (int32_t*)malloc((size_t)totalIntervals * sizeof(int32_t));
    g_rangeUploadLengthsHost =
        (int32_t*)malloc((size_t)totalIntervals * sizeof(int32_t));
    if ((g_rangeUploadStartsHost == NULL) || (g_rangeUploadLengthsHost == NULL)) {
      fprintf(stderr, "Thread %d: range upload host alloc failed\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

    g_rangeUploadCapacity = totalIntervals;
  }

  if (sourceCount > 0) {
    sortedSourceNodes =
        (int32_t*)malloc((size_t)sourceCount * sizeof(int32_t));
    if (sortedSourceNodes == NULL) {
      fprintf(stderr, "Thread %d: range upload source sort alloc failed\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

    memcpy(sortedSourceNodes, Global.theNodesLoadedList,
           (size_t)sourceCount * sizeof(int32_t));
    qsort(sortedSourceNodes, (size_t)sourceCount, sizeof(int32_t),
          compare_int32_asc);
  }

  while ((topoIdx < topoRangeCount) || (sourceIdx < sourceCount)) {
    int32_t nextStart;
    int32_t nextEnd;

    if (topoIdx >= topoRangeCount) {
      nextStart = sortedSourceNodes[sourceIdx];
      nextEnd = nextStart + 1;
      sourceIdx++;
      while ((sourceIdx < sourceCount) &&
             (sortedSourceNodes[sourceIdx] == nextStart)) {
        sourceIdx++;
      }
    } else if (sourceIdx >= sourceCount) {
      nextStart = topoRangeStarts[topoIdx];
      nextEnd = nextStart + topoRangeLengths[topoIdx];
      topoIdx++;
    } else if (sortedSourceNodes[sourceIdx] < topoRangeStarts[topoIdx]) {
      nextStart = sortedSourceNodes[sourceIdx];
      nextEnd = nextStart + 1;
      sourceIdx++;
      while ((sourceIdx < sourceCount) &&
             (sortedSourceNodes[sourceIdx] == nextStart)) {
        sourceIdx++;
      }
    } else {
      nextStart = topoRangeStarts[topoIdx];
      nextEnd = nextStart + topoRangeLengths[topoIdx];
      topoIdx++;
    }

    if (nextEnd <= nextStart) {
      continue;
    }

    if (!haveCurrent) {
      currentStart = nextStart;
      currentEnd = nextEnd;
      haveCurrent = 1;
    } else if (nextStart <= (currentEnd + mergeGapNodes)) {
      if (nextEnd > currentEnd) {
        currentEnd = nextEnd;
      }
    } else {
      g_rangeUploadStartsHost[outCount] = currentStart;
      g_rangeUploadLengthsHost[outCount] = currentEnd - currentStart;
      copiedNodes += (int64_t)(currentEnd - currentStart);
      outCount++;
      currentStart = nextStart;
      currentEnd = nextEnd;
    }
  }

  if (haveCurrent) {
    g_rangeUploadStartsHost[outCount] = currentStart;
    g_rangeUploadLengthsHost[outCount] = currentEnd - currentStart;
    copiedNodes += (int64_t)(currentEnd - currentStart);
    outCount++;
  }

  free(sortedSourceNodes);

  g_rangeUploadCount = outCount;
  g_rangeUploadNodeCount = copiedNodes;
  g_rangeUploadReady = 1;
  g_rangeUploadSourceCountCached = sourceCount;
  g_rangeUploadTopoCountCached = topoRangeCount;
}

static void solver_load_forces_range_gpu(mysolver_t* solver, mesh_t* mesh,
                                         int32_t sourceCount) {
  int32_t i;

  solver_prepare_force_upload_ranges(mesh, sourceCount);

  Timer_Start("Sparse force upload memcpy");
  for (i = 0; i < g_rangeUploadCount; i++) {
    int32_t start = g_rangeUploadStartsHost[i];
    int32_t len = g_rangeUploadLengthsHost[i];
    cudaMemcpyAsync(solver->gpuData->forceDevice + start, solver->force + start,
                    (size_t)len * sizeof(fvector_t), cudaMemcpyHostToDevice,
                    solver->streams[CUDA_STREAM_MAIN]);
  }
  solver->gpu_spec->numbytespci +=
      g_rangeUploadNodeCount * (int64_t)sizeof(fvector_t);
  solver_record_gpu_event(solver, CUDA_EVENT_FORCE_UPLOAD_DONE,
                          "force upload");
  solver_record_force_upload_source_event(solver);
  Timer_Stop("Sparse force upload memcpy");
}

static void solver_load_forces_sparse_gpu(mysolver_t* solver, mesh_t* mesh) {
  int32_t sourceCount = Global.theNodesLoaded;
  int32_t topoCount = topo_get_nodes_count();
  const int32_t* topoNodes = topo_get_nodes_list();
  int32_t totalCount = sourceCount + topoCount;
  int32_t activeCount;
  int32_t useDenseUpload = 0;
  int32_t pathCode;
  int32_t uploadMode;
  int32_t idx;
  int blocksize;
  int gridsize;

  if (totalCount <= 0) {
    g_forceUploadSourceSafePending = 0;
    cudaMemsetAsync(solver->gpuData->forceDevice, 0,
                    mesh->nharbored * sizeof(fvector_t),
                    solver->streams[CUDA_STREAM_MAIN]);
    solver_record_gpu_event(solver, CUDA_EVENT_FORCE_UPLOAD_DONE,
                            "force upload");
    return;
  }

  uploadMode = sparse_force_upload_mode();
  if (uploadMode == 1) {
    if ((Global.myID == 0) && (g_sparseUploadPathReported != 1)) {
      monitor_print(
          "Sparse force upload path: dense (request=dense, source=%d, topo=%d, harbored=%d)\n",
          sourceCount, topoCount, mesh->nharbored);
      g_sparseUploadPathReported = 1;
    }
    Timer_Start("Sparse force upload memcpy");
    cudaMemcpyAsync(solver->gpuData->forceDevice, solver->force,
                    mesh->nharbored * sizeof(fvector_t),
                    cudaMemcpyHostToDevice,
                    solver->streams[CUDA_STREAM_MAIN]);
    solver->gpu_spec->numbytespci += mesh->nharbored * sizeof(fvector_t);
    solver_record_gpu_event(solver, CUDA_EVENT_FORCE_UPLOAD_DONE,
                            "force upload");
    solver_record_force_upload_source_event(solver);
    Timer_Stop("Sparse force upload memcpy");
    return;
  }

  if (uploadMode == 3) {
    solver_prepare_force_upload_ranges(mesh, sourceCount);
    if ((Global.myID == 0) && (g_sparseUploadPathReported != 3)) {
      monitor_print(
          "Sparse force upload path: range (request=range, spans=%d, copied=%lld, source=%d, topo=%d, harbored=%d)\n",
          g_rangeUploadCount, (long long)g_rangeUploadNodeCount, sourceCount,
          topoCount, mesh->nharbored);
      g_sparseUploadPathReported = 3;
    }
    solver_load_forces_range_gpu(solver, mesh, sourceCount);
    return;
  }

  if (totalCount > g_sparseNodeCapacity) {
    if (g_sparseNodeListDevice != NULL) {
      cudaFree(g_sparseNodeListDevice);
      g_sparseNodeListDevice = NULL;
    }
    if (g_sparseForcesDevice != NULL) {
      cudaFree(g_sparseForcesDevice);
      g_sparseForcesDevice = NULL;
    }
    if (g_sparseNodeListHost != NULL) {
      if (g_sparseHostPinned) {
        cudaFreeHost(g_sparseNodeListHost);
      } else {
        free(g_sparseNodeListHost);
      }
      g_sparseNodeListHost = NULL;
    }
    if (g_sparseForcesHost != NULL) {
      if (g_sparseHostPinned) {
        cudaFreeHost(g_sparseForcesHost);
      } else {
        free(g_sparseForcesHost);
      }
      g_sparseForcesHost = NULL;
    }
    g_sparseHostPinned = 0;

    if (cudaMalloc((void**)&g_sparseNodeListDevice,
                   totalCount * sizeof(int32_t)) != cudaSuccess) {
      fprintf(stderr, "Thread %d: sparse node list device alloc failed\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    if (cudaMalloc((void**)&g_sparseForcesDevice,
                   totalCount * sizeof(fvector_t)) != cudaSuccess) {
      fprintf(stderr, "Thread %d: sparse force device alloc failed\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

    if (cudaHostAlloc((void**)&g_sparseNodeListHost,
                      totalCount * sizeof(int32_t), cudaHostAllocDefault) ==
            cudaSuccess &&
        cudaHostAlloc((void**)&g_sparseForcesHost,
                      totalCount * sizeof(fvector_t),
                      cudaHostAllocDefault) == cudaSuccess) {
      g_sparseHostPinned = 1;
    } else {
      if (g_sparseNodeListHost != NULL) {
        cudaFreeHost(g_sparseNodeListHost);
        g_sparseNodeListHost = NULL;
      }
      if (g_sparseForcesHost != NULL) {
        cudaFreeHost(g_sparseForcesHost);
        g_sparseForcesHost = NULL;
      }

      g_sparseNodeListHost = (int32_t*)malloc(totalCount * sizeof(int32_t));
      g_sparseForcesHost = (fvector_t*)malloc(totalCount * sizeof(fvector_t));
      g_sparseHostPinned = 0;
    }

    if ((g_sparseNodeListHost == NULL) || (g_sparseForcesHost == NULL)) {
      fprintf(stderr, "Thread %d: sparse force host alloc failed\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

    g_sparseNodeCapacity = totalCount;
    g_sparseNodeListReady = 0;
    g_sparseSourceCountCached = -1;
    g_sparseTopoCountCached = -1;
    g_sparseCopiedSourceCount = -1;
    g_sparseCopiedTopoCount = -1;
  }

  Timer_Start("Sparse force upload build");
  if ((!g_sparseNodeListReady) || (g_sparseSourceCountCached != sourceCount) ||
      (g_sparseTopoCountCached != topoCount)) {
    int32_t writeIdx = 0;
    int32_t sourceIdx = 0;
    int32_t topoIdx = 0;
    int32_t haveLast = 0;
    int32_t lastNode = 0;
    int32_t* sortedSourceNodes = NULL;

    if (sourceCount > 0) {
      sortedSourceNodes = (int32_t*)malloc((size_t)sourceCount * sizeof(int32_t));
      if (sortedSourceNodes == NULL) {
        fprintf(stderr, "Thread %d: sparse source list scratch alloc failed\n",
                Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }

      memcpy(sortedSourceNodes, Global.theNodesLoadedList,
             (size_t)sourceCount * sizeof(int32_t));
      qsort(sortedSourceNodes, (size_t)sourceCount, sizeof(int32_t),
            compare_int32_asc);
    }

    while ((sourceIdx < sourceCount) || (topoIdx < topoCount)) {
      int32_t candidate;

      if (topoIdx >= topoCount) {
        candidate = sortedSourceNodes[sourceIdx++];
      } else if (sourceIdx >= sourceCount) {
        candidate = topoNodes[topoIdx++];
      } else {
        int32_t sourceNode = sortedSourceNodes[sourceIdx];
        int32_t topoNode = topoNodes[topoIdx];

        if (sourceNode < topoNode) {
          candidate = sourceNode;
          sourceIdx++;
        } else if (topoNode < sourceNode) {
          candidate = topoNode;
          topoIdx++;
        } else {
          candidate = sourceNode;
          sourceIdx++;
          topoIdx++;
        }
      }

      if ((!haveLast) || (candidate != lastNode)) {
        g_sparseNodeListHost[writeIdx++] = candidate;
        lastNode = candidate;
        haveLast = 1;
      }

      while ((sourceIdx < sourceCount) &&
             (sortedSourceNodes[sourceIdx] == lastNode)) {
        sourceIdx++;
      }
      while ((topoIdx < topoCount) && (topoNodes[topoIdx] == lastNode)) {
        topoIdx++;
      }
    }

    if (sortedSourceNodes != NULL) {
      free(sortedSourceNodes);
      sortedSourceNodes = NULL;
    }

    g_sparseActiveCount = writeIdx;
    g_sparseNodeListReady = 1;
    g_sparseSourceCountCached = sourceCount;
    g_sparseTopoCountCached = topoCount;
    g_sparseCopiedSourceCount = -1;
    g_sparseCopiedTopoCount = -1;
  }

  activeCount = g_sparseActiveCount;
  if (uploadMode == 1) {
    useDenseUpload = 1;
  } else if (uploadMode == 2) {
    useDenseUpload = 0;
  } else if ((mesh->nharbored > 0) && ((activeCount * 4) >= mesh->nharbored)) {
    useDenseUpload = 1;
  }

  pathCode = useDenseUpload ? 1 : 2;
  if ((Global.myID == 0) && (g_sparseUploadPathReported != pathCode)) {
    const char* requestName =
        (uploadMode == 1) ? "dense" :
        (uploadMode == 2) ? "sparse" :
        (uploadMode == 3) ? "range" : "auto";
    const char* chosenName = useDenseUpload ? "dense" : "sparse";
    monitor_print(
        "Sparse force upload path: %s (request=%s, active=%d, source=%d, topo=%d, harbored=%d)\n",
        chosenName, requestName, activeCount, sourceCount, topoCount,
        mesh->nharbored);
    g_sparseUploadPathReported = pathCode;
  }

  if (useDenseUpload) {
    Timer_Stop("Sparse force upload build");
    Timer_Start("Sparse force upload memcpy");
    cudaMemcpyAsync(solver->gpuData->forceDevice, solver->force,
                    mesh->nharbored * sizeof(fvector_t),
                    cudaMemcpyHostToDevice,
                    solver->streams[CUDA_STREAM_MAIN]);
    solver->gpu_spec->numbytespci += mesh->nharbored * sizeof(fvector_t);
    solver_record_gpu_event(solver, CUDA_EVENT_FORCE_UPLOAD_DONE,
                            "force upload");
    solver_record_force_upload_source_event(solver);
    Timer_Stop("Sparse force upload memcpy");
    return;
  }

  g_forceUploadSourceSafePending = 0;

  cudaMemsetAsync(solver->gpuData->forceDevice, 0,
                  mesh->nharbored * sizeof(fvector_t),
                  solver->streams[CUDA_STREAM_MAIN]);

  for (idx = 0; idx < activeCount; idx++) {
    int32_t lnid = g_sparseNodeListHost[idx];
    g_sparseForcesHost[idx] = solver->force[lnid];
  }
  Timer_Stop("Sparse force upload build");

  Timer_Start("Sparse force upload memcpy");
  if ((g_sparseSourceCountCached == sourceCount) &&
      (g_sparseTopoCountCached == topoCount) && g_sparseNodeListReady) {
    if ((g_sparseCopiedSourceCount != sourceCount) ||
        (g_sparseCopiedTopoCount != topoCount)) {
      cudaMemcpyAsync(g_sparseNodeListDevice, g_sparseNodeListHost,
                      activeCount * sizeof(int32_t), cudaMemcpyHostToDevice,
                      solver->streams[CUDA_STREAM_MAIN]);
      solver->gpu_spec->numbytespci += activeCount * sizeof(int32_t);
      g_sparseCopiedSourceCount = sourceCount;
      g_sparseCopiedTopoCount = topoCount;
    }
  }

  cudaMemcpyAsync(g_sparseForcesDevice, g_sparseForcesHost,
                  activeCount * sizeof(fvector_t), cudaMemcpyHostToDevice,
                  solver->streams[CUDA_STREAM_MAIN]);
  solver->gpu_spec->numbytespci += activeCount * sizeof(fvector_t);
  Timer_Stop("Sparse force upload memcpy");

  Timer_Start("Sparse force upload scatter");
  blocksize = 256;
  gridsize = (activeCount + blocksize - 1) / blocksize;
  kernelSparseForceScatter<<<gridsize, blocksize, 0,
                             solver->streams[CUDA_STREAM_MAIN]>>>(
      activeCount, g_sparseNodeListDevice, g_sparseForcesDevice,
      solver->gpuData->forceDevice);

  cudaError_t cerror = cudaGetLastError();
  if (cerror != cudaSuccess) {
    fprintf(stderr, "Thread %d: sparse force scatter kernel - %s\n",
            Global.myID, cudaGetErrorString(cerror));
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  solver_record_gpu_event(solver, CUDA_EVENT_FORCE_UPLOAD_DONE,
                          "force upload");
  Timer_Stop("Sparse force upload scatter");
}

static void solver_clear_sparse_host_forces(mysolver_t* solver) {
  int32_t clearMode = solver_sparse_clear_mode();

  Timer_Start("Disp force upload safe wait");
  solver_wait_force_upload_source_event();
  Timer_Stop("Disp force upload safe wait");

  Timer_Start("Disp host force clear only");
  if (clearMode == 3) {
    int32_t fullCount = Global.myMesh->nharbored;
#ifdef _OPENMP
    if (full_clear_use_omp()) {
      int32_t threadCount = omp_get_max_threads();
      int32_t minNodesPerThread = 262144;

      if ((threadCount > 1) &&
          ((int64_t)fullCount >=
           (int64_t)threadCount * (int64_t)minNodesPerThread)) {
#pragma omp parallel for schedule(static)
        for (int32_t i = 0; i < fullCount; i++) {
          solver->force[i].f[0] = 0;
          solver->force[i].f[1] = 0;
          solver->force[i].f[2] = 0;
        }
      } else {
        memset(solver->force, 0, sizeof(fvector_t) * fullCount);
      }
    } else
#endif
    {
      memset(solver->force, 0, sizeof(fvector_t) * fullCount);
    }
  } else if (clearMode == 1) {
    int32_t topoCount = topo_get_nodes_count();
    const int32_t* topoNodes = topo_get_nodes_list();

    for (int32_t i = 0; i < topoCount; i++) {
      int32_t lnid = topoNodes[i];
      solver->force[lnid].f[0] = 0;
      solver->force[lnid].f[1] = 0;
      solver->force[lnid].f[2] = 0;
    }
  } else {
    int32_t rangeCount = topo_get_node_range_count();
    const int32_t* rangeStarts = topo_get_node_range_starts();
    const int32_t* rangeLengths = topo_get_node_range_lengths();

    for (int32_t i = 0; i < rangeCount; i++) {
      int32_t start = rangeStarts[i];
      int32_t len = rangeLengths[i];
      memset(solver->force + start, 0, (size_t)len * sizeof(fvector_t));
    }
  }
  Timer_Stop("Disp host force clear only");
}

static int32_t solver_sparse_clear_mode(void) {
  if (g_sparseClearMode < 0) {
    const char* env = getenv("HERCULES_SPARSE_CLEAR_MODE");
    int32_t requestedMode = 3;
    int32_t harboredCount = Global.myMesh->nharbored;
    int32_t topoCount = topo_get_nodes_count();
    int32_t rangeCount = topo_get_node_range_count();
    int64_t totalRangeSpan = 0;
    const int32_t* rangeLengths = topo_get_node_range_lengths();

    /* 0: range memset, 1: node-list clear, 2: auto-select, 3: full memset */
    if (env != NULL) {
      if (strcasecmp(env, "nodes") == 0 || strcmp(env, "1") == 0) {
        requestedMode = 1;
      } else if (strcasecmp(env, "full") == 0 ||
                 strcasecmp(env, "dense") == 0 || strcmp(env, "3") == 0) {
        requestedMode = 3;
      } else if (strcasecmp(env, "range") == 0 ||
                 strcasecmp(env, "ranges") == 0 || strcmp(env, "0") == 0) {
        requestedMode = 0;
      } else if (strcasecmp(env, "auto") == 0 || strcmp(env, "2") == 0) {
        requestedMode = 2;
      }
    }

    for (int32_t i = 0; i < rangeCount; i++) {
      totalRangeSpan += (int64_t)rangeLengths[i];
    }

    if (requestedMode == 2) {
      int32_t useFull = 0;
      int32_t useRanges = 0;

      /* Auto policy: prefer full clear when topography touches a meaningful
       * fraction of local nodes, otherwise use ranges only when the covered
       * span is compact and the range list is small. */
      if ((harboredCount > 0) && ((topoCount * 4) >= harboredCount)) {
        useFull = 1;
      } else if (topoCount > 0 &&
                 totalRangeSpan <= ((int64_t)topoCount * 12) / 10 &&
                 rangeCount <= 128) {
        useRanges = 1;
      }

      if (useFull) {
        g_sparseClearMode = 3;
      } else {
        g_sparseClearMode = useRanges ? 0 : 1;
      }
    } else {
      g_sparseClearMode = requestedMode;
    }

    if (Global.myID == 0 && g_sparseClearModeReported == 0) {
      const char* modeName =
          (g_sparseClearMode == 0) ? "range" :
          (g_sparseClearMode == 1) ? "nodes" : "full";
      const char* requestName =
          (requestedMode == 2) ? "auto" :
          (requestedMode == 0) ? "range" :
          (requestedMode == 1) ? "nodes" : "full";
      monitor_print("Sparse host force reset mode: %s (request=%s, topo=%d, ranges=%d, span=%lld, harbored=%d)\n",
                    modeName, requestName, topoCount, rangeCount,
                    (long long)totalRangeSpan, harboredCount);
      g_sparseClearModeReported = 1;
    }
  }

  return g_sparseClearMode;
}

static void solver_unload_forces_gpu(mysolver_t* solver, mesh_t* mesh) {
  /* Only copy forces back if Stiffness-Effective or BKT damping */
  if (((Param.theTypeOfDamping != BKT) && (Param.theStiffness == EFFECTIVE)) ||
      (Param.theTypeOfDamping == BKT) || acoustic_pressure_enabled()) {
    /* Copy forces back to host */
    cudaMemcpyAsync(solver->force, solver->gpuData->forceDevice,
                    mesh->nharbored * sizeof(fvector_t), cudaMemcpyDeviceToHost,
                    solver->streams[CUDA_STREAM_MAIN]);
    solver->gpu_spec->numbytespci += mesh->nharbored * sizeof(fvector_t);
  }

  /* Uncomment for timing tests */
  // cudaStreamSynchronize(solver->streams[CUDA_STREAM_MAIN]);

  return;
}

static void solver_record_gpu_event(mysolver_t* solver, int eventID,
                                    const char* eventLabel) {
  cudaError_t cerror =
      cudaEventRecord(solver->events[eventID], solver->streams[CUDA_STREAM_MAIN]);
  if (cerror != cudaSuccess) {
    fprintf(stderr, "Thread %d: record %s event - %s\n", Global.myID,
            eventLabel, cudaGetErrorString(cerror));
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
}

static void solver_record_force_kernel_done_event(mysolver_t* solver) {
  solver_record_gpu_event(solver, CUDA_EVENT_FORCE_KERNEL_DONE,
                          "force kernel");
}

static void solver_unload_disp_gpu(mysolver_t* solver, mesh_t* mesh) {
  if (displacement_sync_mode() != 0) {
    int blocksize;
    int gridsize;
    cudaError_t cerror;

    solver_init_disp_boundary_sync(solver);
    if (g_dispBoundaryNodeCount <= 0) {
      return;
    }

    blocksize = 256;
    gridsize = (g_dispBoundaryNodeCount + blocksize - 1) / blocksize;
    kernelDispBoundaryGather<<<gridsize, blocksize, 0,
                               solver->streams[CUDA_STREAM_MAIN]>>>(
        g_dispBoundaryNodeCount, g_dispBoundaryNodesDevice,
        solver->gpuData->tm2Device, g_dispBoundaryBufferDevice);

    cerror = cudaGetLastError();
    if (cerror != cudaSuccess) {
      fprintf(stderr, "Thread %d: displacement boundary gather - %s\n",
              Global.myID, cudaGetErrorString(cerror));
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

    cudaMemcpyAsync(g_dispBoundaryBufferHost, g_dispBoundaryBufferDevice,
                    sizeof(fvector_t) * g_dispBoundaryNodeCount,
                    cudaMemcpyDeviceToHost,
                    solver->streams[CUDA_STREAM_MAIN]);
    solver->gpu_spec->numbytespci +=
        g_dispBoundaryNodeCount * (int64_t)sizeof(fvector_t);
    return;
  }

  /* Copy updated tm2 displacements from device */
  cudaMemcpyAsync(solver->tm2, solver->gpuData->tm2Device,
                  mesh->nharbored * sizeof(fvector_t), cudaMemcpyDeviceToHost,
                  solver->streams[CUDA_STREAM_MAIN]);

  solver->gpu_spec->numbytespci += mesh->nharbored * sizeof(fvector_t);

  /* Uncomment for timing tests */
  // cudaStreamSynchronize(solver->streams[CUDA_STREAM_MAIN]);

  return;
}

static void solver_finalize_disp_boundary_download(mysolver_t* solver) {
  (void)solver;

  if (displacement_sync_mode() == 0 || g_dispBoundaryNodeCount <= 0) {
    return;
  }

  Timer_Start("Disp boundary host scatter");
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(g_dispBoundaryNodeCount > 100000)
#endif
  for (int32_t i = 0; i < g_dispBoundaryNodeCount; i++) {
    solver->tm2[g_dispBoundaryNodesHost[i]] = g_dispBoundaryBufferHost[i];
  }
  Timer_Stop("Disp boundary host scatter");
}

static void solver_upload_disp_boundary_gpu(mysolver_t* solver) {
  int blocksize;
  int gridsize;
  cudaError_t cerror;

  if (displacement_sync_mode() == 0) {
    return;
  }

  solver_init_disp_boundary_sync(solver);
  if (g_dispBoundaryNodeCount <= 0) {
    return;
  }

  Timer_Start("GPU Patch Disp");
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(g_dispBoundaryNodeCount > 100000)
#endif
  for (int32_t i = 0; i < g_dispBoundaryNodeCount; i++) {
    g_dispBoundaryBufferHost[i] = solver->tm2[g_dispBoundaryNodesHost[i]];
  }

  cudaMemcpyAsync(g_dispBoundaryBufferDevice, g_dispBoundaryBufferHost,
                  sizeof(fvector_t) * g_dispBoundaryNodeCount,
                  cudaMemcpyHostToDevice, solver->streams[CUDA_STREAM_MAIN]);
  solver->gpu_spec->numbytespci +=
      g_dispBoundaryNodeCount * (int64_t)sizeof(fvector_t);

  blocksize = 256;
  gridsize = (g_dispBoundaryNodeCount + blocksize - 1) / blocksize;
  kernelDispBoundaryScatter<<<gridsize, blocksize, 0,
                              solver->streams[CUDA_STREAM_MAIN]>>>(
      g_dispBoundaryNodeCount, g_dispBoundaryNodesDevice,
      g_dispBoundaryBufferDevice, solver->gpuData->tm1Device);

  cerror = cudaGetLastError();
  if (cerror != cudaSuccess) {
    fprintf(stderr, "Thread %d: displacement boundary scatter - %s\n",
            Global.myID, cudaGetErrorString(cerror));
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }
  Timer_Stop("GPU Patch Disp");
}

static void solver_wait_physics_gpu(mysolver_t* solver, mesh_t* mesh) {
  (void)mesh;

  /* Snapshot previous-step station displacement for accelerations. */
  if (Param.printStationAccelerations == YES) {
    Timer_Start("Station Accel History Copy");
    solver_snapshot_station_accel_history(solver);
    Timer_Stop("Station Accel History Copy");
  }

  /* Wait for the appropriate physics calculation to complete */
  if (((Param.theTypeOfDamping != BKT) && (Param.theStiffness == EFFECTIVE)) ||
      (Param.theTypeOfDamping == BKT) || acoustic_pressure_enabled()) {
    Timer_Start("Force Stream Sync");

    Timer_Start("Force kernel wait only");

    Timer_Start("Force stiffness wait only");
    Timer_Start("Force pre-stiffness wait only");
    cudaEventSynchronize(solver->events[CUDA_EVENT_FORCE_UPLOAD_DONE]);
    Timer_Stop("Force pre-stiffness wait only");

    Timer_Start("Force stiffness core wait only");
    cudaEventSynchronize(solver->events[CUDA_EVENT_FORCE_STIFFNESS_DONE]);
    Timer_Stop("Force stiffness core wait only");
    Timer_Stop("Force stiffness wait only");

    Timer_Start("Force damping conv wait only");
    cudaEventSynchronize(solver->events[CUDA_EVENT_FORCE_DAMPING_CONV_DONE]);
    Timer_Stop("Force damping conv wait only");

    Timer_Start("Force damping local wait only");
    cudaEventSynchronize(solver->events[CUDA_EVENT_FORCE_DAMPING_LOCAL_DONE]);
    Timer_Stop("Force damping local wait only");

    cudaEventSynchronize(solver->events[CUDA_EVENT_FORCE_KERNEL_DONE]);
    Timer_Stop("Force kernel wait only");

    Timer_Start("Force postcopy sync wait");
    cudaStreamSynchronize(solver->streams[CUDA_STREAM_MAIN]);
    Timer_Stop("Force postcopy sync wait");

    Timer_Stop("Force Stream Sync");
  }

  return;
}

static void solver_wait_disp_gpu(mysolver_t* solver, mesh_t* mesh) {
  Timer_Start("Compute new displacement");

  Timer_Start("Disp kernel wait only");
  cudaEventSynchronize(solver->events[CUDA_EVENT_DISP_KERNEL_DONE]);
  Timer_Stop("Disp kernel wait only");

  Timer_Start("Disp postcopy sync wait");
  Timer_Start("Disp kernel sync wait");
  cudaStreamSynchronize(solver->streams[CUDA_STREAM_MAIN]);
  Timer_Stop("Disp kernel sync wait");
  Timer_Stop("Disp postcopy sync wait");

  solver_finalize_disp_boundary_download(solver);

  Timer_Stop("Compute new displacement");
  return;
}

/** Compute new displacements of my harbored nodes */
static void solver_compute_displacement_gpu(int32_t myID, mysolver_t* solver,
                                            mesh_t* mesh) {
  Timer_Start("Compute new displacement");

  /* Ensure host forces have been copied to device */
  // cudaStreamSynchronize(solver->streams[CUDA_STREAM_MAIN]);

  int blocksize = solver->gpu_kernel[CUDA_KERNEL_DISPLACEMENT].blocksize;
  int gridsize = solver->gpu_kernel[CUDA_KERNEL_DISPLACEMENT].gridsize;
  int sharedmem = solver->gpu_kernel[CUDA_KERNEL_DISPLACEMENT].sharedmem;

  cudaGetLastError();

  Timer_Start("Disp kernel enqueue");
  kernelDispCalc<<<gridsize, blocksize, sharedmem,
                   solver->streams[CUDA_STREAM_MAIN]>>>(
      mesh->nharbored, solver->gpuDataDevice, Param.printStationAccelerations);

  cudaError_t cerror = cudaGetLastError();
  if (cerror != cudaSuccess) {
    fprintf(stderr, "Thread %d: Calc damping conv kernel - %s\n", myID,
            cudaGetErrorString(cerror));
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  cudaEventRecord(solver->events[CUDA_EVENT_DISP_KERNEL_DONE],
                  solver->streams[CUDA_STREAM_MAIN]);
  Timer_Stop("Disp kernel enqueue");

  solver->gpu_spec->numflops +=
      kernel_flops_per_thread(FLOP_DISPLACEMENT) * mesh->nharbored;

  solver->gpu_spec->numbytes +=
      kernel_mem_per_thread(FLOP_DISPLACEMENT) * mesh->nharbored;

  Timer_Start("Disp host force reset");
  if (use_sparse_force_upload_topo()) {
    solver_clear_sparse_host_forces(solver);
  } else if (need_host_force_upload()) {
    /* zero out the force vector for all nodes */
    memset(solver->force, 0, sizeof(fvector_t) * mesh->nharbored);
  }
  Timer_Stop("Disp host force reset");

  /* Uncomment for timing tests */
  // cudaStreamSynchronize(solver->streams[CUDA_STREAM_MAIN]);

  Timer_Stop("Compute new displacement");
}

static void solver_geostatic_fix(int step) {
  if (Param.includeNonlinearAnalysis == YES) {
    Timer_Start("Compute addforces gravity");
    if (get_geostatic_total_time() > 0) {
      geostatic_displacements_fix(Global.myMesh, Global.mySolver,
                                  Param.theDomainZ, Param.theDeltaT, step);
    }
    Timer_Stop("Compute addforces gravity");
  }
}

/** Share the displacement of anchored nodes with other processors. */
static void solver_send_displacement_anchored(mysolver_t* solver) {
  //    HU_COND_GLOBAL_BARRIER( Param.theTimingBarriersFlag );

  Timer_Start("3rd schedule send data (sharing)");
  schedule_senddata(Global.mySolver->an_sched, Global.mySolver->tm2,
                    sizeof(fvector_t) / sizeof(solver_float), SHARING,
                    AN_DISP_MSG);
  Timer_Stop("3rd schedule send data (sharing)");
}

/** Adjust the displacements of my LOCAL dangling nodes. */
static void solver_adjust_displacement(mysolver_t* solver) {
  Timer_Start("2nd compute adjust (assignment)");
  compute_adjust(Global.mySolver->tm2, sizeof(fvector_t) / sizeof(solver_float),
                 ASSIGNMENT);
  Timer_Stop("2nd compute adjust (assignment)");
}

/** Share the displacement of dangling nodes with other processors. */
static void solver_send_displacement_dangling(mysolver_t* solver) {
  //    HU_COND_GLOBAL_BARRIER( Param.theTimingBarriersFlag );

  Timer_Start("4th schedule send data (sharing)");
  schedule_senddata(Global.mySolver->dn_sched, Global.mySolver->tm2,
                    sizeof(fvector_t) / sizeof(solver_float), SHARING,
                    DN_DISP_MSG);
  Timer_Stop("4th schedule send data (sharing)");
}

/**
 * Hook for instrumentation or other functionality in the solver_run loop.
 */
static void solver_loop_hook_bottom(mysolver_t* solver, mesh_t* mesh,
                                    int step) {
  if (0) {
    solver_debug_overflow(solver, mesh, step);
  }
}

static void solver_output_wavefield_close(void) {
  if (DO_OUTPUT && (Param.FourDOutFp != NULL)) {
    fclose(Param.FourDOutFp); /* close the output file */
  }
}

static void solver_run_collect_timers(void) {
  /* Get min and max of individual cumulative times.
   * this is NOT best and worst single step times
   */
  if (Timer_Exists("Print Planes")) {
    Timer_Reduce("Print Planes", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }

  if (Timer_Exists("Print Stations")) {
    Timer_Reduce("Print Stations", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }

  if (Timer_Exists("Compute Non-linear Entities")) {
    Timer_Reduce("Compute Non-linear Entities",
                 (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }

  if (Timer_Exists("Compute addforces Non-linear")) {
    Timer_Reduce("Compute addforces Non-linear",
                 (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }

  if (Timer_Exists("Compute addforces gravity")) {
    Timer_Reduce("Compute addforces gravity", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }

  if (Timer_Exists("Compute addforces Topography")) {
    Timer_Reduce("Compute addforces Topography",
                 (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }

  if (Timer_Exists("Compute addforces Incident Plane Waves")) {
    Timer_Reduce("Compute addforces Incident Plane Waves",
                 (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }

  if (Timer_Exists("Solver drm output")) {
    Timer_Reduce("Solver drm output", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }

  if (Timer_Exists("Solver drm read displacements")) {
    Timer_Reduce("Solver drm read displacements",
                 (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }

  if (Timer_Exists("Solver drm force compute")) {
    Timer_Reduce("Solver drm force compute", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Sparse force upload build")) {
    Timer_Reduce("Sparse force upload build", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Sparse force upload memcpy")) {
    Timer_Reduce("Sparse force upload memcpy", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Sparse force upload scatter")) {
    Timer_Reduce("Sparse force upload scatter",
                 (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }
  if (Timer_Exists("Load Fixedbase Disps")) {
    Timer_Reduce("Load Fixedbase Disps", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("GPU Memory Copy")) {
    Timer_Reduce("GPU Memory Copy", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("GPU Load Disp")) {
    Timer_Reduce("GPU Load Disp", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("GPU Unload Forces")) {
    Timer_Reduce("GPU Unload Forces", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("GPU Load Forces")) {
    Timer_Reduce("GPU Load Forces", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("GPU Unload Disp")) {
    Timer_Reduce("GPU Unload Disp", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("GPU Patch Disp")) {
    Timer_Reduce("GPU Patch Disp", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("GPU Reset")) {
    Timer_Reduce("GPU Reset", (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }
  if (Timer_Exists("Physics Force Prep")) {
    Timer_Reduce("Physics Force Prep", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Physics Force Kernels Enqueue")) {
    Timer_Reduce("Physics Force Kernels Enqueue",
                 (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }
  if (Timer_Exists("Physics Force Wait/Post")) {
    Timer_Reduce("Physics Force Wait/Post",
                 (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }
  if (Timer_Exists("Station Accel History Copy")) {
    Timer_Reduce("Station Accel History Copy",
                 (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }
  if (Timer_Exists("Force Stream Sync")) {
    Timer_Reduce("Force Stream Sync", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Force stiffness wait only")) {
    Timer_Reduce("Force stiffness wait only",
                 (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }
  if (Timer_Exists("Force pre-stiffness wait only")) {
    Timer_Reduce("Force pre-stiffness wait only",
                 (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }
  if (Timer_Exists("Force stiffness core wait only")) {
    Timer_Reduce("Force stiffness core wait only",
                 (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }
  if (Timer_Exists("Force damping conv wait only")) {
    Timer_Reduce("Force damping conv wait only",
                 (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }
  if (Timer_Exists("Force damping local wait only")) {
    Timer_Reduce("Force damping local wait only",
                 (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }
  if (Timer_Exists("Force kernel wait only")) {
    Timer_Reduce("Force kernel wait only",
                 (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }
  if (Timer_Exists("Acoustic pressure traction")) {
    Timer_Reduce("Acoustic pressure traction",
                 (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }
  if (Timer_Exists("Acoustic Pressure Step")) {
    Timer_Reduce("Acoustic Pressure Step",
                 (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }
  if (Timer_Exists("Force postcopy sync wait")) {
    Timer_Reduce("Force postcopy sync wait",
                 (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }
  if (Timer_Exists("Force Post CPU")) {
    Timer_Reduce("Force Post CPU", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Physics Displacement Enqueue")) {
    Timer_Reduce("Physics Displacement Enqueue",
                 (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }
  if (Timer_Exists("Physics Displacement Wait/Post")) {
    Timer_Reduce("Physics Displacement Wait/Post",
                 (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }
  if (Timer_Exists("Comm Force Exchange")) {
    Timer_Reduce("Comm Force Exchange", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Comm Displacement Exchange")) {
    Timer_Reduce("Comm Displacement Exchange",
                 (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }

  Timer_Reduce("Read My Forces", (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  Timer_Reduce("Compute addforces s", (TimerKind)(MAX | MIN | AVERAGE),
               comm_solver);
  Timer_Reduce("Compute addforces e", (TimerKind)(MAX | MIN | AVERAGE),
               comm_solver);
  Timer_Reduce("Damping addforce", (TimerKind)(MAX | MIN | AVERAGE),
               comm_solver);
  Timer_Reduce("1st schedule send data (contribution)",
               (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  if (Timer_Exists("1st compute adjust (distribution)")) {
    Timer_Reduce("1st compute adjust (distribution)",
                 (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }
  Timer_Reduce("2nd schedule send data (contribution)",
               (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  Timer_Reduce("Compute new displacement", (TimerKind)(MAX | MIN | AVERAGE),
               comm_solver);
  if (Timer_Exists("Disp kernel enqueue")) {
    Timer_Reduce("Disp kernel enqueue", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Disp host force reset")) {
    Timer_Reduce("Disp host force reset", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Disp force upload safe wait")) {
    Timer_Reduce("Disp force upload safe wait",
                 (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }
  if (Timer_Exists("Disp host force clear only")) {
    Timer_Reduce("Disp host force clear only",
                 (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }
  if (Timer_Exists("Disp kernel sync wait")) {
    Timer_Reduce("Disp kernel sync wait", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Disp kernel wait only")) {
    Timer_Reduce("Disp kernel wait only", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Disp postcopy sync wait")) {
    Timer_Reduce("Disp postcopy sync wait", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Disp boundary host scatter")) {
    Timer_Reduce("Disp boundary host scatter",
                 (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }
  Timer_Reduce("3rd schedule send data (sharing)",
               (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  Timer_Reduce("2nd compute adjust (assignment)",
               (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  Timer_Reduce("4th schedule send data (sharing)",
               (TimerKind)(MAX | MIN | AVERAGE), comm_solver);

  if (Timer_Exists("Schedule Post Recv")) {
    Timer_Reduce("Schedule Post Recv", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Schedule Pack Outgoing")) {
    Timer_Reduce("Schedule Pack Outgoing", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Schedule Post Send")) {
    Timer_Reduce("Schedule Post Send", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Schedule Wait Recv")) {
    Timer_Reduce("Schedule Wait Recv", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Schedule Wait Send")) {
    Timer_Reduce("Schedule Wait Send", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Schedule Unpack Incoming")) {
    Timer_Reduce("Schedule Unpack Incoming", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Adjust Cache Build")) {
    Timer_Reduce("Adjust Cache Build", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Adjust Distribution Core")) {
    Timer_Reduce("Adjust Distribution Core", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Adjust Assignment Core")) {
    Timer_Reduce("Adjust Assignment Core", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }

  Timer_Reduce("Solver I/O", (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  Timer_Reduce("Compute Physics", (TimerKind)(MAX | MIN | AVERAGE),
               comm_solver);
  Timer_Reduce("Communication", (TimerKind)(MAX | MIN | AVERAGE), comm_solver);

  solver_diag_reduce();
}

/**
 * March forward in time and output the result whenever necessary.
 */
static void solver_run() {
  int32_t step, startingStep;

  solver_run_init_comm(Global.mySolver);
  solver_diag_reset();
  solver_diag_prepare_local(Global.myMesh);

  /* sets new starting step if loading checkpoint */
  if (Param.theUseCheckPoint == 1) {
    startingStep = checkpoint_read(
        Global.myID, Global.myMesh, Param.theCheckPointingDirOut,
        Global.theGroupSize, Global.mySolver, comm_solver);
  } else {
    startingStep = 0;
  }

  if (Global.myID == 0) {
    /* print header for monitor file */
    monitor_print("solver_run() start\nStarting time step = %d\n\n",
                  startingStep);
    monitor_print(
        "Sim = Simulation time (s), Sol = Solver time (s), WC = Wall Clock "
        "Time (s)\n");
  }

  MPI_Barrier(comm_solver);

  /* march forward in time */
  for (step = startingStep; step < Param.theTotalSteps; step++) {
    fvector_t* tmpvector;

    /* prepare for a new iteration
     * swap displacement vectors for t(n) and t(n-1) */
    tmpvector = Global.mySolver->tm2;
    Global.mySolver->tm2 = Global.mySolver->tm1;
    Global.mySolver->tm1 = tmpvector;
    solver_swap_pressure_buffers(Global.mySolver);

    Timer_Start("GPU Memory Copy");
    Timer_Start("GPU Load Disp");
    solver_load_disp_gpu(Global.mySolver, Global.myMesh);
    solver_load_pressure_gpu(Global.mySolver, Global.myMesh);
    Timer_Stop("GPU Load Disp");
    Timer_Stop("GPU Memory Copy");
    solver_update_acoustic_station_kinematics(step);
    planes_update_acoustic_kinematics(Global.mySolver, Param.theNumberOfPlanes,
                                      step, Param.theDeltaT,
                                      Param.theDeltaTSquared);

    Timer_Start("Solver I/O");
    solver_write_checkpoint(step, startingStep);
    solver_update_status(step, startingStep);
    // Todo Dorian. Ask Ricardo about this file (disp.out), it gets really big
    // at high frequencies, and it is not clear its purpose.
    // solver_output_wavefield( step );
    solver_output_planes(Global.mySolver, Global.myID, step);
    solver_output_stations(step);
    solver_output_drm_nodes(Global.mySolver, step, Param.theTotalSteps);
    solver_read_source_forces(step, Param.theDeltaT);
    solver_read_drm_displacements_v2(step, Param.theDeltaT,
                                     Param.theTotalSteps);
    Timer_Stop("Solver I/O");

    solver_report_force_paths_once();

    Timer_Start("Compute Physics");
    Timer_Start("Physics Force Prep");
    solver_nonlinear_state(Global.mySolver, Global.myMesh, Global.theK1,
                           Global.theK2, step);
    if (use_gpu_source_fastpath()) {
      solver_compute_force_source_gpu(step);
    } else {
      solver_compute_force_source(step);
    }
    // solver_compute_effective_drm_force( Global.mySolver,
    // Global.myMesh,Global.theK1, Global.theK2, step, Param.theDeltaT );
    solver_compute_effective_drm_force_v2(
        Global.mySolver, Global.myMesh, Global.theK1, Global.theK2, step,
        Param.theDeltaT, Param.theTypeOfDamping, Param.theFreq,
        Param.theDeltaTSquared);

    if (topography_has_host_force_path()) {
      solver_compute_force_topography(Global.mySolver, Global.myMesh,
                                      Param.theDeltaTSquared);
    }

    /* Push host-side source/DRM/topography force terms to device before
     * GPU stiffness/damping accumulation. */
    if (need_host_force_upload()) {
      if (use_sparse_force_upload_topo()) {
        solver_load_forces_sparse_gpu(Global.mySolver, Global.myMesh);
      } else {
        solver_load_forces_gpu(Global.mySolver, Global.myMesh);
      }
    }

    if (use_gpu_topography_fastpath()) {
      solver_compute_force_topography(Global.mySolver, Global.myMesh,
                                      Param.theDeltaTSquared);
    }
    Timer_Stop("Physics Force Prep");

    Timer_Start("Physics Force Kernels Enqueue");
    solver_compute_force_stiffness(Global.mySolver, Global.myMesh, Global.theK1,
                     Global.theK2);
    solver_compute_force_damping(Global.mySolver, Global.myMesh, Global.theK1,
                   Global.theK2);
    solver_compute_force_acoustic_pressure(Global.mySolver, Global.myMesh);
    solver_record_force_kernel_done_event(Global.mySolver);
    Timer_Stop("Physics Force Kernels Enqueue");

    Timer_Stop("Compute Physics");

    Timer_Start("GPU Memory Copy");
    Timer_Start("GPU Unload Forces");
    solver_unload_forces_gpu(Global.mySolver, Global.myMesh);
    Timer_Stop("GPU Unload Forces");
    Timer_Stop("GPU Memory Copy");

    Timer_Start("Compute Physics");
    Timer_Start("Physics Force Wait/Post");
    solver_wait_physics_gpu(Global.mySolver, Global.myMesh);

    Timer_Start("Force Post CPU");
    solver_compute_force_gravity(Global.mySolver, Global.myMesh, step);
    solver_compute_force_nonlinear(Global.mySolver, Global.myMesh,
                                   Param.theDeltaTSquared);
    solver_compute_force_planewaves(Global.myMesh, Global.mySolver,
                                    Param.theDeltaT, step, Global.theK1,
                                    Global.theK2);
    Timer_Stop("Force Post CPU");
    Timer_Stop("Physics Force Wait/Post");

    Timer_Stop("Compute Physics");

    Timer_Start("Communication");
    Timer_Start("Comm Force Exchange");
    HU_COND_GLOBAL_BARRIER(Param.theTimingBarriersFlag);
    solver_send_force_dangling(Global.mySolver);
    solver_adjust_forces(Global.mySolver);
    HU_COND_GLOBAL_BARRIER(Param.theTimingBarriersFlag);
    solver_send_force_anchored(Global.mySolver);
    Timer_Stop("Comm Force Exchange");
    Timer_Stop("Communication");

    Timer_Start("GPU Memory Copy");
    Timer_Start("GPU Load Forces");
    solver_load_forces_gpu(Global.mySolver, Global.myMesh);
    Timer_Stop("GPU Load Forces");
    Timer_Stop("GPU Memory Copy");

    Timer_Start("Compute Physics");
    Timer_Start("Physics Displacement Enqueue");
    /* solver_compute_displacement(Global.mySolver, Global.myMesh); */
    solver_compute_displacement_gpu(Global.myID, Global.mySolver,
                                    Global.myMesh);
    Timer_Stop("Physics Displacement Enqueue");
    Timer_Stop("Compute Physics");

    Timer_Start("GPU Memory Copy");
    Timer_Start("GPU Unload Disp");
    solver_unload_disp_gpu(Global.mySolver, Global.myMesh);
    Timer_Stop("GPU Unload Disp");
    Timer_Stop("GPU Memory Copy");

    Timer_Start("Compute Physics");
    Timer_Start("Physics Displacement Wait/Post");
    solver_wait_disp_gpu(Global.mySolver, Global.myMesh);
    solver_run_acoustic_pressure_step(Global.mySolver, Global.myMesh);
    solver_geostatic_fix(step);
    solver_load_fixedbase_displacements(Global.mySolver, step);
    Timer_Stop("Physics Displacement Wait/Post");
    Timer_Stop("Compute Physics");

    Timer_Start("GPU Memory Copy");
    Timer_Start("GPU Reset");
    solver_reset_gpu(Global.mySolver, Global.myMesh);
    Timer_Stop("GPU Reset");
    Timer_Stop("GPU Memory Copy");

    Timer_Start("Communication");
    Timer_Start("Comm Displacement Exchange");
    HU_COND_GLOBAL_BARRIER(Param.theTimingBarriersFlag);
    solver_send_displacement_anchored(Global.mySolver);
    solver_adjust_displacement(Global.mySolver);
    HU_COND_GLOBAL_BARRIER(Param.theTimingBarriersFlag);
    solver_send_displacement_dangling(Global.mySolver);
    Timer_Stop("Comm Displacement Exchange");
    Timer_Stop("Communication");

    solver_upload_disp_boundary_gpu(Global.mySolver);

    solver_loop_hook_bottom(Global.mySolver, Global.myMesh, step);
    g_dnDiagLocal[DN_DIAG_STEPS] += 1.0;
  } /* for (step = ....): all steps */

  solver_drm_close();
  solver_output_wavefield_close();
  solver_run_collect_timers();
}

/**
 * Output the velocity of the mesh nodes for the current timestep. Send
 * the data to Thread 0 who is responsible for dumping the data to disk.
 *
 * \note This is only useful for very small distributed setups where:
 * - There is no parallel file system.
 * - The distributed file system does not have POSIX semantics.
 * - The distributed file system performs poorly with just a few clients
 *   and the effective overall performance is better when only one client
 *   (i.e., PE 0) is writing to the FS.
 */
void solver_output_seq() {
  int32_t nindex;
  int32_t batchlimit, idx;

#ifdef DEBUG
  int64_t gnid_prev, gnid_current;
  int32_t first_counted;
#endif /* DEBUG */

  batchlimit = BATCH * 10;

  /* Allocate a fixed size buffer space if not initiazlied */
  if (Global.myVelocityTable == NULL) {
    Global.myVelocityTable = (fvector_t*)calloc(batchlimit, sizeof(fvector_t));
    if (Global.myVelocityTable == NULL) {
      fprintf(stderr, "Thread %d: solver_output_seq: out of memory\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
  }

  if (Global.myID == 0) {
    int32_t procid;
#ifdef DEBUG
    first_counted = 0;
#endif

    if (Param.FourDOutFp == NULL) {
      out_hdr_t out_hdr;

      /* First output, create the output file */
      Param.FourDOutFp = fopen(Param.FourDOutFile, "w+");
      if (Param.FourDOutFp == NULL) {
        fprintf(stderr, "Thread 0: solver_output_seq: ");
        fprintf(stderr, "cannot create %s\n", Param.FourDOutFile);
        perror("fopen");
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }

      /* Write the header that contains the metadata */
      out_hdr.domain_x = Param.theDomainX;
      out_hdr.domain_y = Param.theDomainY;
      out_hdr.domain_z = Param.theDomainZ;
      out_hdr.total_nodes = Global.theNTotal;
      out_hdr.total_elements = Global.theETotal;
      out_hdr.mesh_ticksize = Global.myMesh->ticksize;
      out_hdr.output_steps = (Param.theTotalSteps - 1) / Param.theRate + 1;

      if (fwrite(&out_hdr, sizeof(out_hdr_t), 1, Param.FourDOutFp) != 1) {
        fprintf(stderr, "Thread 0: solver_output_seq: ");
        fprintf(stderr, "fail to write 4D-out header info\n");
        perror("fwrite");
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }
    }

    /* Output my owned nodes' velocities */
    nindex = 0;
    while (nindex < Global.myMesh->nharbored) {
      fvector_t vel;

      if (Global.myMesh->nodeTable[nindex].ismine) {
        vel.f[0] = (Global.mySolver->tm1[nindex].f[0] -
                    Global.mySolver->tm2[nindex].f[0]) /
                   Param.theDeltaT;
        vel.f[1] = (Global.mySolver->tm1[nindex].f[1] -
                    Global.mySolver->tm2[nindex].f[1]) /
                   Param.theDeltaT;
        vel.f[2] = (Global.mySolver->tm1[nindex].f[2] -
                    Global.mySolver->tm2[nindex].f[2]) /
                   Param.theDeltaT;

        if (fwrite(&vel, sizeof(fvector_t), 1, Param.FourDOutFp) != 1) {
          fprintf(stderr, "Thread 0: solver_output_seq: error\n");
          MPI_Abort(MPI_COMM_WORLD, ERROR);
          exit(1);
        }

        Param.the4DOutSize += sizeof(fvector_t);

#ifdef DEBUG
        gnid_current = Global.myMesh->nodeTable[nindex].gnid;

        if (first_counted) {
          if (gnid_prev != (gnid_current - 1)) {
            fprintf(stderr,
                    "PE 0: uncontinuous gnid\n"
                    "   gnid_prev = %" INT64_FMT ", gnid_current = %" INT64_FMT
                    "\n",
                    gnid_prev, gnid_current);
          }
        } else {
          first_counted = 1;
        }

        gnid_prev = gnid_current;

        if ((vel.f[0] != 0) || (vel.f[1] != 0) || (vel.f[2] != 0)) {
          /*
          fprintf(stderr, "Thread 0: Node %ld  non-zero values\n",
              gnid_current);
          */
        }
#endif /* DEBUG */
      }

      nindex++;
    }

    /* Receive data from other processors */
    for (procid = 1; procid < Global.theGroupSize; procid++) {
      MPI_Status status;
      int32_t rcvbytecount;

      /* Signal the next processor to go ahead */
      MPI_Send(NULL, 0, MPI_CHAR, procid, GOAHEAD_MSG, comm_solver);

      while (1) {
        MPI_Probe(procid, OUT4D_MSG, comm_solver, &status);
        MPI_Get_count(&status, MPI_CHAR, &rcvbytecount);

        /* Receive the data even if rcvbytecount == 0. Otherwise
           the 0-byte message would get stuck in the message queue */
        MPI_Recv(Global.myVelocityTable, rcvbytecount, MPI_CHAR, procid,
                 OUT4D_MSG, comm_solver, &status);

        if (rcvbytecount == 0) {
          /* Done */
          break;
        }

        if (fwrite(Global.myVelocityTable, rcvbytecount, 1, Param.FourDOutFp) !=
            1) {
          fprintf(stderr, "Thread 0: solver_output_seq: error\n");
          MPI_Abort(MPI_COMM_WORLD, ERROR);
          exit(1);
        }

        Param.the4DOutSize += rcvbytecount;

      } /* while there is more data to be received from procid */
    } /* for all the processors */
  } else {
    /* Processors other than 0 needs to send data to 0 */
    int32_t sndbytecount;
    MPI_Status status;

    /* Wait for me turn */
    MPI_Recv(NULL, 0, MPI_CHAR, 0, GOAHEAD_MSG, comm_solver, &status);

#ifdef DEBUG
    first_counted = 0;
#endif

    nindex = 0;
    while (nindex < Global.myMesh->nharbored) {
      fvector_t* velp;

      idx = 0;
      while ((idx < batchlimit) && (nindex < Global.myMesh->nharbored)) {
        if (Global.myMesh->nodeTable[nindex].ismine) {
          velp = &Global.myVelocityTable[idx];

          velp->f[0] = (Global.mySolver->tm1[nindex].f[0] -
                        Global.mySolver->tm2[nindex].f[0]) /
                       Param.theDeltaT;
          velp->f[1] = (Global.mySolver->tm1[nindex].f[1] -
                        Global.mySolver->tm2[nindex].f[1]) /
                       Param.theDeltaT;
          velp->f[2] = (Global.mySolver->tm1[nindex].f[2] -
                        Global.mySolver->tm2[nindex].f[2]) /
                       Param.theDeltaT;

          idx++;

#ifdef DEBUG
          gnid_current = Global.myMesh->nodeTable[nindex].gnid;

          if (first_counted) {
            if (gnid_prev != (gnid_current - 1)) {
              fprintf(stderr,
                      "PE %d uncontinuous gnid\n"
                      "  gnid_prev = %" INT64_FMT ", gnid_current = %" INT64_FMT
                      "\n",
                      Global.myID, gnid_prev, gnid_current);
            }
          } else {
            first_counted = 1;
          }

          gnid_prev = gnid_current;

          /* debug */
          /*
          if ((velp->f[0] != 0) ||
          (velp->f[1] != 0) ||
          (velp->f[2] != 0)) {
          fprintf(stderr,
              "Thread %d: there are non-zero values\n",
              Global.myID);
              }
          */
#endif /* DEBUG */
        }

        nindex++;
      }

      /* Send data to proc 0 */

      if (idx > 0) {
        /* I have some real data to send */
        sndbytecount = idx * sizeof(fvector_t);
        MPI_Send(Global.myVelocityTable, sndbytecount, MPI_CHAR, 0, OUT4D_MSG,
                 comm_solver);
      }
    } /* While there is data left to be sent */

    /* Send an empty message to indicate the end of my transfer */
    MPI_Send(NULL, 0, MPI_CHAR, 0, OUT4D_MSG, comm_solver);
  }

  return;
}

/**
 * Allocate and initialize a scheduler.
 */
static schedule_t* schedule_new() {
  schedule_t* sched;

  sched = (schedule_t*)malloc(sizeof(schedule_t));
  if (sched == NULL) return NULL;

  sched->c_count = 0;
  sched->first_c = NULL;
  sched->messenger_c =
      (messenger_t**)calloc(Global.theGroupSize, sizeof(messenger_t*));
  if (sched->messenger_c == NULL) return NULL;

  sched->s_count = 0;
  sched->first_s = NULL;
  sched->messenger_s =
      (messenger_t**)calloc(Global.theGroupSize, sizeof(messenger_t*));
  if (sched->messenger_s == NULL) return NULL;

  sched->irecvreqs_c = NULL;
  sched->irecvstats_c = NULL;
  sched->waitsome_idxs_c = NULL;
  sched->isendreqs_c = NULL;
  sched->isendstats_c = NULL;
  sched->irecvreqs_s = NULL;
  sched->irecvstats_s = NULL;
  sched->waitsome_idxs_s = NULL;
  sched->isendreqs_s = NULL;
  sched->isendstats_s = NULL;

  return sched;
}

/**
 * Allocate the mapping table for each of my messenger.
 */
static void schedule_allocmapping(schedule_t* sched) {
  messenger_t* messenger;
  int32_t nodecount;

  messenger = sched->first_c;
  while (messenger != NULL) {
    nodecount = messenger->nodecount;

    messenger->mapping = (int32_t*)calloc(nodecount, sizeof(int32_t));

    if (messenger->mapping == NULL) {
      fprintf(stderr, "Thread %d: schedule_allocamapping: ", Global.myID);
      fprintf(stderr, " out of memory\n");
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

    messenger = messenger->next;
  }

  messenger = sched->first_s;
  while (messenger != NULL) {
    nodecount = messenger->nodecount;

    messenger->mapping = (int32_t*)calloc(nodecount, sizeof(int32_t));

    if (messenger->mapping == NULL) {
      fprintf(stderr, "Thread %d: schedule_allocamapping: ", Global.myID);
      fprintf(stderr, " out of memory\n");
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

    messenger = messenger->next;
  }

  return;
}

/**
 * Allocate MPI controls for non-blocing receives.
 */
static void schedule_allocMPIctl(schedule_t* sched) {
  if (sched->c_count != 0) {
    sched->irecvreqs_c =
        (MPI_Request*)malloc(sizeof(MPI_Request) * sched->c_count);
    sched->irecvstats_c =
        (MPI_Status*)malloc(sizeof(MPI_Status) * sched->c_count);
    sched->waitsome_idxs_c = (int*)malloc(sizeof(int) * sched->c_count);
    sched->isendreqs_c =
        (MPI_Request*)malloc(sizeof(MPI_Request) * sched->c_count);
    sched->isendstats_c =
        (MPI_Status*)malloc(sizeof(MPI_Status) * sched->c_count);

    if ((sched->irecvreqs_c == NULL) || (sched->irecvstats_c == NULL) ||
      (sched->waitsome_idxs_c == NULL) ||
        (sched->isendreqs_c == NULL) || (sched->isendstats_c == NULL)) {
      fprintf(stderr, "Thread %d: schedule_allocMPIctl: ", Global.myID);
      fprintf(stderr, "out of memory\n");
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
  } else {
    sched->irecvreqs_c = NULL;
    sched->irecvstats_c = NULL;
    sched->waitsome_idxs_c = NULL;
    sched->isendreqs_c = NULL;
    sched->isendstats_c = NULL;
  }

  if (sched->s_count != 0) {
    sched->irecvreqs_s =
        (MPI_Request*)malloc(sizeof(MPI_Request) * sched->s_count);
    sched->irecvstats_s =
        (MPI_Status*)malloc(sizeof(MPI_Status) * sched->s_count);
    sched->waitsome_idxs_s = (int*)malloc(sizeof(int) * sched->s_count);
    sched->isendreqs_s =
        (MPI_Request*)malloc(sizeof(MPI_Request) * sched->s_count);
    sched->isendstats_s =
        (MPI_Status*)malloc(sizeof(MPI_Status) * sched->s_count);

    if ((sched->irecvreqs_s == NULL) || (sched->irecvstats_s == NULL) ||
      (sched->waitsome_idxs_s == NULL) ||
        (sched->isendreqs_s == NULL) || (sched->isendstats_s == NULL)) {
      fprintf(stderr, "Thread %d: schedule_allocMPIctl: ", Global.myID);
      fprintf(stderr, "out of memory\n");
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
  } else {
    sched->irecvreqs_s = NULL;
    sched->irecvstats_s = NULL;
    sched->waitsome_idxs_s = NULL;
    sched->isendreqs_s = NULL;
    sched->isendstats_s = NULL;
  }

  return;
}

/**
 * Build a communication schedule using local information.
 */
static void schedule_build(mesh_t* mesh, schedule_t* dnsched,
                           schedule_t* ansched) {
  int32_t nindex;
  node_t* nodep;
  messenger_t* messenger;

  for (nindex = 0; nindex < mesh->nharbored; nindex++) {
    nodep = &mesh->nodeTable[nindex];
    int32_t owner, sharer;

    if (!nodep->ismine) {
      /* I do not own this node. Add its owner processor to my c-list */

      owner = nodep->proc.ownerid;

      if (nodep->isanchored) {
        messenger = ansched->messenger_c[owner];
      } else {
        messenger = dnsched->messenger_c[owner];
      }

      if (messenger == NULL) {
        messenger = messenger_new(owner);
        if (messenger == NULL) {
          fprintf(stderr, "Thread %d: schedule_build: ", Global.myID);
          fprintf(stderr, "out of memory.\n");
          MPI_Abort(MPI_COMM_WORLD, ERROR);
          exit(1);
        }

        if (nodep->isanchored) {
          ansched->c_count++;
          ansched->messenger_c[owner] = messenger;
          messenger->next = ansched->first_c;
          ansched->first_c = messenger;
        } else {
          dnsched->c_count++;
          dnsched->messenger_c[owner] = messenger;
          messenger->next = dnsched->first_c;
          dnsched->first_c = messenger;
        }
      }

      /* Update the number of nodecount for the messenger */
      messenger->nodecount++;
    } else {
      /* I own this node. Add any sharing processor to my s-list */

      int32link_t* int32link;

      int32link = nodep->proc.share;
      while (int32link != NULL) {
        sharer = int32link->id;

        if (nodep->isanchored) {
          messenger = ansched->messenger_s[sharer];
        } else {
          messenger = dnsched->messenger_s[sharer];
        }

        if (messenger == NULL) {
          messenger = messenger_new(sharer);
          if (messenger == NULL) {
            fprintf(stderr, "Thread %d: schedule_build: ", Global.myID);
            fprintf(stderr, "out of memory.\n");
            MPI_Abort(MPI_COMM_WORLD, ERROR);
            exit(1);
          }

          if (nodep->isanchored) {
            ansched->s_count++;
            ansched->messenger_s[sharer] = messenger;
            messenger->next = ansched->first_s;
            ansched->first_s = messenger;
          } else {
            dnsched->s_count++;
            dnsched->messenger_s[sharer] = messenger;
            messenger->next = dnsched->first_s;
            dnsched->first_s = messenger;
          }
        }

        /* Update the nodecount */
        messenger->nodecount++;

        /* Move to the next sharing processor */
        int32link = int32link->next;
      }
    }
  }

  /* Allocate MPI controls */
  schedule_allocMPIctl(ansched);
  schedule_allocMPIctl(dnsched);

  /* Allocate localnode table for each of the messegners I have */
  schedule_allocmapping(ansched);
  schedule_allocmapping(dnsched);

  /* Go through the nodes again and fill out the mapping table */
  for (nindex = 0; nindex < mesh->nharbored; nindex++) {
    nodep = &mesh->nodeTable[nindex];
    int32_t owner, sharer;

    if (!nodep->ismine) {
      /* I do not own this node. Add its owner processor to my c-list */
      owner = nodep->proc.ownerid;

      if (nodep->isanchored) {
        messenger = ansched->messenger_c[owner];
      } else {
        messenger = dnsched->messenger_c[owner];
      }

      if (messenger == NULL) {
        fprintf(stderr, "Thread %d: schedule_build: ", Global.myID);
        fprintf(stderr, "encounter NULL messenger.\n");
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }

      /* Fill in the mapping */
      messenger->mapping[messenger->nidx] = nindex;
      messenger->nidx++;
    } else {
      /* I own this node. Add any sharing processor to my s-list */

      int32link_t* int32link;

      int32link = nodep->proc.share;
      while (int32link != NULL) {
        sharer = int32link->id;

        if (nodep->isanchored) {
          messenger = ansched->messenger_s[sharer];
        } else {
          messenger = dnsched->messenger_s[sharer];
        }

        if (messenger == NULL) {
          fprintf(stderr, "Thread %d: schedule_build: ", Global.myID);
          fprintf(stderr, "encounter NULL messenger.\n");
          MPI_Abort(MPI_COMM_WORLD, ERROR);
          exit(1);
        }

        messenger->mapping[messenger->nidx] = nindex;
        messenger->nidx++;

        /* Move to the next sharing processor */
        int32link = int32link->next;
      }
    }
  }

  return;
}

/**
 * Release memory used by a scheduler.
 */
static void schedule_delete(schedule_t* sched) {
  messenger_t *current, *next;

  /* Release messengers overseeing my contributions */
  current = sched->first_c;
  while (current != NULL) {
    next = current->next;

    messenger_delete(current);
    current = next;
  }

  /* Release messengers overseeing shareing with others */
  current = sched->first_s;
  while (current != NULL) {
    next = current->next;

    messenger_delete(current);
    current = next;
  }

  if (sched->irecvreqs_c != NULL) free(sched->irecvreqs_c);
  if (sched->irecvstats_c != NULL) free(sched->irecvstats_c);
  if (sched->waitsome_idxs_c != NULL) free(sched->waitsome_idxs_c);
  if (sched->isendreqs_c != NULL) free(sched->isendreqs_c);
  if (sched->isendstats_c != NULL) free(sched->isendstats_c);
  free(sched->messenger_c);

  if (sched->irecvreqs_s != NULL) free(sched->irecvreqs_s);
  if (sched->irecvstats_s != NULL) free(sched->irecvstats_s);
  if (sched->waitsome_idxs_s != NULL) free(sched->waitsome_idxs_s);
  if (sched->isendreqs_s != NULL) free(sched->isendreqs_s);
  if (sched->isendstats_s != NULL) free(sched->isendstats_s);
  free(sched->messenger_s);

  free(sched);

  return;
}

/**
 * Allocate the memory for data exchange.
 */
static void schedule_prepare(schedule_t* sched, int32_t c_outsize,
                             int32_t c_insize, int32_t s_outsize,
                             int32_t s_insize) {
  messenger_t* messenger;

  messenger = sched->first_c;
  while (messenger != NULL) {
    messenger_set(messenger, c_outsize, c_insize);
    messenger = messenger->next;
  }

  messenger = sched->first_s;
  while (messenger != NULL) {
    messenger_set(messenger, s_outsize, s_insize);
    messenger = messenger->next;
  }

  return;
}

/**
 * Assemble the proper information for the group of messengers managed by
 * a scheduler and send the data.
 *
 * \param direction: CONTRIBUTION or SHARING.
 */
static void schedule_senddata(schedule_t* sched, void* valuetable,
                              int32_t itemsperentry, int32_t direction,
                              int32_t msgtag) {
  messenger_t *send_messenger, *recv_messenger;
  messenger_t *send_head, *recv_head;
  messenger_t** recv_order = NULL;
  int32_t irecvcount, irecvnum, bytesize;
  int32_t isendcount, isendnum;
  double callSendBytes = 0.0;
  double callRecvBytes = 0.0;
  double callSendNodes = 0.0;
  double callRecvNodes = 0.0;
  double callSendMsgs = 0.0;
  double callRecvMsgs = 0.0;
  MPI_Request* irecvreqs;
  MPI_Status* irecvstats;
  MPI_Request* isendreqs;
  MPI_Status* isendstats;
  int32_t waitStage = wait_recv_stage_from_msgtag(msgtag);
  int32_t useWaitsome = 0;

#ifdef DEBUG
  int64_t* gnidp;
#endif /* DEBUG */

  if (direction == CONTRIBUTION) {
    send_messenger = sched->first_c;
    recv_messenger = sched->first_s;
    irecvcount = sched->s_count;
    isendcount = sched->c_count;
    irecvreqs = sched->irecvreqs_s;
    irecvstats = sched->irecvstats_s;
    isendreqs = sched->isendreqs_c;
    isendstats = sched->isendstats_c;
  } else {
    send_messenger = sched->first_s;
    recv_messenger = sched->first_c;
    irecvcount = sched->c_count;
    isendcount = sched->s_count;
    irecvreqs = sched->irecvreqs_c;
    irecvstats = sched->irecvstats_c;
    isendreqs = sched->isendreqs_s;
    isendstats = sched->isendstats_s;
  }

  send_head = send_messenger;
  recv_head = recv_messenger;
  /* Sharing messages are assignments and can stream unpack safely.
   * Contribution messages use += accumulation; in performance mode we stream
   * them too to overlap unpack work with MPI progress, while repro mode keeps
   * the deterministic Waitall/unpack ordering. */
  useWaitsome =
      ((direction == SHARING) && !hercules_repro_mode()) ||
      ((direction == CONTRIBUTION) && contrib_recv_use_waitsome());

  messenger_prepare_offsets(send_head, itemsperentry);
  messenger_prepare_offsets(recv_head, itemsperentry);

  /* Post receives */
  double stageTimerStart = 0.0;
  if (waitStage >= 0) {
    stageTimerStart = MPI_Wtime();
  }
  Timer_Start("Schedule Post Recv");
  irecvnum = 0;
  if (irecvcount > 0) {
    recv_order = (messenger_t**)malloc(sizeof(messenger_t*) * irecvcount);
    if (recv_order == NULL) {
      fprintf(stderr, "Thread %d: schedule_senddata: out of memory\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
  }
  while (recv_messenger != NULL) {
    bytesize = recv_messenger->nodecount * recv_messenger->insize;
    callRecvBytes += (double)bytesize;
    callRecvNodes += (double)recv_messenger->nodecount;
    callRecvMsgs += 1.0;
    recv_order[irecvnum] = recv_messenger;
    MPI_Irecv(recv_messenger->indata, bytesize, MPI_CHAR,
              recv_messenger->procid, msgtag, comm_solver,
              &irecvreqs[irecvnum]);

    irecvnum++;
    recv_messenger = recv_messenger->next;
  }
  Timer_Stop("Schedule Post Recv");
  if (waitStage >= 0) {
    g_stageTimeLocal[stage_time_index(waitStage, STAGE_TIME_POST_RECV)] +=
        MPI_Wtime() - stageTimerStart;
  }

  /* Assemble each outgoing message and post it immediately.  Starting large
   * sends while later messages are still being packed helps hide some of the
   * communication tail on imbalanced dangling-node schedules. */
  isendnum = 0;
  while (send_messenger != NULL) {
    int32_t idx;
    int32_t useOmp = 0;

    if (direction == SHARING) {
      useOmp = 1;
    } else if (direction == CONTRIBUTION) {
      useOmp = contrib_pack_use_omp();
    }

    if (waitStage >= 0) {
      stageTimerStart = MPI_Wtime();
    }
    Timer_Start("Schedule Pack Outgoing");
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if((send_messenger->nidx > 10000) && useOmp)
#endif
    for (idx = 0; idx < send_messenger->nidx; idx++) {
      int32_t lnid, entry;
      solver_float* dvalue;
      solver_float* out;

      lnid = send_messenger->mapping[idx];

      out = (solver_float*)((char*)send_messenger->outdata +
                            send_messenger->outsize * idx);

      dvalue =
          (solver_float*)valuetable + send_messenger->mapping_itemoffset[idx];
      if (itemsperentry == 3) {
        out[0] = dvalue[0];
        out[1] = dvalue[1];
        out[2] = dvalue[2];
      } else {
        for (entry = 0; entry < itemsperentry; entry++)
          *(out + entry) = *(dvalue + entry);
      }

#ifdef DEBUG
      /* For debug, carry the global node id */
      gnidp = (int64_t*)((char*)out + itemsperentry * sizeof(solver_float));
      *gnidp = Global.myMesh->nodeTable[lnid].gnid;
#endif /* DEBUG */
    }
    Timer_Stop("Schedule Pack Outgoing");
    if (waitStage >= 0) {
      g_stageTimeLocal[stage_time_index(waitStage,
                                        STAGE_TIME_PACK_OUTGOING)] +=
          MPI_Wtime() - stageTimerStart;
    }

    if (waitStage >= 0) {
      stageTimerStart = MPI_Wtime();
    }
    Timer_Start("Schedule Post Send");
    bytesize = send_messenger->nodecount * send_messenger->outsize;
    callSendBytes += (double)bytesize;
    callSendNodes += (double)send_messenger->nodecount;
    callSendMsgs += 1.0;
    MPI_Isend(send_messenger->outdata, bytesize, MPI_CHAR,
              send_messenger->procid, msgtag, comm_solver,
              &isendreqs[isendnum]);
    Timer_Stop("Schedule Post Send");
    if (waitStage >= 0) {
      g_stageTimeLocal[stage_time_index(waitStage, STAGE_TIME_POST_SEND)] +=
          MPI_Wtime() - stageTimerStart;
    }

    isendnum++;
    send_messenger = send_messenger->next;
  }

  /* Wait till I receive all the data I want */
  if (irecvcount != 0) {
    int32_t remaining = irecvcount;
    double waitRecvStart = 0.0;
    int* doneidxs = NULL;
    if (waitStage >= 0) {
      waitRecvStart = MPI_Wtime();
    }

    if (useWaitsome) {
      doneidxs = (direction == CONTRIBUTION) ? sched->waitsome_idxs_s
                                             : sched->waitsome_idxs_c;
      if (doneidxs == NULL) {
        fprintf(stderr, "Thread %d: schedule_senddata: out of memory\n",
                Global.myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }
    }

    while (remaining > 0) {
      if (useWaitsome) {
        int outcount;
        int i;

        if (waitStage >= 0) {
          stageTimerStart = MPI_Wtime();
        }
        Timer_Start("Schedule Wait Recv");
        MPI_Waitsome(irecvcount, irecvreqs, &outcount, doneidxs, irecvstats);
        Timer_Stop("Schedule Wait Recv");
        if (waitStage >= 0) {
          g_stageTimeLocal[stage_time_index(waitStage,
                                            STAGE_TIME_WAIT_RECV)] +=
              MPI_Wtime() - stageTimerStart;
        }

        if (outcount == MPI_UNDEFINED) {
          break;
        }

        if (waitStage >= 0) {
          stageTimerStart = MPI_Wtime();
        }
        Timer_Start("Schedule Unpack Incoming");
        for (i = 0; i < outcount; i++) {
          int doneidx = doneidxs[i];
          schedule_unpack_incoming_message(recv_order[doneidx], valuetable,
                                           itemsperentry, direction);
        }
        Timer_Stop("Schedule Unpack Incoming");
        if (waitStage >= 0) {
          g_stageTimeLocal[stage_time_index(waitStage,
                                            STAGE_TIME_UNPACK_INCOMING)] +=
              MPI_Wtime() - stageTimerStart;
        }
        remaining -= outcount;
      } else if ((direction == CONTRIBUTION) || hercules_repro_mode()) {
        int i;

        if (waitStage >= 0) {
          stageTimerStart = MPI_Wtime();
        }
        Timer_Start("Schedule Wait Recv");
        MPI_Waitall(irecvcount, irecvreqs, irecvstats);
        Timer_Stop("Schedule Wait Recv");
        if (waitStage >= 0) {
          g_stageTimeLocal[stage_time_index(waitStage,
                                            STAGE_TIME_WAIT_RECV)] +=
              MPI_Wtime() - stageTimerStart;
        }

        /* Contribution paths use += accumulation, so unpack in a
         * deterministic order independent of receive completion timing. */
        if (waitStage >= 0) {
          stageTimerStart = MPI_Wtime();
        }
        Timer_Start("Schedule Unpack Incoming");
        for (i = 0; i < irecvcount; i++) {
          schedule_unpack_incoming_message(recv_order[i], valuetable,
                                           itemsperentry, direction);
        }
        Timer_Stop("Schedule Unpack Incoming");
        if (waitStage >= 0) {
          g_stageTimeLocal[stage_time_index(waitStage,
                                            STAGE_TIME_UNPACK_INCOMING)] +=
              MPI_Wtime() - stageTimerStart;
        }
        remaining = 0;
      } else {
        int doneidx;
        MPI_Status status;

        if (waitStage >= 0) {
          stageTimerStart = MPI_Wtime();
        }
        Timer_Start("Schedule Wait Recv");
        MPI_Waitany(irecvcount, irecvreqs, &doneidx, &status);
        Timer_Stop("Schedule Wait Recv");
        if (waitStage >= 0) {
          g_stageTimeLocal[stage_time_index(waitStage,
                                            STAGE_TIME_WAIT_RECV)] +=
              MPI_Wtime() - stageTimerStart;
        }

        if (doneidx == MPI_UNDEFINED) {
          break;
        }

        if (waitStage >= 0) {
          stageTimerStart = MPI_Wtime();
        }
        Timer_Start("Schedule Unpack Incoming");
        schedule_unpack_incoming_message(recv_order[doneidx], valuetable,
                                         itemsperentry, direction);
        Timer_Stop("Schedule Unpack Incoming");
        if (waitStage >= 0) {
          g_stageTimeLocal[stage_time_index(waitStage,
                                            STAGE_TIME_UNPACK_INCOMING)] +=
              MPI_Wtime() - stageTimerStart;
        }
        remaining--;
      }
    }

    if (waitStage >= 0) {
      g_waitRecvStageLocal[waitStage] += (MPI_Wtime() - waitRecvStart);
    }
  } else {
    Timer_Start("Schedule Unpack Incoming");
    Timer_Stop("Schedule Unpack Incoming");
  }

  if (isendcount != 0) {
    if (waitStage >= 0) {
      stageTimerStart = MPI_Wtime();
    }
    Timer_Start("Schedule Wait Send");
    MPI_Waitall(isendnum, isendreqs, isendstats);
    Timer_Stop("Schedule Wait Send");
    if (waitStage >= 0) {
      g_stageTimeLocal[stage_time_index(waitStage, STAGE_TIME_WAIT_SEND)] +=
          MPI_Wtime() - stageTimerStart;
    }
  }

  if (recv_order != NULL) {
    free(recv_order);
    recv_order = NULL;
  }

  if (msgtag == DN_FORCE_MSG || msgtag == DN_DISP_MSG) {
    g_dnDiagLocal[DN_DIAG_SEND_BYTES] += callSendBytes;
    g_dnDiagLocal[DN_DIAG_RECV_BYTES] += callRecvBytes;
    g_dnDiagLocal[DN_DIAG_SEND_NODES] += callSendNodes;
    g_dnDiagLocal[DN_DIAG_RECV_NODES] += callRecvNodes;
    g_dnDiagLocal[DN_DIAG_SEND_MSGS] += callSendMsgs;
    g_dnDiagLocal[DN_DIAG_RECV_MSGS] += callRecvMsgs;
  }

    if (waitStage >= 0) {
    g_stageCommLocal[stage_comm_index(waitStage, STAGE_COMM_SEND_BYTES)] +=
      callSendBytes;
    g_stageCommLocal[stage_comm_index(waitStage, STAGE_COMM_RECV_BYTES)] +=
      callRecvBytes;
    g_stageCommLocal[stage_comm_index(waitStage, STAGE_COMM_SEND_MSGS)] +=
      callSendMsgs;
    g_stageCommLocal[stage_comm_index(waitStage, STAGE_COMM_RECV_MSGS)] +=
      callRecvMsgs;
    g_stageCommLocal[stage_comm_index(waitStage, STAGE_COMM_SEND_NODES)] +=
      callSendNodes;
    g_stageCommLocal[stage_comm_index(waitStage, STAGE_COMM_RECV_NODES)] +=
      callRecvNodes;
    }

  /* MPI_Barrier(comm_solver);  */

  return;
}

static void schedule_unpack_incoming_message(messenger_t* recv_messenger,
                                             void* valuetable,
                                             int32_t itemsperentry,
                                             int32_t direction) {
  int32_t idx;

  if (direction == SHARING) {
    const int32_t useOmp = 1;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if((recv_messenger->nidx > 10000) && useOmp && !recv_messenger->mapping_has_duplicates)
#endif
    for (idx = 0; idx < recv_messenger->nidx; idx++) {
      int32_t lnid, entry;
      solver_float* dvalue;
      solver_float* in;

      lnid = recv_messenger->mapping[idx];
      in = (solver_float*)((char*)recv_messenger->indata +
                           recv_messenger->insize * idx);

      dvalue =
          (solver_float*)valuetable + recv_messenger->mapping_itemoffset[idx];
      if (itemsperentry == 3) {
        dvalue[0] = in[0];
        dvalue[1] = in[1];
        dvalue[2] = in[2];
      } else {
        for (entry = 0; entry < itemsperentry; entry++) {
          *(dvalue + entry) = *(in + entry);
        }
      }

#ifdef DEBUG
      {
        int64_t* gnidp;
        gnidp = (int64_t*)((char*)in + itemsperentry * sizeof(solver_float));

        if (*gnidp != Global.myMesh->nodeTable[lnid].gnid) {
          fprintf(stderr, "Thread %d: solver_init: gnids do not match\n",
                  Global.myID);
          MPI_Abort(MPI_COMM_WORLD, ERROR);
          exit(1);
        }
      }
#endif /* DEBUG */
    }

    return;
  }

  {
#ifdef _OPENMP
    const int32_t useOmp = contrib_unpack_use_omp();
#pragma omp parallel for schedule(static) if((recv_messenger->nidx > 10000) && useOmp && !recv_messenger->mapping_has_duplicates)
#endif
    for (idx = 0; idx < recv_messenger->nidx; idx++) {
      int32_t lnid, entry;
      solver_float* dvalue;
      solver_float* in;

      lnid = recv_messenger->mapping[idx];
      in = (solver_float*)((char*)recv_messenger->indata +
                           recv_messenger->insize * idx);

      dvalue = (solver_float*)valuetable + recv_messenger->mapping_itemoffset[idx];
      if (itemsperentry == 3) {
        if (direction == CONTRIBUTION) {
          dvalue[0] += in[0];
          dvalue[1] += in[1];
          dvalue[2] += in[2];
        } else {
          dvalue[0] = in[0];
          dvalue[1] = in[1];
          dvalue[2] = in[2];
        }
      } else {
        for (entry = 0; entry < itemsperentry; entry++) {
          if (direction == CONTRIBUTION) {
            *(dvalue + entry) += *(in + entry);
          } else {
            *(dvalue + entry) = *(in + entry);
          }
        }
      }

#ifdef DEBUG
      {
        int64_t* gnidp;
        gnidp = (int64_t*)((char*)in + itemsperentry * sizeof(solver_float));

        if (*gnidp != Global.myMesh->nodeTable[lnid].gnid) {
          fprintf(stderr, "Thread %d: solver_init: gnids do not match\n",
                  Global.myID);
          MPI_Abort(MPI_COMM_WORLD, ERROR);
          exit(1);
        }
      }
#endif /* DEBUG */
    }
  }

  return;
}

static void solver_diag_reset(void) {
  memset(g_dnDiagLocal, 0, sizeof(g_dnDiagLocal));
  memset(g_waitRecvStageLocal, 0, sizeof(g_waitRecvStageLocal));
  memset(g_waitRecvStageMin, 0, sizeof(g_waitRecvStageMin));
  memset(g_waitRecvStageMax, 0, sizeof(g_waitRecvStageMax));
  memset(g_waitRecvStageAvg, 0, sizeof(g_waitRecvStageAvg));
  memset(g_stageCommLocal, 0, sizeof(g_stageCommLocal));
  memset(g_stageCommMin, 0, sizeof(g_stageCommMin));
  memset(g_stageCommMax, 0, sizeof(g_stageCommMax));
  memset(g_stageCommAvg, 0, sizeof(g_stageCommAvg));
  memset(g_stageTimeLocal, 0, sizeof(g_stageTimeLocal));
  memset(g_stageTimeMin, 0, sizeof(g_stageTimeMin));
  memset(g_stageTimeMax, 0, sizeof(g_stageTimeMax));
  memset(g_stageTimeAvg, 0, sizeof(g_stageTimeAvg));
  g_waitRecvStageMaxRank[0] = -1;
  g_waitRecvStageMaxRank[1] = -1;
  g_waitRecvStageMaxRank[2] = -1;
  g_waitRecvStageMaxRank[3] = -1;
  g_dnDiagReady = 0;
}

static void solver_diag_prepare_local(mesh_t* mesh) {
  int32_t dnindex;
  double depsTotal = 0.0;

  if (mesh == NULL || mesh->dnodeTable == NULL) {
    g_dnDiagLocal[DN_DIAG_ADJUST_NODES] = 0.0;
    g_dnDiagLocal[DN_DIAG_ADJUST_DEPS] = 0.0;
    return;
  }

  g_dnDiagLocal[DN_DIAG_ADJUST_NODES] = (double)mesh->ldnnum;
  for (dnindex = 0; dnindex < mesh->ldnnum; dnindex++) {
    depsTotal += (double)mesh->dnodeTable[dnindex].deps;
  }
  g_dnDiagLocal[DN_DIAG_ADJUST_DEPS] = depsTotal;
}

static void solver_diag_reduce(void) {
  double* gatheredStage = NULL;

  MPI_Reduce(g_dnDiagLocal, g_dnDiagMin, DN_DIAG_METRIC_COUNT, MPI_DOUBLE,
             MPI_MIN, 0, comm_solver);
  MPI_Reduce(g_dnDiagLocal, g_dnDiagMax, DN_DIAG_METRIC_COUNT, MPI_DOUBLE,
             MPI_MAX, 0, comm_solver);
  MPI_Reduce(g_dnDiagLocal, g_dnDiagAvg, DN_DIAG_METRIC_COUNT, MPI_DOUBLE,
             MPI_SUM, 0, comm_solver);

  MPI_Reduce(g_waitRecvStageLocal, g_waitRecvStageMin, WAIT_RECV_STAGE_COUNT,
             MPI_DOUBLE, MPI_MIN, 0, comm_solver);
  MPI_Reduce(g_waitRecvStageLocal, g_waitRecvStageMax, WAIT_RECV_STAGE_COUNT,
             MPI_DOUBLE, MPI_MAX, 0, comm_solver);
  MPI_Reduce(g_waitRecvStageLocal, g_waitRecvStageAvg, WAIT_RECV_STAGE_COUNT,
             MPI_DOUBLE, MPI_SUM, 0, comm_solver);

  MPI_Reduce(g_stageCommLocal, g_stageCommMin, STAGE_COMM_METRIC_COUNT,
             MPI_DOUBLE, MPI_MIN, 0, comm_solver);
  MPI_Reduce(g_stageCommLocal, g_stageCommMax, STAGE_COMM_METRIC_COUNT,
             MPI_DOUBLE, MPI_MAX, 0, comm_solver);
  MPI_Reduce(g_stageCommLocal, g_stageCommAvg, STAGE_COMM_METRIC_COUNT,
             MPI_DOUBLE, MPI_SUM, 0, comm_solver);

  MPI_Reduce(g_stageTimeLocal, g_stageTimeMin, STAGE_TIME_METRIC_COUNT,
             MPI_DOUBLE, MPI_MIN, 0, comm_solver);
  MPI_Reduce(g_stageTimeLocal, g_stageTimeMax, STAGE_TIME_METRIC_COUNT,
             MPI_DOUBLE, MPI_MAX, 0, comm_solver);
  MPI_Reduce(g_stageTimeLocal, g_stageTimeAvg, STAGE_TIME_METRIC_COUNT,
             MPI_DOUBLE, MPI_SUM, 0, comm_solver);

  if (Global.myID == 0) {
    gatheredStage =
        (double*)malloc(sizeof(double) * WAIT_RECV_STAGE_COUNT *
                        (size_t)Global.theGroupSize);
  }

  MPI_Gather(g_waitRecvStageLocal, WAIT_RECV_STAGE_COUNT, MPI_DOUBLE,
             gatheredStage, WAIT_RECV_STAGE_COUNT, MPI_DOUBLE, 0, comm_solver);

  if (Global.myID == 0) {
    int i;
    for (i = 0; i < DN_DIAG_METRIC_COUNT; i++) {
      g_dnDiagAvg[i] /= (double)Global.theGroupSize;
    }

    for (i = 0; i < WAIT_RECV_STAGE_COUNT; i++) {
      int rank;
      double maxVal = -1.0;
      g_waitRecvStageAvg[i] /= (double)Global.theGroupSize;
      g_waitRecvStageMaxRank[i] = 0;
      for (rank = 0; rank < Global.theGroupSize; rank++) {
        double v = gatheredStage[(size_t)rank * WAIT_RECV_STAGE_COUNT + i];
        if (v > maxVal) {
          maxVal = v;
          g_waitRecvStageMaxRank[i] = rank;
        }
      }
    }

    for (i = 0; i < STAGE_COMM_METRIC_COUNT; i++) {
      g_stageCommAvg[i] /= (double)Global.theGroupSize;
    }

    for (i = 0; i < STAGE_TIME_METRIC_COUNT; i++) {
      g_stageTimeAvg[i] /= (double)Global.theGroupSize;
    }

    free(gatheredStage);
    g_dnDiagReady = 1;
  }
}

static int32_t stage_comm_index(int32_t stage, int32_t metric) {
  return (stage * 6) + metric;
}

static int32_t stage_time_index(int32_t stage, int32_t metric) {
  return (stage * STAGE_TIME_METRIC_COUNT_PER_STAGE) + metric;
}

static int32_t wait_recv_stage_from_msgtag(int32_t msgtag) {
  switch (msgtag) {
    case DN_FORCE_MSG:
      return WAIT_RECV_STAGE_1_DN_FORCE;
    case AN_FORCE_MSG:
      return WAIT_RECV_STAGE_2_AN_FORCE;
    case AN_DISP_MSG:
      return WAIT_RECV_STAGE_3_AN_DISP;
    case DN_DISP_MSG:
      return WAIT_RECV_STAGE_4_DN_DISP;
    default:
      return -1;
  }
}

/**
 * Print messenger entry data to the given output stream.
 * The output line contains the following information:
 * "pe_id type c/s rpe nodecount outsize insize"
 *
 * \note Printing is not synchronized accross PEs, so lines can be out of
 * order / interleaved
 */
static int messenger_print(messenger_t* msg, char type, char cs, FILE* out) {
  int ret;
  assert(NULL != msg);

  ret = fprintf(out, "msg_info: %d %c %c %d %d %d %d %d\n", Global.myID, type,
                cs, msg->procid, msg->nodecount, msg->outsize, msg->insize,
                msg->nodecount * msg->outsize);

  if (ret > 0 || !Param.theSchedulePrintErrorCheckFlag) {
    ret = 0;
  } else { /* ret <= 0 && Param.theSchedulePrintErrorCheckFlag */
    perror("Warning! fprintf(...) failed");
    fprintf(stderr, "PE id = %d, ret = %d, out = 0x%p\n", Global.myID, ret,
            out);
    fflush(stderr);
  }

  return ret;
}

/**
 * Print the messenger list for a scheduler
 *
 * \note Printing is not synchronized accross PEs, so lines can be out of
 * order / interleaved
 */
static int schedule_print_messenger_list(schedule_t* sched, messenger_t* msg,
                                         int count, char type, char cs,
                                         FILE* out) {
  messenger_t* m_p;
  int my_count;
  int ret = 0;

  assert(NULL != sched);

  m_p = msg;

  for (my_count = 0; m_p != NULL && ret >= 0; my_count++, m_p = m_p->next) {
    ret += messenger_print(msg, type, cs, out);
  }

  if (ret < 0) {
    return -1;
  }

  if (my_count != count) {
    fprintf(stderr,
            "Warning! schedule list count does not match: "
            "%u %c %c %d %d\n",
            Global.myID, type, cs, count, my_count);
  }

  return ret;
}

/**
 * Print scheduler communication to a given output stream.
 *
 * Print a line per entry in the schedule_t structure, i.e., an line
 * per messenger_t entry in the list.
 *
 * Each line has the following format:
 * "pe_id type c/s rpe nodecount outsize insize"
 */
static int schedule_print_detail(schedule_t* sched, char type, FILE* out) {
  int ret = 0;

  ret += schedule_print_messenger_list(sched, sched->first_c, sched->c_count,
                                       type, 'c', out);
  ret += schedule_print_messenger_list(sched, sched->first_s, sched->s_count,
                                       type, 's', out);

  return ret;
}

/**
 * Print scheduler communication to a given output stream.
 *
 * Print a line per entry in the schedule_t structure, i.e., an line
 * per messenger_t entry in the list.
 *
 * Each line has the following format:
 * "pe_id type s_count c_count"
 *
 * \note Printing is not synchronized accross PEs, so lines can be out of
 * order / interleaved
 */
static int schedule_print(schedule_t* sched, char type, FILE* out) {
  int ret;

  assert(NULL != sched);
  assert(NULL != out);

  ret = fprintf(out, "sch_info: %u %c %d %d\n", Global.myID, type,
                sched->c_count, sched->s_count);

  if (ret > 0 || !Param.theSchedulePrintErrorCheckFlag) {
    ret = 0;
  }

  return ret;
}

/**
 * \note Printing is not synchronized accross PEs, so lines can be out of
 * order / interleaved
 */
static int solver_print_schedules_imp(mysolver_t* solver, FILE* out) {
  int ret = 0;

  assert(solver != NULL);
  assert(out != NULL);

  MPI_Barrier(comm_solver);

  /* print the high-level schedule per PE */
  if (Global.myID == 0) { /* print some header information */
    fputs(
        "# ----------------------------------------------------------\n"
        "# Content: Solver schedule information\n"
        "# pe_id d/a s_count c_count\n",
        out);
  }
  fflush(out);
  fdatasync(fileno(out));

  MPI_Barrier(comm_solver);

  /* this is not synchronized, so it might be out of order */
  ret += schedule_print(solver->an_sched, 'a', out);
  ret += schedule_print(solver->dn_sched, 'd', out);

  fflush(out);
  fdatasync(fileno(out));
  MPI_Barrier(comm_solver);

  /* print the schedule detail */
  if (Global.myID == 0) { /* print some header information */
    fputs(
        "\n\n"
        "# ----------------------------------------------------------\n"
        "# Content: Solver schedule detail\n"
        "# pe_id d/a c/s rpe nodecount outsize insize msgsize\n",
        out);
    fflush(out);
    fdatasync(fileno(out));
  }

  MPI_Barrier(comm_solver);

  /* this is not synchronized, so it might be out of order */
  ret += schedule_print_detail(solver->an_sched, 'a', out);
  ret += schedule_print_detail(solver->dn_sched, 'd', out);

  fflush(out);
  fdatasync(fileno(out));
  MPI_Barrier(comm_solver);

  if (Global.myID == 0) {
    fputs(
        "# ----------------------------------------------------------\n"
        "\n",
        out);
    fflush(out);
    fdatasync(fileno(out));
  }

  MPI_Barrier(comm_solver);

  return ret;
}

/**
 * Wrapper to print the solver schedules to stdout and a file with
 * the given name.
 *
 * \note Printing is not synchronized accross PEs, so lines can be out of
 * order / interleaved
 */
static int solver_print_schedules(mysolver_t* solver) {
  FILE* out;
  int ret = 0;

  if (Param.theSchedulePrintToStdout) {
    /* print schedules to standard output */
    ret += solver_print_schedules_imp(solver, stdout);

    if (ret < 0) {
      fprintf(stderr,
              "Warning! printing schedules to standard output "
              "failed for PE %d\n",
              Global.myID);
    }
  }

  if (Param.theSchedulePrintToFile &&
      (Param.theSchedulePrintFilename != NULL)) {
    /* print schedules to the given file */
    out = fopen(Param.theSchedulePrintFilename, "a");

    if (NULL == out) { /* this is not fatal */
      fprintf(stderr,
              "Warning!, PE# %d failed to open output file for "
              "printing the communication schedule\n",
              Global.myID);
      return -1;
    }

    ret = solver_print_schedules_imp(solver, out);

    if (ret < 0) {
      fprintf(stderr,
              "Warning! PE %d could not print schedules "
              "to file\n",
              Global.myID);
    }

    fclose(out);
  }

  return ret;
}

/**
 * messenger_new: Allocate and initialize a messenger.
 *
 */
static messenger_t* messenger_new(int32_t procid) {
  messenger_t* messenger;

  messenger = (messenger_t*)calloc(1, sizeof(messenger_t));

  if (messenger == NULL) return NULL;

  messenger->procid = procid;

  return messenger;
}

/**
 * messenger_delete: Release memory used by a messenger.
 *
 */
static void messenger_delete(messenger_t* messenger) {
  if (messenger == NULL) return;

  if (messenger->outdata != NULL) free(messenger->outdata);

  if (messenger->indata != NULL) free(messenger->indata);

  if (messenger->mapping_itemoffset != NULL) free(messenger->mapping_itemoffset);

  free(messenger->mapping);

  free(messenger);

  return;
}

/**
 * messenger_set: Free any data memory used by the messenger for
 *                previous communication. And allocate new memory
 *                for the new round of communication.
 *
 */
static void messenger_set(messenger_t* messenger, int32_t outsize,
                          int32_t insize) {
  if (messenger->outdata != NULL) {
    free(messenger->outdata);
    messenger->outdata = NULL;
  }

  if (messenger->indata != NULL) {
    free(messenger->indata);
    messenger->indata = NULL;
  }

  messenger->outsize = outsize;
  messenger->insize = insize;

  if (outsize != 0) {
    messenger->outdata = calloc(messenger->nodecount, outsize);
    if (messenger->outdata == NULL) {
      fprintf(stderr, "Thread %d: messenger_set: out of memory\n", Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
  }

  if (insize != 0) {
    messenger->indata = calloc(messenger->nodecount, insize);
    if (messenger->indata == NULL) {
      fprintf(stderr, "Thread %d: messenger_set: out of memory\n", Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
  }

  return;
}

static void messenger_prepare_offsets(messenger_t* first,
                                      int32_t itemsperentry) {
  messenger_t* messenger = first;
  while (messenger != NULL) {
    if ((messenger->mapping_itemoffset == NULL) ||
        (messenger->mapping_itemsperentry != itemsperentry)) {
      int32_t idx;
      if (messenger->mapping_itemoffset != NULL) {
        free(messenger->mapping_itemoffset);
        messenger->mapping_itemoffset = NULL;
      }

      if (messenger->nodecount > 0) {
        messenger->mapping_itemoffset =
            (int32_t*)malloc(sizeof(int32_t) * messenger->nodecount);
        if (messenger->mapping_itemoffset == NULL) {
          fprintf(stderr, "Thread %d: messenger_prepare_offsets: out of memory\n",
                  Global.myID);
          MPI_Abort(MPI_COMM_WORLD, ERROR);
          exit(1);
        }

        for (idx = 0; idx < messenger->nodecount; idx++) {
          messenger->mapping_itemoffset[idx] =
              messenger->mapping[idx] * itemsperentry;
        }
      }

      messenger->mapping_itemsperentry = itemsperentry;
      messenger->mapping_has_duplicates =
          messenger_detect_duplicate_offsets(messenger);
      if (messenger->mapping_has_duplicates &&
          (messenger->mapping_dup_reported == 0) && (Global.myID == 0)) {
        monitor_print(
            "Messenger duplicate mapping detected for peer %d; forcing "
            "serial unpack for safety\n",
            messenger->procid);
        messenger->mapping_dup_reported = 1;
      }
    }

    messenger = messenger->next;
  }
}

static int compare_int32_asc(const void* a, const void* b) {
  int32_t va = *(const int32_t*)a;
  int32_t vb = *(const int32_t*)b;
  if (va < vb) {
    return -1;
  }
  if (va > vb) {
    return 1;
  }
  return 0;
}

static int32_t messenger_detect_duplicate_offsets(messenger_t* messenger) {
  int32_t count;
  int32_t* copy;

  if ((messenger == NULL) || (messenger->mapping_itemoffset == NULL)) {
    return 0;
  }

  count = messenger->nodecount;
  if (count <= 1) {
    return 0;
  }

  copy = (int32_t*)malloc(sizeof(int32_t) * count);
  if (copy == NULL) {
    fprintf(stderr,
            "Thread %d: messenger_detect_duplicate_offsets: out of memory\n",
            Global.myID);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  memcpy(copy, messenger->mapping_itemoffset, sizeof(int32_t) * count);
  qsort(copy, count, sizeof(int32_t), compare_int32_asc);

  for (int32_t i = 1; i < count; i++) {
    if (copy[i] == copy[i - 1]) {
      free(copy);
      return 1;
    }
  }

  free(copy);
  return 0;
}

/**
 * messenger_countnodes: Count the total number of nodes (that need to
 *                       be communicated) on the messenger link list.
 *
 */
static int32_t messenger_countnodes(messenger_t* first) {
  messenger_t* messenger;
  int32_t total;

  total = 0;
  messenger = first;
  while (messenger != NULL) {
    total += messenger->nodecount;
    messenger = messenger->next;
  }

  return total;
}

/**
 * Compute K1, K2 and K3 for a  linear  elastic,  homogeneous
 * and isotropic cube; and their off-diagonal matrices.
 *
 *            K1 = Integral(   Del(Phi(j))  * T(Del(Phi(i)))     )
 *            K2 = Integral(   Del(Phi(i))  * T(Del(Phi(j)))     )
 *            K3 = Integral( T(Del(Phi(i))) *   Del(Phi(j))  * I )
 *
 *            where:
 *
 *            T(W) = transpose of W.
 *
 * \note K offdiagonal has been commented out because is no longer needed
 * after the changes made on the solution algorith due to the new
 * form to handle the damping involved terms.
 */
static void compute_K() {
  int i, j; /* indices of the block matrices, i for rows, j for columns */
  int k, l; /* indices of 3 x 3 matrices, k for rows, l for columns     */

  double x[3][8] = {{-1, 1, -1, 1, -1, 1, -1, 1},
                    {-1, -1, 1, 1, -1, -1, 1, 1},
                    {-1, -1, -1, -1, 1, 1, 1, 1}};

  /* compute K3 first */
  memset(Global.theK3, 0, 8 * 8 * sizeof(fmatrix_t));
  for (i = 0; i < 8; i++) {
    for (j = 0; j < 8; j++) {
      /* for  each 3 x 3 matrix representing (node i, node j),
             each of these matrices is diagonal, off-diagonal elements
             have been set to 0 in memset() */

      fmatrix_t* matPtr;
      matPtr = &Global.theK3[i][j];

      for (k = 0; k < 3; k++) {
        /* set the diagonal values for the 3 x 3 matrix.
           Note the rotation of the index */
        double I1, I2, I3;

        I1 =
            INTEGRAL_1(x[(k + 0) % 3][i], x[(k + 0) % 3][j], x[(k + 1) % 3][i],
                       x[(k + 1) % 3][j], x[(k + 2) % 3][i], x[(k + 2) % 3][j]);

        I2 =
            INTEGRAL_1(x[(k + 1) % 3][i], x[(k + 1) % 3][j], x[(k + 2) % 3][i],
                       x[(k + 2) % 3][j], x[(k + 0) % 3][i], x[(k + 0) % 3][j]);

        I3 =
            INTEGRAL_1(x[(k + 2) % 3][i], x[(k + 2) % 3][j], x[(k + 0) % 3][i],
                       x[(k + 0) % 3][j], x[(k + 1) % 3][i], x[(k + 1) % 3][j]);

        matPtr->f[k][k] = I1 + I2 + I3;
      } /* k*/
    } /* j*/
  } /* i */

  /* compute K1 and K2. They are not diagonal either, but they are
     indeed symmetric globablly */
  for (i = 0; i < 8; i++) {
    for (j = 0; j < 8; j++) {
      /* for  each 3 x 3 matrix representing (node i, node j)*/

      fmatrix_t *matPtr0, *matPtr1;
      matPtr0 = &Global.theK1[i][j];
      matPtr1 = &Global.theK2[i][j];

      for (k = 0; k < 3; k++)
        for (l = 0; l < 3; l++) {
          if (k == l) {
            /* Initialize the diagnoal elements */
            matPtr0->f[k][k] = INTEGRAL_1(x[(k + 0) % 3][i], x[(k + 0) % 3][j],
                                          x[(k + 1) % 3][i], x[(k + 1) % 3][j],
                                          x[(k + 2) % 3][i], x[(k + 2) % 3][j]);

            matPtr1->f[k][k] = INTEGRAL_1(x[(k + 0) % 3][j], x[(k + 0) % 3][i],
                                          x[(k + 1) % 3][j], x[(k + 1) % 3][i],
                                          x[(k + 2) % 3][j], x[(k + 2) % 3][i]);
          } else {
            /* Initialize off-diagonal elements */

            matPtr0->f[k][l] = INTEGRAL_2(x[k][j], x[l][i], x[3 - (k + l)][j],
                                          x[3 - (k + l)][i]);

            matPtr1->f[k][l] = INTEGRAL_2(x[k][i], x[l][j], x[3 - (k + l)][i],
                                          x[3 - (k + l)][j]);

          } /* else */
        } /* for l */
    } /* for j */
  } /* i */

  /* Old code for damping                                    */
  /* Create the off-diagonal matrices                        */
  /* memcpy(theK1_offdiag, Global.theK1, sizeof(double) * 24 * 24); */
  /* memcpy(theK2_offdiag, Global.theK2, sizeof(double) * 24 * 24); */
  /* memcpy(theK3_offdiag, Global.theK3, sizeof(double) * 24 * 24); */
  /* for (i = 0; i <8; i++) {                                */
  /*     theK1_offdiag[i][i].f[0][0] = 0;                    */
  /*     theK1_offdiag[i][i].f[1][1] = 0;                    */
  /*     theK1_offdiag[i][i].f[2][2] = 0;                    */
  /*     theK2_offdiag[i][i].f[0][0] = 0;                    */
  /*     theK2_offdiag[i][i].f[1][1] = 0;                    */
  /*     theK2_offdiag[i][i].f[2][2] = 0;                    */
  /*     theK3_offdiag[i][i].f[0][0] = 0;                    */
  /*     theK3_offdiag[i][i].f[1][1] = 0;                    */
  /*     theK3_offdiag[i][i].f[2][2] = 0;                    */
  /* }                                                       */

  /* New code:
   * First option to solve double K1 K3 computation in the
   * time steps is to merge both of them in K1 only
   */
  for (i = 0; i < 8; i++) {
    for (j = 0; j < 8; j++) {
      for (k = 0; k < 3; k++) {
        for (l = 0; l < 3; l++) {
          Global.theK1[i][j].f[k][l] += Global.theK3[i][j].f[k][l];
        }
      }
    }
  }

  return;
}

static void constract_Quality_Factor_Table() {
  int i, j;
  /* double local_QTABLE[26][6] = {
          {   5.00,  0.211111102,  0.236842104,  0.032142857,  0.271428571,
  0.1400000 }, {   6.25,  0.188888889,  0.184210526,  0.039893617,  0.336879433,
  0.1015200 }, {   8.33,  0.157777778,  0.139473684,  0.045000000,  0.380000000,
  0.0700000 }, {  10.00,  0.137777765,  0.121052630,  0.032942899,  0.278184480,
  0.0683000 }, {  15.00,  0.097777765,  0.081052630,  0.032942899,  0.278184480,
  0.0450000 }, {  20.00,  0.078139527,  0.060526314,  0.031409788,  0.277574872,
  0.0342250 }, {  25.00,  0.064285708,  0.049999999,  0.031578947,  0.285714286,
  0.0266000 }, {  30.00,  0.053658537,  0.044736842,  0.026640676,  0.246913580,
  0.0230850 }, {  35.00,  0.046341463,  0.038157895,  0.027098480,  0.251156642,
  0.0196690 }, {  40.00,  0.040487805,  0.034210526,  0.025949367,  0.240506329,
  0.0173800 }, {  45.00,  0.036585366,  0.028947368,  0.031393568,  0.290964778,
  0.0143660 }, {  50.00,  0.032926829,  0.026315789,  0.032488114,  0.301109350,
  0.0126200 }, {  60.00,  0.027900000,  0.022300000,  0.027500000,  0.254500000,
  0.0114000 }, {  70.00,  0.024000000,  0.019000000,  0.032488114,  0.301109350,
  0.0083000 }, {  80.00,  0.020700000,  0.017400000,  0.025100000,  0.232600000,
  0.0088000 }, {  90.00,  0.018700000,  0.015400000,  0.024400000,  0.225600000,
  0.0079000 }, { 100.00,  0.017000000,  0.014000000,  0.028021016,  0.288966725,
  0.0062810 }, { 120.00,  0.014200000,  0.011500000,  0.028000000,  0.270000000,
  0.0052000 }, { 150.00,  0.011400000,  0.009400000,  0.024000000,  0.231600000,
  0.0047000 }, { 200.00,  0.008500000,  0.007050000,  0.022603978,  0.226039783,
  0.0035392 }, { 250.00,  0.006900000,  0.005500000,  0.026900000,  0.259600000,
  0.0027000 }, { 300.00,  0.005700000,  0.004700000,  0.027072758,  0.279187817,
  0.0021276 }, { 350.00,  0.004800000,  0.004000000,  0.024200000,  0.233900000,
  0.0020000 }, { 400.00,  0.004300000,  0.003600000,  0.021425572,  0.214255718,
  0.0017935 }, { 450.00,  0.003900000,  0.003000000,  0.028000000,  0.271000000,
  0.0015000 }, { 500.00,  0.003500000,  0.002850000,  0.023408925,  0.241404535,
  0.0013670 }
          //    Qo,      alpha_1,      alpha_2,      gamma_1,      gamma_2, beta
  }; */
  // Doriam's table
  double local_QTABLE[26][6] = {
      {5.0, 0.213824199306997, 0.214975081555006, 0.051565915820727,
       0.392419601076068, 0.114223560617828},
      {6.25, 0.186829389245902, 0.176495993330104, 0.050021147523697,
       0.383250861967562, 0.091379110600293},
      {8.33, 0.153250864931585, 0.135813992388522, 0.048477940850105,
       0.374166449435453, 0.068561729073257},
      {10.00, 0.133541701119348, 0.114543936056555, 0.047705132221359,
       0.369647375406900, 0.057111817639174},
      {15.00, 0.095975244603152, 0.077914986423977, 0.046428268749996,
       0.362223775996977, 0.038074338344773},
      {20.00, 0.074734684302201, 0.059014449564976, 0.045797618339140,
       0.358572526832180, 0.028555662894623},
      {25.00, 0.061146200167406, 0.047488939930242, 0.045423352297652,
       0.356407244204920, 0.022844499869013},
      {30.00, 0.051722724340693, 0.039728560381090, 0.045175923644991,
       0.354975080940478, 0.019037079640288},
      {35.00, 0.044809339288316, 0.034147885877092, 0.045000204762312,
       0.353957211053003, 0.016317504582817},
      {40.00, 0.039523175150440, 0.029941793660133, 0.044868873696694,
       0.353195975615256, 0.014277828479929},
      {45.00, 0.035351165343402, 0.026658136972439, 0.044766892469425,
       0.352604653279352, 0.012691416011863},
      {50.00, 0.031975096967461, 0.024023455052625, 0.044685320595041,
       0.352131663633679, 0.011422286797518},
      {60.00, 0.026846337747289, 0.020058425934387, 0.044562752020458,
       0.351421365568300, 0.009518592500906},
      {70.00, 0.023134990683038, 0.017216647306878, 0.044474817264180,
       0.350912593773531, 0.008158808695107},
      {80.00, 0.020324985209151, 0.015080030053759, 0.044408507155927,
       0.350529798411748, 0.007138968333117},
      {90.00, 0.018123577678620, 0.013415083146655, 0.044356656285723,
       0.350231225564975, 0.006345756901192},
      {100.00, 0.016352382860205, 0.012081157194701, 0.044314985694202,
       0.349991873360687, 0.005711185876928},
      {120.00, 0.013678649106408, 0.010077015436344, 0.044252210323555,
       0.349632434464568, 0.004759325395322},
      {150.00, 0.010984355197906, 0.008069044748514, 0.044189407612328,
       0.349274010595048, 0.003807459231929},
      {200.00, 0.008269251319794, 0.006057391664144, 0.044127549161943,
       0.348920585136250, 0.002855587254772},
      {250.00, 0.006630158272276, 0.004848705344604, 0.044091695658935,
       0.348713672248562, 0.002284462652439},
      {300.00, 0.005533298472639, 0.004042204370751, 0.044068463407196,
       0.348577928429087, 0.001903713971116},
      {350.00, 0.004747858760699, 0.003465747986440, 0.044051784382370,
       0.348480186531356, 0.001631752444839},
      {400.00, 0.004157744966802, 0.003033152859528, 0.044038500649086,
       0.348403654469221, 0.001427783254767},
      {450.00, 0.003698184940471, 0.002696500100933, 0.044026822537802,
       0.348339033919198, 0.001269142412984},
      {500.00, 0.003330190674661, 0.002427023350810, 0.044015678669530,
       0.348280855263629, 0.001142231431869}};
  //    Qo,      alpha_1,      alpha_2,      gamma_1,      gamma_2,       beta

  for (i = 0; i < 26; i++) {
    for (j = 0; j < 6; j++) {
      Global.theQTABLE[i][j] = local_QTABLE[i][j];
    }
  }

  return;
}

/**
 * compute_setflag:
 *
 * - results from the discussion with Leo
 * - set flags as if in the full space.
 * - the main() routine will set the flag properly in case half-space
 *   is desired.
 *
 */
#ifdef BOUNDARY
static char compute_setflag(tick_t ldb[3], tick_t ruf[3], tick_t p1[3],
                            tick_t p2[3]) {
  char flag;

  flag = 13; /* indicate internal element */

  if (ldb[0] == p1[0]) flag = 12;
  if (ldb[1] == p1[1]) flag = 10;
  if (ldb[2] == p1[2]) flag = 4;

  if (ruf[0] == p2[0]) flag = 14;
  if (ruf[1] == p2[1]) flag = 16;
  if (ruf[2] == p2[2]) flag = 22;

  if (ldb[0] == p1[0] && ldb[1] == p1[1]) flag = 9;

  if (ruf[0] == p2[0] && ldb[1] == p1[1]) flag = 11;

  if (ldb[0] == p1[0] && ruf[1] == p2[1]) flag = 15;

  if (ruf[0] == p2[0] && ruf[1] == p2[1]) flag = 17;

  if (ldb[0] == p1[0] && ldb[2] == p1[2]) flag = 3;

  if (ruf[0] == p2[0] && ldb[2] == p1[2]) flag = 5;

  if (ldb[0] == p1[0] && ruf[2] == p2[2]) flag = 21;

  if (ruf[0] == p2[0] && ruf[2] == p2[2]) flag = 23;

  if (ldb[1] == p1[1] && ldb[2] == p1[2]) flag = 1;

  if (ruf[1] == p2[1] && ldb[2] == p1[2]) flag = 7;

  if (ldb[1] == p1[1] && ruf[2] == p2[2]) flag = 19;

  if (ruf[1] == p2[1] && ruf[2] == p2[2]) flag = 25;

  if (ldb[0] == p1[0] && ldb[1] == p1[1] && ldb[2] == p1[2]) flag = 0;

  if (ruf[0] == p2[0] && (ldb[1] == p1[1]) && ldb[2] == p1[2]) flag = 2;

  if (ldb[0] == p1[0] && ruf[1] == p2[1] && ldb[2] == p1[2]) flag = 6;

  if (ruf[0] == p2[0] && ruf[1] == p2[1] && ldb[2] == p1[2]) flag = 8;

  if (ldb[0] == p1[0] && ldb[1] == p1[1] && ruf[2] == p2[2]) flag = 18;

  if (ruf[0] == p2[0] && ldb[1] == p1[1] && ruf[2] == p2[2]) flag = 20;

  if (ldb[0] == p1[0] && ruf[1] == p2[1] && ruf[2] == p2[2]) flag = 24;

  if (ruf[0] == p2[0] && ruf[1] == p2[1] && ruf[2] == p2[2]) flag = 26;

  return flag;
}

static const int theIDBoundaryMatrix[27][8] = {
    {7, 6, 5, 4, 3, 2, 1, 0},
    {6, 6, 4, 4, 2, 2, 0, 0},
    {6, 7, 4, 5, 2, 3, 0, 1},
    {5, 4, 5, 4, 1, 0, 1, 0},
    {4, 4, 4, 4, 0, 0, 0, 0},
    {4, 5, 4, 5, 0, 1, 0, 1},
    {5, 4, 7, 6, 1, 0, 3, 2},
    {4, 4, 6, 6, 0, 0, 2, 2},
    {4, 5, 6, 7, 0, 1, 2, 3},
    {3, 2, 1, 0, 3, 2, 1, 0},
    {2, 2, 0, 0, 2, 2, 0, 0},
    {2, 3, 0, 1, 2, 3, 0, 1},
    {1, 0, 1, 0, 1, 0, 1, 0},
    {0, 0, 0, 0, 0, 0, 0, 0}, /* 13: internal elements */
    {0, 1, 0, 1, 0, 1, 0, 1},
    {1, 0, 3, 2, 1, 0, 3, 2},
    {0, 0, 2, 2, 0, 0, 2, 2},
    {0, 1, 2, 3, 0, 1, 2, 3},
    {3, 2, 1, 0, 7, 6, 5, 4},
    {2, 2, 0, 0, 6, 6, 4, 4},
    {2, 3, 0, 1, 6, 7, 4, 5},
    {1, 0, 1, 0, 5, 4, 5, 4},
    {0, 0, 0, 0, 4, 4, 4, 4},
    {0, 1, 0, 1, 4, 5, 4, 5},
    {1, 0, 3, 2, 5, 4, 7, 6},
    {0, 0, 2, 2, 4, 4, 6, 6},
    {0, 1, 2, 3, 4, 5, 6, 7},
};

static void compute_setboundary(float size, float Vp, float Vs, float rho,
                                int flag, double dashpot[8][3]) {
  int whichNode;
  double scale;

  /* init the damping vector to all zeroes */
  memset(dashpot, 0, sizeof(double) * 8 * 3);

#ifdef HALFSPACE
  flag = (flag < 9) ? flag + 9 : flag;
#endif /* HALFSPACE */

  scale = rho * (size / 2) * (size / 2);

  for (whichNode = 0; whichNode < 8; whichNode++) {
    int bitmark, component;

    bitmark = theIDBoundaryMatrix[flag][whichNode];

    switch (bitmark) {
      case 0:
        break;
      case 7:
        /* Three contributing faces */
        dashpot[whichNode][0] = dashpot[whichNode][1] = dashpot[whichNode][2] =
            (Vp + 2 * Vs) * scale;
        break;
      case 3:
      case 5:
      case 6:
        /* Two contributing faces */
        for (component = 0; component < 3; component++)
          dashpot[whichNode][component] =
              (Vs + ((bitmark & (1 << component)) ? Vp : Vs)) * scale;
        break;
      case 1:
      case 2:
      case 4:
        /* One contributing face */
        for (component = 0; component < 3; component++)
          dashpot[whichNode][component] =
              ((bitmark & (1 << component)) ? Vp : Vs) * scale;
        break;
      default:
        fprintf(stderr, "SetBoundary: Unknown bitmark. Panic!\n");
        exit(1);
    }
  }

  return;
}

/* CPML face-mask: which of the 5 absorbing faces an element touches.
 * The mask is computed from the element centroid's distance to each face;
 * an element belongs to face j's PML zone iff that distance < pml_thickness_m.
 * The z_min face (top, free surface) is never PML.  Corner/edge elements
 * have multiple bits set. */
static uint8_t compute_pml_face_mask(tick_t ldb[3], tick_t ruf[3],
                                     tick_t nearendp[3], tick_t farendp[3],
                                     double ticksize, double pml_thickness_m) {
  double cx = 0.5 * ((double)ldb[0] + (double)ruf[0]) * ticksize;
  double cy = 0.5 * ((double)ldb[1] + (double)ruf[1]) * ticksize;
  double cz = 0.5 * ((double)ldb[2] + (double)ruf[2]) * ticksize;
  double xn = (double)nearendp[0] * ticksize;
  double yn = (double)nearendp[1] * ticksize;
  double xf = (double)farendp[0]  * ticksize;
  double yf = (double)farendp[1]  * ticksize;
  double zf = (double)farendp[2]  * ticksize;

  uint8_t m = 0;
  if (cx - xn < pml_thickness_m) m |= PML_FACE_X_MIN;
  if (xf - cx < pml_thickness_m) m |= PML_FACE_X_MAX;
  if (cy - yn < pml_thickness_m) m |= PML_FACE_Y_MIN;
  if (yf - cy < pml_thickness_m) m |= PML_FACE_Y_MAX;
  if (zf - cz < pml_thickness_m) m |= PML_FACE_Z_MAX;
  return m;
}

/* Compute per-axis CPML ADE recurrence coefficients for one element.
 * Output triples are indexed by stretching axis j in {x=0,y=1,z=2}.
 * For axes outside the mask, returns the no-PML identity:
 *   b_j = 1, a_j = 0, invKappa_j = 1  (memory variable stays at 0).
 *
 * The damping profile uses the local element Vp, so material heterogeneity
 * is handled per-element automatically.
 */
static void compute_pml_profile_coeffs(uint8_t face_mask, tick_t ldb[3],
                                       tick_t ruf[3], tick_t nearendp[3],
                                       tick_t farendp[3], double ticksize,
                                       double L, double Vp, double dt,
                                       int poly_order, double R,
                                       double kappa_max, double alpha_max,
                                       double mpml_cross_fraction,
                                       double alpha_floor_fraction,
                                       double a_out[3], double b_out[3],
                                       double invKappa_out[3],
                                       double d_out[3],
                                       double alpha_out[3],
                                       double kappa_out[3]) {
  const double c[3] = {
      0.5 * ((double)ldb[0] + (double)ruf[0]) * ticksize,
      0.5 * ((double)ldb[1] + (double)ruf[1]) * ticksize,
      0.5 * ((double)ldb[2] + (double)ruf[2]) * ticksize};
  const double bounds_lo[3] = {(double)nearendp[0] * ticksize,
                               (double)nearendp[1] * ticksize,
                               0.0 /* z_min unused, free surface */};
  const double bounds_hi[3] = {(double)farendp[0] * ticksize,
                               (double)farendp[1] * ticksize,
                               (double)farendp[2] * ticksize};
  const uint8_t lo_bits[3] = {PML_FACE_X_MIN, PML_FACE_Y_MIN, 0};
  const uint8_t hi_bits[3] = {PML_FACE_X_MAX, PML_FACE_Y_MAX, PML_FACE_Z_MAX};
  const int p = poly_order;
  const double d_max_pre = (L > 0.0)
                               ? (-(p + 1.0) * Vp * log(R) / (2.0 * L))
                               : 0.0;

  /* First pass: compute the *standard* per-axis CPML profile values
   * d_j, kappa_j, alpha_j. These are the K&M 2007 profiles unchanged.
   * alpha_floor_fraction (default 0.1) prevents alpha from collapsing to
   * zero at xi = L. */
  const double alpha_floor = alpha_floor_fraction * alpha_max;
  double d_raw[3] = {0.0, 0.0, 0.0};
  double poly_arr[3] = {0.0, 0.0, 0.0}; /* per-axis normalized profile (xi/L)^p */
  double alpha_arr[3] = {0.0, 0.0, 0.0};

  for (int j = 0; j < 3; j++) {
    double xi = 0.0;
    if (lo_bits[j] && (face_mask & lo_bits[j])) {
      double d_lo = (bounds_lo[j] + L) - c[j];
      if (d_lo > xi) xi = d_lo;
    }
    if (hi_bits[j] && (face_mask & hi_bits[j])) {
      double d_hi = c[j] - (bounds_hi[j] - L);
      if (d_hi > xi) xi = d_hi;
    }
    if (xi < 0.0) xi = 0.0;
    if (xi > L) xi = L;

    if (xi > 0.0) {
      double xi_norm = xi / L;
      double poly = pow(xi_norm, (double)p);
      d_raw[j] = d_max_pre * poly;
      poly_arr[j] = poly;
      alpha_arr[j] = alpha_max * (1.0 - xi_norm);
      if (alpha_arr[j] < alpha_floor) alpha_arr[j] = alpha_floor;
    }
  }

  /* Second pass: M-PML (Multi-axial PML) cross-coupling.
   *
   * Plain unsplit CPML for the elastic wave equation has a known long-time
   * instability mode at PML corners (K&M 2007 Section 5.4). The alpha
   * floor above prevents the immediate "alpha=0 at the outer edge"
   * variant of the failure, but a slower-growing creeping instability
   * remains, particularly when kappa_max > 1.
   *
   * M-PML (Meza-Fajardo & Papageorgiou 2008) cures this by replacing each
   * axis's damping rate d_j with an effective rate
   *     d_eff_j = d_j + p * (d_{k1} + d_{k2})
   * where p ∈ [0.05, 0.20] is a small cross-coupling fraction. The
   * per-axis psi recurrence then uses d_eff_j in the b and a coefficients,
   * while the constitutive operator (1/kappa_j) is unchanged. Effect:
   * every axis in a PML element gets at least *some* damping rate, even
   * axes that have no own face-mask bit, which restores the symmetry the
   * corner mode exploited.
   *
   * Axes that inherit damping only via cross-coupling have alpha=0 from
   * the first pass; we floor alpha for those axes too so the ADE
   * recurrence has a non-zero decay rate.
   *
   * mpml_cross_fraction is the runtime-tunable p value (default 0.1).
   * Set to 0 to disable M-PML and run plain CPML; values in [0.05, 0.20]
   * are the literature-recommended range. */
  const double mpml_p = mpml_cross_fraction;

  for (int j = 0; j < 3; j++) {
    int k1 = (j + 1) % 3;
    int k2 = (j + 2) % 3;
    double d_total = d_raw[j] + mpml_p * (d_raw[k1] + d_raw[k2]);
    /* Formal-multiaxial kappa (Meza-Fajardo & Papageorgiou 2008; the
     * same profile-cross-mixing the v3 "formal multiaxial PML" uses, and
     * what SPECFEM-class codes ship -- the full Pi=s_x s_y s_z double/
     * triple-delta cross-convolution form is theoretical and not used in
     * production). The normalized profile is cross-mixed exactly like
     * d_total, then kappa is derived from it, so at an edge/corner kappa_j
     * reflects ALL active damping directions (consistent with the
     * cross-coupled d_total and with the matched kappa-inertia). Clamp to
     * [0,1] keeps kappa in [1, kappa_max]. Reduces EXACTLY to the per-axis
     * profile when mpml_p = 0, kappa_max = 1, or on a face element
     * (transverse poly = 0) -> the validated face/kappa=1 results are
     * bit-identical. d_total and alpha_j keep their validated treatment. */
    double prof_mix = poly_arr[j] + mpml_p * (poly_arr[k1] + poly_arr[k2]);
    if (prof_mix > 1.0) prof_mix = 1.0;
    double kappa_j = 1.0 + (kappa_max - 1.0) * prof_mix;
    double alpha_j = alpha_arr[j];

    /* For axes inheriting damping only via M-PML cross-coupling
     * (d_raw[j]==0 but d_total>0), give them a non-zero alpha too. */
    if (d_total > 0.0 && alpha_j <= 0.0) {
      alpha_j = alpha_floor;
    }

    double decay = d_total / kappa_j + alpha_j;
    double b_j = exp(-decay * dt);
    double a_j = 0.0;
    double denom = d_total + kappa_j * alpha_j;
    if (denom > 1.0e-30 && d_total > 0.0) {
      a_j = (d_total / (kappa_j * denom)) * (b_j - 1.0);
    }
    a_out[j] = a_j;
    b_out[j] = b_j;
    invKappa_out[j] = 1.0 / kappa_j;
    /* kappa_j itself, exported so the caller can apply the matched
     * kappa-scaled inertia. The gradient side stretches the strain by
     * 1/kappa_j (invKappa_out); the second-order displacement CFS-CPML
     * requires the inertia to be scaled by the product of kappa over the
     * element's active stretch axes (rho*kappa*d2u/dt2, the leading term
     * of multiplying the momentum eq by s_x s_y s_z). Without it the
     * effective impedance is rho*c0/sqrt(kappa) != rho*c0 -> a spurious
     * reflection at the inner edge for kappa>1. kappa_j = 1 at the inner
     * edge (xi=0) so the matched pair reduces to the native operator. */
    kappa_out[j] = kappa_j;
    /* d_total (M-PML-coupled damping rate, units 1/s) is exported so the
     * caller can build the Matzen second-order divergence-side term:
     * the mass-proportional CPML damping rho * d_j * du/dt injected into
     * the nodal mass arrays (mirroring the Lysmer dashpot path). */
    d_out[j] = d_total;
    /* Floored alpha_a (CFS frequency, rad/s) for this axis. Exported so
     * the caller can build the CFS restoring (-d*alpha*u) and the
     * displacement-convolution memory (+d*alpha^2*zeta) terms. The floor
     * (alpha_floor) already applied above guarantees alpha_a > 0, so the
     * exact zeta integrator gain (1-exp(-alpha dt))/alpha is well-posed. */
    alpha_out[j] = alpha_j;
  }
}
#endif /* BOUNDARY */

/**
 * compute_setab: the base a and b values will be scaled by zeta
 *                specific to each element.
 */
static void compute_setab(double freq, double* aBasePtr, double* bBasePtr) {
  /* old version which caused overflow because of the aproximation in
   * the derivative */

  double w1, w2, lw1, lw2, sw1, sw2, cw1, cw2;
  double numer, denom;

  if (Param.theTypeOfDamping == RAYLEIGH) {
    /* the factors 0.2 and 1 were calibrated heuristically by LEO */
    w1 = 2 * PI * freq * .2;
    w2 = 2 * PI * freq * 1;

    /* logs */
    lw1 = log(w1);
    lw2 = log(w2);

    /* squares */
    sw1 = w1 * w1;
    sw2 = w2 * w2;

    /* cubes */
    cw1 = w1 * w1 * w1;
    cw2 = w2 * w2 * w2;

    /* numerator */
    numer =
        w1 * w2 *
        (-2 * sw1 * lw2 + 2 * sw1 * lw1 - 2 * w1 * w2 * lw2 +
         2 * w1 * w2 * lw1 + 3 * sw2 - 3 * sw1 - 2 * sw2 * lw2 + 2 * sw2 * lw1);

    /* denominator */
    denom = (cw1 - cw2 + 3 * sw2 * w1 - 3 * sw1 * w2);

    /* the a over zeta target is... */
    *aBasePtr = numer / denom;

    /* new numerator */
    numer = 3 * (2 * w1 * w2 * lw2 - 2 * w1 * w2 * lw1 + sw1 - sw2);

    /* the b over zeta target is... */
    *bBasePtr = numer / denom;
  } else if (Param.theTypeOfDamping == MASS) {
    w1 = 2 * PI * freq * .1; /* these .1 and 8 heuristics */
    w2 = 2 * PI * freq * 8;

    numer = 2 * w2 * w1 * log(w2 / w1);
    denom = w2 - w1;

    *aBasePtr = 1.3 * numer / denom; /* this 1.3 comes out from heuristics */
    *bBasePtr = 0;
  } else if (Param.theTypeOfDamping == NONE || Param.theTypeOfDamping >= BKT) {
    *aBasePtr = 0;
    *bBasePtr = 0;
  }

  return;
}

int is_nodeloaded(int32_t iNode, char* onoff) {
  /* \todo use a general bitmap data structure for this */
  /* \todo move the routine declaration to the source generation file */

  int32_t whichByte, whichBit;
  unsigned char mask, test;

  whichByte = iNode / 8;
  whichBit = 7 - iNode % 8;

  mask = (unsigned char)(1u << whichBit);

  test = ((unsigned char)onoff[whichByte]) & mask;

  return (test == mask); /* 1 if equal, 0 otherwise */
}

/**
 * Add the force due to earthquake source.
 *
 * Globals accessed:
 * - Global.mySolver->force (W)
 * - Param.theDeltaTSquared (R)
 * - Global.myForces (R)
 *
 * Iterates over the nodes that are loaded by the source and set the
 * respective forces for those nodes.
 */
static void compute_addforce_s(int32_t timestep) {
  int i; /* index for loaded nodes (from the source) */

  for (i = 0; i < Global.theNodesLoaded; i++) {
    int lnid = Global.theNodesLoadedList[i]; /* local node id */

    /* node's force vector */
    fvector_t* nodalForce = Global.mySolver->force + lnid;

    /* vector-scalar multiply */
    nodalForce->f[0] = (Global.myForces[i].x[0]) * Param.theDeltaTSquared;
    nodalForce->f[1] = (Global.myForces[i].x[1]) * Param.theDeltaTSquared;
    nodalForce->f[2] = (Global.myForces[i].x[2]) * Param.theDeltaTSquared;
  }
}

static void compute_adjust_cache_free(void) {
  free(g_adjustParentOffsets);
  free(g_adjustParentIds);
  free(g_adjustDanglingLnid);
  free(g_adjustParentBase3);
  free(g_adjustDanglingBase3);
  free(g_adjustInvDeps);
  free(g_adjustDnOrder);
  free(g_adjustDnSortKey);
  free(g_adjustDistEdgeDn);
  free(g_adjustDistEdgeParentBase3);
  free(g_adjustDistRunParentBase3);
  free(g_adjustDistRunOffsets);

  g_adjustParentOffsets = NULL;
  g_adjustParentIds = NULL;
  g_adjustDanglingLnid = NULL;
  g_adjustParentBase3 = NULL;
  g_adjustDanglingBase3 = NULL;
  g_adjustInvDeps = NULL;
  g_adjustDnOrder = NULL;
  g_adjustDnSortKey = NULL;
  g_adjustDistEdgeDn = NULL;
  g_adjustDistEdgeParentBase3 = NULL;
  g_adjustDistRunParentBase3 = NULL;
  g_adjustDistRunOffsets = NULL;
  g_adjustDistRunCount = 0;
  g_adjustDnCount = 0;
  g_adjustParentCount = 0;
}

static int compare_adjust_dn_order(const void* a, const void* b) {
  int32_t ia = *(const int32_t*)a;
  int32_t ib = *(const int32_t*)b;
  int32_t ka = g_adjustDnSortKey[ia];
  int32_t kb = g_adjustDnSortKey[ib];
  if (ka < kb) {
    return -1;
  }
  if (ka > kb) {
    return 1;
  }
  if (ia < ib) {
    return -1;
  }
  if (ia > ib) {
    return 1;
  }
  return 0;
}

static int compare_adjust_dist_edge(const void* a, const void* b) {
  int32_t ia = *(const int32_t*)a;
  int32_t ib = *(const int32_t*)b;
  int32_t ka = g_adjustDistEdgeParentBase3[ia];
  int32_t kb = g_adjustDistEdgeParentBase3[ib];
  if (ka < kb) {
    return -1;
  }
  if (ka > kb) {
    return 1;
  }
  if (ia < ib) {
    return -1;
  }
  if (ia > ib) {
    return 1;
  }
  return 0;
}

static void compute_adjust_cache_build(void) {
  int32_t dnindex;
  int32_t totalParents = 0;
  int32_t cursor = 0;

  if ((g_adjustParentOffsets != NULL) &&
      (g_adjustDnCount == Global.myMesh->ldnnum)) {
    return;
  }

  compute_adjust_cache_free();

  g_adjustDnCount = Global.myMesh->ldnnum;
  if (g_adjustDnCount <= 0) {
    return;
  }

  g_adjustParentOffsets =
      (int32_t*)malloc(sizeof(int32_t) * (g_adjustDnCount + 1));
  g_adjustDanglingLnid = (int32_t*)malloc(sizeof(int32_t) * g_adjustDnCount);
  g_adjustDanglingBase3 =
      (int32_t*)malloc(sizeof(int32_t) * g_adjustDnCount);
  g_adjustInvDeps = (solver_float*)malloc(sizeof(solver_float) * g_adjustDnCount);
  g_adjustDnOrder = (int32_t*)malloc(sizeof(int32_t) * g_adjustDnCount);
  g_adjustDnSortKey = (int32_t*)malloc(sizeof(int32_t) * g_adjustDnCount);

  if ((g_adjustParentOffsets == NULL) || (g_adjustDanglingLnid == NULL) ||
      (g_adjustDanglingBase3 == NULL) || (g_adjustInvDeps == NULL) ||
      (g_adjustDnOrder == NULL) || (g_adjustDnSortKey == NULL)) {
    solver_abort(__FUNCTION_NAME, NULL,
                 "Failed allocating compute_adjust cache arrays");
  }

  g_adjustParentOffsets[0] = 0;
  for (dnindex = 0; dnindex < g_adjustDnCount; dnindex++) {
    dnode_t* dnode = &Global.myMesh->dnodeTable[dnindex];
    int32_t deps = (int32_t)dnode->deps;
    if (deps < 0) {
      solver_abort(__FUNCTION_NAME, NULL, "Invalid deps count for dnode %d",
                   dnindex);
    }
    g_adjustDanglingLnid[dnindex] = dnode->ldnid;
    g_adjustDanglingBase3[dnindex] = dnode->ldnid * 3;
    g_adjustInvDeps[dnindex] = (deps > 0) ? ((solver_float)1.0 / deps) : 0;
    g_adjustDnOrder[dnindex] = dnindex;
    totalParents += deps;
    g_adjustParentOffsets[dnindex + 1] = totalParents;
  }

  g_adjustParentIds = (int32_t*)malloc(sizeof(int32_t) * totalParents);
  g_adjustParentBase3 = (int32_t*)malloc(sizeof(int32_t) * totalParents);
  g_adjustDistEdgeDn = (int32_t*)malloc(sizeof(int32_t) * totalParents);
  g_adjustDistEdgeParentBase3 = (int32_t*)malloc(sizeof(int32_t) * totalParents);
  if ((totalParents > 0) &&
      ((g_adjustParentIds == NULL) || (g_adjustParentBase3 == NULL) ||
       (g_adjustDistEdgeDn == NULL) ||
       (g_adjustDistEdgeParentBase3 == NULL))) {
    solver_abort(__FUNCTION_NAME, NULL,
                 "Failed allocating compute_adjust parent id array");
  }

  for (dnindex = 0; dnindex < g_adjustDnCount; dnindex++) {
    int32link_t* int32link = Global.myMesh->dnodeTable[dnindex].lanid;
    while (int32link != NULL) {
      if (cursor >= totalParents) {
        solver_abort(__FUNCTION_NAME, NULL,
                     "compute_adjust cache overflow while building parent ids");
      }
      g_adjustParentIds[cursor++] = int32link->id;
      int32link = int32link->next;
    }
  }

  if (cursor != totalParents) {
    solver_abort(__FUNCTION_NAME, NULL,
                 "compute_adjust cache mismatch: cursor=%d total=%d", cursor,
                 totalParents);
  }

  for (cursor = 0; cursor < totalParents; cursor++) {
    g_adjustParentBase3[cursor] = g_adjustParentIds[cursor] * 3;
    g_adjustDistEdgeParentBase3[cursor] = g_adjustParentBase3[cursor];
  }

  for (dnindex = 0; dnindex < g_adjustDnCount; dnindex++) {
    int32_t start = g_adjustParentOffsets[dnindex];
    int32_t end = g_adjustParentOffsets[dnindex + 1];
    int32_t i, j;
    for (i = start + 1; i < end; i++) {
      int32_t key = g_adjustParentBase3[i];
      j = i - 1;
      while ((j >= start) && (g_adjustParentBase3[j] > key)) {
        g_adjustParentBase3[j + 1] = g_adjustParentBase3[j];
        j--;
      }
      g_adjustParentBase3[j + 1] = key;
    }

    if (start < end) {
      g_adjustDnSortKey[dnindex] = g_adjustParentBase3[start];
    } else {
      g_adjustDnSortKey[dnindex] = g_adjustDanglingBase3[dnindex];
    }
  }

  qsort(g_adjustDnOrder, g_adjustDnCount, sizeof(int32_t),
        compare_adjust_dn_order);

  if (totalParents > 0) {
    int32_t edgeCursor = 0;
    int32_t runCount = 0;

    for (dnindex = 0; dnindex < g_adjustDnCount; dnindex++) {
      int32_t start = g_adjustParentOffsets[dnindex];
      int32_t end = g_adjustParentOffsets[dnindex + 1];
      for (int32_t i = start; i < end; i++) {
        g_adjustDistEdgeDn[edgeCursor++] = dnindex;
      }
    }

    if (edgeCursor != totalParents) {
      solver_abort(__FUNCTION_NAME, NULL,
                   "compute_adjust cache mismatch: edge cursor=%d total=%d",
                   edgeCursor, totalParents);
    }

    {
      int32_t* edgeOrder = (int32_t*)malloc(sizeof(int32_t) * totalParents);
      int32_t* sortedEdgeDn = (int32_t*)malloc(sizeof(int32_t) * totalParents);
      int32_t* sortedEdgeParent =
          (int32_t*)malloc(sizeof(int32_t) * totalParents);
      if ((edgeOrder == NULL) || (sortedEdgeDn == NULL) ||
          (sortedEdgeParent == NULL)) {
        free(edgeOrder);
        free(sortedEdgeDn);
        free(sortedEdgeParent);
        solver_abort(__FUNCTION_NAME, NULL,
                     "Failed allocating compute_adjust distribution edge sort"
                     " arrays");
      }

      for (int32_t i = 0; i < totalParents; i++) {
        edgeOrder[i] = i;
      }

      qsort(edgeOrder, totalParents, sizeof(int32_t), compare_adjust_dist_edge);

      for (int32_t i = 0; i < totalParents; i++) {
        int32_t src = edgeOrder[i];
        sortedEdgeDn[i] = g_adjustDistEdgeDn[src];
        sortedEdgeParent[i] = g_adjustDistEdgeParentBase3[src];
      }

      memcpy(g_adjustDistEdgeDn, sortedEdgeDn, sizeof(int32_t) * totalParents);
      memcpy(g_adjustDistEdgeParentBase3, sortedEdgeParent,
             sizeof(int32_t) * totalParents);

      free(edgeOrder);
      free(sortedEdgeDn);
      free(sortedEdgeParent);
    }

    runCount = 0;
    for (int32_t i = 0; i < totalParents; i++) {
      if ((i == 0) ||
          (g_adjustDistEdgeParentBase3[i] != g_adjustDistEdgeParentBase3[i - 1])) {
        runCount++;
      }
    }

    g_adjustDistRunParentBase3 = (int32_t*)malloc(sizeof(int32_t) * runCount);
    g_adjustDistRunOffsets = (int32_t*)malloc(sizeof(int32_t) * (runCount + 1));
    if ((g_adjustDistRunParentBase3 == NULL) || (g_adjustDistRunOffsets == NULL)) {
      solver_abort(__FUNCTION_NAME, NULL,
                   "Failed allocating compute_adjust distribution run arrays");
    }

    {
      int32_t run = 0;
      g_adjustDistRunOffsets[0] = 0;
      for (int32_t i = 0; i < totalParents; i++) {
        if ((i == 0) ||
            (g_adjustDistEdgeParentBase3[i] !=
             g_adjustDistEdgeParentBase3[i - 1])) {
          if (run > 0) {
            g_adjustDistRunOffsets[run] = i;
          }
          g_adjustDistRunParentBase3[run] = g_adjustDistEdgeParentBase3[i];
          run++;
        }
      }
      g_adjustDistRunOffsets[run] = totalParents;
      g_adjustDistRunCount = run;
    }
  }

  g_adjustParentCount = totalParents;
}

/**
 * compute_adjust: Either distribute the values from LOCAL dangling nodes
 *                 to LOCAL anchored nodes, or assign values from LOCAL
 *                 anchored nodes to LOCAL dangling nodes.
 *
 */
static void compute_adjust(void* valuetable, int32_t itemsperentry,
                           int32_t how) {
  solver_float* vtable = (solver_float*)valuetable;
  int32_t dnindex;

  if (itemsperentry == 3) {
    Timer_Start("Adjust Cache Build");
    compute_adjust_cache_build();
    Timer_Stop("Adjust Cache Build");

    if (how == DISTRIBUTION) {
      Timer_Start("Adjust Distribution Core");
#ifdef _OPENMP
      const int32_t useOmp = adjust_distribution_use_omp();
#pragma omp parallel for if(useOmp) schedule(static)
#endif
      for (int32_t run = 0; run < g_adjustDistRunCount; run++) {
        int32_t start = g_adjustDistRunOffsets[run];
        int32_t end = g_adjustDistRunOffsets[run + 1];
        int32_t parentBase3 = g_adjustDistRunParentBase3[run];
        solver_float s0 = 0;
        solver_float s1 = 0;
        solver_float s2 = 0;

        for (int32_t e = start; e < end; e++) {
          int32_t d = g_adjustDistEdgeDn[e];
          solver_float* myvalue = vtable + g_adjustDanglingBase3[d];
          solver_float invdeps = g_adjustInvDeps[d];
          s0 += myvalue[0] * invdeps;
          s1 += myvalue[1] * invdeps;
          s2 += myvalue[2] * invdeps;
        }

        {
          solver_float* parentvalue = vtable + parentBase3;
          parentvalue[0] += s0;
          parentvalue[1] += s1;
          parentvalue[2] += s2;
        }
      }
      Timer_Stop("Adjust Distribution Core");
    } else {
      Timer_Start("Adjust Assignment Core");
      const int32_t useOmp = adjust_assignment_use_omp();
#ifdef _OPENMP
#pragma omp parallel for if(useOmp) schedule(static)
#endif
      for (int32_t d = 0; d < g_adjustDnCount; d++) {
        dnindex = g_adjustDnOrder[d];
        int32_t start = g_adjustParentOffsets[dnindex];
        int32_t end = g_adjustParentOffsets[dnindex + 1];
        int32_t deps = end - start;
        solver_float* myvalue = vtable + g_adjustDanglingBase3[dnindex];
        const solver_float invdeps = g_adjustInvDeps[dnindex];
        solver_float s0 = 0;
        solver_float s1 = 0;
        solver_float s2 = 0;

        if (deps <= 4) {
          solver_float* p0 = vtable + g_adjustParentBase3[start];
          s0 = p0[0];
          s1 = p0[1];
          s2 = p0[2];

          if (deps >= 2) {
            solver_float* p1 = vtable + g_adjustParentBase3[start + 1];
            s0 += p1[0];
            s1 += p1[1];
            s2 += p1[2];
          }
          if (deps >= 3) {
            solver_float* p2 = vtable + g_adjustParentBase3[start + 2];
            s0 += p2[0];
            s1 += p2[1];
            s2 += p2[2];
          }
          if (deps == 4) {
            solver_float* p3 = vtable + g_adjustParentBase3[start + 3];
            s0 += p3[0];
            s1 += p3[1];
            s2 += p3[2];
          }
        } else if (deps == 5) {
          solver_float* p0 = vtable + g_adjustParentBase3[start];
          solver_float* p1 = vtable + g_adjustParentBase3[start + 1];
          solver_float* p2 = vtable + g_adjustParentBase3[start + 2];
          solver_float* p3 = vtable + g_adjustParentBase3[start + 3];
          solver_float* p4 = vtable + g_adjustParentBase3[start + 4];
          s0 = p0[0] + p1[0] + p2[0] + p3[0] + p4[0];
          s1 = p0[1] + p1[1] + p2[1] + p3[1] + p4[1];
          s2 = p0[2] + p1[2] + p2[2] + p3[2] + p4[2];
        } else if (deps == 6) {
          solver_float* p0 = vtable + g_adjustParentBase3[start];
          solver_float* p1 = vtable + g_adjustParentBase3[start + 1];
          solver_float* p2 = vtable + g_adjustParentBase3[start + 2];
          solver_float* p3 = vtable + g_adjustParentBase3[start + 3];
          solver_float* p4 = vtable + g_adjustParentBase3[start + 4];
          solver_float* p5 = vtable + g_adjustParentBase3[start + 5];
          s0 = p0[0] + p1[0] + p2[0] + p3[0] + p4[0] + p5[0];
          s1 = p0[1] + p1[1] + p2[1] + p3[1] + p4[1] + p5[1];
          s2 = p0[2] + p1[2] + p2[2] + p3[2] + p4[2] + p5[2];
        } else {
          int32_t* parentBase3 = g_adjustParentBase3 + start;
          int32_t parentCount = end - start;
          for (int32_t p = 0; p < parentCount; p++) {
#if defined(__GNUC__) || defined(__clang__)
            if (p + 1 < parentCount) {
              __builtin_prefetch(vtable + parentBase3[p + 1], 0, 1);
            }
#endif
            solver_float* parentvalue = vtable + parentBase3[p];
            s0 += parentvalue[0];
            s1 += parentvalue[1];
            s2 += parentvalue[2];
          }
        }

        myvalue[0] = s0 * invdeps;
        myvalue[1] = s1 * invdeps;
        myvalue[2] = s2 * invdeps;
      }
      Timer_Stop("Adjust Assignment Core");
    }

    return;
  }

  if (how == DISTRIBUTION) {
    for (dnindex = 0; dnindex < Global.myMesh->ldnnum; dnindex++) {
      dnode_t* dnode;
      solver_float *myvalue, *parentvalue;
      solver_float darray[7]; /* A hack to avoid memory allocation */
      int32link_t* int32link;
      int32_t idx, parentlnid;
#ifdef DEBUG
      int32_t deps = 0;
#endif /* DEBUG */

      dnode = &Global.myMesh->dnodeTable[dnindex];
      myvalue = vtable + dnode->ldnid * itemsperentry;

      for (idx = 0; idx < itemsperentry; idx++) {
        darray[idx] = (*(myvalue + idx)) / dnode->deps;
      }

      /* Distribute my darray value to my anchors */
      int32link = dnode->lanid;
      while (int32link != NULL) {
#ifdef DEBUG
        deps++;
#endif

        parentlnid = int32link->id;
        parentvalue = vtable + parentlnid * itemsperentry;

        for (idx = 0; idx < itemsperentry; idx++) {
          /* Accumulation the distributed values */
          *(parentvalue + idx) += darray[idx];
        }

        int32link = int32link->next;
      }

#ifdef DEBUG
      if (deps != (int)dnode->deps) {
        fprintf(stderr, "Thread %d: compute_adjust distri: ", Global.myID);
        fprintf(stderr, "deps don't match\n");
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }
#endif /* DEBUG */
    } /* for all my LOCAL dangling nodes */
  } else {
    /* Assign the value of the anchored parents to the dangling nodes*/

    for (dnindex = 0; dnindex < Global.myMesh->ldnnum; dnindex++) {
      dnode_t* dnode;
      solver_float *myvalue, *parentvalue;
      int32link_t* int32link;
      int32_t idx, parentlnid;
#ifdef DEBUG
      int32_t deps = 0;
#endif /* DEBUG */

      dnode = &Global.myMesh->dnodeTable[dnindex];
      myvalue = vtable + dnode->ldnid * itemsperentry;

      /* Zero out the residual values the dangling node might
         still hold */
      memset(myvalue, 0, sizeof(solver_float) * itemsperentry);

      /* Assign prorated anchored values to a dangling node */
      int32link = dnode->lanid;
      while (int32link != NULL) {
#ifdef DEBUG
        deps++;
#endif

        parentlnid = int32link->id;
        parentvalue = vtable + parentlnid * itemsperentry;

        for (idx = 0; idx < itemsperentry; idx++) {
          *(myvalue + idx) += (*(parentvalue + idx) / dnode->deps);
        }

        int32link = int32link->next;
      }

#ifdef DEBUG
      if (deps != (int)dnode->deps) {
        fprintf(stderr, "Thread %d: compute_adjust assign: ", Global.myID);
        fprintf(stderr, "deps don't match\n");
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
      }
#endif /* DEBUG */

    } /* for all my LOCAL dangling nodes */
  }

  return;
}

static int32_t adjust_distribution_use_omp(void) {
#ifndef _OPENMP
  return 0;
#else
  if (hercules_repro_mode()) {
    return 0;
  }

  if (g_adjustDistributionOmpMode < 0) {
    const char* env = getenv("HERCULES_ADJUST_DIST_OMP");
    if (env != NULL) {
      if ((strcmp(env, "0") == 0) || (strcasecmp(env, "no") == 0) ||
          (strcasecmp(env, "false") == 0)) {
        g_adjustDistributionOmpMode = 0;
      } else {
        g_adjustDistributionOmpMode = 1;
      }
    } else {
      g_adjustDistributionOmpMode = 1;
    }

    if ((g_adjustDistributionOmpMode == 1) &&
        (Global.myID == 0) &&
        (g_adjustDistributionOmpReported == 0)) {
      monitor_print(
          "compute_adjust distribution OpenMP: enabled (threads=%d)\n",
          omp_get_max_threads());
      g_adjustDistributionOmpReported = 1;
    }
  }

  return g_adjustDistributionOmpMode;
#endif
}

static int32_t adjust_assignment_use_omp(void) {
#ifndef _OPENMP
  return 0;
#else
  if (hercules_repro_mode()) {
    return 0;
  }

  if (g_adjustAssignmentOmpMode < 0) {
    const char* env = getenv("HERCULES_ADJUST_ASSIGN_OMP");
    if (env != NULL) {
      if ((strcmp(env, "0") == 0) || (strcasecmp(env, "no") == 0) ||
          (strcasecmp(env, "false") == 0)) {
        g_adjustAssignmentOmpMode = 0;
      } else {
        g_adjustAssignmentOmpMode = 1;
      }
    } else {
      g_adjustAssignmentOmpMode = 1;
    }

    if ((g_adjustAssignmentOmpMode == 1) &&
        (Global.myID == 0) &&
        (g_adjustAssignmentOmpReported == 0)) {
      monitor_print(
          "compute_adjust assignment OpenMP: enabled (threads=%d)\n",
          omp_get_max_threads());
      g_adjustAssignmentOmpReported = 1;
    }
  }

  return g_adjustAssignmentOmpMode;
#endif
}

static int32_t full_clear_use_omp(void) {
#ifndef _OPENMP
  return 0;
#else
  if (hercules_repro_mode()) {
    return 0;
  }

  if (g_fullClearOmpMode < 0) {
    const char* env = getenv("HERCULES_FULL_CLEAR_OMP");

    if (env != NULL) {
      if ((strcmp(env, "0") == 0) || (strcasecmp(env, "no") == 0) ||
          (strcasecmp(env, "false") == 0) || (strcasecmp(env, "off") == 0)) {
        g_fullClearOmpMode = 0;
      } else {
        g_fullClearOmpMode = 1;
      }
    } else {
      g_fullClearOmpMode = 1;
    }

    if ((g_fullClearOmpMode == 1) && (Global.myID == 0) &&
        (g_fullClearOmpReported == 0)) {
      monitor_print("full host force clear OpenMP: enabled (threads=%d)\n",
                    omp_get_max_threads());
      g_fullClearOmpReported = 1;
    }
  }

  return g_fullClearOmpMode;
#endif
}

static void solver_record_force_upload_source_event(mysolver_t* solver) {
  if (g_forceUploadSourceSafeEventInit == 0) {
    if (cudaEventCreateWithFlags(&g_forceUploadSourceSafeEvent,
                                 cudaEventDisableTiming) != cudaSuccess) {
      fprintf(stderr,
              "Thread %d: force upload source-safe event create failed\n",
              Global.myID);
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }
    g_forceUploadSourceSafeEventInit = 1;
  }

  cudaEventRecord(g_forceUploadSourceSafeEvent,
                  solver->streams[CUDA_STREAM_MAIN]);
  g_forceUploadSourceSafePending = 1;
}

static void solver_wait_force_upload_source_event(void) {
  if ((g_forceUploadSourceSafeEventInit == 1) &&
      (g_forceUploadSourceSafePending == 1)) {
    cudaEventSynchronize(g_forceUploadSourceSafeEvent);
    g_forceUploadSourceSafePending = 0;
  }
}

static int32_t contrib_pack_use_omp(void) {
#ifndef _OPENMP
  return 0;
#else
  if (hercules_repro_mode()) {
    return 0;
  }

  if (g_contribPackOmpMode < 0) {
    const char* env = getenv("HERCULES_CONTRIB_PACK_OMP");

    if (env != NULL) {
      if ((strcmp(env, "0") == 0) || (strcasecmp(env, "no") == 0) ||
          (strcasecmp(env, "false") == 0) || (strcasecmp(env, "off") == 0)) {
        g_contribPackOmpMode = 0;
      } else {
        g_contribPackOmpMode = 1;
      }
    } else {
      g_contribPackOmpMode = 1;
    }

    if ((g_contribPackOmpMode == 1) && (Global.myID == 0) &&
        (g_contribPackOmpReported == 0)) {
      monitor_print("contribution pack OpenMP: enabled (threads=%d)\n",
                    omp_get_max_threads());
      g_contribPackOmpReported = 1;
    }
  }

  return g_contribPackOmpMode;
#endif
}

static int32_t contrib_unpack_use_omp(void) {
#ifndef _OPENMP
  return 0;
#else
  if (hercules_repro_mode()) {
    return 0;
  }

  if (g_contribUnpackOmpMode < 0) {
    const char* env = getenv("HERCULES_CONTRIB_UNPACK_OMP");

    if (env != NULL) {
      if ((strcmp(env, "0") == 0) || (strcasecmp(env, "no") == 0) ||
          (strcasecmp(env, "false") == 0) || (strcasecmp(env, "off") == 0)) {
        g_contribUnpackOmpMode = 0;
      } else {
        g_contribUnpackOmpMode = 1;
      }
    } else {
      g_contribUnpackOmpMode = 1;
    }

    if ((g_contribUnpackOmpMode == 1) && (Global.myID == 0) &&
        (g_contribUnpackOmpReported == 0)) {
      monitor_print("contribution unpack OpenMP: enabled (threads=%d)\n",
                    omp_get_max_threads());
      g_contribUnpackOmpReported = 1;
    }
  }

  return g_contribUnpackOmpMode;
#endif
}

static int32_t contrib_recv_use_waitsome(void) {
  if (hercules_repro_mode()) {
    return 0;
  }

  if (g_contribRecvMode < 0) {
    const char* env = getenv("HERCULES_CONTRIB_RECV_MODE");

    if (env != NULL) {
      if ((strcmp(env, "0") == 0) || (strcasecmp(env, "waitall") == 0) ||
          (strcasecmp(env, "deterministic") == 0) ||
          (strcasecmp(env, "ordered") == 0) ||
          (strcasecmp(env, "off") == 0) || (strcasecmp(env, "no") == 0) ||
          (strcasecmp(env, "false") == 0)) {
        g_contribRecvMode = 0;
      } else {
        g_contribRecvMode = 1;
      }
    } else {
      g_contribRecvMode = 1;
    }

    if ((Global.myID == 0) && (g_contribRecvReported == 0)) {
      monitor_print("contribution receive mode: %s\n",
                    g_contribRecvMode ? "waitsome" : "waitall");
      g_contribRecvReported = 1;
    }
  }

  return g_contribRecvMode;
}

static int32_t sparse_force_upload_mode(void) {
  if (g_sparseUploadMode < 0) {
    const char* env = getenv("HERCULES_SPARSE_UPLOAD_MODE");

    /* 0: auto heuristic, 1: force dense upload,
     * 2: force sparse upload, 3: force range upload */
    g_sparseUploadMode = 1;

    if (env != NULL) {
      if ((strcasecmp(env, "dense") == 0) || (strcmp(env, "1") == 0)) {
        g_sparseUploadMode = 1;
      } else if ((strcasecmp(env, "auto") == 0) || (strcmp(env, "0") == 0)) {
        g_sparseUploadMode = 0;
      } else if ((strcasecmp(env, "sparse") == 0) ||
                 (strcasecmp(env, "scatter") == 0) ||
                 (strcmp(env, "2") == 0)) {
        g_sparseUploadMode = 2;
      } else if ((strcasecmp(env, "range") == 0) ||
                 (strcasecmp(env, "ranges") == 0) ||
                 (strcasecmp(env, "batched") == 0) ||
                 (strcmp(env, "3") == 0)) {
        g_sparseUploadMode = 3;
      }
    }

    if ((Global.myID == 0) && (g_sparseUploadModeReported == 0)) {
      const char* modeName =
          (g_sparseUploadMode == 1) ? "dense" :
          (g_sparseUploadMode == 2) ? "sparse" :
          (g_sparseUploadMode == 3) ? "range" : "auto";
      monitor_print("Sparse force upload mode request: %s\n", modeName);
      g_sparseUploadModeReported = 1;
    }
  }

  return g_sparseUploadMode;
}

static void print_timing_stat() {
  double TotalMeshingTime;
  double OtherPhasesTime = 0.0;
  double TotalMeshingTimeMax;
  double OtherPhasesTimeMax = 0.0;
  double BalanceTime = 0.0;
  double BalanceTimeMax = 0.0;
  double CorrectPropsTime = 0.0;
  double CorrectPropsTimeMax = 0.0;
  double SourceInitMax;
  double SourceInitForAccounting;
  double SolverMax;
  double TotalWallClockMax;
  double AccountedMax;
  double UnaccountedMax;
  if (Timer_Exists("Octor Balancetree")) {
    BalanceTime = Timer_Value("Octor Balancetree", (TimerKind)0);
    BalanceTimeMax = Timer_Value("Octor Balancetree", MAX);
  }
  if (Timer_Exists("Mesh correct properties")) {
    CorrectPropsTime = Timer_Value("Mesh correct properties", (TimerKind)0);
    CorrectPropsTimeMax = Timer_Value("Mesh correct properties", MAX);
  }

  TotalMeshingTime = Timer_Value("Octor Newtree", (TimerKind)0) +
                     Timer_Value("Octor Refinetree", (TimerKind)0) +
                     BalanceTime +
                     Timer_Value("Octor Partitiontree", (TimerKind)0) +
                     Timer_Value("Octor Extractmesh", (TimerKind)0) +
                     CorrectPropsTime +
                     Timer_Value("Mesh Stats Print", (TimerKind)0);

  if (Timer_Exists("Slice CVM")) {
    TotalMeshingTime += Timer_Value("Slice CVM", (TimerKind)0);
  }

  if (Timer_Exists("Mesh Out")) {
    TotalMeshingTime += Timer_Value("Mesh Out", (TimerKind)0);
  }

  if (Param.includeBuildings == YES) {
    TotalMeshingTime += Timer_Value("Carve Buildings", (TimerKind)0);
  }

  TotalMeshingTimeMax = Timer_Value("Octor Newtree", MAX) +
                        Timer_Value("Octor Refinetree", MAX) +
                        BalanceTimeMax +
                        Timer_Value("Octor Partitiontree", MAX) +
                        Timer_Value("Octor Extractmesh", MAX) +
                        CorrectPropsTimeMax +
                        Timer_Value("Mesh Stats Print", MAX);

  if (Timer_Exists("Slice CVM")) {
    TotalMeshingTimeMax += Timer_Value("Slice CVM", MAX);
  }

  if (Timer_Exists("Mesh Out")) {
    TotalMeshingTimeMax += Timer_Value("Mesh Out", MAX);
  }

  if (Param.includeBuildings == YES && Timer_Exists("Carve Buildings")) {
    TotalMeshingTimeMax += Timer_Value("Carve Buildings", MAX);
  }

  if (Timer_Exists("Global Init")) {
    OtherPhasesTime += Timer_Value("Global Init", (TimerKind)0);
    OtherPhasesTimeMax += Timer_Value("Global Init", MAX);
  }
  if (Timer_Exists("Pre Solver Setup")) {
    OtherPhasesTime += Timer_Value("Pre Solver Setup", (TimerKind)0);
    OtherPhasesTimeMax += Timer_Value("Pre Solver Setup", MAX);
  }
  if (Timer_Exists("Post Solver Finalize")) {
    OtherPhasesTime += Timer_Value("Post Solver Finalize", (TimerKind)0);
    OtherPhasesTimeMax += Timer_Value("Post Solver Finalize", MAX);
  }

  if (Timer_Exists("Phase Sync/Coordination")) {
    OtherPhasesTime += Timer_Value("Phase Sync/Coordination", (TimerKind)0);
    OtherPhasesTimeMax += Timer_Value("Phase Sync/Coordination", MAX);
  }

  if (Timer_Exists("Init Drm Parameters")) {
    OtherPhasesTime += Timer_Value("Init Drm Parameters", (TimerKind)0);
    OtherPhasesTimeMax += Timer_Value("Init Drm Parameters", MAX);
  }

  if (Timer_Exists("Drm Init")) {
    OtherPhasesTime += Timer_Value("Drm Init", (TimerKind)0);
    OtherPhasesTimeMax += Timer_Value("Drm Init", MAX);
  }

  if (Timer_Exists("Topography Solver Init")) {
    OtherPhasesTime += Timer_Value("Topography Solver Init", (TimerKind)0);
    OtherPhasesTimeMax += Timer_Value("Topography Solver Init", MAX);
  }

  SourceInitMax = Timer_Value("Source Init", MAX);
  SourceInitForAccounting = SourceInitMax;
  if (Timer_Exists("Pre Solver Setup")) {
    SourceInitForAccounting = 0.0;
  }
  SolverMax = Timer_Value("Solver", MAX);
  TotalWallClockMax = Timer_Value("Total Wall Clock", MAX);
  AccountedMax =
      TotalMeshingTimeMax + SourceInitForAccounting + SolverMax +
      OtherPhasesTimeMax;
  UnaccountedMax = TotalWallClockMax - AccountedMax;

  printf(
      "\n\n__________________________Raw Timers__________________________\n\n");
  Timer_PrintAll(comm_solver);
  printf("\n\n\n\n\n");

  printf(
      "\n\n__________________________Timer "
      "Statistics__________________________\n\n");

  printf("\n_____________Summary_____________\n");
  printf("Max Frequency             : %.2f\n", Param.theFreq);
  printf("Vs                        : %.2f\n", Param.theVsCut);
  printf("Total elements            : %" INT64_FMT "\n", Global.theETotal);
  printf("Elements/PE               : %" INT64_FMT "\n",
         Global.theETotal / Global.theGroupSize);
  printf("Simulation duration       : %.2f seconds\n",
         Param.theEndT - Param.theStartT);
  printf("Total steps               : %d\n", Param.theTotalSteps);
  printf("DeltaT used               : %.6f seconds\n", Param.theDeltaT);
  printf("Critical deltaT           : %.6f seconds\n", Global.theCriticalT);
  printf("\n");
  printf("Accounted (MAX basis)     : %.2f seconds\n", AccountedMax);
    printf("Unaccounted vs wall clock (MAX basis): %.2f seconds\n",
      UnaccountedMax);
  printf("\n");
  printf("Total Wall Clock          : %.2f seconds\n",
         Timer_Value("Total Wall Clock", (TimerKind)0));
    printf("Total Wall Clock (MAX)    : %.2f seconds\n", TotalWallClockMax);
  printf("Time/step                 : %.6f seconds\n",
         Timer_Value("Solver", (TimerKind)0) / Param.theTotalSteps);
  printf("Time/step/(elem/PE)       : %.6f millisec\n",
         Timer_Value("Solver", (TimerKind)0) * 1000.0 / Param.theTotalSteps /
             (Global.theETotal * 1.0 / Global.theGroupSize));
  printf(
      "Simulation Rate Variation : %.3f (Average)   %.3f (Min)   %.3f (Max)  "
      "(sec/%d timesteps)\n",
      (Timer_Value("Solver", (TimerKind)0) / Param.theTotalSteps) *
          Param.monitor_stats_rate,
      Global.fastestTimeSteps, Global.slowestTimeSteps,
      Param.monitor_stats_rate);
  printf("\n");

  printf("\n____________Breakdown____________\n");
  printf("TOTAL MESHING                       : %.2f seconds\n",
         TotalMeshingTime);
  printf("    Octor Newtree                   : %.2f seconds\n",
         Timer_Value("Octor Newtree", (TimerKind)0));
  printf("    Octor Refinetree                : %.2f seconds\n",
         Timer_Value("Octor Refinetree", (TimerKind)0));
    if (Timer_Exists("Octor Balancetree")) {
      printf("    Octor Balancetree               : %.2f seconds\n",
        Timer_Value("Octor Balancetree", (TimerKind)0));
    }
  if (Timer_Exists("Carve Buildings"))
    printf("    Octor Carve Buildings           : %.2f seconds\n",
           Timer_Value("Carve Buildings", (TimerKind)0));
  printf("    Octor Partitiontree             : %.2f seconds\n",
         Timer_Value("Octor Partitiontree", (TimerKind)0));
  printf("    Octor Extractmesh               : %.2f seconds\n",
         Timer_Value("Octor Extractmesh", (TimerKind)0));
    if (Timer_Exists("Slice CVM"))
      printf("    Slice CVM                       : %.2f seconds\n",
        Timer_Value("Slice CVM", (TimerKind)0));
  if (Timer_Exists("Mesh correct properties")) {
    printf("    Mesh correct properties         : %.2f seconds\n",
           Timer_Value("Mesh correct properties", (TimerKind)0));
  }
    if (Timer_Exists("Mesh Out"))
      printf("    Mesh Out                        : %.2f seconds\n",
        Timer_Value("Mesh Out", (TimerKind)0));
  printf("    Mesh Stats Print                : %.2f seconds\n",
         Timer_Value("Mesh Stats Print", (TimerKind)0));
  printf("\n");

  if (Param.drmImplement == YES) {
    printf(
        "DRM INIT PARAMETERS                 : %.2f (Max) %.2f (Min) seconds\n",
        Timer_Value("Init Drm Parameters", MAX),
        Timer_Value("Init Drm Parameters", MIN));
    printf("\n");

    printf(
        "DRM INITIALIZATION                  : %.2f (Max) %.2f (Min) seconds\n",
        Timer_Value("Drm Init", MAX), Timer_Value("Drm Init", MIN));

    if (Param.theDrmPart == PART2) {
      printf(
          "    Find Drm File To Readjust       : %.2f (Average)   %.2f (Max) "
          "%.2f (Min) seconds\n",
          Timer_Value("Find Drm File To Readjust", AVERAGE),
          Timer_Value("Find Drm File To Readjust", MAX),
          Timer_Value("Find Drm File To Readjust", MIN));

      printf(
          "    Fill Drm Struct                 : %.2f (Average)   %.2f (Max) "
          "%.2f (Min) seconds\n",
          Timer_Value("Fill Drm Struct", AVERAGE),
          Timer_Value("Fill Drm Struct", MAX),
          Timer_Value("Fill Drm Struct", MIN));

      printf(
          "    Comm of Drm Coordinates         : %.2f (Average)   %.2f (Max) "
          "%.2f (Min) seconds\n",
          Timer_Value("Comm of Drm Coordinates", AVERAGE),
          Timer_Value("Comm of Drm Coordinates", MAX),
          Timer_Value("Comm of Drm Coordinates", MIN));

      printf(
          "    Find Which Drm Files To Print   : %.2f (Average)   %.2f (Max) "
          "%.2f (Min) seconds\n",
          Timer_Value("Find Which Drm Files To Print", AVERAGE),
          Timer_Value("Find Which Drm Files To Print", MAX),
          Timer_Value("Find Which Drm Files To Print", MIN));

      printf(
          "    Read And Rearrange Drm Files    : %.2f (Average)   %.2f (Max) "
          "%.2f (Min) seconds\n",
          Timer_Value("Read And Rearrange Drm Files", AVERAGE),
          Timer_Value("Read And Rearrange Drm Files", MAX),
          Timer_Value("Read And Rearrange Drm Files", MIN));

      printf(
          "    Find Which Drm Files To Read    : %.2f (Average)   %.2f (Max) "
          "%.2f (Min) seconds\n",
          Timer_Value("Find Which Drm Files To Read", AVERAGE),
          Timer_Value("Find Which Drm Files To Read", MAX),
          Timer_Value("Find Which Drm Files To Read", MIN));

      printf(
          "    Locate where I am in file       : %.2f (Average)   %.2f (Max) "
          "%.2f (Min) seconds\n",
          Timer_Value("Locate where I am in file", AVERAGE),
          Timer_Value("Locate where I am in file", MAX),
          Timer_Value("Locate where I am in file", MIN));
    }
    printf("\n");
  }

  printf(
      "SOURCE INITIALIZATION               : %.2f (Max) %.2f (Min) seconds\n",
      Timer_Value("Source Init", MAX), Timer_Value("Source Init", MIN));
  if (Timer_Exists("Global Init") || Timer_Exists("Pre Solver Setup") ||
      Timer_Exists("Post Solver Finalize") ||
      Timer_Exists("Topography Solver Init")) {
    printf("OTHER PHASES (outside mesh/solver/source): %.2f seconds\n",
           OtherPhasesTime);
    if (Timer_Exists("Global Init")) {
      printf("    Global Init                    : %.2f (Average)   %.2f (Max) %.2f "
             "(Min) seconds\n",
             Timer_Value("Global Init", AVERAGE),
             Timer_Value("Global Init", MAX), Timer_Value("Global Init", MIN));
    }
    if (Timer_Exists("Pre Solver Setup")) {
      printf("    Pre Solver Setup               : %.2f (Average)   %.2f (Max) %.2f "
             "(Min) seconds\n",
             Timer_Value("Pre Solver Setup", AVERAGE),
             Timer_Value("Pre Solver Setup", MAX),
             Timer_Value("Pre Solver Setup", MIN));
      if (Timer_Exists("Pre Setup Stations Init")) {
        printf("      Pre Setup Stations Init     : %.2f (Average)   %.2f (Max) %.2f "
               "(Min) seconds\n",
               Timer_Value("Pre Setup Stations Init", AVERAGE),
               Timer_Value("Pre Setup Stations Init", MAX),
               Timer_Value("Pre Setup Stations Init", MIN));
        if (Timer_Exists("Stations Read Input")) {
          printf("        Stations Read Input       : %.2f (Average)   %.2f (Max) %.2f "
                 "(Min) seconds\n",
                 Timer_Value("Stations Read Input", AVERAGE),
                 Timer_Value("Stations Read Input", MAX),
                 Timer_Value("Stations Read Input", MIN));
        }
        if (Timer_Exists("Stations Broadcast")) {
          printf("        Stations Broadcast        : %.2f (Average)   %.2f (Max) %.2f "
                 "(Min) seconds\n",
                 Timer_Value("Stations Broadcast", AVERAGE),
                 Timer_Value("Stations Broadcast", MAX),
                 Timer_Value("Stations Broadcast", MIN));
        }
        if (Timer_Exists("Stations Locate")) {
          printf("        Stations Locate           : %.2f (Average)   %.2f (Max) %.2f "
                 "(Min) seconds\n",
                 Timer_Value("Stations Locate", AVERAGE),
                 Timer_Value("Stations Locate", MAX),
                 Timer_Value("Stations Locate", MIN));
        }
        if (Timer_Exists("Stations Local Init")) {
          printf("        Stations Local Init       : %.2f (Average)   %.2f (Max) %.2f "
                 "(Min) seconds\n",
                 Timer_Value("Stations Local Init", AVERAGE),
                 Timer_Value("Stations Local Init", MAX),
                 Timer_Value("Stations Local Init", MIN));
        }
      }
      if (Timer_Exists("Pre Setup Topo Stations Init")) {
        printf("      Pre Setup Topo Stations Init: %.2f (Average)   %.2f (Max) %.2f "
               "(Min) seconds\n",
               Timer_Value("Pre Setup Topo Stations Init", AVERAGE),
               Timer_Value("Pre Setup Topo Stations Init", MAX),
               Timer_Value("Pre Setup Topo Stations Init", MIN));
      }
      if (Timer_Exists("Pre Setup Planes Init")) {
        printf("      Pre Setup Planes Init       : %.2f (Average)   %.2f (Max) %.2f "
               "(Min) seconds\n",
               Timer_Value("Pre Setup Planes Init", AVERAGE),
               Timer_Value("Pre Setup Planes Init", MAX),
               Timer_Value("Pre Setup Planes Init", MIN));
      }
      if (Timer_Exists("Pre Setup Solver Init")) {
        printf("      Pre Setup Solver Init       : %.2f (Average)   %.2f (Max) %.2f "
               "(Min) seconds\n",
               Timer_Value("Pre Setup Solver Init", AVERAGE),
               Timer_Value("Pre Setup Solver Init", MAX),
               Timer_Value("Pre Setup Solver Init", MIN));
      }
      if (Timer_Exists("Pre Setup Nonlinear Init")) {
        printf("      Pre Setup Nonlinear Init    : %.2f (Average)   %.2f (Max) %.2f "
               "(Min) seconds\n",
               Timer_Value("Pre Setup Nonlinear Init", AVERAGE),
               Timer_Value("Pre Setup Nonlinear Init", MAX),
               Timer_Value("Pre Setup Nonlinear Init", MIN));
      }
      if (Timer_Exists("Pre Setup Topo Stats")) {
        printf("      Pre Setup Topo Stats        : %.2f (Average)   %.2f (Max) %.2f "
               "(Min) seconds\n",
               Timer_Value("Pre Setup Topo Stats", AVERAGE),
               Timer_Value("Pre Setup Topo Stats", MAX),
               Timer_Value("Pre Setup Topo Stats", MIN));
      }
      if (Timer_Exists("Pre Setup PlaneWaves Init")) {
        printf("      Pre Setup PlaneWaves Init   : %.2f (Average)   %.2f (Max) %.2f "
               "(Min) seconds\n",
               Timer_Value("Pre Setup PlaneWaves Init", AVERAGE),
               Timer_Value("Pre Setup PlaneWaves Init", MAX),
               Timer_Value("Pre Setup PlaneWaves Init", MIN));
      }
      if (Timer_Exists("Pre Setup Damping/Stiffness Init")) {
        printf("      Pre Setup Damping/Stiffness : %.2f (Average)   %.2f (Max) %.2f "
               "(Min) seconds\n",
               Timer_Value("Pre Setup Damping/Stiffness Init", AVERAGE),
               Timer_Value("Pre Setup Damping/Stiffness Init", MAX),
               Timer_Value("Pre Setup Damping/Stiffness Init", MIN));
      }
      if (Timer_Exists("Pre Setup Output Init")) {
        printf("      Pre Setup Output Init       : %.2f (Average)   %.2f (Max) %.2f "
               "(Min) seconds\n",
               Timer_Value("Pre Setup Output Init", AVERAGE),
               Timer_Value("Pre Setup Output Init", MAX),
               Timer_Value("Pre Setup Output Init", MIN));
      }
      if (Timer_Exists("Pre Setup Source Init")) {
        printf("      Pre Setup Source Init       : %.2f (Average)   %.2f (Max) %.2f "
               "(Min) seconds\n",
               Timer_Value("Pre Setup Source Init", AVERAGE),
               Timer_Value("Pre Setup Source Init", MAX),
               Timer_Value("Pre Setup Source Init", MIN));
      }
    }
    if (Timer_Exists("Post Solver Finalize")) {
      printf("    Post Solver Finalize           : %.2f (Average)   %.2f (Max) %.2f "
             "(Min) seconds\n",
             Timer_Value("Post Solver Finalize", AVERAGE),
             Timer_Value("Post Solver Finalize", MAX),
             Timer_Value("Post Solver Finalize", MIN));
    }
    if (Timer_Exists("Topography Solver Init")) {
      printf("    Topography Solver Init         : %.2f (Average)   %.2f (Max) %.2f "
             "(Min) seconds\n",
             Timer_Value("Topography Solver Init", AVERAGE),
             Timer_Value("Topography Solver Init", MAX),
             Timer_Value("Topography Solver Init", MIN));
      if (Timer_Exists("Topo Init Select Elements")) {
        printf("      Topo Init Select Elements   : %.2f (Average)   %.2f (Max) %.2f "
               "(Min) seconds\n",
               Timer_Value("Topo Init Select Elements", AVERAGE),
               Timer_Value("Topo Init Select Elements", MAX),
               Timer_Value("Topo Init Select Elements", MIN));
      }
      if (Timer_Exists("Topo Init Build Node Sets")) {
        printf("      Topo Init Build Node Sets   : %.2f (Average)   %.2f (Max) %.2f "
               "(Min) seconds\n",
               Timer_Value("Topo Init Build Node Sets", AVERAGE),
               Timer_Value("Topo Init Build Node Sets", MAX),
               Timer_Value("Topo Init Build Node Sets", MIN));
      }
      if (Timer_Exists("Topo Init Build Constants")) {
        printf("      Topo Init Build Constants   : %.2f (Average)   %.2f (Max) %.2f "
               "(Min) seconds\n",
               Timer_Value("Topo Init Build Constants", AVERAGE),
               Timer_Value("Topo Init Build Constants", MAX),
               Timer_Value("Topo Init Build Constants", MIN));
      }
    }
    if (Timer_Exists("Phase Sync/Coordination")) {
      printf("    Phase Sync/Coordination        : %.2f (Average)   %.2f (Max) %.2f "
             "(Min) seconds\n",
             Timer_Value("Phase Sync/Coordination", AVERAGE),
             Timer_Value("Phase Sync/Coordination", MAX),
             Timer_Value("Phase Sync/Coordination", MIN));
    }
    if (Timer_Exists("Init Drm Parameters")) {
      printf("    Init Drm Parameters            : %.2f (Average)   %.2f (Max) %.2f "
             "(Min) seconds\n",
             Timer_Value("Init Drm Parameters", AVERAGE),
             Timer_Value("Init Drm Parameters", MAX),
             Timer_Value("Init Drm Parameters", MIN));
    }
    if (Timer_Exists("Drm Init")) {
      printf("    Drm Init                       : %.2f (Average)   %.2f (Max) %.2f "
             "(Min) seconds\n",
             Timer_Value("Drm Init", AVERAGE),
             Timer_Value("Drm Init", MAX), Timer_Value("Drm Init", MIN));
    }
  }
  printf("\n");
  printf("TOTAL SOLVER                        : %.2f seconds\n",
         Timer_Value("Solver", (TimerKind)0));
  printf(
      "    Read My Forces                  : %.2f (Average)   %.2f (Max) %.2f "
      "(Min) seconds\n",
      Timer_Value("Read My Forces", AVERAGE),
      Timer_Value("Read My Forces", MAX), Timer_Value("Read My Forces", MIN));
  printf(
      "    Compute addforces s             : %.2f (Average)   %.2f (Max) %.2f "
      "(Min) seconds\n",
      Timer_Value("Compute addforces s", AVERAGE),
      Timer_Value("Compute addforces s", MAX),
      Timer_Value("Compute addforces s", MIN));
  printf(
      "    Compute addforces e             : %.2f (Average)   %.2f (Max) %.2f "
      "(Min) seconds\n",
      Timer_Value("Compute addforces e", AVERAGE),
      Timer_Value("Compute addforces e", MAX),
      Timer_Value("Compute addforces e", MIN));
  printf(
      "    Compute Damping addforce        : %.2f (Average)   %.2f (Max) %.2f "
      "(Min) seconds\n",
      Timer_Value("Damping addforce", AVERAGE),
      Timer_Value("Damping addforce", MAX),
      Timer_Value("Damping addforce", MIN));
  if (Timer_Exists("Compute Non-linear Entities"))
    printf(
        "    Compute Non-linear Entities     : %.2f (Average)   %.2f (Max) "
        "%.2f (Min) seconds\n",
        Timer_Value("Compute Non-linear Entities", AVERAGE),
        Timer_Value("Compute Non-linear Entities", MAX),
        Timer_Value("Compute Non-linear Entities", MIN));
  if (Timer_Exists("Compute addforces Non-linear")) {
    printf(
        "    Compute addforces Non-linear    : %.2f (Average)   %.2f (Max) "
        "%.2f (Min) seconds\n",
        Timer_Value("Compute addforces Non-linear", AVERAGE),
        Timer_Value("Compute addforces Non-linear", MAX),
        Timer_Value("Compute addforces Non-linear", MIN));
  }
  if (Timer_Exists("Compute addforces gravity")) {
    printf(
        "    Compute addforces gravity       : %.2f (Average)   %.2f (Max) "
        "%.2f (Min) seconds\n",
        Timer_Value("Compute addforces gravity", AVERAGE),
        Timer_Value("Compute addforces gravity", MAX),
        Timer_Value("Compute addforces gravity", MIN));
  }

  if (Timer_Exists("Solver drm force compute")) {
    printf(
        "    Solver drm force compute        : %.2f (Average)   %.2f (Max) "
        "%.2f (Min) seconds\n",
        Timer_Value("Solver drm force compute", AVERAGE),
        Timer_Value("Solver drm force compute", MAX),
        Timer_Value("Solver drm force compute", MIN));
  }
  if (Timer_Exists("Compute addforces Topography")) {
    printf(
        "    Compute addforces Topography    : %.2f (Average)   %.2f (Max) "
        "%.2f (Min) seconds\n",
        Timer_Value("Compute addforces Topography", AVERAGE),
        Timer_Value("Compute addforces Topography", MAX),
        Timer_Value("Compute addforces Topography", MIN));
  }
  if (Timer_Exists("Compute addforces Incident Plane Waves")) {
    printf(
        "    Compute addforces Plane Waves   : %.2f (Average)   %.2f (Max) "
        "%.2f (Min) seconds\n",
        Timer_Value("Compute addforces Incident Plane Waves", AVERAGE),
        Timer_Value("Compute addforces Incident Plane Waves", MAX),
        Timer_Value("Compute addforces Incident Plane Waves", MIN));
  }
  if (Timer_Exists("Load Fixedbase Disps")) {
    printf(
        "    Load Fixedbase Disps            : %.2f (Average)   %.2f (Max) "
        "%.2f (Min) seconds\n",
        Timer_Value("Load Fixedbase Disps", AVERAGE),
        Timer_Value("Load Fixedbase Disps", MAX),
        Timer_Value("Load Fixedbase Disps", MIN));
  }
  if (Timer_Exists("Physics Force Prep") ||
      Timer_Exists("Physics Force Kernels Enqueue") ||
      Timer_Exists("Physics Force Wait/Post") ||
      Timer_Exists("Station Accel History Copy") ||
      Timer_Exists("Force Stream Sync") ||
      Timer_Exists("Force stiffness wait only") ||
      Timer_Exists("Force pre-stiffness wait only") ||
      Timer_Exists("Force stiffness core wait only") ||
      Timer_Exists("Force damping conv wait only") ||
      Timer_Exists("Force damping local wait only") ||
      Timer_Exists("Force kernel wait only") ||
      Timer_Exists("Force postcopy sync wait") ||
      Timer_Exists("Force Post CPU") ||
      Timer_Exists("Physics Displacement Enqueue") ||
      Timer_Exists("Physics Displacement Wait/Post")) {
    printf("    Physics phase internals\n");
    if (Timer_Exists("Physics Force Prep")) {
      printf(
          "      Force Prep                    : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("Physics Force Prep", AVERAGE),
          Timer_Value("Physics Force Prep", MAX),
          Timer_Value("Physics Force Prep", MIN));
    }
    if (Timer_Exists("Physics Force Kernels Enqueue")) {
      printf(
          "      Force Kernels Enqueue         : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("Physics Force Kernels Enqueue", AVERAGE),
          Timer_Value("Physics Force Kernels Enqueue", MAX),
          Timer_Value("Physics Force Kernels Enqueue", MIN));
    }
    if (Timer_Exists("Physics Force Wait/Post")) {
      printf(
          "      Force Wait/Post               : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("Physics Force Wait/Post", AVERAGE),
          Timer_Value("Physics Force Wait/Post", MAX),
          Timer_Value("Physics Force Wait/Post", MIN));
    }
    if (Timer_Exists("Station Accel History Copy")) {
      printf(
          "      Station Accel History Copy    : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("Station Accel History Copy", AVERAGE),
          Timer_Value("Station Accel History Copy", MAX),
          Timer_Value("Station Accel History Copy", MIN));
    }
    if (Timer_Exists("Force Stream Sync")) {
      printf(
          "      Force Stream Sync             : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("Force Stream Sync", AVERAGE),
          Timer_Value("Force Stream Sync", MAX),
          Timer_Value("Force Stream Sync", MIN));
    }
    if (Timer_Exists("Force stiffness wait only")) {
      printf(
          "      Force stiffness wait only     : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("Force stiffness wait only", AVERAGE),
          Timer_Value("Force stiffness wait only", MAX),
          Timer_Value("Force stiffness wait only", MIN));
    }
    if (Timer_Exists("Force pre-stiffness wait only")) {
      printf(
          "      Force pre-stiffness wait only : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("Force pre-stiffness wait only", AVERAGE),
          Timer_Value("Force pre-stiffness wait only", MAX),
          Timer_Value("Force pre-stiffness wait only", MIN));
    }
    if (Timer_Exists("Force stiffness core wait only")) {
      printf(
          "      Force stiffness core wait only: %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("Force stiffness core wait only", AVERAGE),
          Timer_Value("Force stiffness core wait only", MAX),
          Timer_Value("Force stiffness core wait only", MIN));
    }
    if (Timer_Exists("Force damping conv wait only")) {
      printf(
          "      Force damping conv wait only  : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("Force damping conv wait only", AVERAGE),
          Timer_Value("Force damping conv wait only", MAX),
          Timer_Value("Force damping conv wait only", MIN));
    }
    if (Timer_Exists("Force damping local wait only")) {
      printf(
          "      Force damping local wait only : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("Force damping local wait only", AVERAGE),
          Timer_Value("Force damping local wait only", MAX),
          Timer_Value("Force damping local wait only", MIN));
    }
    if (Timer_Exists("Force kernel wait only")) {
      printf(
          "      Force kernel wait only        : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("Force kernel wait only", AVERAGE),
          Timer_Value("Force kernel wait only", MAX),
          Timer_Value("Force kernel wait only", MIN));
    }
    if (Timer_Exists("Force postcopy sync wait")) {
      printf(
          "      Force postcopy sync wait      : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("Force postcopy sync wait", AVERAGE),
          Timer_Value("Force postcopy sync wait", MAX),
          Timer_Value("Force postcopy sync wait", MIN));
    }
    if (Timer_Exists("Force Post CPU")) {
      printf(
          "      Force Post CPU                : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("Force Post CPU", AVERAGE),
          Timer_Value("Force Post CPU", MAX),
          Timer_Value("Force Post CPU", MIN));
    }
    if (Timer_Exists("Physics Displacement Enqueue")) {
      printf(
          "      Displacement Enqueue          : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("Physics Displacement Enqueue", AVERAGE),
          Timer_Value("Physics Displacement Enqueue", MAX),
          Timer_Value("Physics Displacement Enqueue", MIN));
    }
    if (Timer_Exists("Physics Displacement Wait/Post")) {
      printf(
          "      Displacement Wait/Post        : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("Physics Displacement Wait/Post", AVERAGE),
          Timer_Value("Physics Displacement Wait/Post", MAX),
          Timer_Value("Physics Displacement Wait/Post", MIN));
    }
  }
  printf(
      "    1st schedule send data          : %.2f (Average)   %.2f (Max) %.2f "
      "(Min) seconds\n",
      Timer_Value("1st schedule send data (contribution)", AVERAGE),
      Timer_Value("1st schedule send data (contribution)", MAX),
      Timer_Value("1st schedule send data (contribution)", MIN));
  printf(
      "    1st compute adjust              : %.2f (Average)   %.2f (Max) %.2f "
      "(Min) seconds\n",
      Timer_Value("1st compute adjust (distribution)", AVERAGE),
      Timer_Value("1st compute adjust (distribution)", MAX),
      Timer_Value("1st compute adjust (distribution)", MIN));
  printf(
      "    2nd schedule send data          : %.2f (Average)   %.2f (Max) %.2f "
      "(Min) seconds\n",
      Timer_Value("2nd schedule send data (contribution)", AVERAGE),
      Timer_Value("2nd schedule send data (contribution)", MAX),
      Timer_Value("2nd schedule send data (contribution)", MIN));
  printf(
      "    Compute new displacement        : %.2f (Average)   %.2f (Max) %.2f "
      "(Min) seconds\n",
      Timer_Value("Compute new displacement", AVERAGE),
      Timer_Value("Compute new displacement", MAX),
      Timer_Value("Compute new displacement", MIN));
  if (Timer_Exists("Disp kernel enqueue")) {
    printf(
        "      Disp kernel enqueue          : %.2f (Average)   %.2f (Max) %.2f "
        "(Min) seconds\n",
        Timer_Value("Disp kernel enqueue", AVERAGE),
        Timer_Value("Disp kernel enqueue", MAX),
        Timer_Value("Disp kernel enqueue", MIN));
  }
  if (Timer_Exists("Disp host force reset")) {
    printf(
        "      Disp host force reset        : %.2f (Average)   %.2f (Max) %.2f "
        "(Min) seconds\n",
        Timer_Value("Disp host force reset", AVERAGE),
        Timer_Value("Disp host force reset", MAX),
        Timer_Value("Disp host force reset", MIN));
  }
  if (Timer_Exists("Disp force upload safe wait")) {
    printf(
        "      Disp force upload safe wait  : %.2f (Average)   %.2f (Max) %.2f "
        "(Min) seconds\n",
        Timer_Value("Disp force upload safe wait", AVERAGE),
        Timer_Value("Disp force upload safe wait", MAX),
        Timer_Value("Disp force upload safe wait", MIN));
  }
  if (Timer_Exists("Disp host force clear only")) {
    printf(
        "      Disp host force clear only   : %.2f (Average)   %.2f (Max) %.2f "
        "(Min) seconds\n",
        Timer_Value("Disp host force clear only", AVERAGE),
        Timer_Value("Disp host force clear only", MAX),
        Timer_Value("Disp host force clear only", MIN));
  }
  if (Timer_Exists("Disp kernel sync wait")) {
    printf(
        "      Disp kernel sync wait        : %.2f (Average)   %.2f (Max) %.2f "
        "(Min) seconds\n",
        Timer_Value("Disp kernel sync wait", AVERAGE),
        Timer_Value("Disp kernel sync wait", MAX),
        Timer_Value("Disp kernel sync wait", MIN));
  }
  if (Timer_Exists("Disp kernel wait only")) {
    printf(
        "      Disp kernel wait only        : %.2f (Average)   %.2f (Max) %.2f "
        "(Min) seconds\n",
        Timer_Value("Disp kernel wait only", AVERAGE),
        Timer_Value("Disp kernel wait only", MAX),
        Timer_Value("Disp kernel wait only", MIN));
  }
  if (Timer_Exists("Disp postcopy sync wait")) {
    printf(
        "      Disp postcopy sync wait      : %.2f (Average)   %.2f (Max) %.2f "
        "(Min) seconds\n",
        Timer_Value("Disp postcopy sync wait", AVERAGE),
        Timer_Value("Disp postcopy sync wait", MAX),
        Timer_Value("Disp postcopy sync wait", MIN));
  }
  if (Timer_Exists("Disp boundary host scatter")) {
    printf(
        "      Disp boundary host scatter   : %.2f (Average)   %.2f (Max) %.2f "
        "(Min) seconds\n",
        Timer_Value("Disp boundary host scatter", AVERAGE),
        Timer_Value("Disp boundary host scatter", MAX),
        Timer_Value("Disp boundary host scatter", MIN));
  }
  printf(
      "    3rd schedule send data          : %.2f (Average)   %.2f (Max) %.2f "
      "(Min) seconds\n",
      Timer_Value("3rd schedule send data (sharing)", AVERAGE),
      Timer_Value("3rd schedule send data (sharing)", MAX),
      Timer_Value("3rd schedule send data (sharing)", MIN));
  printf(
      "    2nd compute adjust              : %.2f (Average)   %.2f (Max) %.2f "
      "(Min) seconds\n",
      Timer_Value("2nd compute adjust (assignment)", AVERAGE),
      Timer_Value("2nd compute adjust (assignment)", MAX),
      Timer_Value("2nd compute adjust (assignment)", MIN));
  printf(
      "    4th schedule send data          : %.2f (Average)   %.2f (Max) %.2f "
      "(Min) seconds\n",
      Timer_Value("4th schedule send data (sharing)", AVERAGE),
      Timer_Value("4th schedule send data (sharing)", MAX),
      Timer_Value("4th schedule send data (sharing)", MIN));
  if (Timer_Exists("Comm Force Exchange") ||
      Timer_Exists("Comm Displacement Exchange")) {
    printf("    Communication phase internals\n");
    if (Timer_Exists("Comm Force Exchange")) {
      printf(
          "      Force Exchange                : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("Comm Force Exchange", AVERAGE),
          Timer_Value("Comm Force Exchange", MAX),
          Timer_Value("Comm Force Exchange", MIN));
    }
    if (Timer_Exists("Comm Displacement Exchange")) {
      printf(
          "      Displacement Exchange         : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("Comm Displacement Exchange", AVERAGE),
          Timer_Value("Comm Displacement Exchange", MAX),
          Timer_Value("Comm Displacement Exchange", MIN));
    }
  }
  if (Timer_Exists("Schedule Post Recv") ||
      Timer_Exists("Schedule Pack Outgoing") ||
      Timer_Exists("Schedule Post Send") ||
      Timer_Exists("Schedule Wait Recv") ||
      Timer_Exists("Schedule Wait Send") ||
      Timer_Exists("Schedule Unpack Incoming")) {
    printf("    Schedule internals (all send phases)\n");
    if (Timer_Exists("Schedule Post Recv")) {
      printf(
          "      Post Recv                     : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("Schedule Post Recv", AVERAGE),
          Timer_Value("Schedule Post Recv", MAX),
          Timer_Value("Schedule Post Recv", MIN));
    }
    if (Timer_Exists("Schedule Pack Outgoing")) {
      printf(
          "      Pack Outgoing                 : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("Schedule Pack Outgoing", AVERAGE),
          Timer_Value("Schedule Pack Outgoing", MAX),
          Timer_Value("Schedule Pack Outgoing", MIN));
    }
    if (Timer_Exists("Schedule Post Send")) {
      printf(
          "      Post Send                     : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("Schedule Post Send", AVERAGE),
          Timer_Value("Schedule Post Send", MAX),
          Timer_Value("Schedule Post Send", MIN));
    }
    if (Timer_Exists("Schedule Wait Recv")) {
      printf(
          "      Wait Recv                     : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("Schedule Wait Recv", AVERAGE),
          Timer_Value("Schedule Wait Recv", MAX),
          Timer_Value("Schedule Wait Recv", MIN));
    }
    if (Timer_Exists("Schedule Wait Send")) {
      printf(
          "      Wait Send                     : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("Schedule Wait Send", AVERAGE),
          Timer_Value("Schedule Wait Send", MAX),
          Timer_Value("Schedule Wait Send", MIN));
    }
    if (Timer_Exists("Schedule Unpack Incoming")) {
      printf(
          "      Unpack Incoming               : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("Schedule Unpack Incoming", AVERAGE),
          Timer_Value("Schedule Unpack Incoming", MAX),
          Timer_Value("Schedule Unpack Incoming", MIN));
    }
      if (g_dnDiagReady && g_dnDiagAvg[DN_DIAG_STEPS] > 0.0) {
        const double gib = 1024.0 * 1024.0 * 1024.0;
        const double mib = 1024.0 * 1024.0;
        const double avgSteps = g_dnDiagAvg[DN_DIAG_STEPS];
        const double sendGiBAvg = g_dnDiagAvg[DN_DIAG_SEND_BYTES] / gib;
        const double sendGiBMax = g_dnDiagMax[DN_DIAG_SEND_BYTES] / gib;
        const double sendGiBMin = g_dnDiagMin[DN_DIAG_SEND_BYTES] / gib;
        const double recvGiBAvg = g_dnDiagAvg[DN_DIAG_RECV_BYTES] / gib;
        const double recvGiBMax = g_dnDiagMax[DN_DIAG_RECV_BYTES] / gib;
        const double recvGiBMin = g_dnDiagMin[DN_DIAG_RECV_BYTES] / gib;
        const double sendSkew =
          (g_dnDiagAvg[DN_DIAG_SEND_BYTES] > 0.0)
            ? (g_dnDiagMax[DN_DIAG_SEND_BYTES] /
             g_dnDiagAvg[DN_DIAG_SEND_BYTES])
            : 0.0;
        const double recvSkew =
          (g_dnDiagAvg[DN_DIAG_RECV_BYTES] > 0.0)
            ? (g_dnDiagMax[DN_DIAG_RECV_BYTES] /
             g_dnDiagAvg[DN_DIAG_RECV_BYTES])
            : 0.0;
        const double depsSkew =
          (g_dnDiagAvg[DN_DIAG_ADJUST_DEPS] > 0.0)
            ? (g_dnDiagMax[DN_DIAG_ADJUST_DEPS] /
             g_dnDiagAvg[DN_DIAG_ADJUST_DEPS])
            : 0.0;

        printf("    Dangling comm/load diagnostics (solver loop cumulative)\n");
        printf(
          "      DN Send Volume                : %.3f (Average)   %.3f (Max) %.3f "
          "(Min) GiB  [skew %.2f]\n",
          sendGiBAvg, sendGiBMax, sendGiBMin, sendSkew);
        printf(
          "      DN Recv Volume                : %.3f (Average)   %.3f (Max) %.3f "
          "(Min) GiB  [skew %.2f]\n",
          recvGiBAvg, recvGiBMax, recvGiBMin, recvSkew);
        printf(
          "      DN Send Msg Count             : %.1f (Average)   %.1f (Max) %.1f "
          "(Min)\n",
          g_dnDiagAvg[DN_DIAG_SEND_MSGS], g_dnDiagMax[DN_DIAG_SEND_MSGS],
          g_dnDiagMin[DN_DIAG_SEND_MSGS]);
        printf(
          "      DN Recv Msg Count             : %.1f (Average)   %.1f (Max) %.1f "
          "(Min)\n",
          g_dnDiagAvg[DN_DIAG_RECV_MSGS], g_dnDiagMax[DN_DIAG_RECV_MSGS],
          g_dnDiagMin[DN_DIAG_RECV_MSGS]);
        printf(
          "      DN Send Per Step              : %.3f MiB/rank (Average)\n",
          (g_dnDiagAvg[DN_DIAG_SEND_BYTES] / avgSteps) / mib);
        printf(
          "      DN Recv Per Step              : %.3f MiB/rank (Average)\n",
          (g_dnDiagAvg[DN_DIAG_RECV_BYTES] / avgSteps) / mib);
        printf(
          "      DN Adjust Deps                : %.1f (Average)   %.1f (Max) %.1f "
          "(Min)  [skew %.2f]\n",
          g_dnDiagAvg[DN_DIAG_ADJUST_DEPS], g_dnDiagMax[DN_DIAG_ADJUST_DEPS],
          g_dnDiagMin[DN_DIAG_ADJUST_DEPS], depsSkew);
        printf(
          "      DN Node Count                 : %.1f (Average)   %.1f (Max) %.1f "
          "(Min)\n",
          g_dnDiagAvg[DN_DIAG_ADJUST_NODES],
          g_dnDiagMax[DN_DIAG_ADJUST_NODES],
          g_dnDiagMin[DN_DIAG_ADJUST_NODES]);

          printf("    Wait Recv by schedule stage (solver loop cumulative)\n");
          for (int i = 0; i < WAIT_RECV_STAGE_COUNT; i++) {
          printf(
            "      %-29s : %.2f (Average)   %.2f (Max, rank %d) %.2f (Min) "
            "seconds\n",
            g_waitRecvStageNames[i], g_waitRecvStageAvg[i],
            g_waitRecvStageMax[i], g_waitRecvStageMaxRank[i],
            g_waitRecvStageMin[i]);
          }

          printf("    Per-stage schedule timing diagnostics\n");
          for (int i = 0; i < WAIT_RECV_STAGE_COUNT; i++) {
            printf("      %s\n", g_waitRecvStageNames[i]);
            for (int metric = 0; metric < STAGE_TIME_METRIC_COUNT_PER_STAGE;
                 metric++) {
              const int32_t idx = stage_time_index(i, metric);
              printf(
                "        %-15s : %.2f (Average)   %.2f (Max) %.2f (Min) seconds\n",
                g_stageTimeMetricNames[metric], g_stageTimeAvg[idx],
                g_stageTimeMax[idx], g_stageTimeMin[idx]);
            }
          }

            printf("    Per-stage comm volume/message diagnostics\n");
            for (int i = 0; i < WAIT_RECV_STAGE_COUNT; i++) {
            const int32_t sendBytesI = stage_comm_index(i, STAGE_COMM_SEND_BYTES);
            const int32_t recvBytesI = stage_comm_index(i, STAGE_COMM_RECV_BYTES);
            const int32_t sendMsgsI = stage_comm_index(i, STAGE_COMM_SEND_MSGS);
            const int32_t recvMsgsI = stage_comm_index(i, STAGE_COMM_RECV_MSGS);
            const int32_t sendNodesI = stage_comm_index(i, STAGE_COMM_SEND_NODES);
            const int32_t recvNodesI = stage_comm_index(i, STAGE_COMM_RECV_NODES);
            const double mib = 1024.0 * 1024.0;
            const double steps = g_dnDiagAvg[DN_DIAG_STEPS] > 0.0
                         ? g_dnDiagAvg[DN_DIAG_STEPS]
                         : 1.0;
            printf(
              "      %-29s : send %.3f MiB/step, recv %.3f MiB/step, msgs %.1f/%.1f per step\n",
              g_waitRecvStageNames[i],
              (g_stageCommAvg[sendBytesI] / steps) / mib,
              (g_stageCommAvg[recvBytesI] / steps) / mib,
              g_stageCommAvg[sendMsgsI] / steps, g_stageCommAvg[recvMsgsI] / steps);
            printf(
              "                                 nodes %.1f/%.1f per step, send skew %.2f, recv skew %.2f\n",
              g_stageCommAvg[sendNodesI] / steps, g_stageCommAvg[recvNodesI] / steps,
              g_stageCommAvg[sendBytesI] > 0.0
                ? (g_stageCommMax[sendBytesI] / g_stageCommAvg[sendBytesI])
                : 0.0,
              g_stageCommAvg[recvBytesI] > 0.0
                ? (g_stageCommMax[recvBytesI] / g_stageCommAvg[recvBytesI])
                : 0.0);
            }
	    }
	  }
  if (Timer_Exists("GPU Memory Copy") || Timer_Exists("GPU Load Disp") ||
      Timer_Exists("GPU Unload Forces") || Timer_Exists("GPU Load Forces") ||
      Timer_Exists("GPU Unload Disp") || Timer_Exists("GPU Patch Disp") ||
      Timer_Exists("GPU Reset")) {
    printf("    GPU copy/enqueue internals\n");
    if (Timer_Exists("GPU Memory Copy")) {
      printf(
          "      GPU Memory Copy               : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("GPU Memory Copy", AVERAGE),
          Timer_Value("GPU Memory Copy", MAX),
          Timer_Value("GPU Memory Copy", MIN));
    }
    if (Timer_Exists("GPU Load Disp")) {
      printf(
          "      Load Disp                     : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("GPU Load Disp", AVERAGE),
          Timer_Value("GPU Load Disp", MAX),
          Timer_Value("GPU Load Disp", MIN));
    }
    if (Timer_Exists("GPU Unload Forces")) {
      printf(
          "      Unload Forces                 : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("GPU Unload Forces", AVERAGE),
          Timer_Value("GPU Unload Forces", MAX),
          Timer_Value("GPU Unload Forces", MIN));
    }
    if (Timer_Exists("GPU Load Forces")) {
      printf(
          "      Load Forces                   : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("GPU Load Forces", AVERAGE),
          Timer_Value("GPU Load Forces", MAX),
          Timer_Value("GPU Load Forces", MIN));
    }
    if (Timer_Exists("GPU Unload Disp")) {
      printf(
          "      Unload Disp                   : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("GPU Unload Disp", AVERAGE),
          Timer_Value("GPU Unload Disp", MAX),
          Timer_Value("GPU Unload Disp", MIN));
    }
    if (Timer_Exists("GPU Patch Disp")) {
      printf(
          "      Patch Disp                    : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("GPU Patch Disp", AVERAGE),
          Timer_Value("GPU Patch Disp", MAX),
          Timer_Value("GPU Patch Disp", MIN));
    }
    if (Timer_Exists("GPU Reset")) {
      printf(
          "      Reset                         : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("GPU Reset", AVERAGE),
          Timer_Value("GPU Reset", MAX),
          Timer_Value("GPU Reset", MIN));
    }
  }
  if (Timer_Exists("Adjust Cache Build") ||
      Timer_Exists("Adjust Distribution Core") ||
      Timer_Exists("Adjust Assignment Core")) {
    printf("    Compute adjust internals\n");
    if (Timer_Exists("Adjust Cache Build")) {
      printf(
          "      Cache Build                   : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("Adjust Cache Build", AVERAGE),
          Timer_Value("Adjust Cache Build", MAX),
          Timer_Value("Adjust Cache Build", MIN));
    }
    if (Timer_Exists("Adjust Distribution Core")) {
      printf(
          "      Distribution Core             : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("Adjust Distribution Core", AVERAGE),
          Timer_Value("Adjust Distribution Core", MAX),
          Timer_Value("Adjust Distribution Core", MIN));
    }
    if (Timer_Exists("Adjust Assignment Core")) {
      printf(
          "      Assignment Core               : %.2f (Average)   %.2f (Max) %.2f "
          "(Min) seconds\n",
          Timer_Value("Adjust Assignment Core", AVERAGE),
          Timer_Value("Adjust Assignment Core", MAX),
          Timer_Value("Adjust Assignment Core", MIN));
    }
  }
  printf("    IO\n");

  if (Timer_Exists("Solver drm output")) {
    printf(
        "        Drm Output                  : %.2f (Average)   %.2f (Max) "
        "%.2f (Min) seconds\n",
        Timer_Value("Solver drm output", AVERAGE),
        Timer_Value("Solver drm output", MAX),
        Timer_Value("Solver drm output", MIN));
  }
  if (Timer_Exists("Solver drm read displacements")) {
    printf(
        "        Solver drm read disp        : %.2f (Average)   %.2f (Max) "
        "%.2f (Min) seconds\n",
        Timer_Value("Solver drm read displacements", AVERAGE),
        Timer_Value("Solver drm read displacements", MAX),
        Timer_Value("Solver drm read displacements", MIN));
  }
  printf(
      "        Solver Stats Print          : %.2f (Average)   %.2f (Max) %.2f "
      "(Min) seconds\n",
      Timer_Value("Solver Stats Print", AVERAGE),
      Timer_Value("Solver Stats Print", MAX),
      Timer_Value("Solver Stats Print", MIN));
  if (Timer_Exists("Print Planes"))
    printf(
        "        Planes                      : %.2f (Average)   %.2f (Max) "
        "%.2f (Min) seconds\n",
        Timer_Value("Print Planes", AVERAGE), Timer_Value("Print Planes", MAX),
        Timer_Value("Print Planes", MIN));
  if (Timer_Exists("Print Stations"))
    printf(
        "        Stations                    : %.2f (Average)   %.2f (Max) "
        "%.2f (Min) seconds\n",
        Timer_Value("Print Stations", AVERAGE),
        Timer_Value("Print Stations", MAX), Timer_Value("Print Stations", MIN));
  if (Timer_Exists("Checkpoint"))
    printf("        Checkpoint                  : %.2f\n",
           Timer_Value("Checkpoint", (TimerKind)0));
  printf("\n");
  printf("TOTAL WALL CLOCK                    : %.2f seconds\n",
         Timer_Value("Total Wall Clock", (TimerKind)0));
  printf("\n");
  printf("\n____________Analysis_____________\n");
  printf(
      "Solver I/O                 : %.2f (Average)   %.2f (Max) %.2f (Min) "
      "seconds\n",
      Timer_Value("Solver I/O", AVERAGE), Timer_Value("Solver I/O", MAX),
      Timer_Value("Solver I/O", MIN));
  printf(
      "Solver Compute             : %.2f (Average)   %.2f (Max) %.2f (Min) "
      "seconds\n",
      Timer_Value("Compute Physics", AVERAGE),
      Timer_Value("Compute Physics", MAX), Timer_Value("Compute Physics", MIN));
  printf(
      "Solver Communicate         : %.2f (Average)   %.2f (Max) %.2f (Min) "
      "seconds\n",
      Timer_Value("Communication", AVERAGE), Timer_Value("Communication", MAX),
      Timer_Value("Communication", MIN));
  printf("Compute/Communicate Ratio  : %.2f \n",
         Timer_Value("Compute Physics", AVERAGE) /
             Timer_Value("Communication", AVERAGE));
  printf("\n\n\n\n");

  fflush(stdout);

  return;
}

/**
 * Prepare data to compute the myForce vector.  It calls compute_force.
 */
static void source_init(const char* physicsin) {
  if ((Param.drmImplement == NO) ||
      (Param.drmImplement == YES && Param.theDrmPart == PART1)) {
    double globalDelayT = 0;
    double surfaceShift = 0;

    /* Load to Global.theMPIInformation */
    Global.theMPIInformation.myid = Global.myID;
    Global.theMPIInformation.groupsize = Global.theGroupSize;

    /* Load to theNumericsInforamation */
    Global.theNumericsInformation.numberoftimesteps = Param.theTotalSteps;
    Global.theNumericsInformation.deltat = Param.theDeltaT;
    Global.theNumericsInformation.validfrequency = Param.theFreq;
    Global.theNumericsInformation.endT = Param.theEndT;
    Global.theNumericsInformation.startT = Param.theStartT;

    Global.theNumericsInformation.xlength = Param.theDomainX;
    Global.theNumericsInformation.ylength = Param.theDomainY;
    Global.theNumericsInformation.zlength = Param.theDomainZ;

    if (Param.includeNonlinearAnalysis == YES) {
      globalDelayT = get_geostatic_total_time();
    }

    if (Param.includeBuildings == YES) {
      surfaceShift = get_surface_shift();
    }

    // TODO: check that this shifting only has effect here
    //  if ( Param.includeTopography == YES ) {
    //      surfaceShift = get_thebase_topo();
    //  }

    /* it will create the files to be read each time step to
   load (force) the mesh */
    if (compute_print_source(
            physicsin, Param.theSourceOutputDir, Global.myOctree, Global.myMesh,
            Global.theNumericsInformation, Global.theMPIInformation,
            globalDelayT, surfaceShift, &Param.utmZone)) {
      monitor_print("Err:cannot create source forces");
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      exit(1);
    }

    Global.theNodesLoaded = source_get_local_loaded_nodes_count();

    if (Global.theNodesLoaded != 0) {
      size_t ret;
      int32_t node_count;

      Global.fpsource = source_open_forces_file("r");

      hu_fread(&node_count, sizeof(int32_t), 1, Global.fpsource);

      /* \todo add assertion for node_count == Global.theNodesLoaded */

      Global.theNodesLoadedList =
          (int32_t*)malloc(sizeof(int32_t) * Global.theNodesLoaded);
      Global.myForces =
          (vector3D_t*)calloc(Global.theNodesLoaded, sizeof(vector3D_t));

      if (Global.myForces == NULL || Global.theNodesLoadedList == NULL) {
        solver_abort("source_init", "memory allocation failed",
                     "Cannot allocate memory for Global.myForces or "
                     "loaded nodes list arrays\n");
      }

      ret = hu_fread(Global.theNodesLoadedList, sizeof(int32_t),
                     Global.theNodesLoaded, Global.fpsource);

      if (ret != Global.theNodesLoaded) {
        solver_abort("source_init(", "fread failed",
                     "Could not read nodal force file");
      }
    }
  }
}

/**
 * Search a point in the domain of the local mesh.
 *
 *   input: coordinates
 *  output: 0 fail 1 success
 */
int32_t search_point(vector3D_t point, octant_t** octant) {
  tick_t xTick, yTick, zTick;

  xTick = point.x[0] / Global.myMesh->ticksize;
  yTick = point.x[1] / Global.myMesh->ticksize;
  zTick = point.x[2] / Global.myMesh->ticksize;

  *octant = octor_searchoctant(Global.myOctree, xTick, yTick, zTick, PIXELLEVEL,
                               AGGREGATE_SEARCH);

  if ((*octant == NULL) || ((*octant)->where == REMOTE)) {
    return 0;
  }

  return 1;
}

static int32_t mesh_element_index_for_octant(octant_t* octant,
                                             const char* context) {
  int32_t eindex = -1;

  if (octant->appdata != NULL) {
    int32_t cached = (int32_t)((intptr_t)octant->appdata - 1);
    if ((cached >= 0) && (cached < Global.myMesh->lenum)) {
      eindex = cached;
    }
  }

  if (eindex < 0) {
    for (eindex = 0; eindex < Global.myMesh->lenum; eindex++) {
      int32_t lnid0 = Global.myMesh->elemTable[eindex].lnid[0];

      if ((Global.myMesh->nodeTable[lnid0].x == octant->lx) &&
          (Global.myMesh->nodeTable[lnid0].y == octant->ly) &&
          (Global.myMesh->nodeTable[lnid0].z == octant->lz)) {
        break;
      }
    }
  }

  if ((eindex < 0) || (eindex == Global.myMesh->lenum)) {
    fprintf(stderr, "Thread %d: %s: ", Global.myID, context);
    fprintf(stderr, "No element matches the containing octant.\n");
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  if (Global.myMesh->elemTable[eindex].level != octant->level) {
    fprintf(stderr, "Thread %d: %s: internal error\n", Global.myID, context);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  return eindex;
}

int32_t solver_get_octant_material(octant_t* octant, edata_t* material) {
  int32_t eindex;

  if ((octant == NULL) || (material == NULL)) {
    return -1;
  }

  eindex = mesh_element_index_for_octant(octant, __FUNCTION_NAME);
  if ((eindex < 0) || (Global.myMesh == NULL) ||
      (Global.myMesh->elemTable[eindex].data == NULL)) {
    return -1;
  }

  memcpy(material, Global.myMesh->elemTable[eindex].data, sizeof(edata_t));
  return eindex;
}

/**
 * \param octant where the point is located.
 * \param pointcoords coordinates of the point.
 * \param localcoords[out] the displacment.
 */
extern int compute_csi_eta_dzeta(octant_t* octant, vector3D_t pointcoords,
                                 vector3D_t* localcoords,
                                 int32_t* localNodeID) {
  tick_t edgeticks;
  int32_t eindex;
  double center_x, center_y, center_z;

  /* various convienient variables */
  double xGlobal = pointcoords.x[0];
  double yGlobal = pointcoords.x[1];
  double zGlobal = pointcoords.x[2];
  double h;

  edgeticks = (tick_t)1 << (PIXELLEVEL - octant->level);
  h = Global.myMesh->ticksize * edgeticks;

  /* Calculate the center coordinate of the element */
  center_x = Global.myMesh->ticksize * (octant->lx + edgeticks / 2);
  center_y = Global.myMesh->ticksize * (octant->ly + edgeticks / 2);
  center_z = Global.myMesh->ticksize * (octant->lz + edgeticks / 2);

  eindex = mesh_element_index_for_octant(octant, __FUNCTION_NAME);

  memcpy(localNodeID, Global.myMesh->elemTable[eindex].lnid,
         sizeof(int32_t) * 8);

  localcoords->x[0] = 2 * (xGlobal - center_x) / h;
  localcoords->x[1] = 2 * (yGlobal - center_y) / h;
  localcoords->x[2] = 2 * (zGlobal - center_z) / h;

  /* redefine local coordinates and nodes to interpolate if VT is on */

  // FILE  *topoinfo = hu_fopen( "topoElem.txt", "w" );

  if ((Param.includeTopography == YES) &&
      (isTopoElement(Global.myMesh, eindex, 0)) && (get_topo_meth() == VT) &&
      (Param.drmImplement == NO)) {
    elem_t* elemp;
    // edata_t *edata;

    elemp = &Global.myMesh->elemTable[eindex];
    // edata = (edata_t *)elemp->data;

    /* Calculate the element's origin */
    double xo = Global.myMesh->ticksize * octant->lx;
    double yo = Global.myMesh->ticksize * octant->ly;
    double zo = Global.myMesh->ticksize * octant->lz;

    double aux[3] = {0};
    int cube_part = get_cube_partition(eindex);

    compute_tetra_localcoord(pointcoords, elemp, localNodeID, aux, xo, yo, zo,
                             h, cube_part);

    localcoords->x[0] = aux[0];
    localcoords->x[1] = aux[1];
    localcoords->x[2] = aux[2];

    *(localNodeID + 4) = 0;
    *(localNodeID + 5) = 0;
    *(localNodeID + 6) = 0;
    *(localNodeID + 7) = 0;
  }

  return 1;
}

/**
 * Read stations info.  This is called by PE 0.
 */
static void read_stations_info(const char* numericalin) {
  static const char* fname = __FUNCTION_NAME;

  int iStation;
  double lon, lat, depth, *auxiliar;
  FILE* fp;

  vector3D_t coords;

  /* obtain the stations specifications */
  if ((fp = fopen(numericalin, "r")) == NULL) {
    solver_abort(fname, numericalin,
                 "Error opening numerical.in configuration file");
  }

  if (parsetext(fp, "output_stations_print_rate", 'i',
                &Param.theStationsPrintRate) != 0) {
    solver_abort(fname, NULL,
                 "Error parsing output_planes_print_rate field from %s\n",
                 numericalin);
  }

  auxiliar = (double*)malloc(sizeof(double) * Param.theNumberOfStations * 3);
  Param.theStationX =
      (double*)malloc(sizeof(double) * Param.theNumberOfStations);
  Param.theStationY =
      (double*)malloc(sizeof(double) * Param.theNumberOfStations);
  Param.theStationZ =
      (double*)malloc(sizeof(double) * Param.theNumberOfStations);

  if (Param.theStationX == NULL || Param.theStationY == NULL ||
      Param.theStationZ == NULL) {
    fprintf(stdout, "Err alloc theStations arrays in output_stations_init");
    fflush(stdout);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
  }

  if (parsedarray(fp, "output_stations", Param.theNumberOfStations * 3,
                  auxiliar) != 0) {
    solver_abort(fname, NULL, "Err parsing output_stations from %s\n",
                 numericalin);
  }

  for (iStation = 0; iStation < Param.theNumberOfStations; iStation++) {
    lat = auxiliar[iStation * 3];
    lon = auxiliar[iStation * 3 + 1];
    depth = auxiliar[iStation * 3 + 2];
    coords = compute_domain_coords_linearinterp(
        lon, lat, Param.theSurfaceCornersLong, Param.theSurfaceCornersLat,
        Param.theDomainY, Param.theDomainX, &Param.utmZone);
    Param.theStationX[iStation] = coords.x[0];
    Param.theStationY[iStation] = coords.x[1];
    Param.theStationZ[iStation] = depth;

    if (Param.includeBuildings == YES) {
      Param.theStationZ[iStation] += get_surface_shift();
    }
    if (Param.includeTopography == YES) {
      //      Param.theStationZ [ iStation ] += get_thebase_topo();
      /* Dorian: place stations on surface */
      Param.theStationZ[iStation] += point_elevation(coords.x[0], coords.x[1]);
    }
  }

  free(auxiliar);

  parsetext(fp, "output_stations_directory", 's', Param.theStationsDirOut);
  createDir(Param.theStationsDirOut, fname);
  Param.theStationMediumsFile[0] = '\0';
  parsetext(fp, "station_mediums_output_file", 's',
            Param.theStationMediumsFile);
  if (Param.theStationMediumsFile[0] == '\0') {
    snprintf(Param.theStationMediumsFile, sizeof(Param.theStationMediumsFile),
             "%s/station_mediums.txt", Param.theStationsDirOut);
  }
  create_parent_dir_for_file(Param.theStationMediumsFile, fname);

  fclose(fp);
  return;
}

/**
 * Broadcast info about the stations.
 */
void broadcast_stations_info() {
  /*initialize the local structures */
  if (Global.myID != 0) {
    Param.theStationX =
        (double*)malloc(sizeof(double) * Param.theNumberOfStations);
    Param.theStationY =
        (double*)malloc(sizeof(double) * Param.theNumberOfStations);
    Param.theStationZ =
        (double*)malloc(sizeof(double) * Param.theNumberOfStations);

    if (Param.theStationX == NULL || Param.theStationY == NULL ||
        Param.theStationZ == NULL) {
      solver_abort("broadcast_stations_info", NULL,
                   "Error: Unable to create stations arrays");
    }
  }

  MPI_Barrier(comm_solver);

  MPI_Bcast(Param.theStationX, Param.theNumberOfStations, MPI_DOUBLE, 0,
            comm_solver);
  MPI_Bcast(Param.theStationY, Param.theNumberOfStations, MPI_DOUBLE, 0,
            comm_solver);
  MPI_Bcast(Param.theStationZ, Param.theNumberOfStations, MPI_DOUBLE, 0,
            comm_solver);
  MPI_Bcast(&Param.theStationsPrintRate, 1, MPI_INT, 0, comm_solver);

  broadcast_char_array(Param.theStationsDirOut, sizeof(Param.theStationsDirOut),
                       0, comm_solver);
  broadcast_char_array(Param.theStationMediumsFile,
                       sizeof(Param.theStationMediumsFile), 0, comm_solver);

  return;
}

static int station_medium_record_compare(const void* a, const void* b) {
  const station_medium_record_t* lhs = (const station_medium_record_t*)a;
  const station_medium_record_t* rhs = (const station_medium_record_t*)b;

  if (lhs->stationId < rhs->stationId) {
    return -1;
  }
  if (lhs->stationId > rhs->stationId) {
    return 1;
  }
  return 0;
}

static void station_write_mediums_file(
    const station_medium_record_t* localRecords) {
  const int32_t localCount = Param.myNumberOfStations;
  const int64_t localBytes64 =
      (int64_t)localCount * (int64_t)sizeof(station_medium_record_t);
  station_medium_record_t emptyRecord;
  int localBytes;
  int* recvCounts = NULL;
  int* displs = NULL;
  station_medium_record_t* records = NULL;
  int totalBytes = 0;
  int totalRecords = 0;
  void* sendBuffer;

  if (localBytes64 > INT_MAX) {
    solver_abort(__FUNCTION_NAME, NULL,
                 "Too many local station medium records to gather");
  }
  localBytes = (int)localBytes64;
  sendBuffer = (localBytes > 0) ? (void*)localRecords : &emptyRecord;

  if (Global.myID == 0) {
    recvCounts = (int*)calloc((size_t)Global.theGroupSize, sizeof(int));
    displs = (int*)calloc((size_t)Global.theGroupSize, sizeof(int));
    if ((recvCounts == NULL) || (displs == NULL)) {
      free(recvCounts);
      free(displs);
      solver_abort(__FUNCTION_NAME, NULL,
                   "Error allocating station medium gather buffers");
    }
  }

  MPI_Gather((void*)&localBytes, 1, MPI_INT, recvCounts, 1, MPI_INT, 0,
             comm_solver);

  if (Global.myID == 0) {
    for (int rank = 0; rank < Global.theGroupSize; rank++) {
      displs[rank] = totalBytes;
      totalBytes += recvCounts[rank];
    }
    totalRecords = totalBytes / (int)sizeof(station_medium_record_t);
    if (totalRecords > 0) {
      records = (station_medium_record_t*)malloc((size_t)totalBytes);
      if (records == NULL) {
        free(recvCounts);
        free(displs);
        solver_abort(__FUNCTION_NAME, NULL,
                     "Error allocating station medium records");
      }
    }
  }

  MPI_Gatherv(sendBuffer, localBytes, MPI_BYTE, records, recvCounts,
              displs, MPI_BYTE, 0, comm_solver);

  if (Global.myID == 0) {
    FILE* fp;

    if (totalRecords > 0) {
      qsort(records, (size_t)totalRecords, sizeof(station_medium_record_t),
            station_medium_record_compare);
    }

    fp = hu_fopen(Param.theStationMediumsFile, "w");
    fprintf(fp,
            "# station_id owner_medium pressure_defined owner_rank element_id "
            "element_level x_m y_m z_m csi eta dzeta h_m Vp_mps Vs_mps "
            "rho_kgm3\n");
    for (int i = 0; i < totalRecords; i++) {
      const station_medium_record_t* rec = &records[i];
      fprintf(fp,
              "%d %s %s %d %" PRId64 " %d %.9e %.9e %.9e %.9e %.9e %.9e "
              "%.9e %.9e %.9e %.9e\n",
              rec->stationId, rec->ownerIsAcoustic ? "acoustic" : "elastic",
              rec->pressureDefined ? "yes" : "no", rec->ownerRank,
              rec->elementId, rec->elementLevel, rec->x, rec->y, rec->z,
              rec->csi, rec->eta, rec->dzeta, rec->h, rec->Vp, rec->Vs,
              rec->rho);
    }
    fclose(fp);

    monitor_print("Station medium metadata:           %s\n",
                  Param.theStationMediumsFile);
    if (totalRecords != Param.theNumberOfStations) {
      monitor_print(
          "WARNING: wrote %d of %d configured station medium records.\n",
          totalRecords, Param.theNumberOfStations);
    }

    free(records);
    free(recvCounts);
    free(displs);
  }
}

/**
 * Prepare all the info every station needs once it is located in a processor
 */
void setup_stations_data() {
  int32_t iStation, iLCStation = 0; /* LC local count */
  int32_t localCapacity = 0;
  int32_t* localStationIds = NULL;
  octant_t** localStationOctants = NULL;
  station_medium_record_t* stationMediumRecords = NULL;
  vector3D_t stationCoords;
  octant_t* octant;

  /* Locate stations once, cache their owning octants, and reuse that for
   * per-station setup to avoid a second full search pass. */
  Timer_Start("Stations Locate");
  for (iStation = 0; iStation < Param.theNumberOfStations; iStation++) {
    stationCoords.x[0] = Param.theStationX[iStation];
    stationCoords.x[1] = Param.theStationY[iStation];
    stationCoords.x[2] = Param.theStationZ[iStation];

    if (!station_point_in_local_interval(stationCoords)) {
      continue;
    }

    if (search_point(stationCoords, &octant) == 1) {
      if (Param.myNumberOfStations == localCapacity) {
        int32_t newCapacity = (localCapacity == 0) ? 1024 : localCapacity * 2;
        int32_t* newIds =
            (int32_t*)realloc(localStationIds, sizeof(int32_t) * newCapacity);
        octant_t** newOctants = (octant_t**)realloc(
            localStationOctants, sizeof(octant_t*) * newCapacity);

        if ((newIds == NULL) || (newOctants == NULL)) {
          free(newIds);
          free(newOctants);
          free(localStationIds);
          free(localStationOctants);
          solver_abort(__FUNCTION_NAME, NULL,
                       "Error allocating local station setup buffers");
        }

        localStationIds = newIds;
        localStationOctants = newOctants;
        localCapacity = newCapacity;
      }

      localStationIds[Param.myNumberOfStations] = iStation;
      localStationOctants[Param.myNumberOfStations] = octant;
      Param.myNumberOfStations++;
    }
  }
  Timer_Stop("Stations Locate");

  /* allocate memory if necessary and generate the list of stations per
   * processor */
  if (Param.myNumberOfStations != 0) {
    //  monitor_print( "PE=%d local_stations=%d total_station_count=%d\n",
    //             Global.myID, Param.myNumberOfStations,
    //             Param.theNumberOfStations );

    XMALLOC_VAR_N(Param.myStations, station_t, Param.myNumberOfStations);
    stationMediumRecords = (station_medium_record_t*)calloc(
        (size_t)Param.myNumberOfStations, sizeof(station_medium_record_t));
    if (stationMediumRecords == NULL) {
      solver_abort(__FUNCTION_NAME, NULL,
                   "Error allocating station medium records");
    }
  }

  Timer_Start("Stations Local Init");
#ifdef _OPENMP
  const int32_t useOmp = station_output_use_omp();
#pragma omp parallel for if(useOmp) schedule(static) \
    private(iStation, octant, stationCoords)
#endif
  for (iLCStation = 0; iLCStation < Param.myNumberOfStations; iLCStation++) {
    char stationFile[256];
    int32_t stationElementIndex;
    edata_t* stationMaterial;

    iStation = localStationIds[iLCStation];
    octant = localStationOctants[iLCStation];
    stationCoords.x[0] = Param.theStationX[iStation];
    stationCoords.x[1] = Param.theStationY[iStation];
    stationCoords.x[2] = Param.theStationZ[iStation];

    Param.myStations[iLCStation].id = iStation;
    Param.myStations[iLCStation].coords = stationCoords;
    stationElementIndex =
        mesh_element_index_for_octant(octant, __FUNCTION_NAME);
    stationMaterial =
        (edata_t*)Global.myMesh->elemTable[stationElementIndex].data;
    Param.myStations[iLCStation].pressure_is_defined =
        (acoustic_pressure_enabled() && element_is_acoustic(stationMaterial))
            ? 1
            : 0;
    Param.myStations[iLCStation].acoustic_kinematics_ready = 0;
    memset(Param.myStations[iLCStation].acoustic_disp, 0,
           sizeof(Param.myStations[iLCStation].acoustic_disp));
    memset(Param.myStations[iLCStation].acoustic_vel, 0,
           sizeof(Param.myStations[iLCStation].acoustic_vel));
    memset(Param.myStations[iLCStation].acoustic_acc, 0,
           sizeof(Param.myStations[iLCStation].acoustic_acc));
    memset(Param.myStations[iLCStation].acoustic_prev_acc, 0,
           sizeof(Param.myStations[iLCStation].acoustic_prev_acc));
    Param.myStations[iLCStation].h =
        (stationMaterial != NULL) ? stationMaterial->edgesize : 0.0;
    Param.myStations[iLCStation].rho =
        (stationMaterial != NULL) ? stationMaterial->rho : 0.0;
    sprintf(stationFile, "%s/station.%d", Param.theStationsDirOut, iStation);
    Param.myStations[iLCStation].fpoutputfile = hu_fopen(stationFile, "w");
    Param.myStations[iLCStation].prev_disp[0] = 0.0;
    Param.myStations[iLCStation].prev_disp[1] = 0.0;
    Param.myStations[iLCStation].prev_disp[2] = 0.0;
    compute_csi_eta_dzeta(octant, Param.myStations[iLCStation].coords,
                          &(Param.myStations[iLCStation].localcoords),
                          Param.myStations[iLCStation].nodestointerpolate);

    stationMediumRecords[iLCStation].stationId = iStation;
    stationMediumRecords[iLCStation].ownerRank = Global.myID;
    stationMediumRecords[iLCStation].ownerIsAcoustic =
        element_is_acoustic(stationMaterial) ? 1 : 0;
    stationMediumRecords[iLCStation].pressureDefined =
        Param.myStations[iLCStation].pressure_is_defined;
    stationMediumRecords[iLCStation].elementLevel =
        Global.myMesh->elemTable[stationElementIndex].level;
    stationMediumRecords[iLCStation].elementId =
        Global.myMesh->elemTable[stationElementIndex].geid;
    stationMediumRecords[iLCStation].x = stationCoords.x[0];
    stationMediumRecords[iLCStation].y = stationCoords.x[1];
    stationMediumRecords[iLCStation].z = stationCoords.x[2];
    stationMediumRecords[iLCStation].csi =
        Param.myStations[iLCStation].localcoords.x[0];
    stationMediumRecords[iLCStation].eta =
        Param.myStations[iLCStation].localcoords.x[1];
    stationMediumRecords[iLCStation].dzeta =
        Param.myStations[iLCStation].localcoords.x[2];
    stationMediumRecords[iLCStation].h =
        (stationMaterial != NULL) ? stationMaterial->edgesize : 0.0;
    stationMediumRecords[iLCStation].Vp =
        (stationMaterial != NULL) ? stationMaterial->Vp : 0.0;
    stationMediumRecords[iLCStation].Vs =
        (stationMaterial != NULL) ? stationMaterial->Vs : 0.0;
    stationMediumRecords[iLCStation].rho =
        (stationMaterial != NULL) ? stationMaterial->rho : 0.0;

    /*
     * This section now enclosed in a DEBUG def because its information
     * is only useful for debugging purposes and otherwise it just
     * introduce noise to the post-processing.
     */
#ifdef DEBUG
    fprintf(Param.myStations[iLCStation].fpoutputfile,
            "# Node identifiers:\n");

    int iNode;

    for (iNode = 0; iNode < 8; iNode++) {
      fprintf(Param.myStations[iLCStation].fpoutputfile,
              "#  %13" INT64_FMT "\n",
              Global.myMesh
                  ->nodeTable[Param.myStations[iLCStation]
                                  .nodestointerpolate[iNode]]
                  .gnid);
    }
#endif /* DEBUG */

    /*
     * This is a one line heading for the station files. Spacing is
     * such that it aligns with data.
     *
     * The lines after X and Y represent the actual orientation of the
     * axes if one looks at the domain's surface on a screen. Zeta's dot
     * means the Z axis goes in(+) and out(-) of the screen.
     */

    fputs("#  Time(s)         X|(m)         Y-(m)         Z.(m)",
          Param.myStations[iLCStation].fpoutputfile);

    if ((Param.printStationVelocities == YES) ||
        (Param.printStationAccelerations == YES)) {
      fputs("       X|(m/s)       Y-(m/s)       Z.(m/s)",
            Param.myStations[iLCStation].fpoutputfile);
    }

    if (Param.printStationAccelerations == YES) {
      fputs("      X|(m/s2)      Y-(m/s2)      Z.(m/s2)",
            Param.myStations[iLCStation].fpoutputfile);
    }

    if (Param.printStationPressures == YES) {
      fputs("      Pressure(Pa)",
            Param.myStations[iLCStation].fpoutputfile);
    }

    /*
     * Additional headings for nonlinear data.
     */
    if (Param.includeNonlinearAnalysis == YES) {
      fputs(
          "    Epsilon_XX      Sigma_XX"
          "    Epsilon_YY      Sigma_YY"
          "    Epsilon_ZZ      Sigma_ZZ"
          "    Epsilon_KK      Sigma_KK"
          "    Epsilon_XY      Sigma_XY"
          "    Epsilon_YZ      Sigma_YZ"
          "    Epsilon_XZ      Sigma_XZ"
          "            Fs         Kappa"
          "       Gamma1d         Tao1d"
          "         GGmax",
          Param.myStations[iLCStation].fpoutputfile);
    }
  }
  Timer_Stop("Stations Local Init");

  station_write_mediums_file(stationMediumRecords);
  free(stationMediumRecords);

  free(localStationIds);
  free(localStationOctants);

  free(Param.theStationX);
  free(Param.theStationY);
  free(Param.theStationZ);
}

static void solver_snapshot_station_accel_history(mysolver_t* solver) {
  int iPhi;
  int32_t iStation;
  double xi[3][8] = {{-1, 1, -1, 1, -1, 1, -1, 1},
                     {-1, -1, 1, 1, -1, -1, 1, 1},
                     {-1, -1, -1, -1, 1, 1, 1, 1}};

  if ((Param.printStationAccelerations != YES) ||
      (Param.myNumberOfStations == 0)) {
    return;
  }

#ifdef _OPENMP
  const int32_t useOmp = station_output_use_omp();
#pragma omp parallel for if(useOmp) schedule(static) private(iPhi)
#endif
  for (iStation = 0; iStation < Param.myNumberOfStations; iStation++) {
    int32_t nodesToInterpolate[8];
    vector3D_t localCoords = Param.myStations[iStation].localcoords;
    double prev_x = 0.0;
    double prev_y = 0.0;
    double prev_z = 0.0;

    for (iPhi = 0; iPhi < 8; iPhi++) {
      nodesToInterpolate[iPhi] =
          Param.myStations[iStation].nodestointerpolate[iPhi];
    }

    for (iPhi = 0; iPhi < 8; iPhi++) {
      double phi = (1 + xi[0][iPhi] * localCoords.x[0]) *
                   (1 + xi[1][iPhi] * localCoords.x[1]) *
                   (1 + xi[2][iPhi] * localCoords.x[2]) / 8;

      prev_x += phi * solver->tm2[nodesToInterpolate[iPhi]].f[0];
      prev_y += phi * solver->tm2[nodesToInterpolate[iPhi]].f[1];
      prev_z += phi * solver->tm2[nodesToInterpolate[iPhi]].f[2];
    }

    Param.myStations[iStation].prev_disp[0] = prev_x;
    Param.myStations[iStation].prev_disp[1] = prev_y;
    Param.myStations[iStation].prev_disp[2] = prev_z;
  }

  topography_snapshot_station_history(solver);
}

static int station_pressure_gradient(int32_t iStation, double grad[3]) {
  int iPhi;
  const double sx[8] = {-1, 1, -1, 1, -1, 1, -1, 1};
  const double sy[8] = {-1, -1, 1, 1, -1, -1, 1, 1};
  const double sz[8] = {-1, -1, -1, -1, 1, 1, 1, 1};
  const double h = Param.myStations[iStation].h;
  const vector3D_t localCoords = Param.myStations[iStation].localcoords;
  const double xi = localCoords.x[0];
  const double eta = localCoords.x[1];
  const double zeta = localCoords.x[2];
  int32_t nodesToInterpolate[8];

  grad[0] = 0.0;
  grad[1] = 0.0;
  grad[2] = 0.0;

  if ((Param.myStations[iStation].pressure_is_defined == 0) ||
      (Global.mySolver == NULL) || (Global.mySolver->pressure_tm1 == NULL) ||
      (h <= 0.0)) {
    return 0;
  }

  for (iPhi = 0; iPhi < 8; iPhi++) {
    nodesToInterpolate[iPhi] =
        Param.myStations[iStation].nodestointerpolate[iPhi];
  }

  /* velocity_potential: the field is phi; p = -dphi/dt
   * = -(phi^n - phi^{n-1})/dt. pressure mode: field is p directly. */
  const int vp_out = (Param.theAcousticVelocityPotential == 1) &&
                     (Param.theDeltaT > 0.0);
  for (iPhi = 0; iPhi < 8; iPhi++) {
    const int nd = nodesToInterpolate[iPhi];
    const double pressure =
        vp_out ? (-(double)(Global.mySolver->pressure_tm1[nd] -
                            Global.mySolver->pressure_tm2[nd]) /
                  Param.theDeltaT)
               : (double)Global.mySolver->pressure_tm1[nd];
    const double dphidx =
        sx[iPhi] * (1.0 + sy[iPhi] * eta) *
        (1.0 + sz[iPhi] * zeta) / (4.0 * h);
    const double dphidy =
        sy[iPhi] * (1.0 + sx[iPhi] * xi) *
        (1.0 + sz[iPhi] * zeta) / (4.0 * h);
    const double dphidz =
        sz[iPhi] * (1.0 + sx[iPhi] * xi) *
        (1.0 + sy[iPhi] * eta) / (4.0 * h);

    grad[0] += dphidx * pressure;
    grad[1] += dphidy * pressure;
    grad[2] += dphidz * pressure;
  }

  return 1;
}

static void solver_update_acoustic_station_kinematics(int32_t step) {
  int32_t iStation;

  if ((!acoustic_pressure_enabled()) || (Param.myNumberOfStations == 0) ||
      (Global.mySolver == NULL) || (Global.mySolver->pressure_tm1 == NULL)) {
    return;
  }

#ifdef _OPENMP
  const int32_t useOmp = station_output_use_omp();
#pragma omp parallel for if(useOmp) schedule(static)
#endif
  for (iStation = 0; iStation < Param.myNumberOfStations; iStation++) {
    double grad[3];
    double accel[3];
    station_t* station = &Param.myStations[iStation];

    if (station->pressure_is_defined == 0) {
      continue;
    }
    if ((station->rho <= 0.0) ||
        (station_pressure_gradient(iStation, grad) == 0)) {
      continue;
    }

    for (int axis = 0; axis < 3; axis++) {
      accel[axis] = -grad[axis] / station->rho;
    }

    if ((station->acoustic_kinematics_ready == 0) || (step == 0)) {
      for (int axis = 0; axis < 3; axis++) {
        station->acoustic_acc[axis] = accel[axis];
        station->acoustic_prev_acc[axis] = accel[axis];
      }
      station->acoustic_kinematics_ready = 1;
      continue;
    }

    for (int axis = 0; axis < 3; axis++) {
      const double oldVel = station->acoustic_vel[axis];
      const double oldAcc = station->acoustic_prev_acc[axis];

      station->acoustic_disp[axis] +=
          oldVel * Param.theDeltaT +
          0.5 * oldAcc * Param.theDeltaTSquared;
      station->acoustic_vel[axis] +=
          0.5 * (oldAcc + accel[axis]) * Param.theDeltaT;
      station->acoustic_acc[axis] = accel[axis];
      station->acoustic_prev_acc[axis] = accel[axis];
    }
  }
}

/**
 * Interpolate the displacements for the stations.
 */
static int interpolate_station_displacements(int32_t step) {
  int iPhi;

  /* Auxiliary array to handle shape functions in a loop */
  double xi[3][8] = {{-1, 1, -1, 1, -1, 1, -1, 1},
                     {-1, -1, 1, 1, -1, -1, 1, 1},
                     {-1, -1, -1, -1, 1, 1, 1, 1}};

  double phi[8];
  double dis_x, dis_y, dis_z;
  double vel_x, vel_y, vel_z;
  double acc_x, acc_y, acc_z;
  int32_t iStation, nodesToInterpolate[8];
  ;
  vector3D_t localCoords; /* convenient renaming */

#ifdef _OPENMP
  const int32_t useOmp = station_output_use_omp();
#pragma omp parallel for if(useOmp) schedule(static) \
    private(iPhi, nodesToInterpolate, localCoords, phi, dis_x, dis_y, dis_z, \
            vel_x, vel_y, vel_z, acc_x, acc_y, acc_z)
#endif
  for (iStation = 0; iStation < Param.myNumberOfStations; iStation++) {
    localCoords = Param.myStations[iStation].localcoords;

    for (iPhi = 0; iPhi < 8; iPhi++) {
      nodesToInterpolate[iPhi] =
          Param.myStations[iStation].nodestointerpolate[iPhi];
    }

    /* Compute interpolation function (phi) for each node and
     * load the displacements
     */
    dis_x = 0;
    dis_y = 0;
    dis_z = 0;

    vel_x = 0;
    vel_y = 0;
    vel_z = 0;

    acc_x = 0;
    acc_y = 0;
    acc_z = 0;

    double time = Param.theDeltaT * step;

    if (Param.myStations[iStation].pressure_is_defined != 0) {
      dis_x = Param.myStations[iStation].acoustic_disp[0];
      dis_y = Param.myStations[iStation].acoustic_disp[1];
      dis_z = Param.myStations[iStation].acoustic_disp[2];
      vel_x = Param.myStations[iStation].acoustic_vel[0];
      vel_y = Param.myStations[iStation].acoustic_vel[1];
      vel_z = Param.myStations[iStation].acoustic_vel[2];
      acc_x = Param.myStations[iStation].acoustic_acc[0];
      acc_y = Param.myStations[iStation].acoustic_acc[1];
      acc_z = Param.myStations[iStation].acoustic_acc[2];

      fprintf(Param.myStations[iStation].fpoutputfile,
              "\n%10.6f % 8e % 8e % 8e", time, dis_x, dis_y, dis_z);

      if ((Param.printStationVelocities == YES) ||
          (Param.printStationAccelerations == YES)) {
        fprintf(Param.myStations[iStation].fpoutputfile, " % 8e % 8e % 8e",
                vel_x, vel_y, vel_z);
      }

      if (Param.printStationAccelerations == YES) {
        fprintf(Param.myStations[iStation].fpoutputfile, " % 8e % 8e % 8e",
                acc_x, acc_y, acc_z);
      }

      if (Param.printStationPressures == YES) {
        fprintf(Param.myStations[iStation].fpoutputfile, " % 8e",
                interpolate_station_pressure_value(iStation));
      }

      if ((step % (Param.theStationsPrintRate * 10)) == 0) {
        fflush(Param.myStations[iStation].fpoutputfile);
      }
    } else if ((Param.includeTopography == YES) &&
               (compute_tetra_displ(&dis_x, &dis_y, &dis_z, &vel_x, &vel_y,
                                    &vel_z, &acc_x, &acc_y, &acc_z,
                                    Param.theDeltaT, Param.theDeltaTSquared,
                                    iStation, Global.mySolver) == 1)) {
      fprintf(Param.myStations[iStation].fpoutputfile,
              "\n%10.6f % 8e % 8e % 8e", time, dis_x, dis_y, dis_z);

      if (Param.printStationVelocities == YES) {
        fprintf(Param.myStations[iStation].fpoutputfile, " % 8e % 8e % 8e",
                vel_x, vel_y, vel_z);
      }

      if (Param.printStationAccelerations == YES) {
        fprintf(Param.myStations[iStation].fpoutputfile, " % 8e % 8e % 8e",
                acc_x, acc_y, acc_z);
      }

      if (Param.printStationPressures == YES) {
        fprintf(Param.myStations[iStation].fpoutputfile, " % 8e",
                interpolate_station_pressure_value(iStation));
      }

      if ((step % (Param.theStationsPrintRate * 10)) == 0) {
        fflush(Param.myStations[iStation].fpoutputfile);
      }
    } else {
      for (iPhi = 0; iPhi < 8; iPhi++) {
        phi[iPhi] = (1 + xi[0][iPhi] * localCoords.x[0]) *
                    (1 + xi[1][iPhi] * localCoords.x[1]) *
                    (1 + xi[2][iPhi] * localCoords.x[2]) / 8;

        dis_x +=
            phi[iPhi] * Global.mySolver->tm1[nodesToInterpolate[iPhi]].f[0];
        dis_y +=
            phi[iPhi] * Global.mySolver->tm1[nodesToInterpolate[iPhi]].f[1];
        dis_z +=
            phi[iPhi] * Global.mySolver->tm1[nodesToInterpolate[iPhi]].f[2];
      }

      /*
       * Please DO NOT CHANGE the format for printing the displacements.
       * It has to be *this* one because it goes in hand with the printing
       * format for the nonlinear information.
       */
      fprintf(Param.myStations[iStation].fpoutputfile,
              "\n%10.6f % 8e % 8e % 8e", time, dis_x, dis_y, dis_z);

      /*
       * Addition for printing velocities on the fly
       */

      if ((Param.printStationVelocities == YES) ||
          (Param.printStationAccelerations == YES)) {
        for (iPhi = 0; iPhi < 8; iPhi++) {
          phi[iPhi] = (1 + xi[0][iPhi] * localCoords.x[0]) *
                      (1 + xi[1][iPhi] * localCoords.x[1]) *
                      (1 + xi[2][iPhi] * localCoords.x[2]) / 8;

          dis_x -=
              phi[iPhi] * Global.mySolver->tm2[nodesToInterpolate[iPhi]].f[0];
          dis_y -=
              phi[iPhi] * Global.mySolver->tm2[nodesToInterpolate[iPhi]].f[1];
          dis_z -=
              phi[iPhi] * Global.mySolver->tm2[nodesToInterpolate[iPhi]].f[2];
        }

        vel_x = dis_x / Param.theDeltaT;
        vel_y = dis_y / Param.theDeltaT;
        vel_z = dis_z / Param.theDeltaT;

        fprintf(Param.myStations[iStation].fpoutputfile, " % 8e % 8e % 8e",
                vel_x, vel_y, vel_z);
      }

      /*
       * Addition for printing accelerations on the fly
       */

      if (Param.printStationAccelerations == YES) {
        for (iPhi = 0; iPhi < 8; iPhi++) {
          phi[iPhi] = (1 + xi[0][iPhi] * localCoords.x[0]) *
                      (1 + xi[1][iPhi] * localCoords.x[1]) *
                      (1 + xi[2][iPhi] * localCoords.x[2]) / 8;

          dis_x -=
              phi[iPhi] * Global.mySolver->tm2[nodesToInterpolate[iPhi]].f[0];
          dis_y -=
              phi[iPhi] * Global.mySolver->tm2[nodesToInterpolate[iPhi]].f[1];
          dis_z -=
              phi[iPhi] * Global.mySolver->tm2[nodesToInterpolate[iPhi]].f[2];
        }

        dis_x += Param.myStations[iStation].prev_disp[0];
        dis_y += Param.myStations[iStation].prev_disp[1];
        dis_z += Param.myStations[iStation].prev_disp[2];

        acc_x = dis_x / Param.theDeltaTSquared;
        acc_y = dis_y / Param.theDeltaTSquared;
        acc_z = dis_z / Param.theDeltaTSquared;

        fprintf(Param.myStations[iStation].fpoutputfile, " % 8e % 8e % 8e",
                acc_x, acc_y, acc_z);
      }

      if (Param.printStationPressures == YES) {
        fprintf(Param.myStations[iStation].fpoutputfile, " % 8e",
                interpolate_station_pressure_value(iStation));
      }

      /* TODO: Have this 10 as a parameter with a default value */
      if ((step % (Param.theStationsPrintRate * 10)) == 0) {
        fflush(Param.myStations[iStation].fpoutputfile);
      }
    }
  }

  return 1;
}

static double interpolate_station_pressure_value(int32_t iStation) {
  int iPhi;
  double xi[3][8] = {{-1, 1, -1, 1, -1, 1, -1, 1},
                     {-1, -1, 1, 1, -1, -1, 1, 1},
                     {-1, -1, -1, -1, 1, 1, 1, 1}};
  int32_t nodesToInterpolate[8];
  vector3D_t localCoords;
  double pressure = 0.0;

  if ((Param.printStationPressures != YES) || !acoustic_pressure_enabled() ||
      (Global.mySolver == NULL) || (Global.mySolver->pressure_tm1 == NULL)) {
    return 0.0;
  }
  if (Param.myStations[iStation].pressure_is_defined == 0) {
    return NAN;
  }

  localCoords = Param.myStations[iStation].localcoords;
  for (iPhi = 0; iPhi < 8; iPhi++) {
    nodesToInterpolate[iPhi] =
        Param.myStations[iStation].nodestointerpolate[iPhi];
  }

  for (iPhi = 0; iPhi < 8; iPhi++) {
    const double sh = (1 + xi[0][iPhi] * localCoords.x[0]) *
                      (1 + xi[1][iPhi] * localCoords.x[1]) *
                      (1 + xi[2][iPhi] * localCoords.x[2]) / 8;
    const int nd = nodesToInterpolate[iPhi];
    /* velocity_potential: p = -dphi/dt = -(phi^n - phi^{n-1})/dt. */
    const double val =
        ((Param.theAcousticVelocityPotential == 1) && (Param.theDeltaT > 0.0))
            ? (-(double)(Global.mySolver->pressure_tm1[nd] -
                         Global.mySolver->pressure_tm2[nd]) /
               Param.theDeltaT)
            : (double)Global.mySolver->pressure_tm1[nd];
    pressure += sh * val;
  }

  return pressure;
}

static int32_t station_output_use_omp(void) {
#ifndef _OPENMP
  return 0;
#else
  if (g_stationOutputOmpMode < 0) {
    const char* env = getenv("HERCULES_STATION_IO_OMP");

    if (env != NULL) {
      if ((strcmp(env, "0") == 0) || (strcasecmp(env, "no") == 0) ||
          (strcasecmp(env, "false") == 0)) {
        g_stationOutputOmpMode = 0;
      } else {
        g_stationOutputOmpMode = 1;
      }
    } else {
      g_stationOutputOmpMode = 1;
    }

    if ((g_stationOutputOmpMode == 1) && (Global.myID == 0) &&
        (g_stationOutputOmpReported == 0)) {
      monitor_print("Station output OpenMP: enabled (threads=%d)\n",
                    omp_get_max_threads());
      g_stationOutputOmpReported = 1;
    }
  }

  return g_stationOutputOmpMode;
#endif
}

static int32_t station_interval_filter_enabled(void) {
  if (g_stationIntervalFilterMode < 0) {
    const char* env = getenv("HERCULES_STATION_ROUTE_MODE");

    if (env != NULL) {
      if ((strcmp(env, "0") == 0) || (strcasecmp(env, "all") == 0) ||
          (strcasecmp(env, "search-all") == 0) ||
          (strcasecmp(env, "off") == 0) || (strcasecmp(env, "no") == 0) ||
          (strcasecmp(env, "false") == 0)) {
        g_stationIntervalFilterMode = 0;
      } else {
        g_stationIntervalFilterMode = 1;
      }
    } else {
      g_stationIntervalFilterMode = 1;
    }

    if ((Global.myID == 0) && (g_stationIntervalFilterReported == 0)) {
      monitor_print("station route mode: %s\n",
                    g_stationIntervalFilterMode ? "interval" : "all");
      g_stationIntervalFilterReported = 1;
    }
  }

  return g_stationIntervalFilterMode;
}

static int32_t station_point_in_local_interval(vector3D_t point) {
  const point_t* intervals;
  point_t stationPoint;

  if (!station_interval_filter_enabled() || Global.theGroupSize <= 1) {
    return 1;
  }

  intervals = octor_getintervaltable(Global.myOctree);
  if (intervals == NULL) {
    return 1;
  }

  stationPoint.x = (tick_t)(point.x[0] / Global.myMesh->ticksize);
  stationPoint.y = (tick_t)(point.x[1] / Global.myMesh->ticksize);
  stationPoint.z = (tick_t)(point.x[2] / Global.myMesh->ticksize);

  if (octor_zcompare(&stationPoint, &intervals[Global.myID]) < 0) {
    return 0;
  }

  if ((Global.myID + 1 < Global.theGroupSize) &&
      (octor_zcompare(&stationPoint, &intervals[Global.myID + 1]) >= 0)) {
    return 0;
  }

  return 1;
}

static int32_t hercules_repro_mode(void) {
  if (g_reproMode < 0) {
    const char* env = getenv("HERCULES_REPRO_MODE");
    if (env != NULL) {
      if ((strcmp(env, "0") == 0) || (strcasecmp(env, "no") == 0) ||
          (strcasecmp(env, "false") == 0) || (strcasecmp(env, "off") == 0)) {
        g_reproMode = 0;
      } else if ((strcmp(env, "1") == 0) || (strcasecmp(env, "on") == 0) ||
                 (strcasecmp(env, "true") == 0)) {
        g_reproMode = 1;
      } else if ((strcmp(env, "2") == 0) ||
                 (strcasecmp(env, "strict") == 0)) {
        g_reproMode = 2;
      } else {
        g_reproMode = 1;
      }
    } else {
      g_reproMode = 0;
    }

    if ((g_reproMode > 0) && (Global.myID == 0)) {
      if (g_reproMode >= 2) {
        monitor_print(
            "Repro mode strict: deterministic MPI unpack + serial CPU OpenMP paths\n");
      } else {
        monitor_print(
            "Repro mode: deterministic MPI unpack path active\n");
      }
    }
  }

  return (g_reproMode > 0) ? 1 : 0;
}

static int32_t hercules_repro_strict(void) {
  if (g_reproMode < 0) {
    (void)hercules_repro_mode();
  }

  return (g_reproMode >= 2) ? 1 : 0;
}

/**
 * Init stations info and data structures
 */
void output_stations_init(const char* numericalin) {
  double t0, t1, setupSeconds, setupSecondsMax;
  int32_t localStations = 0;
  int32_t totalStations = 0;

  t0 = MPI_Wtime();

  Timer_Start("Stations Read Input");
  if (Global.myID == 0) {
    read_stations_info(numericalin);
  }
  Timer_Stop("Stations Read Input");

  Timer_Start("Stations Broadcast");
  broadcast_stations_info();
  Timer_Stop("Stations Broadcast");

  setup_stations_data();

  t1 = MPI_Wtime();
  setupSeconds = t1 - t0;
  localStations = Param.myNumberOfStations;
  MPI_Reduce(&setupSeconds, &setupSecondsMax, 1, MPI_DOUBLE, MPI_MAX, 0,
             comm_solver);
  MPI_Reduce(&localStations, &totalStations, 1, MPI_INT, MPI_SUM, 0,
             comm_solver);

  if (Global.myID == 0) {
    monitor_print(
        "Station setup complete: total=%d max-rank-time=%.2fs\n",
        totalStations, setupSecondsMax);
    if (totalStations != Param.theNumberOfStations) {
      monitor_print(
          "WARNING: located %d of %d configured stations. "
          "Set HERCULES_STATION_ROUTE_MODE=all to compare the full-search path.\n",
          totalStations, Param.theNumberOfStations);
    }
  }

  MPI_Barrier(comm_solver);

  return;
}

/**
 * Reads the profile layers from parameters.in - only done by PE0
 */
static void read_profile(const char* parametersin) {
  static const char* fname = __FUNCTION_NAME;

  int iLayer;
  double* auxtable;
  FILE* fp;

  if ((fp = fopen(parametersin, "r")) == NULL) {
    solver_abort(fname, parametersin, "Error opening parameters.in file\n");
  }

  if (parsetext(fp, "number_profile_layers", 'i', &Param.theNumberOfLayers) !=
      0) {
    solver_abort(fname, NULL, "Error reading the number of layers from %s\n",
                 parametersin);
  }

  auxtable = (double*)malloc(sizeof(double) * Param.theNumberOfLayers * 6);
  Param.theProfileZ = (double*)malloc(sizeof(double) * Param.theNumberOfLayers);
  Param.theProfileVp =
      (double*)malloc(sizeof(double) * Param.theNumberOfLayers);
  Param.theProfileVs =
      (double*)malloc(sizeof(double) * Param.theNumberOfLayers);
  Param.theProfileRho =
      (double*)malloc(sizeof(double) * Param.theNumberOfLayers);
  Param.theProfileQp =
      (double*)malloc(sizeof(double) * Param.theNumberOfLayers);
  Param.theProfileQs =
      (double*)malloc(sizeof(double) * Param.theNumberOfLayers);

  if ((auxtable == NULL) || (Param.theProfileZ == NULL) ||
      (Param.theProfileVp == NULL) || (Param.theProfileVs == NULL) ||
      (Param.theProfileRho == NULL) || (Param.theProfileQp == NULL) ||
      (Param.theProfileQs == NULL)) {
    solver_abort(fname, NULL, "Error allocating memory for the profile\n");
  }

  if (parsedarray(fp, "profile_layers", Param.theNumberOfLayers * 6,
                  auxtable) != 0) {
    solver_abort(fname, NULL, "Error parsing the profile layers from %s\n",
                 parametersin);
  }

  for (iLayer = 0; iLayer < Param.theNumberOfLayers; iLayer++) {
    Param.theProfileZ[iLayer] = auxtable[iLayer * 6];
    Param.theProfileVp[iLayer] = auxtable[iLayer * 6 + 1];
    Param.theProfileVs[iLayer] = auxtable[iLayer * 6 + 2];
    Param.theProfileRho[iLayer] = auxtable[iLayer * 6 + 3];
    Param.theProfileQp[iLayer] = auxtable[iLayer * 6 + 4];
    Param.theProfileQs[iLayer] = auxtable[iLayer * 6 + 5];
  }

  return;
}

/**
 * Broadcasts the profile to all other PEs
 */
void broadcast_profile() {
  static const char* fname = __FUNCTION_NAME;

  MPI_Bcast(&Param.theNumberOfLayers, 1, MPI_INT, 0, comm_solver);

  if (Global.myID != 0) {
    Param.theProfileZ =
        (double*)malloc(sizeof(double) * Param.theNumberOfLayers);
    Param.theProfileVp =
        (double*)malloc(sizeof(double) * Param.theNumberOfLayers);
    Param.theProfileVs =
        (double*)malloc(sizeof(double) * Param.theNumberOfLayers);
    Param.theProfileRho =
        (double*)malloc(sizeof(double) * Param.theNumberOfLayers);
    Param.theProfileQp =
        (double*)malloc(sizeof(double) * Param.theNumberOfLayers);
    Param.theProfileQs =
        (double*)malloc(sizeof(double) * Param.theNumberOfLayers);
  }

  if (Param.theProfileZ == NULL || Param.theProfileVp == NULL ||
      Param.theProfileVs == NULL || Param.theProfileRho == NULL ||
      Param.theProfileQp == NULL || Param.theProfileQs == NULL) {
    solver_abort(fname, NULL, "Error allocating memory for the profile\n");
  }

  MPI_Barrier(comm_solver);

  MPI_Bcast(Param.theProfileZ, Param.theNumberOfLayers, MPI_DOUBLE, 0,
            comm_solver);
  MPI_Bcast(Param.theProfileVp, Param.theNumberOfLayers, MPI_DOUBLE, 0,
            comm_solver);
  MPI_Bcast(Param.theProfileVs, Param.theNumberOfLayers, MPI_DOUBLE, 0,
            comm_solver);
  MPI_Bcast(Param.theProfileRho, Param.theNumberOfLayers, MPI_DOUBLE, 0,
            comm_solver);
  MPI_Bcast(Param.theProfileQp, Param.theNumberOfLayers, MPI_DOUBLE, 0,
            comm_solver);
  MPI_Bcast(Param.theProfileQs, Param.theNumberOfLayers, MPI_DOUBLE, 0,
            comm_solver);

  return;
}

/**
 * Load the velocity profile is given by the user
 */
void load_profile(const char* parametersin) {
  if (Global.myID == 0) {
    read_profile(parametersin);
  }

  broadcast_profile();
  MPI_Barrier(comm_solver);

  return;
}

/**
 * \note This function should only be called by PE with rank 0.
 */
static int load_output_parameters(const char* numericalin,
                                  output_parameters_t* params) {
  FILE* fp;
  int ret, value;
  char filename[LINESIZE];

  assert(NULL != numericalin);
  assert(NULL != params);
  assert(Global.myID == 0);

  /* Read output parameters from numerical.in */
  fp = fopen(numericalin, "r");

  if (NULL == fp) {
    solver_abort("load_output_parameters", "fopen", "numerical.in=\"%s\"",
                 numericalin);
  }

  params->do_output = 0;
  params->parallel_output = 0;
  params->output_displacement = 0;
  params->output_velocity = 0;
  params->output_debug = 0;

  params->displacement_filename = NULL;
  params->velocity_filename = NULL;

  /* read parameters from configuration file */

  value = 0;
  ret = parsetext(fp, "output_parallel", 'i', &value);

  if (0 == ret && 0 != value) {
    value = 0;
    ret = parsetext(fp, "output_displacement", 'i', &value);

    if (0 == ret && 0 != value) { /* output_displacement = 1 in config */
      ret = read_config_string(fp, "output_displacement_file", filename,
                               LINESIZE);

      if (1 == ret && filename[0] != '\0') {
        params->displacement_filename = strdup(filename);
        params->output_displacement = 1;
      } else {
        solver_abort("load_output_parameters", NULL,
                     "Output displacement file name not specified in "
                     "numerical.in=\"%s\"",
                     numericalin);
      }
    }

    value = 0;
    ret = parsetext(fp, "output_velocity", 'i', &value);

    if (0 == ret && 0 != value) { /* output_velocity = 1 in config */
      ret = read_config_string(fp, "output_velocity_file", filename, LINESIZE);

      if (1 == ret && filename[0] != '\0') {
        params->velocity_filename = strdup(filename);
        params->output_velocity = 1;
      } else {
        solver_abort("load_output_parameters", NULL,
                     "Output velocity file name not specified in "
                     "numerical.in=\"%s\"",
                     numericalin);
      }
    }

    params->stats_filename = "output-stats.txt"; /* default value */

    ret = read_config_string(fp, "output_stats_file", filename, LINESIZE);

    if (1 == ret && filename[0] != '\0') {
      params->stats_filename = strdup(filename);
    }

    params->debug_filename = "output-debug.txt"; /* default value */
    ret = read_config_string(fp, "output_debug_file", filename, LINESIZE);

    if (1 == ret && filename[0] != '\0') {
      params->debug_filename = strdup(filename);
    }

    if (params->output_velocity || params->output_displacement) {
      params->do_output = 1;
    }

    params->parallel_output = 1;

    ret = parsetext(fp, "output_debug", 'i', &value);
    if (0 == ret && 0 != value) { /* output_debug = 1 in config */
      params->output_debug = 1;
    }
  }

  ret = 0;

  fclose(fp);

  return ret;
}

/**
 * Initialize output structures, including the opening of 4D output files.
 *
 * \param numericsin Name of the file with the solver and output parameters
 *  i.e., "numerical.in".
 * \param[out] params output parameters, including filenames.
 *
 * \pre The following global variables should be initialized:
 *  - Global.myID.
 *  - Global.theGroupSize.
 *
 * \post On a successful return, the output argument \c params will be
 *  properly initialized.  If the routine fails, the state of the
 *  \c param struct is undefined.
 *
 * \return 0 on success, -1 on error.
 */
static int output_init_parameters(const char* numericalin,
                                  output_parameters_t* params) {
  /* jc: this ugly #define here is because the catamount compiler does not
   * like static const  */
#define VALUES_COUNT 4
  int ret;
  int32_t values[VALUES_COUNT];
  off_t output_steps;

  /* sanity cleanup */
  memset(params, 0, sizeof(output_parameters_t));

  params->do_output = 0;
  params->displacement_filename = NULL;
  params->velocity_filename = NULL;
  params->stats_filename = NULL;

  if (Global.myID == 0) {
    ret = load_output_parameters(numericalin, params);

    if (0 != ret) {
      solver_abort("output_init_parameters", NULL, NULL);
      return -1;
    }
  }

  /* parameters that can be initialized from global variables */
  params->pe_id = Global.myID;
  params->pe_count = Global.theGroupSize;
  params->total_nodes = Global.theNTotal;
  params->total_elements = Global.theETotal;
  params->output_rate = Param.theRate;
  params->domain_x = Param.theDomainX;
  params->domain_y = Param.theDomainY;
  params->domain_z = Param.theDomainZ;
  params->mesh = Global.myMesh;
  params->solver = (solver_t*)Global.mySolver;
  params->delta_t = Param.theDeltaT;
  params->total_time_steps = Param.theTotalSteps;

  output_steps = (Param.theTotalSteps - 1) / Param.theRate + 1;

  values[0] = params->parallel_output;
  values[1] = params->output_displacement;
  values[2] = params->output_velocity;
  values[3] = params->output_debug;

  MPI_Bcast(values, VALUES_COUNT, MPI_INT, 0, comm_solver);
  MPI_Bcast(&params->total_nodes, 1, MPI_INT64_T, 0, comm_solver);

  params->parallel_output = values[0];
  params->output_displacement = values[1];
  params->output_velocity = values[2];
  params->output_debug = values[3];
  params->output_size =
      (output_steps * params->total_nodes * sizeof(fvector_t) +
       sizeof(out_hdr_t));
  Global.theNTotal = params->total_nodes;

  if (params->parallel_output) {
    if (params->output_displacement) {
      broadcast_string(&params->displacement_filename, 0, comm_solver);
    }

    if (params->output_velocity) {
      broadcast_string(&params->velocity_filename, 0, comm_solver);
    }

    if (params->output_debug) {
      broadcast_string(&params->debug_filename, 0, comm_solver);
    }

    /* set the expected file size */
    Param.the4DOutSize = params->output_size;
  }

  return 0;
#undef VALUES_COUNT
}

/**
 * Intialize parallel output structures according to config file
 * parameters.
 *
 * \note params parameter could be a local var instead, it doesn't need
 * to be global, since this is not really used after the initialization
 * except for the stats file name.
 *
 * \return 0 on success, -1 on error (probably aborts instead).
 */
static int output_init(const char* numericalin, output_parameters_t* params) {
  int ret = -1;

  assert(NULL != numericalin);
  assert(NULL != params);

  ret = output_init_parameters(numericalin, params);

  if (ret != 0) {
    return -1;
  }

  /* initialize structures, open files, etc */
  ret = output_init_state(params);

  /* these aren't needed anymore after this point */
  xfree_char(&params->displacement_filename);
  xfree_char(&params->velocity_filename);

  return ret;
}

static int output_get_stats(void) {
  output_stats_t disp_stats, vel_stats;
  int ret = 0;
  // double avg_tput = 0;

  if (Param.theOutputParameters.parallel_output) {
    ret = output_collect_io_stats(
        Param.theOutputParameters.stats_filename, &disp_stats, &vel_stats,
        Timer_Value("Total Wall Clock", (TimerKind)0));

    /* magic trick so print_timing_stat() prints something sensible
     * for the 4D output time in the parallel case.
     *
     * if both displacement and velocity where written out, prefer
     * the displacement stat.
     */
    // if (Param.theOutputParameters.output_velocity)
    // {
    //     avg_tput = vel_stats.tput_avg;
    // }

    // if (Param.theOutputParameters.output_displacement)
    // {
    //     avg_tput = disp_stats.tput_avg;
    // }
  }

  return ret;
}

/**
 * Adjust the values of the properties (Vp, Vs, Rho) for the mesh elements.
 * Initial implementation (by Leo) querying the middle point.1D.
 * Modified implementation (Ricardo) querying the 27 points and averaging.
 *
 * \param cvm the material model database (CVM etree).
 */
static void mesh_correct_properties(etree_t* cvm) {
  elem_t* elemp;
  edata_t* edata;
  int32_t eindex;
  double east_m, north_m, depth_m, VpVsRatio, RhoVpRatio;
  int res, iNorth, iEast, iDepth, numPoints = 3, cnt = 0;
  int acoustic_cnt = 0, elastic_cnt = 0;
  double vs, vp, rho, s_0, qp, qs;
  double acoustic_vp, acoustic_rho, acoustic_qp;
  double elastic_vp, elastic_vs, elastic_rho, elastic_qp, elastic_qs;
  double points[3];
  int32_t lnid0;
  int64_t acoustic_elements = 0;
  int64_t mixed_elements = 0;
  vector3D_t IstVelModel_origin;
  vector3D_t basinVelModel_origin;
  vector3D_t threeDVelModel_origin;

  double Qs, Qp, Qk, L;

  points[0] = 0.005;
  points[1] = 0.5;
  points[2] = 0.995;

  if (Param.IstanbulVelocityModel == YES) {
    IstVelModel_origin = compute_domain_coords_linearinterp(
        28.675, 40.9565, Param.theSurfaceCornersLong,
        Param.theSurfaceCornersLat, Param.theDomainY, Param.theDomainX,
        &Param.utmZone);
  }
  if (Param.basinVelocityModel == YES) {
    basinVelModel_origin = compute_domain_coords_linearinterp(
        Param.theBasinLong, Param.theBasinLat, Param.theSurfaceCornersLong,
        Param.theSurfaceCornersLat, Param.theDomainY, Param.theDomainX,
        &Param.utmZone);
  }
  if (Param.threeDVelocityModel == YES) {
    threeDVelModel_origin = compute_domain_coords_linearinterp(
        Param.the3DVelocityModelLong, Param.the3DVelocityModelLat,
        Param.theSurfaceCornersLong, Param.theSurfaceCornersLat,
        Param.theDomainY, Param.theDomainX, &Param.utmZone);
  }

  if (Global.myID == 0) {
    monitor_print("mesh_correct_properties  ... ");
    fflush(stdout);
  }

  /* iterate over mesh elements */
  for (eindex = 0; eindex < Global.myMesh->lenum; eindex++) {
    elemp = &Global.myMesh->elemTable[eindex];
    edata = (edata_t*)elemp->data;
    lnid0 = Global.myMesh->elemTable[eindex].lnid[0];

    if (Param.includeTopography == YES) {
      if (topo_correctproperties(edata)) {
        continue;
      }
    }

    if (Param.includeBuildings == YES) {
      if (bldgs_correctproperties(Global.myMesh, edata, lnid0)) {
        continue;
      }
    }

    vp = 0;
    vs = 0;
    rho = 0;
    cnt = 0;
    s_0 = 0;
    qp = 0;
    qs = 0;
    acoustic_vp = 0;
    acoustic_rho = 0;
    acoustic_qp = 0;
    acoustic_cnt = 0;
    elastic_vp = 0;
    elastic_vs = 0;
    elastic_rho = 0;
    elastic_qp = 0;
    elastic_qs = 0;
    elastic_cnt = 0;

    for (iNorth = 0; iNorth < numPoints; iNorth++) {
      north_m = (Global.myMesh->ticksize) *
                    (double)Global.myMesh->nodeTable[lnid0].x +
                edata->edgesize * points[iNorth] + Global.theXForMeshOrigin;

      for (iEast = 0; iEast < numPoints; iEast++) {
        east_m = ((Global.myMesh->ticksize) *
                      (double)Global.myMesh->nodeTable[lnid0].y +
                  edata->edgesize * points[iEast] + Global.theYForMeshOrigin);

        for (iDepth = 0; iDepth < numPoints; iDepth++) {
          cvmpayload_t g_props; /* ground properties record */

          depth_m =
              ((Global.myMesh->ticksize) *
                   (double)Global.myMesh->nodeTable[lnid0].z +
               edata->edgesize * points[iDepth] + Global.theZForMeshOrigin);

          /* NOTE: If you want to see the carving process,
           *       activate this and comment the query below */
          if (Param.includeBuildings == YES) {
            depth_m -= get_surface_shift();
          }

          if ((Param.includeTopography == YES) &&
              (get_theetree_type() == SQD)) {
            depth_m -= point_elevation(north_m, east_m);
            if (depth_m < 0) continue;

            if (depth_m > Param.theDomainZ)
              depth_m = Param.theDomainZ - edata->edgesize * 0.005;
          }

          if (belongs2hmgHalfspace(east_m, north_m, depth_m))
            res = get_halfspaceproperties(&g_props);
          else if (Param.IstanbulVelocityModel == YES) {
            res = getIstanbulMaterial(&g_props, IstVelModel_origin, north_m,
                                      east_m, depth_m);
          } else if (Param.basinVelocityModel == YES) {
            res = getBasinMaterial(&g_props, basinVelModel_origin, north_m,
                                   east_m, depth_m);
          } else if (Param.threeDVelocityModel == YES) {
            res = get3DMaterial(&g_props, threeDVelModel_origin, north_m,
                                east_m, depth_m);
          } else if (Param.useProfile == NO) {
            res =
                cvm_query(Global.theCVMEp, east_m, north_m, depth_m, &g_props);
          } else {
            res = profile_query(depth_m, &g_props);
          }

          if (res != 0) {
            fprintf(stderr,
                    "Cannot find the query point: east = %lf, north = %lf, "
                    "depth = %lf \n",
                    east_m, north_m, depth_m);
            exit(1);
          }

          if (material_values_are_acoustic(g_props.Vp, g_props.Vs,
                                           g_props.rho)) {
            acoustic_vp += g_props.Vp;
            acoustic_rho += g_props.rho;
            acoustic_qp += g_props.Qp;
            ++acoustic_cnt;
          } else if (g_props.Vs <= 0.0) {
            solver_abort(__FUNCTION_NAME, NULL,
                         "Encountered Vs <= 0 with include_acoustic = no. "
                         "Set include_acoustic = yes for water/acoustic "
                         "regions.\n");
          } else {
            elastic_vp += g_props.Vp;
            elastic_vs += g_props.Vs;
            elastic_rho += g_props.rho;
            elastic_qp += g_props.Qp;
            elastic_qs += g_props.Qs;
            ++elastic_cnt;
          }

          vp += g_props.Vp;
          vs += g_props.Vs;
          rho += g_props.rho;
          qp += g_props.Qp;
          qs += g_props.Qs;
          ++cnt;

          // get geostatic stress as 1d column
          double nlayers = 100,
                 depth_o = (depth_m - edata->edgesize / 2.0) / nlayers, depth_k;
          int k;
          if (iNorth == 1 && iEast == 1 && iDepth == 1 &&
              Param.includeNonlinearAnalysis == YES) {
            // NOTE from Clifford: This nonlinear section is only applicable for
            // the Istanbul Velocity Model as the nonlinear_id is calculated
            // based on the Istanbul Velocity Model.

            // get nonlinear_id at the elment's center
            double output[4] = {0.0}, x_rel, y_rel;
            x_rel = 4536400.00 + north_m - IstVelModel_origin.x[0];
            y_rel = 388500.00 + east_m - IstVelModel_origin.x[1];

            res = material_property_relative_V10_local(
                y_rel, x_rel, -depth_m, output, IstVelModel_origin.x[1],
                IstVelModel_origin.x[0]);

            edata->nl_id = output[3];

            s_0 = edata->edgesize / 2.0 * edata->rho * 9.81;
            for (k = 0; k < nlayers; k++) {
              depth_k = depth_o * (k + 0.5);

              if (belongs2hmgHalfspace(east_m, north_m, depth_m))
                res = get_halfspaceproperties(&g_props);
              /* NOTE: This section is a little bit different from the previous
              Istanbul Velocity Model. g_props->Qs and g_props->Qp are not
              assigned in the original section. But it should be fine if they
              are still assigned. */
              else if (Param.IstanbulVelocityModel == YES) {
                res = getIstanbulMaterial(&g_props, IstVelModel_origin, north_m,
                                          east_m, depth_k);
              } else if (Param.basinVelocityModel == YES) {
                res = getBasinMaterial(&g_props, basinVelModel_origin, north_m,
                                       east_m, depth_k);
              } else if (Param.threeDVelocityModel == YES) {
                res = get3DMaterial(&g_props, threeDVelModel_origin, north_m,
                                    east_m, depth_k);
              } else if (Param.useProfile == NO) {
                res = cvm_query(Global.theCVMEp, east_m, north_m, depth_k,
                                &g_props);
              } else {
                res = profile_query(depth_k, &g_props);
              }

              s_0 += depth_o * g_props.rho * 9.81;
            }

            edata->sigma_0 = s_0;
          }
        }
      }
    }

    if ((Param.includeAcoustic == YES) && (acoustic_cnt > 0) &&
        (elastic_cnt == 0)) {
      edata->Vp = acoustic_vp / acoustic_cnt;
      edata->Vs = 0.0f;
      edata->rho = acoustic_rho / acoustic_cnt;
      qp = (acoustic_cnt > 0) ? (acoustic_qp / acoustic_cnt) : 1000.0;
      qs = 1000.0;
      acoustic_elements++;
    } else if ((Param.includeAcoustic == YES) && (elastic_cnt > 0) &&
               (acoustic_cnt > 0)) {
      edata->Vp = elastic_vp / elastic_cnt;
      edata->Vs = elastic_vs / elastic_cnt;
      edata->rho = elastic_rho / elastic_cnt;
      qp = elastic_qp / elastic_cnt;
      qs = elastic_qs / elastic_cnt;
      mixed_elements++;
    } else {
      edata->Vp = vp;
      edata->Vs = vs;
      edata->rho = rho;
      // edata->Qp  = qp;
      // edata->Qs  = qs;

      if (cnt != 0) {
        edata->Vp = vp / cnt;
        edata->Vs = vs / cnt;
        edata->rho = rho / cnt;
        qp = qp / cnt;
        qs = qs / cnt;

        // edata->Qp  = qp / cnt;
        // edata->Qs  = qs / cnt;
      }
    }

    if (element_is_acoustic(edata)) {
      acoustic_sanitize_element(edata);
      continue;
    }

    if (edata->Vs <= 0.0f) {
      solver_abort(__FUNCTION_NAME, NULL,
                   "Encountered Vs <= 0 with include_acoustic = no. "
                   "Set include_acoustic = yes for water/acoustic regions.\n");
    }

    /* Auxiliary ratios for adjustments */
    VpVsRatio = edata->Vp / edata->Vs;
    RhoVpRatio = edata->rho / edata->Vp;

    // Doriam says: Adjust minimum value as in La Habra project
    if (VpVsRatio <= 1.453) VpVsRatio = 1.453;  // makes nu=0.05

    if (VpVsRatio > Param.theThresholdVpVs) VpVsRatio = Param.theThresholdVpVs;

    edata->Vp = VpVsRatio * edata->Vs;

    /* Adjust material properties according to the element size and
     * softening factor.
     *
     * A factor of 1 means perfect compliance between the mesh and the
     * elements' material properties resulting in strong changes to the
     * results. A factor of 4 tends to double the simulation delta_t
     * without affecting too much the results. Testing is needed but for
     * now I recommend factors > 4.
     */
    if (Param.theSofteningFactor > 0) {
      double idealVs, factoredVs;

      idealVs = edata->edgesize * Param.theFactor;
      factoredVs = idealVs * Param.theSofteningFactor;

      if (edata->Vs > factoredVs) {
        edata->Vs = factoredVs;
        edata->Vp = factoredVs * VpVsRatio;
        edata->rho = edata->Vp * RhoVpRatio;
      }
    }

    /* Readjust Vs, Vp and Density according to VsCut */
    if (edata->Vs < Param.theVsCut) {
      edata->Vs = Param.theVsCut;
      edata->Vp = Param.theVsCut * VpVsRatio;
      // Clifford's NOTE: The following line is uncommented after discussing
      // with Wenyang
      edata->rho = edata->Vp * RhoVpRatio;
    }

    /*
     * Definition of BTK-Family Damping parameters
     */
    if (Param.theTypeOfDamping >= BKT) {
      double vs_vp_Ratio, vksquared;
      double vs_kms = edata->Vs * 0.001; /* Vs in km/s */

      vksquared = edata->Vp * edata->Vp - 4. / 3. * edata->Vs * edata->Vs;

      vs_vp_Ratio = edata->Vs / edata->Vp;
      L = 4. / 3. * vs_vp_Ratio *
          vs_vp_Ratio; /* As defined in Shearer (2009) */

      if (Param.useParametricQ == YES) {
        /* Use of the formula
         *
         * Qs = C + QAlpha * ( Vs ^ (QBeta) )
         *
         * This formula was introduce by Ricardo and Naeem
         * and it is versatile enough and simpler than the
         * option used in Taborda and Bielak (2013, BSSA)
         */
        Qs = Param.theQConstant + Param.theQAlpha * pow(vs_kms, Param.theQBeta);
      } else {
        /* Use of the formula introduced by Ricardo in the
         * paper Taborda and Bielak (2013, BSSA) which is
         * based on the idea of Brocher (2005)
         */
        if (Param.useProfile == YES)
          Qs = qs;
        else
          Qs =
              10.5 +
              vs_kms *
                  (-16. +
                   vs_kms *
                       (153. +
                        vs_kms * (-103. +
                                  vs_kms * (34.7 + vs_kms * (-5.29 +
                                                             vs_kms * 0.31)))));
      }

      /* Default option for Qp */
      if (Param.useProfile == YES)
        Qp = qp;
      else
        Qp = 2. * Qs;

      // edata->Qs = Qs;  // update Qs and Qp
      // edata->Qp = Qp;

      if (Param.useInfQk == YES) {
        Qk = 1000;
      } else {
        Qk = (1. - L) / (1. / Qp - L / Qs);
      }

      if (Param.theTypeOfDamping == BKT) {
        /* Legacy implementation of original BKT
         * model using a table to set parameters
         */

        int index_Qs, index_Qk;
        int QTable_Size =
            (int)(sizeof(Global.theQTABLE) / (6 * sizeof(double)));

        index_Qs =
            Search_Quality_Table(Qs, &(Global.theQTABLE[0][0]), QTable_Size);
        index_Qk =
            Search_Quality_Table(Qk, &(Global.theQTABLE[0][0]), QTable_Size);

        if ((index_Qs == -2) || (index_Qs >= QTable_Size) || (index_Qk == -2) ||
            (index_Qk >= QTable_Size)) {
          solver_abort(__FUNCTION_NAME, NULL, "Unexpected damping type: %d\n",
                       Param.theTypeOfDamping);
        }

        if (index_Qs == -1) {
          edata->a0_shear = 0;
          edata->a1_shear = 0;
          edata->g0_shear = 0;
          edata->g1_shear = 0;
          edata->b_shear = 0;
        } else {
          edata->a0_shear = Global.theQTABLE[index_Qs][1];
          edata->a1_shear = Global.theQTABLE[index_Qs][2];
          edata->g0_shear = Global.theQTABLE[index_Qs][3];
          edata->g1_shear = Global.theQTABLE[index_Qs][4];
          edata->b_shear = Global.theQTABLE[index_Qs][5];
        }

        if ((Param.useInfQk == YES) || (index_Qk == -1)) {
          edata->a0_kappa = 0;
          edata->a1_kappa = 0;
          edata->g0_kappa = 0;
          edata->g1_kappa = 0;
          edata->b_kappa = 0;
        } else {
          edata->a0_kappa = Global.theQTABLE[index_Qk][1];
          edata->a1_kappa = Global.theQTABLE[index_Qk][2];
          edata->g0_kappa = Global.theQTABLE[index_Qk][3];
          edata->g1_kappa = Global.theQTABLE[index_Qk][4];
          edata->b_kappa = Global.theQTABLE[index_Qk][5];
        }
      } else if (Param.theTypeOfDamping == BKT2) {
        /* Notes:
         * g = normilized_gamma * ( 2 * pi * f_max);
         * b = normilized_beta / (Q * 2 * pi * f_max);
         */
        if (Qs >= 1000) {
          edata->g0_shear = 0;
          edata->g1_shear = 0;
          edata->a0_shear = 0;
          edata->a1_shear = 0;
          edata->b_shear = 0;
        } else {
          // Doriam: Ricardo's functions
          /* edata->g0_shear =   0.0373 * (2. * M_PI * Param.theFreq);
          edata->g1_shear =   0.3082 * (2. * M_PI * Param.theFreq);
          edata->a0_shear = (-2.656  * pow(Qs, -0.8788) + 1.677 ) / Qs;
          edata->a1_shear = (-0.5623 * pow(Qs, -1.0300) + 1.262 ) / Qs;
          edata->b_shear  = ( 0.1876 * pow(Qs, -0.9196) + 0.6137) / (Qs * (2. *
          M_PI * Param.theFreq)); */

          // Doriam: All parameters free
          /* edata->g0_shear = ( 0.008634 * pow(Qs, -0.07081) + 0.04408) * (2. *
          M_PI * Param.theFreq); edata->g1_shear = ( 0.08366 * pow(Qs, -0.1043)
          + 0.3487) * (2. * M_PI * Param.theFreq); edata->a0_shear = ( 1.042  *
          pow(Qs, -0.9077) ) ; edata->a1_shear = ( 1.056 * pow(Qs, -0.9718) ) ;
          edata->b_shear  = ( 0.5711 * pow(Qs, -1.0) ) / ( (2. * M_PI *
          Param.theFreq)); */

          // Doriam: Gamma parameters fixed
          edata->g0_shear = 0.05 * (2. * M_PI * Param.theFreq);
          edata->g1_shear = 0.40 * (2. * M_PI * Param.theFreq);
          edata->a0_shear = (1.056 * pow(Qs, -0.9068));
          edata->a1_shear = (1.078 * pow(Qs, -0.9797));
          edata->b_shear =
              (0.5709 * pow(Qs, -1.021)) / ((2. * M_PI * Param.theFreq));
        }

        if ((Param.useInfQk == YES) || (Qk >= 1000)) {
          edata->g0_kappa = 0;
          edata->g1_kappa = 0;
          edata->a0_kappa = 0;
          edata->a1_kappa = 0;
          edata->b_kappa = 0;
        } else {
          // Doriam: Ricardo's functions
          /* edata->g0_kappa =   0.0373 * (2. * M_PI * Param.theFreq);
          edata->g1_kappa =   0.3082 * (2. * M_PI * Param.theFreq);
          edata->a0_kappa = (-2.656  * pow(Qk, -0.8788) + 1.677 ) / Qk;
          edata->a1_kappa = (-0.5623 * pow(Qk, -1.0300) + 1.262 ) / Qk;
          edata->b_kappa  = ( 0.1876 * pow(Qk, -0.9196) + 0.6137) / (Qk * (2. *
          M_PI * Param.theFreq)); */

          // Doriam: All parameters free
          /* edata->g0_kappa = ( 0.008634 * pow(Qk, -0.07081) + 0.04408 ) * (2.
          * M_PI * Param.theFreq); edata->g1_kappa = ( 0.08366 * pow(Qk,
          -0.1043) + 0.3487 ) * (2. * M_PI * Param.theFreq); edata->a0_kappa =
          ( 1.042  * pow(Qk, -0.9077) ); edata->a1_kappa = ( 1.056 * pow(Qk,
          -0.9718) ); edata->b_kappa  = ( 0.5711 * pow(Qk, -1.0) ) / ( (2. *
          M_PI * Param.theFreq)); */

          // Doriam: Gamma parameters fixed
          edata->g0_kappa = 0.05 * (2. * M_PI * Param.theFreq);
          edata->g1_kappa = 0.40 * (2. * M_PI * Param.theFreq);
          edata->a0_kappa = (1.056 * pow(Qk, -0.9068));
          edata->a1_kappa = (1.078 * pow(Qk, -0.9797));
          edata->b_kappa =
              (0.5709 * pow(Qk, -1.021)) / ((2. * M_PI * Param.theFreq));
        }
      } else if (Param.theTypeOfDamping == BKT3) {
        if (Qs >= 1000) {
          edata->g0_shear = 0;
          edata->g1_shear = 0;
          edata->g2_shear = 0;
          edata->a0_shear = 0;
          edata->a1_shear = 0;
          edata->a2_shear = 0;
          edata->b_shear = 0;
        } else {
          // Doriam: Ricardo's functions
          /* edata->g0_shear =   0.0151 * (2. * M_PI * Param.theFreq);
          edata->g1_shear =   0.1000 * (2. * M_PI * Param.theFreq);
          edata->g2_shear =   0.4814 * (2. * M_PI * Param.theFreq);
          edata->a0_shear = (-2.723  * pow(Qs, -0.8206) + 1.601 ) / Qs;
          edata->a1_shear = (-1.439  * pow(Qs, -0.9668) + 1.04  ) / Qs;
          edata->a2_shear = (-0.3037 * pow(Qs, -0.8911) + 1.032 ) / Qs;
          edata->b_shear  = ( 0.1249 * pow(Qs, -0.804 ) + 0.4782) / (Qs * (2. *
          M_PI * Param.theFreq)); */

          // Doriam: Updated fitting. Gamma Fixed
          edata->g0_shear = 0.0270 * (2. * M_PI * Param.theFreq);
          edata->g1_shear = 0.5900 * (2. * M_PI * Param.theFreq);
          edata->g2_shear = 0.1400 * (2. * M_PI * Param.theFreq);
          edata->a0_shear = (0.8543 * pow(Qs, -0.8953));
          edata->a1_shear = (0.8871 * pow(Qs, -0.9803));
          edata->a2_shear = (0.6718 * pow(Qs, -0.9412));
          edata->b_shear =
              (0.446 * pow(Qs, -1.017)) / ((2. * M_PI * Param.theFreq));
        }

        if ((Param.useInfQk == YES) || (Qk >= 1000)) {
          edata->g0_kappa = 0;
          edata->g1_kappa = 0;
          edata->g2_kappa = 0;
          edata->a0_kappa = 0;
          edata->a1_kappa = 0;
          edata->a2_kappa = 0;
          edata->b_kappa = 0;
        } else {
          // Doriam: Ricardo's functions
          /* edata->g0_kappa =   0.0151 * (2. * M_PI * Param.theFreq);
          edata->g1_kappa =   0.1000 * (2. * M_PI * Param.theFreq);
          edata->g2_kappa =   0.4814 * (2. * M_PI * Param.theFreq);
          edata->a0_kappa = (-2.723  * pow(Qk, -0.8206) + 1.601 ) / Qk;
          edata->a1_kappa = (-1.439  * pow(Qk, -0.9668) + 1.04  ) / Qk;
          edata->a2_kappa = (-0.3037 * pow(Qk, -0.8911) + 1.032 ) / Qk;
          edata->b_kappa  = ( 0.1249 * pow(Qk, -0.804 ) + 0.4782) / (Qk * (2. *
          M_PI * Param.theFreq)); */

          // Doriam: Updated fitting. Gamma Fixed
          edata->g0_kappa = 0.0270 * (2. * M_PI * Param.theFreq);
          edata->g1_kappa = 0.5900 * (2. * M_PI * Param.theFreq);
          edata->g2_kappa = 0.1400 * (2. * M_PI * Param.theFreq);

          edata->a0_kappa = (0.8543 * pow(Qk, -0.8953));
          edata->a1_kappa = (0.8871 * pow(Qk, -0.9803));
          edata->a2_kappa = (0.6718 * pow(Qk, -0.9412));
          edata->b_kappa =
              (0.446 * pow(Qk, -1.017)) / ((2. * M_PI * Param.theFreq));
        }
      } else if (Param.theTypeOfDamping == BKT3F) {
        /*
         * Temporal Implementation of frequency dependent Q with:
         * exponent = 0.8
         * f_max = 10 Hz (fixed)
         * f_o = 0.1 f_max = 1 Hz (fixed)
         */

        if (Qs >= 1000) {
          edata->g0_shear = 0;
          edata->g1_shear = 0;
          edata->g2_shear = 0;
          edata->a0_shear = 0;
          edata->a1_shear = 0;
          edata->a2_shear = 0;
          edata->b_shear = 0;
        } else {
          edata->g0_shear = 0.002 * (2. * M_PI * 10);
          edata->g1_shear = 0.0116 * (2. * M_PI * 10);
          edata->g2_shear = 0.0798 * (2. * M_PI * 10);
          edata->a0_shear = (-2.809 * pow(Qs, -0.7919) + 1.512) / Qs;
          edata->a1_shear = (-1.748 * pow(Qs, -0.882) + 1.064) / Qs;
          edata->a2_shear = (-2.358 * pow(Qs, -1.725) + 1.581) / Qs;
          edata->b_shear =
              (0.09232 * pow(Qs, -0.8876) + 0.006941) / (Qs * (2. * M_PI * 10));
        }

        if ((Param.useInfQk == YES) || (Qk >= 1000)) {
          edata->g0_kappa = 0;
          edata->g1_kappa = 0;
          edata->g2_kappa = 0;
          edata->a0_kappa = 0;
          edata->a1_kappa = 0;
          edata->a2_kappa = 0;
          edata->b_kappa = 0;
        } else {
          edata->g0_kappa = 0.002 * (2. * M_PI * 10);
          edata->g1_kappa = 0.0116 * (2. * M_PI * 10);
          edata->g2_kappa = 0.0798 * (2. * M_PI * 10);
          edata->a0_kappa = (-2.809 * pow(Qk, -0.7919) + 1.512) / Qk;
          edata->a1_kappa = (-1.748 * pow(Qk, -0.882) + 1.064) / Qk;
          edata->a2_kappa = (-2.358 * pow(Qk, -1.725) + 1.581) / Qk;
          edata->b_kappa =
              (0.09232 * pow(Qk, -0.8876) + 0.006941) / (Qk * (2. * M_PI * 10));
        }
      } else {
        /* Should never reach this point */
        solver_abort(__FUNCTION_NAME, NULL, "Unexpected damping type: %d\n",
                     Param.theTypeOfDamping);
      }

      /*
       * Phase velocity adjustment
       */
      if (Param.theFreq_Vel != 0.) {
        double w, w2;
        double shear_vel_corr_factor = 1.0;
        double kappa_vel_corr_factor = 1.0;

        w = Param.theFreq_Vel / Param.theFreq;
        w2 = w * w;

        if ((Param.theTypeOfDamping == BKT) ||
            (Param.theTypeOfDamping == BKT2)) {
          if ((edata->a0_shear != 0) && (edata->a1_shear != 0)) {
            shear_vel_corr_factor =
                sqrt(1. - (edata->a0_shear * edata->g0_shear * edata->g0_shear /
                               (edata->g0_shear * edata->g0_shear + w2) +
                           edata->a1_shear * edata->g1_shear * edata->g1_shear /
                               (edata->g1_shear * edata->g1_shear + w2)));
          }
          if ((edata->a0_kappa != 0) && (edata->a1_kappa != 0)) {
            kappa_vel_corr_factor =
                sqrt(1. - (edata->a0_kappa * edata->g0_kappa * edata->g0_kappa /
                               (edata->g0_kappa * edata->g0_kappa + w2) +
                           edata->a1_kappa * edata->g1_kappa * edata->g1_kappa /
                               (edata->g1_kappa * edata->g1_kappa + w2)));
          }
        } else if ((Param.theTypeOfDamping == BKT3) ||
                   (Param.theTypeOfDamping == BKT3F)) {
          if ((edata->a0_shear != 0) && (edata->a1_shear != 0) &&
              (edata->a2_shear != 0)) {
            shear_vel_corr_factor =
                sqrt(1. - (edata->a0_shear * edata->g0_shear * edata->g0_shear /
                               (edata->g0_shear * edata->g0_shear + w2) +
                           edata->a1_shear * edata->g1_shear * edata->g1_shear /
                               (edata->g1_shear * edata->g1_shear + w2) +
                           edata->a2_shear * edata->g2_shear * edata->g2_shear /
                               (edata->g2_shear * edata->g2_shear + w2)));
          }
          if ((edata->a0_kappa != 0) && (edata->a1_kappa != 0) &&
              (edata->a2_kappa != 0)) {
            kappa_vel_corr_factor =
                sqrt(1. - (edata->a0_kappa * edata->g0_kappa * edata->g0_kappa /
                               (edata->g0_kappa * edata->g0_kappa + w2) +
                           edata->a1_kappa * edata->g1_kappa * edata->g1_kappa /
                               (edata->g1_kappa * edata->g1_kappa + w2) +
                           edata->a2_kappa * edata->g2_kappa * edata->g2_kappa /
                               (edata->g2_kappa * edata->g2_kappa + w2)));
          }
        } else {
          /* Should never reach this point */
          solver_abort(__FUNCTION_NAME, NULL, "Unexpected damping type: %d\n",
                       Param.theTypeOfDamping);
        }

        edata->Vs = shear_vel_corr_factor * edata->Vs;
        edata->Vp =
            sqrt(kappa_vel_corr_factor * kappa_vel_corr_factor * vksquared +
                 4. / 3. * edata->Vs * edata->Vs);
      }
    }
  }

  if (Param.includeAcoustic == YES) {
    int64_t local_counts[2];
    int64_t global_counts[2];

    local_counts[0] = acoustic_elements;
    local_counts[1] = mixed_elements;
    MPI_Reduce(local_counts, global_counts, 2, MPI_INT64_T, MPI_SUM, 0,
               comm_solver);
    if (Global.myID == 0) {
      monitor_print(
          "Acoustic mesh material stats: acoustic elements %.1fM, "
          "mixed interface elements kept elastic %.1fM\n",
          (double)global_counts[0] / 1.0e6,
          (double)global_counts[1] / 1.0e6);
      fflush(stdout);
    }
  }

  mesh_build_acoustic_interface_flags(
      (Param.IstanbulVelocityModel == YES) ? &IstVelModel_origin : NULL,
      (Param.basinVelocityModel == YES) ? &basinVelModel_origin : NULL,
      (Param.threeDVelocityModel == YES) ? &threeDVelModel_origin : NULL);
}

/*** Program's standard entry point. ***/
int main(int argc, char** argv) {
#ifdef DEBUG
  int32_t flag;
  MPI_Status status;
#endif /* DEBUG */
  int need_material_queries = 0;

  /* MPI initialization */
  MPI_Init(&argc, &argv);
  Timer_Start("Total Wall Clock");
  Timer_Start("Phase Sync/Coordination");
  MPI_Barrier(MPI_COMM_WORLD);
  Timer_Stop("Phase Sync/Coordination");
  MPI_Comm_rank(MPI_COMM_WORLD, &Global.myID);
  MPI_Comm_size(MPI_COMM_WORLD, &Global.theGroupSize);

  /* Make sure using correct input arguments */
  if (argc != 2) {
    if (Global.myID == 0) {
      fputs("Usage: psolve <parameters.in>\n", stderr);
    }
    MPI_Finalize();
    exit(1);
  }

  /*Read in and verify IO Pool Configuration */
  char* IO_PES_ENV_VAR;
  IO_PES_ENV_VAR = getenv("IO_PES");
  if (IO_PES_ENV_VAR == NULL)
    Param.IO_pool_pe_count = 0;
  else
    Param.IO_pool_pe_count = atoi(IO_PES_ENV_VAR);
  if (Param.IO_pool_pe_count >= Global.theGroupSize) {
    Param.IO_pool_pe_count = 0;
    if (Global.myID == 0) printf("Warning: IO_PES too large.  Set to 0.\n");
  }

  /* Split off PEs into IO Pool */
  MPI_Comm_dup(MPI_COMM_WORLD, &comm_IO);
  int in_io_pool, color;
  Global.theGroupSize -= Param.IO_pool_pe_count;
  if (Global.myID >= Global.theGroupSize)
    in_io_pool = 1;
  else
    in_io_pool = 0;
  if (in_io_pool)
    color = MPI_UNDEFINED;
  else
    color = 0;
  MPI_Comm_split(MPI_COMM_WORLD, color, Global.myID, &comm_solver);
  if (in_io_pool) {
    planes_IO_PES_main(Global.myID);
    goto IO_PES_REJOIN;
  }

  if (Global.myID == 0) {
    printf("Starting psolve $Revision: 1.166 $ on %d PEs (%d IO Pool PEs).\n\n",
           Global.theGroupSize, Param.IO_pool_pe_count);
    fflush(stdout);
  }

  /* Read input parameters from file */
  Timer_Start("Global Init");
  read_parameters(argc, argv);

  if (Param.IstanbulVelocityModel == YES) {
    g_istanbulVelModelOrigin = compute_domain_coords_linearinterp(
        28.675, 40.9565, Param.theSurfaceCornersLong,
        Param.theSurfaceCornersLat, Param.theDomainY, Param.theDomainX,
        &Param.utmZone);
    g_hasIstanbulVelModelOrigin = 1;
  } else {
    g_hasIstanbulVelModelOrigin = 0;
  }

  need_material_queries =
      (Param.theMeshInputMode == 0) ||
      ((Param.includeAcoustic == YES) &&
       (Param.acousticInterfaceCoupling == YES));

  if (Param.useProfile == YES) {
    /* Read profile to memory */
    load_profile(Param.parameters_input_file);
  } else if (need_material_queries) {
    /* Create and open database */
    open_cvmdb();
  }

  /* Initialize velocity models when meshing or interface sampling uses them. */
  if (Param.IstanbulVelocityModel == YES && need_material_queries) {
    Istanbul_init(Global.myID);
  } else if (Param.threeDVelocityModel == YES && need_material_queries) {
    /* Initialize 3D velocity model */
    general3DVelocityModel_init(Global.myID, Param.the3DVelModelDir);
  } else if (((Param.IstanbulVelocityModel == YES) ||
              (Param.threeDVelocityModel == YES)) &&
             (Global.myID == 0)) {
    monitor_print("Skipping Istanbul/3D velocity model initialization.\n");
  }

  /* Initialize nonlinear parameters */
  if (Param.includeNonlinearAnalysis == YES) {
    nonlinear_init(Global.myID, Param.parameters_input_file, Param.theDeltaT,
                   Param.theEndT);
  }

  if (Param.includeBuildings == YES) {
    bldgs_init(Global.myID, Param.parameters_input_file);
  }

  if (Param.drmImplement == YES) {
    Timer_Start("Init Drm Parameters");
    Param.theDrmPart = drm_init(Global.myID, Param.parameters_input_file,
                                Param.includeBuildings);
    Timer_Stop("Init Drm Parameters");
    Timer_Reduce("Init Drm Parameters", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }

  if (Param.includeTopography == YES) {
    topo_init(Global.myID, Param.parameters_input_file, Param.theTopoDir);
  }

  if (Param.includeIncidentPlaneWaves == YES) {
    /* compute half-space Vs */
    /*      cvmpayload_t props;
            int res = cvm_query( Global.theCVMEp, Param.theDomainY / 2.0,
       Param.theDomainX / 2.0, Param.theDomainZ, &props ); double VsHS =
       props.Vs;*/

    drm_planewaves_init(Global.myID, Param.parameters_input_file);
  }

  if (Param.includeHomogeneousHalfSpace == YES) {
    /* compute half-space Vs */
    hmgHalfspace_init(Global.myID, Param.parameters_input_file);
  }

  // INTRODUCE BKT MODEL
  /* Init Quality Factor Table */
  constract_Quality_Factor_Table();
  Timer_Stop("Global Init");
  Timer_Reduce("Global Init", (TimerKind)(MAX | MIN | AVERAGE), comm_solver);

  /* Generate, partition and output unstructured octree mesh */
  if (Param.theMeshInputMode != 0) {
    mesh_generate_from_etree();
  } else {
    mesh_generate();
  }

  if (Param.includeBuildings == YES) {
    if (get_fixedbase_flag() == YES) {
      bldgs_fixedbase_init(Global.myMesh, Param.theEndT - Param.theStartT);
    }
    bldgs_finalize();
  }

  if (Param.drmImplement == YES) {
    Timer_Start("Drm Init");
    if (Param.theDrmPart == PART0) {
      find_drm_nodes(Global.myMesh, Global.myID, Param.parameters_input_file,
                     Global.myOctree->ticksize, Global.theGroupSize);
    }
    if (Param.theDrmPart == PART1) {
      setup_drm_data(Global.myMesh, Global.myID, Global.theGroupSize);
    }
    if (Param.theDrmPart == PART2) {
      proc_drm_elems(Global.myMesh, Global.myID, Global.theGroupSize,
                     Param.theTotalSteps);
    }
    drm_stats(Global.myID, Global.theGroupSize, Global.theXForMeshOrigin,
              Global.theYForMeshOrigin, Global.theZForMeshOrigin);
    Timer_Stop("Drm Init");
    Timer_Reduce("Drm Init", (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }

  /* Initialize topography solver analysis structures */
  /* must be before solver_init() for proper treatment of the nodal mass */
  if (Param.includeTopography == YES) {
    Timer_Start("Topography Solver Init");
    topo_solver_init(Global.myID, Global.myMesh);
    Timer_Stop("Topography Solver Init");
    Timer_Reduce("Topography Solver Init", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
    if (Timer_Exists("Topo Init Select Elements")) {
      Timer_Reduce("Topo Init Select Elements", (TimerKind)(MAX | MIN | AVERAGE),
                   comm_solver);
    }
    if (Timer_Exists("Topo Init Build Node Sets")) {
      Timer_Reduce("Topo Init Build Node Sets", (TimerKind)(MAX | MIN | AVERAGE),
                   comm_solver);
    }
    if (Timer_Exists("Topo Init Build Constants")) {
      Timer_Reduce("Topo Init Build Constants", (TimerKind)(MAX | MIN | AVERAGE),
                   comm_solver);
    }
  }

  if (Param.theMeshOutFlag && DO_OUTPUT) {
    mesh_output();
  }

  if (Param.storeMeshCoordinatesForMatlab == YES) {
    saveMeshCoordinatesForMatlab(
        Global.myMesh, Global.myID, Param.parameters_input_file,
        Param.theMeshDirMatlab, Global.myOctree->ticksize,
        Param.theTypeOfDamping, Global.theXForMeshOrigin,
        Global.theYForMeshOrigin, Global.theZForMeshOrigin,
        Param.includeBuildings);
  }

  Timer_Start("Mesh Stats Print");
  mesh_print_stat(Global.myOctree, Global.myMesh, Global.myID,
                  Global.theGroupSize, Param.theMeshStatFilename);
  Timer_Stop("Mesh Stats Print");
  Timer_Reduce("Mesh Stats Print", (TimerKind)(MAX | MIN), comm_solver);

  if (Param.theMeshOnlyMode != 0) {
    if (Global.myID == 0) {
      printf("Mesh-only mode enabled: skipping solver phase after meshing.\n");
    }

    Timer_Stop("Total Wall Clock");
    Timer_Reduce("Total Wall Clock", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
    Timer_Reduce("Octor Newtree", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
    Timer_Reduce("Octor Refinetree", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
    if (Timer_Exists("Octor Balancetree")) {
      Timer_Reduce("Octor Balancetree", (TimerKind)(MAX | MIN | AVERAGE),
                   comm_solver);
    }
    Timer_Reduce("Octor Partitiontree", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
    Timer_Reduce("Octor Extractmesh", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
    if (Timer_Exists("Mesh correct properties")) {
      Timer_Reduce("Mesh correct properties", (TimerKind)(MAX | MIN | AVERAGE),
                   comm_solver);
    }

    if (Timer_Exists("Global Init")) {
      Timer_Reduce("Global Init", (TimerKind)(MAX | MIN | AVERAGE),
                   comm_solver);
    }
    if (Timer_Exists("Drm Init")) {
      Timer_Reduce("Drm Init", (TimerKind)(MAX | MIN | AVERAGE),
                   comm_solver);
    }
    if (Timer_Exists("Topography Solver Init")) {
      Timer_Reduce("Topography Solver Init", (TimerKind)(MAX | MIN | AVERAGE),
                   comm_solver);
    }
    if (Timer_Exists("Mesh Out")) {
      Timer_Reduce("Mesh Out", (TimerKind)(MAX | MIN | AVERAGE),
                   comm_solver);
    }

    if (Global.myID == 0) {
      printf("TOTAL WALL CLOCK                    : %.2f seconds\n",
             Timer_Value("Total Wall Clock", MAX));
      printf("Mesh-only mode complete.\n");
    }

    if (Global.myOctree != NULL) {
      octor_deletetree(Global.myOctree);
      Global.myOctree = NULL;
    }
    if (Global.myMesh != NULL) {
      octor_deletemesh(Global.myMesh);
      Global.myMesh = NULL;
    }

    MPI_Finalize();
    return 0;
  }

  Timer_Start("Pre Solver Setup");

  if (Param.theNumberOfStations != 0) {
    Timer_Start("Pre Setup Stations Init");
    output_stations_init(Param.parameters_input_file);
    Timer_Stop("Pre Setup Stations Init");
  }

  /* include additional info for topo stations */

  if ((Param.includeTopography == YES) && (Param.theNumberOfStations != 0)) {
    Timer_Start("Pre Setup Topo Stations Init");
    topography_stations_init(Global.myMesh, Param.myStations,
                             Param.myNumberOfStations);
    Timer_Stop("Pre Setup Topo Stations Init");
  }

  /* Initialize the output planes */
  if (Param.theNumberOfPlanes != 0) {
    double planesSetupStart = MPI_Wtime();
    double planesSetupElapsed;
    double planesSetupElapsedMax;

    Timer_Start("Pre Setup Planes Init");
    planes_setup(Global.myID, &Param.thePlanePrintRate, Param.thePlaneDirOut,
                 Param.IO_pool_pe_count, Param.theNumberOfPlanes,
                 Param.parameters_input_file, get_surface_shift(),
                 Param.theSurfaceCornersLong, Param.theSurfaceCornersLat,
                 Param.theDomainX, Param.theDomainY, Param.theDomainZ,
                 Param.planes_input_file, &Param.utmZone,
                 Param.reconstructAcousticPlaneDisplacements);
    Timer_Stop("Pre Setup Planes Init");

    planesSetupElapsed = MPI_Wtime() - planesSetupStart;
    MPI_Reduce(&planesSetupElapsed, &planesSetupElapsedMax, 1, MPI_DOUBLE,
               MPI_MAX, 0, comm_solver);
    if (Global.myID == 0) {
      monitor_print(
          "Plane setup complete: total=%d max-rank-time=%.2fs\n",
          Param.theNumberOfPlanes, planesSetupElapsedMax);
    }
  }

  /* Initialize the solver, source and output structures */
  Timer_Start("Pre Setup Solver Init");
  solver_init();
  Timer_Stop("Pre Setup Solver Init");
  Timer_Start("Solver Stats Print");
  solver_printstat(Global.mySolver);
  Timer_Stop("Solver Stats Print");
  Timer_Reduce("Solver Stats Print", (TimerKind)(MAX | MIN), comm_solver);

  /* Initialize nonlinear solver analysis structures */
  if (Param.includeNonlinearAnalysis == YES) {
    Timer_Start("Pre Setup Nonlinear Init");
    nonlinear_solver_init(Global.myID, Global.myMesh, Param.theDomainZ);
    if (Param.theNumberOfStations != 0) {
      nonlinear_stations_init(Global.myMesh, Param.myStations,
                              Param.myNumberOfStations);
    }
    nonlinear_stats(Global.myID, Global.theGroupSize);
    Timer_Stop("Pre Setup Nonlinear Init");
  }

  /* include additional info for topo stations */
  if (Param.includeTopography == YES) {
    Timer_Start("Pre Setup Topo Stats");
    topo_stats(Global.myID, Global.theGroupSize, Param.theTopoStatFilename);
    Timer_Stop("Pre Setup Topo Stats");
  }

  if (Param.includeIncidentPlaneWaves == YES) {
    Timer_Start("Pre Setup PlaneWaves Init");
    PlaneWaves_solver_init(Global.myID, Global.myMesh, Global.mySolver);
    Timer_Stop("Pre Setup PlaneWaves Init");
  }

  Timer_Start("Source Init");
  Timer_Start("Pre Setup Source Init");
  source_init(Param.parameters_input_file);
  Timer_Stop("Pre Setup Source Init");
  Timer_Stop("Source Init");
  Timer_Reduce("Source Init", (TimerKind)(MAX | MIN), comm_solver);

  /* Mapping element indices for stiffness
   * This is for compatibility with nonlinear
   * \TODO a more clever way should be possible
   */
  Timer_Start("Pre Setup Damping/Stiffness Init");
  if (Param.theTypeOfDamping < BKT) {
    stiffness_init(
        Global.myID,
        Global.myMesh,
        Global.mySolver);
  }
  // initialize linear elements only for Rayleigh damping
  else {
    damp_init(Global.myID, Global.myMesh);
  }
  Timer_Stop("Pre Setup Damping/Stiffness Init");

  /* this is a little too late to check for output parameters,
   * but let's do this in the mean time
   */
  Timer_Start("Pre Setup Output Init");
  output_init(Param.parameters_input_file, &Param.theOutputParameters);
  Timer_Stop("Pre Setup Output Init");

  Timer_Stop("Pre Solver Setup");
  Timer_Reduce("Pre Solver Setup", (TimerKind)(MAX | MIN | AVERAGE),
               comm_solver);
  if (Timer_Exists("Pre Setup Stations Init")) {
    Timer_Reduce("Pre Setup Stations Init", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Stations Read Input")) {
    Timer_Reduce("Stations Read Input", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Stations Broadcast")) {
    Timer_Reduce("Stations Broadcast", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Stations Locate")) {
    Timer_Reduce("Stations Locate", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Stations Local Init")) {
    Timer_Reduce("Stations Local Init", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Pre Setup Topo Stations Init")) {
    Timer_Reduce("Pre Setup Topo Stations Init",
                 (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }
  if (Timer_Exists("Pre Setup Planes Init")) {
    Timer_Reduce("Pre Setup Planes Init", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Pre Setup Solver Init")) {
    Timer_Reduce("Pre Setup Solver Init", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Pre Setup Nonlinear Init")) {
    Timer_Reduce("Pre Setup Nonlinear Init", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Pre Setup Topo Stats")) {
    Timer_Reduce("Pre Setup Topo Stats", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Pre Setup PlaneWaves Init")) {
    Timer_Reduce("Pre Setup PlaneWaves Init", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Pre Setup Damping/Stiffness Init")) {
    Timer_Reduce("Pre Setup Damping/Stiffness Init",
                 (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  }
  if (Timer_Exists("Pre Setup Output Init")) {
    Timer_Reduce("Pre Setup Output Init", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Pre Setup Source Init")) {
    Timer_Reduce("Pre Setup Source Init", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }

  /* Run the solver and output the results */
  Timer_Start("Phase Sync/Coordination");
  MPI_Barrier(comm_solver);
  Timer_Stop("Phase Sync/Coordination");
  Timer_Start("Solver");
  solver_run();
  Timer_Stop("Solver");
  Timer_Start("Phase Sync/Coordination");
  MPI_Barrier(comm_solver);
  Timer_Stop("Phase Sync/Coordination");

  Timer_Start("Post Solver Finalize");

  if (Param.includeNonlinearAnalysis == YES) {
    nonlinear_yield_stats(Global.myMesh, Global.myID, Param.theTotalSteps,
                          Global.theGroupSize);
  }

  output_fini();

#ifdef DEBUG
  /* Does the OS page out my resident set ? */
  if ((Global.myID % PROCPERNODE) == 0) {
    /* system("ps xl"); */
  }

  /* Are there pending messages that I haven't processed */
  MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm_solver, &flag, &status);
  if (flag != 0) {
    fprintf(stderr, "Thread %d: MPI error: unreceived incoming message\n",
            Global.myID);
  }
#endif /* DEBUG */

  local_finalize();

  output_get_stats();

  Timer_Stop("Post Solver Finalize");
  Timer_Reduce("Post Solver Finalize", (TimerKind)(MAX | MIN | AVERAGE),
               comm_solver);

  Timer_Start("Phase Sync/Coordination");
  MPI_Barrier(comm_solver);
  Timer_Stop("Phase Sync/Coordination");
  Timer_Stop("Total Wall Clock");

  Timer_Reduce("Total Wall Clock", (TimerKind)(MAX | MIN | AVERAGE),
               comm_solver);

  Timer_Reduce("Octor Newtree", (TimerKind)(MAX | MIN | AVERAGE),
               comm_solver);
  Timer_Reduce("Octor Refinetree", (TimerKind)(MAX | MIN | AVERAGE),
               comm_solver);
  if (Timer_Exists("Octor Balancetree")) {
    Timer_Reduce("Octor Balancetree", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  Timer_Reduce("Octor Partitiontree", (TimerKind)(MAX | MIN | AVERAGE),
               comm_solver);
  Timer_Reduce("Octor Extractmesh", (TimerKind)(MAX | MIN | AVERAGE),
               comm_solver);
  if (Timer_Exists("Mesh correct properties")) {
    Timer_Reduce("Mesh correct properties", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  Timer_Reduce("Solver", (TimerKind)(MAX | MIN | AVERAGE), comm_solver);
  if (Timer_Exists("Phase Sync/Coordination")) {
    Timer_Reduce("Phase Sync/Coordination", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }
  if (Timer_Exists("Carve Buildings")) {
    Timer_Reduce("Carve Buildings", (TimerKind)(MAX | MIN | AVERAGE),
                 comm_solver);
  }

  /* Print out the timing stat */
  if (Global.myID == 0) {
    print_timing_stat();
  }

  /* TODO: Think of a better place for this */
  if (Param.includeNonlinearAnalysis == YES) {
    if (get_geostatic_total_time() > 0) {
      check_balance(Global.myID);
    }
  }

  /* Send a message to IO pool PEs to close output files and goto here */
  if (Param.theNumberOfPlanes != 0) {
    planes_close(Global.myID, Param.IO_pool_pe_count, Param.theNumberOfPlanes);
  }
IO_PES_REJOIN:

  MPI_Finalize();

  return 0;
}
