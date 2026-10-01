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

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <strings.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <assert.h>
#include <errno.h>
#include <stdarg.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "io_planes.h"
#include "psolve.h"
#include "geometrics.h"
#include "quake_util.h"
#include "util.h"
#include "cvm.h"
#include "topography.h"
#include "timers.h"


extern int compute_csi_eta_dzeta( octant_t *octant, vector3D_t pointcoords,
				  vector3D_t *localcoords, int32_t *localNodeID );

#define LINESIZE       512

#define STOP_SERVER    10 /* DATA MESSAGE TO STOP IO SERVER */

#define MAX_STRIPS_PER_PLANE 5000
/*Could be determined dynmamically but would cost search_point calls*/


int  New_planes_print(int32_t myID, mysolver_t* mySolver, int theNumberOfPlanes);
void New_planes_setup(int32_t myID, int32_t *thePlanePrintRate, char *thePlaneDirOut, 
	int io_pool_pe_count,
    int theNumberOfPlanes, const char *numericalin, 
    double surfaceShift, double *theSurfaceCornersLong, double *theSurfaceCornersLat,
    double theDomainX, double theDomainY, double theDomainZ, char* planes_input_file,
	UTMZone_t* utmZone, int reconstructAcousticPlaneDisplacements);
void New_planes_close(int32_t myID, int theNumberOfPlanes);
void planes_IO_PES_main(int32_t myID);
static void New_output_planes_construct_strips(int32_t myID, int theNumberOfPlanes);

static void print_planeinfo(int32_t myID, int theNumberOfPlanes);
static void planes_select_output_mode(int32_t myID, int io_pool_pe_count);
static void planes_mpiio_setup(int32_t myID, int theNumberOfPlanes,
                               const char* planeDirOut);
static int planes_mpiio_flush_batch(void);
static void planes_mpiio_close(int theNumberOfPlanes);
static void planes_gpu_root_setup(int32_t myID, int theNumberOfPlanes);
static void planes_gpu_root_close(void);
static int planes_print_gpu_root(int32_t myID, int theNumberOfPlanes,
                                 const double* packedValues,
                                 int32_t pointCount);

int  Old_planes_print(int32_t myID, mysolver_t* mySolver, int theNumberOfPlanes);
void Old_planes_setup(int32_t myID, int32_t *thePlanePrintRate, const char* thePlaneDirOut, 
    int theNumberOfPlanes, const char *numericalin, 
    double surfaceShift, double *theSurfaceCornersLong, double *theSurfaceCornersLat,
    double theDomainX, double theDomainY, double theDomainZ, char* planes_input_file,
	UTMZone_t* utmZone, int reconstructAcousticPlaneDisplacements);
void Old_planes_close(int32_t myID, int theNumberOfPlanes);
static int  Old_print_plane_displacements(int32_t myID, int ThisPlane);
static void Old_output_planes_construct_strips(int32_t myID, int theNumberOfPlanes);

typedef struct plane_strip_element_t {

    int32_t nodestointerpolate[8];
	vector3D_t localcoords;  /* csi, eta, dzeta in (-1,1) x (-1,1) x (-1,1)*/
	int32_t point_is_valid;
	int32_t topostrip;
	int32_t owner_is_acoustic;
	int32_t pressure_is_defined;
	int32_t acoustic_kinematics_ready;
	double h;
	double rho;
	double acoustic_disp[3];
	double acoustic_vel[3];
	double acoustic_acc[3];
	double acoustic_prev_acc[3];

} plane_strip_element_t;


typedef struct plane_t {

    vector3D_t origincoords;              /**< cartesian */

    FILE *fpoutputfile,*fpplanecoordsfile;

    int numberofstepsalongstrike, numberofstepsdowndip;
    double stepalongstrike, stepdowndip, strike, dip, rake;

    int IO_PE_num;
    int numberofstripsthisplane;
    int globalnumberofstripsthisplane, globalnumberofstripsthisplane_recieved; /*valid only on output PE's*/
	int globalnumberofsendersthisplane;

    int stripstart[MAX_STRIPS_PER_PLANE];
    int stripend[MAX_STRIPS_PER_PLANE];
    plane_strip_element_t * strip[MAX_STRIPS_PER_PLANE];

} plane_t;

static plane_t* thePlanes;
static int planes_GlobalLargestStripCount;
static int planes_LocalLargestStripCount;
static int planes_GlobalLargestPlanePayloadDoubles;
static int planes_LocalLargestPlanePayloadDoubles;
static double* planes_output_buffer;
static double* planes_stripMPISendBuffer;
static double* planes_stripMPIRecvBuffer;
static int32_t g_planeIOServerCount = 1;
static int32_t g_planeIOServerFirstRank = 0;
static int32_t g_planeIOSolverRankCount = 0;
static int32_t g_planeOutputOmpMode = -1;
static int32_t g_planeOutputOmpReported = 0;
static int32_t g_planeOutputRankCached = -1;
static int32_t g_planePrintSequence = 0;
static int32_t g_reconstructAcousticPlaneDisplacements = 0;
static int32_t* g_planeDisplacementNodes = NULL;
static int32_t g_planeDisplacementNodeCount = 0;
static int32_t* g_acousticPlanePressureNodes = NULL;
static int32_t g_acousticPlanePressureNodeCount = 0;
static int32_t g_planeGpuOutputEnabled = 0;
static int32_t g_planeGpuRootEnabled = 0;
static int32_t g_planeGpuMpiioEnabled = 0;
static MPI_File* g_planeMpiFiles = NULL;
static MPI_Datatype* g_planeMpiFileTypes = NULL;
static int32_t* g_planeLocalPointOffsets = NULL;
static int32_t* g_planeLocalValueCounts = NULL;
static MPI_Offset* g_planeSnapshotBytes = NULL;
static double** g_planeMpiBatchBuffers = NULL;
static int32_t g_planeMpiBatchCapacity = 64;
static int32_t g_planeMpiBatchCount = 0;
static int32_t g_planeMpiBatchStartSequence = 0;
static int64_t g_planeMpiBatchFlushCount = 0;
static int32_t g_planeLocalPointCount = 0;
static int32_t g_planeMpiPlaneCount = 0;
static int32_t g_planeHostOnlyPointCount = 0;

static void planes_mpi_check(int errorCode, const char* operation) {
    if (errorCode == MPI_SUCCESS) {
        return;
    }

    char errorText[MPI_MAX_ERROR_STRING];
    int errorLength = 0;
    MPI_Error_string(errorCode, errorText, &errorLength);
    fprintf(stderr, "Plane MPI-IO %s failed: %.*s\n", operation,
            errorLength, errorText);
    MPI_Abort(comm_solver, errorCode);
    exit(1);
}

static int32_t planes_mpiio_batch_capacity(int32_t myID) {
    int32_t capacity = 64;

    if (myID == 0) {
        const char* env = getenv("HERCULES_PLANE_MPIIO_BATCH");
        if ((env != NULL) && (env[0] != '\0')) {
            char* end = NULL;
            long parsed;

            errno = 0;
            parsed = strtol(env, &end, 10);
            if ((errno != 0) || (end == env) || (*end != '\0') ||
                (parsed < 1) || (parsed > INT32_MAX)) {
                solver_abort("planes_mpiio_batch_capacity", NULL,
                             "Invalid HERCULES_PLANE_MPIIO_BATCH=%s; "
                             "expected an integer >= 1",
                             env);
            }
            capacity = (int32_t)parsed;
        }
    }

    MPI_Bcast(&capacity, 1, MPI_INT, 0, comm_solver);
    return capacity;
}

static void planes_select_output_mode(int32_t myID, int io_pool_pe_count) {
    const char* env = getenv("HERCULES_PLANE_IO_MODE");
    enum {
        PLANE_BACKEND_LEGACY = 0,
        PLANE_BACKEND_GPU_ROOT = 1,
        PLANE_BACKEND_GPU_MPIIO = 2
    };
    int requestedBackend = PLANE_BACKEND_GPU_MPIIO;

    if ((env == NULL) || (env[0] == '\0') ||
        (strcasecmp(env, "mpiio") == 0) ||
        (strcasecmp(env, "mpi-io") == 0)) {
        requestedBackend = PLANE_BACKEND_GPU_MPIIO;
    } else if ((strcasecmp(env, "hybrid") == 0) ||
               (strcasecmp(env, "gpu-root") == 0) ||
               (strcasecmp(env, "root") == 0) ||
               (strcasecmp(env, "gpu") == 0) ||
               (strcasecmp(env, "on") == 0) ||
               (strcasecmp(env, "yes") == 0) ||
               (strcmp(env, "1") == 0)) {
        requestedBackend = PLANE_BACKEND_GPU_ROOT;
    } else if ((strcasecmp(env, "legacy") == 0) ||
               (strcasecmp(env, "off") == 0) ||
               (strcasecmp(env, "no") == 0) ||
               (strcmp(env, "0") == 0)) {
        requestedBackend = PLANE_BACKEND_LEGACY;
    } else {
        if (myID == 0) {
            fprintf(stderr,
                    "Invalid HERCULES_PLANE_IO_MODE=%s; expected mpiio "
                    "(default), hybrid, or legacy\n",
                    env);
        }
        MPI_Abort(comm_solver, ERROR);
        exit(1);
    }

    if ((requestedBackend != PLANE_BACKEND_LEGACY) &&
        (io_pool_pe_count != 0)) {
        if (myID == 0) {
            fprintf(stderr,
                    "GPU plane-output backends are incompatible with "
                    "dedicated plane I/O ranks; using the existing I/O-pool "
                    "backend\n");
        }
        requestedBackend = PLANE_BACKEND_LEGACY;
    }

    g_planeGpuOutputEnabled =
        (requestedBackend != PLANE_BACKEND_LEGACY) ? 1 : 0;
    g_planeGpuRootEnabled =
        (requestedBackend == PLANE_BACKEND_GPU_ROOT) ? 1 : 0;
    g_planeGpuMpiioEnabled =
        (requestedBackend == PLANE_BACKEND_GPU_MPIIO) ? 1 : 0;
    g_planePrintSequence = 0;
    g_planeMpiBatchCount = 0;
    g_planeMpiBatchStartSequence = 0;
    g_planeMpiBatchFlushCount = 0;
    g_planeHostOnlyPointCount = 0;
    if (myID == 0) {
        const char* backendName = "legacy CPU interpolation + rank-0/I/O-pool";
        if (g_planeGpuRootEnabled) {
            backendName = "GPU interpolation + rank-0 gather/fwrite";
        } else if (g_planeGpuMpiioEnabled) {
            backendName = "GPU interpolation + MPI-IO";
        }
        printf("Plane output backend: %s\n", backendName);
    }
}

static void planes_mpiio_setup(int32_t myID, int theNumberOfPlanes,
                               const char* planeDirOut) {
    int64_t localPointCount = 0;
    int64_t localBatchBytes = 0;
    long long globalPointCount = 0;
    long long maxRankPointCount = 0;
    long long maxRankBatchBytes = 0;

    if (!g_planeGpuMpiioEnabled || (theNumberOfPlanes <= 0)) {
        return;
    }

    g_planeMpiBatchCapacity = planes_mpiio_batch_capacity(myID);
    g_planeMpiBatchCount = 0;
    g_planeMpiBatchStartSequence = 0;
    g_planeMpiBatchFlushCount = 0;

    g_planeMpiFiles = (MPI_File*)calloc((size_t)theNumberOfPlanes,
                                        sizeof(MPI_File));
    g_planeMpiFileTypes = (MPI_Datatype*)calloc((size_t)theNumberOfPlanes,
                                                sizeof(MPI_Datatype));
    g_planeLocalPointOffsets = (int32_t*)calloc(
        (size_t)theNumberOfPlanes + 1, sizeof(int32_t));
    g_planeLocalValueCounts = (int32_t*)calloc(
        (size_t)theNumberOfPlanes, sizeof(int32_t));
    g_planeSnapshotBytes = (MPI_Offset*)calloc(
        (size_t)theNumberOfPlanes, sizeof(MPI_Offset));
    g_planeMpiBatchBuffers = (double**)calloc(
        (size_t)theNumberOfPlanes, sizeof(double*));
    if ((g_planeMpiFiles == NULL) || (g_planeMpiFileTypes == NULL) ||
        (g_planeLocalPointOffsets == NULL) ||
        (g_planeLocalValueCounts == NULL) ||
        (g_planeSnapshotBytes == NULL) ||
        (g_planeMpiBatchBuffers == NULL)) {
        solver_abort("planes_mpiio_setup", NULL,
                     "Unable to allocate plane MPI-IO state");
    }

    for (int iPlane = 0; iPlane < theNumberOfPlanes; iPlane++) {
        char path[1024];
        const int stripCount = thePlanes[iPlane].numberofstripsthisplane;
        int* blockLengths = NULL;
        MPI_Aint* displacements = NULL;
        MPI_Datatype rawType = MPI_DATATYPE_NULL;
        MPI_Datatype fileType = MPI_DATATYPE_NULL;
        int64_t localValues64 = 0;
        const int64_t planePoints =
            (int64_t)thePlanes[iPlane].numberofstepsalongstrike *
            (int64_t)thePlanes[iPlane].numberofstepsdowndip;
        int64_t snapshotBytes64;

        if ((planePoints <= 0) ||
            (planePoints > INT64_MAX / (3 * (int64_t)sizeof(double))) ||
            (planePoints > INTPTR_MAX / (3 * (int64_t)sizeof(double)))) {
            solver_abort("planes_mpiio_setup", NULL,
                         "Plane %d dimensions exceed MPI-IO offset limits",
                         iPlane);
        }
        snapshotBytes64 = planePoints * 3 * (int64_t)sizeof(double);

        if (stripCount > 0) {
            blockLengths = (int*)malloc((size_t)stripCount * sizeof(int));
            displacements =
                (MPI_Aint*)malloc((size_t)stripCount * sizeof(MPI_Aint));
            if ((blockLengths == NULL) || (displacements == NULL)) {
                free(blockLengths);
                free(displacements);
                solver_abort("planes_mpiio_setup", NULL,
                             "Unable to allocate plane MPI file view");
            }

            for (int iStrip = 0; iStrip < stripCount; iStrip++) {
                const int stripPoints = thePlanes[iPlane].stripend[iStrip] -
                                        thePlanes[iPlane].stripstart[iStrip] +
                                        1;
                if ((stripPoints <= 0) || (stripPoints > (INT_MAX / 3)) ||
                    (localValues64 > (int64_t)INT_MAX -
                                         (int64_t)stripPoints * 3)) {
                    free(blockLengths);
                    free(displacements);
                    solver_abort("planes_mpiio_setup", NULL,
                                 "Plane %d local MPI-IO value count exceeds "
                                 "the MPI int-count limit",
                                 iPlane);
                }
                blockLengths[iStrip] = stripPoints * 3;
                displacements[iStrip] =
                    (MPI_Aint)thePlanes[iPlane].stripstart[iStrip] * 3 *
                    (MPI_Aint)sizeof(double);
                localValues64 += blockLengths[iStrip];
            }

            planes_mpi_check(
                MPI_Type_create_hindexed(stripCount, blockLengths,
                                         displacements, MPI_DOUBLE, &rawType),
                "filetype creation");
            planes_mpi_check(
                MPI_Type_create_resized(rawType, 0,
                                        (MPI_Aint)snapshotBytes64, &fileType),
                "filetype resize");
            planes_mpi_check(MPI_Type_commit(&fileType), "filetype commit");
            MPI_Type_free(&rawType);
            free(blockLengths);
            free(displacements);
        } else {
            /* This rank has no plane points. It still participates in every
             * collective write with count zero, using a harmless full view. */
            planes_mpi_check(MPI_Type_contiguous(1, MPI_DOUBLE, &fileType),
                             "empty filetype creation");
            planes_mpi_check(MPI_Type_commit(&fileType),
                             "empty filetype commit");
        }

        g_planeMpiFileTypes[iPlane] = fileType;
        g_planeLocalPointOffsets[iPlane] = (int32_t)localPointCount;
        g_planeLocalValueCounts[iPlane] = (int32_t)localValues64;
        g_planeSnapshotBytes[iPlane] = (MPI_Offset)snapshotBytes64;
        if ((localValues64 > 0) &&
            ((localValues64 > INT_MAX / g_planeMpiBatchCapacity) ||
             ((uint64_t)localValues64 >
              SIZE_MAX / sizeof(double) /
                  (uint64_t)g_planeMpiBatchCapacity))) {
            solver_abort("planes_mpiio_setup", NULL,
                         "Plane %d MPI-IO batch size %d exceeds count or "
                         "address-space limits on this rank",
                         iPlane, g_planeMpiBatchCapacity);
        }
        if (localValues64 > 0) {
            const size_t batchValues =
                (size_t)localValues64 * (size_t)g_planeMpiBatchCapacity;
            const size_t batchBytes = batchValues * sizeof(double);

            g_planeMpiBatchBuffers[iPlane] = (double*)malloc(batchBytes);
            if (g_planeMpiBatchBuffers[iPlane] == NULL) {
                solver_abort("planes_mpiio_setup", NULL,
                             "Unable to allocate %zu-byte MPI-IO batch "
                             "buffer for plane %d",
                             batchBytes, iPlane);
            }
            if ((uint64_t)batchBytes >
                (uint64_t)INT64_MAX - (uint64_t)localBatchBytes) {
                solver_abort("planes_mpiio_setup", NULL,
                             "Total local plane MPI-IO batch memory "
                             "exceeds INT64_MAX");
            }
            localBatchBytes += (int64_t)batchBytes;
        }
        localPointCount += localValues64 / 3;
        if (localPointCount > INT32_MAX) {
            solver_abort("planes_mpiio_setup", NULL,
                         "Local plane interpolation point count exceeds "
                         "INT32_MAX");
        }

        const int pathLength =
            snprintf(path, sizeof(path), "%s/planedisplacements.%d",
                     planeDirOut, iPlane);
        if ((pathLength < 0) || (pathLength >= (int)sizeof(path))) {
            solver_abort("planes_mpiio_setup", NULL,
                         "Plane output path is too long for plane %d", iPlane);
        }
        planes_mpi_check(
            MPI_File_open(comm_solver, path, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                          MPI_INFO_NULL, &g_planeMpiFiles[iPlane]),
            "file open");
        planes_mpi_check(
            MPI_File_set_errhandler(g_planeMpiFiles[iPlane], MPI_ERRORS_RETURN),
            "error-handler setup");
        planes_mpi_check(MPI_File_set_size(g_planeMpiFiles[iPlane], 0),
                         "file truncation");
        planes_mpi_check(
            MPI_File_set_view(g_planeMpiFiles[iPlane], 0, MPI_DOUBLE, fileType,
                              "native", MPI_INFO_NULL),
            "file-view setup");
    }
    g_planeLocalPointOffsets[theNumberOfPlanes] = (int32_t)localPointCount;
    g_planeLocalPointCount = (int32_t)localPointCount;
    g_planeMpiPlaneCount = theNumberOfPlanes;

    {
        const long long local = (long long)localPointCount;
        MPI_Reduce(&local, &globalPointCount, 1, MPI_LONG_LONG, MPI_SUM, 0,
                   comm_solver);
        MPI_Reduce(&local, &maxRankPointCount, 1, MPI_LONG_LONG, MPI_MAX, 0,
                   comm_solver);
        {
            const long long localBytes = (long long)localBatchBytes;
            MPI_Reduce(&localBytes, &maxRankBatchBytes, 1, MPI_LONG_LONG,
                       MPI_MAX, 0, comm_solver);
        }
    }
    if (myID == 0) {
        printf("Plane GPU/MPI-IO layout: %lld interpolation points globally "
               "(max/rank=%lld)\n",
               globalPointCount, maxRankPointCount);
        printf("Plane MPI-IO batching: %d snapshots/collective "
               "(max buffer/rank=%.2f MiB)\n",
               g_planeMpiBatchCapacity,
               (double)maxRankBatchBytes / (1024.0 * 1024.0));
    }
}

static int planes_mpiio_flush_batch(void) {
    static double emptyBuffer = 0.0;

    if (!g_planeGpuMpiioEnabled || (g_planeMpiFiles == NULL) ||
        (g_planeMpiBatchBuffers == NULL) || (g_planeMpiBatchCount < 0) ||
        (g_planeMpiBatchCount > g_planeMpiBatchCapacity)) {
        return 0;
    }
    if (g_planeMpiBatchCount == 0) {
        return 1;
    }
    if ((g_planeMpiBatchStartSequence + g_planeMpiBatchCount) !=
        g_planePrintSequence) {
        return 0;
    }

    for (int iPlane = 0; iPlane < g_planeMpiPlaneCount; iPlane++) {
        MPI_Status status;
        const int32_t localValues = g_planeLocalValueCounts[iPlane];
        const int64_t writeCount64 =
            (int64_t)localValues * (int64_t)g_planeMpiBatchCount;
        int writeCount;
        const MPI_Offset fileOffset =
            (MPI_Offset)g_planeMpiBatchStartSequence *
            (MPI_Offset)localValues;
        const double* buffer =
            (localValues > 0) ? g_planeMpiBatchBuffers[iPlane]
                              : &emptyBuffer;

        if ((writeCount64 < 0) || (writeCount64 > INT_MAX) ||
            ((localValues > 0) && (buffer == NULL))) {
            return 0;
        }
        writeCount = (int)writeCount64;
        planes_mpi_check(
            MPI_File_write_at_all(g_planeMpiFiles[iPlane], fileOffset,
                                  (void*)buffer, writeCount, MPI_DOUBLE,
                                  &status),
            "batched collective write");
        if (writeCount > 0) {
            int written = 0;
            planes_mpi_check(MPI_Get_count(&status, MPI_DOUBLE, &written),
                             "batched write count");
            if (written != writeCount) {
                solver_abort("planes_mpiio_flush_batch", NULL,
                             "Short plane MPI-IO batch write: expected %d "
                             "doubles, wrote %d",
                             writeCount, written);
            }
        }
    }

    g_planeMpiBatchCount = 0;
    g_planeMpiBatchStartSequence = g_planePrintSequence;
    g_planeMpiBatchFlushCount++;
    return 1;
}

int planes_flush_gpu_mpiio(void) {
    int result;

    if (!g_planeGpuMpiioEnabled || (g_planeMpiBatchCount == 0)) {
        return 1;
    }
    Timer_Start("Plane MPI-IO Write");
    result = planes_mpiio_flush_batch();
    Timer_Stop("Plane MPI-IO Write");
    return result;
}

static void planes_mpiio_close(int theNumberOfPlanes) {
    int solverRank = -1;

    if (!g_planeGpuMpiioEnabled || (g_planeMpiFiles == NULL)) {
        return;
    }

    if (!planes_mpiio_flush_batch()) {
        solver_abort("planes_mpiio_close", NULL,
                     "Unable to flush final plane MPI-IO batch");
    }
    MPI_Comm_rank(comm_solver, &solverRank);
    if (solverRank == 0) {
        printf("Plane MPI-IO completed: %" PRId64
               " collective batches for %d snapshots\n",
               g_planeMpiBatchFlushCount, g_planePrintSequence);
    }

    for (int iPlane = 0; iPlane < theNumberOfPlanes; iPlane++) {
        const MPI_Offset finalSize =
            (MPI_Offset)g_planePrintSequence * g_planeSnapshotBytes[iPlane];
        planes_mpi_check(MPI_File_set_size(g_planeMpiFiles[iPlane], finalSize),
                         "final file sizing");
        planes_mpi_check(MPI_File_close(&g_planeMpiFiles[iPlane]),
                         "file close");
        if (g_planeMpiFileTypes[iPlane] != MPI_DATATYPE_NULL) {
            MPI_Type_free(&g_planeMpiFileTypes[iPlane]);
        }
        free(g_planeMpiBatchBuffers[iPlane]);
    }

    free(g_planeMpiFiles);
    free(g_planeMpiFileTypes);
    free(g_planeLocalPointOffsets);
    free(g_planeLocalValueCounts);
    free(g_planeSnapshotBytes);
    free(g_planeMpiBatchBuffers);
    g_planeMpiFiles = NULL;
    g_planeMpiFileTypes = NULL;
    g_planeLocalPointOffsets = NULL;
    g_planeLocalValueCounts = NULL;
    g_planeSnapshotBytes = NULL;
    g_planeMpiBatchBuffers = NULL;
    g_planeMpiBatchCapacity = 64;
    g_planeMpiBatchCount = 0;
    g_planeMpiBatchStartSequence = 0;
    g_planeMpiBatchFlushCount = 0;
    g_planeLocalPointCount = 0;
    g_planeMpiPlaneCount = 0;
    g_planeHostOnlyPointCount = 0;
}

static void planes_gpu_root_setup(int32_t myID, int theNumberOfPlanes) {
    int64_t localPointCount = 0;
    long long globalPointCount = 0;
    long long maxRankPointCount = 0;

    if (!g_planeGpuRootEnabled || (theNumberOfPlanes <= 0)) {
        return;
    }

    g_planeLocalPointOffsets = (int32_t*)calloc(
        (size_t)theNumberOfPlanes + 1, sizeof(int32_t));
    g_planeLocalValueCounts = (int32_t*)calloc(
        (size_t)theNumberOfPlanes, sizeof(int32_t));
    if ((g_planeLocalPointOffsets == NULL) ||
        (g_planeLocalValueCounts == NULL)) {
        solver_abort("planes_gpu_root_setup", NULL,
                     "Unable to allocate hybrid plane-output layout");
    }

    for (int iPlane = 0; iPlane < theNumberOfPlanes; iPlane++) {
        int64_t localPlanePoints = 0;

        g_planeLocalPointOffsets[iPlane] = (int32_t)localPointCount;
        for (int iStrip = 0;
             iStrip < thePlanes[iPlane].numberofstripsthisplane; iStrip++) {
            const int stripPoints = thePlanes[iPlane].stripend[iStrip] -
                                    thePlanes[iPlane].stripstart[iStrip] + 1;
            if ((stripPoints <= 0) ||
                (localPlanePoints > INT32_MAX - (int64_t)stripPoints)) {
                solver_abort("planes_gpu_root_setup", NULL,
                             "Plane %d local point count exceeds INT32_MAX",
                             iPlane);
            }
            localPlanePoints += stripPoints;
        }
        if (localPointCount > INT32_MAX - localPlanePoints) {
            solver_abort("planes_gpu_root_setup", NULL,
                         "Local plane interpolation point count exceeds "
                         "INT32_MAX");
        }
        localPointCount += localPlanePoints;
        if (localPlanePoints > INT_MAX / 3) {
            solver_abort("planes_gpu_root_setup", NULL,
                         "Plane %d local packed-value count exceeds the "
                         "MPI int-count limit",
                         iPlane);
        }
        g_planeLocalValueCounts[iPlane] = (int32_t)(localPlanePoints * 3);
    }

    g_planeLocalPointOffsets[theNumberOfPlanes] = (int32_t)localPointCount;
    g_planeLocalPointCount = (int32_t)localPointCount;
    g_planeMpiPlaneCount = theNumberOfPlanes;

    {
        const long long local = (long long)localPointCount;
        MPI_Reduce(&local, &globalPointCount, 1, MPI_LONG_LONG, MPI_SUM, 0,
                   comm_solver);
        MPI_Reduce(&local, &maxRankPointCount, 1, MPI_LONG_LONG, MPI_MAX, 0,
                   comm_solver);
    }
    if (myID == 0) {
        printf("Plane GPU/root layout: %lld interpolation points globally "
               "(max/rank=%lld)\n",
               globalPointCount, maxRankPointCount);
    }
}

static void planes_gpu_root_close(void) {
    free(g_planeLocalPointOffsets);
    free(g_planeLocalValueCounts);
    g_planeLocalPointOffsets = NULL;
    g_planeLocalValueCounts = NULL;
    g_planeLocalPointCount = 0;
    g_planeMpiPlaneCount = 0;
    g_planeHostOnlyPointCount = 0;
}

static int compare_int32_plane_nodes(const void* lhs, const void* rhs) {
	const int32_t a = *(const int32_t*)lhs;
	const int32_t b = *(const int32_t*)rhs;
	return (a > b) - (a < b);
}

static void planes_reset_acoustic_pressure_nodes(void) {
	free(g_acousticPlanePressureNodes);
	g_acousticPlanePressureNodes = NULL;
	g_acousticPlanePressureNodeCount = 0;
}

static void planes_reset_displacement_nodes(void) {
	free(g_planeDisplacementNodes);
	g_planeDisplacementNodes = NULL;
	g_planeDisplacementNodeCount = 0;
}

static void planes_collect_displacement_nodes(int theNumberOfPlanes) {
	int capacity = 0;
	int count = 0;
	int uniqueCount = 0;
	int32_t* nodes = NULL;

	planes_reset_displacement_nodes();

	if ((theNumberOfPlanes <= 0) || (thePlanes == NULL)) {
		return;
	}

	for (int iPlane = 0; iPlane < theNumberOfPlanes; iPlane++) {
		for (int iStrip = 0; iStrip < thePlanes[iPlane].numberofstripsthisplane;
		     iStrip++) {
			const int elemCount =
				thePlanes[iPlane].stripend[iStrip] -
				thePlanes[iPlane].stripstart[iStrip] + 1;
			for (int elemnum = 0; elemnum < elemCount; elemnum++) {
				const plane_strip_element_t* element =
					&(thePlanes[iPlane].strip[iStrip][elemnum]);
				if (element->point_is_valid != 0) {
					capacity += 8;
				}
			}
		}
	}

	if (capacity <= 0) {
		return;
	}

	nodes = (int32_t*)malloc((size_t)capacity * sizeof(int32_t));
	if (nodes == NULL) {
		solver_abort("planes_collect_displacement_nodes", NULL,
		             "Failed allocating plane displacement-node list\n");
	}

	for (int iPlane = 0; iPlane < theNumberOfPlanes; iPlane++) {
		for (int iStrip = 0; iStrip < thePlanes[iPlane].numberofstripsthisplane;
		     iStrip++) {
			const int elemCount =
				thePlanes[iPlane].stripend[iStrip] -
				thePlanes[iPlane].stripstart[iStrip] + 1;
			for (int elemnum = 0; elemnum < elemCount; elemnum++) {
				const plane_strip_element_t* element =
					&(thePlanes[iPlane].strip[iStrip][elemnum]);
				if (element->point_is_valid == 0) {
					continue;
				}
				for (int corner = 0; corner < 8; corner++) {
					const int32_t node = element->nodestointerpolate[corner];
					if (node >= 0) {
						nodes[count++] = node;
					}
				}
			}
		}
	}

	if (count <= 0) {
		free(nodes);
		return;
	}

	qsort(nodes, (size_t)count, sizeof(int32_t), compare_int32_plane_nodes);
	for (int i = 0; i < count; i++) {
		if ((i == 0) || (nodes[i] != nodes[i - 1])) {
			nodes[uniqueCount++] = nodes[i];
		}
	}

	g_planeDisplacementNodes =
		(int32_t*)malloc((size_t)uniqueCount * sizeof(int32_t));
	if (g_planeDisplacementNodes == NULL) {
		free(nodes);
		solver_abort("planes_collect_displacement_nodes", NULL,
		             "Failed allocating plane displacement-node copy\n");
	}
	memcpy(g_planeDisplacementNodes, nodes,
	       (size_t)uniqueCount * sizeof(int32_t));
	g_planeDisplacementNodeCount = uniqueCount;
	free(nodes);
}

static void planes_collect_acoustic_pressure_nodes(int theNumberOfPlanes) {
	int capacity = 0;
	int count = 0;
	int uniqueCount = 0;
	int32_t* nodes = NULL;

	planes_reset_acoustic_pressure_nodes();

	if ((g_reconstructAcousticPlaneDisplacements == 0) ||
	    (theNumberOfPlanes <= 0) || (thePlanes == NULL)) {
		return;
	}

	for (int iPlane = 0; iPlane < theNumberOfPlanes; iPlane++) {
		for (int iStrip = 0; iStrip < thePlanes[iPlane].numberofstripsthisplane;
		     iStrip++) {
			const int elemCount =
				thePlanes[iPlane].stripend[iStrip] -
				thePlanes[iPlane].stripstart[iStrip] + 1;
			for (int elemnum = 0; elemnum < elemCount; elemnum++) {
				const plane_strip_element_t* element =
					&(thePlanes[iPlane].strip[iStrip][elemnum]);
				if (element->pressure_is_defined != 0) {
					capacity += 8;
				}
			}
		}
	}

	if (capacity <= 0) {
		return;
	}

	nodes = (int32_t*)malloc((size_t)capacity * sizeof(int32_t));
	if (nodes == NULL) {
		solver_abort("planes_collect_acoustic_pressure_nodes", NULL,
		             "Failed allocating acoustic plane pressure-node list\n");
	}

	for (int iPlane = 0; iPlane < theNumberOfPlanes; iPlane++) {
		for (int iStrip = 0; iStrip < thePlanes[iPlane].numberofstripsthisplane;
		     iStrip++) {
			const int elemCount =
				thePlanes[iPlane].stripend[iStrip] -
				thePlanes[iPlane].stripstart[iStrip] + 1;
			for (int elemnum = 0; elemnum < elemCount; elemnum++) {
				const plane_strip_element_t* element =
					&(thePlanes[iPlane].strip[iStrip][elemnum]);
				if (element->pressure_is_defined == 0) {
					continue;
				}
				for (int corner = 0; corner < 8; corner++) {
					const int32_t node = element->nodestointerpolate[corner];
					if (node >= 0) {
						nodes[count++] = node;
					}
				}
			}
		}
	}

	if (count <= 0) {
		free(nodes);
		return;
	}

	qsort(nodes, (size_t)count, sizeof(int32_t), compare_int32_plane_nodes);
	for (int i = 0; i < count; i++) {
		if ((i == 0) || (nodes[i] != nodes[i - 1])) {
			nodes[uniqueCount++] = nodes[i];
		}
	}

	g_acousticPlanePressureNodes =
		(int32_t*)malloc((size_t)uniqueCount * sizeof(int32_t));
	if (g_acousticPlanePressureNodes == NULL) {
		free(nodes);
		solver_abort("planes_collect_acoustic_pressure_nodes", NULL,
		             "Failed allocating acoustic plane pressure-node copy\n");
	}
	memcpy(g_acousticPlanePressureNodes, nodes,
	       (size_t)uniqueCount * sizeof(int32_t));
	g_acousticPlanePressureNodeCount = uniqueCount;
	free(nodes);
}

int32_t planes_displacement_node_count(void) {
	return g_planeDisplacementNodeCount;
}

const int32_t* planes_displacement_nodes(void) {
	return g_planeDisplacementNodes;
}

int32_t planes_acoustic_pressure_node_count(void) {
	return g_acousticPlanePressureNodeCount;
}

const int32_t* planes_acoustic_pressure_nodes(void) {
	return g_acousticPlanePressureNodes;
}

int32_t planes_gpu_output_enabled(void) {
	return g_planeGpuOutputEnabled;
}

int32_t planes_gpu_mpiio_enabled(void) {
	return g_planeGpuMpiioEnabled;
}

int32_t planes_local_interpolation_point_count(void) {
	return g_planeGpuOutputEnabled ? g_planeLocalPointCount : 0;
}

int32_t planes_host_only_point_count(void) {
	return g_planeGpuOutputEnabled ? g_planeHostOnlyPointCount : 0;
}

int planes_export_gpu_interpolation_metadata(
    int32_t pointCount, int32_t* nodes, double* coords, int32_t* modes) {
	int32_t packed = 0;
	g_planeHostOnlyPointCount = 0;

	if (!g_planeGpuOutputEnabled || (pointCount != g_planeLocalPointCount) ||
	    ((pointCount > 0) &&
	     ((nodes == NULL) || (coords == NULL) || (modes == NULL)))) {
		return 0;
	}

	for (int iPlane = 0; iPlane < g_planeMpiPlaneCount; iPlane++) {
		const int32_t planeEnd = g_planeLocalPointOffsets[iPlane + 1];
		for (int iStrip = 0;
		     iStrip < thePlanes[iPlane].numberofstripsthisplane;
		     iStrip++) {
			const int stripLength = thePlanes[iPlane].stripend[iStrip] -
			                        thePlanes[iPlane].stripstart[iStrip] + 1;
			for (int elemnum = 0; elemnum < stripLength; elemnum++) {
				const plane_strip_element_t* element =
					&(thePlanes[iPlane].strip[iStrip][elemnum]);
				int32_t mode = PLANE_GPU_POINT_INVALID;

				for (int corner = 0; corner < 8; corner++) {
					nodes[(size_t)packed * 8 + corner] =
						element->nodestointerpolate[corner];
				}
				coords[(size_t)packed * 3] = element->localcoords.x[0];
				coords[(size_t)packed * 3 + 1] = element->localcoords.x[1];
				coords[(size_t)packed * 3 + 2] = element->localcoords.x[2];

				if (element->point_is_valid != 0) {
					if ((g_reconstructAcousticPlaneDisplacements != 0) &&
					    (element->owner_is_acoustic != 0)) {
						mode = PLANE_GPU_POINT_HOST_ACOUSTIC;
						g_planeHostOnlyPointCount++;
					} else if (element->topostrip != 0) {
						mode = PLANE_GPU_POINT_TOPO;
					} else {
						mode = PLANE_GPU_POINT_TRILINEAR;
					}
				}
				modes[packed++] = mode;
			}
		}
		if (packed != planeEnd) {
			return 0;
		}
	}

	return packed == pointCount;
}

int planes_patch_host_only_values(double* packedValues,
                                  int32_t pointCount) {
	int32_t packed = 0;

	if (!g_planeGpuOutputEnabled || (pointCount != g_planeLocalPointCount) ||
	    ((pointCount > 0) && (packedValues == NULL))) {
		return 0;
	}

	for (int iPlane = 0; iPlane < g_planeMpiPlaneCount; iPlane++) {
		const int32_t planeEnd = g_planeLocalPointOffsets[iPlane + 1];
		for (int iStrip = 0;
		     iStrip < thePlanes[iPlane].numberofstripsthisplane;
		     iStrip++) {
			const int stripLength = thePlanes[iPlane].stripend[iStrip] -
			                        thePlanes[iPlane].stripstart[iStrip] + 1;
			for (int elemnum = 0; elemnum < stripLength; elemnum++) {
				const plane_strip_element_t* element =
					&(thePlanes[iPlane].strip[iStrip][elemnum]);
				double* const result = packedValues + (size_t)packed * 3;
				if (element->point_is_valid == 0) {
					result[0] = result[1] = result[2] = 0.0;
				} else if ((g_reconstructAcousticPlaneDisplacements != 0) &&
				           (element->owner_is_acoustic != 0)) {
					result[0] = element->acoustic_disp[0];
					result[1] = element->acoustic_disp[1];
					result[2] = element->acoustic_disp[2];
				}
				packed++;
			}
		}
		if (packed != planeEnd) {
			return 0;
		}
	}

	return packed == pointCount;
}

static int planes_print_gpu_mpiio(int32_t myID, int theNumberOfPlanes,
                                  const double* packedValues,
                                  int32_t pointCount) {
	(void)myID;

	if (!g_planeGpuMpiioEnabled || (g_planeMpiFiles == NULL) ||
	    (g_planeMpiBatchBuffers == NULL) ||
	    (g_planeMpiBatchCount < 0) ||
	    (g_planeMpiBatchCount >= g_planeMpiBatchCapacity) ||
	    (theNumberOfPlanes != g_planeMpiPlaneCount) ||
	    (pointCount != g_planeLocalPointCount) ||
	    ((pointCount > 0) && (packedValues == NULL))) {
		return 0;
	}

	Timer_Start("Plane MPI-IO Write");
	if (g_planeMpiBatchCount == 0) {
		g_planeMpiBatchStartSequence = g_planePrintSequence;
	}
	if ((g_planeMpiBatchStartSequence + g_planeMpiBatchCount) !=
	    g_planePrintSequence) {
		Timer_Stop("Plane MPI-IO Write");
		return 0;
	}
	for (int iPlane = 0; iPlane < theNumberOfPlanes; iPlane++) {
		const int32_t localValues = g_planeLocalValueCounts[iPlane];
		const int32_t pointOffset = g_planeLocalPointOffsets[iPlane];
		if (localValues > 0) {
			double* const batchDestination =
				g_planeMpiBatchBuffers[iPlane] +
				(size_t)g_planeMpiBatchCount * (size_t)localValues;
			const double* const planeValues =
				packedValues + (size_t)pointOffset * 3;
			memcpy(batchDestination, planeValues,
			       (size_t)localValues * sizeof(double));
		}
	}
	g_planeMpiBatchCount++;
	g_planePrintSequence++;
	if ((g_planeMpiBatchCount == g_planeMpiBatchCapacity) &&
	    !planes_mpiio_flush_batch()) {
		Timer_Stop("Plane MPI-IO Write");
		return 0;
	}
	Timer_Stop("Plane MPI-IO Write");
	return 1;
}

static int planes_print_gpu_root(int32_t myID, int theNumberOfPlanes,
                                 const double* packedValues,
                                 int32_t pointCount) {
    int32_t packedPoint = 0;

    if (!g_planeGpuRootEnabled ||
        (theNumberOfPlanes != g_planeMpiPlaneCount) ||
        (pointCount != g_planeLocalPointCount) ||
        ((pointCount > 0) && (packedValues == NULL))) {
        return 0;
    }

    for (int iPlane = 0; iPlane < theNumberOfPlanes; iPlane++) {
        MPI_Status status;

        if (packedPoint != g_planeLocalPointOffsets[iPlane]) {
            return 0;
        }
        if (myID == 0) {
            const int64_t planePointCount =
                (int64_t)thePlanes[iPlane].numberofstepsalongstrike *
                (int64_t)thePlanes[iPlane].numberofstepsdowndip;
            memset(planes_output_buffer, 0,
                   (size_t)planePointCount * 3 * sizeof(double));
        }

        Timer_Start("Plane Root Gather");
        for (int iStrip = 0;
             iStrip < thePlanes[iPlane].numberofstripsthisplane; iStrip++) {
            const int stripLength = thePlanes[iPlane].stripend[iStrip] -
                                    thePlanes[iPlane].stripstart[iStrip] + 1;
            const int valueCount = stripLength * 3;
            const double* stripValues =
                packedValues + (size_t)packedPoint * 3;

            if (myID == 0) {
                memcpy(planes_output_buffer +
                           (size_t)thePlanes[iPlane].stripstart[iStrip] * 3,
                       stripValues, (size_t)valueCount * sizeof(double));
            } else {
                memcpy(planes_stripMPISendBuffer, stripValues,
                       (size_t)valueCount * sizeof(double));
                planes_stripMPISendBuffer[valueCount] =
                    (double)thePlanes[iPlane].stripstart[iStrip] + 0.1;
                MPI_Send(planes_stripMPISendBuffer, valueCount + 1,
                         MPI_DOUBLE, 0, iPlane, comm_solver);
            }
            packedPoint += stripLength;
        }

        if (myID == 0) {
            const int remoteStripCount =
                thePlanes[iPlane].globalnumberofstripsthisplane -
                thePlanes[iPlane].numberofstripsthisplane;
            for (int rStrip = 0; rStrip < remoteStripCount; rStrip++) {
                int recvCount = 0;
                int startLocation;

                MPI_Recv(planes_stripMPIRecvBuffer,
                         planes_GlobalLargestStripCount * 3 + 1, MPI_DOUBLE,
                         MPI_ANY_SOURCE, iPlane, comm_solver, &status);
                MPI_Get_count(&status, MPI_DOUBLE, &recvCount);
                if ((recvCount <= 1) || (((recvCount - 1) % 3) != 0)) {
                    solver_abort("planes_print_gpu_root", NULL,
                                 "Invalid plane strip message length %d",
                                 recvCount);
                }
                startLocation =
                    (int)planes_stripMPIRecvBuffer[recvCount - 1];
                memcpy(planes_output_buffer + (size_t)startLocation * 3,
                       planes_stripMPIRecvBuffer,
                       (size_t)(recvCount - 1) * sizeof(double));
            }
        }
        Timer_Stop("Plane Root Gather");

        if (myID == 0) {
            Timer_Start("Plane Root fwrite");
            Old_print_plane_displacements(myID, iPlane);
            Timer_Stop("Plane Root fwrite");
        }

        /* Retain the synchronization semantics of the proven legacy writer. */
        MPI_Barrier(comm_solver);
        if (packedPoint != g_planeLocalPointOffsets[iPlane + 1]) {
            return 0;
        }
    }

    g_planePrintSequence++;
    return packedPoint == pointCount;
}

int planes_print_gpu(int32_t myID, int theNumberOfPlanes,
                     const double* packedValues, int32_t pointCount) {
    if (g_planeGpuRootEnabled) {
        return planes_print_gpu_root(myID, theNumberOfPlanes, packedValues,
                                     pointCount);
    }
    if (g_planeGpuMpiioEnabled) {
        return planes_print_gpu_mpiio(myID, theNumberOfPlanes, packedValues,
                                      pointCount);
    }
    return 0;
}

static void plane_unpack_message(const double* messageBuffer, int recvCount) {
	int offset = 0;

	while (offset < recvCount) {
		int startLocation;
		int stripLength;
		int valueCount;

		if ((offset + 2) > recvCount) {
			solver_abort("plane_unpack_message", NULL,
					 "Malformed plane message: recvCount=%d offset=%d",
					 recvCount, offset);
		}

		startLocation = (int)messageBuffer[offset++];
		stripLength = (int)messageBuffer[offset++];
		valueCount = stripLength * 3;

		if ((stripLength < 0) || ((offset + valueCount) > recvCount)) {
			solver_abort("plane_unpack_message", NULL,
					 "Malformed plane payload: stripLength=%d recvCount=%d offset=%d",
					 stripLength, recvCount, offset);
		}

		memcpy(&(planes_output_buffer[startLocation * 3]),
		       messageBuffer + offset,
		       (size_t)valueCount * sizeof(double));
		offset += valueCount;
	}
}

static int32_t plane_element_is_acoustic(const edata_t* edata) {
	return (edata != NULL) && (edata->Vp > 0.0) && (edata->rho > 0.0) &&
	       (edata->Vs <= ACOUSTIC_VS_EPS);
}

/* Optional setup-only audit of the actual containing mesh elements. A material
 * model lookup or distance to a domain face cannot replace CPML membership. */
void planes_write_point_metadata(int32_t myID, int theNumberOfPlanes,
                                 const char* planeDirOut) {
    const char* env = getenv("HERCULES_PLANE_METADATA");
    if ((env == NULL) || (env[0] == '\0') || (strcmp(env, "0") == 0)) {
        return;
    }
    if (strcmp(env, "1") != 0) {
        solver_abort("planes_write_point_metadata", NULL,
                     "HERCULES_PLANE_METADATA must be 0 or 1");
    }
    if ((theNumberOfPlanes == 0) || (thePlanes == NULL)) return;

    FILE* fp = NULL;
    for (int p = 0; p < theNumberOfPlanes; ++p) {
        const plane_t* plane = &thePlanes[p];
        for (int strip = 0; strip < plane->numberofstripsthisplane; ++strip) {
            for (int point = plane->stripstart[strip];
                 point <= plane->stripend[strip]; ++point) {
                const plane_strip_element_t* element =
                    &plane->strip[strip][point - plane->stripstart[strip]];
                const int i = point / plane->numberofstepsdowndip;
                const int j = point % plane->numberofstepsdowndip;
                vector3D_t local = {{i * plane->stepalongstrike,
                                     j * plane->stepdowndip, 0.0}};
                vector3D_t xyz = compute_global_coords(
                    plane->origincoords, local, plane->dip, 0, plane->strike);
                if (get_thebase_topo() != 0.0) {
                    xyz.x[2] += point_elevation(xyz.x[0], xyz.x[1]);
                }
                edata_t material;
                memset(&material, 0, sizeof(material));
                material.Vp = material.Vs = material.rho = material.edgesize = NAN;
                const char* medium = "unknown";
                int32_t eindex = -1;
                int cpml = -1;
                if (element->point_is_valid) {
                    octant_t* octant = NULL;
                    if (search_point(xyz, &octant) != 1 ||
                        (eindex = solver_get_octant_material(octant, &material)) < 0) {
                        solver_abort("planes_write_point_metadata", NULL,
                                     "Cannot resolve plane %d point (%d,%d)", p, i, j);
                    }
                    medium = plane_element_is_acoustic(&material) ? "acoustic" :
                        ((material.Vp > 0.0 && material.rho > 0.0) ? "solid" : "inactive");
                    cpml = cpml_pml_element(eindex);
                }
                if (fp == NULL) {
                    char path[1024];
                    int length = snprintf(path, sizeof(path),
                                          "%s/plane_point_metadata.%d.csv", planeDirOut, myID);
                    if (length < 0 || (size_t)length >= sizeof(path)) {
                        solver_abort("planes_write_point_metadata", NULL,
                                     "Plane metadata output path is too long");
                    }
                    fp = hu_fopen(path, "w");
                    fprintf(fp, "plane,i,j,rank,element,x_m,y_m,z_m,medium,cpml,"
                                "vp_m_s,vs_m_s,rho_kg_m3,h_m,valid,topostrip,"
                                "acoustic_reconstruction,pressure_defined\n");
                }
                fprintf(fp, "%d,%d,%d,%d,%d,%.17g,%.17g,%.17g,%s,%d,"
                            "%.17g,%.17g,%.17g,%.17g,%d,%d,%d,%d\n",
                        p, i, j, myID, eindex, xyz.x[0], xyz.x[1], xyz.x[2], medium, cpml,
                        (double)material.Vp, (double)material.Vs, (double)material.rho,
                        (double)material.edgesize, element->point_is_valid, element->topostrip,
                        g_reconstructAcousticPlaneDisplacements && element->owner_is_acoustic,
                        element->pressure_is_defined);
            }
        }
    }
    if (fp != NULL) {
        const int writeError = ferror(fp);
        const int closeError = fclose(fp);
        if (writeError || closeError) {
            solver_abort("planes_write_point_metadata", NULL,
                         "Failed writing plane point metadata on rank %d", myID);
        }
    }
    if (myID == 0) {
        printf("Plane point metadata: %s/plane_point_metadata.RANK.csv "
               "(setup only; containing material and CPML membership)\n", planeDirOut);
    }
}

static void plane_init_acoustic_metadata(plane_strip_element_t* element,
		octant_t* octant) {
	edata_t material;

	if (element == NULL) {
		return;
	}

	element->owner_is_acoustic = 0;
	element->pressure_is_defined = 0;
	element->acoustic_kinematics_ready = 0;
	element->h = 0.0;
	element->rho = 0.0;
	memset(element->acoustic_disp, 0, sizeof(element->acoustic_disp));
	memset(element->acoustic_vel, 0, sizeof(element->acoustic_vel));
	memset(element->acoustic_acc, 0, sizeof(element->acoustic_acc));
	memset(element->acoustic_prev_acc, 0, sizeof(element->acoustic_prev_acc));

	if ((g_reconstructAcousticPlaneDisplacements == 0) ||
	    (solver_get_octant_material(octant, &material) < 0) ||
	    !plane_element_is_acoustic(&material)) {
		return;
	}

	element->owner_is_acoustic = 1;

	/* Topography tetrahedra do not have the full hex pressure stencil needed
	 * by plane_pressure_gradient(); keep them from falling back to tm1. */
	if (element->topostrip != 0) {
		return;
	}

	element->pressure_is_defined = 1;
	element->h = material.edgesize;
	element->rho = material.rho;
}

static int plane_pressure_gradient(const mysolver_t* mySolver,
		const plane_strip_element_t* element,
		double grad[3]) {
	int iPhi;
	const double sx[8] = {-1, 1, -1, 1, -1, 1, -1, 1};
	const double sy[8] = {-1, -1, 1, 1, -1, -1, 1, 1};
	const double sz[8] = {-1, -1, -1, -1, 1, 1, 1, 1};
	const double h = element->h;
	const vector3D_t localCoords = element->localcoords;
	const double xi = localCoords.x[0];
	const double eta = localCoords.x[1];
	const double zeta = localCoords.x[2];

	grad[0] = 0.0;
	grad[1] = 0.0;
	grad[2] = 0.0;

	if ((element->pressure_is_defined == 0) || (mySolver == NULL) ||
	    (mySolver->pressure_tm1 == NULL) || (h <= 0.0)) {
		return 0;
	}

	for (iPhi = 0; iPhi < 8; iPhi++) {
		const int32_t node = element->nodestointerpolate[iPhi];
		const double pressure = mySolver->pressure_tm1[node];
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

static void plane_update_acoustic_point(const mysolver_t* mySolver,
		plane_strip_element_t* element,
		int32_t step,
		double deltaT,
		double deltaTSquared) {
	double grad[3];
	double accel[3];

	if ((element == NULL) || (element->pressure_is_defined == 0) ||
	    (element->rho <= 0.0) ||
	    (plane_pressure_gradient(mySolver, element, grad) == 0)) {
		return;
	}

	for (int axis = 0; axis < 3; axis++) {
		accel[axis] = -grad[axis] / element->rho;
	}

	if ((element->acoustic_kinematics_ready == 0) || (step == 0)) {
		for (int axis = 0; axis < 3; axis++) {
			element->acoustic_acc[axis] = accel[axis];
			element->acoustic_prev_acc[axis] = accel[axis];
		}
		element->acoustic_kinematics_ready = 1;
		return;
	}

	for (int axis = 0; axis < 3; axis++) {
		const double oldVel = element->acoustic_vel[axis];
		const double oldAcc = element->acoustic_prev_acc[axis];

		element->acoustic_disp[axis] +=
			oldVel * deltaT + 0.5 * oldAcc * deltaTSquared;
		element->acoustic_vel[axis] +=
			0.5 * (oldAcc + accel[axis]) * deltaT;
		element->acoustic_acc[axis] = accel[axis];
		element->acoustic_prev_acc[axis] = accel[axis];
	}
}

static void plane_interpolate_point(const mysolver_t* mySolver,
		const plane_strip_element_t* element,
		double* outX,
		double* outY,
		double* outZ) {
	int iPhi;
	int32_t nodesToInterpolate[8];
	vector3D_t localCoords;
	double phi[8];
	static const double xi[3][8] = {
		{-1,  1, -1,  1, -1,  1, -1, 1},
		{-1, -1,  1,  1, -1, -1,  1, 1},
		{-1, -1, -1, -1,  1,  1,  1, 1}
	};
	double displacementsX = 0.0;
	double displacementsY = 0.0;
	double displacementsZ = 0.0;

	if ((element == NULL) || (element->point_is_valid == 0)) {
		*outX = 0.0;
		*outY = 0.0;
		*outZ = 0.0;
		return;
	}

	if ((g_reconstructAcousticPlaneDisplacements != 0) &&
	    (element->owner_is_acoustic != 0) &&
	    (mySolver != NULL) && (mySolver->pressure_tm1 != NULL)) {
		*outX = element->acoustic_disp[0];
		*outY = element->acoustic_disp[1];
		*outZ = element->acoustic_disp[2];
		return;
	}

	localCoords = element->localcoords;
	for (iPhi = 0; iPhi < 8; iPhi++) {
		nodesToInterpolate[iPhi] = element->nodestointerpolate[iPhi];
	}

	if (element->topostrip == 1) {
		double ux_0 = mySolver->tm1[nodesToInterpolate[0]].f[0];
		double ux_10 = mySolver->tm1[nodesToInterpolate[1]].f[0] - ux_0;
		double ux_20 = mySolver->tm1[nodesToInterpolate[2]].f[0] - ux_0;
		double ux_30 = mySolver->tm1[nodesToInterpolate[3]].f[0] - ux_0;
		double uy_0 = mySolver->tm1[nodesToInterpolate[0]].f[1];
		double uy_10 = mySolver->tm1[nodesToInterpolate[1]].f[1] - uy_0;
		double uy_20 = mySolver->tm1[nodesToInterpolate[2]].f[1] - uy_0;
		double uy_30 = mySolver->tm1[nodesToInterpolate[3]].f[1] - uy_0;
		double uz_0 = mySolver->tm1[nodesToInterpolate[0]].f[2];
		double uz_10 = mySolver->tm1[nodesToInterpolate[1]].f[2] - uz_0;
		double uz_20 = mySolver->tm1[nodesToInterpolate[2]].f[2] - uz_0;
		double uz_30 = mySolver->tm1[nodesToInterpolate[3]].f[2] - uz_0;

		displacementsX = ux_0 + localCoords.x[0] * ux_20
						  + localCoords.x[1] * ux_10
						  + localCoords.x[2] * ux_30;
		displacementsY = uy_0 + localCoords.x[0] * uy_20
						  + localCoords.x[1] * uy_10
						  + localCoords.x[2] * uy_30;
		displacementsZ = uz_0 + localCoords.x[0] * uz_20
						  + localCoords.x[1] * uz_10
						  + localCoords.x[2] * uz_30;
	} else {
		for (iPhi = 0; iPhi < 8; iPhi++) {
			phi[iPhi] = (1 + xi[0][iPhi] * localCoords.x[0])
					  * (1 + xi[1][iPhi] * localCoords.x[1])
					  * (1 + xi[2][iPhi] * localCoords.x[2]) / 8;
			displacementsX += phi[iPhi] * mySolver->tm1[nodesToInterpolate[iPhi]].f[0];
			displacementsY += phi[iPhi] * mySolver->tm1[nodesToInterpolate[iPhi]].f[1];
			displacementsZ += phi[iPhi] * mySolver->tm1[nodesToInterpolate[iPhi]].f[2];
		}
	}

	*outX = displacementsX;
	*outY = displacementsY;
	*outZ = displacementsZ;
}

static int32_t plane_output_use_omp(void) {
#ifndef _OPENMP
	return 0;
#else
	if (g_planeOutputOmpMode < 0) {
		const char* env = getenv("HERCULES_PLANE_IO_OMP");
		if ((env != NULL) && (strcmp(env, "0") == 0)) {
			g_planeOutputOmpMode = 0;
		} else if (env != NULL) {
			g_planeOutputOmpMode = 1;
		} else {
			g_planeOutputOmpMode = 1;
		}

		if ((g_planeOutputOmpMode == 1) && (g_planeOutputOmpReported == 0)) {
			if (g_planeOutputRankCached < 0) {
				MPI_Comm_rank(MPI_COMM_WORLD, &g_planeOutputRankCached);
			}
			if (g_planeOutputRankCached == 0) {
				fprintf(stderr,
						"Plane output OpenMP: enabled (threads=%d)\n",
						omp_get_max_threads());
			}
			g_planeOutputOmpReported = 1;
		}
	}

	return g_planeOutputOmpMode;
#endif
}

/******* Wrapper on routines to maintain backwards compatiblity for single core routines ******/

int planes_print(int32_t myID, int IO_pool_pe_count, int theNumberOfPlanes, mysolver_t* mySolver){
    if (IO_pool_pe_count)
         New_planes_print(myID, mySolver, theNumberOfPlanes);
    else
         Old_planes_print(myID, mySolver, theNumberOfPlanes);
    return 1;
}

void planes_update_acoustic_kinematics(mysolver_t* mySolver,
    int theNumberOfPlanes, int32_t step, double deltaT,
    double deltaTSquared) {
	int iPlane, iStrip;

	if ((g_reconstructAcousticPlaneDisplacements == 0) ||
	    (theNumberOfPlanes == 0) || (thePlanes == NULL) ||
	    (mySolver == NULL) || (mySolver->pressure_tm1 == NULL)) {
		return;
	}

	for (iPlane = 0; iPlane < theNumberOfPlanes; iPlane++) {
		for (iStrip = 0; iStrip < thePlanes[iPlane].numberofstripsthisplane;
		     iStrip++) {
			int elemnum;
			const int elemCount =
				thePlanes[iPlane].stripend[iStrip] -
				thePlanes[iPlane].stripstart[iStrip] + 1;

	#ifdef _OPENMP
			const int32_t useOmp = plane_output_use_omp();
	#pragma omp parallel for if(useOmp) schedule(static)
	#endif
			for (elemnum = 0; elemnum < elemCount; elemnum++) {
				plane_update_acoustic_point(mySolver,
					&(thePlanes[iPlane].strip[iStrip][elemnum]),
					step, deltaT, deltaTSquared);
			}
		}
	}
}

void planes_setup(int32_t myID, int32_t *thePlanePrintRate, char *thePlaneDirOut, 
    int IO_pool_pe_count, int theNumberOfPlanes, const char *numericalin, 
    double surfaceShift, double *theSurfaceCornersLong, double *theSurfaceCornersLat,
    double theDomainX, double theDomainY, double theDomainZ, char* planes_input_file,
	UTMZone_t* utmZone, int reconstructAcousticPlaneDisplacements){

//	IO_pool_pe_count = 1;
	planes_select_output_mode(myID, IO_pool_pe_count);
	g_reconstructAcousticPlaneDisplacements =
		reconstructAcousticPlaneDisplacements ? 1 : 0;
	planes_reset_displacement_nodes();
	planes_reset_acoustic_pressure_nodes();
    if (IO_pool_pe_count)
	New_planes_setup(myID, thePlanePrintRate, thePlaneDirOut,
		IO_pool_pe_count,
        theNumberOfPlanes, numericalin, 
        surfaceShift, theSurfaceCornersLong, theSurfaceCornersLat,
        theDomainX, theDomainY, theDomainZ, planes_input_file, utmZone,
        reconstructAcousticPlaneDisplacements);
    else
    Old_planes_setup(myID, thePlanePrintRate, thePlaneDirOut, 
        theNumberOfPlanes, numericalin, 
        surfaceShift, theSurfaceCornersLong, theSurfaceCornersLat,
        theDomainX, theDomainY, theDomainZ, planes_input_file, utmZone,
        reconstructAcousticPlaneDisplacements);

	if (!g_planeGpuOutputEnabled) {
		planes_collect_displacement_nodes(theNumberOfPlanes);
	}
	planes_collect_acoustic_pressure_nodes(theNumberOfPlanes);
}

void planes_close(int32_t myID, int IO_pool_pe_count, int theNumberOfPlanes){
    if (IO_pool_pe_count)
	New_planes_close(myID, theNumberOfPlanes);
	else
	Old_planes_close(myID, theNumberOfPlanes);
	planes_reset_displacement_nodes();
	planes_reset_acoustic_pressure_nodes();
}


/********************* Old Version To Maintain Backwards Compatbility For One Core ******/

/**
 * Interpolate the displacenents and communicate strips to printing PE.
 * Bigben version has no control flow, just chugs through planes.
 */
int Old_planes_print(int32_t myID, mysolver_t* mySolver, int theNumberOfPlanes)
{
	int iPlane;
    int stripLength;
    int iStrip, elemnum, rStrip;
    //    int32_t     howManyDisplacements, iDisplacement;
    //    vector3D_t* localCoords; /* convinient renaming */
    //    int32_t*    nodesToInterpolate[8];
    MPI_Status  status;
    //    MPI_Request sendstat;
    int recvStripCount, StartLocation;
#ifdef _OPENMP
	const int32_t useOmp = plane_output_use_omp();
#endif

    for (iPlane = 0; iPlane < theNumberOfPlanes ;iPlane++) {
	if (myID == 0) {
	    int planeValues = 3 * thePlanes[iPlane].numberofstepsalongstrike
		* thePlanes[iPlane].numberofstepsdowndip;
	    memset(planes_output_buffer, 0, (size_t)planeValues * sizeof(double));
	}
	for (iStrip = 0; iStrip < thePlanes[iPlane].numberofstripsthisplane;
	     iStrip++) {
	    int elemCount = thePlanes[iPlane].stripend[iStrip]
	        - thePlanes[iPlane].stripstart[iStrip] + 1;
	    double* stripBuffer = planes_stripMPISendBuffer;

	    /* interpolate points directly into send buffer for this strip */
	#ifdef _OPENMP
	#pragma omp parallel for if(useOmp) schedule(static)
	#endif
	    for (elemnum = 0; elemnum < elemCount; elemnum++)
	    {
	    	double displacementsX;
	    	double displacementsY;
	    	double displacementsZ;

	    	plane_interpolate_point(mySolver,
	    			&(thePlanes[iPlane].strip[iStrip][elemnum]),
	    			&displacementsX,
	    			&displacementsY,
	    			&displacementsZ);

	    	stripBuffer[elemnum*3] = displacementsX;
	    	stripBuffer[elemnum*3 + 1] = displacementsY;
	    	stripBuffer[elemnum*3 + 2] = displacementsZ;

	    }

	    /* add start location as last element, add .1 to insure
	       int-float conversion (yes, bad form)*/
            stripLength = thePlanes[iPlane].stripend[iStrip]
		- thePlanes[iPlane].stripstart[iStrip] + 1;

		    stripBuffer[stripLength*3]
		= (double) (thePlanes[iPlane].stripstart[iStrip] + 0.1);

            if (myID==0) { /* don't try to send to same PE, just memcopy */
			memcpy( &(planes_output_buffer[ thePlanes[iPlane].stripstart[iStrip]*3 ] ), stripBuffer,
			stripLength*3*sizeof(double) );
            }
            else {
			MPI_Send( stripBuffer, (stripLength*3)+1,
			  MPI_DOUBLE, 0 , iPlane, comm_solver );
            }
	} /*for iStrips*/

	/* IO PE recieves for this plane */
	if(myID == 0){
	    /* Recv all strips not directly memory copied above */
	    for (rStrip = 0;
		 rStrip < (thePlanes[iPlane].globalnumberofstripsthisplane
			   - thePlanes[iPlane].numberofstripsthisplane);
		 rStrip++) {

		MPI_Recv( planes_stripMPIRecvBuffer,
			  planes_GlobalLargestStripCount * 3 + 1, MPI_DOUBLE,
			  MPI_ANY_SOURCE, iPlane, comm_solver, &status );
		MPI_Get_count(&status, MPI_DOUBLE, &recvStripCount);
		StartLocation = (int)planes_stripMPIRecvBuffer[recvStripCount-1];
		memcpy( &(planes_output_buffer[StartLocation * 3]),
			planes_stripMPIRecvBuffer,
			(recvStripCount-1)*sizeof(double) );
	    }
	}

	/* do print for just this plane */
	Old_print_plane_displacements(myID, iPlane );

	/* may not be required */
	MPI_Barrier( comm_solver );
    } /* for iPlane */

    return 1;
}


static int Old_print_plane_displacements(int32_t myID, int ThisPlane)
{
    int ret;

    if  (myID == 0){
	ret = fwrite( planes_output_buffer, sizeof(double),
                      3 * thePlanes[ThisPlane].numberofstepsalongstrike
		      * thePlanes[ThisPlane].numberofstepsdowndip,
		      thePlanes[ThisPlane].fpoutputfile );
	if (ret != (3*thePlanes[ThisPlane].numberofstepsalongstrike
		    * thePlanes[ThisPlane].numberofstepsdowndip) ) {
	    fprintf( stderr, "Error writing all values in planes file %d.\n",
		     ThisPlane );
	    exit(1);
	}
    }

    return 1;
}



/**
 * This is the stripped BigBen version.  It is simplified:
 * 1. It assumes all output is PE0.
 * 2. It does not have a strip size limit.
 */

void Old_planes_setup(int32_t myID, int32_t *thePlanePrintRate, const char *thePlaneDirOut,
    int theNumberOfPlanes, const char *numericalin,
    double surfaceShift, double *theSurfaceCornersLong, double *theSurfaceCornersLat,
    double theDomainX, double theDomainY, double theDomainZ, char* planes_input_file,
	UTMZone_t* utmZone, int reconstructAcousticPlaneDisplacements) {

    static const char* fname = "output_planes_setup()";
    char       planedisplacementsout[1024]; 
	// char       planecoordsfile[1024];
    int        iPlane; 
    vector3D_t originPlaneCoords;
    int        largestplanesize;
	FILE*      fp;
	FILE*      fp_planes;

    g_reconstructAcousticPlaneDisplacements =
        reconstructAcousticPlaneDisplacements ? 1 : 0;

    if (myID == 0) {

	/* obtain the general planes specifications in parameter file */
	if ( (fp = fopen ( numericalin, "r")) == NULL ) {
	    solver_abort (fname, numericalin,
			  "Error opening parameters configuration file");
	}
	if ( (parsetext(fp, "output_planes_print_rate", 'i',
			thePlanePrintRate) != 0) )
	    {
		solver_abort (fname, NULL, "Error parsing output_planes_print_rate field from %s\n",
			      numericalin);
	    }

	thePlanes = (plane_t *) malloc ( sizeof( plane_t ) * theNumberOfPlanes);
	if ( thePlanes == NULL ) {
	    solver_abort( fname, "Allocating memory for planes array",
			  "Unable to create plane information arrays" );
	}
	fclose(fp);


	/* read in the individual planes coords in planes file */
	if ( (fp_planes = fopen ( planes_input_file, "r")) == NULL ) {
	    solver_abort (fname, planes_input_file,
			  "Error opening planes configuration file");
	}

	int found=0;
	while (!found) {
	    char line[LINESIZE];
	    char *name;
	    char delimiters[] = " =\n";
	    char querystring[] = "output_planes";
	    /* Read in one line */
	    if (fgets (line, LINESIZE, fp_planes) == NULL) {
		break;
	    }
	    name = strtok(line, delimiters);
	    if ( (name != NULL) && (strcmp(name, querystring) == 0) ) {
		largestplanesize = 0;
		found = 1;
		for ( iPlane = 0; iPlane < theNumberOfPlanes; iPlane++ ) {
		    int fret0, fret1;
		    FILE* fp0;
		    // FILE* fp1;
		    fret0 = fscanf( fp_planes," %lf %lf %lf",
				    &thePlanes[iPlane].origincoords.x[0],
				    &thePlanes[iPlane].origincoords.x[1],
				    &thePlanes[iPlane].origincoords.x[2] );

		    /* RICARDO: Correction for buildings case */
		    thePlanes[iPlane].origincoords.x[2] += surfaceShift;

		    /* origin coordinates must be in Lat, Long and depth */
		    if (fret0 == 0) {
			solver_abort (fname, NULL,
				      "Unable to read planes origin in %s",
				      planes_input_file);
		    }
		    /* convert to cartesian refered to the mesh */
                    originPlaneCoords =
			compute_domain_coords_linearinterp(
				thePlanes[ iPlane ].origincoords.x[1], thePlanes[ iPlane ].origincoords.x[0],
				theSurfaceCornersLong, theSurfaceCornersLat,
				theDomainY, theDomainX, utmZone);

                    thePlanes[iPlane].origincoords.x[0]=originPlaneCoords.x[0];
                    thePlanes[iPlane].origincoords.x[1]=originPlaneCoords.x[1];

                    fret1 = fscanf(fp_planes," %lf %d %lf %d %lf %lf",
                                   &thePlanes[iPlane].stepalongstrike,
                                   &thePlanes[iPlane].numberofstepsalongstrike,
                                   &thePlanes[iPlane].stepdowndip,
				   &thePlanes[iPlane].numberofstepsdowndip,
                                   &thePlanes[iPlane].strike,
                                   &thePlanes[iPlane].dip);
                    /*Find largest plane for output buffer allocation */
                    if ( (thePlanes[iPlane].numberofstepsdowndip
                          * thePlanes[iPlane].numberofstepsalongstrike)
                         > largestplanesize)
			{
			    largestplanesize
				= thePlanes[iPlane].numberofstepsdowndip
				* thePlanes[iPlane].numberofstepsalongstrike;
			}
                    if (fret1 == 0) {
			solver_abort(fname, NULL,
				     "Unable to read plane specification in %s",
				     planes_input_file);
                    }

                    /* The scalable path opens these files collectively after
                     * strip ownership is known. The legacy path retains the
                     * rank-0 stdio writer and its byte-for-byte file layout. */
                    if (!g_planeGpuMpiioEnabled) {
                        sprintf(planedisplacementsout,
                                "%s/planedisplacements.%d", thePlaneDirOut,
                                iPlane);
                        fp0 = hu_fopen(planedisplacementsout, "w");
                        thePlanes[iPlane].fpoutputfile = fp0;
                    } else {
                        thePlanes[iPlane].fpoutputfile = NULL;
                    }

                    // NOTE: Commented out by Clifford: Since planecoords files are not implemented, the files should not be created in the first place
					// sprintf (planecoordsfile, "%s/planecoords.%d",
                    //          thePlaneDirOut, iPlane);
                    // fp1 = hu_fopen (planecoordsfile, "w");
                    // thePlanes[iPlane].fpplanecoordsfile = fp1;
		} /* for */
	    } /* if ( (name != NULL) ... ) */
	} /* while */

	fclose(fp_planes);



    } /* if (myID == 0) */

    /* broadcast plane info */
    MPI_Bcast( &theNumberOfPlanes, 1, MPI_INT, 0, comm_solver );
    MPI_Bcast( thePlanePrintRate, 1, MPI_INT, 0, comm_solver );

    /* initialize the local structures */
    if (myID != 0) {
	thePlanes = (plane_t*)malloc( sizeof( plane_t ) * theNumberOfPlanes );
	if (thePlanes == NULL) {
	    solver_abort( "broadcast_planeinfo", NULL,
			  "Error: Unable to create plane information arrays" );
	}
    }

    for ( iPlane = 0; iPlane < theNumberOfPlanes; iPlane++ ) {
	MPI_Bcast( &(thePlanes[iPlane].numberofstepsalongstrike), 1, MPI_INT, 0, comm_solver);
	MPI_Bcast( &(thePlanes[iPlane].numberofstepsdowndip), 1, MPI_INT, 0, comm_solver);
	MPI_Bcast( &(thePlanes[iPlane].stepalongstrike), 1, MPI_DOUBLE,0, comm_solver);
	MPI_Bcast( &(thePlanes[iPlane].stepdowndip), 1, MPI_DOUBLE,0, comm_solver);
	MPI_Bcast( &(thePlanes[iPlane].strike), 1, MPI_DOUBLE,0, comm_solver);
	MPI_Bcast( &(thePlanes[iPlane].dip), 1, MPI_DOUBLE,0, comm_solver);
	MPI_Bcast( &(thePlanes[iPlane].rake), 1, MPI_DOUBLE,0, comm_solver);
	MPI_Bcast( &(thePlanes[iPlane].origincoords.x[0]), 3, MPI_DOUBLE,0, comm_solver);
	MPI_Barrier(comm_solver);
    }

    Old_output_planes_construct_strips(myID, theNumberOfPlanes);

    if (g_planeGpuMpiioEnabled) {
        planes_mpiio_setup(myID, theNumberOfPlanes, thePlaneDirOut);
    } else if (g_planeGpuRootEnabled) {
        planes_gpu_root_setup(myID, theNumberOfPlanes);
    }

    /* master allocates largest plane */
    if ((myID == 0) && !g_planeGpuMpiioEnabled) {
	planes_output_buffer = (double*)malloc( 3 * sizeof(double)
						* largestplanesize );
	if (planes_output_buffer == NULL) {
	    fprintf(stderr, "Error creating buffer for plane output\n");
	    exit(1);
	}
    }

    print_planeinfo(myID, theNumberOfPlanes);

    return;
}

/**
 * This builds all of the strips used in the planes output. BigBen's
 * version does not have strip size limits
 */

static void Old_output_planes_construct_strips(int32_t myID, int theNumberOfPlanes)
{

    //    int iNode, iPlane, iStrike, iDownDip;
    int iPlane, iStrike, iDownDip;
    double xLocal, yLocal;
    vector3D_t origin, pointLocal, pointGlobal;
    octant_t *octant;
    int onstrip ;
#ifdef _OPENMP
	const int32_t useOmp = plane_output_use_omp();
#endif
    //    int32_t nodesToInterpolate[8];

    planes_LocalLargestStripCount = 0;

    for ( iPlane = 0; iPlane < theNumberOfPlanes ;iPlane++ ){
    	thePlanes[ iPlane ].numberofstripsthisplane = 0;
    	/* Compute the coordinates of the nodes of the plane */
    	origin.x[0] = thePlanes[ iPlane ].origincoords.x[0];
    	origin.x[1] = thePlanes[ iPlane ].origincoords.x[1];
    	origin.x[2] = thePlanes[ iPlane ].origincoords.x[2];

//    	/*Dorian. sanity check for topo-planes "Only horizontal planes are supported" */
//    	if ( ( get_thebase_topo() != 0.0  ) &&
//    			( ( thePlanes[ iPlane ].dip > 0 ) && ( thePlanes[ iPlane ].dip < 180 ) ) ) {
//    		fprintf(stderr,
//    				"Old_output_planes_construct_strips: "
//    				"Topo-plane must be horizontal IPlane=%d, dip %f\n",iPlane, thePlanes[ iPlane ].dip);
//    		MPI_Abort(MPI_COMM_WORLD, ERROR);
//    		exit(1);
//    	}


    	/*Find limits of consecutive strips*/
    	onstrip = 0;
    	for ( iStrike = 0;
    			iStrike < thePlanes[iPlane].numberofstepsalongstrike;
    			iStrike++ )
    	{
    		xLocal = iStrike*thePlanes[iPlane].stepalongstrike;
    		for ( iDownDip = 0;
    				iDownDip < thePlanes[iPlane].numberofstepsdowndip;
    				iDownDip++ )
    		{
    			yLocal = iDownDip*thePlanes[ iPlane ].stepdowndip;
    			pointLocal.x[0] = xLocal;
    			pointLocal.x[1] = yLocal;
    			pointLocal.x[2] = 0;
    			pointGlobal
    			= compute_global_coords( origin, pointLocal,
    					thePlanes[ iPlane ].dip, 0,
    					thePlanes[ iPlane ].strike );

    			/* Dorian: correct z coordinate if topo-plane */
    			if ( get_thebase_topo() != 0.0  ) {
    				pointGlobal.x[2] +=  point_elevation ( pointGlobal.x[0], pointGlobal.x[1] );
    			}

    			if (search_point( pointGlobal, &octant ) == 1) {
    				if (!onstrip) { /* start new strip */
    					thePlanes[iPlane].stripstart[thePlanes[iPlane].numberofstripsthisplane]
    					                             = (iStrike * thePlanes[iPlane].numberofstepsdowndip
    					                            		 + iDownDip);
    					onstrip = 1;
    				}

    			}

    			else {
    				if (onstrip) { /* close strip */
    					onstrip = 0;
    					thePlanes[iPlane].stripend[thePlanes[iPlane]
    					                                     .numberofstripsthisplane]
    					                                     = (iStrike * thePlanes[iPlane].numberofstepsdowndip
    					                                    		 + iDownDip - 1);
    					thePlanes[ iPlane ].numberofstripsthisplane++;
    					if (thePlanes[ iPlane ].numberofstripsthisplane
    							> MAX_STRIPS_PER_PLANE)
    					{
    						fprintf( stderr,
    								"Number of strips on plane exceeds "
    								"MAX_STRIPS_PER_PLANE\n");
    						exit(1);
    					}
    				}
    			}
    		}
    	}

    	if (onstrip){  /*if on a strip at end of loop, close strip */
    		thePlanes[iPlane].stripend[thePlanes[iPlane].numberofstripsthisplane]
    		                           = (thePlanes[iPlane].numberofstepsdowndip
    		                        		   * thePlanes[iPlane].numberofstepsalongstrike) - 1;
    		thePlanes[ iPlane ].numberofstripsthisplane++;
    	}

    	/* get strip counts to IO PEs */
    	MPI_Reduce( &thePlanes[ iPlane ].numberofstripsthisplane,
    			&thePlanes[iPlane].globalnumberofstripsthisplane,
    			1, MPI_INT, MPI_SUM, 0, comm_solver );

    	/* allocate strips */
    	int stripnum, stripLength;

    	for (stripnum = 0;
    			stripnum < thePlanes[ iPlane ].numberofstripsthisplane;
    			stripnum++)
    	{
    		stripLength = thePlanes[iPlane].stripend[stripnum]
    		                                         - thePlanes[iPlane].stripstart[stripnum] + 1;

    		thePlanes[iPlane].strip[stripnum]
    		                        = (plane_strip_element_t*)calloc( (size_t)stripLength,
								 sizeof(plane_strip_element_t) );
    		if (thePlanes[ iPlane ].strip[stripnum] == NULL) {
    			fprintf( stderr,
    					"Error malloc'ing array strips for plane output\n"
    					"PE: %d  Plane: %d  Strip: %d  Size: %zu\n",
    					myID, iPlane, stripnum,
    					sizeof(plane_strip_element_t) * stripLength );
    			exit(1);
    		}

    		if (stripLength>planes_LocalLargestStripCount) {
    			planes_LocalLargestStripCount = stripLength;
    		}
    	}

    	/* fill strips */
    	origin.x[0] = thePlanes[ iPlane ].origincoords.x[0];
    	origin.x[1] = thePlanes[ iPlane ].origincoords.x[1];
    	origin.x[2] = thePlanes[ iPlane ].origincoords.x[2];

    	int elemnum;
    	for (stripnum = 0;
    			stripnum < thePlanes[ iPlane ].numberofstripsthisplane;
    			stripnum++)
    	{
			int elemCount = thePlanes[iPlane].stripend[stripnum]
			    - thePlanes[iPlane].stripstart[stripnum] + 1;
		#ifdef _OPENMP
		#pragma omp parallel for if(useOmp) schedule(static) \
		    private(iStrike, iDownDip, xLocal, yLocal, pointLocal, pointGlobal, octant)
		#endif
    		for (elemnum = 0;
					elemnum < elemCount;
    				elemnum++)
    		{
			thePlanes[iPlane].strip[stripnum][elemnum].topostrip = 0;
			thePlanes[iPlane].strip[stripnum][elemnum].point_is_valid = 0;
    			iStrike  = (elemnum+thePlanes[iPlane].stripstart[stripnum])
			    						/ thePlanes[iPlane].numberofstepsdowndip;
    			iDownDip = (elemnum+thePlanes[iPlane].stripstart[stripnum])
			    						% thePlanes[iPlane].numberofstepsdowndip;
    			xLocal = iStrike*thePlanes[iPlane].stepalongstrike;
    			yLocal = iDownDip*thePlanes[ iPlane ].stepdowndip;
    			pointLocal.x[0] = xLocal;
    			pointLocal.x[1] = yLocal;
    			pointLocal.x[2] = 0;
    			pointGlobal
    			=  compute_global_coords( origin, pointLocal,
    					thePlanes[ iPlane ].dip, 0,
    					thePlanes[ iPlane ].strike );

    			/* Dorian: correct z coordinate if topo-plane */
    			if ( get_thebase_topo() != 0.0 ) {
    				pointGlobal.x[2] +=  point_elevation ( pointGlobal.x[0], pointGlobal.x[1] );
    			}

    			if (search_point( pointGlobal, &octant) == 1) {
					thePlanes[iPlane].strip[stripnum][elemnum].point_is_valid = 1;
    				compute_csi_eta_dzeta( octant, pointGlobal,
							&(thePlanes[ iPlane ].strip[stripnum][elemnum].localcoords),
							thePlanes[ iPlane ].strip[stripnum][elemnum].nodestointerpolate);

					if ( ( thePlanes[ iPlane ].strip[stripnum][elemnum].nodestointerpolate[4] == 0 ) &&
						 ( thePlanes[ iPlane ].strip[stripnum][elemnum].nodestointerpolate[5] == 0 ) &&
						 ( thePlanes[ iPlane ].strip[stripnum][elemnum].nodestointerpolate[6] == 0 ) &&
						 ( thePlanes[ iPlane ].strip[stripnum][elemnum].nodestointerpolate[7] == 0 ) ) {
						thePlanes[ iPlane ].strip[stripnum][elemnum].topostrip = 1;}
					plane_init_acoustic_metadata(
						&(thePlanes[iPlane].strip[stripnum][elemnum]),
						octant);
    			}
    		}
    	}

    } /* end of plane loop: for (stripnum = 0; ...) */

    MPI_Reduce( &planes_LocalLargestStripCount,
		&planes_GlobalLargestStripCount,
		1, MPI_INT, MPI_MAX, 0, comm_solver );

    /* Legacy gather buffers are unnecessary when every rank writes its
     * owned strips directly through MPI-IO. */
    if (g_planeGpuMpiioEnabled) {
        return;
    }

    /* allocate MPI recv buffers for IO PE */
    if (myID == 0) {
    	int buf_len0 = (3 * planes_GlobalLargestStripCount) + 1;
    	planes_stripMPIRecvBuffer = (double*)malloc(sizeof(double) * buf_len0);

    	/* This is large enough for 3 double componants of biggest strip
	   plus location number */
    	if (planes_stripMPIRecvBuffer == NULL) {
    		perror("Error creating MPI master buffer strip for plane output\n");
    		exit(1);
    	}
    }

    /* allocate MPI send buffers for all PEs */
    int buf_len = (3 * planes_LocalLargestStripCount) + 1;
    planes_stripMPISendBuffer = (double*)malloc( sizeof(double) * buf_len );

    /* This is large enough for 3 double componants of biggest strip plus
       location number */
    if (planes_stripMPISendBuffer == NULL) {
	fprintf(stderr,
		"Error creating MPI send buffer strip for plane output\n");
	exit(1);
    }
}

/** Close planes files. */
void Old_planes_close(int32_t myID, int theNumberOfPlanes)
{
    int iPlane;

    if (g_planeGpuMpiioEnabled) {
        planes_mpiio_close(theNumberOfPlanes);
        return;
    }

    if (myID==0){
	for (iPlane = 0; iPlane < theNumberOfPlanes ;iPlane++) {
	    fclose( thePlanes[iPlane].fpoutputfile );
	    // NOTE: Commented out by Clifford: Since planecoords files are not implemented, the files should not be created in the first place
		// fclose( thePlanes[iPlane].fpplanecoordsfile );
	}
    }

    if (g_planeGpuRootEnabled) {
        planes_gpu_root_close();
    }
}

/** Print Plane Stats, and possibly old version plane coords file. */
void print_planeinfo(int32_t myID, int theNumberOfPlanes)
{
    int iPlane;
    int MaxStrip;

    if (myID == 0) {
	fprintf( stderr, "\nPlane Output Stats___________________________\n");

	for (iPlane = 0; iPlane < theNumberOfPlanes ;iPlane++) {
	    fprintf( stderr, "Plane dimensions for plane #%d: %d,%d\n", iPlane,
		     thePlanes[iPlane].numberofstepsalongstrike,
		     thePlanes[iPlane].numberofstepsdowndip );
	}

	for (iPlane = 0; iPlane < theNumberOfPlanes ;iPlane++) {
	    fprintf( stderr, "Total strips for plane #%d: %d\n", iPlane,
		     thePlanes[iPlane].globalnumberofstripsthisplane );
	}
    }

    for (iPlane = 0; iPlane < theNumberOfPlanes ;iPlane++) {
	MPI_Reduce( &thePlanes[ iPlane ].numberofstripsthisplane, &MaxStrip,
		    1, MPI_INT, MPI_MAX, 0, comm_solver );
	if (myID==0) {
	    fprintf(stderr, "Max strips for a PE in plane #%d: %d\n",
		    iPlane, MaxStrip);
	}
    }

    if (myID==0) {
	fprintf(stderr,
		"\nDone Plane Output Stats___________________________\n");
	fflush( stderr );
    }
}


/**************************************************************************/
/*New IO Server Version Of Code Below Here*********************************/
/**************************************************************************/


void New_planes_setup(int32_t PENum, int32_t *thePlanePrintRate, char *thePlaneDirOut,
	int io_pool_pe_count,
    int theNumberOfPlanes, const char *numericalin,
    double surfaceShift, double *theSurfaceCornersLong, double *theSurfaceCornersLat,
    double theDomainX, double theDomainY, double theDomainZ, char* planes_input_file,
	UTMZone_t* utmZone, int reconstructAcousticPlaneDisplacements) {

    static const char* fname = "new_output_planes_setup()";
    int        iPlane;
    vector3D_t originPlaneCoords;
    int largestplanesize;
    FILE*      fp;
    FILE*      fp_planes;

    g_reconstructAcousticPlaneDisplacements =
        reconstructAcousticPlaneDisplacements ? 1 : 0;

    if (PENum == 0) {

        /* obtain the general planes specifications in parameter file */
	if ( (fp = fopen ( numericalin, "r")) == NULL ) {
	    solver_abort (fname, numericalin,
			  "Error opening parameters configuration file");
	}
	if ( (parsetext(fp, "output_planes_print_rate", 'i',
			thePlanePrintRate) != 0) )
	    {
		solver_abort (fname, NULL, "Error parsing output_planes_print_rate field from %s\n",
			      numericalin);
	    }

	thePlanes = (plane_t *) malloc ( sizeof( plane_t ) * theNumberOfPlanes);
	if ( thePlanes == NULL ) {
	    solver_abort( fname, "Allocating memory for planes array",
			  "Unable to create plane information arrays" );
	}
	fclose(fp);

        /* read in the individual planes coords in planes file */
        if ( (fp_planes = fopen ( planes_input_file, "r")) == NULL ) {
            solver_abort (fname, planes_input_file,
                          "Error opening planes configuration file");
        }

	int found=0;
	while (!found) {
	    char line[LINESIZE];
	    char *name;
	    char delimiters[] = " =\n";
	    char querystring[] = "output_planes";
	    /* Read in one line */
	    if (fgets (line, LINESIZE, fp_planes) == NULL) {
		break;
	    }
	    name = strtok(line, delimiters);
	    if ( (name != NULL) && (strcmp(name, querystring) == 0) ) {
		largestplanesize = 0;
		found = 1;
		for ( iPlane = 0; iPlane < theNumberOfPlanes; iPlane++ ) {
		    int fret0, fret1;
		    fret0 = fscanf( fp_planes," %lf %lf %lf",
				    &thePlanes[iPlane].origincoords.x[0],
				    &thePlanes[iPlane].origincoords.x[1],
				    &thePlanes[iPlane].origincoords.x[2] );

		    /* RICARDO: For buildings correction */
		    thePlanes[iPlane].origincoords.x[2] += surfaceShift;

		    /* origin coordinates must be in Lat, Long and depth */
		    if (fret0 == 0) {
			solver_abort (fname, NULL,
				      "Unable to read planes origin in %s",
				      planes_input_file);
		    }
		    /* convert to cartesian refered to the mesh */
                    originPlaneCoords =
			compute_domain_coords_linearinterp(
				thePlanes[ iPlane ].origincoords.x[1], thePlanes[ iPlane ].origincoords.x[0],
				theSurfaceCornersLong, theSurfaceCornersLat,
				theDomainY, theDomainX, utmZone);

                    thePlanes[iPlane].origincoords.x[0]=originPlaneCoords.x[0];
                    thePlanes[iPlane].origincoords.x[1]=originPlaneCoords.x[1];

                    fret1 = fscanf(fp_planes," %lf %d %lf %d %lf %lf",
                                   &thePlanes[iPlane].stepalongstrike,
                                   &thePlanes[iPlane].numberofstepsalongstrike,
                                   &thePlanes[iPlane].stepdowndip,
				   &thePlanes[iPlane].numberofstepsdowndip,
                                   &thePlanes[iPlane].strike,
                                   &thePlanes[iPlane].dip);
                    /*Find largest plane for output buffer allocation */
                    if ( (thePlanes[iPlane].numberofstepsdowndip
                          * thePlanes[iPlane].numberofstepsalongstrike)
                         > largestplanesize)
			{
			    largestplanesize
				= thePlanes[iPlane].numberofstepsdowndip
				* thePlanes[iPlane].numberofstepsalongstrike;
			}
                    if (fret1 == 0) {
			solver_abort(fname, NULL,
				     "Unable to read plane specification in %s",
				     planes_input_file);
                    }

		} /* for */
	    } /* if ( (name != NULL) ... ) */
	} /* while */

	fclose(fp_planes);

    } /* if (PENum == 0) */

	/* Determine IO-server rank range in MPI_COMM_WORLD. */
	if (PENum == 0) {
		int world_size;
		int server_count;

		MPI_Comm_size(MPI_COMM_WORLD, &world_size);
		server_count = io_pool_pe_count;
		if (server_count < 1) {
			server_count = 1;
		}
		if (server_count > world_size) {
			server_count = world_size;
		}

		g_planeIOServerCount = server_count;
		g_planeIOServerFirstRank = world_size - server_count;
		g_planeIOSolverRankCount = world_size - io_pool_pe_count;
		if (g_planeIOSolverRankCount < 1) {
			g_planeIOSolverRankCount = g_planeIOServerFirstRank;
		}
	}

	/* broadcast plane info to whole IO group */
    MPI_Bcast( &theNumberOfPlanes, 1, MPI_INT, 0, comm_IO );
    MPI_Bcast( thePlanePrintRate, 1, MPI_INT, 0, comm_IO );
    MPI_Bcast( &largestplanesize, 1, MPI_INT, 0, comm_IO );
    MPI_Bcast( thePlaneDirOut, 256, MPI_CHAR, 0, comm_IO );
	MPI_Bcast( &g_planeIOServerCount, 1, MPI_INT, 0, comm_IO );
	MPI_Bcast( &g_planeIOServerFirstRank, 1, MPI_INT, 0, comm_IO );
	MPI_Bcast( &g_planeIOSolverRankCount, 1, MPI_INT, 0, comm_IO );


    /* initialize the local structures */
    if (PENum != 0) {
	thePlanes = (plane_t*)malloc( sizeof( plane_t ) * theNumberOfPlanes );
	if (thePlanes == NULL) {
	    solver_abort( "broadcast_planeinfo", NULL,
			  "Error: Unable to create plane information arrays" );
	}
    }

    /* broadcast plane params to whole IO group */
    for ( iPlane = 0; iPlane < theNumberOfPlanes; iPlane++ ) {
	MPI_Bcast( &(thePlanes[iPlane].numberofstepsalongstrike), 1, MPI_INT, 0, comm_IO);
	MPI_Bcast( &(thePlanes[iPlane].numberofstepsdowndip), 1, MPI_INT, 0, comm_IO);
	MPI_Bcast( &(thePlanes[iPlane].stepalongstrike), 1, MPI_DOUBLE,0, comm_IO);
	MPI_Bcast( &(thePlanes[iPlane].stepdowndip), 1, MPI_DOUBLE,0, comm_IO);
	MPI_Bcast( &(thePlanes[iPlane].strike), 1, MPI_DOUBLE,0, comm_IO);
	MPI_Bcast( &(thePlanes[iPlane].dip), 1, MPI_DOUBLE,0, comm_IO);
	MPI_Bcast( &(thePlanes[iPlane].rake), 1, MPI_DOUBLE,0, comm_IO);
	MPI_Bcast( &(thePlanes[iPlane].origincoords.x[0]), 3, MPI_DOUBLE,0, comm_IO);
	MPI_Barrier(comm_IO);
    }

    /* Assign IO server rank to each plane (round-robin over IO pool ranks). */
    for ( iPlane = 0; iPlane < theNumberOfPlanes; iPlane++ ) {
	thePlanes[iPlane].IO_PE_num = g_planeIOServerFirstRank +
	    (iPlane % g_planeIOServerCount);
    }

    New_output_planes_construct_strips(PENum, theNumberOfPlanes);

    print_planeinfo(PENum, theNumberOfPlanes);

    return;
}


/* This builds all of the strips used in the planes output. XT5
   version does not have strip size limits */

static void New_output_planes_construct_strips(int32_t myID, int theNumberOfPlanes)
{

    int iPlane, iStrike, iDownDip;
    double xLocal, yLocal;
    vector3D_t origin, pointLocal, pointGlobal;
    octant_t *octant;
    int onstrip ;
#ifdef _OPENMP
	const int32_t useOmp = plane_output_use_omp();
#endif

    planes_LocalLargestStripCount = 0;
	planes_LocalLargestPlanePayloadDoubles = 0;

    for ( iPlane = 0; iPlane < theNumberOfPlanes ;iPlane++ ){
	int localPlanePayloadDoubles = 0;
	thePlanes[ iPlane ].numberofstripsthisplane = 0;
	/* Compute the coordinates of the nodes of the plane */
	origin.x[0] = thePlanes[ iPlane ].origincoords.x[0];
	origin.x[1] = thePlanes[ iPlane ].origincoords.x[1];
	origin.x[2] = thePlanes[ iPlane ].origincoords.x[2];

//	/*Dorian. sanity check for topo-planes "Only horizontal planes are supported" */
//    if ( ( get_thebase_topo() != 0.0  ) &&
//    	 ( ( thePlanes[ iPlane ].dip > 0 ) && ( thePlanes[ iPlane ].dip < 180 ) ) ) {
//        fprintf(stderr,
//                "Old_output_planes_construct_strips: "
//                "Topo-plane must be horizontal IPlane=%d, dip %f\n",iPlane, thePlanes[ iPlane ].dip);
//        MPI_Abort(MPI_COMM_WORLD, ERROR);
//        exit(1);
//    }

	/*Find limits of consecutive strips*/
	onstrip = 0;
	for ( iStrike = 0;
	      iStrike < thePlanes[iPlane].numberofstepsalongstrike;
	      iStrike++ )
	    {
		xLocal = iStrike*thePlanes[iPlane].stepalongstrike;
		for ( iDownDip = 0;
		      iDownDip < thePlanes[iPlane].numberofstepsdowndip;
		      iDownDip++ )
		    {
			yLocal = iDownDip*thePlanes[ iPlane ].stepdowndip;
			pointLocal.x[0] = xLocal;
			pointLocal.x[1] = yLocal;
			pointLocal.x[2] = 0;
			pointGlobal
			    = compute_global_coords( origin, pointLocal,
						     thePlanes[ iPlane ].dip, 0,
						     thePlanes[ iPlane ].strike );

			/* Dorian: correct z coordinate if topo-plane */
			if ( get_thebase_topo() != 0.0  ) {
				pointGlobal.x[2] +=  point_elevation ( pointGlobal.x[0], pointGlobal.x[1] );
			}

			if (search_point( pointGlobal, &octant ) == 1) {
			    if (!onstrip) { /* start new strip */
				thePlanes[iPlane].stripstart[thePlanes[iPlane].numberofstripsthisplane]
				    = (iStrike * thePlanes[iPlane].numberofstepsdowndip
				       + iDownDip);
				onstrip = 1;
			    }

			}

			else {
			    if (onstrip) { /* close strip */
				onstrip = 0;
				thePlanes[iPlane].stripend[thePlanes[iPlane]
							   .numberofstripsthisplane]
				    = (iStrike * thePlanes[iPlane].numberofstepsdowndip
				       + iDownDip - 1);
				thePlanes[ iPlane ].numberofstripsthisplane++;
				if (thePlanes[ iPlane ].numberofstripsthisplane
				    > MAX_STRIPS_PER_PLANE)
				    {
					fprintf( stderr,
						 "Number of strips on plane exceeds "
						 "MAX_STRIPS_PER_PLANE\n");
					exit(1);
				    }
			    }
			}
		    }
	    }

	if (onstrip){  /*if on a strip at end of loop, close strip */
            thePlanes[iPlane].stripend[thePlanes[iPlane].numberofstripsthisplane]
		= (thePlanes[iPlane].numberofstepsdowndip
		   * thePlanes[iPlane].numberofstepsalongstrike) - 1;
            thePlanes[ iPlane ].numberofstripsthisplane++;
	}

	/* get strip counts to IO PEs */
	MPI_Reduce( &thePlanes[ iPlane ].numberofstripsthisplane,
		    &thePlanes[iPlane].globalnumberofstripsthisplane,
		    1, MPI_INT, MPI_SUM, 0, comm_solver );
	MPI_Bcast( &thePlanes[iPlane].globalnumberofstripsthisplane, 1, MPI_INT, 0, comm_IO );
	{
		int localSenderCount =
			(thePlanes[iPlane].numberofstripsthisplane > 0) ? 1 : 0;
		MPI_Reduce(&localSenderCount,
			   &thePlanes[iPlane].globalnumberofsendersthisplane,
			   1, MPI_INT, MPI_SUM, 0, comm_solver);
		MPI_Bcast(&thePlanes[iPlane].globalnumberofsendersthisplane,
			  1, MPI_INT, 0, comm_IO);
	}

	/* allocate strips */
	int stripnum, stripLength;

	for (stripnum = 0;
	     stripnum < thePlanes[ iPlane ].numberofstripsthisplane;
	     stripnum++)
	    {
		stripLength = thePlanes[iPlane].stripend[stripnum]
		    - thePlanes[iPlane].stripstart[stripnum] + 1;

    		thePlanes[iPlane].strip[stripnum]
    		                        = (plane_strip_element_t*)calloc( (size_t)stripLength,
								 sizeof(plane_strip_element_t) );
    		if (thePlanes[ iPlane ].strip[stripnum] == NULL) {
    			fprintf( stderr,
    					"Error malloc'ing array strips for plane output\n"
			     "PE: %d  Plane: %d  Strip: %d  Size: %zu\n",
			     myID, iPlane, stripnum,
			     sizeof(plane_strip_element_t) * stripLength );
		    exit(1);
		}

		if (stripLength>planes_LocalLargestStripCount) {
		    planes_LocalLargestStripCount = stripLength;
		}
		localPlanePayloadDoubles += (stripLength * 3) + 2;
	    }
	if (localPlanePayloadDoubles > planes_LocalLargestPlanePayloadDoubles) {
		planes_LocalLargestPlanePayloadDoubles = localPlanePayloadDoubles;
	}

	/* fill strips */
	origin.x[0] = thePlanes[ iPlane ].origincoords.x[0];
	origin.x[1] = thePlanes[ iPlane ].origincoords.x[1];
	origin.x[2] = thePlanes[ iPlane ].origincoords.x[2];

	int elemnum;
	for (stripnum = 0;
	     stripnum < thePlanes[ iPlane ].numberofstripsthisplane;
	     stripnum++)
	    {
		int elemCount = thePlanes[iPlane].stripend[stripnum]
		    - thePlanes[iPlane].stripstart[stripnum] + 1;
	#ifdef _OPENMP
	#pragma omp parallel for if(useOmp) schedule(static) \
	    private(iStrike, iDownDip, xLocal, yLocal, pointLocal, pointGlobal, octant)
	#endif
    		for (elemnum = 0;
		     elemnum < elemCount;
    		     elemnum++)
    		    {
			thePlanes[iPlane].strip[stripnum][elemnum].topostrip = 0;
			thePlanes[iPlane].strip[stripnum][elemnum].point_is_valid = 0;
			iStrike  = (elemnum+thePlanes[iPlane].stripstart[stripnum])
			    / thePlanes[iPlane].numberofstepsdowndip;
			iDownDip = (elemnum+thePlanes[iPlane].stripstart[stripnum])
			    % thePlanes[iPlane].numberofstepsdowndip;
			xLocal = iStrike*thePlanes[iPlane].stepalongstrike;
			yLocal = iDownDip*thePlanes[ iPlane ].stepdowndip;
			pointLocal.x[0] = xLocal;
			pointLocal.x[1] = yLocal;
			pointLocal.x[2] = 0;
			pointGlobal
			    =  compute_global_coords( origin, pointLocal,
						      thePlanes[ iPlane ].dip, 0,
						      thePlanes[ iPlane ].strike );

			/* Dorian: correct z coordinate if topo-plane */
			if ( get_thebase_topo() != 0.0 ) {
				pointGlobal.x[2] +=  point_elevation ( pointGlobal.x[0], pointGlobal.x[1] );
			}

			if (search_point( pointGlobal, &octant) == 1) {
				thePlanes[iPlane].strip[stripnum][elemnum].point_is_valid = 1;
				compute_csi_eta_dzeta( octant, pointGlobal,
						&(thePlanes[ iPlane ].strip[stripnum][elemnum].localcoords),
						thePlanes[ iPlane ].strip[stripnum][elemnum].nodestointerpolate);

				if ( ( thePlanes[ iPlane ].strip[stripnum][elemnum].nodestointerpolate[4] == 0 ) &&
					 ( thePlanes[ iPlane ].strip[stripnum][elemnum].nodestointerpolate[5] == 0 ) &&
					 ( thePlanes[ iPlane ].strip[stripnum][elemnum].nodestointerpolate[6] == 0 ) &&
					 ( thePlanes[ iPlane ].strip[stripnum][elemnum].nodestointerpolate[7] == 0 ) ) {
					thePlanes[ iPlane ].strip[stripnum][elemnum].topostrip = 1;}
				plane_init_acoustic_metadata(
					&(thePlanes[iPlane].strip[stripnum][elemnum]),
					octant);
			}
		    }
	    }
    } /* end of plane loop: for (stripnum = 0; ...) */

    MPI_Reduce( &planes_LocalLargestStripCount,
		&planes_GlobalLargestStripCount,
		1, MPI_INT, MPI_MAX, 0, comm_solver );
    MPI_Bcast( &planes_GlobalLargestStripCount, 1, MPI_INT, 0, comm_IO );
	MPI_Reduce(&planes_LocalLargestPlanePayloadDoubles,
		   &planes_GlobalLargestPlanePayloadDoubles,
		   1, MPI_INT, MPI_MAX, 0, comm_solver);
	MPI_Bcast(&planes_GlobalLargestPlanePayloadDoubles, 1, MPI_INT, 0,
		  comm_IO);

    /* allocate MPI send buffers for all PEs */
    int buf_len = (planes_LocalLargestPlanePayloadDoubles > 0)
	? planes_LocalLargestPlanePayloadDoubles : 1;
    planes_stripMPISendBuffer = (double*)malloc( sizeof(double) * buf_len );

    /* This is large enough for one packed per-plane message from this rank. */
    if (planes_stripMPISendBuffer == NULL) {
	fprintf(stderr,
		"Error creating MPI send buffer strip for plane output\n");
	exit(1);
    }
}

/**
 * Interpolate the displacenents and communicate strips to printing PE.
 * Bigben version has no control flow, just chugs through planes.
 */
int New_planes_print(int32_t PENum, mysolver_t* mySolver, int theNumberOfPlanes)
{
	int iPlane;
    int iStrip, elemnum;
#ifdef _OPENMP
    const int32_t useOmp = plane_output_use_omp();
#endif

    for (iPlane = 0; iPlane < theNumberOfPlanes ;iPlane++) {
	int planeOffset = 0;
	double* planeBuffer = planes_stripMPISendBuffer;
	for (iStrip = 0; iStrip < thePlanes[iPlane].numberofstripsthisplane;
	     iStrip++) {
	    int elemCount = thePlanes[iPlane].stripend[iStrip]
	        - thePlanes[iPlane].stripstart[iStrip] + 1;
	    int stripStart = thePlanes[iPlane].stripstart[iStrip];
	    int stripLength = elemCount;
	    double* stripBuffer = planeBuffer + planeOffset + 2;

	    /* interpolate points directly into send buffer for this strip */
	#ifdef _OPENMP
	#pragma omp parallel for if(useOmp) schedule(static)
	#endif
	    for (elemnum = 0; elemnum < elemCount; elemnum++){
	    	double displacementsX;
	    	double displacementsY;
	    	double displacementsZ;

	    	plane_interpolate_point(mySolver,
	    			&(thePlanes[iPlane].strip[iStrip][elemnum]),
	    			&displacementsX,
	    			&displacementsY,
	    			&displacementsZ);

	    	stripBuffer[elemnum*3] = displacementsX;
	    	stripBuffer[elemnum*3 + 1] = displacementsY;
	    	stripBuffer[elemnum*3 + 2] = displacementsZ;

	    } /*for elemnum*/

	    planeBuffer[planeOffset] = (double)stripStart;
	    planeBuffer[planeOffset + 1] = (double)stripLength;
	    planeOffset += (stripLength * 3) + 2;

	} /*for iStrips*/

	if (planeOffset > 0) {
	    int planeTag = g_planePrintSequence * theNumberOfPlanes + iPlane;
	    MPI_Send(planeBuffer, planeOffset, MPI_DOUBLE,
		     thePlanes[iPlane].IO_PE_num, planeTag, comm_IO);
	}


	MPI_Barrier( comm_solver ); /* may not be required */


    } /* for iPlane */

    g_planePrintSequence++;

    return 1;
}


/** Close planes files. */
void New_planes_close(int32_t PEnum, int theNumberOfPlanes){

    int iServer;

    /* Every solver rank sends STOP after its own plane data. The IO server
     * waits for all STOPs, which uses MPI's per-source message ordering to
     * guarantee that all earlier plane messages have been received first. */
    if (PEnum < g_planeIOSolverRankCount){
	planes_stripMPISendBuffer[0]= (STOP_SERVER + 0.1);
	for (iServer = 0; iServer < g_planeIOServerCount; iServer++) {
	    int serverRank = g_planeIOServerFirstRank + iServer;
	    MPI_Send( planes_stripMPISendBuffer, 1, MPI_DOUBLE,
	              serverRank, 0, comm_IO );
	}
    }
}



/* Server IO PEs Main Loop */
void planes_IO_PES_main(int32_t PEnum){

    int        largestplanesize, comm_size, iPlane;
    MPI_Status status;
    char       thePlaneDirOut[256], planedisplacementsout[1024];
    FILE*      fp0;
    int32_t    thePlanePrintRate;
    int        theNumberOfPlanes;

    /* Get basic planes parameters */
    MPI_Bcast( &theNumberOfPlanes, 1, MPI_INT, 0, comm_IO );
    MPI_Bcast( &thePlanePrintRate, 1, MPI_INT, 0, comm_IO );
    MPI_Bcast( &largestplanesize, 1, MPI_INT, 0, comm_IO );
    MPI_Bcast( thePlaneDirOut, 256, MPI_CHAR, 0, comm_IO );
	MPI_Bcast( &g_planeIOServerCount, 1, MPI_INT, 0, comm_IO );
	MPI_Bcast( &g_planeIOServerFirstRank, 1, MPI_INT, 0, comm_IO );
	MPI_Bcast( &g_planeIOSolverRankCount, 1, MPI_INT, 0, comm_IO );

    if(theNumberOfPlanes==0) return;

    /* Create Planes info structure */
    thePlanes = (plane_t*)malloc( sizeof( plane_t ) * theNumberOfPlanes );
    if (thePlanes == NULL) {
	solver_abort( "broadcast_planeinfo", NULL,
		      "Error: Unable to create plane information arrays" );
    }

    /* Get plane params */
    for ( iPlane = 0; iPlane < theNumberOfPlanes; iPlane++ ) {
	MPI_Bcast( &(thePlanes[iPlane].numberofstepsalongstrike), 1, MPI_INT, 0, comm_IO);
	MPI_Bcast( &(thePlanes[iPlane].numberofstepsdowndip), 1, MPI_INT, 0, comm_IO);
	MPI_Bcast( &(thePlanes[iPlane].stepalongstrike), 1, MPI_DOUBLE,0, comm_IO);
	MPI_Bcast( &(thePlanes[iPlane].stepdowndip), 1, MPI_DOUBLE,0, comm_IO);
	MPI_Bcast( &(thePlanes[iPlane].strike), 1, MPI_DOUBLE,0, comm_IO);
	MPI_Bcast( &(thePlanes[iPlane].dip), 1, MPI_DOUBLE,0, comm_IO);
	MPI_Bcast( &(thePlanes[iPlane].rake), 1, MPI_DOUBLE,0, comm_IO);
	MPI_Bcast( &(thePlanes[iPlane].origincoords.x[0]), 3, MPI_DOUBLE,0, comm_IO);
	MPI_Barrier(comm_IO);
    }

    /* Get strip count for each plane */
    for (iPlane = 0; iPlane < theNumberOfPlanes; iPlane++)
	MPI_Bcast( &thePlanes[iPlane].globalnumberofstripsthisplane, 1, MPI_INT, 0, comm_IO );
	for (iPlane = 0; iPlane < theNumberOfPlanes; iPlane++)
	MPI_Bcast( &thePlanes[iPlane].globalnumberofsendersthisplane, 1, MPI_INT, 0, comm_IO );

    /* Find size of largest strip */
    MPI_Bcast( &planes_GlobalLargestStripCount, 1, MPI_INT, 0, comm_IO );
	MPI_Bcast( &planes_GlobalLargestPlanePayloadDoubles, 1, MPI_INT, 0, comm_IO );

	/* Idle IO pool PEs can leave, don't use collectives past here! */
	MPI_Comm_size(comm_IO, &comm_size);
	if ((PEnum < g_planeIOServerFirstRank) ||
		(PEnum >= (g_planeIOServerFirstRank + g_planeIOServerCount))) {
		return;
	}

    /* Allocate largest required strip recieve buffer */
	int buf_len0 = (planes_GlobalLargestPlanePayloadDoubles > 0)
	? planes_GlobalLargestPlanePayloadDoubles : 1;
    planes_stripMPIRecvBuffer = (double*)malloc(sizeof(double) * buf_len0);
	/* large enough for one packed per-plane message from any rank */
    if (planes_stripMPIRecvBuffer == NULL) {
	perror("Error creating MPI master buffer strip for plane output\n");
	exit(1);
    }

    /* Allocate largest plane output buffer */
    planes_output_buffer = (double*)malloc( 3 * sizeof(double)
					    * largestplanesize );
    if (planes_output_buffer == NULL) {
	fprintf(stderr, "Error creating buffer for plane output\n");
	exit(1);
    }

    /* Open displacement files only for planes owned by this IO server rank. */
    for (iPlane = 0; iPlane < theNumberOfPlanes ;iPlane++){
	if (iPlane % g_planeIOServerCount !=
	    (PEnum - g_planeIOServerFirstRank)) {
	    thePlanes[iPlane].fpoutputfile = NULL;
	    continue;
	}

	sprintf( planedisplacementsout, "%s/planedisplacements.%d",
		 thePlaneDirOut, iPlane );
	fp0 = fopen( planedisplacementsout, "w" );
	thePlanes[iPlane].fpoutputfile = fp0;
	/* Not activating planecoored files yet
	   sprintf (planecoordsfile, "%s/planecoords.%d",
	   thePlaneDirOut, iPlane);
	   fp1 = fopen (planecoordsfile, "w");
	   thePlanes[iPlane].fpplanecoordsfile = fp1;
	   Turn this on when and if it is needed */
    }

    int stopCount = 0;
    while(1){

	int recvStripCount, rStrip;

	/* Get any message from anywhere */
	MPI_Recv( planes_stripMPIRecvBuffer, buf_len0, MPI_DOUBLE,
		  MPI_ANY_SOURCE, MPI_ANY_TAG, comm_IO, &status );
	MPI_Get_count(&status, MPI_DOUBLE, &recvStripCount);

	/* Check for a stop message */
	if( recvStripCount == 1){
	    if( (int)planes_stripMPIRecvBuffer[0] == STOP_SERVER) {
		stopCount++;
		if (stopCount >= g_planeIOSolverRankCount)
		    break;
		continue;
	    }
	    solver_abort( "IO_Server_Main_Loop", NULL,
	                  "Error: malformed one-value plane IO message" );
	}
	/* else pick this plane and fill in first strip */
	else{
	    int planeTag = status.MPI_TAG;
	    iPlane = planeTag % theNumberOfPlanes;
	    {
		int planeValues = 3 * thePlanes[iPlane].numberofstepsalongstrike
		    * thePlanes[iPlane].numberofstepsdowndip;
		memset(planes_output_buffer, 0, (size_t)planeValues * sizeof(double));
	    }
	    plane_unpack_message(planes_stripMPIRecvBuffer, recvStripCount);

	    /* Fill in rest of this plane for the same print sequence tag. */
	    for (rStrip = 1; rStrip < thePlanes[iPlane].globalnumberofsendersthisplane; rStrip++) {
		MPI_Recv( planes_stripMPIRecvBuffer,
			  buf_len0, MPI_DOUBLE,
			  MPI_ANY_SOURCE, planeTag, comm_IO, &status );
		MPI_Get_count(&status, MPI_DOUBLE, &recvStripCount);
		plane_unpack_message(planes_stripMPIRecvBuffer, recvStripCount);
	    }
	}

	/* Print just this one plane */
	int ret, items_to_write;
	items_to_write = 3 * thePlanes[iPlane].numberofstepsalongstrike * thePlanes[iPlane].numberofstepsdowndip;
	if (thePlanes[iPlane].fpoutputfile == NULL) {
	    solver_abort( "IO_Server_Main_Loop", NULL,
	                 "Error: server rank %d received unassigned plane %d",
	                 PEnum, iPlane );
	}
	ret = fwrite( planes_output_buffer, sizeof(double), items_to_write, thePlanes[iPlane].fpoutputfile );
	if (ret != items_to_write)
	    solver_abort( "IO_Server_Main_Loop", NULL, "Error: Unable to output displacement file to disk" );

    }/*while*/

    /* Close all planes.  Will need to make sure these are seperated on multiple PEs*/
    // int ret;
    for (iPlane = 0; iPlane < theNumberOfPlanes ;iPlane++){
		if (thePlanes[iPlane].fpoutputfile != NULL) {
			fclose( thePlanes[iPlane].fpoutputfile );
		}
		// ret = fclose( thePlanes[iPlane].fpoutputfile );
		// if (!ret) printf("Error closing plane file %d!!!\n", iPlane);  /* Sort this out later */
    }

    return;

}
