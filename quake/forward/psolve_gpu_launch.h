/* -*- C++ -*- */

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
 * @copyright_notice_end */

#ifndef PSOLVE_GPU_LAUNCH_H_
#define PSOLVE_GPU_LAUNCH_H_

#include <stddef.h>
#include <stdint.h>

#include <cuda_runtime.h>

#include "psolve.h"
#include "quake_util.h"
#include "topography.h"

/* Device-only launch details live behind this host-callable interface.  This
 * keeps the large psolve translation unit free of CUDA launch syntax so the
 * MPI C++ compiler can build it in parallel with this much smaller NVCC unit. */
typedef struct topo_bkt_const_t {
  double tetraVol[5];
  double edgesize;
  double mu_bkt;
  double kappa_bkt;
  int cube_part;
} topo_bkt_const_t;

typedef enum psolve_gpu_kernel_id_t {
  PSOLVE_GPU_KERNEL_STIFFNESS = 0,
  PSOLVE_GPU_KERNEL_DAMPING_BKT2,
  PSOLVE_GPU_KERNEL_DAMPING_BKT3,
  PSOLVE_GPU_KERNEL_DAMPING_GENERIC,
  PSOLVE_GPU_KERNEL_DAMPING_FORCE,
  PSOLVE_GPU_KERNEL_DISPLACEMENT,
  PSOLVE_GPU_KERNEL_STATION_TENSOR
} psolve_gpu_kernel_id_t;

char* psolve_gpu_kernel_handle(psolve_gpu_kernel_id_t kernel);
int32_t psolve_gpu_kernel_blocksize(gpu_spec_t* gpuSpecs,
                                   psolve_gpu_kernel_id_t kernel,
                                   int32_t memPerThread);

void psolve_launch_station_tensor(
    int grid, int block, size_t shared, cudaStream_t stream,
    int32_t stationCount, int32_t outputStride, int32_t strainOffset,
    int32_t stressOffset, gpu_data_t* gpuData);
void psolve_launch_source_force_scatter(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t count,
    const int32_t* lnidList, const vector3D_t* srcForces,
    fvector_t* forceDevice, double dt2);
void psolve_launch_sparse_force_scatter(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t count,
    const int32_t* lnidList, const fvector_t* srcForces,
    fvector_t* forceDevice);
void psolve_launch_disp_boundary_gather(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t count,
    const int32_t* lnidList, const fvector_t* dispDevice,
    fvector_t* bufferDevice);
void psolve_launch_disp_boundary_scatter(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t count,
    const int32_t* lnidList, const fvector_t* bufferDevice,
    fvector_t* dispDevice);
void psolve_launch_scalar_boundary_gather(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t count,
    const int32_t* lnidList, const solver_float* srcDevice,
    solver_float* bufferDevice);
void psolve_launch_scalar_boundary_scatter(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t count,
    const int32_t* lnidList, const solver_float* bufferDevice,
    solver_float* dstDevice);
void psolve_launch_scalar_boundary_scatter_serial(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t count,
    const int32_t* lnidList, const solver_float* bufferDevice,
    solver_float* dstDevice);
void psolve_launch_scalar_boundary_accumulate(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t count,
    const int32_t* lnidList, const solver_float* bufferDevice,
    solver_float* dstDevice);
void psolve_launch_scalar_adjust_distribution(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t runCount,
    const int32_t* runOffsets, const int32_t* runParentBase3,
    const int32_t* edgeDn, const int32_t* danglingBase3,
    const solver_float* invDeps, solver_float* values);
void psolve_launch_scalar_adjust_assignment(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t dnCount,
    const int32_t* dnOrder, const int32_t* parentOffsets,
    const int32_t* parentBase3, const int32_t* danglingBase3,
    const solver_float* invDeps, solver_float* values);
void psolve_launch_vector_boundary_gather(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t count,
    const int32_t* lnidList, const fvector_t* srcDevice,
    fvector_t* bufferDevice);
void psolve_launch_vector_boundary_scatter(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t count,
    const int32_t* lnidList, const fvector_t* bufferDevice,
    fvector_t* dstDevice);
void psolve_launch_vector_boundary_scatter_serial(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t count,
    const int32_t* lnidList, const fvector_t* bufferDevice,
    fvector_t* dstDevice);
void psolve_launch_vector_boundary_accumulate(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t count,
    const int32_t* lnidList, const fvector_t* bufferDevice,
    fvector_t* dstDevice);
void psolve_launch_vector_adjust_distribution(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t runCount,
    const int32_t* runOffsets, const int32_t* runParentBase3,
    const int32_t* edgeDn, const int32_t* danglingBase3,
    const solver_float* invDeps, fvector_t* values);
void psolve_launch_vector_adjust_assignment(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t dnCount,
    const int32_t* dnOrder, const int32_t* parentOffsets,
    const int32_t* parentBase3, const int32_t* danglingBase3,
    const solver_float* invDeps, fvector_t* values);

void psolve_launch_stiffness(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t lenum,
    gpu_data_t* gpuData, int32_t linearElementCount,
    int32_t* linearElementMapper);
void psolve_launch_pml_elastic(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t lenum,
    gpu_data_t* gpuData, int32_t pmlElasticCount);
void psolve_launch_damping_bkt2(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t lenum,
    int32_t numcomp, gpu_data_t* gpuData, double rmax,
    double shearCoef1, double shearCoef2, double shearExpCoef0,
    double shearCoef3, double shearCoef4, double shearExpCoef1,
    double kappaCoef1, double kappaCoef2, double kappaExpCoef0,
    double kappaCoef3, double kappaCoef4, double kappaExpCoef1,
    int32_t elementCount, int32_t* elementMapper);
void psolve_launch_damping_bkt3(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t lenum,
    int32_t numcomp, gpu_data_t* gpuData, double rmax,
    double shearCoef1, double shearCoef2, double shearExpCoef0,
    double shearCoef3, double shearCoef4, double shearExpCoef1,
    double shearCoef5, double shearCoef6, double shearExpCoef2,
    double kappaCoef1, double kappaCoef2, double kappaExpCoef0,
    double kappaCoef3, double kappaCoef4, double kappaExpCoef1,
    double kappaCoef5, double kappaCoef6, double kappaExpCoef2,
    int32_t elementCount, int32_t* elementMapper);
void psolve_launch_damping_generic(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t lenum,
    int32_t numcomp, gpu_data_t* gpuData, double rmax, int32_t dampingType,
    int32_t elementCount, int32_t* elementMapper);
void psolve_launch_damping_force(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t lenum,
    gpu_data_t* gpuData, int32_t elementCount, int32_t* elementMapper);
void psolve_launch_topography_bkt(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t lenum,
    gpu_data_t* gpuData, int32_t elementCount, const int32_t* elementMapper,
    const topo_bkt_const_t* constants, double dt2);
void psolve_launch_topography_effective(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t lenum,
    gpu_data_t* gpuData, int32_t elementCount, const int32_t* elementMapper,
    const topoconstants_t* constants, double dt2);

void psolve_launch_acoustic_displacement_coupling(
    int grid, int block, size_t shared, cudaStream_t stream,
    int32_t faceCount, gpu_data_t* gpuData, double dt2, double dt);
void psolve_launch_acoustic_pressure_force(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t lenum,
    gpu_data_t* gpuData, double dt2);
void psolve_launch_pml_acoustic(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t lenum,
    gpu_data_t* gpuData, int32_t pmlAcousticCount);
void psolve_launch_pml_acoustic_gravity_surface(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t lenum,
    gpu_data_t* gpuData, int32_t surfaceCount);
void psolve_launch_acoustic_interface_force(
    int grid, int block, size_t shared, cudaStream_t stream,
    int32_t faceCount, gpu_data_t* gpuData, double dt2, double dt);
void psolve_launch_acoustic_pressure_update(
    int grid, int block, size_t shared, cudaStream_t stream,
    int32_t nharbored, gpu_data_t* gpuData);
void psolve_launch_poroelastic_pressure_traction(
    int grid, int block, size_t shared, cudaStream_t stream,
    int32_t elementCount, gpu_data_t* gpuData, double dt2,
    double biotAlpha);
void psolve_launch_acoustic_open_pore_target(
    int grid, int block, size_t shared, cudaStream_t stream,
    int32_t faceCount, gpu_data_t* gpuData, double dt);
void psolve_launch_poroelastic_pressure_force(
    int grid, int block, size_t shared, cudaStream_t stream,
    int32_t elementCount, gpu_data_t* gpuData, double dt, double mobilityX,
    double mobilityY, double mobilityZ, int haveDiffusion, double biotAlpha);
void psolve_launch_poroelastic_pressure_update(
    int grid, int block, size_t shared, cudaStream_t stream,
    int32_t nharbored, gpu_data_t* gpuData, int useDrainMass,
    int useOpenPore);

void psolve_launch_plane_interpolate(
    int grid, int block, size_t shared, cudaStream_t stream,
    int32_t pointCount, const int32_t* pointNodes, const double* pointCoords,
    const int32_t* pointModes, const fvector_t* tm1, double* pointOutput);
void psolve_launch_acoustic_station_interpolate(
    int grid, int block, size_t shared, cudaStream_t stream,
    int32_t stationCount, const int32_t* stationNodes,
    const double* stationCoords, const double* stationH,
    const solver_float* pressureTm1, const solver_float* pressureTm2,
    int32_t velocityPotential, double inverseDt, double* stationOutput);
void psolve_launch_pore_pressure_station_interpolate(
    int grid, int block, size_t shared, cudaStream_t stream,
    int32_t stationCount, const int32_t* stationNodes,
    const double* stationCoords, const solver_float* porePressureTm1,
    const solver_float* porePressureMass, double* stationOutput);
void psolve_launch_displacement(
    int grid, int block, size_t shared, cudaStream_t stream,
    int32_t nharbored, gpu_data_t* gpuData, noyesflag_t printAccel);

#endif /* PSOLVE_GPU_LAUNCH_H_ */
