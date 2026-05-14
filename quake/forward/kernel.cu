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

#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "psolve.h"
#include "kernel.h"
#include "quake_util.h"
#include "damping.h"
#include <cuda.h>
#include <cuda_runtime.h>


/**
 * For regular simulations use 1e-20.  For performance and scalability,
 * uncomment the immediate return in vector_is_zero() and vector_is_all_zero().
 * Alternatives with each platform underflow precision limit for 'double'
 * may also work, though this have not been thoroughly tested.
 *
 * Known underflow limits:
 * - NICS' Kraken: Limit is ~2.2e-308 --> use 1e-200
 *
 */
#define UNDERFLOW_CAP_STIFFNESS 1e-20


/* Constant memory symbols */
__constant__ fmatrix_t theK1Device[8][8];
__constant__ fmatrix_t theK2Device[8][8];

__host__ __device__ static inline uint8_t acoustic_node_axis_mask(int node,
                                                                  uint8_t flags)
{
  uint8_t mask = 0;

  if ((flags & ACOUSTIC_FACE_X_MINUS) &&
      ((node == 0) || (node == 2) || (node == 4) || (node == 6))) {
    mask |= 0x01u;
  }
  if ((flags & ACOUSTIC_FACE_X_PLUS) &&
      ((node == 1) || (node == 3) || (node == 5) || (node == 7))) {
    mask |= 0x01u;
  }
  if ((flags & ACOUSTIC_FACE_Y_MINUS) &&
      ((node == 0) || (node == 1) || (node == 4) || (node == 5))) {
    mask |= 0x02u;
  }
  if ((flags & ACOUSTIC_FACE_Y_PLUS) &&
      ((node == 2) || (node == 3) || (node == 6) || (node == 7))) {
    mask |= 0x02u;
  }
  if ((flags & ACOUSTIC_FACE_Z_MINUS) &&
      ((node == 0) || (node == 1) || (node == 2) || (node == 3))) {
    mask |= 0x04u;
  }
  if ((flags & ACOUSTIC_FACE_Z_PLUS) &&
      ((node == 4) || (node == 5) || (node == 6) || (node == 7))) {
    mask |= 0x04u;
  }

  return mask;
}

__host__ __device__ static inline void acoustic_project_interface_nodes(
    fvector_t localVec[8], uint8_t flags)
{
  for (int node = 0; node < 8; node++) {
    uint8_t mask = acoustic_node_axis_mask(node, flags);

    if (mask == 0) {
      continue;
    }

    if ((mask & 0x01u) == 0) {
      localVec[node].f[0] = 0.0;
    }
    if ((mask & 0x02u) == 0) {
      localVec[node].f[1] = 0.0;
    }
    if ((mask & 0x04u) == 0) {
      localVec[node].f[2] = 0.0;
    }
  }
}

__host__ __device__ static inline int acoustic_face_has_node(uint8_t face,
                                                             int node)
{
  switch (face) {
    case ACOUSTIC_FACE_X_MINUS:
      return (node == 0) || (node == 2) || (node == 4) || (node == 6);
    case ACOUSTIC_FACE_X_PLUS:
      return (node == 1) || (node == 3) || (node == 5) || (node == 7);
    case ACOUSTIC_FACE_Y_MINUS:
      return (node == 0) || (node == 1) || (node == 4) || (node == 5);
    case ACOUSTIC_FACE_Y_PLUS:
      return (node == 2) || (node == 3) || (node == 6) || (node == 7);
    case ACOUSTIC_FACE_Z_MINUS:
      return (node == 0) || (node == 1) || (node == 2) || (node == 3);
    case ACOUSTIC_FACE_Z_PLUS:
      return (node == 4) || (node == 5) || (node == 6) || (node == 7);
    default:
      return 0;
  }
}

__host__ __device__ static inline double acoustic_scalar_stiffness_shape(int i,
                                                                         int j)
{
  const double sx[8] = {-1, 1, -1, 1, -1, 1, -1, 1};
  const double sy[8] = {-1, -1, 1, 1, -1, -1, 1, 1};
  const double sz[8] = {-1, -1, -1, -1, 1, 1, 1, 1};

  const double x = sx[i] * sx[j] * (1.0 + sy[i] * sy[j] / 3.0) *
                   (1.0 + sz[i] * sz[j] / 3.0);
  const double y = sy[i] * sy[j] * (1.0 + sx[i] * sx[j] / 3.0) *
                   (1.0 + sz[i] * sz[j] / 3.0);
  const double z = sz[i] * sz[j] * (1.0 + sx[i] * sx[j] / 3.0) *
                   (1.0 + sy[i] * sy[j] / 3.0);

  return x + y + z;
}

__host__ __device__ static inline double acoustic_shape_value_at(int node,
                                                                 double csi,
                                                                 double eta,
                                                                 double dzeta)
{
  const double sx[8] = {-1, 1, -1, 1, -1, 1, -1, 1};
  const double sy[8] = {-1, -1, 1, 1, -1, -1, 1, 1};
  const double sz[8] = {-1, -1, -1, -1, 1, 1, 1, 1};

  return 0.125 * (1.0 + sx[node] * csi) * (1.0 + sy[node] * eta) *
         (1.0 + sz[node] * dzeta);
}


/* Flag denoting if kernel FLOP count cache initialized */
static int kernel_cached_flag = 0;

/* Operation type array indices */
#define FLOP_MULT_INDEX  0
#define FLOP_ADD_INDEX   1
#define FLOP_TRANS_INDEX 2
#define FLOP_MEM_INDEX   3

/* Operations counts for each kernel by thread (MULT/DIV, ADD/SUB, TRANS, MEM).
   These values were manually tabulated by code review. */
static int64_t kernel_ops[4][4] = {52, 373, 0, 416,
				   63,  54, 4, 716,
				   85, 547, 0, 592,
                   6,   9, 0, 128};

/* Operation coefficents (MULT/DIV, ADD/SUB, TRANS, MEM). */
static int64_t kernel_coef[4] = {1, 1, 1, 1};

/* Cached FLOP/memory storage */
static int64_t kernel_ops_cached[4];
static int64_t kernel_mem_cached[4];


/* Initialize the kernel FLOP count cache */
int kernel_init_cache()
{
  int i, j;

  for (i = 0; i < FLOP_MAX_KERNEL; i++) {
    kernel_ops_cached[i] = 0;
    for (j = 0; j < 3; j++) {
      kernel_ops_cached[i] += kernel_coef[j] * kernel_ops[i][j];
    }
    kernel_mem_cached[i] = kernel_ops[i][3];
  }

  return (0);
}


/* Return FLOP count for the specified kernel */
int64_t kernel_flops_per_thread(int kid)
{
  if ((kid < 0) || (kid >= FLOP_MAX_KERNEL)) {
    return (0);
  }

  /* Initialize the kernel FLOP count cache */
  if (kernel_cached_flag == 0) {
    kernel_init_cache();
    kernel_cached_flag = 1;
  }

  return(kernel_ops_cached[kid]);
}


/* Return memory count for the specified kernel */
int64_t kernel_mem_per_thread(int kid)
{
  if ((kid < 0) || (kid >= FLOP_MAX_KERNEL)) {
    return (0);
  }

  /* Initialize the kernel FLOP count cache */
  if (kernel_cached_flag == 0) {
    kernel_init_cache();
    kernel_cached_flag = 1;
  }

  return(kernel_mem_cached[kid]);
}



/* Get the register count for the specified kernel function handle */
int32_t gpu_get_reg_count(char* kernel)
{
    /* Get kernel attributes */
    cudaFuncAttributes attributes;
    cudaFuncGetAttributes(&attributes, kernel);

    return(attributes.numRegs);
}


/* Get the recommended block size for a kernel given the GPU specifications
   and shared memory per thread */
int32_t gpu_get_blocksize(gpu_spec_t *gpuSpecs,
			  char* kernel,
			  int32_t memPerThread)
{
    int32_t regs_per_warp;
    int32_t max_threads;
    int32_t computed;

    /* Get kernel attributes */
    cudaFuncAttributes attributes;
    cudaFuncGetAttributes(&attributes, kernel);

    /* Maximum theads possible based on register use */
    if (attributes.numRegs > 0) {
      //computed = gpuSpecs->regs_per_block / attributes.numRegs;
      //computed = 1 << (int)floor(log(computed)/log(2));
      regs_per_warp = (int)ceil((attributes.numRegs * gpuSpecs->warp_size)
				/ (float)gpuSpecs->register_allocation_size)
	* gpuSpecs->register_allocation_size;
      computed = (int)floor(gpuSpecs->regs_per_block
			    / (float)(regs_per_warp *
				      gpuSpecs->warp_allocation_size)) *
	gpuSpecs->warp_allocation_size * gpuSpecs->warp_size;
      max_threads = imin(computed, gpuSpecs->max_threads);
    } else {
      max_threads = gpuSpecs->max_threads;
    }

    /* Maximum threads possible based on shared memory use */
    if (memPerThread > 0) {
      computed = gpuSpecs->shared_per_block / memPerThread;
      computed = 1 << (int)floor(log(computed)/log(2));
      max_threads = imin(computed, max_threads);
    }

    return(max_threads);
}


int gpu_copy_constant_symbols(int myID, fmatrix_t (*theK1)[8],
			      fmatrix_t (*theK2)[8])
{
  cudaError_t cerror;

  cudaMemcpyToSymbol(theK1Device, &(theK1[0][0]),
		     8 * 8 * sizeof(fmatrix_t));
  cudaMemcpyToSymbol(theK2Device, &(theK2[0][0]), 
		     8 * 8 * sizeof(fmatrix_t));

  cerror = cudaGetLastError();
  if (cerror != cudaSuccess) {
    fprintf(stderr, "Thread %d: Failed to copy symbols - %s\n", myID,
	    cudaGetErrorString(cerror));
    return(1);
  }

  return(0);
}


/* Stiffness Calc-Force Kernel
   FLOPs: 4M + 25A + ((147A) + (48M + 33A) + (168A)) = 52M + 373A
   MEM  : 4 + 4 + 8*(24 + 1) + 2*8 + 8 * (24) = 416
*/
__global__ 
//__launch_bounds__(768, 1)
void kernelStiffnessCalcLocal(int32_t lenum,
			      gpu_data_t *gpuData,
			      int32_t   myLinearElementsCount,
			      int32_t*  myLinearElementsMapperDevice)
{
    int       i;
    int32_t   eindex;
    int32_t   lin_eindex = (blockIdx.x * blockDim.x) + threadIdx.x;
    fvector_t localForce[8];
    int       do_contrib = 0;

    /* Extra threads and threads for non-linear elements exit here */
    if (lin_eindex >= myLinearElementsCount) {
      return;
    }

    eindex = myLinearElementsMapperDevice[lin_eindex];

    /* CPML: PML elements are handled by kernelPMLElasticStep, so skip them
     * here when abc_type=cpml. pmlFlagArrayDevice is NULL in dashpot mode. */
    if (gpuData->pmlFlagArrayDevice != NULL &&
        gpuData->pmlFlagArrayDevice[eindex] != 0) {
      return;
    }

    int32_t *lnid = &(gpuData->lnidArrayDevice[eindex]);

    for (i = 0; i < 8; i++) {
      const __restrict__ fvector_t *tm1Disp = gpuData->tm1Device + *lnid;

      localForce[i].f[0] = tm1Disp->f[0];
      do_contrib = ( fabs(localForce[i].f[0]) > UNDERFLOW_CAP_STIFFNESS ) ? 1 : do_contrib;
      localForce[i].f[1] = tm1Disp->f[1];
      do_contrib = ( fabs(localForce[i].f[1]) > UNDERFLOW_CAP_STIFFNESS ) ? 1 : do_contrib;
      localForce[i].f[2] = tm1Disp->f[2];
      do_contrib = ( fabs(localForce[i].f[2]) > UNDERFLOW_CAP_STIFFNESS ) ? 1 : do_contrib;

      lnid += lenum;
    }

    if (!do_contrib) {
      return;
    }

    uint8_t acousticFaceFlags = gpuData->acousticFaceFlagsArrayDevice[eindex];
    if (acousticFaceFlags != 0) {
      acoustic_project_interface_nodes(localForce, acousticFaceFlags);
    }

    /* Coefficients for new stiffness matrix calculation */
    //if (vector_is_zero( localForce ) != 0) {
    double first_coeff  = -0.5625 * (gpuData->c2ArrayDevice[eindex] +
				     2 * gpuData->c1ArrayDevice[eindex]);
    double second_coeff = -0.5625 * (gpuData->c2ArrayDevice[eindex]);
    double third_coeff  = -0.5625 * (gpuData->c1ArrayDevice[eindex]);

    double atu[24];
    double firstVec[24];

    aTransposeU( localForce, atu );
    firstVector( atu, firstVec, first_coeff, second_coeff, third_coeff );
    au( localForce, firstVec );
    //} else {
    //  /* Returning saves us from memset of localForce buffer */
    //  return;
    //}

    if (acousticFaceFlags != 0) {
      acoustic_project_interface_nodes(localForce, acousticFaceFlags);
    }

    /* Sum the nodal forces */
    lnid = &(gpuData->lnidArrayDevice[eindex]);
    for (i = 0; i < 8; i++) {
      fvector_t *nodalForce = gpuData->forceDevice + *lnid;

      atomicAdd((solver_float *)&(nodalForce->f[0]),
		(solver_float)localForce[i].f[0]);
      atomicAdd((solver_float *)&(nodalForce->f[1]),
		(solver_float)localForce[i].f[1]);
      atomicAdd((solver_float *)&(nodalForce->f[2]),
		(solver_float)localForce[i].f[2]);

      lnid += lenum;
    }

    return;
}

__global__ void kernelAcousticPressureDisplacementCoupling(int32_t faceCount,
                                                           gpu_data_t* gpuData,
                                                           double deltaT2)
{
  int32_t faceIndex = (blockIdx.x * blockDim.x) + threadIdx.x;

  if (faceIndex >= faceCount) {
    return;
  }

  const int32_t eindex =
      gpuData->acousticInterfaceElementIdsDevice[faceIndex];
  const double face_coeff =
      deltaT2 * gpuData->acousticInterfaceFaceAreasDevice[faceIndex];
  const double center[3] = {
      gpuData->acousticInterfaceFaceCentersDevice[3 * faceIndex + 0],
      gpuData->acousticInterfaceFaceCentersDevice[3 * faceIndex + 1],
      gpuData->acousticInterfaceFaceCentersDevice[3 * faceIndex + 2]};
  const double normal[3] = {
      gpuData->acousticInterfaceFaceNormalsDevice[3 * faceIndex + 0],
      gpuData->acousticInterfaceFaceNormalsDevice[3 * faceIndex + 1],
      gpuData->acousticInterfaceFaceNormalsDevice[3 * faceIndex + 2]};
  const int32_t adjacent_eindex =
      (gpuData->acousticInterfaceAdjacentElementIndexesDevice != NULL)
          ? gpuData->acousticInterfaceAdjacentElementIndexesDevice[faceIndex]
          : -1;
  const double adjacent_center[3] = {
      (gpuData->acousticInterfaceAdjacentCentersDevice != NULL)
          ? gpuData->acousticInterfaceAdjacentCentersDevice[3 * faceIndex + 0]
          : center[0],
      (gpuData->acousticInterfaceAdjacentCentersDevice != NULL)
          ? gpuData->acousticInterfaceAdjacentCentersDevice[3 * faceIndex + 1]
          : center[1],
      (gpuData->acousticInterfaceAdjacentCentersDevice != NULL)
          ? gpuData->acousticInterfaceAdjacentCentersDevice[3 * faceIndex + 2]
          : center[2]};
  const int32_t force_eindex = (adjacent_eindex >= 0) ? adjacent_eindex : eindex;
  const double* force_center = (adjacent_eindex >= 0) ? adjacent_center : center;
  double pressure = 0.0;

  for (int i = 0; i < 8; i++) {
    const int32_t lnid =
        gpuData->lnidArrayDevice[eindex + i * gpuData->lenum];
    const double weight =
        acoustic_shape_value_at(i, center[0], center[1], center[2]);

    pressure += weight * (double)gpuData->pressureTm1Device[lnid];
  }

  for (int i = 0; i < 8; i++) {
    const double weight =
        acoustic_shape_value_at(i, force_center[0], force_center[1],
                                force_center[2]);
    if (weight == 0.0) {
      continue;
    }

    const int32_t lnid =
        gpuData->lnidArrayDevice[force_eindex + i * gpuData->lenum];
    fvector_t* nodalForce = gpuData->forceDevice + lnid;
    atomicAdd((solver_float*)&nodalForce->f[0],
              (solver_float)(face_coeff * weight * pressure * normal[0]));
    atomicAdd((solver_float*)&nodalForce->f[1],
              (solver_float)(face_coeff * weight * pressure * normal[1]));
    atomicAdd((solver_float*)&nodalForce->f[2],
              (solver_float)(face_coeff * weight * pressure * normal[2]));
  }
}

__global__ void kernelAcousticPressureForce(int32_t lenum,
                                            gpu_data_t* gpuData,
                                            double deltaT2)
{
  int32_t eindex = (blockIdx.x * blockDim.x) + threadIdx.x;
  solver_float pressure[8];
  (void)deltaT2;

  if (eindex >= lenum) {
    return;
  }
  if (gpuData->acousticElementFlagsArrayDevice[eindex] == 0) {
    return;
  }

  /* CPML: acoustic PML elements are handled by kernelPMLAcousticStep. */
  if (gpuData->pmlFlagArrayDevice != NULL &&
      gpuData->pmlFlagArrayDevice[eindex] != 0) {
    return;
  }

  for (int i = 0; i < 8; i++) {
    const int32_t lnid = gpuData->lnidArrayDevice[eindex + i * lenum];
    pressure[i] = gpuData->pressureTm1Device[lnid];
  }

  const double stiffness_coeff =
      gpuData->acousticPressureCoeffArrayDevice[eindex];

  for (int i = 0; i < 8; i++) {
    double local = 0.0;
    const int32_t lnid = gpuData->lnidArrayDevice[eindex + i * lenum];

    for (int j = 0; j < 8; j++) {
      local += acoustic_scalar_stiffness_shape(i, j) * pressure[j];
    }

    atomicAdd((solver_float*)(gpuData->pressureForceDevice + lnid),
              (solver_float)(-stiffness_coeff * local));
  }
}

__global__ void kernelAcousticPressureInterfaceForce(int32_t faceCount,
                                                     gpu_data_t* gpuData,
                                                     double deltaT2)
{
  int32_t faceIndex = (blockIdx.x * blockDim.x) + threadIdx.x;
  const double invDeltaT2 = 1.0 / deltaT2;

  if ((faceIndex >= faceCount) || (gpuData->tm3Device == NULL)) {
    return;
  }

  const int32_t eindex =
      gpuData->acousticInterfaceElementIdsDevice[faceIndex];
  const double face_coeff =
      gpuData->acousticInterfaceFaceAreasDevice[faceIndex];
  const double center[3] = {
      gpuData->acousticInterfaceFaceCentersDevice[3 * faceIndex + 0],
      gpuData->acousticInterfaceFaceCentersDevice[3 * faceIndex + 1],
      gpuData->acousticInterfaceFaceCentersDevice[3 * faceIndex + 2]};
  const double normal[3] = {
      gpuData->acousticInterfaceFaceNormalsDevice[3 * faceIndex + 0],
      gpuData->acousticInterfaceFaceNormalsDevice[3 * faceIndex + 1],
      gpuData->acousticInterfaceFaceNormalsDevice[3 * faceIndex + 2]};
  const int32_t adjacent_eindex =
      (gpuData->acousticInterfaceAdjacentElementIndexesDevice != NULL)
          ? gpuData->acousticInterfaceAdjacentElementIndexesDevice[faceIndex]
          : -1;
  const double adjacent_center[3] = {
      (gpuData->acousticInterfaceAdjacentCentersDevice != NULL)
          ? gpuData->acousticInterfaceAdjacentCentersDevice[3 * faceIndex + 0]
          : center[0],
      (gpuData->acousticInterfaceAdjacentCentersDevice != NULL)
          ? gpuData->acousticInterfaceAdjacentCentersDevice[3 * faceIndex + 1]
          : center[1],
      (gpuData->acousticInterfaceAdjacentCentersDevice != NULL)
          ? gpuData->acousticInterfaceAdjacentCentersDevice[3 * faceIndex + 2]
          : center[2]};
  const int32_t accel_eindex =
      (adjacent_eindex >= 0) ? adjacent_eindex : eindex;
  const double* accel_center =
      (adjacent_eindex >= 0) ? adjacent_center : center;
  double accel_n = 0.0;

  for (int i = 0; i < 8; i++) {
    const double weight =
        acoustic_shape_value_at(i, accel_center[0], accel_center[1],
                                accel_center[2]);
    if (weight == 0.0) {
      continue;
    }

    const int32_t lnid =
        gpuData->lnidArrayDevice[accel_eindex + i * gpuData->lenum];
    const fvector_t u_next = gpuData->tm2Device[lnid];
    const fvector_t u_curr = gpuData->tm1Device[lnid];
    const fvector_t u_prev = gpuData->tm3Device[lnid];

    accel_n += weight *
               ((((double)u_next.f[0] - 2.0 * (double)u_curr.f[0] +
                  (double)u_prev.f[0]) *
                     invDeltaT2 * normal[0]) +
                (((double)u_next.f[1] - 2.0 * (double)u_curr.f[1] +
                  (double)u_prev.f[1]) *
                     invDeltaT2 * normal[1]) +
                (((double)u_next.f[2] - 2.0 * (double)u_curr.f[2] +
                  (double)u_prev.f[2]) *
                     invDeltaT2 * normal[2]));
  }

  for (int i = 0; i < 8; i++) {
    const double weight =
        acoustic_shape_value_at(i, center[0], center[1], center[2]);
    if (weight == 0.0) {
      continue;
    }

    const int32_t lnid =
        gpuData->lnidArrayDevice[eindex + i * gpuData->lenum];
    atomicAdd((solver_float*)(gpuData->pressureForceDevice + lnid),
              (solver_float)(-face_coeff * weight * deltaT2 * accel_n));
  }
}

__global__ void kernelAcousticPressureUpdate(int32_t nharbored,
                                             gpu_data_t* gpuData)
{
  lnid_t nindex = (blockIdx.x * blockDim.x) + threadIdx.x;

  if (nindex >= nharbored) {
    return;
  }

  if ((gpuData->pressureFreeSurfaceNodeFlagsDevice != NULL) &&
      (gpuData->pressureFreeSurfaceNodeFlagsDevice[nindex] != 0)) {
    gpuData->pressureTm2Device[nindex] = 0.0;
    gpuData->pressureForceDevice[nindex] = 0.0;
    return;
  }

  const solver_float mass = gpuData->pressureMassDevice[nindex];
  solver_float p_new = 0.0;

  if (mass > 0.0) {
    p_new = (gpuData->pressureForceDevice[nindex] +
             (solver_float)2.0 * mass * gpuData->pressureTm1Device[nindex] -
             mass * gpuData->pressureTm2Device[nindex]) /
            mass;
  }

  gpuData->pressureTm2Device[nindex] = p_new;
  gpuData->pressureForceDevice[nindex] = 0.0;
}


/* Damping Calc-Conv Kernel
   FLOPs: 63M + 54A + 4T
   MEM  : 4 + 8*(48 + 1) + 2*(8*5 + 8*(6*2 + 3)) = 716
*/
__global__  void kernelDampingCalcConvBkt2(int32_t lenum,
                                           int32_t numcomp,
                                           gpu_data_t *gpuData,
                                           double rmax,
                                           double shear_coef_1,
                                           double shear_coef_2,
                                           double shear_exp_coef_0,
                                           double shear_coef_3,
                                           double shear_coef_4,
                                           double shear_exp_coef_1,
                                           double kappa_coef_1,
                                           double kappa_coef_2,
                                           double kappa_exp_coef_0,
                                           double kappa_coef_3,
                                           double kappa_coef_4,
                                           double kappa_exp_coef_1,
                                           int32_t dampingElemCount,
                                           int32_t* dampingElemMapperDevice)
{
  int i;
  int32_t lin_eindex = (blockIdx.x * blockDim.x) + threadIdx.x;
  int32_t eindex;

  if (lin_eindex >= dampingElemCount) {
    return;
  }

  eindex = dampingElemMapperDevice[lin_eindex];

  /* CPML: skip BKT damping on PML elements. */
  if (gpuData->pmlFlagArrayDevice != NULL &&
      gpuData->pmlFlagArrayDevice[eindex] != 0) {
    return;
  }

  const double b_shear = gpuData->b_shearArrayDevice[eindex];
  const double a0_shear = gpuData->a0_shearArrayDevice[eindex];
  const double a1_shear = gpuData->a1_shearArrayDevice[eindex];
  const double coef_shear = (rmax != 0.0) ? (b_shear / rmax) : 0.0;
  const int need_conv_history_shear = ((a0_shear != 0.0) || (a1_shear != 0.0));
  const int do_conv_update_shear =
      (need_conv_history_shear && (shear_coef_1 != 0.0) &&
       (shear_coef_3 != 0.0));
  const int do_damping_calc_shear =
      (need_conv_history_shear || (b_shear != 0.0));

  const double b_kappa = gpuData->b_kappaArrayDevice[eindex];
  const double a0_kappa = gpuData->a0_kappaArrayDevice[eindex];
  const double a1_kappa = gpuData->a1_kappaArrayDevice[eindex];
  const double coef_kappa = (rmax != 0.0) ? (b_kappa / rmax) : 0.0;
  const int need_conv_history_kappa =
      ((a0_kappa != 0.0) || (a1_kappa != 0.0));
  const int do_conv_update_kappa =
      (need_conv_history_kappa && (kappa_coef_1 != 0.0) &&
       (kappa_coef_3 != 0.0));
  const int do_damping_calc_kappa =
      (need_conv_history_kappa || (b_kappa != 0.0));

  for (i = 0; i < 8; i++) {
    const int32_t offset0 = eindex + (i * lenum);
    const int32_t offset1 = offset0 + numcomp;
    const int32_t offset2 = offset1 + numcomp;
    const int32_t lnid = *(gpuData->lnidArrayDevice + offset0);
    const __restrict__ fvector_t *disp = gpuData->tm1Device + lnid;
    const solver_float tm1x = disp->f[0];
    const solver_float tm1y = disp->f[1];
    const solver_float tm1z = disp->f[2];

    disp = gpuData->tm2Device + lnid;
    const solver_float tm2x = disp->f[0];
    const solver_float tm2y = disp->f[1];
    const solver_float tm2z = disp->f[2];

    if (need_conv_history_shear) {
      solver_float *shear0_x =
          (solver_float *)gpuData->conv_shear_1Device + offset0;
      solver_float *shear0_y =
          (solver_float *)gpuData->conv_shear_1Device + offset1;
      solver_float *shear0_z =
          (solver_float *)gpuData->conv_shear_1Device + offset2;
      solver_float *shear1_x =
          (solver_float *)gpuData->conv_shear_2Device + offset0;
      solver_float *shear1_y =
          (solver_float *)gpuData->conv_shear_2Device + offset1;
      solver_float *shear1_z =
          (solver_float *)gpuData->conv_shear_2Device + offset2;
      solver_float conv0_x = *shear0_x;
      solver_float conv0_y = *shear0_y;
      solver_float conv0_z = *shear0_z;
      solver_float conv1_x = *shear1_x;
      solver_float conv1_y = *shear1_y;
      solver_float conv1_z = *shear1_z;

      if (do_conv_update_shear) {
        conv0_x = shear_coef_2 * tm1x + shear_coef_1 * tm2x +
                  shear_exp_coef_0 * conv0_x;
        conv0_y = shear_coef_2 * tm1y + shear_coef_1 * tm2y +
                  shear_exp_coef_0 * conv0_y;
        conv0_z = shear_coef_2 * tm1z + shear_coef_1 * tm2z +
                  shear_exp_coef_0 * conv0_z;
        conv1_x = shear_coef_4 * tm1x + shear_coef_3 * tm2x +
                  shear_exp_coef_1 * conv1_x;
        conv1_y = shear_coef_4 * tm1y + shear_coef_3 * tm2y +
                  shear_exp_coef_1 * conv1_y;
        conv1_z = shear_coef_4 * tm1z + shear_coef_3 * tm2z +
                  shear_exp_coef_1 * conv1_z;

        *shear0_x = conv0_x;
        *shear0_y = conv0_y;
        *shear0_z = conv0_z;
        *shear1_x = conv1_x;
        *shear1_y = conv1_y;
        *shear1_z = conv1_z;
      }

      solver_float *convVec0 =
          (solver_float *)gpuData->shearVectorDevice + offset0;
      solver_float *convVec1 =
          (solver_float *)gpuData->shearVectorDevice + offset1;
      solver_float *convVec2 =
          (solver_float *)gpuData->shearVectorDevice + offset2;
      *convVec0 = coef_shear * (tm1x - tm2x) -
                  (a0_shear * conv0_x + a1_shear * conv1_x) + tm1x;
      *convVec1 = coef_shear * (tm1y - tm2y) -
                  (a0_shear * conv0_y + a1_shear * conv1_y) + tm1y;
      *convVec2 = coef_shear * (tm1z - tm2z) -
                  (a0_shear * conv0_z + a1_shear * conv1_z) + tm1z;
    } else {
      solver_float *convVec0 =
          (solver_float *)gpuData->shearVectorDevice + offset0;
      solver_float *convVec1 =
          (solver_float *)gpuData->shearVectorDevice + offset1;
      solver_float *convVec2 =
          (solver_float *)gpuData->shearVectorDevice + offset2;
      if (do_damping_calc_shear) {
        *convVec0 = coef_shear * (tm1x - tm2x) + tm1x;
        *convVec1 = coef_shear * (tm1y - tm2y) + tm1y;
        *convVec2 = coef_shear * (tm1z - tm2z) + tm1z;
      } else {
        *convVec0 = tm1x;
        *convVec1 = tm1y;
        *convVec2 = tm1z;
      }
    }

    if (need_conv_history_kappa) {
      solver_float *kappa0_x =
          (solver_float *)gpuData->conv_kappa_1Device + offset0;
      solver_float *kappa0_y =
          (solver_float *)gpuData->conv_kappa_1Device + offset1;
      solver_float *kappa0_z =
          (solver_float *)gpuData->conv_kappa_1Device + offset2;
      solver_float *kappa1_x =
          (solver_float *)gpuData->conv_kappa_2Device + offset0;
      solver_float *kappa1_y =
          (solver_float *)gpuData->conv_kappa_2Device + offset1;
      solver_float *kappa1_z =
          (solver_float *)gpuData->conv_kappa_2Device + offset2;
      solver_float conv0_x = *kappa0_x;
      solver_float conv0_y = *kappa0_y;
      solver_float conv0_z = *kappa0_z;
      solver_float conv1_x = *kappa1_x;
      solver_float conv1_y = *kappa1_y;
      solver_float conv1_z = *kappa1_z;

      if (do_conv_update_kappa) {
        conv0_x = kappa_coef_2 * tm1x + kappa_coef_1 * tm2x +
                  kappa_exp_coef_0 * conv0_x;
        conv0_y = kappa_coef_2 * tm1y + kappa_coef_1 * tm2y +
                  kappa_exp_coef_0 * conv0_y;
        conv0_z = kappa_coef_2 * tm1z + kappa_coef_1 * tm2z +
                  kappa_exp_coef_0 * conv0_z;
        conv1_x = kappa_coef_4 * tm1x + kappa_coef_3 * tm2x +
                  kappa_exp_coef_1 * conv1_x;
        conv1_y = kappa_coef_4 * tm1y + kappa_coef_3 * tm2y +
                  kappa_exp_coef_1 * conv1_y;
        conv1_z = kappa_coef_4 * tm1z + kappa_coef_3 * tm2z +
                  kappa_exp_coef_1 * conv1_z;

        *kappa0_x = conv0_x;
        *kappa0_y = conv0_y;
        *kappa0_z = conv0_z;
        *kappa1_x = conv1_x;
        *kappa1_y = conv1_y;
        *kappa1_z = conv1_z;
      }

      solver_float *convVec0 =
          (solver_float *)gpuData->kappaVectorDevice + offset0;
      solver_float *convVec1 =
          (solver_float *)gpuData->kappaVectorDevice + offset1;
      solver_float *convVec2 =
          (solver_float *)gpuData->kappaVectorDevice + offset2;
      *convVec0 = coef_kappa * (tm1x - tm2x) -
                  (a0_kappa * conv0_x + a1_kappa * conv1_x) + tm1x;
      *convVec1 = coef_kappa * (tm1y - tm2y) -
                  (a0_kappa * conv0_y + a1_kappa * conv1_y) + tm1y;
      *convVec2 = coef_kappa * (tm1z - tm2z) -
                  (a0_kappa * conv0_z + a1_kappa * conv1_z) + tm1z;
    } else {
      solver_float *convVec0 =
          (solver_float *)gpuData->kappaVectorDevice + offset0;
      solver_float *convVec1 =
          (solver_float *)gpuData->kappaVectorDevice + offset1;
      solver_float *convVec2 =
          (solver_float *)gpuData->kappaVectorDevice + offset2;
      if (do_damping_calc_kappa) {
        *convVec0 = coef_kappa * (tm1x - tm2x) + tm1x;
        *convVec1 = coef_kappa * (tm1y - tm2y) + tm1y;
        *convVec2 = coef_kappa * (tm1z - tm2z) + tm1z;
      } else {
        *convVec0 = tm1x;
        *convVec1 = tm1y;
        *convVec2 = tm1z;
      }
    }
  }
}

__global__  void kernelDampingCalcConv(int32_t lenum,
				       int32_t numcomp,
				       gpu_data_t *gpuData,
               double rmax,
               int32_t dampingType,
               int32_t dampingElemCount,
               int32_t* dampingElemMapperDevice)
{
  int i;
  int32_t   lin_eindex = (blockIdx.x * blockDim.x) + threadIdx.x;
  int32_t   eindex;

  /* Since number of elements may not be exactly divisible by block size,
     check that we are not off the end of the element array */
  if (lin_eindex >= dampingElemCount) {
    return;
  }

  eindex = dampingElemMapperDevice[lin_eindex];

  /* CPML: PML elements run their own elastic/acoustic memory recurrence;
   * the BKT (Q-attenuation) memory variables are skipped on PML elements
   * to keep the two memory-variable systems decoupled. */
  if (gpuData->pmlFlagArrayDevice != NULL &&
      gpuData->pmlFlagArrayDevice[eindex] != 0) {
    return;
  }

  double g0, g1, g2;
  double coef_1_shear, coef_2_shear, coef_3_shear, coef_4_shear;
  double coef_5_shear, coef_6_shear;
  double exp_coef_0_shear, exp_coef_1_shear, exp_coef_2_shear;
  double coef_1_kappa, coef_2_kappa, coef_3_kappa, coef_4_kappa;
  double coef_5_kappa, coef_6_kappa;
  double exp_coef_0_kappa, exp_coef_1_kappa, exp_coef_2_kappa;
  fvector_t     tm1Disp, tm2Disp;
  solver_float *f0_tm1[3], *f1_tm1[3], *f2_tm1[3], *convVec[3];

  double coef_shear, a0_shear, a1_shear, a2_shear;
  double coef_kappa, a0_kappa, a1_kappa, a2_kappa;
  int do_conv_update_shear, do_damping_calc_shear;
  int do_conv_update_kappa, do_damping_calc_kappa;
  int use_bkt3;
  int32_t offset0;

  use_bkt3 = (dampingType >= BKT3);

  // SHEAR RELATED CONVOLUTION  
  g0 = gpuData->g0_shearArrayDevice[eindex];
  coef_1_shear = g0 / 2.;
  coef_2_shear = coef_1_shear * ( 1. - g0 );
  exp_coef_0_shear = exp( -g0 );

  g1 = gpuData->g1_shearArrayDevice[eindex];
  coef_3_shear = g1 / 2.;
  coef_4_shear = coef_3_shear * ( 1. - g1 );
  exp_coef_1_shear = exp( -g1 );

  if (use_bkt3) {
    g2 = gpuData->g2_shearArrayDevice[eindex];
    coef_5_shear = g2 / 2.;
    coef_6_shear = coef_5_shear * (1. - g2);
    exp_coef_2_shear = exp(-g2);
  }

  a0_shear = gpuData->a0_shearArrayDevice[eindex];
  a1_shear = gpuData->a1_shearArrayDevice[eindex];
  a2_shear = use_bkt3 ? gpuData->a2_shearArrayDevice[eindex] : 0.0;
  coef_shear = (rmax != 0.0) ? (gpuData->b_shearArrayDevice[eindex] / rmax) : 0.0;

  do_conv_update_shear = ((g0 != 0) && (g1 != 0));
  do_damping_calc_shear =
      ((a0_shear + a1_shear + a2_shear + gpuData->b_shearArrayDevice[eindex]) !=
       0);

  // DILATION RELATED CONVOLUTION
  g0 = gpuData->g0_kappaArrayDevice[eindex];
  coef_1_kappa = g0 / 2.;
  coef_2_kappa = coef_1_kappa * ( 1. - g0 );
  exp_coef_0_kappa = exp( -g0 );

  g1 = gpuData->g1_kappaArrayDevice[eindex];
  coef_3_kappa = g1 / 2.;
  coef_4_kappa = coef_3_kappa * ( 1. - g1 );
  exp_coef_1_kappa = exp( -g1 );

  if (use_bkt3) {
    g2 = gpuData->g2_kappaArrayDevice[eindex];
    coef_5_kappa = g2 / 2.;
    coef_6_kappa = coef_5_kappa * (1. - g2);
    exp_coef_2_kappa = exp(-g2);
  }

  a0_kappa = gpuData->a0_kappaArrayDevice[eindex];
  a1_kappa = gpuData->a1_kappaArrayDevice[eindex];
  a2_kappa = use_bkt3 ? gpuData->a2_kappaArrayDevice[eindex] : 0.0;
  coef_kappa = (rmax != 0.0) ? (gpuData->b_kappaArrayDevice[eindex] / rmax) : 0.0;

  do_conv_update_kappa = ((g0 != 0) && (g1 != 0));
  do_damping_calc_kappa =
      ((a0_kappa + a1_kappa + a2_kappa + gpuData->b_kappaArrayDevice[eindex]) !=
       0);

  for (i = 0, offset0 = eindex; i < 8; i++, offset0 += lenum) {
    int32_t offset1 = offset0 + numcomp;
    int32_t offset2 = offset1 + numcomp;
    int32_t lnid = *(gpuData->lnidArrayDevice + offset0);

    const __restrict__ fvector_t *disp = gpuData->tm1Device + lnid;
    tm1Disp.f[0] = disp->f[0];
    tm1Disp.f[1] = disp->f[1];
    tm1Disp.f[2] = disp->f[2];

    disp = gpuData->tm2Device + lnid;
    tm2Disp.f[0] = disp->f[0];
    tm2Disp.f[1] = disp->f[1];
    tm2Disp.f[2] = disp->f[2];

    /* Read aligned tm2disp from the shear vector storage */
    //tm2Disp.f[0] = *((solver_float *)gpuData->shearVectorDevice + offset0);
    //tm2Disp.f[1] = *((solver_float *)gpuData->shearVectorDevice + offset1);
    //tm2Disp.f[2] = *((solver_float *)gpuData->shearVectorDevice + offset2);

    f0_tm1[0] = (solver_float *)gpuData->conv_shear_1Device + offset0;
    f0_tm1[1] = (solver_float *)gpuData->conv_shear_1Device + offset1;
    f0_tm1[2] = (solver_float *)gpuData->conv_shear_1Device + offset2;

    f1_tm1[0] = (solver_float *)gpuData->conv_shear_2Device + offset0;
    f1_tm1[1] = (solver_float *)gpuData->conv_shear_2Device + offset1;
    f1_tm1[2] = (solver_float *)gpuData->conv_shear_2Device + offset2;

    if (use_bkt3) {
      f2_tm1[0] = (solver_float *)gpuData->conv_shear_3Device + offset0;
      f2_tm1[1] = (solver_float *)gpuData->conv_shear_3Device + offset1;
      f2_tm1[2] = (solver_float *)gpuData->conv_shear_3Device + offset2;
    }
    
    /* Update the convolution */
    if (do_conv_update_shear) {

      *f0_tm1[0] = coef_2_shear * tm1Disp.f[0] + coef_1_shear * tm2Disp.f[0] + exp_coef_0_shear * (*f0_tm1[0]);
      *f1_tm1[0] = coef_4_shear * tm1Disp.f[0] + coef_3_shear * tm2Disp.f[0] + exp_coef_1_shear * (*f1_tm1[0]);
      
      *f0_tm1[1] = coef_2_shear * tm1Disp.f[1] + coef_1_shear * tm2Disp.f[1] + exp_coef_0_shear * (*f0_tm1[1]);
      *f1_tm1[1] = coef_4_shear * tm1Disp.f[1] + coef_3_shear * tm2Disp.f[1] + exp_coef_1_shear * (*f1_tm1[1]);
      
      *f0_tm1[2] = coef_2_shear * tm1Disp.f[2] + coef_1_shear * tm2Disp.f[2] + exp_coef_0_shear * (*f0_tm1[2]);
      *f1_tm1[2] = coef_4_shear * tm1Disp.f[2] + coef_3_shear * tm2Disp.f[2] + exp_coef_1_shear * (*f1_tm1[2]);

      if (use_bkt3) {
        *f2_tm1[0] = coef_6_shear * tm1Disp.f[0] + coef_5_shear * tm2Disp.f[0] + exp_coef_2_shear * (*f2_tm1[0]);
        *f2_tm1[1] = coef_6_shear * tm1Disp.f[1] + coef_5_shear * tm2Disp.f[1] + exp_coef_2_shear * (*f2_tm1[1]);
        *f2_tm1[2] = coef_6_shear * tm1Disp.f[2] + coef_5_shear * tm2Disp.f[2] + exp_coef_2_shear * (*f2_tm1[2]);
      }
      
    }

    convVec[0] = (solver_float *)gpuData->shearVectorDevice + offset0;
    convVec[1] = (solver_float *)gpuData->shearVectorDevice + offset1;
    convVec[2] = (solver_float *)gpuData->shearVectorDevice + offset2;
    
    /* Construct the damping vector */
    if (do_damping_calc_shear ) {

      *convVec[0] = coef_shear * (tm1Disp.f[0] - tm2Disp.f[0])
	      - (a0_shear * (*f0_tm1[0]) + a1_shear * (*f1_tm1[0])) + tm1Disp.f[0];

      *convVec[1] = coef_shear * (tm1Disp.f[1] - tm2Disp.f[1])
	      - (a0_shear * (*f0_tm1[1]) + a1_shear * (*f1_tm1[1])) + tm1Disp.f[1];

      *convVec[2] = coef_shear * (tm1Disp.f[2] - tm2Disp.f[2])
	      - (a0_shear * (*f0_tm1[2]) + a1_shear * (*f1_tm1[2])) + tm1Disp.f[2];

      if (use_bkt3) {
        *convVec[0] -= a2_shear * (*f2_tm1[0]);
        *convVec[1] -= a2_shear * (*f2_tm1[1]);
        *convVec[2] -= a2_shear * (*f2_tm1[2]);
      }

    } else {

      *convVec[0] = tm1Disp.f[0];
      *convVec[1] = tm1Disp.f[1];
      *convVec[2] = tm1Disp.f[2];

    }

    f0_tm1[0] = (solver_float *)gpuData->conv_kappa_1Device + offset0;
    f0_tm1[1] = (solver_float *)gpuData->conv_kappa_1Device + offset1;
    f0_tm1[2] = (solver_float *)gpuData->conv_kappa_1Device + offset2;

    f1_tm1[0] = (solver_float *)gpuData->conv_kappa_2Device + offset0;
    f1_tm1[1] = (solver_float *)gpuData->conv_kappa_2Device + offset1;
    f1_tm1[2] = (solver_float *)gpuData->conv_kappa_2Device + offset2;

    if (use_bkt3) {
      f2_tm1[0] = (solver_float *)gpuData->conv_kappa_3Device + offset0;
      f2_tm1[1] = (solver_float *)gpuData->conv_kappa_3Device + offset1;
      f2_tm1[2] = (solver_float *)gpuData->conv_kappa_3Device + offset2;
    }

    /* Update the convolution */
    if (do_conv_update_kappa) {

      *f0_tm1[0] = coef_2_kappa * tm1Disp.f[0] + coef_1_kappa * tm2Disp.f[0] + exp_coef_0_kappa * (*f0_tm1[0]);
      *f1_tm1[0] = coef_4_kappa * tm1Disp.f[0] + coef_3_kappa * tm2Disp.f[0] + exp_coef_1_kappa * (*f1_tm1[0]);

      *f0_tm1[1] = coef_2_kappa * tm1Disp.f[1] + coef_1_kappa * tm2Disp.f[1] + exp_coef_0_kappa * (*f0_tm1[1]);
      *f1_tm1[1] = coef_4_kappa * tm1Disp.f[1] + coef_3_kappa * tm2Disp.f[1] + exp_coef_1_kappa * (*f1_tm1[1]);

      *f0_tm1[2] = coef_2_kappa * tm1Disp.f[2] + coef_1_kappa * tm2Disp.f[2] + exp_coef_0_kappa * (*f0_tm1[2]);
      *f1_tm1[2] = coef_4_kappa * tm1Disp.f[2] + coef_3_kappa * tm2Disp.f[2] + exp_coef_1_kappa * (*f1_tm1[2]);

      if (use_bkt3) {
        *f2_tm1[0] = coef_6_kappa * tm1Disp.f[0] + coef_5_kappa * tm2Disp.f[0] + exp_coef_2_kappa * (*f2_tm1[0]);
        *f2_tm1[1] = coef_6_kappa * tm1Disp.f[1] + coef_5_kappa * tm2Disp.f[1] + exp_coef_2_kappa * (*f2_tm1[1]);
        *f2_tm1[2] = coef_6_kappa * tm1Disp.f[2] + coef_5_kappa * tm2Disp.f[2] + exp_coef_2_kappa * (*f2_tm1[2]);
      }

    }

    convVec[0] = (solver_float *)gpuData->kappaVectorDevice + offset0;
    convVec[1] = (solver_float *)gpuData->kappaVectorDevice + offset1;
    convVec[2] = (solver_float *)gpuData->kappaVectorDevice + offset2;
    
    /* Construct the damping vector */
    if (do_damping_calc_kappa) {

      *convVec[0] = coef_kappa * (tm1Disp.f[0] - tm2Disp.f[0])
	      - (a0_kappa * (*f0_tm1[0]) + a1_kappa * (*f1_tm1[0])) + tm1Disp.f[0];

      *convVec[1] = coef_kappa * (tm1Disp.f[1] - tm2Disp.f[1])
	      - (a0_kappa * (*f0_tm1[1]) + a1_kappa * (*f1_tm1[1])) + tm1Disp.f[1];

      *convVec[2] = coef_kappa * (tm1Disp.f[2] - tm2Disp.f[2])
	      - (a0_kappa * (*f0_tm1[2]) + a1_kappa * (*f1_tm1[2])) + tm1Disp.f[2];

      if (use_bkt3) {
        *convVec[0] -= a2_kappa * (*f2_tm1[0]);
        *convVec[1] -= a2_kappa * (*f2_tm1[1]);
        *convVec[2] -= a2_kappa * (*f2_tm1[2]);
      }

    } else {

      *convVec[0] = tm1Disp.f[0];
      *convVec[1] = tm1Disp.f[1];
      *convVec[2] = tm1Disp.f[2];

    }

  }

  return;
}



/* Damping Calc-force Kernel
   FLOPs: (4M + 1A) + 2*(147A) + (60M + 48A) + (21M + 36A) + 168A = 85M + 547A
   MEM  : 2*(8 * 24 + 8) + 8(24 + 1) = 592
*/
__global__
//__launch_bounds__(512, 1)
void kernelDampingCalcLocal(int32_t lenum,
			    gpu_data_t *gpuData,
			    int32_t dampingElemCount,
			    int32_t* dampingElemMapperDevice)
{
  int       i;
  int32_t   lin_eindex = (blockIdx.x * blockDim.x) + threadIdx.x;
  int32_t   eindex;
  fvector_t damping_vector[8];

  /* Since number of elements may not be exactly divisible by block size,
     check that we are not off the end of the element array */
  if (lin_eindex >= dampingElemCount) {
    return;
  }

  eindex = dampingElemMapperDevice[lin_eindex];

  /* CPML: skip BKT damping force on PML elements. */
  if (gpuData->pmlFlagArrayDevice != NULL &&
      gpuData->pmlFlagArrayDevice[eindex] != 0) {
    return;
  }

  int         numcomp = lenum * 8;
  int32_t     *lnid;
  int         do_write = 0;
  double      atu[24], firstVec[24];
  double      mu, kappa;
  solver_float *cVec;
  int         do_contrib;
  uint8_t acousticFaceFlags = gpuData->acousticFaceFlagsArrayDevice[eindex];

  /* SHEAR CONTRIBUTION */
  cVec = (solver_float *)gpuData->shearVectorDevice + eindex;
  do_contrib = 0;
  for (i = 0; i < 8; i++) {
    damping_vector[i].f[0] = *cVec;
    do_contrib = (fabs( damping_vector[i].f[0] ) > UNDERFLOW_CAP_STIFFNESS) ? 1 : do_contrib;
    firstVec[i] = 0.0;
    damping_vector[i].f[1] = *(cVec + numcomp);
    do_contrib = (fabs( damping_vector[i].f[1] ) > UNDERFLOW_CAP_STIFFNESS) ? 1 : do_contrib;
    firstVec[i+8] = 0.0;
    damping_vector[i].f[2] = *(cVec + 2*numcomp);
    do_contrib = (fabs( damping_vector[i].f[2] ) > UNDERFLOW_CAP_STIFFNESS) ? 1 : do_contrib;
    firstVec[i+16] = 0.0;
    cVec += lenum;
  }

  mu = -0.5625 * gpuData->c1ArrayDevice[eindex];  
  if (do_contrib) {
    //if(vector_is_zero( damping_vector ) != 0) {
    aTransposeU( damping_vector, atu );
    firstVector_mu( atu, firstVec, mu);
    do_write = 1;
  }

  /* DILATION CONTRIBUTION */
  cVec = (solver_float *)gpuData->kappaVectorDevice + eindex;
  do_contrib = 0;
  for (i = 0; i < 8; i++) {
    damping_vector[i].f[0] = *cVec;
    do_contrib = (fabs( damping_vector[i].f[0] ) > UNDERFLOW_CAP_STIFFNESS) ? 1 : do_contrib;
    damping_vector[i].f[1] = *(cVec + numcomp);
    do_contrib = (fabs( damping_vector[i].f[1] ) > UNDERFLOW_CAP_STIFFNESS) ? 1 : do_contrib;
    damping_vector[i].f[2] = *(cVec + 2*numcomp);
    do_contrib = (fabs( damping_vector[i].f[2] ) > UNDERFLOW_CAP_STIFFNESS) ? 1 : do_contrib;
    cVec += lenum;
  }

  kappa = -0.5625 * (gpuData->c2ArrayDevice[eindex] + 
		     2. / 3. * gpuData->c1ArrayDevice[eindex]);
  if (do_contrib) {
    if (acousticFaceFlags != 0) {
      acoustic_project_interface_nodes(damping_vector, acousticFaceFlags);
    }
    //if(vector_is_zero( damping_vector ) != 0) {
    aTransposeU( damping_vector, atu );
    firstVector_kappa( atu, firstVec, kappa);
    do_write = 1;
  }

  if (!do_write) {
    return;
  }

  /* Re-use the damping vector for the local force */
  au( damping_vector, firstVec );

  if (acousticFaceFlags != 0) {
    acoustic_project_interface_nodes(damping_vector, acousticFaceFlags);
  }
  
  /* Sum the nodal forces */
  lnid = &(gpuData->lnidArrayDevice[eindex]);
  for (i = 0; i < 8; i++) {
    fvector_t *nodalForce = gpuData->forceDevice + *lnid;

    atomicAdd((solver_float *)&(nodalForce->f[0]),
	      (solver_float)damping_vector[i].f[0]);
    atomicAdd((solver_float *)&(nodalForce->f[1]),
	      (solver_float)damping_vector[i].f[1]);
    atomicAdd((solver_float *)&(nodalForce->f[2]),
	      (solver_float)damping_vector[i].f[2]);

    lnid += lenum;
  }

  return;
}


/* Displacement Calculation Kernel
   FLOPs: 3*(2A + 3M) = 6A + 9M
   MEM  : 8(3*3) + 8 + 8*6 + 8*3 = 128
*/
__global__  void kernelDispCalc(int32_t nharbored,
				gpu_data_t* gpuData,
				noyesflag_t printAccel)
{
    (void)nharbored;
    (void)printAccel;

    lnid_t nindex = (blockIdx.x * blockDim.x) + threadIdx.x;
    solver_float mass2_minusaM, mass_minusaM, mass_simple;
    fvector_t tm1Disp, tm2Disp, nodalForce;

    if (nindex >= gpuData->nharbored) {
      return;
    }

    tm1Disp = gpuData->tm1Device[nindex];
    tm2Disp = gpuData->tm2Device[nindex];
    nodalForce = gpuData->forceDevice[nindex];

    if (gpuData->tm3Device != NULL) {
      gpuData->tm3Device[nindex] = tm2Disp;
    }

    mass_simple = gpuData->mass_simpleArrayDevice[nindex];
    if (mass_simple > 0) {
      solver_float inv_mass_simple = (solver_float)1.0 / mass_simple;

      mass2_minusaM = gpuData->mass2_minusaMArrayDevice[0][nindex];
      mass_minusaM = gpuData->mass_minusaMArrayDevice[0][nindex];
      tm2Disp.f[0] =
          (nodalForce.f[0] + mass2_minusaM * tm1Disp.f[0] -
           mass_minusaM * tm2Disp.f[0]) *
          inv_mass_simple;

      mass2_minusaM = gpuData->mass2_minusaMArrayDevice[1][nindex];
      mass_minusaM = gpuData->mass_minusaMArrayDevice[1][nindex];
      tm2Disp.f[1] =
          (nodalForce.f[1] + mass2_minusaM * tm1Disp.f[1] -
           mass_minusaM * tm2Disp.f[1]) *
          inv_mass_simple;

      mass2_minusaM = gpuData->mass2_minusaMArrayDevice[2][nindex];
      mass_minusaM = gpuData->mass_minusaMArrayDevice[2][nindex];
      tm2Disp.f[2] =
          (nodalForce.f[2] + mass2_minusaM * tm1Disp.f[2] -
           mass_minusaM * tm2Disp.f[2]) *
          inv_mass_simple;
    } else {
      tm2Disp.f[0] = 0;
      tm2Disp.f[1] = 0;
      tm2Disp.f[2] = 0;
    }

    gpuData->tm2Device[nindex] = tm2Disp;
}



/* tm2Disp Alignment Kernel
   FLOPs: 0 (overhead kernel)
   MEM  : 0 (overhead kernel)
*/
__global__  void kernelAlignTm2Disp(int32_t lenum,
				    int32_t numcomp,
				    gpu_data_t *gpuData)
{

  int i;
  int32_t   eindex = (blockIdx.x * blockDim.x) + threadIdx.x; // Element number
  int offset;

  /* Since number of elements may not be exactly divisible by block size,
     check that we are not off the end of the element array */
  if (eindex >= lenum) {
    return;
  }

  /* Align tm2Disp */
  for (i = 0, offset = eindex; i < 8; i++, offset += lenum) {
    int32_t lnid = *(gpuData->lnidArrayDevice + offset);
    const __restrict__ fvector_t *disp = gpuData->tm1Device + lnid;
    solver_float *dampVec = (solver_float *)gpuData->shearVectorDevice + offset;

    *dampVec = disp->f[0];
    *(dampVec + numcomp) = disp->f[1];
    *(dampVec + 2*numcomp) = disp->f[2];

  }

  return;
}


/* -------------------------------------------------------------------------- */
/*                         Efficient Method Utilities                         */
/* -------------------------------------------------------------------------- */


/**
 * For effective stiffness method:
 *
 * Check whether all components of a 3D vector are close to zero,
 * i.e., less than a small threshold.
 *
 * \return 1 when there is at least one "non-zero" component;
 *         0 when all the components are "zero".
 */
__host__ __device__ int vector_is_zero( const fvector_t* v )
{
    /*
     * For scalability studies, uncomment the immediate return.
     */

    /* return 1; */

    int i;
    const solver_float *tmp = (solver_float *)v;

    for (i = 0; i < 24; i++) {
      if (fabs( tmp[i] ) > UNDERFLOW_CAP_STIFFNESS) {
	return 1;
      }
    }

    return(0);

}


__host__ __device__ int vector_is_all_zero( const fvector_t* v )
{
    /*
     * For scalability studies, uncomment the immediate return.
     */

    /* return 1; */

  if ((fabs( v->f[0] ) > UNDERFLOW_CAP_STIFFNESS) ||
      (fabs( v->f[1] ) > UNDERFLOW_CAP_STIFFNESS) ||
      (fabs( v->f[2] ) > UNDERFLOW_CAP_STIFFNESS)) {
    return 1;
  }

    return 0;
}


/* FLOPs: 147A */
__host__ __device__ void aTransposeU( fvector_t* un, double* atu )
{
#ifdef  SINGLE_PRECISION_SOLVER
    double temp[24];
    double u[24];
    int    i, j;

    /* arrange displacement values in an array */
    for (i=0; i<8; i++) {
        for(j=0; j<3; j++) {
            temp[i*3 + j] = un[i].f[j];     /* u1 u2 u3 .... v1 v2 v3 ... z1 z2 
z3 */
        }
    }

    reformU( temp, u );
#else
    double u[24];

    reformU( (double *)un, u );
#endif

    /* atu[0] = u[0] + u[1] + u[2] + u[3] + u[4] + u[5] + u[6] + u[7]; */
    atu[0]  = 0;
    atu[1]  = -u[0] - u[1] - u[2] - u[3] + u[4] + u[5] + u[6] + u[7];
    atu[2]  = -u[0] - u[1] + u[2] + u[3] - u[4] - u[5] + u[6] + u[7];
    atu[3]  = -u[0] + u[1] - u[2] + u[3] - u[4] + u[5] - u[6] + u[7];
    atu[4]  =  u[0] + u[1] - u[2] - u[3] - u[4] - u[5] + u[6] + u[7];
    atu[5]  =  u[0] - u[1] + u[2] - u[3] - u[4] + u[5] - u[6] + u[7];
    atu[6]  =  u[0] - u[1] - u[2] + u[3] + u[4] - u[5] - u[6] + u[7];
    atu[7]  = -u[0] + u[1] + u[2] - u[3] + u[4] - u[5] - u[6] + u[7];

    /* atu[8] = u[8] + u[9] + u[10] + u[11] + u[12] + u[13] + u[14] + u[15]; */
    atu[8]  = 0;
    atu[9]  = -u[8] - u[9] - u[10] - u[11] + u[12] + u[13] + u[14] + u[15];
    atu[10] = -u[8] - u[9] + u[10] + u[11] - u[12] - u[13] + u[14] + u[15];
    atu[11] = -u[8] + u[9] - u[10] + u[11] - u[12] + u[13] - u[14] + u[15];
    atu[12] =  u[8] + u[9] - u[10] - u[11] - u[12] - u[13] + u[14] + u[15];
    atu[13] =  u[8] - u[9] + u[10] - u[11] - u[12] + u[13] - u[14] + u[15];
    atu[14] =  u[8] - u[9] - u[10] + u[11] + u[12] - u[13] - u[14] + u[15];
    atu[15] = -u[8] + u[9] + u[10] - u[11] + u[12] - u[13] - u[14] + u[15];

    /* atu[16] = u[16] + u[17] + u[18] + u[19] + u[20] + u[21] + u[22] + u[23]; */
    atu[16] = 0;
    atu[17] = -u[16] - u[17] - u[18] - u[19] + u[20] + u[21] + u[22] + u[23];
    atu[18] = -u[16] - u[17] + u[18] + u[19] - u[20] - u[21] + u[22] + u[23];
    atu[19] = -u[16] + u[17] - u[18] + u[19] - u[20] + u[21] - u[22] + u[23];
    atu[20] =  u[16] + u[17] - u[18] - u[19] - u[20] - u[21] + u[22] + u[23];
    atu[21] =  u[16] - u[17] + u[18] - u[19] - u[20] + u[21] - u[22] + u[23];
    atu[22] =  u[16] - u[17] - u[18] + u[19] + u[20] - u[21] - u[22] + u[23];
    atu[23] = -u[16] + u[17] + u[18] - u[19] + u[20] - u[21] - u[22] + u[23];
}


/* FLOPs: 48M + 33A */
__host__ __device__ void firstVector( const double* atu,
				      double* finalVector,
				      double a,
				      double c,
				      double b )
{
    finalVector[0] = 0;
    finalVector[1] = b * (atu[19] + atu[1]);
    finalVector[2] = b * (atu[11] + atu[2]);
    finalVector[3] = a * atu[3] + c * (atu[10] + atu[17]);
    finalVector[4] = b * (atu[13] + atu[22] + 2. * atu[4]) / 3.;
    finalVector[5] = ( (a + b) * atu[5] + c * atu[12] ) /3.;
    finalVector[6] = ( (a + b) * atu[6] + c * atu[20] ) /3.;
    finalVector[7] = ( (a + 2.*b) * atu[7] ) / 9.;

    finalVector[8] = 0;
    finalVector[9] = b * (atu[18] + atu[9]);
    finalVector[10] = a * atu[10] + c * (atu[3] + atu[17]);
    finalVector[11] = b * (atu[11] + atu[2]);
    finalVector[12] = ( (a + b) * atu[12] + c * atu[5] ) / 3.;
    finalVector[13] = b * (atu[4] + atu[22] + 2. * atu[13]) / 3.;
    finalVector[14] = ( (a + b) * atu[14] + c * atu[21] ) /3.;
    finalVector[15] = (a + 2. * b) * atu[15] / 9.;

    finalVector[16] = 0;
    finalVector[17] = a * atu[17] + c * (atu[3] + atu[10]);
    finalVector[18] = b * (atu[18] + atu[9]);
    finalVector[19] = b * (atu[19] + atu[1]);
    finalVector[20] = ( (a + b) * atu[20] + c * atu[6] ) / 3.;
    finalVector[21] = ( (a + b) * atu[21] + c * atu[14] ) / 3.;
    finalVector[22] = b * ( atu[4] + atu[13] + 2. * atu[22]) / 3.;
    finalVector[23] = (a + 2. * b) * atu[23] / 9.;
}


/* FLOPs: 168A */
__host__ __device__ void au( fvector_t* resVec, const double* u )
{
#ifdef  SINGLE_PRECISION_SOLVER
    int    i, j;
    double finVec[24];
    double temp[24];
#else
    double finVec[24];
#endif

    finVec[0]  = u[0]  - u[1] - u[2] - u[3] + u[4] + u[5] + u[6] - u[7];
    finVec[1]  = u[0]  - u[1] - u[2] + u[3] + u[4] - u[5] - u[6] + u[7];
    finVec[2]  = u[0]  - u[1] + u[2] - u[3] - u[4] + u[5] - u[6] + u[7];
    finVec[3]  = u[0]  - u[1] + u[2] + u[3] - u[4] - u[5] + u[6] - u[7];
    finVec[4]  = u[0]  + u[1] - u[2] - u[3] - u[4] - u[5] + u[6] + u[7];
    finVec[5]  = u[0]  + u[1] - u[2] + u[3] - u[4] + u[5] - u[6] - u[7];
    finVec[6]  = u[0]  + u[1] + u[2] - u[3] + u[4] - u[5] - u[6] - u[7];
    finVec[7]  = u[0]  + u[1] + u[2] + u[3] + u[4] + u[5] + u[6] + u[7];

    finVec[8]  = u[8]  - u[9] - u[10] - u[11] + u[12] + u[13] + u[14] - u[15];
    finVec[9]  = u[8]  - u[9] - u[10] + u[11] + u[12] - u[13] - u[14] + u[15];
    finVec[10] = u[8]  - u[9] + u[10] - u[11] - u[12] + u[13] - u[14] + u[15];
    finVec[11] = u[8]  - u[9] + u[10] + u[11] - u[12] - u[13] + u[14] - u[15];
    finVec[12] = u[8]  + u[9] - u[10] - u[11] - u[12] - u[13] + u[14] + u[15];
    finVec[13] = u[8]  + u[9] - u[10] + u[11] - u[12] + u[13] - u[14] - u[15];
    finVec[14] = u[8]  + u[9] + u[10] - u[11] + u[12] - u[13] - u[14] - u[15];
    finVec[15] = u[8]  + u[9] + u[10] + u[11] + u[12] + u[13] + u[14] + u[15];

    finVec[16] = u[16] - u[17] - u[18] - u[19] + u[20] + u[21] + u[22] - u[23];
    finVec[17] = u[16] - u[17] - u[18] + u[19] + u[20] - u[21] - u[22] + u[23];
    finVec[18] = u[16] - u[17] + u[18] - u[19] - u[20] + u[21] - u[22] + u[23];
    finVec[19] = u[16] - u[17] + u[18] + u[19] - u[20] - u[21] + u[22] - u[23];
    finVec[20] = u[16] + u[17] - u[18] - u[19] - u[20] - u[21] + u[22] + u[23];
    finVec[21] = u[16] + u[17] - u[18] + u[19] - u[20] + u[21] - u[22] - u[23];
    finVec[22] = u[16] + u[17] + u[18] - u[19] + u[20] - u[21] - u[22] - u[23];
    finVec[23] = u[16] + u[17] + u[18] + u[19] + u[20] + u[21] + u[22] + u[23];

#ifdef  SINGLE_PRECISION_SOLVER
    reformF( finVec, temp );

    for (j = 0; j<8; j++)
    {
        for (i = 0; i<3; i++)
        {
            resVec[j].f[i] += temp[j*3 + i];
        }
    }
#else
    reformF( finVec, (double *)resVec );
#endif

}


__host__ __device__ void reformF( const double* u, double* newU )
{
    newU[0]  = u[0];
    newU[1]  = u[8];
    newU[2]  = u[16];
    newU[3]  = u[1];
    newU[4]  = u[9];
    newU[5]  = u[17];
    newU[6]  = u[2];
    newU[7]  = u[10];
    newU[8]  = u[18];
    newU[9]  = u[3];
    newU[10] = u[11];
    newU[11] = u[19];
    newU[12] = u[4];
    newU[13] = u[12];
    newU[14] = u[20];
    newU[15] = u[5];
    newU[16] = u[13];
    newU[17] = u[21];
    newU[18] = u[6];
    newU[19] = u[14];
    newU[20] = u[22];
    newU[21] = u[7];
    newU[22] = u[15];
    newU[23] = u[23];
}


__host__ __device__ void reformU( const double* u, double* newU )
{
    newU[0]  = u[0];
    newU[1]  = u[3];
    newU[2]  = u[6];
    newU[3]  = u[9];
    newU[4]  = u[12];
    newU[5]  = u[15];
    newU[6]  = u[18];
    newU[7]  = u[21];
    newU[8]  = u[1];
    newU[9]  = u[4];
    newU[10] = u[7];
    newU[11] = u[10];
    newU[12] = u[13];
    newU[13] = u[16];
    newU[14] = u[19];
    newU[15] = u[22];
    newU[16] = u[2];
    newU[17] = u[5];
    newU[18] = u[8];
    newU[19] = u[11];
    newU[20] = u[14];
    newU[21] = u[17];
    newU[22] = u[20];
    newU[23] = u[23];
}


/* FLOPs: 21M + 36A */
__host__ __device__ void firstVector_kappa( const double* atu, double* finalVector, double kappa)
{
    finalVector[0] += 0.;
    finalVector[1] += 0.;
    finalVector[2] += 0.;
    finalVector[3] += kappa * (atu[3] + atu[10] + atu[17]);
    finalVector[4] += 0.;
    finalVector[5] += kappa * ( atu[5] + atu[12] ) /3.;
    finalVector[6] += kappa * ( atu[6] + atu[20] ) /3.;
    finalVector[7] += kappa * atu[7] / 9.;

    finalVector[8] += 0.;
    finalVector[9] += 0.;
    finalVector[10] += kappa * ( atu[10] + atu[3] + atu[17]);
    finalVector[11] += 0.;
    finalVector[12] += kappa * ( atu[12] + atu[5] ) / 3.;
    finalVector[13] += 0.;
    finalVector[14] += kappa * ( atu[14] + atu[21] ) /3.;
    finalVector[15] += kappa * atu[15] / 9.;

    finalVector[16] += 0.;
    finalVector[17] += kappa * ( atu[17] + atu[3] + atu[10] );
    finalVector[18] += 0.;
    finalVector[19] += 0.;
    finalVector[20] += kappa * ( atu[20] + atu[6] ) / 3.;
    finalVector[21] += kappa * ( atu[21] + atu[14] ) / 3.;
    finalVector[22] += 0.;
    finalVector[23] += kappa * atu[23] / 9.;
}


/* FLOPs: 60M + 48A */
__host__ __device__ void firstVector_mu( const double* atu, double* finalVector, double b )
{
    finalVector[0] += 0;
    finalVector[1] += b * (atu[19] + atu[1]);
    finalVector[2] += b * (atu[11] + atu[2]);
    finalVector[3] += b * ( 4. * atu[3] - 2. * (atu[10] + atu[17]) ) / 3.;
    finalVector[4] += b * (atu[13] + atu[22] + 2. * atu[4]) / 3.;
    finalVector[5] += b * ( 7. * atu[5] - 2. * atu[12] ) / 9.;
    finalVector[6] += b * ( 7. * atu[6] - 2. * atu[20] ) / 9.;
    finalVector[7] += ( 10. * b * atu[7] ) / 27.;

    finalVector[8]  += 0;
    finalVector[9]  += b * (atu[18] + atu[9]);
    finalVector[10] += b * ( 4. * atu[10] - 2. * (atu[3] + atu[17]) ) / 3.;
    finalVector[11] += b * (atu[11] + atu[2]);
    finalVector[12] += b * ( 7. * atu[12] - 2. * atu[5] ) / 9.;
    finalVector[13] += b * (atu[4] + atu[22] + 2. * atu[13]) / 3.;
    finalVector[14] += b * ( 7. * atu[14] - 2. * atu[21] ) / 9.;
    finalVector[15] += (10. * b * atu[15] ) / 27.;

    finalVector[16] += 0;
    finalVector[17] += b * ( 4. * atu[17] - 2. * (atu[3] + atu[10]) ) / 3.;
    finalVector[18] += b * (atu[18] + atu[9]);
    finalVector[19] += b * (atu[19] + atu[1]);
    finalVector[20] += b * ( 7. * atu[20] - 2. * atu[6] ) / 9.;
    finalVector[21] += b * ( 7. * atu[21] - 2. * atu[14] ) / 9.;
    finalVector[22] += b * ( atu[4] + atu[13] + 2. * atu[22]) / 3.;
    finalVector[23] += (10. * b * atu[23] ) / 27.;
}

/**
#ifndef __CUDA_ARCH__ || (__CUDA_ARCH__ < 600)
__device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
                                          (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed, 
                        __double_as_longlong(val + 
                        __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}
#endif
*/


/**
 * MultAddMatVec: Multiply a 3 x 3 Matrix (M) with a 3 x 1 vector (V1)
 *                and a constant (c). Then add the result to the same
 *                target vector (V2)
 *
 *  V2 = V2 + M * V1 * c
 *
 */
__host__ __device__ void MultAddMatVecGPU( fmatrix_t* M, fvector_t* V1, 
					   double c, fvector_t* V2 )
{
  if (1) {
  V2->f[0] += c*(M->f[0][0] * V1->f[0] + 
		 M->f[0][1] * V1->f[1] + M->f[0][2] * V1->f[2]);

  V2->f[1] += c*(M->f[1][0] * V1->f[0] + 
		 M->f[1][1] * V1->f[1] + M->f[1][2] * V1->f[2]);

  V2->f[2] += c*(M->f[2][0] * V1->f[0] + 
		 M->f[2][1] * V1->f[1] + M->f[2][2] * V1->f[2]);
  }

  if (0) {
  V2->f[0] += (M->f[0][0] * V1->f[0] + 
		 M->f[0][1] * V1->f[1] + M->f[0][2] * V1->f[2]);

  V2->f[1] += (M->f[1][0] * V1->f[0] + 
		 M->f[1][1] * V1->f[1] + M->f[1][2] * V1->f[2]);

  V2->f[2] += (M->f[2][0] * V1->f[0] +
		 M->f[2][1] * V1->f[1] + M->f[2][2] * V1->f[2]);
  }

  return;
}


/* ------------------------------------------------------------------------- */
/*                  C-PML (unsplit ADE) GPU kernels                          */
/* ------------------------------------------------------------------------- */
/*
 * Hex8 1-pt-rule node ordering: node n in [0..7] has corner bits
 *   (n & 1) -> +x face when 1, -x when 0       => s_x(n) = (n & 1) ? +1 : -1
 *   (n & 2) -> +y face when 1, -y when 0       => s_y(n) = ...
 *   (n & 4) -> +z face when 1, -z when 0       => s_z(n) = ...
 * At the element center (xi = eta = zeta = 0):
 *   dN_n/dx_j = s_j(n) / (4h).
 * The (1/(4h)) factor is absorbed into the stored per-element a and
 * 1/kappa coefficients (see solver_init), so the kernels work with the
 * raw signed displacement (or pressure) sums S_j = sum_n s_j(n) u_i^(n).
 * The dt^2 * h^2 / 4 volume/time factor is folded into muCoeff,
 * lambdaCoeff (elastic) and forceCoeff (acoustic).
 *
 * Force sign: the central-difference update in kernelDispCalc reads
 *   u^{n+1} = (force + 2M u^n - M u^{n-1}) / M
 * so "force" already encodes dt^2 * F_stiff; an elastic restoring force
 * F_stiff = -K u gives a negative sign in the assembly below.
 */

__global__ void kernelPMLElasticStep(int32_t lenum,
                                     gpu_data_t* gpuData,
                                     int32_t pmlElasticCount)
{
    int32_t pml_idx = (blockIdx.x * blockDim.x) + threadIdx.x;
    if (pml_idx >= pmlElasticCount) return;

    const int32_t eindex = gpuData->pmlElasticElemMapperDevice[pml_idx];
    const int32_t* lnid_base = &gpuData->lnidArrayDevice[eindex];

    /* Read 8 nodal displacements (u_alpha at each node n in [0..7]). */
    double u[8][3];
    for (int n = 0; n < 8; n++) {
        const int32_t lnid = lnid_base[n * lenum];
        u[n][0] = gpuData->tm1Device[lnid].f[0];
        u[n][1] = gpuData->tm1Device[lnid].f[1];
        u[n][2] = gpuData->tm1Device[lnid].f[2];
    }

    /* Signed sums S[j][alpha] = sum_n s_j(n) * u_alpha^(n).
     * du_alpha/dx_j = S[j][alpha] / (4h); the 1/(4h) factor is absorbed
     * into the per-element a and 1/kappa coefficients. */
    double S[3][3];
    for (int j = 0; j < 3; j++) {
        for (int a = 0; a < 3; a++) {
            double s = 0.0;
            for (int n = 0; n < 8; n++) {
                const double sign = ((n >> j) & 1) ? 1.0 : -1.0;
                s += sign * u[n][a];
            }
            S[j][a] = s;
        }
    }

    /* Per-axis ADE coefficients (a and 1/kappa already pre-divided by 4h). */
    const double aJ[3] = { gpuData->pml_elastic_a_xArrayDevice[pml_idx],
                           gpuData->pml_elastic_a_yArrayDevice[pml_idx],
                           gpuData->pml_elastic_a_zArrayDevice[pml_idx] };
    const double bJ[3] = { gpuData->pml_elastic_b_xArrayDevice[pml_idx],
                           gpuData->pml_elastic_b_yArrayDevice[pml_idx],
                           gpuData->pml_elastic_b_zArrayDevice[pml_idx] };
    const double invKJ[3] = {
        gpuData->pml_elastic_invKappa_xArrayDevice[pml_idx],
        gpuData->pml_elastic_invKappa_yArrayDevice[pml_idx],
        gpuData->pml_elastic_invKappa_zArrayDevice[pml_idx] };

    /* ADE recurrence + modified partial derivative.
     *   psi_new        = b_j * psi_old + a_j * S[j][i]
     *   partial[j][i]  = invK_j * S[j][i] + psi_new
     * Index convention: psi[i][j] holds the j-axis convolution of d/dx_j u_i.
     * partial[j][i] is the CPML-modified du_i/dx_j. */
    double partial[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            const double psi_old =
                (double)gpuData->pml_elastic_psiArrayDevice[i][j][pml_idx];
            const double psi_new = bJ[j] * psi_old + aJ[j] * S[j][i];
            gpuData->pml_elastic_psiArrayDevice[i][j][pml_idx] =
                (solver_float)psi_new;
            partial[j][i] = invKJ[j] * S[j][i] + psi_new;
        }
    }

    /* Strain eps_ab = 1/2 ( du_a/dx_b + du_b/dx_a ) using CPML-modified partials.
     * (partial[j][i] = du_i/dx_j, so du_a/dx_b = partial[b][a].) */
    double eps[3][3];
    double trace = 0.0;
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            eps[a][b] = 0.5 * (partial[b][a] + partial[a][b]);
        }
        trace += eps[a][a];
    }

    /* Pre-scaled stress (sigma * dt^2 * h^2 / 4). */
    const double muCoeff =
        gpuData->pml_elastic_muCoeffArrayDevice[pml_idx];
    const double lambdaCoeff =
        gpuData->pml_elastic_lambdaCoeffArrayDevice[pml_idx];
    double sigma[3][3];
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            sigma[a][b] = 2.0 * muCoeff * eps[a][b];
        }
        sigma[a][a] += lambdaCoeff * trace;
    }

    /* Nodal force: f_alpha^(n) = - sum_beta s_beta(n) * sigma_{alpha,beta}. */
    for (int n = 0; n < 8; n++) {
        const double sx = (n & 1) ? 1.0 : -1.0;
        const double sy = (n & 2) ? 1.0 : -1.0;
        const double sz = (n & 4) ? 1.0 : -1.0;
        const int32_t lnid = lnid_base[n * lenum];
        fvector_t* nodalForce = gpuData->forceDevice + lnid;
        for (int a = 0; a < 3; a++) {
            const double f = -(sx * sigma[a][0] +
                               sy * sigma[a][1] +
                               sz * sigma[a][2]);
            atomicAdd((solver_float*)&(nodalForce->f[a]), (solver_float)f);
        }
    }
}


__global__ void kernelPMLAcousticStep(int32_t lenum,
                                      gpu_data_t* gpuData,
                                      int32_t pmlAcousticCount)
{
    int32_t pml_idx = (blockIdx.x * blockDim.x) + threadIdx.x;
    if (pml_idx >= pmlAcousticCount) return;

    const int32_t eindex = gpuData->pmlAcousticElemMapperDevice[pml_idx];
    const int32_t* lnid_base = &gpuData->lnidArrayDevice[eindex];

    /* Read 8 nodal pressures. */
    double p[8];
    for (int n = 0; n < 8; n++) {
        const int32_t lnid = lnid_base[n * lenum];
        p[n] = (double)gpuData->pressureTm1Device[lnid];
    }

    /* Signed sums S[j] = sum_n s_j(n) * p^(n).
     * dp/dx_j = S[j] / (4h); 1/(4h) absorbed into a and 1/kappa. */
    double S[3];
    for (int j = 0; j < 3; j++) {
        double s = 0.0;
        for (int n = 0; n < 8; n++) {
            const double sign = ((n >> j) & 1) ? 1.0 : -1.0;
            s += sign * p[n];
        }
        S[j] = s;
    }

    const double aJ[3] = { gpuData->pml_acoustic_a_xArrayDevice[pml_idx],
                           gpuData->pml_acoustic_a_yArrayDevice[pml_idx],
                           gpuData->pml_acoustic_a_zArrayDevice[pml_idx] };
    const double bJ[3] = { gpuData->pml_acoustic_b_xArrayDevice[pml_idx],
                           gpuData->pml_acoustic_b_yArrayDevice[pml_idx],
                           gpuData->pml_acoustic_b_zArrayDevice[pml_idx] };
    const double invKJ[3] = {
        gpuData->pml_acoustic_invKappa_xArrayDevice[pml_idx],
        gpuData->pml_acoustic_invKappa_yArrayDevice[pml_idx],
        gpuData->pml_acoustic_invKappa_zArrayDevice[pml_idx] };

    /* ADE update and modified gradient: ~dp/dx_j = invK_j*S[j] + psi_new. */
    double partial[3];
    for (int j = 0; j < 3; j++) {
        const double psi_old =
            (double)gpuData->pml_acoustic_psiArrayDevice[j][pml_idx];
        const double psi_new = bJ[j] * psi_old + aJ[j] * S[j];
        gpuData->pml_acoustic_psiArrayDevice[j][pml_idx] =
            (solver_float)psi_new;
        partial[j] = invKJ[j] * S[j] + psi_new;
    }

    /* Pressure force: f^(n) = - forceCoeff * sum_j s_j(n) * partial[j].
     * forceCoeff = dt^2 * h^2 / (4 rho). */
    const double forceCoeff =
        gpuData->pml_acoustic_forceCoeffArrayDevice[pml_idx];
    for (int n = 0; n < 8; n++) {
        const double sx = (n & 1) ? 1.0 : -1.0;
        const double sy = (n & 2) ? 1.0 : -1.0;
        const double sz = (n & 4) ? 1.0 : -1.0;
        const int32_t lnid = lnid_base[n * lenum];
        const double f = -forceCoeff *
                         (sx * partial[0] + sy * partial[1] + sz * partial[2]);
        atomicAdd((solver_float*)(gpuData->pressureForceDevice + lnid),
                  (solver_float)f);
    }
}
