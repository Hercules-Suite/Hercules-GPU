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

__host__ __device__ static inline double acoustic_node_sign(int node, int axis)
{
  return (node & (1 << axis)) ? 1.0 : -1.0;
}

__device__ static inline void acoustic_scalar_stiffness_sums(
    const solver_float p[8], double sums[7])
{
  for (int k = 0; k < 7; k++) {
    sums[k] = 0.0;
  }

  for (int node = 0; node < 8; node++) {
    const double sx = acoustic_node_sign(node, 0);
    const double sy = acoustic_node_sign(node, 1);
    const double sz = acoustic_node_sign(node, 2);
    const double pv = (double)p[node];

    sums[0] += sx * pv;
    sums[1] += sy * pv;
    sums[2] += sz * pv;
    sums[3] += sx * sy * pv;
    sums[4] += sx * sz * pv;
    sums[5] += sy * sz * pv;
    sums[6] += sx * sy * sz * pv;
  }
}

__device__ static inline void acoustic_scalar_stiffness_delta_sums(
    const solver_float p[8], const solver_float p_prev[8], double sums[7])
{
  for (int k = 0; k < 7; k++) {
    sums[k] = 0.0;
  }

  for (int node = 0; node < 8; node++) {
    const double sx = acoustic_node_sign(node, 0);
    const double sy = acoustic_node_sign(node, 1);
    const double sz = acoustic_node_sign(node, 2);
    const double pv = (double)p[node] - (double)p_prev[node];

    sums[0] += sx * pv;
    sums[1] += sy * pv;
    sums[2] += sz * pv;
    sums[3] += sx * sy * pv;
    sums[4] += sx * sz * pv;
    sums[5] += sy * sz * pv;
    sums[6] += sx * sy * sz * pv;
  }
}

__device__ static inline double acoustic_scalar_stiffness_from_sums(
    int node, const double sums[7])
{
  const double sx = acoustic_node_sign(node, 0);
  const double sy = acoustic_node_sign(node, 1);
  const double sz = acoustic_node_sign(node, 2);

  return sx * sums[0] + sy * sums[1] + sz * sums[2] +
         (2.0 / 3.0) * sx * sy * sums[3] +
         (2.0 / 3.0) * sx * sz * sums[4] +
         (2.0 / 3.0) * sy * sz * sums[5] +
         (1.0 / 3.0) * sx * sy * sz * sums[6];
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

__host__ __device__ static inline int acoustic_interface_face_axes_from_face(
    uint8_t face, int tangent_axes[2])
{
  int normal_axis = 2;

  if ((face == ACOUSTIC_FACE_X_MINUS) || (face == ACOUSTIC_FACE_X_PLUS)) {
    normal_axis = 0;
  } else if ((face == ACOUSTIC_FACE_Y_MINUS) ||
             (face == ACOUSTIC_FACE_Y_PLUS)) {
    normal_axis = 1;
  }

  int cursor = 0;
  for (int dim = 0; dim < 3; dim++) {
    if (dim != normal_axis) {
      tangent_axes[cursor++] = dim;
    }
  }

  return normal_axis;
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

    if (attributes.maxThreadsPerBlock > 0) {
      max_threads = imin(max_threads, attributes.maxThreadsPerBlock);
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

    /* PML elements are integrated completely by their own modal kernel. */
    if (gpuData->pmlFlagArrayDevice != NULL &&
        gpuData->pmlFlagArrayDevice[eindex] != 0) return;

    int32_t *lnid = &(gpuData->lnidArrayDevice[eindex]);

    /* Rayleigh stiffness-proportional damping b*K, ported to the main
     * GPU stiffness path. It was previously only in the CPU conventional
     * path (stiffness.cu) and the GPU topo kernel -- absent here, so
     * type_of_damping=rayleigh GPU (non-topo) runs silently lost b*K.
     * Effective displacement  u_eff = (1 + b/dt) u^n - (b/dt) u^{n-1},
     * b/dt = c3/c1, fed to the SAME operator => K u^n + b K u_dot.
     * c3 == 0 for every non-Rayleigh damping (compute_setab zeroes bBase
     * for NONE/BKT/MASS) => sd_b_over_dt == 0 => u_eff == tm1 EXACTLY
     * (the else branch is the byte-identical original) => bit-identical
     * for all non-Rayleigh runs incl. every BKT/none CPML validation.
     * PML elements have already returned and use their own complete operator. */
    double sd_b_over_dt = 0.0;
    {
      const double c1v = gpuData->c1ArrayDevice[eindex];
      if (c1v != 0.0) sd_b_over_dt = gpuData->c3ArrayDevice[eindex] / c1v;

    }

    for (i = 0; i < 8; i++) {
      const __restrict__ fvector_t *tm1Disp = gpuData->tm1Device + *lnid;

      if (sd_b_over_dt != 0.0) {
        const __restrict__ fvector_t *tm2Disp = gpuData->tm2Device + *lnid;
        localForce[i].f[0] = tm1Disp->f[0] * (1.0 + sd_b_over_dt) -
                             sd_b_over_dt * tm2Disp->f[0];
        localForce[i].f[1] = tm1Disp->f[1] * (1.0 + sd_b_over_dt) -
                             sd_b_over_dt * tm2Disp->f[1];
        localForce[i].f[2] = tm1Disp->f[2] * (1.0 + sd_b_over_dt) -
                             sd_b_over_dt * tm2Disp->f[2];
      } else {
        localForce[i].f[0] = tm1Disp->f[0];
        localForce[i].f[1] = tm1Disp->f[1];
        localForce[i].f[2] = tm1Disp->f[2];
      }
      do_contrib = ( fabs(localForce[i].f[0]) > UNDERFLOW_CAP_STIFFNESS ) ? 1 : do_contrib;
      do_contrib = ( fabs(localForce[i].f[1]) > UNDERFLOW_CAP_STIFFNESS ) ? 1 : do_contrib;
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

/* Stretch the interface flux by the cofactor J/S_j. */
__device__ static inline void acoustic_cpml_interface_flux(gpu_data_t* gpu,
    int face, int q, int side, double value[3]) {
  if (gpu->acousticInterfacePmlSideMaskDevice == NULL ||
      !(gpu->acousticInterfacePmlSideMaskDevice[face] & side)) return;
  const size_t stride = (size_t)gpu->acousticInterfaceFaceCount;
  const int solid = side == ACOUSTIC_INTERFACE_PML_ELASTIC_SIDE;
  const cpml_interface_t* coeff = solid ? gpu->acousticInterfacePmlElasticDevice
                                       : gpu->acousticInterfacePmlAcousticDevice;
  solver_float* state = solid ? gpu->acousticInterfacePsiPressureDevice
                              : gpu->acousticInterfacePsiAccelDevice;
  cpml_interface_apply(coeff[face], value, state+face+6*q*stride, stride);
}

__global__ void kernelAcousticPressureDisplacementCoupling(int32_t faceCount,
                                                           gpu_data_t* gpuData,
                                                           double deltaT2,
                                                           double deltaT)
{
  const int vp = (gpuData->acousticVelocityPotential == 1);
  const double invDt = (deltaT > 0.0) ? (1.0 / deltaT) : 0.0;
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
  const uint8_t face =
      (gpuData->acousticInterfaceFaceIdsDevice != NULL)
          ? gpuData->acousticInterfaceFaceIdsDevice[faceIndex]
          : ACOUSTIC_FACE_Z_PLUS;
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
  const int32_t* solid_nodes = gpuData->acousticInterfaceSolidNodesDevice;
  const bool solid_trace = (solid_nodes != NULL) || (adjacent_eindex >= 0);
  const int32_t force_eindex = (adjacent_eindex >= 0) ? adjacent_eindex : eindex;
  const double* force_center = solid_trace ? adjacent_center : center;

  if ((gpuData->acousticInterfaceProjectionMode ==
       ACOUSTIC_INTERFACE_PROJECTION_FACE_MASS) &&
      solid_trace &&
      (gpuData->acousticInterfaceFaceHalfWidthsDevice != NULL) &&
      (gpuData->acousticInterfaceAdjacentHalfWidthsDevice != NULL)) {
    const double* half =
        gpuData->acousticInterfaceFaceHalfWidthsDevice + 3 * faceIndex;
    const double* adjacent_half =
        gpuData->acousticInterfaceAdjacentHalfWidthsDevice + 3 * faceIndex;
    const double gauss[2] = {-0.57735026918962576451,
                              0.57735026918962576451};
    const double q_coeff = face_coeff * 0.25;
    int tangent_axes[2];

    (void)acoustic_interface_face_axes_from_face(face, tangent_axes);

    for (int g0 = 0; g0 < 2; g0++) {
      for (int g1 = 0; g1 < 2; g1++) {
        double q_center[3] = {center[0], center[1], center[2]};
        double q_force_center[3] = {adjacent_center[0], adjacent_center[1],
                                    adjacent_center[2]};
        double q_pressure = 0.0;

        q_center[tangent_axes[0]] += half[tangent_axes[0]] * gauss[g0];
        q_center[tangent_axes[1]] += half[tangent_axes[1]] * gauss[g1];
        q_force_center[tangent_axes[0]] +=
            adjacent_half[tangent_axes[0]] * gauss[g0];
        q_force_center[tangent_axes[1]] +=
            adjacent_half[tangent_axes[1]] * gauss[g1];

        for (int i = 0; i < 8; i++) {
          const int32_t lnid =
              gpuData->lnidArrayDevice[eindex + i * gpuData->lenum];
          const double weight =
              acoustic_shape_value_at(i, q_center[0], q_center[1],
                                      q_center[2]);

          if (vp) {
            q_pressure +=
                weight *
                (-((double)gpuData->pressureTm1Device[lnid] -
                   (double)gpuData->pressureTm2Device[lnid]) * invDt);
          } else {
            q_pressure +=
                weight * (double)gpuData->pressureTm1Device[lnid];
          }
        }

        double traction[3] = {q_pressure*normal[0], q_pressure*normal[1],
                              q_pressure*normal[2]};
        acoustic_cpml_interface_flux(gpuData, faceIndex, 2*g0+g1,
                                     ACOUSTIC_INTERFACE_PML_ELASTIC_SIDE, traction);

        for (int i = 0; i < 8; i++) {
          const double weight =
              acoustic_shape_value_at(i, q_force_center[0],
                                      q_force_center[1], q_force_center[2]);
          if (weight == 0.0) {
            continue;
          }

          const int32_t lnid =
              solid_nodes ? solid_nodes[8 * faceIndex + i] :
              gpuData->lnidArrayDevice[force_eindex + i * gpuData->lenum];
          fvector_t* nodalForce = gpuData->forceDevice + lnid;
          atomicAdd((solver_float*)&nodalForce->f[0],
                    (solver_float)(q_coeff * weight * traction[0]));
          atomicAdd((solver_float*)&nodalForce->f[1],
                    (solver_float)(q_coeff * weight * traction[1]));
          atomicAdd((solver_float*)&nodalForce->f[2],
                    (solver_float)(q_coeff * weight * traction[2]));
        }
      }
    }
    return;
  }

  double pressure = 0.0;

  for (int i = 0; i < 8; i++) {
    const int32_t lnid =
        gpuData->lnidArrayDevice[eindex + i * gpuData->lenum];
    const double weight =
        acoustic_shape_value_at(i, center[0], center[1], center[2]);

    /* velocity_potential: the solid feels p = -dphi/dt
     * = -(phi^n - phi^{n-1})/dt. pressure mode: field IS p. */
    if (vp) {
      pressure += weight *
                  (-((double)gpuData->pressureTm1Device[lnid] -
                     (double)gpuData->pressureTm2Device[lnid]) * invDt);
    } else {
      pressure += weight * (double)gpuData->pressureTm1Device[lnid];
    }
  }

  double traction[3] = {pressure*normal[0], pressure*normal[1], pressure*normal[2]};
  acoustic_cpml_interface_flux(gpuData, faceIndex, 0,
                               ACOUSTIC_INTERFACE_PML_ELASTIC_SIDE, traction);

  for (int i = 0; i < 8; i++) {
    const double weight =
        acoustic_shape_value_at(i, force_center[0], force_center[1],
                                force_center[2]);
    if (weight == 0.0) {
      continue;
    }

    const int32_t lnid =
        solid_nodes ? solid_nodes[8 * faceIndex + i] :
        gpuData->lnidArrayDevice[force_eindex + i * gpuData->lenum];
    fvector_t* nodalForce = gpuData->forceDevice + lnid;
    atomicAdd((solver_float*)&nodalForce->f[0],
              (solver_float)(face_coeff * weight * traction[0]));
    atomicAdd((solver_float*)&nodalForce->f[1],
              (solver_float)(face_coeff * weight * traction[1]));
    atomicAdd((solver_float*)&nodalForce->f[2],
              (solver_float)(face_coeff * weight * traction[2]));
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

  /* The PML kernel supplies the complete selected element operator. */
  if (gpuData->pmlFlagArrayDevice != NULL &&
      gpuData->pmlFlagArrayDevice[eindex] != 0) return;

  for (int i = 0; i < 8; i++) {
    const int32_t lnid = gpuData->lnidArrayDevice[eindex + i * lenum];
    pressure[i] = gpuData->pressureTm1Device[lnid];
  }

  const double stiffness_coeff =
      gpuData->acousticPressureCoeffArrayDevice[eindex];

  /* Optional stiffness-proportional (Kelvin-Voigt) acoustic damping.
   * Adds -coeff*beta*K*(p^n - p^{n-1}) to the acoustic force: the same
   * sign-mode stiffness operator K applied to the per-step field increment,
   * so the damping rate ~ omega^2
   * (Q ~ 1/omega: high-frequency selective -- attenuates the lossless-
   * cavity c/4h grid ring while barely touching the source-band
   * signal). zeta(w) = beta*w*dt/2, so beta = 1/(Q*2*pi*theFreq*dt)
   * (host-side, damping-type independent; see solver_init). When
   * beta == 0 (Q unset/<=0), the lossless path still uses the same
   * algebraically equivalent sign-mode stiffness evaluation below, avoiding
   * the older dense 8x8 loop. Skipped on CPML PML elements so the PML zone
   * stays pure CPML (consistent with the material-damping skip in commit
   * 5d91b09). Non-SPECFEM numerical regularizer; see README. */
  const double damp_beta = gpuData->acousticStiffnessDampingBeta;
  const int damp_on =
      (damp_beta != 0.0) &&
      !(gpuData->pmlFlagArrayDevice != NULL &&
        gpuData->pmlFlagArrayDevice[eindex] != 0);
  solver_float pressure_prev[8];
  if (damp_on) {
    for (int i = 0; i < 8; i++) {
      const int32_t lnid = gpuData->lnidArrayDevice[eindex + i * lenum];
      pressure_prev[i] = gpuData->pressureTm2Device[lnid];
    }

    /* Fast algebraically equivalent application of the same consistent
     * acoustic stiffness matrix.  The shape matrix is diagonal in the
     * seven non-constant trilinear sign modes (x, y, z, xy, xz, yz,
     * xyz), so Q-active damping can evaluate K*p and K*(p-p_prev)
     * without two dense 8x8 loops. */
    double pressure_sums[7];
    double delta_sums[7];
    acoustic_scalar_stiffness_sums(pressure, pressure_sums);
    acoustic_scalar_stiffness_delta_sums(pressure, pressure_prev, delta_sums);

    for (int i = 0; i < 8; i++) {
      const int32_t lnid = gpuData->lnidArrayDevice[eindex + i * lenum];
      const double local =
          acoustic_scalar_stiffness_from_sums(i, pressure_sums);
      const double dlocal =
          acoustic_scalar_stiffness_from_sums(i, delta_sums);
      atomicAdd((solver_float*)(gpuData->pressureForceDevice + lnid),
                (solver_float)(-stiffness_coeff * local -
                               stiffness_coeff * damp_beta * dlocal));
    }
    return;
  }

  double pressure_sums[7];
  acoustic_scalar_stiffness_sums(pressure, pressure_sums);

  for (int i = 0; i < 8; i++) {
    const int32_t lnid = gpuData->lnidArrayDevice[eindex + i * lenum];
    const double local = acoustic_scalar_stiffness_from_sums(i, pressure_sums);

    atomicAdd((solver_float*)(gpuData->pressureForceDevice + lnid),
              (solver_float)(-stiffness_coeff * local));
  }
}

__global__ void kernelAcousticPressureInterfaceForce(int32_t faceCount,
                                                     gpu_data_t* gpuData,
                                                     double deltaT2,
                                                     double deltaT)
{
  int32_t faceIndex = (blockIdx.x * blockDim.x) + threadIdx.x;
  const double invDeltaT2 = 1.0 / deltaT2;
  const int vp = (gpuData->acousticVelocityPotential == 1);
  const double invDt = (deltaT > 0.0) ? (1.0 / deltaT) : 0.0;

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
  const uint8_t face =
      (gpuData->acousticInterfaceFaceIdsDevice != NULL)
          ? gpuData->acousticInterfaceFaceIdsDevice[faceIndex]
          : ACOUSTIC_FACE_Z_PLUS;
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
  const int32_t* solid_nodes = gpuData->acousticInterfaceSolidNodesDevice;
  const bool solid_trace = (solid_nodes != NULL) || (adjacent_eindex >= 0);
  const double* accel_center =
      solid_trace ? adjacent_center : center;

  if ((gpuData->acousticInterfaceProjectionMode ==
       ACOUSTIC_INTERFACE_PROJECTION_FACE_MASS) &&
      solid_trace &&
      (gpuData->acousticInterfaceFaceHalfWidthsDevice != NULL) &&
      (gpuData->acousticInterfaceAdjacentHalfWidthsDevice != NULL)) {
    const double* half =
        gpuData->acousticInterfaceFaceHalfWidthsDevice + 3 * faceIndex;
    const double* adjacent_half =
        gpuData->acousticInterfaceAdjacentHalfWidthsDevice + 3 * faceIndex;
    const double gauss[2] = {-0.57735026918962576451,
                              0.57735026918962576451};
    const double sgn = vp ? 1.0 : -1.0;
    const double q_coeff = face_coeff * deltaT2 * 0.25;
    int tangent_axes[2];

    (void)acoustic_interface_face_axes_from_face(face, tangent_axes);

    for (int g0 = 0; g0 < 2; g0++) {
      for (int g1 = 0; g1 < 2; g1++) {
        double q_center[3] = {center[0], center[1], center[2]};
        double q_accel_center[3] = {adjacent_center[0], adjacent_center[1],
                                    adjacent_center[2]};
        double q_source[3] = {0,0,0};

        q_center[tangent_axes[0]] += half[tangent_axes[0]] * gauss[g0];
        q_center[tangent_axes[1]] += half[tangent_axes[1]] * gauss[g1];
        q_accel_center[tangent_axes[0]] +=
            adjacent_half[tangent_axes[0]] * gauss[g0];
        q_accel_center[tangent_axes[1]] +=
            adjacent_half[tangent_axes[1]] * gauss[g1];

        for (int i = 0; i < 8; i++) {
          const double weight =
              acoustic_shape_value_at(i, q_accel_center[0],
                                      q_accel_center[1],
                                      q_accel_center[2]);
          if (weight == 0.0) {
            continue;
          }

          const int32_t lnid =
              solid_nodes ? solid_nodes[8 * faceIndex + i] :
              gpuData->lnidArrayDevice[accel_eindex + i * gpuData->lenum];
          const fvector_t u_next = gpuData->tm2Device[lnid];
          const fvector_t u_curr = gpuData->tm1Device[lnid];
          const fvector_t u_prev = vp ? u_curr : gpuData->tm3Device[lnid];
          for (int j=0; j<3; ++j) {
            const double motion = vp ? ((double)u_next.f[j]-u_curr.f[j])*invDt
                : ((double)u_next.f[j]-2.0*u_curr.f[j]+u_prev.f[j])*invDeltaT2;
            q_source[j] += weight*motion*normal[j];
          }
        }
        acoustic_cpml_interface_flux(gpuData, faceIndex, 2*g0+g1,
                                     ACOUSTIC_INTERFACE_PML_ACOUSTIC_SIDE, q_source);
        const double q_source_n = q_source[0]+q_source[1]+q_source[2];

        for (int i = 0; i < 8; i++) {
          const double weight =
              acoustic_shape_value_at(i, q_center[0], q_center[1],
                                      q_center[2]);
          if (weight == 0.0) {
            continue;
          }

          const int32_t lnid =
              gpuData->lnidArrayDevice[eindex + i * gpuData->lenum];
          atomicAdd((solver_float*)(gpuData->pressureForceDevice + lnid),
                    (solver_float)(sgn * q_coeff * weight * q_source_n));
        }
      }
    }
    return;
  }

  /* source_n is the scalar exchanged into the acoustic RHS:
   *  - pressure mode: solid normal ACCELERATION (2nd time-difference,
   *    ~omega^2) -- the validated path (and the grid pump);
   *  - velocity_potential: solid normal VELOCITY (uⁿ⁺¹-uⁿ)/dt . n
   *    (1st time-difference, ~omega) -- the milder, stable coupling. */
  double source[3] = {0,0,0};

  for (int i = 0; i < 8; i++) {
    const double weight =
        acoustic_shape_value_at(i, accel_center[0], accel_center[1],
                                accel_center[2]);
    if (weight == 0.0) {
      continue;
    }

    const int32_t lnid =
        solid_nodes ? solid_nodes[8 * faceIndex + i] :
        gpuData->lnidArrayDevice[accel_eindex + i * gpuData->lenum];
    const fvector_t u_next = gpuData->tm2Device[lnid];
    const fvector_t u_curr = gpuData->tm1Device[lnid];
    const fvector_t u_prev = vp ? u_curr : gpuData->tm3Device[lnid];
    for (int j=0; j<3; ++j) {
      const double motion = vp ? ((double)u_next.f[j]-u_curr.f[j])*invDt
          : ((double)u_next.f[j]-2.0*u_curr.f[j]+u_prev.f[j])*invDeltaT2;
      source[j] += weight*motion*normal[j];
    }
  }
  acoustic_cpml_interface_flux(gpuData, faceIndex, 0,
                               ACOUSTIC_INTERFACE_PML_ACOUSTIC_SIDE, source);
  const double source_n = source[0]+source[1]+source[2];

  /* pressure: acoustic RHS += -C u_s_ddot ; velocity_potential:
   * phi-RHS += +G v_s (sign flip; same |coeff|). */
  const double sgn = vp ? 1.0 : -1.0;
  for (int i = 0; i < 8; i++) {
    const double weight =
        acoustic_shape_value_at(i, center[0], center[1], center[2]);
    if (weight == 0.0) {
      continue;
    }

    const int32_t lnid =
        gpuData->lnidArrayDevice[eindex + i * gpuData->lenum];
    atomicAdd((solver_float*)(gpuData->pressureForceDevice + lnid),
              (solver_float)(sgn * face_coeff * weight * deltaT2 * source_n));
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
      (gpuData->pressureFreeSurfaceNodeFlagsDevice[nindex] != 0) &&
      (gpuData->acousticFreeSurfaceGravity == 0)) {
    /* Pressure-release free surface, p = 0. Pressure mode: field=p, set
     * 0. velocity_potential: p = -dphi/dt = 0  ->  phi^{n+1} = phi^n. */
    gpuData->pressureTm2Device[nindex] =
        (gpuData->acousticVelocityPotential == 1)
            ? gpuData->pressureTm1Device[nindex]
            : (solver_float)0.0;
    gpuData->pressureForceDevice[nindex] = 0.0;
    return;
  }

  /* Gravity free-surface nodes carry their weak boundary contribution in
   * pressureMassDevice and therefore use the ordinary scalar update below. */
  const solver_float mass = gpuData->pressureMassDevice[nindex];
  /* Step B: acoustic CPML mass-proportional damping. cd is the
   * assembled dt*C term (dt baked in at build, dt*C summed across
   * elements and ranks like pressure_mass). cd == 0 everywhere unless
   * an acoustic PML element contributes, in which case the update is
   * the damped central difference
   *   mass*(p^{n+1}-2p^n+p^{n-1}) + dt*C*(p^n - p^{n-1}) = F
   * mirroring the validated elastic v9 dashpot discretization. When
   * cd == 0 this reduces EXACTLY to the legacy update (bit-identical
   * for any run without acoustic PML). */
  const solver_float cd = (gpuData->pressureDampDevice != NULL)
                              ? gpuData->pressureDampDevice[nindex]
                              : (solver_float)0.0;
  solver_float p_new = 0.0;

  if (mass > 0.0) {
    p_new = (gpuData->pressureForceDevice[nindex] +
             ((solver_float)2.0 * mass - cd) *
                 gpuData->pressureTm1Device[nindex] -
             (mass - cd) * gpuData->pressureTm2Device[nindex]) /
            mass;
  }

  /* CPML outer-face Dirichlet clamp -- acoustic field. The clamp mask
   * is geometric (any node on an x_min/x_max/y_min/y_max/z_max domain-
   * extent face), so acoustic-PML outer-face nodes are already flagged;
   * kernelDispCalc zeros u there. Here the clamped scalar is pressure or
   * density-scaled velocity potential, depending on the formulation.
   * Absorption is provided by the CPML layer; this outer termination can
   * reflect residual waves. No-op when cpml_clamp_outer=no (mask NULL).
   * Pressure-release top nodes return earlier; gravity top nodes reach here. */
  if (gpuData->pmlClampNodeMaskDevice != NULL &&
      gpuData->pmlClampNodeMaskDevice[nindex] != 0) {
    p_new = 0.0;
  }

  gpuData->pressureTm2Device[nindex] = p_new;
  gpuData->pressureForceDevice[nindex] = 0.0;
}

__device__ static double poro_node_sign_device(int node, int axis)
{
  return (node & (1 << axis)) ? 1.0 : -1.0;
}

__global__ void kernelPoroelasticPressureTraction(int32_t elemCount,
                                                  gpu_data_t* gpuData,
                                                  double deltaT2,
                                                  double biotAlpha)
{
  const int32_t pindex = (blockIdx.x * blockDim.x) + threadIdx.x;
  const double gp = 0.57735026918962576451;
  const double gcoord[2] = {-gp, gp};

  if (pindex >= elemCount) {
    return;
  }

  const int32_t eindex = gpuData->poroelasticElemMapperDevice[pindex];
  const double h = gpuData->poroelasticElemHDevice[pindex];
  const double detJ = gpuData->poroelasticElemDetJDevice[pindex];
  if ((h <= 0.0) || (detJ <= 0.0)) {
    return;
  }

  double qForce[8][3];
  for (int j = 0; j < 8; j++) {
    qForce[j][0] = 0.0;
    qForce[j][1] = 0.0;
    qForce[j][2] = 0.0;
  }

  for (int gx = 0; gx < 2; gx++) {
    const double xi = gcoord[gx];
    for (int gy = 0; gy < 2; gy++) {
      const double eta = gcoord[gy];
      for (int gz = 0; gz < 2; gz++) {
        const double zeta = gcoord[gz];
        double dNdx[8][3];
        double pressureAtGp = 0.0;

        for (int j = 0; j < 8; j++) {
          const double sx = poro_node_sign_device(j, 0);
          const double sy = poro_node_sign_device(j, 1);
          const double sz = poro_node_sign_device(j, 2);
          const double ax = 1.0 + sx * xi;
          const double ay = 1.0 + sy * eta;
          const double az = 1.0 + sz * zeta;
          const int32_t lnid =
              gpuData->lnidArrayDevice[eindex + j * gpuData->lenum];
          const double N = 0.125 * ax * ay * az;
          const double shapeScale =
              gpuData->poroelasticPressureShapeScaleDevice[j * elemCount +
                                                           pindex];
          const double Np = N * shapeScale;

          dNdx[j][0] = sx * ay * az / (4.0 * h);
          dNdx[j][1] = sy * ax * az / (4.0 * h);
          dNdx[j][2] = sz * ax * ay / (4.0 * h);
          pressureAtGp +=
              Np * (double)gpuData->porePressureTm1Device[lnid];
        }

        const double scale = biotAlpha * pressureAtGp * detJ;
        for (int j = 0; j < 8; j++) {
          qForce[j][0] += scale * dNdx[j][0];
          qForce[j][1] += scale * dNdx[j][1];
          qForce[j][2] += scale * dNdx[j][2];
        }
      }
    }
  }

  for (int j = 0; j < 8; j++) {
    const int32_t lnid =
        gpuData->lnidArrayDevice[eindex + j * gpuData->lenum];
    fvector_t* nodalForce = gpuData->forceDevice + lnid;
    atomicAdd((solver_float*)&nodalForce->f[0],
              (solver_float)(deltaT2 * qForce[j][0]));
    atomicAdd((solver_float*)&nodalForce->f[1],
              (solver_float)(deltaT2 * qForce[j][1]));
    atomicAdd((solver_float*)&nodalForce->f[2],
              (solver_float)(deltaT2 * qForce[j][2]));
  }
}

__global__ void kernelPoroelasticPressureForce(int32_t elemCount,
                                               gpu_data_t* gpuData,
                                               double deltaT,
                                               double mobilityX,
                                               double mobilityY,
                                               double mobilityZ,
                                               int haveDiffusion,
                                               double biotAlpha)
{
  const int32_t pindex = (blockIdx.x * blockDim.x) + threadIdx.x;
  const double gp = 0.57735026918962576451;
  const double gcoord[2] = {-gp, gp};

  if ((pindex >= elemCount) || ((biotAlpha == 0.0) && !haveDiffusion)) {
    return;
  }

  const int32_t eindex = gpuData->poroelasticElemMapperDevice[pindex];
  const double h = gpuData->poroelasticElemHDevice[pindex];
  const double detJ = gpuData->poroelasticElemDetJDevice[pindex];
  const double storageM = gpuData->poroelasticStorageMDevice[pindex];
  if ((h <= 0.0) || (detJ <= 0.0) || (storageM <= 0.0)) {
    return;
  }

  double qtDeltaU[8];
  for (int j = 0; j < 8; j++) {
    qtDeltaU[j] = 0.0;
  }

  for (int gx = 0; gx < 2; gx++) {
    const double xi = gcoord[gx];
    for (int gy = 0; gy < 2; gy++) {
      const double eta = gcoord[gy];
      for (int gz = 0; gz < 2; gz++) {
        const double zeta = gcoord[gz];
        double Np[8];
        double divDeltaU = 0.0;

        for (int j = 0; j < 8; j++) {
          const double sx = poro_node_sign_device(j, 0);
          const double sy = poro_node_sign_device(j, 1);
          const double sz = poro_node_sign_device(j, 2);
          const double ax = 1.0 + sx * xi;
          const double ay = 1.0 + sy * eta;
          const double az = 1.0 + sz * zeta;
          const int32_t lnid =
              gpuData->lnidArrayDevice[eindex + j * gpuData->lenum];
          const fvector_t unew = gpuData->tm2Device[lnid];
          const fvector_t uold = gpuData->tm1Device[lnid];
          const double N = 0.125 * ax * ay * az;
          const double shapeScale =
              gpuData->poroelasticPressureShapeScaleDevice[j * elemCount +
                                                           pindex];
          const double dNdx0 = sx * ay * az / (4.0 * h);
          const double dNdx1 = sy * ax * az / (4.0 * h);
          const double dNdx2 = sz * ax * ay / (4.0 * h);

          Np[j] = N * shapeScale;
          divDeltaU += dNdx0 * ((double)unew.f[0] - (double)uold.f[0]);
          divDeltaU += dNdx1 * ((double)unew.f[1] - (double)uold.f[1]);
          divDeltaU += dNdx2 * ((double)unew.f[2] - (double)uold.f[2]);
        }

        const double scale = biotAlpha * divDeltaU * detJ;
        for (int j = 0; j < 8; j++) {
          qtDeltaU[j] += scale * Np[j];
        }
      }
    }
  }

  if (haveDiffusion) {
    double diffusionWeak[8];
    const double mobility[3] = {mobilityX, mobilityY, mobilityZ};

    for (int j = 0; j < 8; j++) {
      diffusionWeak[j] = 0.0;
    }

    for (int gx = 0; gx < 2; gx++) {
      const double xi = gcoord[gx];
      for (int gy = 0; gy < 2; gy++) {
        const double eta = gcoord[gy];
        for (int gz = 0; gz < 2; gz++) {
          const double zeta = gcoord[gz];
          double gradp[3] = {0.0, 0.0, 0.0};

          for (int j = 0; j < 8; j++) {
            const double sx = poro_node_sign_device(j, 0);
            const double sy = poro_node_sign_device(j, 1);
            const double sz = poro_node_sign_device(j, 2);
            const double ax = 1.0 + sx * xi;
            const double ay = 1.0 + sy * eta;
            const double az = 1.0 + sz * zeta;
            const int32_t lnid =
                gpuData->lnidArrayDevice[eindex + j * gpuData->lenum];
            const double pold = (double)gpuData->porePressureTm1Device[lnid];
            const double shapeScale =
                gpuData->poroelasticPressureShapeScaleDevice[j * elemCount +
                                                             pindex];
            const double dNpdx0 = shapeScale * sx * ay * az / (4.0 * h);
            const double dNpdx1 = shapeScale * sy * ax * az / (4.0 * h);
            const double dNpdx2 = shapeScale * sz * ax * ay / (4.0 * h);

            gradp[0] += dNpdx0 * pold;
            gradp[1] += dNpdx1 * pold;
            gradp[2] += dNpdx2 * pold;
          }

          for (int j = 0; j < 8; j++) {
            const double sx = poro_node_sign_device(j, 0);
            const double sy = poro_node_sign_device(j, 1);
            const double sz = poro_node_sign_device(j, 2);
            const double ax = 1.0 + sx * xi;
            const double ay = 1.0 + sy * eta;
            const double az = 1.0 + sz * zeta;
            const double shapeScale =
                gpuData->poroelasticPressureShapeScaleDevice[j * elemCount +
                                                             pindex];
            const double dNpdx0 = shapeScale * sx * ay * az / (4.0 * h);
            const double dNpdx1 = shapeScale * sy * ax * az / (4.0 * h);
            const double dNpdx2 = shapeScale * sz * ax * ay / (4.0 * h);
            double weakDiffusion = 0.0;

            weakDiffusion += mobility[0] * gradp[0] * dNpdx0;
            weakDiffusion += mobility[1] * gradp[1] * dNpdx1;
            weakDiffusion += mobility[2] * gradp[2] * dNpdx2;
            diffusionWeak[j] += detJ * weakDiffusion;
          }
        }
      }
    }

    for (int j = 0; j < 8; j++) {
      const int32_t lnid =
          gpuData->lnidArrayDevice[eindex + j * gpuData->lenum];
      double rhs = -qtDeltaU[j] - deltaT * diffusionWeak[j];

      atomicAdd((solver_float*)(gpuData->porePressureForceDevice + lnid),
                (solver_float)rhs);
    }
    return;
  }

  for (int j = 0; j < 8; j++) {
    const int32_t lnid =
        gpuData->lnidArrayDevice[eindex + j * gpuData->lenum];
    const double rhs = -qtDeltaU[j];
    atomicAdd((solver_float*)(gpuData->porePressureForceDevice + lnid),
              (solver_float)rhs);
  }
}

__global__ void kernelAcousticPoroelasticOpenPoreTarget(int32_t faceCount,
                                                        gpu_data_t* gpuData,
                                                        double deltaT)
{
  const int32_t faceIndex = (blockIdx.x * blockDim.x) + threadIdx.x;
  const int vp = (gpuData->acousticVelocityPotential == 1);
  const double invDt = (deltaT > 0.0) ? (1.0 / deltaT) : 0.0;

  if ((faceIndex >= faceCount) ||
      (gpuData->acousticInterfaceAdjacentPoroelasticDevice == NULL) ||
      (gpuData->acousticInterfacePoroelasticBoundaryScalesDevice == NULL) ||
      (gpuData->porePressureOpenPoreTargetDevice == NULL)) {
    return;
  }
  if (gpuData->acousticInterfaceAdjacentPoroelasticDevice[faceIndex] == 0) {
    return;
  }

  const double scale =
      gpuData->acousticInterfacePoroelasticBoundaryScalesDevice[faceIndex];
  if (scale <= 0.0) {
    return;
  }

  const int32_t eindex =
      gpuData->acousticInterfaceElementIdsDevice[faceIndex];
  const int32_t adjacent_eindex =
      (gpuData->acousticInterfaceAdjacentElementIndexesDevice != NULL)
          ? gpuData->acousticInterfaceAdjacentElementIndexesDevice[faceIndex]
          : -1;
  const int32_t target_eindex =
      (adjacent_eindex >= 0) ? adjacent_eindex : eindex;
  const int32_t* solid_nodes = gpuData->acousticInterfaceSolidNodesDevice;
  const double center[3] = {
      gpuData->acousticInterfaceFaceCentersDevice[3 * faceIndex + 0],
      gpuData->acousticInterfaceFaceCentersDevice[3 * faceIndex + 1],
      gpuData->acousticInterfaceFaceCentersDevice[3 * faceIndex + 2]};
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
  const double* target_center =
      (solid_nodes != NULL || adjacent_eindex >= 0) ? adjacent_center : center;
  double pressure = 0.0;

  for (int i = 0; i < 8; i++) {
    const int32_t lnid =
        gpuData->lnidArrayDevice[eindex + i * gpuData->lenum];
    const double weight =
        acoustic_shape_value_at(i, center[0], center[1], center[2]);

    if (vp) {
      pressure +=
          weight *
          (-((double)gpuData->pressureTm2Device[lnid] -
             (double)gpuData->pressureTm1Device[lnid]) * invDt);
    } else {
      pressure += weight * (double)gpuData->pressureTm2Device[lnid];
    }
  }

  for (int j = 0; j < 8; j++) {
    const double weight =
        acoustic_shape_value_at(j, target_center[0], target_center[1],
                                target_center[2]);
    if (weight == 0.0) {
      continue;
    }

    const int32_t lnid =
        solid_nodes ? solid_nodes[8 * faceIndex + j]
                    : gpuData->lnidArrayDevice[target_eindex + j * gpuData->lenum];
    atomicAdd((solver_float*)(gpuData->porePressureOpenPoreTargetDevice +
                              lnid),
              (solver_float)(scale * weight * pressure));
  }
}

__global__ void kernelPoroelasticPressureUpdate(int32_t nharbored,
                                                gpu_data_t* gpuData,
                                                int useDrainMass,
                                                int useOpenPore)
{
  const int32_t nindex = (blockIdx.x * blockDim.x) + threadIdx.x;

  if (nindex >= nharbored) {
    return;
  }

  const solver_float mass = gpuData->porePressureMassDevice[nindex];
  solver_float pnew = 0.0;
  if (mass > 0.0) {
    const solver_float pfree =
        gpuData->porePressureTm1Device[nindex] +
        gpuData->porePressureForceDevice[nindex] / mass;
    const solver_float drainMass =
        useDrainMass ? gpuData->porePressureSurfaceDevice[nindex]
                     : (solver_float)0.0;
    const solver_float openMass =
        useOpenPore ? gpuData->porePressureOpenPoreMassDevice[nindex]
                    : (solver_float)0.0;
    const solver_float openTarget =
        useOpenPore ? gpuData->porePressureOpenPoreTargetDevice[nindex]
                    : (solver_float)0.0;
    if ((drainMass > (solver_float)0.0) ||
        (openMass > (solver_float)0.0)) {
      pnew = (mass * pfree + openTarget) /
             (mass + drainMass + openMass);
    } else {
      pnew = pfree;
    }
  }

  gpuData->porePressureTm2Device[nindex] = pnew;
  gpuData->porePressureForceDevice[nindex] = (solver_float)0.0;
}

int poroelastic_gpu_uniform_element_history(
    const double* strain_history, const double* initial_strain,
    int32_t steps, double storage_modulus, double biot_alpha,
    double initial_pressure, double* pressure_history,
    double* maximum_nodal_spread)
{
  gpu_data_t host_data = {};
  gpu_data_t* data_device = NULL;
  fvector_t* tm1_device = NULL;
  fvector_t* tm2_device = NULL;
  solver_float* pressure_tm1_device = NULL;
  solver_float* pressure_tm2_device = NULL;
  solver_float* pressure_force_device = NULL;
  solver_float* pressure_mass_device = NULL;
  int32_t* element_mapper_device = NULL;
  int32_t* lnid_device = NULL;
  double* element_h_device = NULL;
  double* element_detj_device = NULL;
  double* storage_modulus_device = NULL;
  double* shape_scale_device = NULL;
  fvector_t previous_displacement[8] = {};
  fvector_t current_displacement[8] = {};
  solver_float nodal_pressure[8] = {};
  solver_float initial_nodal_pressure[8] = {};
  solver_float nodal_mass[8] = {};
  int32_t lnid[8] = {0, 1, 2, 3, 4, 5, 6, 7};
  double shape_scale[8] = {1.0, 1.0, 1.0, 1.0,
                           1.0, 1.0, 1.0, 1.0};
  int32_t element_mapper = 0;
  double element_h = 1.0;
  double element_detj = 0.125;
  double previous_strain[6] = {};
  double maximum_spread = 0.0;
  cudaError_t error = cudaSuccess;

  if (strain_history == NULL || initial_strain == NULL ||
      pressure_history == NULL || maximum_nodal_spread == NULL ||
      steps <= 0 || !isfinite(storage_modulus) || storage_modulus <= 0.0 ||
      !isfinite(biot_alpha) || biot_alpha <= 0.0 ||
      !isfinite(initial_pressure)) {
    return -1;
  }

  for (int component = 0; component < 6; component++) {
    previous_strain[component] = initial_strain[component];
  }
  for (int node = 0; node < 8; node++) {
    initial_nodal_pressure[node] = (solver_float)initial_pressure;
    nodal_mass[node] = (solver_float)(0.125 / storage_modulus);
  }

#define PORO_TEST_ALLOC(pointer, count) do {                              \
    error = cudaMalloc((void**)&(pointer), (count) * sizeof(*(pointer))); \
    if (error != cudaSuccess) goto cleanup;                               \
  } while (0)
#define PORO_TEST_COPY_TO_DEVICE(pointer, source, count) do {              \
    error = cudaMemcpy((pointer), (source),                                \
                       (count) * sizeof(*(pointer)), cudaMemcpyHostToDevice); \
    if (error != cudaSuccess) goto cleanup;                                \
  } while (0)

  PORO_TEST_ALLOC(data_device, 1);
  PORO_TEST_ALLOC(tm1_device, 8);
  PORO_TEST_ALLOC(tm2_device, 8);
  PORO_TEST_ALLOC(pressure_tm1_device, 8);
  PORO_TEST_ALLOC(pressure_tm2_device, 8);
  PORO_TEST_ALLOC(pressure_force_device, 8);
  PORO_TEST_ALLOC(pressure_mass_device, 8);
  PORO_TEST_ALLOC(element_mapper_device, 1);
  PORO_TEST_ALLOC(lnid_device, 8);
  PORO_TEST_ALLOC(element_h_device, 1);
  PORO_TEST_ALLOC(element_detj_device, 1);
  PORO_TEST_ALLOC(storage_modulus_device, 1);
  PORO_TEST_ALLOC(shape_scale_device, 8);

  PORO_TEST_COPY_TO_DEVICE(pressure_tm1_device, initial_nodal_pressure, 8);
  PORO_TEST_COPY_TO_DEVICE(pressure_tm2_device, initial_nodal_pressure, 8);
  PORO_TEST_COPY_TO_DEVICE(pressure_mass_device, nodal_mass, 8);
  PORO_TEST_COPY_TO_DEVICE(element_mapper_device, &element_mapper, 1);
  PORO_TEST_COPY_TO_DEVICE(lnid_device, lnid, 8);
  PORO_TEST_COPY_TO_DEVICE(element_h_device, &element_h, 1);
  PORO_TEST_COPY_TO_DEVICE(element_detj_device, &element_detj, 1);
  PORO_TEST_COPY_TO_DEVICE(storage_modulus_device, &storage_modulus, 1);
  PORO_TEST_COPY_TO_DEVICE(shape_scale_device, shape_scale, 8);

  host_data.nharbored = 8;
  host_data.lenum = 1;
  host_data.tm1Device = tm1_device;
  host_data.tm2Device = tm2_device;
  host_data.porePressureTm1Device = pressure_tm1_device;
  host_data.porePressureTm2Device = pressure_tm2_device;
  host_data.porePressureForceDevice = pressure_force_device;
  host_data.porePressureMassDevice = pressure_mass_device;
  host_data.poroelasticElemMapperDevice = element_mapper_device;
  host_data.lnidArrayDevice = lnid_device;
  host_data.poroelasticElemHDevice = element_h_device;
  host_data.poroelasticElemDetJDevice = element_detj_device;
  host_data.poroelasticStorageMDevice = storage_modulus_device;
  host_data.poroelasticPressureShapeScaleDevice = shape_scale_device;
  error = cudaMemcpy(data_device, &host_data, sizeof(host_data),
                     cudaMemcpyHostToDevice);
  if (error != cudaSuccess) goto cleanup;

  for (int32_t step = 0; step < steps; step++) {
    const double* current_strain = strain_history + 6 * step;

    for (int node = 0; node < 8; node++) {
      const double x = (node & 1) ? 1.0 : 0.0;
      const double y = (node & 2) ? 1.0 : 0.0;
      const double z = (node & 4) ? 1.0 : 0.0;

      previous_displacement[node].f[0] = (solver_float)(
          previous_strain[0] * x + previous_strain[3] * y +
          previous_strain[5] * z);
      previous_displacement[node].f[1] = (solver_float)(
          previous_strain[3] * x + previous_strain[1] * y +
          previous_strain[4] * z);
      previous_displacement[node].f[2] = (solver_float)(
          previous_strain[5] * x + previous_strain[4] * y +
          previous_strain[2] * z);
      current_displacement[node].f[0] = (solver_float)(
          current_strain[0] * x + current_strain[3] * y +
          current_strain[5] * z);
      current_displacement[node].f[1] = (solver_float)(
          current_strain[3] * x + current_strain[1] * y +
          current_strain[4] * z);
      current_displacement[node].f[2] = (solver_float)(
          current_strain[5] * x + current_strain[4] * y +
          current_strain[2] * z);
    }

    PORO_TEST_COPY_TO_DEVICE(tm1_device, previous_displacement, 8);
    PORO_TEST_COPY_TO_DEVICE(tm2_device, current_displacement, 8);
    error = cudaMemset(pressure_force_device, 0,
                       8 * sizeof(*pressure_force_device));
    if (error != cudaSuccess) goto cleanup;

    kernelPoroelasticPressureForce<<<1, 1>>>(
        1, data_device, 1.0, 0.0, 0.0, 0.0, 0, biot_alpha);
    error = cudaGetLastError();
    if (error != cudaSuccess) goto cleanup;
    kernelPoroelasticPressureUpdate<<<1, 8>>>(8, data_device, 0, 0);
    error = cudaGetLastError();
    if (error != cudaSuccess) goto cleanup;
    error = cudaDeviceSynchronize();
    if (error != cudaSuccess) goto cleanup;
    error = cudaMemcpy(nodal_pressure, pressure_tm2_device,
                       sizeof(nodal_pressure), cudaMemcpyDeviceToHost);
    if (error != cudaSuccess) goto cleanup;

    double pressure_mean = 0.0;
    for (int node = 0; node < 8; node++) {
      pressure_mean += (double)nodal_pressure[node] / 8.0;
    }
    for (int node = 0; node < 8; node++) {
      const double spread = fabs((double)nodal_pressure[node] - pressure_mean);
      if (spread > maximum_spread) maximum_spread = spread;
    }
    pressure_history[step] = pressure_mean;

    error = cudaMemcpy(pressure_tm1_device, pressure_tm2_device,
                       8 * sizeof(*pressure_tm1_device),
                       cudaMemcpyDeviceToDevice);
    if (error != cudaSuccess) goto cleanup;
    for (int component = 0; component < 6; component++) {
      previous_strain[component] = current_strain[component];
    }
  }

  *maximum_nodal_spread = maximum_spread;

cleanup:
  if (shape_scale_device != NULL) cudaFree(shape_scale_device);
  if (storage_modulus_device != NULL) cudaFree(storage_modulus_device);
  if (element_detj_device != NULL) cudaFree(element_detj_device);
  if (element_h_device != NULL) cudaFree(element_h_device);
  if (lnid_device != NULL) cudaFree(lnid_device);
  if (element_mapper_device != NULL) cudaFree(element_mapper_device);
  if (pressure_mass_device != NULL) cudaFree(pressure_mass_device);
  if (pressure_force_device != NULL) cudaFree(pressure_force_device);
  if (pressure_tm2_device != NULL) cudaFree(pressure_tm2_device);
  if (pressure_tm1_device != NULL) cudaFree(pressure_tm1_device);
  if (tm2_device != NULL) cudaFree(tm2_device);
  if (tm1_device != NULL) cudaFree(tm1_device);
  if (data_device != NULL) cudaFree(data_device);
#undef PORO_TEST_COPY_TO_DEVICE
#undef PORO_TEST_ALLOC

  if (error != cudaSuccess) {
    fprintf(stderr, "BrickUP uniform-element GPU regression failed: %s\n",
            cudaGetErrorString(error));
    return -1;
  }
  return 0;
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

  /* CPML: keep BKT material damping out of PML elements. Inside the PML
   * zone, BKT-family runs use native elastic stiffness plus the CPML
   * correction; stretching the BKT memory variables over-damped large
   * heterogeneous runs to numerical zero. */
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

  /* CPML: keep BKT material damping out of PML elements; see
   * kernelDampingCalcConvBkt2 for the policy. */
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


/* Damping Calc-Conv Kernel -- BKT3/BKT3F fast path (3 Maxwell
 * mechanisms). Numerically a faithful, bit-identical specialization of
 * the generic kernelDampingCalcConv for the spatially-uniform-g case:
 * the 18 spectral coefficients are hoisted to scalar kernel args
 * (computed once on the host -> no per-element/-node exp() and no
 * g0/g1/g2 array traffic), while the per-element weights a*,b* are
 * still read per element. Gates and arithmetic order match the generic
 * kernel exactly. Dispatched only when g_dampingGUniform==1; the
 * generic kernel is the heterogeneous-g fallback. */
__global__ void kernelDampingCalcConvBkt3(int32_t lenum,
                                          int32_t numcomp,
                                          gpu_data_t *gpuData,
                                          double rmax,
                                          double shear_coef_1,
                                          double shear_coef_2,
                                          double shear_exp_coef_0,
                                          double shear_coef_3,
                                          double shear_coef_4,
                                          double shear_exp_coef_1,
                                          double shear_coef_5,
                                          double shear_coef_6,
                                          double shear_exp_coef_2,
                                          double kappa_coef_1,
                                          double kappa_coef_2,
                                          double kappa_exp_coef_0,
                                          double kappa_coef_3,
                                          double kappa_coef_4,
                                          double kappa_exp_coef_1,
                                          double kappa_coef_5,
                                          double kappa_coef_6,
                                          double kappa_exp_coef_2,
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

  /* CPML: keep BKT material damping out of PML elements; see
   * kernelDampingCalcConvBkt2 for the policy. */
  if (gpuData->pmlFlagArrayDevice != NULL &&
      gpuData->pmlFlagArrayDevice[eindex] != 0) {
    return;
  }

  const double a0_shear = gpuData->a0_shearArrayDevice[eindex];
  const double a1_shear = gpuData->a1_shearArrayDevice[eindex];
  const double a2_shear = gpuData->a2_shearArrayDevice[eindex];
  const double b_shear = gpuData->b_shearArrayDevice[eindex];
  const double coef_shear = (rmax != 0.0) ? (b_shear / rmax) : 0.0;
  const int do_conv_update_shear =
      ((shear_coef_1 != 0.0) && (shear_coef_3 != 0.0));
  const int do_damping_calc_shear =
      ((a0_shear + a1_shear + a2_shear + b_shear) != 0.0);

  const double a0_kappa = gpuData->a0_kappaArrayDevice[eindex];
  const double a1_kappa = gpuData->a1_kappaArrayDevice[eindex];
  const double a2_kappa = gpuData->a2_kappaArrayDevice[eindex];
  const double b_kappa = gpuData->b_kappaArrayDevice[eindex];
  const double coef_kappa = (rmax != 0.0) ? (b_kappa / rmax) : 0.0;
  const int do_conv_update_kappa =
      ((kappa_coef_1 != 0.0) && (kappa_coef_3 != 0.0));
  const int do_damping_calc_kappa =
      ((a0_kappa + a1_kappa + a2_kappa + b_kappa) != 0.0);

  fvector_t tm1Disp, tm2Disp;
  solver_float *f0_tm1[3], *f1_tm1[3], *f2_tm1[3], *convVec[3];
  int32_t offset0;

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

    /* SHEAR convolution (3 mechanisms) */
    f0_tm1[0] = (solver_float *)gpuData->conv_shear_1Device + offset0;
    f0_tm1[1] = (solver_float *)gpuData->conv_shear_1Device + offset1;
    f0_tm1[2] = (solver_float *)gpuData->conv_shear_1Device + offset2;
    f1_tm1[0] = (solver_float *)gpuData->conv_shear_2Device + offset0;
    f1_tm1[1] = (solver_float *)gpuData->conv_shear_2Device + offset1;
    f1_tm1[2] = (solver_float *)gpuData->conv_shear_2Device + offset2;
    f2_tm1[0] = (solver_float *)gpuData->conv_shear_3Device + offset0;
    f2_tm1[1] = (solver_float *)gpuData->conv_shear_3Device + offset1;
    f2_tm1[2] = (solver_float *)gpuData->conv_shear_3Device + offset2;

    if (do_conv_update_shear) {
      *f0_tm1[0] = shear_coef_2 * tm1Disp.f[0] + shear_coef_1 * tm2Disp.f[0] + shear_exp_coef_0 * (*f0_tm1[0]);
      *f1_tm1[0] = shear_coef_4 * tm1Disp.f[0] + shear_coef_3 * tm2Disp.f[0] + shear_exp_coef_1 * (*f1_tm1[0]);

      *f0_tm1[1] = shear_coef_2 * tm1Disp.f[1] + shear_coef_1 * tm2Disp.f[1] + shear_exp_coef_0 * (*f0_tm1[1]);
      *f1_tm1[1] = shear_coef_4 * tm1Disp.f[1] + shear_coef_3 * tm2Disp.f[1] + shear_exp_coef_1 * (*f1_tm1[1]);

      *f0_tm1[2] = shear_coef_2 * tm1Disp.f[2] + shear_coef_1 * tm2Disp.f[2] + shear_exp_coef_0 * (*f0_tm1[2]);
      *f1_tm1[2] = shear_coef_4 * tm1Disp.f[2] + shear_coef_3 * tm2Disp.f[2] + shear_exp_coef_1 * (*f1_tm1[2]);

      *f2_tm1[0] = shear_coef_6 * tm1Disp.f[0] + shear_coef_5 * tm2Disp.f[0] + shear_exp_coef_2 * (*f2_tm1[0]);
      *f2_tm1[1] = shear_coef_6 * tm1Disp.f[1] + shear_coef_5 * tm2Disp.f[1] + shear_exp_coef_2 * (*f2_tm1[1]);
      *f2_tm1[2] = shear_coef_6 * tm1Disp.f[2] + shear_coef_5 * tm2Disp.f[2] + shear_exp_coef_2 * (*f2_tm1[2]);
    }

    convVec[0] = (solver_float *)gpuData->shearVectorDevice + offset0;
    convVec[1] = (solver_float *)gpuData->shearVectorDevice + offset1;
    convVec[2] = (solver_float *)gpuData->shearVectorDevice + offset2;

    if (do_damping_calc_shear) {
      *convVec[0] = coef_shear * (tm1Disp.f[0] - tm2Disp.f[0])
          - (a0_shear * (*f0_tm1[0]) + a1_shear * (*f1_tm1[0])) + tm1Disp.f[0];
      *convVec[1] = coef_shear * (tm1Disp.f[1] - tm2Disp.f[1])
          - (a0_shear * (*f0_tm1[1]) + a1_shear * (*f1_tm1[1])) + tm1Disp.f[1];
      *convVec[2] = coef_shear * (tm1Disp.f[2] - tm2Disp.f[2])
          - (a0_shear * (*f0_tm1[2]) + a1_shear * (*f1_tm1[2])) + tm1Disp.f[2];

      *convVec[0] -= a2_shear * (*f2_tm1[0]);
      *convVec[1] -= a2_shear * (*f2_tm1[1]);
      *convVec[2] -= a2_shear * (*f2_tm1[2]);
    } else {
      *convVec[0] = tm1Disp.f[0];
      *convVec[1] = tm1Disp.f[1];
      *convVec[2] = tm1Disp.f[2];
    }

    /* DILATATION (kappa) convolution (3 mechanisms) */
    f0_tm1[0] = (solver_float *)gpuData->conv_kappa_1Device + offset0;
    f0_tm1[1] = (solver_float *)gpuData->conv_kappa_1Device + offset1;
    f0_tm1[2] = (solver_float *)gpuData->conv_kappa_1Device + offset2;
    f1_tm1[0] = (solver_float *)gpuData->conv_kappa_2Device + offset0;
    f1_tm1[1] = (solver_float *)gpuData->conv_kappa_2Device + offset1;
    f1_tm1[2] = (solver_float *)gpuData->conv_kappa_2Device + offset2;
    f2_tm1[0] = (solver_float *)gpuData->conv_kappa_3Device + offset0;
    f2_tm1[1] = (solver_float *)gpuData->conv_kappa_3Device + offset1;
    f2_tm1[2] = (solver_float *)gpuData->conv_kappa_3Device + offset2;

    if (do_conv_update_kappa) {
      *f0_tm1[0] = kappa_coef_2 * tm1Disp.f[0] + kappa_coef_1 * tm2Disp.f[0] + kappa_exp_coef_0 * (*f0_tm1[0]);
      *f1_tm1[0] = kappa_coef_4 * tm1Disp.f[0] + kappa_coef_3 * tm2Disp.f[0] + kappa_exp_coef_1 * (*f1_tm1[0]);

      *f0_tm1[1] = kappa_coef_2 * tm1Disp.f[1] + kappa_coef_1 * tm2Disp.f[1] + kappa_exp_coef_0 * (*f0_tm1[1]);
      *f1_tm1[1] = kappa_coef_4 * tm1Disp.f[1] + kappa_coef_3 * tm2Disp.f[1] + kappa_exp_coef_1 * (*f1_tm1[1]);

      *f0_tm1[2] = kappa_coef_2 * tm1Disp.f[2] + kappa_coef_1 * tm2Disp.f[2] + kappa_exp_coef_0 * (*f0_tm1[2]);
      *f1_tm1[2] = kappa_coef_4 * tm1Disp.f[2] + kappa_coef_3 * tm2Disp.f[2] + kappa_exp_coef_1 * (*f1_tm1[2]);

      *f2_tm1[0] = kappa_coef_6 * tm1Disp.f[0] + kappa_coef_5 * tm2Disp.f[0] + kappa_exp_coef_2 * (*f2_tm1[0]);
      *f2_tm1[1] = kappa_coef_6 * tm1Disp.f[1] + kappa_coef_5 * tm2Disp.f[1] + kappa_exp_coef_2 * (*f2_tm1[1]);
      *f2_tm1[2] = kappa_coef_6 * tm1Disp.f[2] + kappa_coef_5 * tm2Disp.f[2] + kappa_exp_coef_2 * (*f2_tm1[2]);
    }

    convVec[0] = (solver_float *)gpuData->kappaVectorDevice + offset0;
    convVec[1] = (solver_float *)gpuData->kappaVectorDevice + offset1;
    convVec[2] = (solver_float *)gpuData->kappaVectorDevice + offset2;

    if (do_damping_calc_kappa) {
      *convVec[0] = coef_kappa * (tm1Disp.f[0] - tm2Disp.f[0])
          - (a0_kappa * (*f0_tm1[0]) + a1_kappa * (*f1_tm1[0])) + tm1Disp.f[0];
      *convVec[1] = coef_kappa * (tm1Disp.f[1] - tm2Disp.f[1])
          - (a0_kappa * (*f0_tm1[1]) + a1_kappa * (*f1_tm1[1])) + tm1Disp.f[1];
      *convVec[2] = coef_kappa * (tm1Disp.f[2] - tm2Disp.f[2])
          - (a0_kappa * (*f0_tm1[2]) + a1_kappa * (*f1_tm1[2])) + tm1Disp.f[2];

      *convVec[0] -= a2_kappa * (*f2_tm1[0]);
      *convVec[1] -= a2_kappa * (*f2_tm1[1]);
      *convVec[2] -= a2_kappa * (*f2_tm1[2]);
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

  /* CPML: keep BKT material damping force out of PML elements. The solver
   * adds native elastic stiffness for PML elements in BKT-family runs before
   * applying the CPML correction. */
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


/* The common BKT2 spectrum has identical displacement-history recurrences
 * for shear and bulk, independent of their Q-dependent weights. With zero
 * initial histories, one pair in conv_shear_{1,2} therefore suffices, including
 * bulk-only elements. Do not alias the legacy arrays: VT and fallback kernels
 * still use their independent histories. No precision or constitutive change.
 *
 * Keep solver_float assignment boundaries and the force operator's evaluation
 * order identical to kernelDampingCalcConvBkt2 + kernelDampingCalcLocal.
 * Nodal atomic summation order remains nondeterministic, as on the old path.
 */
/* The unconstrained fused kernel's register demand limits occupancy. This
 * bound trades some compiler spills for more resident warps; measured against
 * the separate kernels on GB200. Keep the launcher at 128 threads or fewer. */
__global__ __launch_bounds__(128, 6) void kernelDampingBkt2Fused(
    int32_t lenum, gpu_data_t* gpuData, double rmax,
    double coef1, double coef2, double expCoef0,
    double coef3, double coef4, double expCoef1,
    int32_t elementCount, const int32_t* elementMapper, int saveVectors) {
  const int32_t idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= elementCount) return;
  const int32_t eindex = elementMapper[idx];
  if (gpuData->pmlFlagArrayDevice != NULL &&
      gpuData->pmlFlagArrayDevice[eindex] != 0) return;

  const int32_t numcomp = lenum * 8;
  const double a0s = gpuData->a0_shearArrayDevice[eindex];
  const double a1s = gpuData->a1_shearArrayDevice[eindex];
  const double bs = gpuData->b_shearArrayDevice[eindex];
  const double a0k = gpuData->a0_kappaArrayDevice[eindex];
  const double a1k = gpuData->a1_kappaArrayDevice[eindex];
  const double bk = gpuData->b_kappaArrayDevice[eindex];
  const double cs = (rmax != 0.0) ? bs / rmax : 0.0;
  const double ck = (rmax != 0.0) ? bk / rmax : 0.0;
  const int historyS = (a0s != 0.0 || a1s != 0.0);
  const int historyK = (a0k != 0.0 || a1k != 0.0);
  const int update = (historyS || historyK) && coef1 != 0.0 && coef3 != 0.0;
  fvector_t shear[8], bulk[8];
  int contributeS = 0, contributeK = 0;

  for (int i = 0; i < 8; ++i) {
    const int32_t offset = eindex + i * lenum;
    const int32_t node = gpuData->lnidArrayDevice[offset];
    const fvector_t u = gpuData->tm1Device[node];
    const fvector_t prev = gpuData->tm2Device[node];
    for (int d = 0; d < 3; ++d) {
      const int32_t pos = offset + d * numcomp;
      solver_float h0 = 0.0, h1 = 0.0;
      if (historyS || historyK) {
        h0 = ((solver_float*)gpuData->conv_shear_1Device)[pos];
        h1 = ((solver_float*)gpuData->conv_shear_2Device)[pos];
        if (update) {
          h0 = coef2 * u.f[d] + coef1 * prev.f[d] + expCoef0 * h0;
          h1 = coef4 * u.f[d] + coef3 * prev.f[d] + expCoef1 * h1;
          ((solver_float*)gpuData->conv_shear_1Device)[pos] = h0;
          ((solver_float*)gpuData->conv_shear_2Device)[pos] = h1;
        }
      }
      if (historyS) {
        shear[i].f[d] = cs * (u.f[d] - prev.f[d]) -
                        (a0s * h0 + a1s * h1) + u.f[d];
      } else if (bs != 0.0) {
        shear[i].f[d] = cs * (u.f[d] - prev.f[d]) + u.f[d];
      } else {
        shear[i].f[d] = u.f[d];
      }
      if (historyK) {
        bulk[i].f[d] = ck * (u.f[d] - prev.f[d]) -
                       (a0k * h0 + a1k * h1) + u.f[d];
      } else if (bk != 0.0) {
        bulk[i].f[d] = ck * (u.f[d] - prev.f[d]) + u.f[d];
      } else {
        bulk[i].f[d] = u.f[d];
      }
      contributeS |= fabs(shear[i].f[d]) > UNDERFLOW_CAP_STIFFNESS;
      contributeK |= fabs(bulk[i].f[d]) > UNDERFLOW_CAP_STIFFNESS;
      /* Station stress gathers consume the previous step's vectors. Save
       * every step when requested, even if station output is decimated. */
      if (saveVectors) {
        ((solver_float*)gpuData->shearVectorDevice)[pos] = shear[i].f[d];
        ((solver_float*)gpuData->kappaVectorDevice)[pos] = bulk[i].f[d];
      }
    }
  }
  if (!contributeS && !contributeK) return;

  double atu[24], firstVec[24];
  for (int i = 0; i < 24; ++i) firstVec[i] = 0.0;
  const uint8_t flags = gpuData->acousticFaceFlagsArrayDevice[eindex];
  if (contributeS) {
    const double mu = -0.5625 * gpuData->c1ArrayDevice[eindex];
    aTransposeU(shear, atu);
    firstVector_mu(atu, firstVec, mu);
  }
  if (contributeK) {
    const double kappa = -0.5625 * (gpuData->c2ArrayDevice[eindex] +
                                  2. / 3. * gpuData->c1ArrayDevice[eindex]);
    if (flags != 0) acoustic_project_interface_nodes(bulk, flags);
    aTransposeU(bulk, atu);
    firstVector_kappa(atu, firstVec, kappa);
  }
  /* Match the legacy output buffer too: au adds to its input buffer in a
   * SINGLE_PRECISION_SOLVER build, but overwrites it in double precision. */
  au(bulk, firstVec);
  if (flags != 0) acoustic_project_interface_nodes(bulk, flags);
  for (int i = 0; i < 8; ++i) {
    const int32_t node = gpuData->lnidArrayDevice[eindex + i * lenum];
    for (int d = 0; d < 3; ++d)
      atomicAdd(&gpuData->forceDevice[node].f[d], bulk[i].f[d]);
  }
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

    /* Optional CPML outer-face Dirichlet clamp. pmlClampNodeMaskDevice is
     * non-NULL only when cpml_clamp_outer=yes; it marks nodes lying on
     * x_min/x_max/y_min/y_max/z_max domain-extent faces touched by a PML
     * element. Force u=0 at those nodes; the CPML layer attenuates waves
     * before and after any reflection at this outer termination. The
     * top face is not selected by Z, but its lateral rim can be flagged. */
    if (gpuData->pmlClampNodeMaskDevice != NULL &&
        gpuData->pmlClampNodeMaskDevice[nindex] != 0) {
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


__global__ void kernelPMLElasticStep(int32_t lenum, gpu_data_t* gpuData,
                                    int32_t count) {
    const int32_t index=blockIdx.x*blockDim.x+threadIdx.x;
    if (index>=count) return;
    const int32_t e=gpuData->pmlElasticElemMapperDevice[index];
    const cpml_element_t& c=gpuData->pmlElasticCoeffDevice[index];
    const int native=gpuData->pmlNativeStiffness;
    solver_float* state=gpuData->pmlElasticStateDevice+index;
    double u[8][3],force[8][3];
    fvector_t effective[8],stiff_force[8];
    for (int n=0;n<8;++n) {
        const int32_t ln=gpuData->lnidArrayDevice[e+n*lenum];
        for (int i=0;i<3;++i) {
            u[n][i]=(double)gpuData->tm1Device[ln].f[i];
            effective[n].f[i]=(solver_float)(u[n][i]+c.b_over_dt*
                (u[n][i]-(double)gpuData->tm2Device[ln].f[i]));
        }
    }
    const uint8_t faces=gpuData->acousticFaceFlagsArrayDevice[e];
    if (faces) acoustic_project_interface_nodes(effective,faces);
    double ue[8][3];
    for (int n=0;n<8;++n)
        for (int i=0;i<3;++i) ue[n][i]=(double)effective[n].f[i];
    cpml_stiffness(c,3,native,ue,force,state,(size_t)count);
    for (int n=0;n<8;++n)
        for (int i=0;i<3;++i) stiff_force[n].f[i]=(solver_float)force[n][i];
    if (faces) acoustic_project_interface_nodes(stiff_force,faces);
    const int offset=cpml_stiffness_states(3,native);
    for (int n=0;n<8;++n) {
        const int32_t ln=gpuData->lnidArrayDevice[e+n*lenum];
        for (int i=0;i<3;++i) {
            const size_t channel=(size_t)offset+3*(3*n+i);
            const double inertia=c.nodal_mass_dt2*cpml_inertia_force(
                c.inertia,u[n][i],state+channel*count,(size_t)count);
            atomicAdd(&gpuData->forceDevice[ln].f[i],
                (solver_float)((double)stiff_force[n].f[i]+inertia));
        }
    }
}

__global__ void kernelPMLAcousticStep(int32_t lenum, gpu_data_t* gpuData,
                                     int32_t count) {
    const int32_t index=blockIdx.x*blockDim.x+threadIdx.x;
    if (index>=count) return;
    const int32_t e=gpuData->pmlAcousticElemMapperDevice[index];
    const cpml_element_t& c=gpuData->pmlAcousticCoeffDevice[index];
    const int native=gpuData->pmlNativeStiffnessAcoustic;
    solver_float* state=gpuData->pmlAcousticStateDevice+index;
    double u[8][3]={},force[8][3];
    for (int n=0;n<8;++n) {
        const int32_t ln=gpuData->lnidArrayDevice[e+n*lenum];
        u[n][0]=(double)gpuData->pressureTm1Device[ln];
    }
    cpml_stiffness(c,1,native,u,force,state,(size_t)count);
    const int offset=cpml_stiffness_states(1,native);
    for (int n=0;n<8;++n) {
        const int32_t ln=gpuData->lnidArrayDevice[e+n*lenum];
        const size_t channel=(size_t)offset+3*n;
        const double inertia=c.nodal_mass_dt2*cpml_inertia_force(
            c.inertia,u[n][0],state+channel*count,(size_t)count);
        atomicAdd(gpuData->pressureForceDevice+ln,
                  (solver_float)(force[n][0]+inertia));
    }
}

__global__ void kernelPMLAcousticGravitySurfaceStep(
    int32_t lenum,gpu_data_t* gpuData,int32_t count) {
    const int32_t index=blockIdx.x*blockDim.x+threadIdx.x;
    if (index>=count) return;
    const int32_t pml=gpuData->pml_acoustic_gravitySurfacePmlIndexArrayDevice[index];
    const int32_t e=gpuData->pmlAcousticElemMapperDevice[pml];
    const cpml_surface_t& c=gpuData->pmlGravityCoeffDevice[index];
    /* Surface Jacobian S_x*S_y, with one history per lumped top node.
     * Cross terms and coincident poles use the same inertia builder as bulk. */
    for (int n=0;n<4;++n) {
        const int32_t ln=gpuData->lnidArrayDevice[e+n*lenum];
        const double force=c.nodal_mass_dt2*cpml_inertia_force(
            c.inertia,(double)gpuData->pressureTm1Device[ln],
            gpuData->pmlGravityStateDevice+index+(size_t)(2*n)*count,
            (size_t)count);
        atomicAdd(gpuData->pressureForceDevice+ln,(solver_float)force);
    }
}


/* ---------------------------------------------------------------------
 * GPU station strain / viscoelastic-stress gather.
 *
 * One thread per station. For kinematic strain it reads the current
 * displacement at the owning element's eight nodes. For viscoelastic
 * stress it reads the matching slots of shearVectorDevice (deviatoric
 * channel) and kappaVectorDevice (volumetric channel), written by the
 * damping conv kernel each timestep for BKT-family runs, and uses the
 * same dev/vol split as the damping force kernel:
 *
 *   sigma_visco_ij = 2*mu * (eps^s_ij - 1/3 tr(eps^s) delta_ij)
 *                  + K   * tr(eps^k) delta_ij     (K = lambda + 2*mu/3)
 *
 * where eps^s, eps^k are the strains derived from shearVector and
 * kappaVector respectively. This reduces to the elastic constitutive
 * when shearVector = kappaVector = u (no SLS memory), so the value
 * the host writes for the Sig_* columns is the same as v2's elastic
 * point_stress in the no-damping limit; under BKT-family it carries
 * the full SLS-memory contribution.
 *
 * For regular (non-cut) cube stations the strain at the station's
 * local coords is the trilinear-hex point_dxi + point_strain (same as
 * nonlinear.cc's point_strain). For VT cut-cell stations (teth_i>=0)
 * the constant-per-tet strain is computed via the same algebra as
 * point_strain_tetrah (nonlinear.cc), ported here inline so the kernel
 * needs no host call.
 *
 * Requested tensors are packed in gpuData->stationTensorArrayDevice.
 * outputStride is 6 for a single tensor or 12 for both; a negative
 * strainOffset/stressOffset disables that path. Every tensor uses the
 * XX,YY,ZZ,XY,YZ,XZ component order.
 *
 * Acoustic stations may be processed too, but the host row writer ignores
 * gathered tensors for them (the acoustic branch returns before the tensor
 * columns). Their shearVector/kappaVector slots remain zero because acoustic
 * elements are excluded from the BKT damping mapper.
 * --------------------------------------------------------------------- */

__device__ static inline void v3_cube_strain(
    const double *ux, const double *uy, const double *uz,
    double csi, double eta, double dze, double h,
    double *exx, double *eyy, double *ezz,
    double *exy, double *eyz, double *exz)
{
    /* mirror nonlinear.cc point_dxi + point_strain (trilinear hex8). */
    const int xi0[8] = {-1, 1,-1, 1,-1, 1,-1, 1};
    const int xi1[8] = {-1,-1, 1, 1,-1,-1, 1, 1};
    const int xi2[8] = {-1,-1,-1,-1, 1, 1, 1, 1};
    const double Jij = 0.25 / h;
    double Exx=0, Eyy=0, Ezz=0, Exy=0, Eyz=0, Exz=0;
    #pragma unroll
    for (int i = 0; i < 8; i++) {
        double dx = Jij * (double)xi0[i] * (1.0 + xi1[i]*eta) * (1.0 + xi2[i]*dze);
        double dy = Jij * (1.0 + xi0[i]*csi) * (double)xi1[i] * (1.0 + xi2[i]*dze);
        double dz = Jij * (1.0 + xi0[i]*csi) * (1.0 + xi1[i]*eta) * (double)xi2[i];
        Exx += dx * ux[i];
        Eyy += dy * uy[i];
        Ezz += dz * uz[i];
        Exy += 0.5 * (dy * ux[i] + dx * uy[i]);
        Eyz += 0.5 * (dz * uy[i] + dy * uz[i]);
        Exz += 0.5 * (dz * ux[i] + dx * uz[i]);
    }
    *exx = Exx; *eyy = Eyy; *ezz = Ezz;
    *exy = Exy; *eyz = Eyz; *exz = Exz;
}

__device__ static inline void v3_tet_strain(
    const double *ux, const double *uy, const double *uz,
    double h, int teth_i, int cube_part,
    double *exx, double *eyy, double *ezz,
    double *exy, double *eyz, double *exz)
{
    /* Faithful port of nonlinear.cc point_strain_tetrah. 5 tets x 2
     * cube_parts = 10 cases. The N0..N3 picks specific cube nodes for
     * each tet (the same indexing the CPU function uses). All cases
     * yield the constant-per-tet strain (linear shape functions over
     * the tet). Cube node ordering is the standard Hercules hex8
     * (xi[3][8] sign-table). */
    double a_xx=0, a_yy=0, a_zz=0, a_xy=0, a_yz=0, a_xz=0;
    int N0=0, N1=0, N2=0, N3=0;
    if (cube_part == 1) {
        switch (teth_i) {
            case 0:
                N0=0; N1=2; N2=1; N3=4;
                a_xx = -(ux[N0] - ux[N2]) / h;
                a_yy = -(uy[N0] - uy[N1]) / h;
                a_zz = -(uz[N0] - uz[N3]) / h;
                a_xy = -0.5*(ux[N0]-ux[N1] + uy[N0]-uy[N2]) / h;
                a_yz = -0.5*(uy[N0]-uy[N3] + uz[N0]-uz[N1]) / h;
                a_xz = -0.5*(ux[N0]-ux[N3] + uz[N0]-uz[N2]) / h;
                break;
            case 1:
                N0=3; N1=1; N2=2; N3=7;
                a_xx =  (ux[N0] - ux[N2]) / h;
                a_yy =  (uy[N0] - uy[N1]) / h;
                a_zz = -(uz[N0] - uz[N3]) / h;
                a_xy =  0.5*(ux[N0]-ux[N1] + uy[N0]-uy[N2]) / h;
                a_yz = -0.5*(uy[N0]-uy[N3] - uz[N0]+uz[N1]) / h;
                a_xz = -0.5*(ux[N0]-ux[N3] - uz[N0]+uz[N2]) / h;
                break;
            case 2:
                N0=6; N1=4; N2=7; N3=2;
                a_xx = -(ux[N0] - ux[N2]) / h;
                a_yy =  (uy[N0] - uy[N1]) / h;
                a_zz =  (uz[N0] - uz[N3]) / h;
                a_xy =  0.5*(ux[N0]-ux[N1] - uy[N0]+uy[N2]) / h;
                a_yz =  0.5*(uy[N0]-uy[N3] + uz[N0]-uz[N1]) / h;
                a_xz =  0.5*(ux[N0]-ux[N3] - uz[N0]+uz[N2]) / h;
                break;
            case 3:
                N0=5; N1=7; N2=4; N3=1;
                a_xx =  (ux[N0] - ux[N2]) / h;
                a_yy = -(uy[N0] - uy[N1]) / h;
                a_zz =  (uz[N0] - uz[N3]) / h;
                a_xy = -0.5*(ux[N0]-ux[N1] - uy[N0]+uy[N2]) / h;
                a_yz =  0.5*(uy[N0]-uy[N3] - uz[N0]+uz[N1]) / h;
                a_xz =  0.5*(ux[N0]-ux[N3] + uz[N0]-uz[N2]) / h;
                break;
            case 4:
                N0=2; N1=4; N2=7; N3=1;
                a_xx = -0.5*(ux[N0]+ux[N1]-ux[N2]-ux[N3]) / h;
                a_yy =  0.5*(uy[N0]-uy[N1]+uy[N2]-uy[N3]) / h;
                a_zz = -0.5*(uz[N0]-uz[N1]-uz[N2]+uz[N3]) / h;
                a_xy =  0.25*(ux[N0]-ux[N1]+ux[N2]-ux[N3] - uy[N0]-uy[N1]+uy[N2]+uy[N3]) / h;
                a_yz = -0.25*(uy[N0]-uy[N1]-uy[N2]+uy[N3] - uz[N0]+uz[N1]-uz[N2]+uz[N3]) / h;
                a_xz = -0.25*(ux[N0]-ux[N1]-ux[N2]+ux[N3] + uz[N0]+uz[N1]-uz[N2]-uz[N3]) / h;
                break;
            default: break;
        }
    } else {
        switch (teth_i) {
            case 0:
                N0=0; N1=3; N2=1; N3=5;
                a_xx = -(ux[N0] - ux[N2]) / h;
                a_yy =  (uy[N1] - uy[N2]) / h;
                a_zz = -(uz[N2] - uz[N3]) / h;
                a_xy =  0.5*(ux[N1]-ux[N2] - uy[N0]+uy[N2]) / h;
                a_yz = -0.5*(uy[N2]-uy[N3] - uz[N1]+uz[N2]) / h;
                a_xz = -0.5*(ux[N2]-ux[N3] + uz[N0]-uz[N2]) / h;
                break;
            case 1:
                N0=0; N1=2; N2=3; N3=6;
                a_xx = -(ux[N1] - ux[N2]) / h;
                a_yy = -(uy[N0] - uy[N1]) / h;
                a_zz = -(uz[N1] - uz[N3]) / h;
                a_xy = -0.5*(ux[N0]-ux[N1] + uy[N1]-uy[N2]) / h;
                a_yz = -0.5*(uy[N1]-uy[N3] + uz[N0]-uz[N1]) / h;
                a_xz = -0.5*(ux[N1]-ux[N3] + uz[N1]-uz[N2]) / h;
                break;
            case 2:
                N0=4; N1=5; N2=6; N3=0;
                a_xx = -(ux[N0] - ux[N1]) / h;
                a_yy = -(uy[N0] - uy[N2]) / h;
                a_zz =  (uz[N0] - uz[N3]) / h;
                a_xy = -0.5*(ux[N0]-ux[N2] + uy[N0]-uy[N1]) / h;
                a_yz =  0.5*(uy[N0]-uy[N3] - uz[N0]+uz[N2]) / h;
                a_xz =  0.5*(ux[N0]-ux[N3] - uz[N0]+uz[N1]) / h;
                break;
            case 3:
                N0=6; N1=5; N2=7; N3=3;
                a_xx = -(ux[N0] - ux[N2]) / h;
                a_yy = -(uy[N1] - uy[N2]) / h;
                a_zz =  (uz[N2] - uz[N3]) / h;
                a_xy = -0.5*(ux[N1]-ux[N2] + uy[N0]-uy[N2]) / h;
                a_yz =  0.5*(uy[N2]-uy[N3] - uz[N1]+uz[N2]) / h;
                a_xz =  0.5*(ux[N2]-ux[N3] - uz[N0]+uz[N2]) / h;
                break;
            case 4:
                N0=0; N1=6; N2=3; N3=5;
                a_xx = -0.5*(ux[N0]+ux[N1]-ux[N2]-ux[N3]) / h;
                a_yy = -0.5*(uy[N0]-uy[N1]-uy[N2]+uy[N3]) / h;
                a_zz = -0.5*(uz[N0]-uz[N1]+uz[N2]-uz[N3]) / h;
                a_xy = -0.25*(ux[N0]-ux[N1]-ux[N2]+ux[N3] + uy[N0]+uy[N1]-uy[N2]-uy[N3]) / h;
                a_yz = -0.25*(uy[N0]-uy[N1]+uy[N2]-uy[N3] + uz[N0]-uz[N1]-uz[N2]+uz[N3]) / h;
                a_xz = -0.25*(ux[N0]-ux[N1]+ux[N2]-ux[N3] + uz[N0]+uz[N1]-uz[N2]-uz[N3]) / h;
                break;
            default: break;
        }
    }
    *exx=a_xx; *eyy=a_yy; *ezz=a_zz; *exy=a_xy; *eyz=a_yz; *exz=a_xz;
}

__global__ void kernelStationVEStressGather(int32_t stationCount,
                                            int32_t outputStride,
                                            int32_t strainOffset,
                                            int32_t stressOffset,
                                            gpu_data_t* gpuData)
{
    int32_t lin = blockIdx.x * blockDim.x + threadIdx.x;
    if (lin >= stationCount) return;

    const int32_t eindex = gpuData->stationGatherEindexDevice[lin];
    const double  csi    = gpuData->stationGatherCsiDevice[lin];
    const double  eta    = gpuData->stationGatherEtaDevice[lin];
    const double  dze    = gpuData->stationGatherDzetaDevice[lin];
    const double  h      = gpuData->stationGatherHDevice[lin];
    const int     teth_i = gpuData->stationGatherTethIDevice[lin];
    const int     cp     = gpuData->stationGatherCubePartDevice[lin];

    const int32_t lenum   = gpuData->lenum;
    double* const stationOut =
        gpuData->stationTensorArrayDevice + (size_t)lin * outputStride;

    if (strainOffset >= 0) {
        double ux[8], uy[8], uz[8];
        #pragma unroll
        for (int i = 0; i < 8; i++) {
            const int32_t lnid =
                gpuData->lnidArrayDevice[eindex + i * lenum];
            const fvector_t u = gpuData->tm1Device[lnid];
            ux[i] = (double)u.f[0];
            uy[i] = (double)u.f[1];
            uz[i] = (double)u.f[2];
        }

        double exx, eyy, ezz, exy, eyz, exz;
        if (teth_i < 0) {
            v3_cube_strain(ux, uy, uz, csi, eta, dze, h,
                           &exx, &eyy, &ezz, &exy, &eyz, &exz);
        } else {
            v3_tet_strain(ux, uy, uz, h, teth_i, cp,
                          &exx, &eyy, &ezz, &exy, &eyz, &exz);
        }

        double* const out = stationOut + strainOffset;
        out[0] = exx;
        out[1] = eyy;
        out[2] = ezz;
        out[3] = exy;
        out[4] = eyz;
        out[5] = exz;
    }

    if (stressOffset >= 0) {
        const int32_t numcomp = lenum * 8;
        const double  mu = gpuData->stationGatherMuDevice[lin];
        const double  lambda = gpuData->stationGatherLambdaDevice[lin];

        /* Gather 8 cube nodes x 3 components of shearVector and
         * kappaVector at the station's owning element. */
        double sx[8], sy[8], sz[8], kx[8], ky[8], kz[8];
        const solver_float* sh =
            (const solver_float*)gpuData->shearVectorDevice;
        const solver_float* ka =
            (const solver_float*)gpuData->kappaVectorDevice;
        #pragma unroll
        for (int i = 0; i < 8; i++) {
            const int32_t o0 = eindex + i * lenum;
            const int32_t o1 = o0 + numcomp;
            const int32_t o2 = o1 + numcomp;
            sx[i] = (double)sh[o0];
            sy[i] = (double)sh[o1];
            sz[i] = (double)sh[o2];
            kx[i] = (double)ka[o0];
            ky[i] = (double)ka[o1];
            kz[i] = (double)ka[o2];
        }

        double esxx, esyy, eszz, esxy, esyz, esxz;
        double ekxx, ekyy, ekzz, ekxy, ekyz, ekxz;
        if (teth_i < 0) {
            v3_cube_strain(sx, sy, sz, csi, eta, dze, h,
                           &esxx, &esyy, &eszz, &esxy, &esyz, &esxz);
            v3_cube_strain(kx, ky, kz, csi, eta, dze, h,
                           &ekxx, &ekyy, &ekzz, &ekxy, &ekyz, &ekxz);
        } else {
            v3_tet_strain(sx, sy, sz, h, teth_i, cp,
                          &esxx, &esyy, &eszz, &esxy, &esyz, &esxz);
            v3_tet_strain(kx, ky, kz, h, teth_i, cp,
                          &ekxx, &ekyy, &ekzz, &ekxy, &ekyz, &ekxz);
        }

        /* Viscoelastic stress (deviatoric/volumetric split, consistent
         * with Hercules' damping force operator). */
        const double tr_s = esxx + esyy + eszz;
        const double tr_k = ekxx + ekyy + ekzz;
        const double K_bulk = lambda + (2.0/3.0) * mu;
        const double Ktr_k = K_bulk * tr_k;
        const double dev_off = (1.0/3.0) * tr_s;
        const double two_mu = 2.0 * mu;

        double* const out = stationOut + stressOffset;
        out[0] = two_mu * (esxx - dev_off) + Ktr_k;
        out[1] = two_mu * (esyy - dev_off) + Ktr_k;
        out[2] = two_mu * (eszz - dev_off) + Ktr_k;
        out[3] = two_mu * esxy;
        out[4] = two_mu * esyz;
        out[5] = two_mu * esxz;
    }
}

/* Interpolate acoustic pressure and pressure gradient directly from the two
 * device-resident history fields. The shape-function algebra and component
 * ordering mirror interpolate_station_pressure_value and
 * station_pressure_gradient in psolve.cc. */
__global__ void kernelAcousticStationInterpolate(
    int32_t stationCount, const int32_t* stationNodes,
    const double* stationCoords, const double* stationH,
    const solver_float* pressureTm1, const solver_float* pressureTm2,
    int32_t velocityPotential, double inverseDt, double* stationOutput)
{
    const int32_t station = blockIdx.x * blockDim.x + threadIdx.x;
    if (station >= stationCount) return;

    const double sx[8] = {-1, 1, -1, 1, -1, 1, -1, 1};
    const double sy[8] = {-1, -1, 1, 1, -1, -1, 1, 1};
    const double sz[8] = {-1, -1, -1, -1, 1, 1, 1, 1};
    const double xi = stationCoords[(size_t)station * 3];
    const double eta = stationCoords[(size_t)station * 3 + 1];
    const double zeta = stationCoords[(size_t)station * 3 + 2];
    const double inv4h = 0.25 / stationH[station];
    double pressure = 0.0;
    double gradX = 0.0;
    double gradY = 0.0;
    double gradZ = 0.0;

    #pragma unroll
    for (int corner = 0; corner < 8; corner++) {
        const int32_t node = stationNodes[(size_t)station * 8 + corner];
        const solver_float currentField = pressureTm1[node];
        const double current = (double)currentField;
        const double value =
            velocityPotential
                ? -(double)(currentField - pressureTm2[node]) * inverseDt
                : current;
        const double oneX = 1.0 + sx[corner] * xi;
        const double oneY = 1.0 + sy[corner] * eta;
        const double oneZ = 1.0 + sz[corner] * zeta;
        pressure += 0.125 * oneX * oneY * oneZ * value;
        gradX += sx[corner] * oneY * oneZ * inv4h * value;
        gradY += sy[corner] * oneX * oneZ * inv4h * value;
        gradZ += sz[corner] * oneX * oneY * inv4h * value;
    }

    double* const out = stationOutput + (size_t)station * 4;
    out[0] = pressure;
    out[1] = gradX;
    out[2] = gradY;
    out[3] = gradZ;
}

/* Match the host poroelastic station interpolation exactly: inactive nodal
 * pressure masses are skipped and the surviving trilinear weights are
 * renormalized. */
__global__ void kernelPorePressureStationInterpolate(
    int32_t stationCount, const int32_t* stationNodes,
    const double* stationCoords, const solver_float* porePressureTm1,
    const solver_float* porePressureMass, double* stationOutput)
{
    const int32_t station = blockIdx.x * blockDim.x + threadIdx.x;
    if (station >= stationCount) return;

    const double sx[8] = {-1, 1, -1, 1, -1, 1, -1, 1};
    const double sy[8] = {-1, -1, 1, 1, -1, -1, 1, 1};
    const double sz[8] = {-1, -1, -1, -1, 1, 1, 1, 1};
    const double xi = stationCoords[(size_t)station * 3];
    const double eta = stationCoords[(size_t)station * 3 + 1];
    const double zeta = stationCoords[(size_t)station * 3 + 2];
    double pressure = 0.0;
    double weight = 0.0;

    #pragma unroll
    for (int corner = 0; corner < 8; corner++) {
        const int32_t node = stationNodes[(size_t)station * 8 + corner];
        if ((double)porePressureMass[node] > 0.0) {
            const double shape =
                0.125 * (1.0 + sx[corner] * xi) *
                (1.0 + sy[corner] * eta) *
                (1.0 + sz[corner] * zeta);
            pressure += shape * (double)porePressureTm1[node];
            weight += shape;
        }
    }

    stationOutput[station] =
        (weight > 0.0) ? (pressure / weight) : nan("");
}

/* Interpolate one XYZ displacement per locally owned plane point. The output
 * order exactly matches io_planes.cc's plane/strip/point traversal, allowing
 * it to be copied once and written through precomputed MPI file views. */
__global__ void kernelPlaneDisplacementInterpolate(
    int32_t pointCount, const int32_t* pointNodes,
    const double* pointCoords, const int32_t* pointModes,
    const fvector_t* tm1, double* pointOutput)
{
    const int32_t point = blockIdx.x * blockDim.x + threadIdx.x;
    if (point >= pointCount) return;

    const int32_t mode = pointModes[point];
    double* const out = pointOutput + (size_t)point * 3;
    if ((mode == PLANE_GPU_POINT_INVALID) ||
        (mode == PLANE_GPU_POINT_HOST_ACOUSTIC)) {
        out[0] = 0.0;
        out[1] = 0.0;
        out[2] = 0.0;
        return;
    }

    const int32_t* const nodes = pointNodes + (size_t)point * 8;
    const double xi = pointCoords[(size_t)point * 3];
    const double eta = pointCoords[(size_t)point * 3 + 1];
    const double zeta = pointCoords[(size_t)point * 3 + 2];

    if (mode == PLANE_GPU_POINT_TOPO) {
        const fvector_t u0 = tm1[nodes[0]];
        const fvector_t u1 = tm1[nodes[1]];
        const fvector_t u2 = tm1[nodes[2]];
        const fvector_t u3 = tm1[nodes[3]];
        #pragma unroll
        for (int axis = 0; axis < 3; axis++) {
            const double base = (double)u0.f[axis];
            out[axis] = base + xi * ((double)u2.f[axis] - base) +
                        eta * ((double)u1.f[axis] - base) +
                        zeta * ((double)u3.f[axis] - base);
        }
        return;
    }

    const double sx[8] = {-1, 1, -1, 1, -1, 1, -1, 1};
    const double sy[8] = {-1, -1, 1, 1, -1, -1, 1, 1};
    const double sz[8] = {-1, -1, -1, -1, 1, 1, 1, 1};
    double ux = 0.0;
    double uy = 0.0;
    double uz = 0.0;
    #pragma unroll
    for (int corner = 0; corner < 8; corner++) {
        const double shape =
            0.125 * (1.0 + sx[corner] * xi) *
            (1.0 + sy[corner] * eta) *
            (1.0 + sz[corner] * zeta);
        const fvector_t value = tm1[nodes[corner]];
        ux += shape * (double)value.f[0];
        uy += shape * (double)value.f[1];
        uz += shape * (double)value.f[2];
    }
    out[0] = ux;
    out[1] = uy;
    out[2] = uz;
}
