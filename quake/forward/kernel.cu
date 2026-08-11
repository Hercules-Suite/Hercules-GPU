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

    /* CPML elastic PML stiffness mode (gpuData->pmlNativeStiffness):
     *  - native (1, default): do NOT skip -- PML elements get the FULL
     *    native effective stiffness here, byte-identical to interior;
     *    kernelPMLElasticStep then adds only the CPML correction
     *    -(sigma_pml - sigma_raw) + CFS (zero where d=0, so the
     *    PML<->interior operator is continuous and PML elements keep the
     *    native higher / hourglass-resisting modes).
     *  - reduced (0): SKIP PML elements here; kernelPMLElasticStep then
     *    supplies the full standalone 1-point reduced CPML stiffness.
     * pmlFlagArrayDevice is NULL in dashpot mode (no skip either way),
     * and is still used unconditionally by the acoustic-pressure and BKT
     * kernels (their skips are unchanged). */
    if (gpuData->pmlFlagArrayDevice != NULL &&
        gpuData->pmlFlagArrayDevice[eindex] != 0 &&
        gpuData->pmlNativeStiffness == 0) {
      return;
    }

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
     * For CPML, b*K is skipped on PML elements (5d91b09 material-damping-
     * on-PML policy) UNLESS the viscoelastic Rayleigh PML is requested
     * (cpmlPmlMaterialDamping==2) -- then PML gets the full native b*K
     * and kernelPMLElasticStep adds the CPML stretch delta. */
    double sd_b_over_dt = 0.0;
    {
      const double c1v = gpuData->c1ArrayDevice[eindex];
      if (c1v != 0.0) sd_b_over_dt = gpuData->c3ArrayDevice[eindex] / c1v;
      if (sd_b_over_dt != 0.0 &&
          gpuData->pmlFlagArrayDevice != NULL &&
          gpuData->pmlFlagArrayDevice[eindex] != 0 &&
          gpuData->cpmlPmlMaterialDamping != 2) {
        sd_b_over_dt = 0.0;
      }
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
  const int32_t force_eindex = (adjacent_eindex >= 0) ? adjacent_eindex : eindex;
  const double* force_center = (adjacent_eindex >= 0) ? adjacent_center : center;
  const int pml_face =
      (gpuData->acousticInterfacePmlSideMaskDevice != NULL) &&
      (gpuData->acousticInterfacePmlSideMaskDevice[faceIndex] != 0);

  if ((gpuData->acousticInterfaceProjectionMode ==
       ACOUSTIC_INTERFACE_PROJECTION_FACE_MASS) &&
      (adjacent_eindex >= 0) && !pml_face &&
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

        for (int i = 0; i < 8; i++) {
          const double weight =
              acoustic_shape_value_at(i, q_force_center[0],
                                      q_force_center[1], q_force_center[2]);
          if (weight == 0.0) {
            continue;
          }

          const int32_t lnid =
              gpuData->lnidArrayDevice[force_eindex + i * gpuData->lenum];
          fvector_t* nodalForce = gpuData->forceDevice + lnid;
          atomicAdd((solver_float*)&nodalForce->f[0],
                    (solver_float)(q_coeff * weight * q_pressure *
                                   normal[0]));
          atomicAdd((solver_float*)&nodalForce->f[1],
                    (solver_float)(q_coeff * weight * q_pressure *
                                   normal[1]));
          atomicAdd((solver_float*)&nodalForce->f[2],
                    (solver_float)(q_coeff * weight * q_pressure *
                                   normal[2]));
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

  /* CPML-aware interface correction for the pressure traction in the
   * elastic equation. Use the elastic-side stretch only; using the
   * acoustic-side profile here would apply water Vp/sigma coefficients to
   * the solid PML interface flux. Identity on non-elastic-PML faces and
   * at the PML inner edge. */
  if (gpuData->acousticInterfacePmlSideMaskDevice != NULL &&
      (gpuData->acousticInterfacePmlSideMaskDevice[faceIndex] &
       ACOUSTIC_INTERFACE_PML_ELASTIC_SIDE) != 0) {
    const double a =
        gpuData->acousticInterfacePmlElastic_aDevice[faceIndex];
    const double b =
        gpuData->acousticInterfacePmlElastic_bDevice[faceIndex];
    const double iK =
        gpuData->acousticInterfacePmlElastic_invKappaDevice[faceIndex];
    const double psi_old =
        (double)gpuData->acousticInterfacePsiPressureDevice[faceIndex];
    const double psi_new = b * psi_old + a * pressure;
    gpuData->acousticInterfacePsiPressureDevice[faceIndex] =
        (solver_float)psi_new;
    pressure = iK * pressure + psi_new;
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

  /* CPML acoustic PML stiffness mode (mirror of the validated elastic
   * native+correction in kernelStiffnessCalcLocal):
   *  - native (1, default): do NOT skip -- PML acoustic elements get the
   *    FULL consistent operator here, byte-identical to the interior (the
   *    operator validated against CPS, commit 20180b6).
   *    kernelPMLAcousticStep then adds only the CPML delta
   *    -(P_pml - P_raw) + CFS, which is exactly zero where d=0, so the
   *    PML<->interior operator is continuous and PML elements keep the
   *    consistent, hourglass-free modes (the bare 1-pt reduced operator
   *    is rank-deficient and rings undamped in the lossless fluid).
   *  - reduced (0): SKIP here; kernelPMLAcousticStep then supplies the
   *    full standalone 1-point reduced CPML stiffness (pre-restructure
   *    behaviour, bit-identical). */
  if (gpuData->pmlFlagArrayDevice != NULL &&
      gpuData->pmlFlagArrayDevice[eindex] != 0 &&
      gpuData->pmlNativeStiffnessAcoustic == 0) {
    return;
  }

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
  const double* accel_center =
      (adjacent_eindex >= 0) ? adjacent_center : center;
  const int pml_face =
      (gpuData->acousticInterfacePmlSideMaskDevice != NULL) &&
      (gpuData->acousticInterfacePmlSideMaskDevice[faceIndex] != 0);

  if ((gpuData->acousticInterfaceProjectionMode ==
       ACOUSTIC_INTERFACE_PROJECTION_FACE_MASS) &&
      (adjacent_eindex >= 0) && !pml_face &&
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
        double q_source_n = 0.0;

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
              gpuData->lnidArrayDevice[accel_eindex + i * gpuData->lenum];
          const fvector_t u_next = gpuData->tm2Device[lnid];
          const fvector_t u_curr = gpuData->tm1Device[lnid];
          if (vp) {
            q_source_n +=
                weight *
                (((double)u_next.f[0] - (double)u_curr.f[0]) * invDt *
                     normal[0] +
                 ((double)u_next.f[1] - (double)u_curr.f[1]) * invDt *
                     normal[1] +
                 ((double)u_next.f[2] - (double)u_curr.f[2]) * invDt *
                     normal[2]);
          } else {
            const fvector_t u_prev = gpuData->tm3Device[lnid];
            q_source_n +=
                weight *
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
        }

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
  double source_n = 0.0;

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
    if (vp) {
      source_n += weight *
                  (((double)u_next.f[0] - (double)u_curr.f[0]) * invDt *
                       normal[0] +
                   ((double)u_next.f[1] - (double)u_curr.f[1]) * invDt *
                       normal[1] +
                   ((double)u_next.f[2] - (double)u_curr.f[2]) * invDt *
                       normal[2]);
    } else {
      const fvector_t u_prev = gpuData->tm3Device[lnid];
      source_n += weight *
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
  }

  /* CPML-aware interface correction for the elastic normal-motion source in
   * the acoustic equation. Use the acoustic-side stretch only; the acoustic
   * RHS is the one being modified. Identity on non-acoustic-PML faces / at
   * the PML inner edge. */
  if (gpuData->acousticInterfacePmlSideMaskDevice != NULL &&
      (gpuData->acousticInterfacePmlSideMaskDevice[faceIndex] &
       ACOUSTIC_INTERFACE_PML_ACOUSTIC_SIDE) != 0) {
    const double a =
        gpuData->acousticInterfacePmlAcoustic_aDevice[faceIndex];
    const double b =
        gpuData->acousticInterfacePmlAcoustic_bDevice[faceIndex];
    const double iK =
        gpuData->acousticInterfacePmlAcoustic_invKappaDevice[faceIndex];
    const double psi_old =
        (double)gpuData->acousticInterfacePsiAccelDevice[faceIndex];
    const double psi_new = b * psi_old + a * source_n;
    gpuData->acousticInterfacePsiAccelDevice[faceIndex] =
        (solver_float)psi_new;
    source_n = iK * source_n + psi_new;
  }

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
      (gpuData->pressureFreeSurfaceNodeFlagsDevice[nindex] != 0)) {
    /* Pressure-release free surface, p = 0. Pressure mode: field=p, set
     * 0. velocity_potential: p = -dphi/dt = 0  ->  phi^{n+1} = phi^n. */
    gpuData->pressureTm2Device[nindex] =
        (gpuData->acousticVelocityPotential == 1)
            ? gpuData->pressureTm1Device[nindex]
            : (solver_float)0.0;
    gpuData->pressureForceDevice[nindex] = 0.0;
    return;
  }

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
   * kernelDispCalc zeros u there, and we must likewise zero the acoustic
   * pressure (the primary field in the fluid) so the outgoing wave is
   * terminated single-pass instead of reflecting off the outer face.
   * Mirrors the elastic clamp exactly; no-op when cpml_clamp_outer=no
   * (mask NULL). */
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
      (adjacent_eindex >= 0) ? adjacent_center : center;
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
        gpuData->lnidArrayDevice[target_eindex + j * gpuData->lenum];
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
     * element. Force u=0 at those nodes so the wave cannot reflect off
     * the outer free surface of the PML zone. z_min (free surface) is
     * never clamped. */
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

    /* === CPML stiffness CORRECTION (native + delta) ===
     * The FULL native effective stiffness for this element is added
     * separately by kernelStiffnessCalcLocal (its PML skip was removed),
     * byte-identical to the interior. Here we add only the CPML delta
     *   -(sigma_pml - sigma_raw)
     * evaluated at the 1-point element centre in the SAME native
     * normalization (mu = 9*c1, lambda = 9*c2, dN = 0.25*sign), so
     * sigma_raw is exactly the constant-strain projection the native
     * operator contains. Where d=0 the psi memory stays 0 and
     * sigma_pml == sigma_raw, so the delta is exactly zero: no operator
     * discontinuity at the PML inner edge, and PML elements keep the
     * native higher (hourglass-resisting) modes. This mirrors the codex
     * v3 reduced-correction (gradient-side psi only -- no stress
     * convolution). */
    const double mu     = 9.0 * gpuData->c1ArrayDevice[eindex];
    const double lambda = 9.0 * gpuData->c2ArrayDevice[eindex];

    /* Natural-coord gradient at the centre: dN_n,axis = 0.25 * sgn.
     * grad[i][j] = sum_n u_i^(n) * dN_n,j  (= du_i/dx_j in this norm). */
    double grad[3][3] = {{0.0, 0.0, 0.0},
                         {0.0, 0.0, 0.0},
                         {0.0, 0.0, 0.0}};
    for (int n = 0; n < 8; n++) {
        const double gN[3] = { (n & 1) ? 0.25 : -0.25,
                               (n & 2) ? 0.25 : -0.25,
                               (n & 4) ? 0.25 : -0.25 };
        for (int i = 0; i < 3; i++) {
            grad[i][0] += u[n][i] * gN[0];
            grad[i][1] += u[n][i] * gN[1];
            grad[i][2] += u[n][i] * gN[2];
        }
    }

    /* Raw (un-stretched) stress. */
    const double traceRaw = grad[0][0] + grad[1][1] + grad[2][2];
    double sigmaRaw[3][3];
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            sigmaRaw[a][b] = mu * (grad[a][b] + grad[b][a]);
        }
        sigmaRaw[a][a] += lambda * traceRaw;
    }

    /* Gradient-side CPML psi memory (raw a / 1/kappa; b = exp decay).
     * psi[i][j] convolves the j-axis gradient of u_i. */
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
    double gradPml[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            const double psi_old =
                (double)gpuData->pml_elastic_psiArrayDevice[i][j][pml_idx];
            const double psi_new = bJ[j] * psi_old + aJ[j] * grad[i][j];
            gpuData->pml_elastic_psiArrayDevice[i][j][pml_idx] =
                (solver_float)psi_new;
            gradPml[i][j] = invKJ[j] * grad[i][j] + psi_new;
        }
    }
    const double tracePml = gradPml[0][0] + gradPml[1][1] + gradPml[2][2];
    double sigmaPml[3][3];
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            sigmaPml[a][b] = mu * (gradPml[a][b] + gradPml[b][a]);
        }
        sigmaPml[a][a] += lambda * tracePml;
    }

    /* BKT-family material damping is intentionally not stretched inside the
     * elastic CPML. The PML zone uses native elastic stiffness plus this
     * CPML correction; BKT memory/force kernels skip PML elements. */
    /* Rayleigh stiffness-proportional (b*K) viscoelastic CPML, gradient
     * side. Single effective displacement:
     *   u_eff = (1 + b/dt) u^n - (b/dt) u^{n-1},   b/dt = c3/c1
     * relaxed constant-strain stress is the plain elastic form on u_eff
     *   sigma_re = mu*(g_e+g_e^T) + lambda*tr(g_e)*I,  g_e = grad(u_eff)
     * which reduces EXACTLY to the elastic stress when b/dt==0
     * (u_eff==u^n), so the d=0 delta cancellation against the native
     * (now b*K-inclusive, kernelStiffnessCalcLocal) modal force holds ->
     * inner edge reflectionless (v10-safe). Reuses the psiShear memory
     * family.
     * b/dt==0 (any non-Rayleigh, or zeta=0) -> block skipped, sigma
     * stays elastic -> byte-identical. */
    if (gpuData->cpmlPmlMaterialDamping == 2) {
        double bod = 0.0;
        const double c1v = gpuData->c1ArrayDevice[eindex];
        if (c1v != 0.0) bod = gpuData->c3ArrayDevice[eindex] / c1v;
        if (bod != 0.0) {
            double gE[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
            for (int n = 0; n < 8; n++) {
                const int32_t ln = lnid_base[n * lenum];
                const double ue[3] = {
                    (double)gpuData->tm1Device[ln].f[0]*(1.0+bod) -
                        bod*(double)gpuData->tm2Device[ln].f[0],
                    (double)gpuData->tm1Device[ln].f[1]*(1.0+bod) -
                        bod*(double)gpuData->tm2Device[ln].f[1],
                    (double)gpuData->tm1Device[ln].f[2]*(1.0+bod) -
                        bod*(double)gpuData->tm2Device[ln].f[2] };
                const double gN[3] = { (n & 1) ? 0.25 : -0.25,
                                       (n & 2) ? 0.25 : -0.25,
                                       (n & 4) ? 0.25 : -0.25 };
                for (int i = 0; i < 3; i++) {
                    gE[i][0] += ue[i]*gN[0]; gE[i][1] += ue[i]*gN[1];
                    gE[i][2] += ue[i]*gN[2];
                }
            }
            const double trEr = gE[0][0] + gE[1][1] + gE[2][2];
            for (int a = 0; a < 3; a++) {
                for (int c = 0; c < 3; c++)
                    sigmaRaw[a][c] = mu * (gE[a][c] + gE[c][a]);
                sigmaRaw[a][a] += lambda * trEr;
            }
            double gEp[3][3];
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    const double p_old = (double)
                        gpuData->pml_elastic_psiShearArrayDevice[i][j][pml_idx];
                    const double p_new = bJ[j]*p_old + aJ[j]*gE[i][j];
                    gpuData->pml_elastic_psiShearArrayDevice[i][j][pml_idx] =
                        (solver_float)p_new;
                    gEp[i][j] = invKJ[j]*gE[i][j] + p_new;
                }
            }
            const double trEp = gEp[0][0] + gEp[1][1] + gEp[2][2];
            for (int a = 0; a < 3; a++) {
                for (int c = 0; c < 3; c++)
                    sigmaPml[a][c] = mu * (gEp[a][c] + gEp[c][a]);
                sigmaPml[a][a] += lambda * trEp;
            }
        }
    }

    /* Divergence-side CPML, correct 2nd-order displacement form (Matzen
     * 2011; Xie et al. 2014, the SPECFEM3D FE-PML). Multiplying the
     * momentum equation by s_j turns the divergence-side stretch into,
     * per damping axis a:
     *   rho*d_a du_i/dt           -> mass-prop. damping (mass arrays)
     *   +rho*d_a*alpha_a u_i      -> CFS restoring     (HERE)
     *   -rho*d_a*alpha_a^2 zeta_i -> CFS memory force  (HERE)
     * The damping is injected at the node in solver_init (validated,
     * stable in v9). The CFS restoring + memory are a pair that must
     * cancel at DC; v10 applied the restoring per-node via mass2_minusaM
     * while zeta was element-centroid based -- the discrete cancellation
     * failed and PML corners grew exponentially. Both halves are now
     * evaluated HERE from the SAME element-centroid displacement
     *   ubar_i = (1/8) sum_n u_i^(n)   (exact for trilinear hex8)
     * so the per-element net CFS force
     *   cfs_i = restoreCoeff*ubar_i - sum_a memForce_a*zeta_i^(a)
     * is an exact discrete null at DC (zeta_i^(a) -> ubar_i/alpha_a =>
     * sum_a memForce_a*zeta_i -> restoreCoeff*ubar_i). zeta advances by
     *   zeta_i^(a) = zetaDecay_a*zeta_i^(a) + zetaGain_a*ubar_i
     * restoreCoeff = dt^2*(rho h^3/8)*sum_a d_a*alpha_a and
     * memForce_a   = dt^2*(rho h^3/8)*d_a*alpha_a^2 are PER-NODE shares,
     * so the same cfs_i is added to each of the 8 nodes. Strict no-op
     * where d_a=0 (restoreCoeff=0, memForce_a=0). zeta state reuses
     * pml_elastic_divPsiArrayDevice[i][a]. */
    double ubar[3] = {0.0, 0.0, 0.0};
    for (int n = 0; n < 8; n++) {
        ubar[0] += u[n][0];
        ubar[1] += u[n][1];
        ubar[2] += u[n][2];
    }
    ubar[0] *= 0.125; ubar[1] *= 0.125; ubar[2] *= 0.125;

    const double zDec[3] = {
        gpuData->pml_elastic_zetaDecayArrayDevice[0][pml_idx],
        gpuData->pml_elastic_zetaDecayArrayDevice[1][pml_idx],
        gpuData->pml_elastic_zetaDecayArrayDevice[2][pml_idx] };
    const double zGain[3] = {
        gpuData->pml_elastic_zetaGainArrayDevice[0][pml_idx],
        gpuData->pml_elastic_zetaGainArrayDevice[1][pml_idx],
        gpuData->pml_elastic_zetaGainArrayDevice[2][pml_idx] };
    const double memF[3] = {
        gpuData->pml_elastic_memForceArrayDevice[0][pml_idx],
        gpuData->pml_elastic_memForceArrayDevice[1][pml_idx],
        gpuData->pml_elastic_memForceArrayDevice[2][pml_idx] };
    const double restoreCoeff =
        gpuData->pml_elastic_restoreCoeffArrayDevice[pml_idx];

    /* Advance zeta_i^(a); net per-node CFS force
     * cfs_i = restoreCoeff*ubar_i - sum_a memForce_a * zeta_i^(a). */
    double cfs[3];
    for (int i = 0; i < 3; i++) {
        double fmem = 0.0;
        for (int a = 0; a < 3; a++) {
            const double z_old =
                (double)gpuData->pml_elastic_divPsiArrayDevice[i][a][pml_idx];
            const double z_new = zDec[a] * z_old + zGain[a] * ubar[i];
            gpuData->pml_elastic_divPsiArrayDevice[i][a][pml_idx] =
                (solver_float)z_new;
            fmem += memF[a] * z_new;
        }
        cfs[i] = restoreCoeff * ubar[i] - fmem;
    }

    /* Assemble (centre derivatives dN_n,axis = 0.25*sgn) + CFS pair.
     * Mode (gpuData->pmlNativeStiffness):
     *  - native (1): add the CPML delta -(nodalPml - nodalRaw); the FULL
     *    native stiffness was added by kernelStiffnessCalcLocal. At d=0
     *    sigmaPml==sigmaRaw and cfs==0, so only the native operator
     *    remains -> exact PML<->interior match.
     *  - reduced (0): kernelStiffnessCalcLocal skipped this element, so
     *    add the FULL standalone 1-point reduced CPML stiffness -nodalPml
     *    (no native higher modes; the pre-a1c0fb8 behaviour). */
    const int native = gpuData->pmlNativeStiffness;
    for (int n = 0; n < 8; n++) {
        const double gN[3] = { (n & 1) ? 0.25 : -0.25,
                               (n & 2) ? 0.25 : -0.25,
                               (n & 4) ? 0.25 : -0.25 };
        const int32_t lnid = lnid_base[n * lenum];
        fvector_t* nodalForce = gpuData->forceDevice + lnid;
        for (int a = 0; a < 3; a++) {
            const double nodalPml = gN[0] * sigmaPml[a][0] +
                                    gN[1] * sigmaPml[a][1] +
                                    gN[2] * sigmaPml[a][2];
            double f;
            if (native) {
                const double nodalRaw = gN[0] * sigmaRaw[a][0] +
                                        gN[1] * sigmaRaw[a][1] +
                                        gN[2] * sigmaRaw[a][2];
                f = -(nodalPml - nodalRaw) + cfs[a];
            } else {
                f = -nodalPml + cfs[a];
            }
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

    /* Step B -- divergence-side CFS (scalar mirror of the elastic
     * kernelPMLElasticStep CFS pair). The gradient side is partial[j];
     * the divergence-side mass-prop damping is applied to the nodal
     * pressure mass (pressure_damp), and the CFS restoring + memory pair
     * is evaluated here from the element-centroid pressure
     *   pbar = (1/8) sum_n p^(n)
     *   cfs  = restoreCoeff*pbar - sum_a memForce_a * zeta_p^(a)
     *   zeta_p^(a) = zetaDecay_a*zeta_p^(a) + zetaGain_a*pbar
     * restoreCoeff/memForce already pack dt^2 * nodalPressureMass *
     * (d*alpha, d*alpha^2), i.e. the PER-NODE share, so the same cfs is
     * added to each of the 8 nodes. At d=0 (inner edge / non-PML axis)
     * restoreCoeff = memForce = 0 -> cfs = 0 (strict no-op), and the
     * pair cancels at DC (zeta_p -> pbar/alpha). */
    double pbar = 0.0;
    for (int n = 0; n < 8; n++) pbar += p[n];
    pbar *= 0.125;

    const double zDec[3] = {
        gpuData->pml_acoustic_zetaDecayArrayDevice[0][pml_idx],
        gpuData->pml_acoustic_zetaDecayArrayDevice[1][pml_idx],
        gpuData->pml_acoustic_zetaDecayArrayDevice[2][pml_idx] };
    const double zGain[3] = {
        gpuData->pml_acoustic_zetaGainArrayDevice[0][pml_idx],
        gpuData->pml_acoustic_zetaGainArrayDevice[1][pml_idx],
        gpuData->pml_acoustic_zetaGainArrayDevice[2][pml_idx] };
    const double memF[3] = {
        gpuData->pml_acoustic_memForceArrayDevice[0][pml_idx],
        gpuData->pml_acoustic_memForceArrayDevice[1][pml_idx],
        gpuData->pml_acoustic_memForceArrayDevice[2][pml_idx] };
    const double restoreCoeff =
        gpuData->pml_acoustic_restoreCoeffArrayDevice[pml_idx];

    double fmem = 0.0;
    for (int a = 0; a < 3; a++) {
        const double z_old =
            (double)gpuData->pml_acoustic_zetaPArrayDevice[a][pml_idx];
        const double z_new = zDec[a] * z_old + zGain[a] * pbar;
        gpuData->pml_acoustic_zetaPArrayDevice[a][pml_idx] =
            (solver_float)z_new;
        fmem += memF[a] * z_new;
    }
    const double cfs = restoreCoeff * pbar - fmem;

    /* Pressure force -- native+correction mirror of kernelPMLElasticStep.
     * Cstiff = dt^2 * h / (16 rho) is the corrected consistent interior
     * acoustic stiffness coefficient (commit 20180b6), stored (raw a/b/
     * invKappa, natural signed-sum normalization) exactly as the elastic
     * path stores raw a/b/invKappa and folds the physical scaling into
     * mu=9*c1. With  S[j] = sum_n s_j(n) p_n  (constant-gradient signed
     * sum) and  partial[j] = invKappa*S[j] + psi  (CPML-stretched):
     *   P_raw(n) = Cstiff * sum_j s_j(n) * S[j]        -- the constant-
     *     gradient projection, exactly the content the consistent
     *     interior operator already carries (1-pt Gauss is exact on the
     *     linear subspace), and
     *   P_pml(n) = Cstiff * sum_j s_j(n) * partial[j]  -- CPML-stretched.
     * native (default): kernelAcousticPressureForce added the FULL
     *   consistent operator for this element; here add only the CPML
     *   delta  -(P_pml - P_raw) + cfs.  At d=0 partial==S so the delta
     *   is exactly 0 and cfs==0 -> the element is purely the consistent
     *   interior operator (validated, hourglass-free, inner-edge exact).
     * reduced: interior skipped this element -> add the full standalone
     *   1-point reduced CPML stiffness  -P_pml + cfs  (bit-identical to
     *   the pre-restructure kernel; the raw/Cstiff rescaling is an exact
     *   algebraic identity with the old inv4h / dt^2 h^2/4 form). */
    const double Cstiff =
        gpuData->pml_acoustic_forceCoeffArrayDevice[pml_idx];
    const int native = gpuData->pmlNativeStiffnessAcoustic;
    for (int n = 0; n < 8; n++) {
        const double sx = (n & 1) ? 1.0 : -1.0;
        const double sy = (n & 2) ? 1.0 : -1.0;
        const double sz = (n & 4) ? 1.0 : -1.0;
        const int32_t lnid = lnid_base[n * lenum];
        const double Ppml = Cstiff *
            (sx * partial[0] + sy * partial[1] + sz * partial[2]);
        double f;
        if (native) {
            const double Praw = Cstiff *
                (sx * S[0] + sy * S[1] + sz * S[2]);
            f = -(Ppml - Praw) + cfs;
        } else {
            f = -Ppml + cfs;
        }
        atomicAdd((solver_float*)(gpuData->pressureForceDevice + lnid),
                  (solver_float)f);
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
