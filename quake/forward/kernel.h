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

#ifndef KERNEL_H_
#define KERNEL_H_

#include "psolve.h"
#include "quake_util.h"
#include <cuda_runtime.h>


/* Kernel identifiers for FLOP counts */
#define FLOP_MAX_KERNEL        4
#define FLOP_STIFFNESS_KERNEL  0
#define FLOP_CALCCONV_KERNEL   1
#define FLOP_DAMPING_KERNEL    2
#define FLOP_DISPLACEMENT      3


/* Debug and utility functions */
int64_t kernel_flops_per_thread(int kid);
int64_t kernel_mem_per_thread(int kid);
int32_t gpu_get_reg_count(char* kernel);
int32_t gpu_get_blocksize(gpu_spec_t *gpuSpecs, char* kernel,
			  int32_t memPerThread);
int gpu_copy_constant_symbols(int myID, fmatrix_t (*theK1)[8],
			      fmatrix_t (*theK2)[8]);

/* Kernels */
__global__  void kernelStiffnessCalcLocal(int32_t lenum,
					  gpu_data_t *gpuData,
					  int32_t   myLinearElementsCount,
					  int32_t*  myLinearElementsMapperDevice);

__global__  void kernelDampingCalcConv(int32_t lenum,
				       int32_t numcomp,
				       gpu_data_t *gpuData,
				       double rmax,
				       int32_t dampingType,
				       int32_t dampingElemCount,
				       int32_t* dampingElemMapperDevice);

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
					   int32_t* dampingElemMapperDevice);

__global__  void kernelDampingCalcConvBkt3(int32_t lenum,
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
	int32_t* dampingElemMapperDevice);

__global__  void kernelDampingCalcLocal(int32_t lenum,
					gpu_data_t *gpuData,
					int32_t dampingElemCount,
					int32_t* dampingElemMapperDevice);

__global__  void kernelDispCalc(int32_t nharbored,
					gpu_data_t* gpuData,
	                                noyesflag_t printAccel);

/* v3 station viscoelastic stress gather: one thread per station; reads
 * the 8-node shearVector/kappaVector slot of the station's owning
 * element, computes the deviatoric (shear) and volumetric (kappa)
 * modified strains via inline trilinear-hex or VT-tet shape-function
 * derivatives, and writes [eps_s_xx..xz, sig_visco_xx..xz] (12 doubles
 * per station) into gpuData->stationVEStressArrayDevice. Launched only
 * when v3 is active (print_station_strains||print_station_stresses,
 * type_of_damping >= BKT). */
__global__ void kernelStationVEStressGather(int32_t stationCount,
                                            gpu_data_t* gpuData);

__global__ void kernelAcousticPressureDisplacementCoupling(int32_t faceCount,
							   gpu_data_t* gpuData,
							   double deltaT2,
								   double deltaT);

__global__ void kernelAcousticPressureForce(int32_t lenum,
					    gpu_data_t* gpuData,
					    double deltaT2);

__global__ void kernelAcousticPressureInterfaceForce(int32_t faceCount,
						     gpu_data_t* gpuData,
						     double deltaT2,
								   double deltaT);

__global__ void kernelAcousticPressureUpdate(int32_t nharbored,
					     gpu_data_t* gpuData);

__global__ void kernelPoroelasticPressureTraction(int32_t elemCount,
                                                  gpu_data_t* gpuData,
                                                  double deltaT2,
                                                  double biotAlpha);

__global__ void kernelPoroelasticPressureForce(int32_t elemCount,
                                               gpu_data_t* gpuData,
                                               double deltaT,
                                               double mobilityX,
                                               double mobilityY,
                                               double mobilityZ,
                                               int haveDiffusion,
                                               double biotAlpha);

__global__ void kernelPoroelasticPressureUpdate(int32_t nharbored,
                                                gpu_data_t* gpuData);

__global__  void kernelAlignTm2Disp(int32_t lenum,
				    int32_t numcomp,
				    gpu_data_t *gpuData);

/* C-PML kernels (unsplit ADE form). Each kernel processes the PML
 * elements of one material type (elastic or acoustic) via its dedicated
 * mapper. The kernel fuses the ADE memory-variable update and the
 * modified-strain/gradient force assembly so the per-element nodal data
 * (tm1 disp / pressureTm1) is loaded only once. */
__global__ void kernelPMLElasticStep(int32_t lenum,
				     gpu_data_t* gpuData,
				     int32_t pmlElasticCount);

__global__ void kernelPMLAcousticStep(int32_t lenum,
				      gpu_data_t* gpuData,
				      int32_t pmlAcousticCount);

/* Physics functions used on both host and device */
#ifdef __cplusplus
extern "C" {
#endif
__host__ __device__ int vector_is_zero( const fvector_t* v );
__host__ __device__ int vector_is_all_zero( const fvector_t* v );
__host__ __device__ void aTransposeU( fvector_t* un, double* atu );
__host__ __device__ void firstVector( const double* atu, double* finalVector, double a, double c, double b );
__host__ __device__ void au( fvector_t* resVec, const double* u );
__host__ __device__ void reformF( const double* u, double* newU );
__host__ __device__ void reformU( const double* u, double* newU );
__host__ __device__ void firstVector_kappa( const double* atu, double* finalVector, double kappa);
__host__ __device__ void firstVector_mu( const double* atu, double* finalVector, double b);

__host__ __device__ void MultAddMatVecGPU( fmatrix_t* M, fvector_t* V1, 
					   double c, fvector_t* V2 );
#ifdef __cplusplus
}
#endif

__device__ double atomicAdd(double* address, double val);

#endif /* KERNEL_H_ */
