/* Dedicated production CUDA translation unit for the frozen PJ-Liq V8
 * constitutive model. Keeping it separate prevents its large explicit state
 * update from increasing register pressure in the other nonlinear kernels. */
#define NL_GPU_PRODUCTION_ONLY
#define NL_GPU_PJ_LIQ_V8_ONLY
#include "nonlinear_gpu_implementation.cuh"
