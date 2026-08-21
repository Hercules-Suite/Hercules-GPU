/* Dedicated production CUDA translation unit for the frozen RIVA-Sand
 * constitutive model. Keeping it separate prevents its large explicit state
 * update from increasing register pressure in the other nonlinear kernels. */
#define NL_GPU_PRODUCTION_ONLY
#define NL_GPU_RIVA_SAND_ONLY
#include "nonlinear_gpu_implementation.cuh"
