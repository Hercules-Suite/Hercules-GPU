/* Dedicated solver instantiation for the private intermediate-bias RIVA-Sand
 * research successor. Keeping this separate from the production
 * RIVA-Sand translation prevents either model from inheriting the other's
 * compile-time constitutive path. */
#define NL_GPU_PRODUCTION_ONLY
#define NL_GPU_RIVA_INTERMEDIATE_BIAS_RESEARCH_ONLY
#include "nonlinear_gpu_implementation.cuh"
