/* Standalone host/CUDA material-point replay for the private intermediate-bias
 * RIVA-Sand research successor.  Keep this translation unit independent of the
 * generic material-point kernel: the research state and constitutive core are
 * substantially larger and should not increase compilation cost or register
 * pressure for the established nonlinear fixtures. */

#include <cuda_runtime.h>

#include <math.h>
#include <stdint.h>
#include <stdio.h>

/* nonlinear_gpu.h exposes solver-facing declarations whose mesh, matrix, and
 * MPI types are defined by psolve.h.  The shared nonlinear CUDA translation
 * unit includes this prerequisite before nonlinear_gpu.h; preserve the same
 * order in this isolated material-point unit. */
#include "psolve.h"
#include "nonlinear_gpu.h"
#include "nonlinear_riva_sand_intermediate_bias_research.h"

static_assert(NL_RIVA_INTERMEDIATE_BIAS_STATE_VALUE_COUNT ==
                  RIVA_IB_STATE_VALUE_COUNT,
              "intermediate-bias material-point state width mismatch");

__host__ __device__ static inline int riva_ib_zero_increment(tensor_t value)
{
    return value.xx == 0.0 && value.yy == 0.0 && value.zz == 0.0 &&
           value.xy == 0.0 && value.yz == 0.0 && value.xz == 0.0;
}

__host__ __device__ static inline void riva_ib_material_point_result(
    const riva_ib_state_t *state,
    const riva_update_info_t *info,
    int ok,
    int zero_increment,
    nonlinear_material_point_result_t *result)
{
    *result = {};
    result->stress = state->base.stress;
    result->backstress = state->base.alpha;
    result->ep_bar = state->base.lambda_total;
    result->fs = zero_increment
        ? 0.0 : fabs(riva_compatibility_residual(&state->base));
    result->psi = state->base.D;
    result->load_unload = (double)info->accepted_substeps;
    result->kappa = riva_pressure(state->base.stress);
    result->stress_reference = state->base.alpha01;
    result->gamma1d = state->base.eps_v_irreversible;
    result->tao1d = state->base.eps_v_reversible;
    result->ggmax = state->base.amplitude_factor;
    result->model_error = result->fs;
    result->fabric = state->base.fabric;
    result->alpha_in = state->base.alpha0;
    result->void_ratio = state->base.void_ratio;
    result->plastic_multiplier = state->base.lambda_total;
    result->accepted_substeps = info->accepted_substeps;
    riva_ib_state_values(state,
        result->riva_intermediate_bias_state_values);
    result->riva_geostatic_admitted = state->base.geostatic_admitted;
    result->riva_reversal_registered = info->reversal_registered;
    result->status = ok ? 0 : 1;
}

__host__ __device__ static inline int riva_ib_material_point_replay(
    const riva_material_parameters_t *material,
    double void_ratio,
    const nonlinear_gpu_config_t *config,
    const nonlinear_material_point_initial_state_t *initial_state,
    const tensor_t *strain_history,
    int32_t steps,
    nonlinear_material_point_result_t *results)
{
    if (material == NULL || config == NULL || initial_state == NULL ||
        strain_history == NULL || results == NULL || steps <= 0 ||
        config->material_model != RIVA_SAND_INTERMEDIATE_BIAS_RESEARCH ||
        !initial_state->provided || config->no_substeps < 1) return 0;

    const riva_ib_parameters_t parameters =
        riva_ib_reference_parameters(1000.0);
    tensor_t strain_previous = initial_state->strain;
    tensor_t reference_stress = initial_state->stress;
    reference_stress.xy = 0.0;
    reference_stress.yz = 0.0;
    reference_stress.xz = 0.0;
    riva_ib_state_t state = {};
    int32_t admitted = 0;
    int ok = riva_ib_initialize_material(
        &parameters, material, initial_state->stress, void_ratio, &state);
    if (ok)
        ok = riva_ib_admit_geostatic_state(
            &parameters, material, &state, &admitted);
    if (ok)
        ok = riva_ib_begin_dynamic_phase(
            &parameters, material, &reference_stress, &state);
    if (!ok) return 0;

    for (int32_t step = 0; step < steps; step++) {
        riva_update_info_t info = {};
        const tensor_t increment = riva_sub(
            strain_history[step], strain_previous);
        tensor_t stress = state.base.stress;
        const int zero_increment = riva_ib_zero_increment(increment);
        if (!zero_increment) {
            ok = riva_ib_update_material(
                &parameters, material, increment, config->no_substeps,
                &state, &stress, &info);
        } else {
            ok = 1;
        }
        riva_ib_material_point_result(
            &state, &info, ok, zero_increment, &results[step]);
        strain_previous = strain_history[step];
        if (!ok) return 1;
    }
    return 1;
}

__global__ static void kernel_riva_ib_material_point_history(
    riva_material_parameters_t material,
    double void_ratio,
    nonlinear_gpu_config_t config,
    nonlinear_material_point_initial_state_t initial_state,
    const tensor_t *strain_history,
    int32_t steps,
    nonlinear_material_point_result_t *results)
{
    (void)riva_ib_material_point_replay(
        &material, void_ratio, &config, &initial_state, strain_history, steps,
        results);
}

int nonlinear_riva_intermediate_bias_research_host_history(
    const riva_material_parameters_t *riva_material,
    double riva_void_ratio,
    const nonlinear_gpu_config_t *config,
    const nonlinear_material_point_initial_state_t *initial_state,
    const tensor_t *strain_history,
    int32_t steps,
    nonlinear_material_point_result_t *results)
{
    return riva_ib_material_point_replay(
        riva_material, riva_void_ratio, config, initial_state, strain_history,
        steps, results) ? 0 : -1;
}

int nonlinear_gpu_riva_intermediate_bias_research_material_point_history(
    const riva_material_parameters_t *riva_material,
    double riva_void_ratio,
    const nonlinear_gpu_config_t *config,
    const nonlinear_material_point_initial_state_t *initial_state,
    const tensor_t *strain_history,
    int32_t steps,
    nonlinear_material_point_result_t *results)
{
    if (riva_material == NULL || config == NULL || initial_state == NULL ||
        strain_history == NULL || results == NULL || steps <= 0) return -1;

    tensor_t *strain_device = NULL;
    nonlinear_material_point_result_t *result_device = NULL;
    cudaError_t error = cudaMalloc(
        (void **)&strain_device, (size_t)steps*sizeof(*strain_device));
    if (error == cudaSuccess)
        error = cudaMalloc(
            (void **)&result_device, (size_t)steps*sizeof(*result_device));
    if (error == cudaSuccess)
        error = cudaMemcpy(
            strain_device, strain_history,
            (size_t)steps*sizeof(*strain_device), cudaMemcpyHostToDevice);
    if (error == cudaSuccess) {
        kernel_riva_ib_material_point_history<<<1, 1>>>(
            *riva_material, riva_void_ratio, *config, *initial_state,
            strain_device, steps, result_device);
        error = cudaGetLastError();
    }
    if (error == cudaSuccess) error = cudaDeviceSynchronize();
    if (error == cudaSuccess)
        error = cudaMemcpy(
            results, result_device, (size_t)steps*sizeof(*result_device),
            cudaMemcpyDeviceToHost);

    if (result_device != NULL) cudaFree(result_device);
    if (strain_device != NULL) cudaFree(strain_device);
    if (error != cudaSuccess) {
        fprintf(stderr,
                "RIVA-Sand intermediate-bias research GPU material-point "
                "test failed: %s\n",
                cudaGetErrorString(error));
        return -1;
    }
    return 0;
}
