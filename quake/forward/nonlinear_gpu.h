/* -*- C++ -*- */
#ifndef HERCULES_NONLINEAR_GPU_H
#define HERCULES_NONLINEAR_GPU_H

#include "nonlinear.h"

typedef struct nonlinear_gpu_config_t {
    materialmodel_t material_model;
    noyesflag_t approximate_geostatic_state;
    noyesflag_t tension_cutoff;
    int32_t geostatic_final_step;
    int32_t no_substeps;
    double error_tolerance;
    double backbone_error_tolerance;
    double stiff_damp;
    double stiff_damp_frequency;
} nonlinear_gpu_config_t;

/* Complete state captured after one prescribed material-point strain step.
 * This is intentionally independent of the mesh solver so the CPU and CUDA
 * constitutive implementations can be regression-tested directly. */
typedef struct nonlinear_material_point_result_t {
    tensor_t stress;
    tensor_t plastic_strain;
    tensor_t backstress;
    double ep_bar;
    double fs;
    double psi;
    double load_unload;
    double tao;
    double tao_max;
    double kappa;
    tensor_t stress_reference;
    double gamma1d;
    double tao1d;
    double ggmax;
    double model_error;
    int32_t status;
} nonlinear_material_point_result_t;

/* Backend selection is deliberately runtime-only: it changes execution, not
 * the physical input model. CPU remains the default/reference path. */
int nonlinear_gpu_requested(void);
int nonlinear_gpu_active(void);

int nonlinear_gpu_initialize(int32_t myID, mesh_t *mesh, mysolver_t *solver,
                             int32_t nonlinear_count,
                             const int32_t *nonlinear_mapping,
                             nlsolver_t *host_state,
                             const nonlinear_gpu_config_t *config);

void nonlinear_gpu_enqueue(mysolver_t *solver, int32_t step, double dt,
                           double dt_squared);

/* Copy device state into the existing host structures. This is used only on
 * nonlinear station output steps and at finalization/statistics. */
void nonlinear_gpu_download_state(mysolver_t *solver, nlsolver_t *host_state);

/* Called after the force-stream wait. Reports the first numerical failure
 * recorded by the kernel without adding another time-step synchronization. */
void nonlinear_gpu_check_error(int32_t myID);

int nonlinear_gpu_material_point_history(
    const nlconstants_t *constants,
    const nonlinear_gpu_config_t *config,
    const tensor_t *strain_history,
    int32_t steps,
    nonlinear_material_point_result_t *results);

void nonlinear_gpu_finalize(void);

#endif /* HERCULES_NONLINEAR_GPU_H */
