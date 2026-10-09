/* -*- C++ -*- */
#ifndef HERCULES_NONLINEAR_GPU_H
#define HERCULES_NONLINEAR_GPU_H

#include <sys/types.h>

#include "nonlinear.h"
#include "quake_util.h"

#define NL_RIVA_STATE_VALUE_COUNT 93
#define NL_RIVA_INTERMEDIATE_BIAS_STATE_VALUE_COUNT 138

typedef enum nl_md_recovery_mode_t {
    NL_MD_RECOVERY_STRICT = 0,
    NL_MD_RECOVERY_FORCE_ACCEPT = 1,
    NL_MD_RECOVERY_UNLIMITED_SUBSTEPS = 2,
    NL_MD_RECOVERY_OPENSEES = 3
} nl_md_recovery_mode_t;

typedef enum nl_md_integration_scheme_t {
    NL_MD_INTEGRATION_MODIFIED_EULER = 0,
    NL_MD_INTEGRATION_FORWARD_EULER = 1
} nl_md_integration_scheme_t;

typedef enum nl_md_forward_euler_substeps_t {
    NL_MD_FORWARD_EULER_SUBSTEPS_FIXED = 0,
    NL_MD_FORWARD_EULER_SUBSTEPS_DEMAND = 1
} nl_md_forward_euler_substeps_t;

typedef struct nonlinear_gpu_config_t {
    materialmodel_t material_model;
    noyesflag_t approximate_geostatic_state;
    noyesflag_t tension_cutoff;
    int32_t geostatic_final_step;
    int32_t no_substeps;
    int32_t riva_reversal_type; /* RIVASand02/Custom02: 1 (default), 2 or 3 */
    nl_md_integration_scheme_t md_integration_scheme;
    nl_md_forward_euler_substeps_t md_forward_euler_substeps;
    nl_md_recovery_mode_t md_recovery_mode;
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
    tensor_t fabric;
    tensor_t alpha_in;
    double void_ratio;
    double plastic_multiplier;
    double integration_error;
    int32_t accepted_substeps;
    int32_t pdmy03_active_surface;
    int32_t pdmy03_requested_substeps;
    int32_t pdmy03_crossed_surfaces;
    int32_t pdmy03_substeps_clamped;
    double pdmy03_pressure_d;
    double pdmy03_center_norm;
    double riva_state_values[NL_RIVA_STATE_VALUE_COUNT];
    double riva_intermediate_bias_state_values[
        NL_RIVA_INTERMEDIATE_BIAS_STATE_VALUE_COUNT];
    int32_t riva_geostatic_admitted;
    int32_t riva_reversal_registered;
    int32_t status;
} nonlinear_material_point_result_t;

/* Optional committed state from which a prescribed material-point history is
 * replayed. This is used by independent reference comparisons that perform a
 * separate elastic consolidation stage before enabling plasticity. */
typedef struct nonlinear_material_point_initial_state_t {
    tensor_t strain;
    tensor_t stress;
    tensor_t plastic_strain;
    tensor_t alpha;
    tensor_t alpha_in;
    tensor_t fabric;
    int32_t provided;
} nonlinear_material_point_initial_state_t;

/* Backend selection is deliberately runtime-only: it changes execution, not
 * the physical input model. GPU is the default; CPU remains available as an
 * explicit reference/diagnostic override. */
int nonlinear_gpu_requested(void);
int nonlinear_gpu_active(void);

int nonlinear_gpu_initialize(int32_t myID, mesh_t *mesh, mysolver_t *solver,
                             int32_t nonlinear_count,
                             const int32_t *nonlinear_mapping,
                             nlsolver_t *host_state,
                             const nonlinear_gpu_config_t *config);

void nonlinear_gpu_enqueue(mysolver_t *solver, int32_t step, double dt,
                           double dt_squared);

/* Configure and download compact nonlinear-station records. The supplied
 * indices address the nonlinear-element arrays, not the full mesh arrays. */
int nonlinear_gpu_configure_station_output(
    int32_t myID, int32_t station_count,
    const int32_t *station_nonlinear_indices);
void nonlinear_gpu_download_station_state(mysolver_t *solver,
                                          nlsolver_t *host_state);

/* Copy the complete device state into the existing host structures. This is
 * reserved for finalization/statistics; station output uses the sparse gather
 * above. */
void nonlinear_gpu_download_state(mysolver_t *solver, nlsolver_t *host_state);

/* Called after the force-stream wait. Reports the first numerical failure
 * recorded by the kernel without adding another time-step synchronization. */
void nonlinear_gpu_check_error(int32_t myID);

int nonlinear_gpu_material_point_history(
    const nlconstants_t *constants,
    const manzari_dafalias_constants_t *md_material,
    const pdmy03_constants_t *pdmy03_material,
    const riva_material_parameters_t *riva_material,
    double riva_void_ratio,
    const nonlinear_gpu_config_t *config,
    const nonlinear_material_point_initial_state_t *initial_state,
    const tensor_t *strain_history,
    int32_t steps,
    nonlinear_material_point_result_t *results);

/* Host execution of the same allocation-free Manzari--Dafalias integrator is
 * used only by the material-point regression path. Production MD updates are
 * CUDA-only. */
int nonlinear_manzari_dafalias_host_history(
    const manzari_dafalias_constants_t *md_material,
    const nonlinear_gpu_config_t *config,
    const nonlinear_material_point_initial_state_t *initial_state,
    const tensor_t *strain_history,
    int32_t steps,
    nonlinear_material_point_result_t *results);

int nonlinear_pdmy03_host_history(
    const pdmy03_constants_t *pdmy03_material,
    const nonlinear_gpu_config_t *config,
    const nonlinear_material_point_initial_state_t *initial_state,
    const tensor_t *strain_history,
    int32_t steps,
    nonlinear_material_point_result_t *results);

int nonlinear_riva_host_history(
    const nlconstants_t *constants,
    const riva_material_parameters_t *riva_material,
    double riva_void_ratio,
    const nonlinear_gpu_config_t *config,
    const nonlinear_material_point_initial_state_t *initial_state,
    const tensor_t *strain_history,
    int32_t steps,
    nonlinear_material_point_result_t *results);

/* The intermediate-bias research material has a distinct 138-value restart
 * contract.  Its material-point replay is compiled separately so adding the
 * larger research kernel does not increase the existing all-model regression
 * translation unit. */
int nonlinear_riva_intermediate_bias_research_host_history(
    const riva_material_parameters_t *riva_material,
    double riva_void_ratio,
    const nonlinear_gpu_config_t *config,
    const nonlinear_material_point_initial_state_t *initial_state,
    const tensor_t *strain_history,
    int32_t steps,
    nonlinear_material_point_result_t *results);

int nonlinear_gpu_riva_intermediate_bias_research_material_point_history(
    const riva_material_parameters_t *riva_material,
    double riva_void_ratio,
    const nonlinear_gpu_config_t *config,
    const nonlinear_material_point_initial_state_t *initial_state,
    const tensor_t *strain_history,
    int32_t steps,
    nonlinear_material_point_result_t *results);

/* Append/read the versioned nonlinear/BrickUP extension of a displacement
 * checkpoint. The legacy displacement prefix is unchanged; extension_offset
 * is the first byte after its reserved per-rank records. */
int nonlinear_gpu_checkpoint_write(
    const char *filename, off_t extension_offset, int32_t myID,
    int32_t group_size, int32_t node_count, mysolver_t *solver,
    MPI_Comm comm_solver);
int nonlinear_gpu_checkpoint_read(
    const char *filename, off_t extension_offset, int32_t myID,
    int32_t group_size, int32_t node_count, mysolver_t *solver,
    MPI_Comm comm_solver);

void nonlinear_gpu_finalize(void);

#endif /* HERCULES_NONLINEAR_GPU_H */
