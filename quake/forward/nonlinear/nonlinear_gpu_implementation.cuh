/* -*- C++ -*- */
/*
 * Shared implementation included by the independent nonlinear CUDA build
 * units.  Each production model family instantiates its own kernel so ptxas
 * can optimize the constitutive paths independently and Make can compile them
 * concurrently.
 *
 * The legacy nonlinear.cc implementation remains the CPU/GSL reference.
 * These CUDA build units intentionally contain no GSL calls and perform no
 * dynamic allocation from device code.
 */

#include <cuda.h>
#include <cuda_runtime.h>

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "util.h"
#include "psolve.h"
#include "kernel.h"
#include "nonlinear_gpu.h"
#include "nonlinear_pdmy03.h"
#if defined(NL_GPU_RUNTIME_ONLY) || defined(NL_GPU_PJ_LIQ_V8_ONLY) || \
    defined(NL_GPU_MATERIAL_POINT_ONLY)
#include "nonlinear_pj_liq_v8.h"
static_assert(NL_PJ8_STATE_VALUE_COUNT == PJ8_STATE_VALUE_COUNT,
              "PJ-Liq V8 material-point state width mismatch");
#endif

#define NL_GPU_QP 8
#define NL_GPU_UPDATE_BLOCK_SIZE 128
#define NL_GPU_KAPPA_MIN 1.0e-10
#define NL_GPU_PI 3.141592653589793238462643383279502884

#if (defined(NL_GPU_RUNTIME_ONLY) + defined(NL_GPU_CLASSIC_ONLY) + \
     defined(NL_GPU_MANZARI_DAFALIAS_ONLY) + \
     defined(NL_GPU_PDMY03_ONLY) + \
     defined(NL_GPU_PJ_LIQ_V8_ONLY) + \
     defined(NL_GPU_MATERIAL_POINT_ONLY)) != 1
#error "select exactly one nonlinear CUDA compilation unit"
#endif

#ifndef NL_GPU_MATERIAL_POINT_ONLY
typedef struct nonlinear_gpu_element_t {
    int32_t eindex;
    double h;
} nonlinear_gpu_element_t;

typedef struct nonlinear_gpu_error_t {
    int32_t code;
    int32_t eindex;
    int32_t qp;
    int32_t md_reason;
    int32_t md_attempts;
    int32_t md_accepted_substeps;
    int32_t md_correction_retries;
    int32_t md_correction_iterations;
    int32_t md_projection_fallbacks;
    int32_t md_force_accepts;
    double residual;
    double value;
    double md_substep_fraction;
    double md_progress_fraction;
    double md_effective_mean_stress;
    double md_strain_increment_norm;
    double md_first_dgamma;
    double md_second_dgamma;
    double md_first_raw_dgamma;
    double md_second_raw_dgamma;
    double md_first_denominator;
    double md_second_denominator;
    double md_first_h;
    double md_second_h;
    double md_first_hardening_gap;
    double md_second_hardening_gap;
    tensor_t md_stress;
    tensor_t md_strain_increment;
} nonlinear_gpu_error_t;

/* Only these quadrature-point values are consumed by
 * print_nonlinear_stations(). Keeping the record compact prevents station
 * output from downloading every per-element history array. */
typedef struct nonlinear_gpu_station_state_t {
    int32_t nonlinear_index;
    tensor_t stress;
    tensor_t strain;
    double fs;
    double kappa;
    double gamma1d;
    double tao1d;
    double ggmax;
    double sigma_z_static;
    double phi;
} nonlinear_gpu_station_state_t;

typedef struct nonlinear_gpu_station_buffers_t {
    int32_t count;
    int32_t *indices;
    nonlinear_gpu_station_state_t *device;
    nonlinear_gpu_station_state_t *host;
} nonlinear_gpu_station_buffers_t;

typedef struct nonlinear_gpu_workload_stats_t {
    unsigned long long zero_or_inactive_points;
    unsigned long long one_step_points;
    unsigned long long multi_step_points;
    unsigned long long accepted_substeps;
    unsigned long long maximum_substeps;
} nonlinear_gpu_workload_stats_t;

typedef struct nonlinear_gpu_workload_buffers_t {
    int32_t interval;
    int32_t pending;
    int32_t step;
    nonlinear_gpu_workload_stats_t *device;
    nonlinear_gpu_workload_stats_t host;
} nonlinear_gpu_workload_buffers_t;

typedef struct nonlinear_gpu_state_t {
    int32_t count;
    int32_t pdmy03_surface_stride;
    nonlinear_gpu_element_t *elements;
    nlconstants_t *constants;
    qptensors_t *stresses;
    qptensors_t *strains;
    qptensors_t *strains1;
    qptensors_t *pstrains1;
    qptensors_t *pstrains2;
    qptensors_t *alphastress1;
    qptensors_t *alphastress2;
    qpvectors_t *ep1;
    qpvectors_t *ep2;
    qpvectors_t *psi_n;
    qpvectors_t *load_unload_n;
    qpvectors_t *sv_n;
    qpvectors_t *sv_max;
    qpvectors_t *kappa;
    qptensors_t *sref;
    qptensors_t *fabric;
    qpvectors_t *gamma1d;
    qpvectors_t *tao1d;
    qpvectors_t *ggmax;
    pdmy03_state_t *pdmy03_states;
    tensor_t *pdmy03_centers;
    void *pj8_parameters;
    void *pj8_states;
    nonlinear_gpu_error_t *error;
} nonlinear_gpu_state_t;

extern "C" void nonlinear_gpu_launch_classic(
    nonlinear_gpu_state_t state, nonlinear_gpu_config_t config,
    gpu_data_t *gpu, int32_t grid_size, int32_t block_size,
    int32_t step, double delta_t, double delta_t_squared,
    cudaStream_t stream);
extern "C" void nonlinear_gpu_launch_manzari_dafalias(
    nonlinear_gpu_state_t state, nonlinear_gpu_config_t config,
    gpu_data_t *gpu, int32_t grid_size, int32_t block_size,
    int32_t step, double delta_t, double delta_t_squared,
    cudaStream_t stream);
extern "C" void nonlinear_gpu_launch_pdmy03(
    nonlinear_gpu_state_t state, nonlinear_gpu_config_t config,
    gpu_data_t *gpu, int32_t grid_size, int32_t block_size,
    int32_t step, double delta_t, double delta_t_squared,
    cudaStream_t stream);
extern "C" void nonlinear_gpu_launch_pj_liq_v8(
    nonlinear_gpu_state_t state, nonlinear_gpu_config_t config,
    gpu_data_t *gpu, int32_t grid_size, int32_t block_size,
    int32_t step, double delta_t, double delta_t_squared,
    cudaStream_t stream);

#ifdef NL_GPU_RUNTIME_ONLY
static nonlinear_gpu_state_t g_nl = {};
static nonlinear_gpu_station_buffers_t g_station = {};
static nonlinear_gpu_workload_buffers_t g_workload = {};
static nonlinear_gpu_config_t g_config = {};
static nonlinear_gpu_error_t g_error_host = {};
static mesh_t *g_mesh_host = NULL;
static int g_backend_request = -1;
static int g_active = 0;

static int parse_boolean_gpu(const char *value)
{
    return (strcasecmp(value, "gpu") == 0 ||
            strcasecmp(value, "device") == 0 ||
            strcasecmp(value, "yes") == 0 ||
            strcasecmp(value, "true") == 0 || strcmp(value, "1") == 0);
}

static int parse_boolean_cpu(const char *value)
{
    return (strcasecmp(value, "cpu") == 0 ||
            strcasecmp(value, "host") == 0 ||
            strcasecmp(value, "no") == 0 ||
            strcasecmp(value, "false") == 0 || strcmp(value, "0") == 0);
}

static int nonlinear_gpu_workload_interval(int32_t myID)
{
    const char *value = getenv("HERCULES_NONLINEAR_STATS_INTERVAL");
    if (value == NULL || value[0] == '\0' || strcmp(value, "0") == 0)
        return 0;

    char *end = NULL;
    const long interval = strtol(value, &end, 10);
    if (end == value || *end != '\0' || interval < 0 ||
        interval > INT32_MAX) {
        fprintf(stderr,
                "Thread %d: HERCULES_NONLINEAR_STATS_INTERVAL must be "
                "zero or a positive integer; got '%s'\n",
                myID, value);
        return -1;
    }
    return (int)interval;
}

int nonlinear_gpu_requested(void)
{
    if (g_backend_request < 0) {
        const char *value = getenv("HERCULES_NONLINEAR_BACKEND");
        /* This CUDA build defaults nonlinear constitutive updates and force
         * assembly to the device path. Keep an explicit CPU override for
         * reference comparisons and diagnosis. */
        g_backend_request = 1;
        if (value != NULL) {
            if (parse_boolean_gpu(value)) {
                g_backend_request = 1;
            } else if (parse_boolean_cpu(value)) {
                g_backend_request = 0;
            } else {
                fprintf(stderr,
                        "HERCULES_NONLINEAR_BACKEND must be cpu or gpu; got '%s'\n",
                        value);
                exit(EXIT_FAILURE);
            }
        }
    }
    return g_backend_request;
}

int nonlinear_gpu_active(void)
{
    return g_active;
}
#endif
#endif

__device__ static inline double nl_max(double a, double b)
{
    return a > b ? a : b;
}

__device__ static inline double nl_min(double a, double b)
{
    return a < b ? a : b;
}

__device__ static inline tensor_t nl_zero(void)
{
    tensor_t a = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    return a;
}

__device__ static inline tensor_t nl_add(tensor_t a, tensor_t b)
{
    tensor_t c = {a.xx + b.xx, a.yy + b.yy, a.zz + b.zz,
                  a.xy + b.xy, a.yz + b.yz, a.xz + b.xz};
    return c;
}

__device__ static inline tensor_t nl_sub(tensor_t a, tensor_t b)
{
    tensor_t c = {a.xx - b.xx, a.yy - b.yy, a.zz - b.zz,
                  a.xy - b.xy, a.yz - b.yz, a.xz - b.xz};
    return c;
}

__device__ static inline tensor_t nl_scale(tensor_t a, double s)
{
    tensor_t c = {a.xx * s, a.yy * s, a.zz * s,
                  a.xy * s, a.yz * s, a.xz * s};
    return c;
}

__device__ static inline tensor_t nl_iso(double s)
{
    tensor_t a = {s, s, s, 0.0, 0.0, 0.0};
    return a;
}

__device__ static inline double nl_i1(tensor_t a)
{
    return a.xx + a.yy + a.zz;
}

__device__ static inline tensor_t nl_dev(tensor_t a)
{
    const double p = nl_i1(a) / 3.0;
    a.xx -= p;
    a.yy -= p;
    a.zz -= p;
    return a;
}

__device__ static inline double nl_j2(tensor_t a)
{
    return 0.5 * (a.xx * a.xx + a.yy * a.yy + a.zz * a.zz) +
           a.xy * a.xy + a.yz * a.yz + a.xz * a.xz;
}

__device__ static inline double nl_ddot(tensor_t a, tensor_t b)
{
    return a.xx * b.xx + a.yy * b.yy + a.zz * b.zz +
           2.0 * (a.xy * b.xy + a.yz * b.yz + a.xz * b.xz);
}

__device__ static inline double nl_comb_j2(tensor_t a, tensor_t b)
{
    return 0.5 * (a.xx * b.xx + a.yy * b.yy + a.zz * b.zz) +
           a.xy * b.xy + a.yz * b.yz + a.xz * b.xz;
}

__device__ static inline tensor_t nl_stress(tensor_t strain,
                                             double mu, double lambda)
{
    const double lkk = lambda * nl_i1(strain);
    const double mu2 = 2.0 * mu;
    tensor_t stress = {mu2 * strain.xx + lkk,
                       mu2 * strain.yy + lkk,
                       mu2 * strain.zz + lkk,
                       mu2 * strain.xy,
                       mu2 * strain.yz,
                       mu2 * strain.xz};
    return stress;
}

__device__ static inline tensor_t nl_elastic_strain(tensor_t stress,
                                                     double mu, double bulk)
{
    const double skk = nl_i1(stress);
    tensor_t dev = nl_dev(stress);
    tensor_t strain = nl_scale(dev, 1.0 / (2.0 * mu));
    const double ev = skk / (9.0 * bulk);
    strain.xx += ev;
    strain.yy += ev;
    strain.zz += ev;
    return strain;
}

/* -------------------------------------------------------------------------- */
/* Manzari--Dafalias explicit Modified-Euler integration.                      */
/*                                                                            */
/* The core is host/device so the material-point command can run the exact     */
/* same arithmetic on the host and GPU. Production execution remains GPU-only. */
/* Hercules stores tension-positive stress/strain; the MD equations below use  */
/* compression-positive tensors, converted only at the integration boundary.   */
/* -------------------------------------------------------------------------- */

#define NL_HD __host__ __device__

typedef enum nl_md_failure_reason_t {
    NL_MD_FAILURE_NONE = 0,
    NL_MD_FAILURE_INITIAL_PRESSURE = 1,
    NL_MD_FAILURE_INVALID_CONTROLS = 2,
    NL_MD_FAILURE_ELASTIC_MODULI = 3,
    NL_MD_FAILURE_YIELD_INTERSECTION = 4,
    NL_MD_FAILURE_ATTEMPT_LIMIT = 5,
    NL_MD_FAILURE_FIRST_STAGE = 6,
    NL_MD_FAILURE_TENSILE_PREDICTOR = 7,
    NL_MD_FAILURE_SECOND_STAGE = 8,
    NL_MD_FAILURE_NONFINITE_LOCAL_ERROR = 9,
    NL_MD_FAILURE_LOCAL_ERROR_MIN_STEP = 10,
    NL_MD_FAILURE_SUBSTEP_LIMIT = 11,
    NL_MD_FAILURE_STRESS_CORRECTION_MIN_STEP = 12,
    NL_MD_FAILURE_FINAL_STATE = 13
} nl_md_failure_reason_t;

typedef struct nl_md_diagnostics_t {
    int32_t failure_reason;
    int32_t attempts;
    int32_t correction_retries;
    int32_t correction_iterations;
    int32_t projection_fallbacks;
    int32_t force_accepts;
    double substep_fraction;
    double progress_fraction;
    double effective_mean_stress;
    double local_error;
    double correction_residual;
    double first_dgamma;
    double second_dgamma;
    double first_raw_dgamma;
    double second_raw_dgamma;
    double first_denominator;
    double second_denominator;
    double first_h;
    double second_h;
    double first_hardening_gap;
    double second_hardening_gap;
} nl_md_diagnostics_t;

typedef struct nl_md_stage_t {
    tensor_t dsigma;
    tensor_t dalpha;
    tensor_t dfabric;
    tensor_t dpstrain;
    double dgamma;
    double raw_dgamma;
    double loading_numerator;
    double denominator;
    double h;
    double hardening_gap;
    int neutral;
    int ok;
} nl_md_stage_t;

typedef struct nl_md_dependent_t {
    tensor_t n;
    tensor_t b;
    tensor_t R;
    double h;
    double D;
    double B;
    double C;
    double hardening_gap;
    int ok;
} nl_md_dependent_t;

NL_HD static inline int md_force_accept_enabled(
    nl_md_recovery_mode_t mode)
{
    return mode == NL_MD_RECOVERY_FORCE_ACCEPT ||
           mode == NL_MD_RECOVERY_OPENSEES;
}

NL_HD static inline int md_unlimited_substeps_enabled(
    nl_md_recovery_mode_t mode)
{
    return mode == NL_MD_RECOVERY_UNLIMITED_SUBSTEPS ||
           mode == NL_MD_RECOVERY_OPENSEES;
}

NL_HD static inline int md_valid_recovery_mode(nl_md_recovery_mode_t mode)
{
    return mode >= NL_MD_RECOVERY_STRICT &&
           mode <= NL_MD_RECOVERY_OPENSEES;
}

NL_HD static inline int md_valid_integration_scheme(
    nl_md_integration_scheme_t scheme)
{
    return scheme == NL_MD_INTEGRATION_MODIFIED_EULER ||
           scheme == NL_MD_INTEGRATION_FORWARD_EULER;
}

NL_HD static inline int md_valid_forward_euler_substeps(
    nl_md_forward_euler_substeps_t mode)
{
    return mode == NL_MD_FORWARD_EULER_SUBSTEPS_FIXED ||
           mode == NL_MD_FORWARD_EULER_SUBSTEPS_DEMAND;
}

NL_HD static inline tensor_t md_zero(void)
{
    tensor_t a = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    return a;
}

NL_HD static inline tensor_t md_add(tensor_t a, tensor_t b)
{
    tensor_t c = {a.xx+b.xx, a.yy+b.yy, a.zz+b.zz,
                  a.xy+b.xy, a.yz+b.yz, a.xz+b.xz};
    return c;
}

NL_HD static inline tensor_t md_sub(tensor_t a, tensor_t b)
{
    tensor_t c = {a.xx-b.xx, a.yy-b.yy, a.zz-b.zz,
                  a.xy-b.xy, a.yz-b.yz, a.xz-b.xz};
    return c;
}

NL_HD static inline tensor_t md_scale(tensor_t a, double s)
{
    tensor_t c = {a.xx*s, a.yy*s, a.zz*s, a.xy*s, a.yz*s, a.xz*s};
    return c;
}

NL_HD static inline tensor_t md_iso(double s)
{
    tensor_t a = {s, s, s, 0.0, 0.0, 0.0};
    return a;
}

NL_HD static inline double md_trace(tensor_t a)
{
    return a.xx+a.yy+a.zz;
}

NL_HD static inline tensor_t md_dev(tensor_t a)
{
    const double mean = md_trace(a)/3.0;
    a.xx -= mean; a.yy -= mean; a.zz -= mean;
    return a;
}

NL_HD static inline double md_ddot(tensor_t a, tensor_t b)
{
    return a.xx*b.xx+a.yy*b.yy+a.zz*b.zz+
           2.0*(a.xy*b.xy+a.yz*b.yz+a.xz*b.xz);
}

NL_HD static inline double md_norm(tensor_t a)
{
    return sqrt(fmax(md_ddot(a, a), 0.0));
}

NL_HD static inline int md_finite(tensor_t a)
{
    return isfinite(a.xx) && isfinite(a.yy) && isfinite(a.zz) &&
           isfinite(a.xy) && isfinite(a.yz) && isfinite(a.xz);
}

NL_HD static inline tensor_t md_square(tensor_t a)
{
    tensor_t b;
    b.xx = a.xx*a.xx+a.xy*a.xy+a.xz*a.xz;
    b.yy = a.xy*a.xy+a.yy*a.yy+a.yz*a.yz;
    b.zz = a.xz*a.xz+a.yz*a.yz+a.zz*a.zz;
    b.xy = a.xx*a.xy+a.xy*a.yy+a.xz*a.yz;
    b.yz = a.xy*a.xz+a.yy*a.yz+a.yz*a.zz;
    b.xz = a.xx*a.xz+a.xy*a.yz+a.xz*a.zz;
    return b;
}

NL_HD static inline tensor_t md_elastic_increment(tensor_t de,
                                                   double K, double G)
{
    return md_add(md_scale(md_dev(de), 2.0*G),
                  md_iso(K*md_trace(de)));
}

NL_HD static inline tensor_t md_elastic_compliance(tensor_t dsigma,
                                                    double K, double G)
{
    return md_add(md_scale(md_dev(dsigma), 1.0/(2.0*G)),
                  md_iso(md_trace(dsigma)/(9.0*K)));
}

NL_HD static inline double md_minimum_mean_stress(
    const manzari_dafalias_constants_t *m)
{
    return 1.0e-4*m->P_atm;
}

/* Optional demand-driven Forward-Euler subdivision targets. This is a
 * Manzari--Dafalias-specific heuristic: it limits stress-ratio excursion and
 * octahedral engineering deviatoric strain per subincrement. It is not the
 * PDMY03 surface-spacing criterion. Fixed subdivision remains the default for
 * backward-compatible and reproducible integration-sensitivity studies. */
#define NL_MD_FE_RATIO_STEP  0.05
#define NL_MD_FE_STRAIN_STEP 1.0e-5

/* Convert a real-valued subdivision demand into a substep count. Clamping
 * before the cast is deliberate: converting a non-finite or out-of-range
 * double to int32_t is undefined behavior, and on the device it saturates to
 * INT32_MAX, which is how a diverged point reports a nonsense substep
 * request instead of a usable diagnostic. */
NL_HD static inline int32_t md_substep_request(double demand, int32_t cap)
{
    if (cap < 1) return 1;
    if (!(demand > 0.0)) return 1;                /* also catches NaN */
    if (!(demand < (double)cap)) return cap;      /* also catches +inf */
    return (int32_t)demand + 1;
}

/* Complete the unresolved remainder on the low-effective-pressure cutoff.
 * Sand cannot carry tensile effective mean stress. Collapse the skeleton to
 * the isotropic OpenSees p_min state and assign the remaining total strain,
 * less the elastic stress change, to plastic strain. This preserves the
 * additive strain split instead of returning a partially integrated state. */
NL_HD static int md_complete_at_pressure_floor(
    const manzari_dafalias_constants_t *m, tensor_t current_stress,
    tensor_t remaining_strain, tensor_t *stress, tensor_t *plastic_strain,
    tensor_t *alpha, tensor_t *alpha_in, double K, double G)
{
    const tensor_t floor_stress = md_iso(md_minimum_mean_stress(m));
    const tensor_t elastic_change = md_elastic_compliance(
        md_sub(floor_stress, current_stress), K, G);
    *plastic_strain = md_add(
        *plastic_strain, md_sub(remaining_strain, elastic_change));
    *stress = floor_stress;
    *alpha = md_zero();
    *alpha_in = md_zero();
    return md_finite(*stress) && md_finite(*plastic_strain);
}

/* Replace an already formed candidate by the pressure-floor state while
 * holding its total strain fixed. */
NL_HD static int md_project_candidate_to_pressure_floor(
    const manzari_dafalias_constants_t *m, tensor_t *stress,
    tensor_t *plastic_strain, tensor_t *alpha, tensor_t *alpha_in,
    double K, double G)
{
    const tensor_t old_stress = *stress;
    const tensor_t floor_stress = md_iso(md_minimum_mean_stress(m));
    *plastic_strain = md_add(
        *plastic_strain,
        md_elastic_compliance(md_sub(old_stress, floor_stress), K, G));
    *stress = floor_stress;
    *alpha = md_zero();
    *alpha_in = md_zero();
    return md_finite(*stress) && md_finite(*plastic_strain) &&
           md_finite(*alpha) && md_finite(*alpha_in);
}

NL_HD static inline void md_elastic_moduli(
    const manzari_dafalias_constants_t *m, tensor_t stress,
    double *K, double *G)
{
    const double pmin = 1.0e-4*m->P_atm;
    const double p = fmax(md_trace(stress)/3.0, pmin);
    const double void_factor = (2.97-m->e_init)*(2.97-m->e_init)/
                               (1.0+m->e_init);
    *G = m->G0*m->P_atm*void_factor*sqrt(p/m->P_atm);
    *K = (2.0/3.0)*(1.0+m->nu)/(1.0-2.0*m->nu)*(*G);
}

NL_HD static inline double md_yield(
    const manzari_dafalias_constants_t *m, tensor_t stress, tensor_t alpha)
{
    const double p = md_trace(stress)/3.0+1.0e-2*m->P_atm;
    return md_norm(md_sub(md_dev(stress), md_scale(alpha, p)))-
           sqrt(2.0/3.0)*m->m*p;
}

/* Commit a finite candidate after the adaptive Modified-Euler recovery has
 * reached its configured minimum fraction. This follows the intent of the
 * OpenSees force-accept branch: cap the stress ratio and complete the local
 * increment instead of aborting. Hercules additionally preserves the
 * additive elastic/plastic strain split after the cap and recenters the small
 * kinematic yield surface so the committed state is finite and admissible. */
NL_HD static int md_force_accept_candidate(
    const manzari_dafalias_constants_t *m, tensor_t proposed_stress,
    tensor_t proposed_plastic_strain, tensor_t proposed_fabric,
    tensor_t *stress, tensor_t *plastic_strain, tensor_t *alpha,
    tensor_t *alpha_in, tensor_t *fabric, double K, double G,
    double *residual)
{
    if (!md_finite(proposed_stress) ||
        !md_finite(proposed_plastic_strain) ||
        !md_finite(proposed_fabric)) return 0;

    tensor_t accepted_stress = proposed_stress;
    tensor_t accepted_pstrain = proposed_plastic_strain;
    const double mean_stress = md_trace(accepted_stress)/3.0;
    if (mean_stress < md_minimum_mean_stress(m)) {
        tensor_t accepted_alpha = *alpha;
        tensor_t accepted_alpha_in = *alpha_in;
        if (!md_project_candidate_to_pressure_floor(
                m, &accepted_stress, &accepted_pstrain,
                &accepted_alpha, &accepted_alpha_in, K, G)) return 0;
        *stress = accepted_stress;
        *plastic_strain = accepted_pstrain;
        *alpha = accepted_alpha;
        *alpha_in = accepted_alpha_in;
        *fabric = proposed_fabric;
        *residual = 0.0;
        return md_finite(*fabric);
    }

    const tensor_t proposed_dev = md_dev(accepted_stress);
    const double trace = md_trace(accepted_stress);
    const double eta = sqrt(13.5)*md_norm(proposed_dev)/trace;
    if (isfinite(eta) && eta > m->Mc) {
        const tensor_t capped_stress = md_add(
            md_iso(mean_stress), md_scale(proposed_dev, m->Mc/eta));
        accepted_pstrain = md_add(
            accepted_pstrain,
            md_elastic_compliance(
                md_sub(accepted_stress, capped_stress), K, G));
        accepted_stress = capped_stress;
    }

    const double shifted_p = md_trace(accepted_stress)/3.0+
                             1.0e-2*m->P_atm;
    if (!(shifted_p > 0.0) || !isfinite(shifted_p)) return 0;
    const tensor_t accepted_alpha = md_scale(
        md_dev(accepted_stress), 1.0/shifted_p);
    const double f = md_yield(m, accepted_stress, accepted_alpha);
    if (!isfinite(f) || !md_finite(accepted_pstrain) ||
        !md_finite(accepted_alpha)) return 0;

    *stress = accepted_stress;
    *plastic_strain = accepted_pstrain;
    *alpha = accepted_alpha;
    *alpha_in = accepted_alpha;
    *fabric = proposed_fabric;
    *residual = fmax(f, 0.0);
    return md_finite(*fabric);
}

/* Locate the first yield crossing along an elastic stress increment. The
 * current state must be strictly inside and the full trial state outside.
 * Bisection is inexpensive here and gives identical host/device behavior. */
NL_HD static int md_elastic_yield_intersection(
    const manzari_dafalias_constants_t *m, tensor_t stress,
    tensor_t elastic_increment, tensor_t alpha, double tolerance,
    double *ratio)
{
    double lower = 0.0;
    double upper = 1.0;
    const double flower = md_yield(m, stress, alpha);
    const double fupper = md_yield(
        m, md_add(stress, elastic_increment), alpha);
    if (!(flower < -tolerance) || !(fupper > tolerance) ||
        !isfinite(flower) || !isfinite(fupper)) return 0;

    for (int iteration = 0; iteration < 80; iteration++) {
        const double middle = 0.5*(lower+upper);
        const double fmiddle = md_yield(
            m, md_add(stress, md_scale(elastic_increment, middle)), alpha);
        if (!isfinite(fmiddle)) return 0;
        if (fabs(fmiddle) <= tolerance) {
            *ratio = middle;
            return 1;
        }
        if (fmiddle > 0.0) {
            upper = middle;
        } else {
            lower = middle;
        }
    }
    *ratio = 0.5*(lower+upper);
    return isfinite(*ratio) && *ratio >= 0.0 && *ratio <= 1.0;
}

NL_HD static nl_md_dependent_t md_state_dependent(
    const manzari_dafalias_constants_t *m, tensor_t stress, tensor_t alpha,
    tensor_t fabric, double void_ratio, tensor_t alpha_in)
{
    nl_md_dependent_t out = {};
    const double small = 1.0e-14;
    const double p = md_trace(stress)/3.0+1.0e-2*m->P_atm;
    if (!(p > 0.0) || !(m->c > 0.0) || !(m->P_atm > 0.0)) return out;

    tensor_t shifted = md_sub(md_dev(stress), md_scale(alpha, p));
    const double shifted_norm = md_norm(shifted);
    out.n = (shifted_norm > small) ? md_scale(shifted, 1.0/shifted_norm)
                                   : md_zero();
    const tensor_t n2 = md_square(out.n);
    double cos3theta = sqrt(6.0)*md_ddot(out.n, n2);
    cos3theta = fmin(1.0, fmax(-1.0, cos3theta));
    const double gden = (1.0+m->c)-(1.0-m->c)*cos3theta;
    if (fabs(gden) <= DBL_MIN) return out;
    const double g = 2.0*m->c/gden;
    const double psi = void_ratio-
        (m->e0-m->lambda_c*pow(p/m->P_atm, m->ksi));
    const double alpha_b = g*m->Mc*exp(-m->nb*psi)-m->m;
    const double alpha_d = g*m->Mc*exp( m->nd*psi)-m->m;
    out.b = md_sub(md_scale(out.n, sqrt(2.0/3.0)*alpha_b), alpha);
    const tensor_t d = md_sub(
        md_scale(out.n, sqrt(2.0/3.0)*alpha_d), alpha);
    const double denom_h = md_ddot(md_sub(alpha, alpha_in), out.n);
    out.hardening_gap = denom_h;
    const double b0 = m->G0*m->h0*(1.0-m->ch*void_ratio)/
                      sqrt(p/m->P_atm);
    out.h = (fabs(denom_h) < small) ? 1.0e10 : b0/denom_h;
    const double A = m->A0*(1.0+fmax(md_ddot(fabric, out.n), 0.0));
    out.D = A*md_ddot(d, out.n);
    if (p < 0.05*m->P_atm) {
        /* The published/UW low-pressure regularization was calibrated with
         * pressure expressed in kPa. Hercules input is SI, so convert here
         * instead of making the response depend on the user's unit system. */
        const double p_kpa = p*1.0e-3;
        out.D /= 1.0+exp(7.6349-7.2713*p_kpa);
    }
    out.B = 1.0+1.5*(1.0-m->c)/m->c*g*cos3theta;
    out.C = 3.0*sqrt(1.5)*(1.0-m->c)/m->c*g;
    out.R = md_sub(md_scale(out.n, out.B),
                   md_scale(md_sub(n2, md_iso(1.0/3.0)), out.C));
    out.R = md_add(out.R, md_iso(out.D/3.0));
    out.ok = md_finite(out.n) && md_finite(out.b) && md_finite(out.R) &&
             isfinite(out.h) && isfinite(out.D) && isfinite(out.B) &&
             isfinite(out.C);
    return out;
}

NL_HD static nl_md_stage_t md_stage(
    const manzari_dafalias_constants_t *m, tensor_t stress, tensor_t alpha,
    tensor_t fabric, tensor_t alpha_in, double void_ratio, tensor_t de,
    double K, double G)
{
    nl_md_stage_t out = {};
    const nl_md_dependent_t dep = md_state_dependent(
        m, stress, alpha, fabric, void_ratio, alpha_in);
    if (!dep.ok) return out;
    const double p = md_trace(stress)/3.0+1.0e-2*m->P_atm;
    const tensor_t r = md_scale(md_dev(stress), 1.0/p);
    const tensor_t n2 = md_square(dep.n);
    const double Kp = (2.0/3.0)*p*dep.h*md_ddot(dep.b, dep.n);
    const double denominator = Kp+
        2.0*G*(dep.B-dep.C*md_ddot(dep.n, n2))-
        K*dep.D*md_ddot(dep.n, r);
    const double dev_work = md_ddot(dep.n, md_dev(de));
    const double numerator = 2.0*G*dev_work-
                             K*md_trace(de)*md_ddot(dep.n, r);
    out.loading_numerator = numerator;
    out.denominator = denominator;
    out.h = dep.h;
    out.hardening_gap = dep.hardening_gap;

    if (fabs(denominator) < 1.0e-10) {
        /* Match the UW/OpenSees Modified-Euler neutral-loading branch. When
         * the consistency denominator vanishes, the increment is assigned
         * to plastic strain while stress and directional history remain
         * unchanged. Returning an invalid stage here made Hercules abort a
         * trial that OpenSees handles explicitly. */
        out.raw_dgamma = 0.0;
        out.dgamma = 0.0;
        out.neutral = 1;
        out.dsigma = md_zero();
        out.dalpha = md_zero();
        out.dfabric = md_zero();
        out.dpstrain = de;
        out.ok = md_finite(out.dpstrain) &&
                 isfinite(out.loading_numerator) &&
                 isfinite(out.denominator) && isfinite(out.h) &&
                 isfinite(out.hardening_gap);
        return out;
    }
    out.raw_dgamma = numerator/denominator;
    out.dgamma = out.raw_dgamma;
    if (out.dgamma <= 0.0) {
        out.dgamma = 0.0;
        out.dsigma = md_elastic_increment(de, K, G);
        const tensor_t next = md_add(stress, out.dsigma);
        const double pnext = md_trace(next)/3.0+1.0e-2*m->P_atm;
        if (!(pnext > 0.0)) return out;
        /* Alpha, fabric, and plastic strain are plastic history variables.
         * They must remain fixed during an elastic loading/unloading stage.
         * Moving alpha with the elastic stress ratio artificially translates
         * the yield surface and can hide a subsequent yield crossing. */
        out.dalpha = md_zero();
    } else {
        tensor_t plastic_direction = md_sub(
            md_scale(dep.n, dep.B),
            md_scale(md_sub(n2, md_iso(1.0/3.0)), dep.C));
        tensor_t correction = md_add(md_scale(plastic_direction, 2.0*G),
                                     md_iso(K*dep.D));
        out.dsigma = md_sub(md_elastic_increment(de, K, G),
                            md_scale(correction, out.dgamma));
        out.dalpha = md_scale(dep.b, out.dgamma*(2.0/3.0)*dep.h);
        const double fabric_scale = -out.dgamma*m->cz*fmax(-dep.D, 0.0);
        out.dfabric = md_scale(
            md_add(md_scale(dep.n, m->z_max), fabric), fabric_scale);
        out.dpstrain = md_scale(dep.R, out.dgamma);
    }
    out.ok = md_finite(out.dsigma) && md_finite(out.dalpha) &&
             md_finite(out.dfabric) && md_finite(out.dpstrain) &&
             isfinite(out.dgamma) && isfinite(out.raw_dgamma) &&
             isfinite(out.loading_numerator) &&
             isfinite(out.denominator) && isfinite(out.h) &&
             isfinite(out.hardening_gap);
    return out;
}

/* OpenSees 3.8 integration scheme 5 evaluates one explicit derivative at the
 * beginning of its increment. Hercules applies that stage at the beginning
 * of each fixed Forward-Euler substep. Keep this separate from md_stage: the
 * reference ForwardEuler routine clamps a near-zero denominator instead of
 * taking Modified Euler's neutral-loading branch, applies the Macauley
 * bracket to stress/backstress/fabric evolution, and retains the raw plastic
 * multiplier for its plastic-strain increment. Its local stress-ratio vector
 * is zero initialized by the reference routine, so the two consistency terms
 * involving r vanish; reproducing that detail is required for a meaningful
 * OpenSees scheme-5 comparison. */
NL_HD static nl_md_stage_t md_forward_euler_stage(
    const manzari_dafalias_constants_t *m, tensor_t stress, tensor_t alpha,
    tensor_t fabric, tensor_t alpha_in, double void_ratio, tensor_t de,
    double K, double G)
{
    nl_md_stage_t out = {};
    const nl_md_dependent_t dep = md_state_dependent(
        m, stress, alpha, fabric, void_ratio, alpha_in);
    if (!dep.ok) return out;
    const double p = md_trace(stress)/3.0+1.0e-2*m->P_atm;
    const tensor_t r = md_zero();
    const tensor_t n2 = md_square(dep.n);
    const double Kp = (2.0/3.0)*p*dep.h*md_ddot(dep.b, dep.n);
    double denominator = Kp+
        2.0*G*(dep.B-dep.C*md_ddot(dep.n, n2))-
        K*dep.D*md_ddot(dep.n, r);
    if (fabs(denominator) < 1.0e-10) denominator = 1.0e-10;
    const double numerator =
        2.0*G*md_ddot(dep.n, md_dev(de))-
        K*md_trace(de)*md_ddot(dep.n, r);
    out.loading_numerator = numerator;
    out.denominator = denominator;
    out.h = dep.h;
    out.hardening_gap = dep.hardening_gap;
    out.raw_dgamma = numerator/denominator;
    out.dgamma = fmax(out.raw_dgamma, 0.0);

    tensor_t plastic_direction = md_sub(
        md_scale(dep.n, dep.B),
        md_scale(md_sub(n2, md_iso(1.0/3.0)), dep.C));
    const tensor_t correction = md_add(
        md_scale(plastic_direction, 2.0*G), md_iso(K*dep.D));
    out.dsigma = md_sub(
        md_elastic_increment(de, K, G),
        md_scale(correction, out.dgamma));
    out.dalpha = md_scale(dep.b, out.dgamma*(2.0/3.0)*dep.h);
    const double fabric_scale = -out.dgamma*m->cz*fmax(-dep.D, 0.0);
    out.dfabric = md_scale(
        md_add(md_scale(dep.n, m->z_max), fabric), fabric_scale);
    out.dpstrain = md_scale(dep.R, out.raw_dgamma);
    out.ok = md_finite(out.dsigma) && md_finite(out.dalpha) &&
             md_finite(out.dfabric) && md_finite(out.dpstrain) &&
             isfinite(out.dgamma) && isfinite(out.raw_dgamma) &&
             isfinite(out.loading_numerator) &&
             isfinite(out.denominator) && isfinite(out.h) &&
             isfinite(out.hardening_gap);
    return out;
}

/* OpenSees falls back to a segment search when its local Newton correction
 * exhausts 50 iterations. Keep the same admissible-state idea, but require a
 * genuine yield-function bracket and preserve total strain by mapping the
 * removed stress through the elastic compliance. */
NL_HD static int md_stress_segment_projection(
    const manzari_dafalias_constants_t *m, tensor_t anchor_stress,
    tensor_t proposed_stress, tensor_t alpha,
    tensor_t proposed_plastic_strain, tensor_t *stress,
    tensor_t *plastic_strain,
    double K, double G, double tolerance, double *residual)
{
    const tensor_t stress_delta = md_sub(proposed_stress, anchor_stress);
    double lower = 0.0;
    double upper = 1.0;
    double flower = md_yield(m, anchor_stress, alpha);
    double fupper = md_yield(m, proposed_stress, alpha);
    if (!isfinite(flower) || !isfinite(fupper) ||
        flower*fupper > 0.0) return 0;

    tensor_t projected = proposed_stress;
    double fprojected = fupper;
    for (int iteration = 0; iteration < 80; iteration++) {
        const double middle = 0.5*(lower+upper);
        projected = md_add(anchor_stress, md_scale(stress_delta, middle));
        fprojected = md_yield(m, projected, alpha);
        if (!isfinite(fprojected)) return 0;
        if (fabs(fprojected) < tolerance) break;
        if ((flower <= 0.0 && fprojected <= 0.0) ||
            (flower >= 0.0 && fprojected >= 0.0)) {
            lower = middle;
            flower = fprojected;
        } else {
            upper = middle;
        }
    }
    if (!(fabs(fprojected) < tolerance)) return 0;

    *stress = projected;
    *plastic_strain = md_add(
        proposed_plastic_strain,
        md_elastic_compliance(md_sub(proposed_stress, projected), K, G));
    *residual = fabs(fprojected);
    return 1;
}

NL_HD static int md_stress_correction(
    const manzari_dafalias_constants_t *m, tensor_t anchor_stress,
    tensor_t *stress,
    tensor_t *alpha, tensor_t fabric, tensor_t *alpha_in,
    double void_ratio, tensor_t direction_hint, tensor_t *plastic_strain,
    double K, double G, double tolerance, int maximum_iterations,
    double *residual, int32_t *iterations_used, int32_t *projection_used)
{
    *iterations_used = 0;
    *projection_used = 0;
    const tensor_t proposed_stress = *stress;
    const tensor_t proposed_alpha = *alpha;
    const tensor_t proposed_plastic_strain = *plastic_strain;
    if (md_trace(*stress)/3.0 < md_minimum_mean_stress(m)) {
        if (!md_project_candidate_to_pressure_floor(
                m, stress, plastic_strain, alpha, alpha_in, K, G)) return 0;
        *residual = 0.0;
        return 1;
    }
    double f = md_yield(m, *stress, *alpha);
    *residual = isfinite(f) ? fabs(f) : DBL_MAX;
    if (!isfinite(f)) return 0;
    if (fabs(f) < tolerance) {
        return 1;
    }

    const int iterations = maximum_iterations < 50 ? maximum_iterations : 50;
    for (int iteration = 0; iteration < iterations; iteration++) {
        *iterations_used = iteration+1;
        const double p = md_trace(*stress)/3.0+1.0e-2*m->P_atm;
        if (!(p > 0.0)) return 0;
        const nl_md_dependent_t dep = md_state_dependent(
            m, *stress, *alpha, fabric, void_ratio, *alpha_in);
        if (!dep.ok) return 0;
        const tensor_t r = md_scale(md_dev(*stress), 1.0/p);
        const tensor_t dsigma_p = md_add(md_scale(md_dev(dep.R), 2.0*G),
                                         md_iso(K*md_trace(dep.R)));
        const tensor_t abar = md_scale(dep.b, (2.0/3.0)*dep.h);
        tensor_t dfds = md_sub(dep.n,
            md_iso(md_ddot(dep.n, r)/3.0));
        const tensor_t dfdalpha = md_scale(dep.n, -p);
        const double denominator = md_ddot(dfds, dsigma_p)-
                                   md_ddot(dfdalpha, abar);
        if (fabs(denominator) <= DBL_MIN) {
            /* At the exact center of the small kinematic yield surface,
             * n and the Newton denominator are undefined/zero. If Modified
             * Euler crossed into that center, recover the loading direction
             * from the elastic trial and project radially to |f|=0. */
            const tensor_t hint = md_dev(direction_hint);
            const double hint_norm = md_norm(hint);
            if (!(f < -tolerance) || hint_norm <= DBL_MIN) return 0;
            const tensor_t old_stress = *stress;
            const double radius = sqrt(2.0/3.0)*m->m*p;
            const tensor_t target_dev = md_add(
                md_scale(*alpha, p), md_scale(hint, radius/hint_norm));
            const tensor_t next_stress = md_add(
                target_dev, md_iso(md_trace(old_stress)/3.0));
            const double next_f = md_yield(m, next_stress, *alpha);
            if (!isfinite(next_f) || fabs(next_f) >= tolerance) return 0;
            *stress = next_stress;
            *plastic_strain = md_add(
                *plastic_strain,
                md_elastic_compliance(
                    md_sub(old_stress, next_stress), K, G));
            *residual = fabs(next_f);
            return 1;
        }
        const double dlambda = f/denominator;
        const tensor_t old_stress = *stress;
        const tensor_t old_alpha = *alpha;
        tensor_t next_stress = old_stress;
        tensor_t next_alpha = old_alpha;
        double next_f = f;
        int improved = 0;

        /* A full Newton correction can cross the small kinematic yield
         * surface. Backtrack until |f| actually decreases instead of
         * accepting a large negative residual or oscillating across f=0. */
        double line_scale = 1.0;
        for (int line = 0; line < 20; line++) {
            next_stress = md_sub(
                old_stress, md_scale(dsigma_p, line_scale*dlambda));
            next_alpha = md_add(
                old_alpha, md_scale(abar, line_scale*dlambda));
            next_f = md_yield(m, next_stress, next_alpha);
            if (isfinite(next_f) && fabs(next_f) < fabs(f)) {
                improved = 1;
                break;
            }
            line_scale *= 0.5;
        }

        if (!improved) {
            const double grad_norm = md_ddot(dfds, dfds);
            if (grad_norm <= DBL_MIN) return 0;
            const double gradient_step = f/grad_norm;
            line_scale = 1.0;
            for (int line = 0; line < 20; line++) {
                next_stress = md_sub(
                    old_stress,
                    md_scale(dfds, line_scale*gradient_step));
                next_alpha = old_alpha;
                next_f = md_yield(m, next_stress, next_alpha);
                if (isfinite(next_f) && fabs(next_f) < fabs(f)) {
                    improved = 1;
                    break;
                }
                line_scale *= 0.5;
            }
            if (!improved) return 0;
        }
        if (md_trace(next_stress)/3.0 < md_minimum_mean_stress(m)) {
            tensor_t floor_stress = old_stress;
            tensor_t floor_alpha = old_alpha;
            if (!md_project_candidate_to_pressure_floor(
                    m, &floor_stress, plastic_strain, &floor_alpha,
                    alpha_in, K, G))
                return 0;
            *stress = floor_stress;
            *alpha = floor_alpha;
            *residual = 0.0;
            return 1;
        }
        *stress = next_stress;
        *alpha = next_alpha;
        /* Total strain is fixed during projection. Therefore the plastic
         * correction is the elastic-compliance image of the stress removed,
         * including when the robust yield-gradient fallback was selected. */
        *plastic_strain = md_add(
            *plastic_strain,
            md_elastic_compliance(md_sub(old_stress, next_stress), K, G));
        f = next_f;
        *residual = fabs(f);
        if (fabs(f) < tolerance) return 1;
    }
    if (!md_stress_segment_projection(
            m, anchor_stress, proposed_stress, proposed_alpha,
            proposed_plastic_strain, stress, plastic_strain,
            K, G, tolerance, residual)) return 0;
    *alpha = proposed_alpha;
    *projection_used = 1;
    return 1;
}

/* The integration scheme is a compile-time parameter so each production kernel
 * contains exactly one integrator. Modified Euler holds two nl_md_stage_t
 * estimates plus the stress-correction working set and allocates all 255
 * registers with heavy spilling; Forward Euler needs a fraction of that. When
 * the scheme was a runtime branch, ptxas had to allocate for the worse of the
 * two, so a Forward-Euler run paid Modified Euler's spill traffic. The host
 * material-point driver keeps a runtime-dispatching wrapper below. */
template <nl_md_integration_scheme_t MD_SCHEME,
          nl_md_forward_euler_substeps_t MD_FE_SUBSTEPS>
NL_HD static int nl_update_manzari_dafalias_scheme(
    const nlconstants_t *__restrict__ c,
    const nonlinear_gpu_config_t *__restrict__ cfg,
    tensor_t strain, tensor_t strain_prev, tensor_t static_stress,
    tensor_t *__restrict__ dynamic_stress,
    tensor_t *__restrict__ plastic_strain,
    tensor_t *__restrict__ alpha, tensor_t *__restrict__ alpha_in,
    tensor_t *__restrict__ fabric, double *__restrict__ initialized,
    double *__restrict__ void_ratio,
    double *__restrict__ plastic_multiplier,
    double *__restrict__ yield_error,
    double *__restrict__ integration_error,
    int32_t *__restrict__ accepted_substeps,
    nl_md_diagnostics_t *__restrict__ diagnostics)
{
    const manzari_dafalias_constants_t *m = &c->md;
    const tensor_t de_h = md_sub(strain, strain_prev);
    const tensor_t de = md_scale(de_h, -1.0);
    tensor_t stress = md_scale(md_add(*dynamic_stress, static_stress), -1.0);
    tensor_t pstrain = md_scale(*plastic_strain, -1.0);
    tensor_t alpha_work = *alpha;
    tensor_t alpha_in_work = *alpha_in;
    tensor_t fabric_work = *fabric;
    double initialized_work = *initialized;
    *diagnostics = nl_md_diagnostics_t{};
    diagnostics->effective_mean_stress = md_trace(stress)/3.0;
    *void_ratio = m->e_init-(1.0+m->e_init)*md_trace(md_scale(strain, -1.0));
    *plastic_multiplier = 0.0;
    *yield_error = 0.0;
    *integration_error = 0.0;
    *accepted_substeps = 0;

    if (!(initialized_work > 0.5)) {
        const double p = md_trace(stress)/3.0+1.0e-2*m->P_atm;
        if (!(p > 0.0)) {
            diagnostics->failure_reason = NL_MD_FAILURE_INITIAL_PRESSURE;
            return 0;
        }
        /* Treat a supplied geostatic tensor as a pre-existing consolidated
         * state: center the small kinematic yield surface on its stress ratio
         * so self weight is not mistaken for a new plastic load increment.
         * For the usual isotropic material-point start this is exactly zero. */
        alpha_work = md_scale(md_dev(stress), 1.0/p);
        alpha_in_work = alpha_work;
        fabric_work = md_zero();
        initialized_work = 1.0;
    }
    if (md_norm(de) == 0.0) {
        *alpha = alpha_work;
        *alpha_in = alpha_in_work;
        *fabric = fabric_work;
        *initialized = initialized_work;
        return 1;
    }
    if (!(*void_ratio > -0.99) || cfg->no_substeps <= 0 ||
        !md_valid_integration_scheme(cfg->md_integration_scheme) ||
        !md_valid_forward_euler_substeps(
            cfg->md_forward_euler_substeps) ||
        !md_valid_recovery_mode(cfg->md_recovery_mode) ||
        /* Guards a launch that dispatched to the wrong specialization. */
        cfg->md_integration_scheme != MD_SCHEME ||
        cfg->md_forward_euler_substeps != MD_FE_SUBSTEPS ||
        (MD_SCHEME == NL_MD_INTEGRATION_MODIFIED_EULER &&
         MD_FE_SUBSTEPS != NL_MD_FORWARD_EULER_SUBSTEPS_FIXED) ||
        (MD_SCHEME == NL_MD_INTEGRATION_FORWARD_EULER &&
         cfg->md_recovery_mode != NL_MD_RECOVERY_STRICT) ||
        !(cfg->backbone_error_tolerance > 0.0) ||
        !(cfg->error_tolerance > 0.0)) {
        diagnostics->failure_reason = NL_MD_FAILURE_INVALID_CONTROLS;
        return 0;
    }

    double K, G;
    md_elastic_moduli(m, stress, &K, &G);
    if (!(K > 0.0) || !(G > 0.0) || !isfinite(K) || !isfinite(G)) {
        diagnostics->failure_reason = NL_MD_FAILURE_ELASTIC_MODULI;
        return 0;
    }
    const double minimum_mean_stress = md_minimum_mean_stress(m);
    if (md_trace(stress)/3.0 < minimum_mean_stress) {
        const tensor_t old_stress = stress;
        stress = md_add(md_dev(stress), md_iso(minimum_mean_stress));
        pstrain = md_add(
            pstrain,
            md_elastic_compliance(md_sub(old_stress, stress), K, G));
        if (!md_finite(stress) || !md_finite(pstrain)) {
            diagnostics->failure_reason = NL_MD_FAILURE_INITIAL_PRESSURE;
            return 0;
        }
    }
    const tensor_t elastic_trial = md_elastic_increment(de, K, G);
    if (md_ddot(md_sub(alpha_work, alpha_in_work), elastic_trial) < 0.0)
        alpha_in_work = alpha_work;

    const tensor_t trial = md_add(stress, elastic_trial);
    const double trial_mean_stress = md_trace(trial)/3.0;
    const double trial_f = md_yield(m, trial, alpha_work);
    if (trial_mean_stress >= minimum_mean_stress &&
        trial_f <= cfg->error_tolerance) {
        stress = trial;
        const tensor_t dynamic_next = md_sub(
            md_scale(stress, -1.0), static_stress);
        if (!md_finite(dynamic_next)) {
            diagnostics->failure_reason = NL_MD_FAILURE_FINAL_STATE;
            return 0;
        }
        *dynamic_stress = dynamic_next;
        *alpha = alpha_work;
        *alpha_in = alpha_in_work;
        *fabric = fabric_work;
        *initialized = initialized_work;
        *yield_error = fmax(trial_f, 0.0);
        *accepted_substeps = 1;
        diagnostics->effective_mean_stress = md_trace(stress)/3.0;
        return 1;
    }

    double T = 0.0;
    double dT = 1.0;
    const double current_f = md_yield(m, stress, alpha_work);
    if (current_f < -cfg->error_tolerance &&
        trial_f > cfg->error_tolerance) {
        double elastic_ratio = 0.0;
        if (!md_elastic_yield_intersection(
                m, stress, elastic_trial, alpha_work, cfg->error_tolerance,
                &elastic_ratio)) {
            diagnostics->failure_reason = NL_MD_FAILURE_YIELD_INTERSECTION;
            return 0;
        }
        stress = md_add(stress, md_scale(elastic_trial, elastic_ratio));
        T = elastic_ratio;
        dT = 1.0-elastic_ratio;
    }

    if (MD_SCHEME == NL_MD_INTEGRATION_FORWARD_EULER) {
        /* OpenSees integration scheme 5 supplies the first-order update used
         * by every substep. The default fixed policy divides the plastic
         * remainder after the elastic/yield intersection into exactly
         * no_substeps equal parts, preserving the historical and canonical
         * sensitivity-study behavior. The opt-in demand policy instead sizes
         * the count from an MD-specific stress-ratio/deviatoric-strain
         * heuristic and treats no_substeps as a cap. There is no embedded-error
         * estimate, adaptive retry, or yield-surface stress correction. */
        const double plastic_fraction = dT;
        int32_t requested = cfg->no_substeps;
        if (MD_FE_SUBSTEPS == NL_MD_FORWARD_EULER_SUBSTEPS_DEMAND) {
            const tensor_t de_plastic = md_scale(de, plastic_fraction);
            const tensor_t trial_plastic = md_add(
                stress, md_elastic_increment(de_plastic, K, G));
            const double p_now = md_trace(stress)/3.0+1.0e-2*m->P_atm;
            const double p_trial =
                md_trace(trial_plastic)/3.0+1.0e-2*m->P_atm;
            requested = 1;
            if (p_now > 0.0 && p_trial > 0.0) {
                const tensor_t ratio_change = md_sub(
                    md_scale(md_dev(trial_plastic), 1.0/p_trial),
                    md_scale(md_dev(stress), 1.0/p_now));
                requested = md_substep_request(
                    sqrt(1.5*md_ddot(ratio_change, ratio_change))/
                        NL_MD_FE_RATIO_STEP,
                    cfg->no_substeps);
            }
            const double dev_octahedral =
                2.0*md_norm(md_dev(de_plastic))/sqrt(3.0);
            const int32_t strain_request = md_substep_request(
                dev_octahedral/NL_MD_FE_STRAIN_STEP, cfg->no_substeps);
            if (strain_request > requested) requested = strain_request;
        }
        const double fixed_fraction = plastic_fraction/(double)requested;
        for (int32_t substep = 0; substep < requested; substep++) {
            const double step_fraction = substep == requested-1
                ? fmax(1.0-T, 0.0) : fixed_fraction;
            diagnostics->attempts = substep+1;
            diagnostics->substep_fraction = step_fraction;
            diagnostics->progress_fraction = T;
            diagnostics->effective_mean_stress = md_trace(stress)/3.0;
            double step_K = 0.0;
            double step_G = 0.0;
            md_elastic_moduli(m, stress, &step_K, &step_G);
            if (!(step_K > 0.0) || !(step_G > 0.0) ||
                !isfinite(step_K) || !isfinite(step_G)) {
                diagnostics->failure_reason = NL_MD_FAILURE_ELASTIC_MODULI;
                return 0;
            }
            const tensor_t de_step = md_scale(de, step_fraction);
            const tensor_t elastic_step = md_elastic_increment(
                de_step, step_K, step_G);
            if (md_ddot(
                    md_sub(alpha_work, alpha_in_work), elastic_step) < 0.0)
                alpha_in_work = alpha_work;
            const tensor_t strain_at_end = md_add(
                md_scale(strain_prev, -1.0),
                md_scale(de, T+step_fraction));
            const double e_step = m->e_init-
                (1.0+m->e_init)*md_trace(strain_at_end);
            const nl_md_stage_t first = md_forward_euler_stage(
                m, stress, alpha_work, fabric_work, alpha_in_work,
                e_step, de_step, step_K, step_G);
            diagnostics->first_dgamma = first.dgamma;
            diagnostics->first_raw_dgamma = first.raw_dgamma;
            diagnostics->first_denominator = first.denominator;
            diagnostics->first_h = first.h;
            diagnostics->first_hardening_gap = first.hardening_gap;
            if (!first.ok) {
                diagnostics->failure_reason = NL_MD_FAILURE_FIRST_STAGE;
                return 0;
            }

            const tensor_t candidate_stress = md_add(stress, first.dsigma);
            if (md_trace(candidate_stress)/3.0 < minimum_mean_stress) {
                if (!md_complete_at_pressure_floor(
                        m, stress, de_step, &stress, &pstrain,
                        &alpha_work, &alpha_in_work, step_K, step_G)) {
                    diagnostics->failure_reason = NL_MD_FAILURE_FINAL_STATE;
                    return 0;
                }
            } else {
                stress = candidate_stress;
                pstrain = md_add(pstrain, first.dpstrain);
                alpha_work = md_add(alpha_work, first.dalpha);
                fabric_work = md_add(fabric_work, first.dfabric);
            }
            *plastic_multiplier += first.dgamma;
            (*accepted_substeps)++;
            T += step_fraction;
        }
        *integration_error = 0.0;

        const double final_f = md_yield(m, stress, alpha_work);
        const tensor_t dynamic_next = md_sub(
            md_scale(stress, -1.0), static_stress);
        const tensor_t plastic_next = md_scale(pstrain, -1.0);
        *yield_error = isfinite(final_f) ? fabs(final_f) : DBL_MAX;
        diagnostics->correction_residual = *yield_error;
        diagnostics->effective_mean_stress = md_trace(stress)/3.0;
        if (!md_finite(dynamic_next) || !md_finite(plastic_next) ||
            !md_finite(alpha_work) || !md_finite(alpha_in_work) ||
            !md_finite(fabric_work) || !isfinite(*yield_error)) {
            diagnostics->failure_reason = NL_MD_FAILURE_FINAL_STATE;
            return 0;
        }
        *dynamic_stress = dynamic_next;
        *plastic_strain = plastic_next;
        *alpha = alpha_work;
        *alpha_in = alpha_in_work;
        *fabric = fabric_work;
        *initialized = initialized_work;
        return 1;
    }

    /* OpenSees stops reducing the normalized constitutive fraction at 1e-6
     * and then accepts a special first-order state when its embedded error is
     * too large. The opensees policy uses that floor; all other policies allow
     * up to three additional decades when the imposed strain increment can
     * represent them accurately. */
    const double strain_roundoff =
        64.0*DBL_EPSILON*(1.0+md_norm(strain_prev));
    const double dT_min = cfg->md_recovery_mode == NL_MD_RECOVERY_OPENSEES
        ? 1.0e-6
        : fmin(1.0e-6,
               fmax(1.0e-9,
                    strain_roundoff/fmax(md_norm(de), DBL_MIN)));
    const int force_accept = md_force_accept_enabled(cfg->md_recovery_mode);
    const int unlimited_substeps =
        md_unlimited_substeps_enabled(cfg->md_recovery_mode);
    /* no_substeps caps accepted work. Rejected error-control trials are
     * bounded separately so one difficult point cannot spin forever. */
    int64_t attempts = 0;
    while (T < 1.0-DBL_EPSILON) {
        attempts++;
        if (!unlimited_substeps && attempts > 4LL*cfg->no_substeps) {
            diagnostics->failure_reason = NL_MD_FAILURE_ATTEMPT_LIMIT;
            diagnostics->attempts = attempts > INT32_MAX
                ? INT32_MAX : (int32_t)attempts;
            return 0;
        }
        diagnostics->attempts = attempts > INT32_MAX
            ? INT32_MAX : (int32_t)attempts;
        if (!unlimited_substeps &&
            *accepted_substeps >= cfg->no_substeps) {
            diagnostics->failure_reason = NL_MD_FAILURE_SUBSTEP_LIMIT;
            return 0;
        }
        dT = fmin(dT, 1.0-T);
        diagnostics->substep_fraction = dT;
        diagnostics->progress_fraction = T;
        diagnostics->effective_mean_stress = md_trace(stress)/3.0;
        const tensor_t de_step = md_scale(de, dT);
        const tensor_t elastic_step = md_elastic_increment(de_step, K, G);
        /* alpha_in is the backstress ratio at the most recent reversal.
         * The global predictor applies the OpenSees reversal test once before
         * integration, but alpha evolves during accepted Modified-Euler
         * substeps. Reapply the same test to the current accepted state so a
         * multiaxial loading-direction switch inside the global increment
         * cannot leave (alpha-alpha_in):n negative during outward loading.
         * Re-anchoring before the stage restores the limiting h -> infinity
         * response at reversal and avoids a false consistency-denominator
         * pole; rejected candidates still cannot commit this local state. */
        if (md_ddot(md_sub(alpha_work, alpha_in_work), elastic_step) < 0.0)
            alpha_in_work = alpha_work;
        const tensor_t strain_at_end = md_add(
            md_scale(strain_prev, -1.0), md_scale(de, T+dT));
        const double e_step = m->e_init-
            (1.0+m->e_init)*md_trace(strain_at_end);
        nl_md_stage_t first = md_stage(
            m, stress, alpha_work, fabric_work, alpha_in_work,
            e_step, de_step, K, G);
        if (!first.ok) {
            if (dT > dT_min) {
                dT = fmax(0.5*dT, dT_min);
                continue;
            }
            if (force_accept) {
                tensor_t forced_stress = stress;
                tensor_t forced_pstrain = md_add(pstrain, de_step);
                double forced_residual = 0.0;
                if (md_force_accept_candidate(
                        m, forced_stress, forced_pstrain, fabric_work,
                        &stress, &pstrain, &alpha_work, &alpha_in_work,
                        &fabric_work, K, G, &forced_residual)) {
                    diagnostics->force_accepts++;
                    diagnostics->local_error = 1.0;
                    diagnostics->correction_residual = forced_residual;
                    *integration_error = fmax(*integration_error, 1.0);
                    *yield_error = fmax(*yield_error, forced_residual);
                    (*accepted_substeps)++;
                    T += dT;
                    if (T >= 1.0-DBL_EPSILON) break;
                    dT = fmin(fmax(2.0*dT, dT_min), 1.0-T);
                    continue;
                }
            }
            diagnostics->failure_reason = NL_MD_FAILURE_FIRST_STAGE;
            return 0;
        }
        /* A positive consistency numerator is outward loading. Its
         * loading-side hardening gap must not be negative: that combination
         * means the evolving multiaxial direction has passed the stored
         * reversal anchor even if the broader trial-direction dot product is
         * numerically inconclusive. Re-anchor and recompute the stage instead
         * of letting a negative h create a denominator pole and a spurious
         * negative plastic multiplier. */
        if (first.hardening_gap < 0.0 &&
            first.loading_numerator > 0.0) {
            alpha_in_work = alpha_work;
            first = md_stage(
                m, stress, alpha_work, fabric_work, alpha_in_work,
                e_step, de_step, K, G);
            if (!first.ok) {
                if (dT > dT_min) {
                    dT = fmax(0.5*dT, dT_min);
                    continue;
                }
                if (force_accept) {
                    tensor_t forced_stress = stress;
                    tensor_t forced_pstrain = md_add(pstrain, de_step);
                    double forced_residual = 0.0;
                    if (md_force_accept_candidate(
                            m, forced_stress, forced_pstrain, fabric_work,
                            &stress, &pstrain, &alpha_work,
                            &alpha_in_work, &fabric_work, K, G,
                            &forced_residual)) {
                        diagnostics->force_accepts++;
                        diagnostics->local_error = 1.0;
                        diagnostics->correction_residual = forced_residual;
                        *integration_error = fmax(*integration_error, 1.0);
                        *yield_error = fmax(*yield_error, forced_residual);
                        (*accepted_substeps)++;
                        T += dT;
                        if (T >= 1.0-DBL_EPSILON) break;
                        dT = fmin(
                            fmax(2.0*dT, dT_min), 1.0-T);
                        continue;
                    }
                }
                diagnostics->failure_reason = NL_MD_FAILURE_FIRST_STAGE;
                return 0;
            }
        }
        diagnostics->first_dgamma = first.dgamma;
        diagnostics->first_raw_dgamma = first.raw_dgamma;
        diagnostics->first_denominator = first.denominator;
        diagnostics->first_h = first.h;
        diagnostics->first_hardening_gap = first.hardening_gap;
        diagnostics->second_dgamma = 0.0;
        diagnostics->second_raw_dgamma = 0.0;
        diagnostics->second_denominator = 0.0;
        diagnostics->second_h = 0.0;
        diagnostics->second_hardening_gap = 0.0;
        const tensor_t predictor_stress = md_add(stress, first.dsigma);
        const double predictor_mean_stress = md_trace(predictor_stress)/3.0;
        if (predictor_mean_stress < minimum_mean_stress) {
            diagnostics->effective_mean_stress =
                predictor_mean_stress;
            if (dT <= dT_min) {
                if (!md_complete_at_pressure_floor(
                        m, stress, md_scale(de, 1.0-T), &stress, &pstrain,
                        &alpha_work, &alpha_in_work, K, G)) {
                    diagnostics->failure_reason = NL_MD_FAILURE_FINAL_STATE;
                    return 0;
                }
                diagnostics->correction_residual = 0.0;
                diagnostics->effective_mean_stress = minimum_mean_stress;
                (*accepted_substeps)++;
                T = 1.0;
                break;
            }
            const double current_mean_stress = md_trace(stress)/3.0;
            const double pressure_ratio =
                (current_mean_stress-minimum_mean_stress)/
                fmax(current_mean_stress-predictor_mean_stress, DBL_MIN);
            dT = fmax(fmax(0.8*pressure_ratio, 0.1)*dT, dT_min);
            continue;
        }
        const double predictor_f = md_yield(
            m, predictor_stress, alpha_work);
        if (!isfinite(predictor_f)) {
            if (dT > dT_min) {
                dT = fmax(0.5*dT, dT_min);
                continue;
            }
            if (force_accept) {
                tensor_t forced_stress = predictor_stress;
                tensor_t forced_pstrain = md_add(pstrain, first.dpstrain);
                tensor_t forced_fabric = md_add(
                    fabric_work, first.dfabric);
                double forced_residual = 0.0;
                if (md_force_accept_candidate(
                        m, forced_stress, forced_pstrain, forced_fabric,
                        &stress, &pstrain, &alpha_work, &alpha_in_work,
                        &fabric_work, K, G, &forced_residual)) {
                    diagnostics->force_accepts++;
                    diagnostics->local_error = 1.0;
                    diagnostics->correction_residual = forced_residual;
                    *plastic_multiplier += first.dgamma;
                    *integration_error = fmax(*integration_error, 1.0);
                    *yield_error = fmax(*yield_error, forced_residual);
                    (*accepted_substeps)++;
                    T += dT;
                    if (T >= 1.0-DBL_EPSILON) break;
                    dT = fmin(fmax(2.0*dT, dT_min), 1.0-T);
                    continue;
                }
            }
            diagnostics->failure_reason = NL_MD_FAILURE_FIRST_STAGE;
            return 0;
        }

        /* Modified Euler is a hybrid elastic/plastic update. When the first
         * stage detects unloading (a nonpositive plastic multiplier), its
         * elastic predictor is already the complete substep solution if that
         * predictor remains admissible. Do not evaluate a plastic second
         * stage in that case. Near a loading switch, the consistency
         * denominator can cross zero over an arbitrarily small elastic
         * increment; comparing that discontinuous plastic stage with the
         * elastic first stage creates an O(1) embedded-error estimate that
         * cannot be cured by reducing dT. The next accepted substep starts
         * from the other side of the switch and dispatches plastic flow when
         * its own first-stage multiplier is positive. */
        if (first.dgamma == 0.0 && !first.neutral &&
            predictor_f <= cfg->error_tolerance) {
            stress = predictor_stress;
            diagnostics->local_error = 0.0;
            diagnostics->correction_residual = 0.0;
            *yield_error = fmax(*yield_error, fmax(predictor_f, 0.0));
            (*accepted_substeps)++;
            T += dT;
            if (T >= 1.0-DBL_EPSILON) break;
            dT = fmin(fmax(2.0*dT, dT_min), 1.0-T);
            continue;
        }
        const nl_md_stage_t second = md_stage(
            m, predictor_stress, md_add(alpha_work, first.dalpha),
            md_add(fabric_work, first.dfabric), alpha_in_work,
            e_step, de_step, K, G);
        if (!second.ok) {
            if (dT > dT_min) {
                dT = fmax(0.5*dT, dT_min);
                continue;
            }
            if (force_accept) {
                tensor_t forced_stress = predictor_stress;
                tensor_t forced_pstrain = md_add(pstrain, first.dpstrain);
                tensor_t forced_fabric = md_add(
                    fabric_work, first.dfabric);
                double forced_residual = 0.0;
                if (md_force_accept_candidate(
                        m, forced_stress, forced_pstrain, forced_fabric,
                        &stress, &pstrain, &alpha_work, &alpha_in_work,
                        &fabric_work, K, G, &forced_residual)) {
                    diagnostics->force_accepts++;
                    diagnostics->local_error = 1.0;
                    diagnostics->correction_residual = forced_residual;
                    *plastic_multiplier += first.dgamma;
                    *integration_error = fmax(*integration_error, 1.0);
                    *yield_error = fmax(*yield_error, forced_residual);
                    (*accepted_substeps)++;
                    T += dT;
                    if (T >= 1.0-DBL_EPSILON) break;
                    dT = fmin(fmax(2.0*dT, dT_min), 1.0-T);
                    continue;
                }
            }
            diagnostics->failure_reason = NL_MD_FAILURE_SECOND_STAGE;
            return 0;
        }
        diagnostics->second_dgamma = second.dgamma;
        diagnostics->second_raw_dgamma = second.raw_dgamma;
        diagnostics->second_denominator = second.denominator;
        diagnostics->second_h = second.h;
        diagnostics->second_hardening_gap = second.hardening_gap;

        /* OpenSees checks the averaged Modified-Euler candidate before its
         * embedded-error test. A candidate below the physical pressure floor
         * is a cutoff event, not a smooth integration-error event. Reduce the
         * fraction until the crossing is resolved, then complete the
         * unresolved strain on the cutoff state at the minimum fraction. */
        tensor_t candidate_stress = md_add(stress, md_scale(
            md_add(first.dsigma, second.dsigma), 0.5));
        tensor_t candidate_alpha = md_add(alpha_work, md_scale(
            md_add(first.dalpha, second.dalpha), 0.5));
        tensor_t candidate_alpha_in = alpha_in_work;
        tensor_t candidate_fabric = md_add(fabric_work, md_scale(
            md_add(first.dfabric, second.dfabric), 0.5));
        tensor_t candidate_pstrain = md_add(pstrain, md_scale(
            md_add(first.dpstrain, second.dpstrain), 0.5));
        const double candidate_multiplier = *plastic_multiplier+
            0.5*(first.dgamma+second.dgamma);
        const double candidate_mean_stress = md_trace(candidate_stress)/3.0;
        if (candidate_mean_stress < minimum_mean_stress) {
            diagnostics->effective_mean_stress = candidate_mean_stress;
            if (dT <= dT_min) {
                if (!md_complete_at_pressure_floor(
                        m, stress, md_scale(de, 1.0-T), &stress, &pstrain,
                        &alpha_work, &alpha_in_work, K, G)) {
                    diagnostics->failure_reason = NL_MD_FAILURE_FINAL_STATE;
                    return 0;
                }
                diagnostics->correction_residual = 0.0;
                diagnostics->effective_mean_stress = minimum_mean_stress;
                (*accepted_substeps)++;
                T = 1.0;
                break;
            }
            const double current_mean_stress = md_trace(stress)/3.0;
            const double pressure_ratio =
                (current_mean_stress-minimum_mean_stress)/
                fmax(current_mean_stress-candidate_mean_stress, DBL_MIN);
            dT = fmax(fmax(0.8*pressure_ratio, 0.1)*dT, dT_min);
            continue;
        }

        const double stress_norm = md_norm(stress);
        const double error_denominator = stress_norm < 0.5
            ? 1.0 : 2.0*stress_norm;
        const double step_error =
            md_norm(md_sub(second.dsigma, first.dsigma))/error_denominator;
        diagnostics->local_error = step_error;
        if (!isfinite(step_error)) {
            if (dT > dT_min) {
                dT = fmax(0.5*dT, dT_min);
                continue;
            }
            if (force_accept) {
                double forced_residual = 0.0;
                if (md_force_accept_candidate(
                        m, candidate_stress, candidate_pstrain,
                        candidate_fabric, &stress, &pstrain, &alpha_work,
                        &alpha_in_work, &fabric_work, K, G,
                        &forced_residual)) {
                    diagnostics->force_accepts++;
                    diagnostics->local_error = 1.0;
                    diagnostics->correction_residual = forced_residual;
                    *plastic_multiplier = candidate_multiplier;
                    *integration_error = fmax(*integration_error, 1.0);
                    *yield_error = fmax(*yield_error, forced_residual);
                    (*accepted_substeps)++;
                    T += dT;
                    if (T >= 1.0-DBL_EPSILON) break;
                    dT = fmin(fmax(2.0*dT, dT_min), 1.0-T);
                    continue;
                }
            }
            diagnostics->failure_reason = NL_MD_FAILURE_NONFINITE_LOCAL_ERROR;
            return 0;
        }
        if (step_error > cfg->backbone_error_tolerance) {
            if (dT <= dT_min) {
                *integration_error = step_error;
                if (force_accept) {
                    double forced_residual = 0.0;
                    if (md_force_accept_candidate(
                            m, candidate_stress, candidate_pstrain,
                            candidate_fabric, &stress, &pstrain,
                            &alpha_work, &alpha_in_work, &fabric_work,
                            K, G, &forced_residual)) {
                        diagnostics->force_accepts++;
                        diagnostics->correction_residual = forced_residual;
                        *plastic_multiplier = candidate_multiplier;
                        *yield_error = fmax(
                            *yield_error, forced_residual);
                        (*accepted_substeps)++;
                        T += dT;
                        if (T >= 1.0-DBL_EPSILON) break;
                        dT = fmin(
                            fmax(2.0*dT, dT_min), 1.0-T);
                        continue;
                    }
                }
                diagnostics->failure_reason =
                    NL_MD_FAILURE_LOCAL_ERROR_MIN_STEP;
                return 0;
            }
            const double q = fmax(0.8*sqrt(
                cfg->backbone_error_tolerance/step_error), 0.1);
            dT = fmax(q*dT, dT_min);
            continue;
        }

        /* A failed yield correction must not modify committed history;
         * reject it, reduce the constitutive fraction, and retry from the
         * same accepted state. */
        double correction_residual = 0.0;
        int32_t correction_iterations = 0;
        int32_t projection_used = 0;
        const tensor_t force_candidate_stress = candidate_stress;
        const tensor_t force_candidate_pstrain = candidate_pstrain;
        const tensor_t force_candidate_fabric = candidate_fabric;
        if (!md_stress_correction(
                m, stress, &candidate_stress, &candidate_alpha,
                candidate_fabric, &candidate_alpha_in, e_step,
                elastic_trial, &candidate_pstrain, K, G,
                cfg->error_tolerance, cfg->no_substeps,
                &correction_residual, &correction_iterations,
                &projection_used)) {
            diagnostics->correction_residual = correction_residual;
            diagnostics->correction_iterations =
                correction_iterations > diagnostics->correction_iterations
                    ? correction_iterations
                    : diagnostics->correction_iterations;
            diagnostics->correction_retries++;
            diagnostics->effective_mean_stress =
                md_trace(candidate_stress)/3.0;
            if (dT <= dT_min) {
                *yield_error = correction_residual;
                if (force_accept) {
                    double forced_residual = 0.0;
                    if (md_force_accept_candidate(
                            m, force_candidate_stress,
                            force_candidate_pstrain,
                            force_candidate_fabric, &stress, &pstrain,
                            &alpha_work, &alpha_in_work, &fabric_work,
                            K, G, &forced_residual)) {
                        diagnostics->force_accepts++;
                        diagnostics->correction_residual = forced_residual;
                        *plastic_multiplier = candidate_multiplier;
                        *integration_error = fmax(
                            *integration_error, step_error);
                        *yield_error = forced_residual;
                        (*accepted_substeps)++;
                        T += dT;
                        if (T >= 1.0-DBL_EPSILON) break;
                        dT = fmin(
                            fmax(2.0*dT, dT_min), 1.0-T);
                        continue;
                    }
                }
                diagnostics->failure_reason =
                    NL_MD_FAILURE_STRESS_CORRECTION_MIN_STEP;
                return 0;
            }
            dT = fmax(0.5*dT, dT_min);
            continue;
        }
        diagnostics->correction_iterations =
            correction_iterations > diagnostics->correction_iterations
                ? correction_iterations : diagnostics->correction_iterations;
        diagnostics->correction_residual = correction_residual;
        diagnostics->projection_fallbacks += projection_used;

        stress = candidate_stress;
        alpha_work = candidate_alpha;
        alpha_in_work = candidate_alpha_in;
        fabric_work = candidate_fabric;
        pstrain = candidate_pstrain;
        *plastic_multiplier = candidate_multiplier;
        *integration_error = fmax(*integration_error, step_error);
        (*accepted_substeps)++;
        *yield_error = fmax(*yield_error, correction_residual);
        T += dT;
        if (T >= 1.0-DBL_EPSILON) break;
        const double q = (step_error <= DBL_MIN) ? 2.0 : fmax(
            0.8*sqrt(cfg->backbone_error_tolerance/step_error), 0.5);
        dT = fmin(fmax(q*dT, dT_min), 1.0-T);
    }

    const tensor_t dynamic_next = md_sub(md_scale(stress, -1.0), static_stress);
    const tensor_t plastic_next = md_scale(pstrain, -1.0);
    const double final_f = md_yield(m, stress, alpha_work);
    if (!isfinite(final_f)) {
        *yield_error = DBL_MAX;
        diagnostics->failure_reason = NL_MD_FAILURE_FINAL_STATE;
        return 0;
    }
    /* error_tolerance is absolute only for a plastic stress correction. The
     * final state can instead be an elastic unloading state strictly inside
     * the surface (f < 0), which remains admissible by the same one-sided
     * criterion used by the elastic predictor. Keep prior absolute plastic
     * correction residuals, but count only a positive final violation. */
    *yield_error = fmax(*yield_error, fmax(final_f, 0.0));
    if (!md_finite(dynamic_next) || !md_finite(plastic_next) ||
        !md_finite(alpha_work) || !md_finite(alpha_in_work) ||
        !md_finite(fabric_work) || !isfinite(*yield_error) ||
        !(*yield_error < cfg->error_tolerance)) {
        diagnostics->failure_reason = NL_MD_FAILURE_FINAL_STATE;
        return 0;
    }
    *dynamic_stress = dynamic_next;
    *plastic_strain = plastic_next;
    *alpha = alpha_work;
    *alpha_in = alpha_in_work;
    *fabric = fabric_work;
    *initialized = initialized_work;
    diagnostics->effective_mean_stress = md_trace(stress)/3.0;
    return 1;
}

/* Runtime-dispatching wrapper. The production kernels call the specialization
 * directly so each one carries a single integrator; this entry point exists
 * for the host material-point driver, which must honor whichever scheme the
 * fixture selected and is not performance critical. */
NL_HD static int nl_update_manzari_dafalias(
    const nlconstants_t *__restrict__ c,
    const nonlinear_gpu_config_t *__restrict__ cfg,
    tensor_t strain, tensor_t strain_prev, tensor_t static_stress,
    tensor_t *__restrict__ dynamic_stress,
    tensor_t *__restrict__ plastic_strain,
    tensor_t *__restrict__ alpha, tensor_t *__restrict__ alpha_in,
    tensor_t *__restrict__ fabric, double *__restrict__ initialized,
    double *__restrict__ void_ratio,
    double *__restrict__ plastic_multiplier,
    double *__restrict__ yield_error,
    double *__restrict__ integration_error,
    int32_t *__restrict__ accepted_substeps,
    nl_md_diagnostics_t *__restrict__ diagnostics)
{
    if (cfg->md_integration_scheme == NL_MD_INTEGRATION_FORWARD_EULER) {
        if (cfg->md_forward_euler_substeps ==
                NL_MD_FORWARD_EULER_SUBSTEPS_DEMAND) {
            return nl_update_manzari_dafalias_scheme<
                NL_MD_INTEGRATION_FORWARD_EULER,
                NL_MD_FORWARD_EULER_SUBSTEPS_DEMAND>(
                c, cfg, strain, strain_prev, static_stress, dynamic_stress,
                plastic_strain, alpha, alpha_in, fabric, initialized,
                void_ratio, plastic_multiplier, yield_error,
                integration_error, accepted_substeps, diagnostics);
        }
        return nl_update_manzari_dafalias_scheme<
            NL_MD_INTEGRATION_FORWARD_EULER,
            NL_MD_FORWARD_EULER_SUBSTEPS_FIXED>(
            c, cfg, strain, strain_prev, static_stress, dynamic_stress,
            plastic_strain, alpha, alpha_in, fabric, initialized, void_ratio,
            plastic_multiplier, yield_error, integration_error,
            accepted_substeps, diagnostics);
    }
    return nl_update_manzari_dafalias_scheme<
        NL_MD_INTEGRATION_MODIFIED_EULER,
        NL_MD_FORWARD_EULER_SUBSTEPS_FIXED>(
        c, cfg, strain, strain_prev, static_stress, dynamic_stress,
        plastic_strain, alpha, alpha_in, fabric, initialized, void_ratio,
        plastic_multiplier, yield_error, integration_error,
        accepted_substeps, diagnostics);
}

#undef NL_HD

#ifndef NL_GPU_MATERIAL_POINT_ONLY
__device__ static inline void nl_record_error(nonlinear_gpu_error_t *error,
                                               int32_t code,
                                               int32_t eindex, int32_t qp,
                                               double residual, double value)
{
    if (atomicCAS(&error->code, 0, code) == 0) {
        error->eindex = eindex;
        error->qp = qp;
        error->residual = residual;
        error->value = value;
    }
}

__device__ static inline void nl_record_md_error(
    nonlinear_gpu_error_t *error, int32_t eindex, int32_t qp,
    double residual, double integration_error,
    int32_t accepted_substeps, const nl_md_diagnostics_t *diagnostics,
    tensor_t stress,
    tensor_t strain_increment)
{
    if (atomicCAS(&error->code, 0, 3) == 0) {
        error->eindex = eindex;
        error->qp = qp;
        error->md_reason = diagnostics->failure_reason;
        error->md_attempts = diagnostics->attempts;
        error->md_accepted_substeps = accepted_substeps;
        error->md_correction_retries = diagnostics->correction_retries;
        error->md_correction_iterations = diagnostics->correction_iterations;
        error->md_projection_fallbacks = diagnostics->projection_fallbacks;
        error->md_force_accepts = diagnostics->force_accepts;
        error->residual = residual;
        error->value = integration_error;
        error->md_substep_fraction = diagnostics->substep_fraction;
        error->md_progress_fraction = diagnostics->progress_fraction;
        error->md_effective_mean_stress =
            diagnostics->effective_mean_stress;
        error->md_strain_increment_norm = md_norm(strain_increment);
        error->md_first_dgamma = diagnostics->first_dgamma;
        error->md_second_dgamma = diagnostics->second_dgamma;
        error->md_first_raw_dgamma = diagnostics->first_raw_dgamma;
        error->md_second_raw_dgamma = diagnostics->second_raw_dgamma;
        error->md_first_denominator = diagnostics->first_denominator;
        error->md_second_denominator = diagnostics->second_denominator;
        error->md_first_h = diagnostics->first_h;
        error->md_second_h = diagnostics->second_h;
        error->md_first_hardening_gap = diagnostics->first_hardening_gap;
        error->md_second_hardening_gap = diagnostics->second_hardening_gap;
        error->md_stress = stress;
        error->md_strain_increment = strain_increment;
    }
}

__device__ static inline void nl_record_pdmy03_error(
    nonlinear_gpu_error_t *error, int32_t eindex, int32_t qp,
    const pdmy03_diagnostics_t *diagnostics, tensor_t stress,
    tensor_t strain_increment)
{
    if (atomicCAS(&error->code, 0, 4) == 0) {
        error->eindex = eindex;
        error->qp = qp;
        error->md_reason = diagnostics->reason;
        error->md_attempts = diagnostics->requested_substeps;
        error->md_accepted_substeps = diagnostics->accepted_substeps;
        error->md_projection_fallbacks = diagnostics->crossed_surfaces;
        error->md_force_accepts = diagnostics->substeps_clamped;
        error->residual = diagnostics->yield_residual;
        error->value = diagnostics->mean_stress;
        error->md_effective_mean_stress = -diagnostics->mean_stress;
        error->md_strain_increment_norm = diagnostics->strain_increment_norm;
        error->md_stress = stress;
        error->md_strain_increment = strain_increment;
    }
}

__device__ static inline void nl_point_dxi(double *dx, double *dy, double *dz,
                                            double lx, double ly, double lz,
                                            double h, int node)
{
    const int sx = (node & 1) ? 1 : -1;
    const int sy = (node & 2) ? 1 : -1;
    const int sz = (node & 4) ? 1 : -1;
    const double jac = 0.25 / h;
    *dx = jac * sx * (1.0 + sy * ly) * (1.0 + sz * lz);
    *dy = jac * (1.0 + sx * lx) * sy * (1.0 + sz * lz);
    *dz = jac * (1.0 + sx * lx) * (1.0 + sy * ly) * sz;
}

__device__ static tensor_t nl_hex_strain(const fvector_t u[8], int qp,
                                          double h)
{
    const double q = 0.577350269189;
    const double lx = (qp & 1) ? q : -q;
    const double ly = (qp & 2) ? q : -q;
    const double lz = (qp & 4) ? q : -q;
    tensor_t strain = nl_zero();
    for (int node = 0; node < 8; node++) {
        double dx, dy, dz;
        nl_point_dxi(&dx, &dy, &dz, lx, ly, lz, h, node);
        strain.xx += dx * (double)u[node].f[0];
        strain.yy += dy * (double)u[node].f[1];
        strain.zz += dz * (double)u[node].f[2];
        strain.xy += 0.5 * (dy * (double)u[node].f[0] +
                            dx * (double)u[node].f[1]);
        strain.yz += 0.5 * (dz * (double)u[node].f[1] +
                            dy * (double)u[node].f[2]);
        strain.xz += 0.5 * (dz * (double)u[node].f[0] +
                            dx * (double)u[node].f[2]);
    }
    return strain;
}

__device__ static inline void nl_tetra_nodes(int part, int tet, int n[4])
{
    const int part1[5][4] = {{0,2,1,4}, {3,1,2,7}, {6,4,7,2},
                             {5,7,4,1}, {2,4,7,1}};
    const int part2[5][4] = {{0,3,1,5}, {0,2,3,6}, {4,5,6,0},
                             {6,5,7,3}, {0,6,3,5}};
    for (int i = 0; i < 4; i++) n[i] = (part == 1) ? part1[tet][i]
                                                    : part2[tet][i];
}

__device__ static inline double nl_det3(const double a[3][3])
{
    return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
           a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
           a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
}

__device__ static int nl_tetra_gradients(int part, int tet, double h,
                                         int nodes[4], double grad[4][3],
                                         double *volume)
{
    nl_tetra_nodes(part, tet, nodes);
    double xyz[4][3];
    for (int a = 0; a < 4; a++) {
        xyz[a][0] = (nodes[a] & 1) ? h : 0.0;
        xyz[a][1] = (nodes[a] & 2) ? h : 0.0;
        xyz[a][2] = (nodes[a] & 4) ? h : 0.0;
    }
    double j[3][3];
    for (int r = 0; r < 3; r++) {
        j[r][0] = xyz[1][r] - xyz[0][r];
        j[r][1] = xyz[2][r] - xyz[0][r];
        j[r][2] = xyz[3][r] - xyz[0][r];
    }
    const double det = nl_det3(j);
    if (fabs(det) < DBL_MIN) return 0;
    const double invdet = 1.0 / det;
    double inv[3][3];
    inv[0][0] =  (j[1][1]*j[2][2]-j[1][2]*j[2][1])*invdet;
    inv[0][1] = -(j[0][1]*j[2][2]-j[0][2]*j[2][1])*invdet;
    inv[0][2] =  (j[0][1]*j[1][2]-j[0][2]*j[1][1])*invdet;
    inv[1][0] = -(j[1][0]*j[2][2]-j[1][2]*j[2][0])*invdet;
    inv[1][1] =  (j[0][0]*j[2][2]-j[0][2]*j[2][0])*invdet;
    inv[1][2] = -(j[0][0]*j[1][2]-j[0][2]*j[1][0])*invdet;
    inv[2][0] =  (j[1][0]*j[2][1]-j[1][1]*j[2][0])*invdet;
    inv[2][1] = -(j[0][0]*j[2][1]-j[0][1]*j[2][0])*invdet;
    inv[2][2] =  (j[0][0]*j[1][1]-j[0][1]*j[1][0])*invdet;

    /* grad(N1..N3) are rows of J^-1; grad(N0) is their negative sum. */
    for (int d = 0; d < 3; d++) {
        grad[1][d] = inv[0][d];
        grad[2][d] = inv[1][d];
        grad[3][d] = inv[2][d];
        grad[0][d] = -grad[1][d] - grad[2][d] - grad[3][d];
    }
    *volume = fabs(det) / 6.0;
    return 1;
}

__device__ static tensor_t nl_tetra_strain(const fvector_t u[8], int part,
                                            int tet, double h)
{
    int nodes[4];
    double grad[4][3], volume;
    tensor_t strain = nl_zero();
    if (!nl_tetra_gradients(part, tet, h, nodes, grad, &volume)) return strain;
    for (int a = 0; a < 4; a++) {
        const fvector_t v = u[nodes[a]];
        strain.xx += grad[a][0] * (double)v.f[0];
        strain.yy += grad[a][1] * (double)v.f[1];
        strain.zz += grad[a][2] * (double)v.f[2];
        strain.xy += 0.5 * (grad[a][1] * (double)v.f[0] +
                            grad[a][0] * (double)v.f[1]);
        strain.yz += 0.5 * (grad[a][2] * (double)v.f[1] +
                            grad[a][1] * (double)v.f[2]);
        strain.xz += 0.5 * (grad[a][2] * (double)v.f[0] +
                            grad[a][0] * (double)v.f[2]);
    }
    return strain;
}
#endif

__device__ static inline double nl_gqh_value(const nlconstants_t *c,
                                              double gamma)
{
    const double t1 = c->thetaGQH[0], t2 = c->thetaGQH[1];
    const double t3 = c->thetaGQH[2], t4 = c->thetaGQH[3];
    const double t5 = c->thetaGQH[4];
    const double gp = pow(nl_max(gamma, 0.0), t5);
    const double theta = t1 + t2 * (t4 * gp) / (pow(t3, t5) + t4 * gp);
    const double rad = nl_max((1.0 + gamma) * (1.0 + gamma) -
                              4.0 * theta * gamma, 0.0);
    return 2.0 * gamma / (1.0 + gamma + sqrt(rad));
}

__device__ static inline double nl_mkz_value(const nlconstants_t *c,
                                              double gamma)
{
    return gamma / (1.0 + c->beta_MKZ * pow(nl_max(gamma, 0.0),
                                             c->s_MKZ));
}

__device__ static inline double nl_gqh_derivative(const nlconstants_t *c,
                                                   double gamma)
{
    const double t1 = c->thetaGQH[0], t2 = c->thetaGQH[1];
    const double t3 = c->thetaGQH[2], t4 = c->thetaGQH[3];
    const double t5 = c->thetaGQH[4];
    const double gp = pow(nl_max(gamma, DBL_MIN), t5);
    const double t3p = pow(t3, t5);
    const double denominator = t3p + t4 * gp;
    const double theta = t1 + t2 * t4 * gp / denominator;
    const double dtheta = t2 * t3p * t4 * t5 *
        pow(nl_max(gamma, DBL_MIN), t5 - 1.0) /
        (denominator * denominator);
    const double root = sqrt(nl_max(
        (1.0 + gamma) * (1.0 + gamma) - 4.0 * theta * gamma,
        DBL_MIN));
    const double a = 1.0 + gamma + root;
    const double b = 1.0 + (1.0 + gamma -
        2.0 * (theta + gamma * dtheta)) / root;
    return 2.0 * (a - gamma * b) / (a * a);
}

__device__ static inline double nl_mkz_derivative(const nlconstants_t *c,
                                                   double gamma)
{
    const double p = pow(nl_max(gamma, DBL_MIN), c->s_MKZ);
    return (1.0 + c->beta_MKZ * (1.0 - c->s_MKZ) * p) /
           ((1.0 + c->beta_MKZ * p) *
            (1.0 + c->beta_MKZ * p));
}

__device__ static double nl_backbone_gamma(const nlconstants_t *c,
                                            materialmodel_t model,
                                            double target, double tolerance)
{
    double lo = 0.0, hi = 1.0;
    double flo = (model == VONMISES_GQH) ? nl_gqh_value(c, lo)
                                          : nl_mkz_value(c, lo);
    double fhi = (model == VONMISES_GQH) ? nl_gqh_value(c, hi)
                                          : nl_mkz_value(c, hi);
    for (int i = 0; i < 200 && fhi < target; i++) {
        hi *= 2.0;
        fhi = (model == VONMISES_GQH) ? nl_gqh_value(c, hi)
                                       : nl_mkz_value(c, hi);
    }

    /* This inversion is nested inside the Borja--Amies stress-residual
     * solve. Stopping at the user-facing dimensionless backbone tolerance
     * quantizes H(kappa), so the outer residual can stall well above its
     * Pa-valued tolerance. Refine to roundoff here; this remains stricter
     * than the configured backbone tolerance and makes the outer function
     * continuous at the accuracy it requires. */
    const double solve_tolerance = nl_min(
        tolerance, DBL_EPSILON * nl_max(1.0, fabs(target)));
    double gamma = 0.5 * (lo + hi);
    for (int i = 0; i < 200; i++) {
        const double value = (model == VONMISES_GQH)
                                 ? nl_gqh_value(c, gamma)
                                 : nl_mkz_value(c, gamma);
        const double residual = value - target;
        if (fabs(value - target) <= solve_tolerance) return gamma;
        if (residual < 0.0) {
            lo = gamma;
            flo = value;
        } else {
            hi = gamma;
            fhi = value;
        }

        const double derivative = (model == VONMISES_GQH)
                                      ? nl_gqh_derivative(c, gamma)
                                      : nl_mkz_derivative(c, gamma);
        double candidate = gamma - residual / derivative;
        if (!isfinite(candidate) || derivative <= 0.0 || candidate <= lo ||
            candidate >= hi) {
            candidate = 0.5 * (lo + hi);
        }
        if (candidate == gamma)
            return (fabs(flo - target) <= fabs(fhi - target)) ? lo : hi;
        gamma = candidate;
    }
    return (fabs(flo - target) <= fabs(fhi - target)) ? lo : hi;
}

__device__ static double nl_hardening(const nlconstants_t *c,
                                      materialmodel_t model, double kappa,
                                      double gamma_n, double tolerance)
{
    const double g = c->mu;
    if (kappa <= 0.0 && model != VONMISES_MKZ) return 0.0;
    if (model == VONMISES_BAE)
        return c->psi0 * g * pow(nl_max(kappa, 0.0), c->m);
    if (model == VONMISES_BAH)
        return 3.0 * g * kappa * kappa / (1.0 + 2.0 * kappa);
    if (model == VONMISES_RO)
        return 3.0 * g / (c->alpha_RO * c->eta_RO) *
               pow(c->phi_RO * (1.0 + kappa), c->eta_RO - 1.0);

    if (model == VONMISES_GQH || model == VONMISES_MKZ) {
        double target = 1.0 / (1.0 + kappa);
        if (model == VONMISES_MKZ) target /= c->phi_MKZ;
        const double gamma = nl_backbone_gamma(c, model, target, tolerance);
        double derivative = (model == VONMISES_MKZ)
                                ? nl_mkz_derivative(c, gamma)
                                : nl_gqh_derivative(c, gamma);
        derivative = nl_min(derivative, 1.0 - 1.0e-14);
        return 3.0 * g * derivative / nl_max(1.0 - derivative, DBL_MIN);
    }
    (void)gamma_n;
    return 0.0;
}

typedef struct nl_ba_params_t {
    double a, b, c, d, e, f;
    double gamma_n;
    nlconstants_t constants;
    materialmodel_t model;
    double backbone_tolerance;
} nl_ba_params_t;

__device__ static double nl_ba_residual(double kappa,
                                         const nl_ba_params_t *p)
{
    const double h = nl_hardening(&p->constants, p->model, kappa,
                                  p->gamma_n, p->backbone_tolerance);
    const double mu = p->constants.mu;
    const double psi = (h + 3.0 * mu != 0.0)
                           ? 2.0 * mu * h / (h + 3.0 * mu) : 0.0;
    const double phi = psi * (1.0 + kappa);
    const double rad = p->a + 2.0 * p->b * kappa +
                       kappa * kappa * p->c + 2.0 * phi * p->d +
                       2.0 * kappa * phi * p->e + phi * phi * p->f;
    const double r = p->constants.c * sqrt(8.0 / 3.0);
    return sqrt(nl_max(rad, 0.0)) - r;
}

__device__ static int nl_bracketed_bisection(const nl_ba_params_t *p,
                                              double *lo, double *hi,
                                              double tolerance,
                                              double *root,
                                              double *residual)
{
    double flo = nl_ba_residual(*lo, p);
    double fhi = nl_ba_residual(*hi, p);
    if (!isfinite(flo) || !isfinite(fhi) || flo * fhi > 0.0) return 0;
    for (int i = 0; i < 1000; i++) {
        const double mid = 0.5 * (*lo + *hi);
        const double fm = nl_ba_residual(mid, p);
        if (!isfinite(fm)) return 0;
        if (fabs(fm) <= tolerance) {
            *root = mid;
            *residual = fm;
            return 1;
        }
        if (fabs(*hi - *lo) <= DBL_EPSILON * nl_max(1.0, fabs(mid))) {
            *root = mid;
            *residual = fm;
            return 0;
        }
        if (flo * fm <= 0.0) {
            *hi = mid;
            fhi = fm;
        } else {
            *lo = mid;
            flo = fm;
        }
    }
    *root = 0.5 * (*lo + *hi);
    *residual = nl_ba_residual(*root, p);
    return fabs(*residual) <= tolerance;
}

__device__ static void nl_von_mises_boundary(tensor_t sigma_n,
                                              tensor_t de_dev,
                                              tensor_t *sigma,
                                              double g, double bulk,
                                              double de_vol, double radius)
{
    tensor_t s_trial = nl_add(nl_dev(sigma_n), nl_scale(de_dev, 2.0 * g));
    const double mag = sqrt(nl_max(nl_ddot(s_trial, s_trial), DBL_MIN));
    tensor_t p = nl_iso(nl_i1(sigma_n) / 3.0 + bulk * de_vol);
    *sigma = nl_add(p, nl_scale(s_trial, radius / mag));
}

__device__ static void nl_backbone_diagnostics(const nlconstants_t *c,
                                                materialmodel_t model,
                                                double kappa,
                                                double gamma_n,
                                                double tolerance,
                                                double *gamma, double *tao,
                                                double *ggmax)
{
    const double tao_max = 2.0 * c->c / sqrt(3.0);
    const double tao_bar = 1.0 / (1.0 + kappa);
    if (model == VONMISES_MKZ) {
        const double target = tao_bar / c->phi_MKZ;
        const double gb = nl_backbone_gamma(c, model, target, tolerance);
        *gamma = gb * c->gammaref_MKZ;
        *tao = tao_bar * tao_max;
        *ggmax = 1.0 / (1.0 + c->beta_MKZ * pow(gb, c->s_MKZ));
    } else if (model == VONMISES_BAH) {
        const double gamma_ref = tao_max / c->mu;
        const double gb = (tao_bar == 1.0) ? 0.10 / gamma_ref
                                           : tao_bar / (1.0 - tao_bar);
        *gamma = gb * gamma_ref;
        *tao = tao_bar * tao_max;
        *ggmax = 1.0 / (1.0 + gb);
    } else if (model == VONMISES_BAE) {
        const double points[10] = {
             0.1488743389816312108848260,  0.4333953941292471907992659,
             0.6794095682990244062343274,  0.8650633666889845107320967,
             0.9739065285171717200779640, -0.1488743389816312108848260,
            -0.4333953941292471907992659, -0.6794095682990244062343274,
            -0.8650633666889845107320967, -0.9739065285171717200779640};
        const double weights[10] = {
            0.2955242247147528701738930, 0.2692667193099963550912269,
            0.2190863625159820439955349, 0.1494513491505805931457763,
            0.0666713443086881375935688, 0.2955242247147528701738930,
            0.2692667193099963550912269, 0.2190863625159820439955349,
            0.1494513491505805931457763, 0.0666713443086881375935688};
        const double h = c->psi0*c->mu;
        const double tau = tao_bar*tao_max;
        double gamma_value = tau/c->mu;
        double reduction = 0.0;
        for (int i = 0; i < 10; i++) {
            const double xi = 0.5*tau*(points[i] + 1.0);
            gamma_value += 3.0/h*pow(xi/(tao_max-xi), c->m)*
                           weights[i]*tau*0.5;
            reduction += pow(2.0*xi/(tao_max-2.0*xi+tau), c->m)*
                         weights[i]*tau;
        }
        *gamma = gamma_value;
        *tao = tau;
        *ggmax = 1.0/(1.0 + 3.0*c->mu*reduction/(2.0*h*tau));
    }
    /* The CPU reference currently leaves the GQH and Ramberg--Osgood
     * station-only backbone diagnostic values unchanged. Preserve that
     * behavior while still updating their constitutive stress and kappa. */
    (void)gamma_n;
    (void)tolerance;
}

__device__ static int nl_update_backbone(const nlconstants_t *c,
                                          materialmodel_t model,
                                          tensor_t strain,
                                          tensor_t strain_prev,
                                          tensor_t *sref,
                                          tensor_t *sigma,
                                          double *kappa,
                                          double *gamma1d,
                                          double *tao1d,
                                          double *ggmax,
                                          const nonlinear_gpu_config_t *cfg,
                                          double *error)
{
    const double g = c->mu;
    const double bulk = c->lambda + 2.0 * g / 3.0;
    const tensor_t sigma_n = *sigma;
    const double kappa_n = *kappa;
    const double gamma_n = *gamma1d;
    tensor_t sdev_n = nl_dev(sigma_n);
    tensor_t de = nl_sub(strain, strain_prev);
    const double de_vol = nl_i1(de);
    tensor_t de_dev = nl_dev(de);

    /* The model is incremental: a true zero increment preserves the previous
     * stress and every backbone/history variable without a reversal test or
     * scalar root solve. */
    if (de.xx == 0.0 && de.yy == 0.0 && de.zz == 0.0 &&
        de.xy == 0.0 && de.yz == 0.0 && de.xz == 0.0) {
        *error = 0.0;
        return 1;
    }

    /* Pressure-insensitive von Mises backbones have an analytic volumetric
     * update. Do not send a zero or roundoff-only deviatoric increment to a
     * scalar equation that has no finite kappa root. */
    const tensor_t dsdev_elastic = nl_scale(de_dev, 2.0 * g);
    const double dsdev_measure =
        sqrt(nl_max(nl_j2(dsdev_elastic), 0.0));
    double stress_scale = nl_max(1.0, fabs(bulk * de_vol));
    stress_scale = nl_max(stress_scale,
                          sqrt(nl_max(nl_j2(sdev_n), 0.0)));
    if (dsdev_measure <= 64.0 * DBL_EPSILON * stress_scale) {
        *sigma = nl_add(sigma_n, nl_iso(bulk * de_vol));
        *error = 0.0;
        return 1;
    }

    /* The kappa=100 sentinel denotes virgin state. Near zero deviatoric
     * stress, use the elastic asymptote and accumulate it incrementally until
     * the boundary root becomes numerically resolvable. This is deliberately
     * not applied to established cyclic history. */
    const tensor_t sdev_trial = nl_add(sdev_n, dsdev_elastic);
    if (kappa_n == 100.0 &&
        sqrt(nl_max(nl_j2(sdev_trial), 0.0)) <
            NONLINEAR_BACKBONE_STRESS_FLOOR_PA) {
        *sigma = nl_add(nl_add(sigma_n, nl_iso(bulk * de_vol)),
                        dsdev_elastic);
        *error = 0.0;
        return 1;
    }

    tensor_t sno = nl_sub(sdev_n, *sref);
    const double den = nl_ddot(sdev_n, sno) + kappa_n * nl_ddot(sno, sno);
    const tensor_t num = nl_add(nl_scale(sdev_n, 1.0 + kappa_n),
                                nl_scale(sno, kappa_n * (1.0 + kappa_n)));
    const double load_unload = (fabs(den) > DBL_MIN)
                                   ? -nl_ddot(num, de_dev) / den : 0.0;
    if (load_unload > 0.0) *sref = sdev_n;

    sno = nl_sub(sdev_n, *sref);
    nl_ba_params_t p;
    p.a = nl_ddot(sdev_n, sdev_n);
    p.b = nl_ddot(sdev_n, sno);
    p.c = nl_ddot(sno, sno);
    p.d = nl_ddot(sdev_n, de_dev);
    p.e = nl_ddot(sno, de_dev);
    p.f = nl_ddot(de_dev, de_dev);
    p.gamma_n = gamma_n;
    p.constants = *c;
    p.model = model;
    p.backbone_tolerance = cfg->backbone_error_tolerance;

    double f0 = nl_ba_residual(0.0, &p);
    if (fabs(f0) <= cfg->error_tolerance) {
        *kappa = 0.0;
        *error = f0;
        nl_von_mises_boundary(sigma_n, de_dev, sigma, g, bulk, de_vol,
                              c->c * sqrt(8.0 / 3.0));
        nl_backbone_diagnostics(c, model, *kappa, gamma_n,
                                cfg->backbone_error_tolerance,
                                gamma1d, tao1d, ggmax);
        return 1;
    }

    double lo, hi;
    if (kappa_n == 100.0 || load_unload > 0.0) {
        lo = 0.0;
        hi = NL_GPU_KAPPA_MIN;
    } else {
        hi = nl_max(kappa_n, NL_GPU_KAPPA_MIN);
        lo = 0.5 * hi;
    }
    double flo = nl_ba_residual(lo, &p);
    double fhi = nl_ba_residual(hi, &p);
    for (int i = 0; i < 1000 && flo * fhi > 0.0; i++) {
        if (kappa_n != 100.0 && load_unload <= 0.0 && lo > NL_GPU_KAPPA_MIN) {
            hi = lo;
            fhi = flo;
            lo *= 0.5;
            flo = nl_ba_residual(lo, &p);
        } else {
            lo = hi;
            flo = fhi;
            hi *= 2.0;
            fhi = nl_ba_residual(hi, &p);
        }
    }
    if (flo * fhi > 0.0) {
        lo = 0.0;
        hi = NL_GPU_KAPPA_MIN;
        flo = nl_ba_residual(lo, &p);
        fhi = nl_ba_residual(hi, &p);
        for (int i = 0; i < 1000 && flo * fhi > 0.0; i++) {
            hi *= 2.0;
            fhi = nl_ba_residual(hi, &p);
        }
    }
    if (hi <= NL_GPU_KAPPA_MIN) {
        *kappa = 0.0;
        *error = f0;
        nl_von_mises_boundary(sigma_n, de_dev, sigma, g, bulk, de_vol,
                              c->c * sqrt(8.0 / 3.0));
        nl_backbone_diagnostics(c, model, *kappa, gamma_n,
                                cfg->backbone_error_tolerance,
                                gamma1d, tao1d, ggmax);
        return 1;
    }

    double root = 0.0, residual = DBL_MAX;
    if (!nl_bracketed_bisection(&p, &lo, &hi, cfg->error_tolerance,
                                &root, &residual)) {
        *error = residual;
        return 0;
    }
    *kappa = root;
    *error = fabs(residual);
    const double hmod = nl_hardening(c, model, root, gamma_n,
                                     cfg->backbone_error_tolerance);
    const double psi = 2.0 * g * hmod / nl_max(hmod + 3.0 * g, DBL_MIN);
    *sigma = nl_add(nl_add(sigma_n, nl_iso(bulk * de_vol)),
                    nl_scale(de_dev, psi));
    nl_backbone_diagnostics(c, model, root, gamma_n,
                            cfg->backbone_error_tolerance,
                            gamma1d, tao1d, ggmax);
    return 1;
}

__device__ static double nl_positive_quartic_root(double c1, double c2,
                                                   double c3, double c4,
                                                   double c5,
                                                   double tolerance,
                                                   int *ok)
{
    typedef struct nl_complex_t { double real, imag; } nl_complex_t;
    if (fabs(c1) > DBL_MIN) {
        const double a = c2/c1, b = c3/c1, c = c4/c1, d = c5/c1;
        const double radius = 1.0+nl_max(nl_max(fabs(a), fabs(b)),
                                         nl_max(fabs(c), fabs(d)));
        nl_complex_t roots[4] = {{radius, 0.0},
                                 {0.4*radius, 0.9*radius},
                                 {-0.65*radius, 0.72*radius},
                                 {-0.7*radius, -0.8*radius}};
        for (int iteration = 0; iteration < 160; iteration++) {
            double maximum_update = 0.0;
            for (int i = 0; i < 4; i++) {
                const nl_complex_t z = roots[i];
                nl_complex_t value = {1.0, 0.0};
                const double coefficients[4] = {a, b, c, d};
                for (int k = 0; k < 4; k++) {
                    const nl_complex_t next = {
                        value.real*z.real-value.imag*z.imag+coefficients[k],
                        value.real*z.imag+value.imag*z.real};
                    value = next;
                }
                nl_complex_t product = {1.0, 0.0};
                for (int j = 0; j < 4; j++) {
                    if (j == i) continue;
                    const nl_complex_t difference = {
                        z.real-roots[j].real, z.imag-roots[j].imag};
                    const nl_complex_t next = {
                        product.real*difference.real-
                            product.imag*difference.imag,
                        product.real*difference.imag+
                            product.imag*difference.real};
                    product = next;
                }
                const double denominator = product.real*product.real+
                                           product.imag*product.imag;
                if (denominator <= DBL_MIN) continue;
                const nl_complex_t update = {
                    (value.real*product.real+value.imag*product.imag)/denominator,
                    (value.imag*product.real-value.real*product.imag)/denominator};
                roots[i].real -= update.real;
                roots[i].imag -= update.imag;
                maximum_update = nl_max(maximum_update,
                    sqrt(update.real*update.real+update.imag*update.imag));
            }
            if (maximum_update <= DBL_EPSILON*radius*16.0) break;
        }
        double selected = DBL_MAX;
        for (int i = 0; i < 4; i++) {
            if (roots[i].real >= 0.0 &&
                fabs(roots[i].imag) <= 1.0e-9*
                    nl_max(1.0, fabs(roots[i].real))) {
                selected = nl_min(selected, roots[i].real);
            }
        }
        if (selected != DBL_MAX) {
            const double residual = ((((c1*selected+c2)*selected+c3)*
                                       selected+c4)*selected+c5);
            const double scale = fabs(c5)+fabs(c4*selected)+
                fabs(c3*selected*selected)+
                fabs(c2*selected*selected*selected)+
                fabs(c1*selected*selected*selected*selected)+1.0;
            if (isfinite(residual) && fabs(residual) <=
                nl_max(tolerance, 1.0e-10)*scale) {
                *ok = 1;
                return selected;
            }
        }
    }

    /* Degenerate or poorly conditioned fallback: locate a sign-changing
     * interval on the non-negative real axis and bisect it. */
    double lo = 0.0;
    double flo = c5;
    if (fabs(flo) <= tolerance) {
        *ok = 1;
        return 0.0;
    }
    double hi = 1.0e-16;
    double fhi = ((((c1 * hi + c2) * hi + c3) * hi + c4) * hi + c5);
    for (int i = 0; i < 256 && flo * fhi > 0.0; i++) {
        lo = hi;
        flo = fhi;
        hi *= 2.0;
        fhi = ((((c1 * hi + c2) * hi + c3) * hi + c4) * hi + c5);
    }
    if (!isfinite(fhi) || flo * fhi > 0.0) {
        *ok = 0;
        return 0.0;
    }
    for (int i = 0; i < 256; i++) {
        const double mid = 0.5 * (lo + hi);
        const double fm = ((((c1 * mid + c2) * mid + c3) * mid + c4) * mid + c5);
        const double scale = fabs(c5) + fabs(c4 * mid) + fabs(c3 * mid * mid) +
                             fabs(c2 * mid * mid * mid) +
                             fabs(c1 * mid * mid * mid * mid) + 1.0;
        if (fabs(fm) <= tolerance * scale ||
            fabs(hi - lo) <= DBL_EPSILON * nl_max(1.0, mid)) {
            *ok = 1;
            return mid;
        }
        if (flo * fm <= 0.0) {
            hi = mid;
            fhi = fm;
        } else {
            lo = mid;
            flo = fm;
        }
    }
    *ok = 1;
    return 0.5 * (lo + hi);
}

__device__ static int nl_update_fa(const nlconstants_t *c,
                                    materialmodel_t model,
                                    tensor_t dev_predictor,
                                    tensor_t strain_prev,
                                    tensor_t pstrain_prev,
                                    tensor_t backstress_prev,
                                    double epbar_prev,
                                    tensor_t elastic_stress,
                                    tensor_t *pstrain,
                                    tensor_t *backstress,
                                    tensor_t *sigma,
                                    double *epbar,
                                    double *fs,
                                    double *psi_state,
                                    double *load_unload_state,
                                    double *tao_state,
                                    double *tao_max_state,
                                    double tolerance)
{
    double psi = c->psi0;
    double su = c->c;
    double sy = c->Sstrain0 * c->mu * sqrt(2.0);
    if (model == VONMISES_FAM && *psi_state == 0.0) *psi_state = psi;

    const double s_ss = 2.0 * nl_j2(dev_predictor);
    const double s_aa = 2.0 * nl_j2(backstress_prev);
    const double s_sa = 2.0 * nl_comb_j2(backstress_prev, dev_predictor);
    double hkin = ((model == VONMISES_FAM) ? *psi_state : psi) * c->mu;
    double hnlin = hkin / nl_max(sqrt(8.0 / 3.0) * su - sy, DBL_MIN);
    double g1 = c->mu + 0.5 * hkin;

    double c1 = (2.0 * c->mu * hnlin) * (2.0 * c->mu * hnlin);
    double c2 = (4.0 * sy * c->mu * hnlin + 8.0 * c->mu * g1) * hnlin;
    double c3 = hnlin * hnlin * (sy * sy - s_ss) + 4.0 * g1 * g1 +
                4.0 * hnlin * sy * (c->mu + g1);
    double c4 = 2.0 * hnlin * (sy * sy + s_sa - s_ss) + 4.0 * sy * g1;
    double c5 = sy * sy - s_aa + 2.0 * s_sa - s_ss;
    int ok = 0;
    double dl = nl_positive_quartic_root(c1, c2, c3, c4, c5,
                                         tolerance, &ok);
    if (!ok) return 0;

    double tlambda = 1.0 / (1.0 + hnlin * dl);
    tensor_t z = nl_sub(dev_predictor, nl_scale(backstress_prev, tlambda));
    const double magz = sqrt(nl_max(2.0 * nl_j2(z), DBL_MIN));
    tensor_t n = nl_scale(z, 1.0 / magz);
    *pstrain = nl_add(pstrain_prev, nl_scale(n, dl));
    *epbar = epbar_prev + dl;
    *backstress = nl_add(nl_scale(backstress_prev, tlambda),
                         nl_scale(n, hkin * tlambda * dl));
    *sigma = nl_sub(elastic_stress, nl_scale(n, 2.0 * c->mu * dl));
    tensor_t sdev = nl_sub(dev_predictor, nl_scale(n, 2.0 * c->mu * dl));
    *fs = sqrt(nl_max(2.0 * nl_j2(nl_sub(sdev, *backstress)), 0.0)) - sy;

    if (model == VONMISES_FAM) {
        double load_unload = 2.0 * nl_comb_j2(n, *backstress);
        double tao = sqrt(nl_max(nl_j2(sdev), 0.0));
        if ((*load_unload_state) * load_unload < 0.0 &&
            *tao_state >= *tao_max_state) {
            *tao_max_state = *tao_state;
            tensor_t stress_prev = nl_stress(strain_prev, c->mu, c->lambda);
            const double tao_elastic = sqrt(nl_max(nl_j2(nl_dev(stress_prev)),
                                                   DBL_MIN));
            su *= 2.0 / sqrt(3.0);
            const double ratio = (su + *tao_max_state) /
                                 nl_max(su - *tao_max_state, DBL_MIN);
            *psi_state = log(nl_max(ratio, DBL_MIN)) * su /
                         nl_max(tao_elastic - *tao_max_state, DBL_MIN);
            hkin = *psi_state * c->mu;
            hnlin = hkin / nl_max(sqrt(2.0) * su - sy, DBL_MIN);
            g1 = c->mu + 0.5 * hkin;
            c1 = (2.0 * c->mu * hnlin) * (2.0 * c->mu * hnlin);
            c2 = (4.0 * sy * c->mu * hnlin + 8.0 * c->mu * g1) * hnlin;
            c3 = hnlin * hnlin * (sy * sy - s_ss) + 4.0 * g1 * g1 +
                 4.0 * hnlin * sy * (c->mu + g1);
            c4 = 2.0 * hnlin * (sy * sy + s_sa - s_ss) + 4.0 * sy * g1;
            c5 = sy * sy - s_aa + 2.0 * s_sa - s_ss;
            dl = nl_positive_quartic_root(c1, c2, c3, c4, c5,
                                          tolerance, &ok);
            if (!ok) return 0;
            tlambda = 1.0 / (1.0 + hnlin * dl);
            z = nl_sub(dev_predictor, nl_scale(backstress_prev, tlambda));
            n = nl_scale(z, 1.0 / sqrt(nl_max(2.0 * nl_j2(z), DBL_MIN)));
            *pstrain = nl_add(pstrain_prev, nl_scale(n, dl));
            *epbar = epbar_prev + dl;
            *backstress = nl_add(nl_scale(backstress_prev, tlambda),
                                 nl_scale(n, hkin * tlambda * dl));
            *sigma = nl_sub(elastic_stress, nl_scale(n, 2.0 * c->mu * dl));
            sdev = nl_sub(dev_predictor, nl_scale(n, 2.0 * c->mu * dl));
            *fs = sqrt(nl_max(2.0 * nl_j2(nl_sub(sdev, *backstress)), 0.0)) - sy;
            load_unload = 2.0 * nl_comb_j2(n, *backstress);
            tao = sqrt(nl_max(nl_j2(sdev), 0.0));
            *tao_max_state = nl_max(tao, *tao_state);
        } else {
            *tao_max_state = nl_max(tao, *tao_max_state);
        }
        *load_unload_state = load_unload;
        *tao_state = tao;
    }
    return 1;
}

__device__ static void nl_jacobi_eigen(tensor_t a, double eig[3],
                                        double vec[3][3])
{
    double m[3][3] = {{a.xx, a.xy, a.xz},
                      {a.xy, a.yy, a.yz},
                      {a.xz, a.yz, a.zz}};
    vec[0][0] = vec[1][1] = vec[2][2] = 1.0;
    vec[0][1] = vec[0][2] = vec[1][0] = 0.0;
    vec[1][2] = vec[2][0] = vec[2][1] = 0.0;
    for (int sweep = 0; sweep < 16; sweep++) {
        for (int pair = 0; pair < 3; pair++) {
            const int p = (pair == 0) ? 0 : (pair == 1 ? 0 : 1);
            const int q = (pair == 0) ? 1 : (pair == 1 ? 2 : 2);
            const double apq = m[p][q];
            if (fabs(apq) <= DBL_EPSILON *
                             (fabs(m[p][p]) + fabs(m[q][q]) + 1.0)) continue;
            const double tau = (m[q][q] - m[p][p]) / (2.0 * apq);
            const double t = copysign(1.0, tau) /
                             (fabs(tau) + sqrt(1.0 + tau * tau));
            const double cs = 1.0 / sqrt(1.0 + t * t);
            const double sn = t * cs;
            const double app = m[p][p], aqq = m[q][q];
            m[p][p] = app - t * apq;
            m[q][q] = aqq + t * apq;
            m[p][q] = m[q][p] = 0.0;
            for (int r = 0; r < 3; r++) {
                if (r != p && r != q) {
                    const double mrp = m[r][p], mrq = m[r][q];
                    m[r][p] = m[p][r] = cs * mrp - sn * mrq;
                    m[r][q] = m[q][r] = sn * mrp + cs * mrq;
                }
                const double vrp = vec[r][p], vrq = vec[r][q];
                vec[r][p] = cs * vrp - sn * vrq;
                vec[r][q] = sn * vrp + cs * vrq;
            }
        }
    }
    eig[0] = m[0][0]; eig[1] = m[1][1]; eig[2] = m[2][2];
    for (int i = 0; i < 2; i++) {
        int best = i;
        for (int j = i + 1; j < 3; j++) if (eig[j] > eig[best]) best = j;
        if (best != i) {
            const double e = eig[i]; eig[i] = eig[best]; eig[best] = e;
            for (int r = 0; r < 3; r++) {
                const double v = vec[r][i]; vec[r][i] = vec[r][best];
                vec[r][best] = v;
            }
        }
    }
}

__device__ static tensor_t nl_spectral_recompose(const double eig[3],
                                                  const double v[3][3])
{
    tensor_t s = nl_zero();
    for (int a = 0; a < 3; a++) {
        s.xx += eig[a] * v[0][a] * v[0][a];
        s.yy += eig[a] * v[1][a] * v[1][a];
        s.zz += eig[a] * v[2][a] * v[2][a];
        s.xy += eig[a] * v[0][a] * v[1][a];
        s.yz += eig[a] * v[1][a] * v[2][a];
        s.xz += eig[a] * v[0][a] * v[2][a];
    }
    return s;
}

__device__ static void nl_mc_plane_return(double epbar_prev,
                                           const double trial[3],
                                           double phi, double dil,
                                           double hmod, double cohesion,
                                           double bulk, double g,
                                           double out[3], double *epbar)
{
    const double denom = 4.0 * hmod * cos(phi) * cos(phi) +
        4.0 * g * (1.0 + sin(dil) * sin(phi) / 3.0) +
        4.0 * bulk * sin(phi) * sin(dil);
    const double dg = (trial[0] - trial[2] + (trial[0] + trial[2]) * sin(phi) -
                       2.0 * (cohesion + hmod * epbar_prev) * cos(phi)) /
                      nl_max(denom, DBL_MIN);
    out[0] = trial[0] - (2.0*g*(1.0 + sin(dil)/3.0) +
                          2.0*bulk*sin(dil)) * dg;
    out[1] = trial[1] + (4.0*g/3.0 - 2.0*bulk) * sin(dil) * dg;
    out[2] = trial[2] + (2.0*g*(1.0 - sin(dil)/3.0) -
                          2.0*bulk*sin(dil)) * dg;
    *epbar = epbar_prev + 2.0 * cos(phi) * dg;
}

__device__ static void nl_mc_edge_return(double epbar_prev,
                                          const double trial[3],
                                          double phi, double dil,
                                          double hmod, double cohesion,
                                          double bulk, double g, int edge,
                                          double out[3], double *epbar)
{
    const double aux = 2.0 * cos(phi) * (cohesion + hmod * epbar_prev);
    const double fa = trial[0] - trial[2] + (trial[0] + trial[2])*sin(phi) - aux;
    const double a = 4.0*g*(1.0 + sin(phi)*sin(dil)/3.0) +
                     4.0*bulk*sin(phi)*sin(dil);
    double b, fb;
    if (edge == 1) {
        b = 2.0*g*(1.0 + sin(phi) + sin(dil) -
                   sin(phi)*sin(dil)/3.0) + 4.0*bulk*sin(phi)*sin(dil);
        fb = trial[0] - trial[1] + (trial[0] + trial[1])*sin(phi) - aux;
    } else {
        b = 2.0*g*(1.0 - sin(phi) - sin(dil) -
                   sin(phi)*sin(dil)/3.0) + 4.0*bulk*sin(phi)*sin(dil);
        fb = trial[1] - trial[2] + (trial[1] + trial[2])*sin(phi) - aux;
    }
    const double aa = a + 4.0*hmod*cos(phi)*cos(phi);
    const double bb = b + 4.0*hmod*cos(phi)*cos(phi);
    const double det = nl_max(aa*aa - bb*bb, DBL_MIN);
    const double d0 = (aa*fa - bb*fb)/det;
    const double d1 = (aa*fb - bb*fa)/det;
    const double sum = d0 + d1;
    *epbar = epbar_prev + 2.0*cos(phi)*sum;
    if (edge == 1) {
        out[0] = trial[0] - (2.0*g*(1.0+sin(dil)/3.0)+2.0*bulk*sin(dil))*sum;
        out[1] = trial[1] + (4.0*g/3.0-2.0*bulk)*sin(dil)*d0 +
                 (2.0*g*(1.0-sin(dil)/3.0)-2.0*bulk*sin(dil))*d1;
        out[2] = trial[2] + (2.0*g*(1.0-sin(dil)/3.0)-2.0*bulk*sin(dil))*d0 +
                 (4.0*g/3.0-2.0*bulk)*sin(dil)*d1;
    } else {
        out[0] = trial[0] - (2.0*g*(1.0+sin(dil)/3.0)+2.0*bulk*sin(dil))*d0 +
                 (4.0*g/3.0-2.0*bulk)*sin(dil)*d1;
        out[1] = trial[1] + (4.0*g/3.0-2.0*bulk)*sin(dil)*d0 -
                 (2.0*g*(1.0+sin(dil)/3.0)+2.0*bulk*sin(dil))*d1;
        out[2] = trial[2] + (2.0*g*(1.0-sin(dil)/3.0)-2.0*bulk*sin(dil))*sum;
    }
}

__device__ static void nl_mc_apex_return(double epbar_prev, double p_trial,
                                          double phi, double dil,
                                          double hmod, double cohesion,
                                          double bulk, double out[3],
                                          double *epbar)
{
    const double cot = cos(phi) / nl_max(sin(phi), DBL_MIN);
    const double omega = sin(dil) / nl_max(cos(phi), DBL_MIN);
    const double dep = (p_trial - (cohesion + hmod*epbar_prev)*cot) /
                       nl_max(hmod*cot + omega*bulk, DBL_MIN);
    *epbar = epbar_prev + dep;
    out[0] = out[1] = out[2] = p_trial - bulk*omega*dep;
}

__device__ static int nl_mc_corner_zone(const double sigma[3], double sp,
                                         double updated[3], double phi[5])
{
    const double tolerance = -1.0e-10;
    phi[0] = sigma[0];
    phi[1] = sp-sigma[2];
    phi[2] = sigma[1];
    phi[3] = sigma[2];
    phi[4] = sp-sigma[1];
    int active = 0;
    for (int i = 0; i < 5; i++) active += (phi[i] > 0.0);
    if (active > 3) return -1;
    if (phi[2] > tolerance && phi[3] > tolerance) {
        updated[0] = updated[1] = updated[2] = 0.0;
        return 1;
    }
    if (phi[3] > tolerance && phi[4] > tolerance) {
        updated[0] = 0.0; updated[1] = sp; updated[2] = 0.0;
        return 3;
    }
    if (phi[1] > tolerance && phi[4] > tolerance) {
        updated[0] = 0.0; updated[1] = sp; updated[2] = sp;
        return 5;
    }
    if (phi[1] > tolerance && phi[2] > tolerance) {
        updated[0] = 0.0; updated[1] = 0.0; updated[2] = sp;
        return 7;
    }
    updated[0] = updated[1] = updated[2] = 0.0;
    return 0;
}

__device__ static int nl_mc_tension_return(double bulk, double g,
                                            double phi_angle,
                                            double cohesion,
                                            const double trial[3],
                                            double updated[3])
{
    const double nphi = (1.0+sin(phi_angle)) /
                        nl_max(1.0-sin(phi_angle), DBL_MIN);
    const double sp = -2.0*cohesion*sqrt(nphi);
    const double alpha1 = bulk+4.0*g/3.0;
    const double alpha2 = bulk-2.0*g/3.0;
    double first[3], second[3], planes[5];
    int zone = nl_mc_corner_zone(trial, sp, first, planes);
    if (zone < 0) return 0;
    if (zone != 0) {
        updated[0] = first[0]; updated[1] = first[1]; updated[2] = first[2];
        return 1;
    }

    const double plane1 = planes[0];
    double plane2 = 0.0;
    int position = 0;
    for (int i = 1; i < 5; i++) {
        if (planes[i] > 0.0) {
            plane2 = planes[i];
            position = i;
            break;
        }
    }
    if (plane2 != 0.0) {
        double bb = alpha2;
        if (position == 2 || position == 5) bb = -bb;
        const double determinant = alpha1*alpha1-bb*bb;
        const double dl1 = (alpha1*plane1-bb*plane2)/determinant;
        const double dl2 = (-bb*plane1+alpha1*plane2)/determinant;
        if (position == 2) {
            second[0] = trial[0]-(alpha1*dl1-alpha2*dl2);
            second[1] = trial[1]-(alpha2*dl1-alpha2*dl2);
            second[2] = trial[2]-(alpha2*dl1-alpha1*dl2);
        } else if (position == 3) {
            second[0] = trial[0]-(alpha1*dl1+alpha2*dl2);
            second[1] = trial[1]-(alpha2*dl1+alpha1*dl2);
            second[2] = trial[2]-(alpha2*dl1+alpha2*dl2);
        } else if (position == 4) {
            second[0] = trial[0]-(alpha1*dl1+alpha2*dl2);
            second[1] = trial[1]-(alpha2*dl1+alpha2*dl2);
            second[2] = trial[2]-(alpha2*dl1+alpha1*dl2);
        } else {
            second[0] = trial[0]-(alpha1*dl1-alpha2*dl2);
            second[1] = trial[1]-(alpha2*dl1-alpha1*dl2);
            second[2] = trial[2]-(alpha2*dl1-alpha2*dl2);
        }
    } else {
        second[0] = trial[0]-plane1;
        second[1] = trial[1]-alpha2*plane1/alpha1;
        second[2] = trial[2]-alpha2*plane1/alpha1;
    }

    zone = nl_mc_corner_zone(second, sp, updated, planes);
    if (zone < 0) return 0;
    if (zone == 0) {
        updated[0] = second[0]; updated[1] = second[1]; updated[2] = second[2];
        if (updated[2] < sp) updated[2] = sp;
        if (updated[1] < sp) updated[1] = sp;
    }
    return 1;
}

__device__ static int nl_update_mohr_coulomb(const nlconstants_t *c,
                                              tensor_t strain,
                                              tensor_t sigma_trial,
                                              tensor_t static_stress,
                                              tensor_t *pstrain,
                                              tensor_t *sigma,
                                              double epbar_prev,
                                              double *epbar,
                                              double *fs,
                                              noyesflag_t tension_cutoff)
{
    const double bulk = c->lambda + 2.0*c->mu/3.0;
    double eig[3], vec[3][3];
    nl_jacobi_eigen(sigma_trial, eig, vec);
    if (tension_cutoff == YES) {
        const double nphi = (1.0 + sin(c->phi)) /
                            nl_max(1.0 - sin(c->phi), DBL_MIN);
        const double alpha_p = sqrt(1.0+nphi*nphi)+nphi;
        const double sp = -2.0*c->c*sqrt(nphi);
        const double shear_tension = eig[0]+alpha_p*(eig[2]-sp);
        if (shear_tension > 0.0 && eig[0] > 0.0) {
            double corrected[3];
            if (!nl_mc_tension_return(bulk, c->mu, c->phi, c->c,
                                      eig, corrected)) return 0;
            const tensor_t total_stress =
                nl_spectral_recompose(corrected, vec);
            *sigma = nl_sub(total_stress, static_stress);
            *pstrain = nl_sub(strain,
                              nl_elastic_strain(*sigma, c->mu, bulk));
            *epbar = 0.0;
            *fs = corrected[0];
            return 1;
        }
        if (shear_tension > 0.0 && eig[0] <= 0.0) {
            *sigma = nl_sub(sigma_trial, static_stress);
            *epbar = epbar_prev;
            *fs = eig[0];
            return 1;
        }
    }
    double yield = eig[0] - eig[2] + (eig[0] + eig[2])*sin(c->phi) -
                   2.0*(c->c + c->h*epbar_prev)*cos(c->phi);
    if (yield < 0.0) {
        *sigma = nl_sub(sigma_trial, static_stress);
        *epbar = epbar_prev;
        *fs = yield;
        return 1;
    }
    double corrected[3];
    nl_mc_plane_return(epbar_prev, eig, c->phi, c->dil_angle, c->h, c->c,
                       bulk, c->mu, corrected, epbar);
    if (corrected[0] <= corrected[1] || corrected[1] <= corrected[2]) {
        const int edge = ((1.0-sin(c->dil_angle))*eig[0] - 2.0*eig[1] +
                          (1.0+sin(c->dil_angle))*eig[2] > 0.0) ? 1 : -1;
        nl_mc_edge_return(epbar_prev, eig, c->phi, c->dil_angle, c->h, c->c,
                          bulk, c->mu, edge, corrected, epbar);
        const double cond1 = corrected[0]-corrected[1];
        const double cond2 = corrected[1]-corrected[2];
        if ((cond1 <= 0.0 && fabs(cond1) >= 1.0e-3) ||
            (cond2 <= 0.0 && fabs(cond2) >= 1.0e-3)) {
            nl_mc_apex_return(epbar_prev, (eig[0]+eig[1]+eig[2])/3.0,
                              c->phi, c->dil_angle, c->h, c->c, bulk,
                              corrected, epbar);
        }
    }
    const tensor_t total_stress = nl_spectral_recompose(corrected, vec);
    *sigma = nl_sub(total_stress, static_stress);
    *pstrain = nl_sub(strain, nl_elastic_strain(*sigma, c->mu, bulk));
    *fs = corrected[0]-corrected[2] +
          (corrected[0]+corrected[2])*sin(c->phi) -
          2.0*(c->c+c->h*(*epbar))*cos(c->phi);
    return isfinite(*fs);
}

__device__ static int nl_update_material(const nlconstants_t *c,
                                          const nonlinear_gpu_config_t *cfg,
                                          tensor_t strain,
                                          tensor_t strain_prev,
                                          tensor_t static_stress,
                                          tensor_t pstrain_prev,
                                          tensor_t backstress_prev,
                                          double epbar_prev,
                                          tensor_t *pstrain,
                                          tensor_t *backstress,
                                          tensor_t *sigma,
                                          double *epbar,
                                          double *fs,
                                          double *psi_state,
                                          double *load_unload_state,
                                          double *tao_state,
                                          double *tao_max_state,
                                          double *kappa,
                                          tensor_t *sref,
                                          double *gamma1d,
                                          double *tao1d,
                                          double *ggmax,
                                          double *model_error,
                                          tensor_t *fabric,
                                          double *md_initialized,
                                          double *md_integration_error,
                                          int32_t *md_substeps,
                                          nl_md_diagnostics_t *md_diagnostics)
{
    const materialmodel_t model = cfg->material_model;

#ifndef NL_GPU_CLASSIC_ONLY
    if (model == MANZARI_DAFALIAS) {
        *pstrain = pstrain_prev;
        *backstress = backstress_prev;
        int32_t accepted = 0;
        const int ok = nl_update_manzari_dafalias(
            c, cfg, strain, strain_prev, static_stress, sigma, pstrain,
            backstress, sref, fabric, md_initialized, epbar, psi_state, fs,
            md_integration_error, &accepted, md_diagnostics);
        *load_unload_state = (double)accepted;
        *md_substeps = accepted;
        *model_error = *fs;
        return ok;
    }
#endif

    /* Cyclic bounding/backbone models decide loading and reversal from the
     * strain increment and their stored stress/history. They must not pass
     * through the conventional total-strain yield predictor. */
    if (MATERIAL_MODEL_IS_BACKBONE(model)) {
        return nl_update_backbone(c, model, strain, strain_prev, sref, sigma,
                                  kappa, gamma1d, tao1d, ggmax, cfg,
                                  model_error);
    }

    tensor_t elastic_strain = nl_sub(strain, pstrain_prev);
    tensor_t stress_predictor = nl_stress(elastic_strain, c->mu, c->lambda);
    const tensor_t total_predictor = nl_add(stress_predictor, static_stress);
    tensor_t relative_dev = nl_sub(nl_dev(total_predictor), backstress_prev);
    const double j2 = nl_j2(relative_dev);
    const double i1 = nl_i1(total_predictor);
    double hardening;
    if (model == VONMISES_EP) hardening = c->c * 2.0/sqrt(3.0);
    else if (model == DRUCKERPRAGER)
        hardening = c->gamma*(c->c+c->h*epbar_prev);
    else if (model == VONMISES_FA) hardening = c->Sstrain0*c->mu;
    else if (model == MOHR_COULOMB) {
        double eig[3], vec[3][3];
        nl_jacobi_eigen(total_predictor, eig, vec);
        hardening = 2.0*(c->c+c->h*epbar_prev)*cos(c->phi);
        const double pred = eig[0]-eig[2]+(eig[0]+eig[2])*sin(c->phi)-hardening;
        if (pred < 1.0e-3) {
            *pstrain = pstrain_prev; *sigma = stress_predictor;
            *epbar = epbar_prev; *fs = pred; return 1;
        }
        return nl_update_mohr_coulomb(c, strain, total_predictor,
                                      static_stress, pstrain, sigma,
                                      epbar_prev, epbar, fs,
                                      cfg->tension_cutoff);
    } else hardening = 0.0;

    const double surface = ((model == DRUCKERPRAGER) ? c->alpha*i1 : 0.0) +
                           sqrt(nl_max(j2, 0.0));
    const double predictor = surface-hardening;
    if (model == LINEAR || predictor < 1.0e-3) {
        *pstrain = pstrain_prev;
        *backstress = backstress_prev;
        *sigma = stress_predictor;
        *epbar = epbar_prev;
        *fs = predictor;
        return 1;
    }

    if (model == VONMISES_FA || model == VONMISES_FAM) {
        const tensor_t restored_dev = nl_add(relative_dev, backstress_prev);
        return nl_update_fa(c, model, restored_dev, strain_prev,
                            pstrain_prev, backstress_prev, epbar_prev,
                            stress_predictor, pstrain, backstress, sigma,
                            epbar, fs, psi_state, load_unload_state,
                            tao_state, tao_max_state,
                            cfg->error_tolerance);
    }

    const double bulk = c->lambda + 2.0*c->mu/3.0;
    const double rootj2 = sqrt(nl_max(j2, DBL_MIN));
    tensor_t direction = nl_scale(relative_dev, 1.0/(2.0*rootj2));
    direction.xx += c->beta; direction.yy += c->beta; direction.zz += c->beta;
    const double dl = predictor /
        nl_max(c->mu + 9.0*bulk*c->alpha*c->beta +
               c->h*c->gamma*c->gamma, DBL_MIN);
    *pstrain = nl_add(pstrain_prev, nl_scale(direction, dl));
    *sigma = nl_stress(nl_sub(strain, *pstrain), c->mu, c->lambda);
    *epbar = epbar_prev + dl*c->gamma;
    const tensor_t updated_total = nl_add(*sigma, static_stress);
    tensor_t updated_dev = nl_dev(updated_total);
    double updated_j2 = nl_j2(updated_dev);
    *fs = ((model == DRUCKERPRAGER) ? c->alpha*nl_i1(updated_total) : 0.0) +
          sqrt(nl_max(updated_j2, 0.0));
    if (model == VONMISES_EP) *fs -= c->c*2.0/sqrt(3.0);
    else *fs -= c->gamma*(c->c+c->h*(*epbar));

    if (model == DRUCKERPRAGER && rootj2-c->mu*dl < 0.0) {
        const double dl_apex = (c->alpha*i1-
            c->gamma*(c->c+c->h*epbar_prev)) /
            nl_max(9.0*bulk*c->alpha*c->beta+
                   c->h*c->gamma*c->gamma, DBL_MIN);
        *epbar = epbar_prev+dl_apex*c->gamma;
        const double skk = i1-9.0*bulk*c->beta*dl_apex;
        *sigma = nl_sub(nl_iso(skk/3.0), static_stress);
        *pstrain = nl_sub(strain, nl_elastic_strain(*sigma, c->mu, bulk));
        *fs = c->alpha*skk-c->gamma*(c->c+c->h*(*epbar));
    }
    *model_error = fabs(*fs);
    return isfinite(*fs);
}

#ifndef NL_GPU_PRODUCTION_ONLY
__global__ static void kernel_nonlinear_material_point_history(
    nlconstants_t constants,
    nonlinear_gpu_config_t config,
    nonlinear_material_point_initial_state_t initial_state,
    const tensor_t *strain_history,
    int32_t steps,
    nonlinear_material_point_result_t *results)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) return;

    const int seeded = (config.material_model == MANZARI_DAFALIAS ||
                        config.material_model == PRESSURE_DEPEND_MULTIYIELD_03 ||
                        config.material_model == PJ_LIQ_V8) &&
                       initial_state.provided != 0;
    tensor_t strain_prev = seeded ? initial_state.strain : nl_zero();
    tensor_t plastic_strain = seeded ? initial_state.plastic_strain : nl_zero();
    tensor_t backstress = seeded ? initial_state.alpha : nl_zero();
    tensor_t stress_reference = seeded ? initial_state.alpha_in : nl_zero();
    tensor_t fabric = seeded ? initial_state.fabric : nl_zero();
    tensor_t stress = seeded ? initial_state.stress : nl_zero();
    double ep_bar = 0.0;
    double psi = 0.0;
    double load_unload = 0.0;
    double tao = 0.0;
    double tao_max = 0.0;
    double kappa = MATERIAL_MODEL_IS_BACKBONE(config.material_model)
        ? 100.0 : 0.0;
    double gamma1d = MATERIAL_MODEL_IS_BACKBONE(config.material_model)
        ? 1.0e-10 : 0.0;
    double tao1d = 0.0;
    double ggmax = MATERIAL_MODEL_IS_BACKBONE(config.material_model)
        ? 1.0 : 0.0;
    double md_initialized = seeded ? 1.0 : 0.0;
    double md_integration_error = 0.0;
    int32_t md_substeps = 0;
    pdmy03_state_t pdmy03_state = {};
    tensor_t pdmy03_centers[PDMY03_MAX_SURFACES] = {};
    pj8_state_t pj8_state = {};
    pj8_parameters_t pj8_parameters = pj8_reference_parameters(1000.0);
    if (config.material_model == PJ_LIQ_V8 && seeded) {
        tensor_t reference = stress;
        reference.xy = reference.yz = reference.xz = 0.0;
        if (!pj8_initialize(&pj8_parameters, stress,
                            constants.pj8_void_ratio, &pj8_state) ||
            !pj8_begin_dynamic_phase(&pj8_parameters, reference,
                                     &pj8_state)) return;
    }

    for (int32_t step = 0; step < steps; step++) {
        tensor_t plastic_next = plastic_strain;
        tensor_t backstress_next = backstress;
        double ep_bar_next = ep_bar;
        double fs = 0.0;
        double model_error = 0.0;
        nl_md_diagnostics_t md_diagnostics = {};
        pdmy03_diagnostics_t pdmy03_diagnostics = {};
        pj8_update_info_t pj8_info = {};

        int ok;
        if (config.material_model == PJ_LIQ_V8) {
            ok = pj8_update(&pj8_parameters,
                pj8_sub(strain_history[step], strain_prev),
                config.no_substeps, &pj8_state, &stress, &pj8_info);
            fs = fabs(pj8_compatibility_residual(&pj8_state));
            model_error = fs;
        } else if (config.material_model == PRESSURE_DEPEND_MULTIYIELD_03) {
            ok = pdmy03_update(
                &constants.pdmy03, config.no_substeps,
                strain_history[step], strain_prev, nl_zero(), &stress,
                &plastic_next, &pdmy03_state, pdmy03_centers,
                &pdmy03_diagnostics);
            fs = pdmy03_diagnostics.yield_residual;
            model_error = fs;
        } else {
            ok = nl_update_material(
                &constants, &config, strain_history[step], strain_prev, nl_zero(),
                plastic_strain, backstress, ep_bar,
                &plastic_next, &backstress_next, &stress, &ep_bar_next, &fs,
                &psi, &load_unload, &tao, &tao_max, &kappa,
                &stress_reference, &gamma1d, &tao1d, &ggmax, &model_error,
                &fabric, &md_initialized, &md_integration_error, &md_substeps,
                &md_diagnostics);
        }

        nonlinear_material_point_result_t result = {};
        result.stress = stress;
        result.plastic_strain = plastic_next;
        result.backstress = backstress_next;
        result.ep_bar = ep_bar_next;
        result.fs = MATERIAL_MODEL_IS_BACKBONE(config.material_model)
            ? model_error : fs;
        result.psi = psi;
        result.load_unload = load_unload;
        result.tao = tao;
        result.tao_max = tao_max;
        result.kappa = kappa;
        result.stress_reference = stress_reference;
        result.gamma1d = gamma1d;
        result.tao1d = tao1d;
        result.ggmax = ggmax;
        result.model_error = MATERIAL_MODEL_IS_BACKBONE(config.material_model)
            ? model_error : 0.0;
        if (config.material_model == MANZARI_DAFALIAS) {
            result.fs = fs;
            result.model_error = model_error;
            result.fabric = fabric;
            result.alpha_in = stress_reference;
            result.void_ratio = ep_bar_next;
            result.plastic_multiplier = psi;
            result.integration_error = md_integration_error;
            result.accepted_substeps = md_substeps;
        } else if (config.material_model == PRESSURE_DEPEND_MULTIYIELD_03) {
            result.fs = pdmy03_diagnostics.yield_residual;
            result.model_error = pdmy03_diagnostics.yield_residual;
            result.accepted_substeps = pdmy03_diagnostics.accepted_substeps;
            result.pdmy03_active_surface = pdmy03_state.active_surface;
            result.pdmy03_requested_substeps =
                pdmy03_diagnostics.requested_substeps;
            result.pdmy03_crossed_surfaces =
                pdmy03_diagnostics.crossed_surfaces;
            result.pdmy03_substeps_clamped =
                pdmy03_diagnostics.substeps_clamped;
            result.pdmy03_pressure_d = pdmy03_state.pressure_d;
            double center_norm = 0.0;
            for (int i = 0; i < constants.pdmy03.surface_count; i++)
                center_norm += pdmy03_ddot(pdmy03_centers[i],
                                           pdmy03_centers[i]);
            result.pdmy03_center_norm = sqrt(center_norm);
        } else if (config.material_model == PJ_LIQ_V8) {
            result.fs = fabs(pj8_compatibility_residual(&pj8_state));
            result.model_error = result.fs;
            result.backstress = pj8_state.alpha;
            result.ep_bar = pj8_state.lambda_total;
            result.psi = pj8_state.D;
            result.load_unload = (double)pj8_info.accepted_substeps;
            result.kappa = pj8_pressure(pj8_state.stress);
            result.stress_reference = pj8_state.alpha01;
            result.gamma1d = pj8_state.eps_v_irreversible;
            result.tao1d = pj8_state.eps_v_reversible;
            result.ggmax = pj8_state.amplitude_factor;
            result.fabric = pj8_state.fabric;
            result.alpha_in = pj8_state.alpha0;
            result.accepted_substeps = pj8_info.accepted_substeps;
            result.void_ratio = pj8_state.void_ratio;
            result.plastic_multiplier = pj8_state.lambda_total;
            pj8_state_values(&pj8_state, result.pj8_state_values);
            result.pj8_reversal_registered =
                pj8_info.reversal_registered;
        }
        result.status = ok ? 0 : 1;
        results[step] = result;

        /* Every history variable returned by step n is the previous state for
         * step n+1. The production element loop must obey the same rule. */
        strain_prev = strain_history[step];
        plastic_strain = plastic_next;
        backstress = backstress_next;
        ep_bar = ep_bar_next;
    }
}
#endif

#if defined(NL_GPU_CLASSIC_ONLY) || \
    defined(NL_GPU_MANZARI_DAFALIAS_ONLY) || defined(NL_GPU_PDMY03_ONLY) || \
    defined(NL_GPU_PJ_LIQ_V8_ONLY)
__device__ static inline void nl_accumulate_stress_force(
    fvector_t local_force[8], int node, tensor_t stress,
    double dx, double dy, double dz, double weight)
{
    local_force[node].f[0] += (solver_float)(
        (dx*stress.xx + dy*stress.xy + dz*stress.xz)*weight);
    local_force[node].f[1] += (solver_float)(
        (dy*stress.yy + dx*stress.xy + dz*stress.yz)*weight);
    local_force[node].f[2] += (solver_float)(
        (dz*stress.zz + dy*stress.yz + dx*stress.xz)*weight);
}

__device__ static void nl_hex_internal_force(const qptensors_t *stresses,
                                              double h,
                                              fvector_t local_force[8])
{
    const double q = 0.577350269189;
    const double weight = h*h*h*0.125;
    for (int node = 0; node < 8; node++) {
        for (int qp = 0; qp < 8; qp++) {
            const double lx = (qp & 1) ? q : -q;
            const double ly = (qp & 2) ? q : -q;
            const double lz = (qp & 4) ? q : -q;
            double dx, dy, dz;
            nl_point_dxi(&dx, &dy, &dz, lx, ly, lz, h, node);
            nl_accumulate_stress_force(local_force, node, stresses->qp[qp],
                                       dx, dy, dz, weight);
        }
    }
}

__device__ static void nl_tetra_internal_force(const qptensors_t *stresses,
                                                const nlconstants_t *c,
                                                double h,
                                                fvector_t local_force[8])
{
    for (int tet = 0; tet < 5; tet++) {
        if (c->tetraVol[tet] == 0.0) continue;
        int nodes[4];
        double grad[4][3], volume;
        if (!nl_tetra_gradients(c->topoPart, tet, h, nodes, grad, &volume))
            continue;
        const double weight = volume*c->tetraVol[tet];
        for (int a = 0; a < 4; a++) {
            nl_accumulate_stress_force(local_force, nodes[a],
                                       stresses->qp[tet],
                                       grad[a][0], grad[a][1], grad[a][2],
                                       weight);
        }
    }
}

__device__ static void nl_elastic_internal_force(const fvector_t u[8],
                                                  const nlconstants_t *c,
                                                  double h,
                                                  fvector_t local_force[8])
{
    qptensors_t stresses = {};
    if (c->isTopoNonlin == 0) {
        for (int qp = 0; qp < 8; qp++)
            stresses.qp[qp] = nl_stress(nl_hex_strain(u, qp, h),
                                        c->mu, c->lambda);
        nl_hex_internal_force(&stresses, h, local_force);
    } else {
        for (int tet = 0; tet < 5; tet++) {
            if (c->tetraVol[tet] != 0.0)
                stresses.qp[tet] = nl_stress(
                    nl_tetra_strain(u, c->topoPart, tet, h),
                    c->mu, c->lambda);
        }
        nl_tetra_internal_force(&stresses, c, h, local_force);
    }
}

__device__ static inline void nl_atomic_add(solver_float *address,
                                             solver_float value)
{
    atomicAdd(address, value);
}

__device__ static tensor_t nl_approx_static_stress(
    const nlconstants_t *c, const nonlinear_gpu_config_t *cfg,
    double h, int qp)
{
    if (cfg->approximate_geostatic_state != YES || c->sigmaZ_st <= 0.0)
        return nl_zero();
    const double lz = c->isTopoNonlin ? 0.0 :
        ((qp & 4) ? 0.577350269189 : -0.577350269189);
    const double vertical = nl_max(
        c->sigmaZ_st+c->sigmaZ_gradient*h*lz*0.5, 0.0);
    const double k0 = 1.0-sin(c->phi);
    tensor_t stress = {-k0*vertical, -k0*vertical, -vertical,
                       0.0, 0.0, 0.0};
    return stress;
}

#ifdef NL_GPU_MANZARI_DAFALIAS_ONLY
__device__ static tensor_t nl_approx_static_stress_k0(
    const nlconstants_t *c, const nonlinear_gpu_config_t *cfg,
    double geostatic_k0, double h, int qp)
{
    if (cfg->approximate_geostatic_state != YES || c->sigmaZ_st <= 0.0)
        return nl_zero();
    const double lz = c->isTopoNonlin ? 0.0 :
        ((qp & 4) ? 0.577350269189 : -0.577350269189);
    const double vertical = nl_max(
        c->sigmaZ_st+c->sigmaZ_gradient*h*lz*0.5, 0.0);
    tensor_t stress = {-geostatic_k0*vertical,
                       -geostatic_k0*vertical, -vertical,
                       0.0, 0.0, 0.0};
    return stress;
}
#endif

#if defined(NL_GPU_MANZARI_DAFALIAS_ONLY)
/* Diagnostics are consumed only after a failed constitutive update. Keeping
 * one record per thread in shared memory prevents that failure-only record
 * from enlarging every production thread's local stack. The extra double
 * gives consecutive records a 152-byte stride, avoiding a power-of-two
 * shared-memory bank stride for warp-wide member stores. */
typedef struct nl_md_shared_diagnostics_t {
    nl_md_diagnostics_t value;
    double bank_padding;
} nl_md_shared_diagnostics_t;
#endif

#if defined(NL_GPU_MANZARI_DAFALIAS_ONLY)
/* Minimum resident blocks per SM requested from ptxas, per integration scheme.
 *
 * Without a bound ptxas allocates up to the 255-register ceiling, which on
 * sm_100 admits only 65536/255 = 257 threads, i.e. 2 blocks of
 * NL_GPU_UPDATE_BLOCK_SIZE and 8 of 64 warps. Raising the requested residency
 * lowers the per-thread budget (65536/(BLOCK*N)) and buys latency hiding at
 * the cost of more spilling, so the useful setting is empirical: the measured
 * demand Forward-Euler entry spills only a few hundred bytes at 255 registers
 * and may have room to trade, while Modified Euler already spills tens of
 * kilobytes and is left alone. Recheck every entry after changing the bound.
 *
 * Override per build to sweep, e.g. -DNL_MD_FE_MIN_BLOCKS_PER_SM=3. The
 * neutral default of 2 reproduces the previous allocation on the measured
 * sm_100 build until target-system timing justifies a higher residency. */
#ifndef NL_MD_FE_MIN_BLOCKS_PER_SM
#define NL_MD_FE_MIN_BLOCKS_PER_SM 2
#endif
#ifndef NL_MD_ME_MIN_BLOCKS_PER_SM
#define NL_MD_ME_MIN_BLOCKS_PER_SM 2
#endif

/* One entry function per integration scheme and Forward-Euler subdivision
 * policy. ptxas allocates registers per entry function, so specializing here
 * keeps both Forward-Euler paths from inheriting Modified Euler's allocation
 * and spill traffic. */
template <nl_md_integration_scheme_t MD_SCHEME,
          nl_md_forward_euler_substeps_t MD_FE_SUBSTEPS>
__global__ static void
__launch_bounds__(NL_GPU_UPDATE_BLOCK_SIZE,
                  MD_SCHEME == NL_MD_INTEGRATION_FORWARD_EULER
                      ? NL_MD_FE_MIN_BLOCKS_PER_SM
                      : NL_MD_ME_MIN_BLOCKS_PER_SM)
kernel_nonlinear_update_force(
    nonlinear_gpu_state_t state, nonlinear_gpu_config_t cfg,
    gpu_data_t *gpu, int32_t step, double delta_t, double delta_t_squared)
#else
__global__ static void kernel_nonlinear_update_force(
    nonlinear_gpu_state_t state, nonlinear_gpu_config_t cfg,
    gpu_data_t *gpu, int32_t step, double delta_t, double delta_t_squared)
#endif
{
#if defined(NL_GPU_MANZARI_DAFALIAS_ONLY)
    __shared__ nl_md_shared_diagnostics_t
        md_diagnostics_shared[NL_GPU_UPDATE_BLOCK_SIZE];
#endif
    const int32_t nl_eindex = (int32_t)(blockIdx.x*blockDim.x + threadIdx.x);
    if (nl_eindex >= state.count) return;

    const nonlinear_gpu_element_t element = state.elements[nl_eindex];
    const int32_t eindex = element.eindex;
    nlconstants_t *c = state.constants + nl_eindex;
    fvector_t u[8], du[8];
    int any_displacement = 0;
    for (int node = 0; node < 8; node++) {
        const int32_t lnid = gpu->lnidArrayDevice[eindex + node*gpu->lenum];
        u[node] = gpu->tm1Device[lnid];
        for (int d = 0; d < 3; d++) {
            any_displacement |= (u[node].f[d] != (solver_float)0);
            du[node].f[d] = 0;
        }
    }

#ifdef NL_GPU_CLASSIC_ONLY
    if (step == 0 && MATERIAL_MODEL_IS_BACKBONE(cfg.material_model)) {
        for (int qp = 0; qp < NL_GPU_QP; qp++) {
            state.kappa[nl_eindex].qv[qp] = 100.0;
            state.gamma1d[nl_eindex].qv[qp] = 1.0e-10;
            state.ggmax[nl_eindex].qv[qp] = 1.0;
        }
    }
#endif

    int update_constitutive_state = any_displacement;
#ifdef NL_GPU_PJ_LIQ_V8_ONLY
    /* A cyclic element can return to exactly zero absolute displacement from a
     * nonzero committed strain.  The old absolute-displacement shortcut would
     * skip that real unloading increment and assemble the previous V8 stress.
     * Inspect committed strain only after the point has initialized, retaining
     * the inexpensive untouched-wavefront path. */
    if (!update_constitutive_state) {
        const pj8_state_t *pj8_states =
            (const pj8_state_t *)state.pj8_states;
        int history_qp = 0;
        if (c->isTopoNonlin) {
            while (history_qp < 5 && c->tetraVol[history_qp] == 0.0)
                history_qp++;
        }
        if (history_qp < (c->isTopoNonlin ? 5 : NL_GPU_QP) &&
            pj8_states[nl_eindex*NL_GPU_QP+history_qp].initialized) {
            const int ngp = c->isTopoNonlin ? 5 : NL_GPU_QP;
            for (int qp = 0; qp < ngp; qp++) {
                if (c->isTopoNonlin && c->tetraVol[qp] == 0.0) continue;
                const tensor_t previous = state.strains[nl_eindex].qp[qp];
                if (previous.xx != 0.0 || previous.yy != 0.0 ||
                    previous.zz != 0.0 || previous.xy != 0.0 ||
                    previous.yz != 0.0 || previous.xz != 0.0) {
                    update_constitutive_state = 1;
                    break;
                }
            }
        }
    }
#endif

    if (update_constitutive_state) {
#ifdef NL_GPU_MANZARI_DAFALIAS_ONLY
        const double geostatic_k0 =
            step > cfg.geostatic_final_step &&
            cfg.approximate_geostatic_state == YES && c->sigmaZ_st > 0.0
                ? 1.0-sin(c->phi) : 0.0;
#endif
        const int ngp = c->isTopoNonlin ? 5 : 8;
        for (int qp = 0; qp < ngp; qp++) {
            if (c->isTopoNonlin && c->tetraVol[qp] == 0.0) continue;

            const tensor_t strain_prev = state.strains[nl_eindex].qp[qp];
#if defined(NL_GPU_CLASSIC_ONLY) || defined(NL_GPU_PDMY03_ONLY)
            state.strains1[nl_eindex].qp[qp] = strain_prev;
#endif
            const tensor_t strain = c->isTopoNonlin
                ? nl_tetra_strain(u, c->topoPart, qp, element.h)
                : nl_hex_strain(u, qp, element.h);
            state.strains[nl_eindex].qp[qp] = strain;
#if defined(NL_GPU_CLASSIC_ONLY) || defined(NL_GPU_PDMY03_ONLY)
            state.pstrains1[nl_eindex].qp[qp] =
                state.pstrains2[nl_eindex].qp[qp];
            state.alphastress1[nl_eindex].qp[qp] =
                state.alphastress2[nl_eindex].qp[qp];
            state.ep1[nl_eindex].qv[qp] = state.ep2[nl_eindex].qv[qp];
#endif

            if (cfg.material_model == LINEAR ||
                step <= cfg.geostatic_final_step) {
                state.stresses[nl_eindex].qp[qp] =
                    nl_stress(strain, c->mu, c->lambda);
                continue;
            }

            double model_error = 0.0;
#ifdef NL_GPU_MANZARI_DAFALIAS_ONLY
            const tensor_t static_stress = nl_approx_static_stress_k0(
                c, &cfg, geostatic_k0, element.h, qp);
#else
            const tensor_t static_stress =
                nl_approx_static_stress(c, &cfg, element.h, qp);
#endif
            int ok;
#ifdef NL_GPU_PDMY03_ONLY
            pdmy03_diagnostics_t pdmy03_diagnostics = {};
            const int32_t point = nl_eindex*NL_GPU_QP+qp;
            tensor_t *centers = state.pdmy03_centers+
                (size_t)point*state.pdmy03_surface_stride;
            ok = pdmy03_update(
                &c->pdmy03, cfg.no_substeps, strain, strain_prev,
                static_stress, &state.stresses[nl_eindex].qp[qp],
                &state.pstrains2[nl_eindex].qp[qp],
                state.pdmy03_states+point, centers,
                &pdmy03_diagnostics);
            model_error = pdmy03_diagnostics.yield_residual;
            c->fs[qp] = model_error;
            if (!ok || !isfinite(model_error)) {
                nl_record_pdmy03_error(
                    state.error, eindex, qp, &pdmy03_diagnostics,
                    pdmy03_scale(pdmy03_add(
                        state.stresses[nl_eindex].qp[qp], static_stress),
                        -1.0),
                    pdmy03_sub(strain, strain_prev));
            }
#elif defined(NL_GPU_PJ_LIQ_V8_ONLY)
            const int32_t point = nl_eindex*NL_GPU_QP+qp;
            pj8_state_t *pj8_states = (pj8_state_t *)state.pj8_states;
            const pj8_parameters_t *pj8_parameters =
                (const pj8_parameters_t *)state.pj8_parameters;
            pj8_state_t *pj8_state = pj8_states+point;
            if (!pj8_state->initialized) {
                const tensor_t equilibrated = pj8_add(
                    state.stresses[nl_eindex].qp[qp], static_stress);
                tensor_t zero_static_shear = equilibrated;
                zero_static_shear.xy = 0.0;
                zero_static_shear.yz = 0.0;
                zero_static_shear.xz = 0.0;
                ok = pj8_initialize(pj8_parameters, equilibrated,
                                    c->pj8_void_ratio, pj8_state) &&
                     pj8_begin_dynamic_phase(pj8_parameters,
                                             zero_static_shear, pj8_state);
            } else {
                ok = 1;
            }
            pj8_update_info_t pj8_info = {};
            tensor_t full_effective_stress = pj8_state->stress;
            if (ok) {
                ok = pj8_update(pj8_parameters,
                    pj8_sub(strain, strain_prev), cfg.no_substeps,
                    pj8_state, &full_effective_stress, &pj8_info);
            }
            state.stresses[nl_eindex].qp[qp] =
                pj8_sub(full_effective_stress, static_stress);
            model_error = fabs(pj8_compatibility_residual(pj8_state));
            c->fs[qp] = model_error;
            state.load_unload_n[nl_eindex].qv[qp] =
                (double)pj8_info.accepted_substeps;
            if (!ok || !isfinite(model_error) ||
                !pj8_finite_tensor(full_effective_stress)) {
                nl_record_error(state.error, 5, eindex, qp, model_error,
                                pj8_pressure(full_effective_stress));
            }
#elif defined(NL_GPU_MANZARI_DAFALIAS_ONLY)
            int32_t md_substeps = 0;
            nl_md_diagnostics_t *md_diagnostics =
                &md_diagnostics_shared[threadIdx.x].value;
            ok = nl_update_manzari_dafalias_scheme<MD_SCHEME,
                                                    MD_FE_SUBSTEPS>(
                c, &cfg, strain, strain_prev, static_stress,
                &state.stresses[nl_eindex].qp[qp],
                &state.pstrains2[nl_eindex].qp[qp],
                &state.alphastress2[nl_eindex].qp[qp],
                &state.sref[nl_eindex].qp[qp],
                &state.fabric[nl_eindex].qp[qp],
                &state.sv_max[nl_eindex].qv[qp],
                &state.ep2[nl_eindex].qv[qp],
                &state.psi_n[nl_eindex].qv[qp], &c->fs[qp],
                &state.kappa[nl_eindex].qv[qp], &md_substeps,
                md_diagnostics);
            state.load_unload_n[nl_eindex].qv[qp] = (double)md_substeps;
            model_error = c->fs[qp];
            /* Modified Euler commits only an admissible stress-corrected
             * state, so its residual remains a production failure criterion.
             * OpenSees scheme-5 Forward Euler intentionally performs no
             * stress correction: a finite residual is a reported drift, not
             * constitutive nonconvergence.  Keep the finite-state/integrator
             * failure checks for both schemes. */
            const int enforce_md_yield_tolerance =
                MD_SCHEME == NL_MD_INTEGRATION_MODIFIED_EULER;
            if (!ok || !isfinite(model_error) ||
                (enforce_md_yield_tolerance &&
                 model_error > cfg.error_tolerance)) {
                const tensor_t md_stress = md_scale(md_add(
                    state.stresses[nl_eindex].qp[qp], static_stress), -1.0);
                nl_record_md_error(
                    state.error, eindex, qp,
                    fmax(model_error, md_diagnostics->correction_residual),
                    fmax(state.kappa[nl_eindex].qv[qp],
                         md_diagnostics->local_error),
                    md_substeps, md_diagnostics, md_stress,
                    md_sub(strain, strain_prev));
            }
#else
            int32_t unused_substeps = 0;
            nl_md_diagnostics_t unused_diagnostics = {};
            ok = nl_update_material(
                c, &cfg, strain, strain_prev, static_stress,
                state.pstrains1[nl_eindex].qp[qp],
                state.alphastress1[nl_eindex].qp[qp],
                state.ep1[nl_eindex].qv[qp],
                &state.pstrains2[nl_eindex].qp[qp],
                &state.alphastress2[nl_eindex].qp[qp],
                &state.stresses[nl_eindex].qp[qp],
                &state.ep2[nl_eindex].qv[qp], &c->fs[qp],
                &state.psi_n[nl_eindex].qv[qp],
                &state.load_unload_n[nl_eindex].qv[qp],
                &state.sv_n[nl_eindex].qv[qp],
                &state.sv_max[nl_eindex].qv[qp],
                &state.kappa[nl_eindex].qv[qp],
                &state.sref[nl_eindex].qp[qp],
                &state.gamma1d[nl_eindex].qv[qp],
                &state.tao1d[nl_eindex].qv[qp],
                &state.ggmax[nl_eindex].qv[qp], &model_error, NULL,
                &state.sv_max[nl_eindex].qv[qp],
                &state.kappa[nl_eindex].qv[qp], &unused_substeps,
                &unused_diagnostics);
            if (MATERIAL_MODEL_IS_BACKBONE(cfg.material_model)) {
                c->fs[qp] = model_error;
                if (!ok || !isfinite(model_error) ||
                    model_error > cfg.error_tolerance) {
                    nl_record_error(state.error, 2, eindex, qp,
                                    model_error,
                                    state.stresses[nl_eindex].qp[qp].xx);
                }
            } else if (!ok || !isfinite(c->fs[qp])) {
                nl_record_error(state.error, 1, eindex, qp, c->fs[qp],
                                state.stresses[nl_eindex].qp[qp].xx);
            }
#endif
        }
    }

    fvector_t local_force[8] = {};
    if (c->isTopoNonlin)
        nl_tetra_internal_force(state.stresses + nl_eindex, c, element.h,
                                local_force);
    else
        nl_hex_internal_force(state.stresses + nl_eindex, element.h,
                              local_force);

    fvector_t local_damping[8] = {};
    double damping_scale = 0.0;
    if (cfg.stiff_damp != 0.0 && cfg.stiff_damp_frequency != 0.0 &&
        delta_t > 0.0) {
        const double b_over_dt = cfg.stiff_damp /
            (cfg.stiff_damp_frequency*NL_GPU_PI*delta_t);
        int any_delta = 0;
        for (int node = 0; node < 8; node++) {
            const int32_t lnid = gpu->lnidArrayDevice[eindex + node*gpu->lenum];
            for (int d = 0; d < 3; d++) {
                du[node].f[d] = (solver_float)(
                    ((double)gpu->tm1Device[lnid].f[d] -
                     (double)gpu->tm2Device[lnid].f[d])*b_over_dt);
                any_delta |= (du[node].f[d] != (solver_float)0);
            }
        }
        if (any_delta) {
            nl_elastic_internal_force(du, c, element.h, local_damping);
            /* Preserve the legacy topographic damping convention. */
            damping_scale = c->isTopoNonlin ? 1.0 : delta_t_squared;
        }
    }

    for (int node = 0; node < 8; node++) {
        const int32_t lnid = gpu->lnidArrayDevice[eindex + node*gpu->lenum];
        for (int d = 0; d < 3; d++) {
            const solver_float value = (solver_float)(
                -(double)local_force[node].f[d]*delta_t_squared -
                 (double)local_damping[node].f[d]*damping_scale);
            nl_atomic_add(&gpu->forceDevice[lnid].f[d], value);
        }
    }
}

#if defined(NL_GPU_CLASSIC_ONLY)
extern "C" void nonlinear_gpu_launch_classic(
    nonlinear_gpu_state_t state, nonlinear_gpu_config_t config,
    gpu_data_t *gpu, int32_t grid_size, int32_t block_size,
    int32_t step, double delta_t, double delta_t_squared,
    cudaStream_t stream)
{
    kernel_nonlinear_update_force<<<grid_size, block_size, 0, stream>>>(
        state, config, gpu, step, delta_t, delta_t_squared);
}
#elif defined(NL_GPU_MANZARI_DAFALIAS_ONLY)
extern "C" void nonlinear_gpu_launch_manzari_dafalias(
    nonlinear_gpu_state_t state, nonlinear_gpu_config_t config,
    gpu_data_t *gpu, int32_t grid_size, int32_t block_size,
    int32_t step, double delta_t, double delta_t_squared,
    cudaStream_t stream)
{
    /* Dispatch to the specialization that matches the configured scheme and
     * Forward-Euler subdivision policy. Each is a separate ptxas entry
     * function, so neither Forward-Euler path inherits Modified Euler's
     * register allocation. */
    if (config.md_integration_scheme == NL_MD_INTEGRATION_FORWARD_EULER) {
        if (config.md_forward_euler_substeps ==
                NL_MD_FORWARD_EULER_SUBSTEPS_DEMAND) {
            kernel_nonlinear_update_force<
                NL_MD_INTEGRATION_FORWARD_EULER,
                NL_MD_FORWARD_EULER_SUBSTEPS_DEMAND>
                <<<grid_size, block_size, 0, stream>>>(
                    state, config, gpu, step, delta_t, delta_t_squared);
        } else {
            kernel_nonlinear_update_force<
                NL_MD_INTEGRATION_FORWARD_EULER,
                NL_MD_FORWARD_EULER_SUBSTEPS_FIXED>
                <<<grid_size, block_size, 0, stream>>>(
                    state, config, gpu, step, delta_t, delta_t_squared);
        }
    } else {
        kernel_nonlinear_update_force<
            NL_MD_INTEGRATION_MODIFIED_EULER,
            NL_MD_FORWARD_EULER_SUBSTEPS_FIXED>
            <<<grid_size, block_size, 0, stream>>>(
                state, config, gpu, step, delta_t, delta_t_squared);
    }
}
#elif defined(NL_GPU_PDMY03_ONLY)
extern "C" void nonlinear_gpu_launch_pdmy03(
    nonlinear_gpu_state_t state, nonlinear_gpu_config_t config,
    gpu_data_t *gpu, int32_t grid_size, int32_t block_size,
    int32_t step, double delta_t, double delta_t_squared,
    cudaStream_t stream)
{
    kernel_nonlinear_update_force<<<grid_size, block_size, 0, stream>>>(
        state, config, gpu, step, delta_t, delta_t_squared);
}
#elif defined(NL_GPU_PJ_LIQ_V8_ONLY)
extern "C" void nonlinear_gpu_launch_pj_liq_v8(
    nonlinear_gpu_state_t state, nonlinear_gpu_config_t config,
    gpu_data_t *gpu, int32_t grid_size, int32_t block_size,
    int32_t step, double delta_t, double delta_t_squared,
    cudaStream_t stream)
{
    kernel_nonlinear_update_force<<<grid_size, block_size, 0, stream>>>(
        state, config, gpu, step, delta_t, delta_t_squared);
}
#endif
#endif

#ifdef NL_GPU_RUNTIME_ONLY
#define NL_GPU_WORKLOAD_BLOCK_SIZE 256

__global__ static void kernel_nonlinear_workload_stats(
    const qpvectors_t *substeps, size_t point_count,
    nonlinear_gpu_workload_stats_t *stats)
{
    __shared__ unsigned long long zero_or_inactive[NL_GPU_WORKLOAD_BLOCK_SIZE];
    __shared__ unsigned long long one_step[NL_GPU_WORKLOAD_BLOCK_SIZE];
    __shared__ unsigned long long multi_step[NL_GPU_WORKLOAD_BLOCK_SIZE];
    __shared__ unsigned long long accepted[NL_GPU_WORKLOAD_BLOCK_SIZE];
    __shared__ unsigned long long maximum[NL_GPU_WORKLOAD_BLOCK_SIZE];

    const int lane = (int)threadIdx.x;
    const size_t first = (size_t)blockIdx.x*blockDim.x+threadIdx.x;
    const size_t stride = (size_t)blockDim.x*gridDim.x;
    const double *values = (const double *)substeps;
    unsigned long long local_zero = 0;
    unsigned long long local_one = 0;
    unsigned long long local_multi = 0;
    unsigned long long local_accepted = 0;
    unsigned long long local_maximum = 0;

    for (size_t point = first; point < point_count; point += stride) {
        const double value = values[point];
        const unsigned long long count = value > 0.0
            ? (unsigned long long)value : 0ULL;
        if (count == 0ULL)
            local_zero++;
        else if (count == 1ULL)
            local_one++;
        else
            local_multi++;
        local_accepted += count;
        if (count > local_maximum) local_maximum = count;
    }

    zero_or_inactive[lane] = local_zero;
    one_step[lane] = local_one;
    multi_step[lane] = local_multi;
    accepted[lane] = local_accepted;
    maximum[lane] = local_maximum;
    __syncthreads();

    for (int offset = NL_GPU_WORKLOAD_BLOCK_SIZE/2;
         offset > 0; offset /= 2) {
        if (lane < offset) {
            zero_or_inactive[lane] += zero_or_inactive[lane+offset];
            one_step[lane] += one_step[lane+offset];
            multi_step[lane] += multi_step[lane+offset];
            accepted[lane] += accepted[lane+offset];
            if (maximum[lane+offset] > maximum[lane])
                maximum[lane] = maximum[lane+offset];
        }
        __syncthreads();
    }

    if (lane == 0) {
        atomicAdd(&stats->zero_or_inactive_points, zero_or_inactive[0]);
        atomicAdd(&stats->one_step_points, one_step[0]);
        atomicAdd(&stats->multi_step_points, multi_step[0]);
        atomicAdd(&stats->accepted_substeps, accepted[0]);
        atomicMax(&stats->maximum_substeps, maximum[0]);
    }
}

__global__ static void kernel_nonlinear_gather_station_state(
    nonlinear_gpu_state_t state, int32_t station_count,
    const int32_t *station_indices,
    nonlinear_gpu_station_state_t *station_state)
{
    const int32_t station =
        (int32_t)(blockIdx.x*blockDim.x + threadIdx.x);
    if (station >= station_count) return;

    const int32_t nl_eindex = station_indices[station];
    const int qp = 4;
    nonlinear_gpu_station_state_t result = {};
    result.nonlinear_index = nl_eindex;
    result.stress = state.stresses[nl_eindex].qp[qp];
    result.strain = state.strains[nl_eindex].qp[qp];
    result.fs = state.constants[nl_eindex].fs[qp];
    result.kappa = state.kappa[nl_eindex].qv[qp];
    result.gamma1d = state.gamma1d[nl_eindex].qv[qp];
    result.tao1d = state.tao1d[nl_eindex].qv[qp];
    result.ggmax = state.ggmax[nl_eindex].qv[qp];
    result.sigma_z_static = state.constants[nl_eindex].sigmaZ_st;
    result.phi = state.constants[nl_eindex].phi;
    station_state[station] = result;
}

#undef NL_GPU_WORKLOAD_BLOCK_SIZE

static int nonlinear_cuda_failure(int32_t myID, cudaError_t error,
                                  const char *operation)
{
    if (error == cudaSuccess) return 0;
    fprintf(stderr, "Thread %d: nonlinear GPU %s failed: %s\n",
            myID, operation, cudaGetErrorString(error));
    return -1;
}

template <typename T>
static int nonlinear_alloc_copy(int32_t myID, T **device, const T *host,
                                size_t count, const char *name)
{
    *device = NULL;
    if (count == 0) return 0;
    cudaError_t error = cudaMalloc((void **)device, count*sizeof(T));
    if (nonlinear_cuda_failure(myID, error, name) != 0) return -1;
    error = cudaMemcpy(*device, host, count*sizeof(T), cudaMemcpyHostToDevice);
    if (nonlinear_cuda_failure(myID, error, name) != 0) return -1;
    return 0;
}

#define NL_ALLOC_COPY(field, host_field)                                      \
    do {                                                                       \
        if (nonlinear_alloc_copy(myID, &g_nl.field, host_state->host_field,    \
                                 count, #field) != 0) return -1;               \
    } while (0)

int nonlinear_gpu_initialize(int32_t myID, mesh_t *mesh, mysolver_t *solver,
                             int32_t count, const int32_t *mapping,
                             nlsolver_t *host_state,
                             const nonlinear_gpu_config_t *config)
{
    if (!nonlinear_gpu_requested()) return 0;
    if (mesh == NULL || solver == NULL || host_state == NULL || config == NULL ||
        count < 0 ||
        !md_valid_integration_scheme(config->md_integration_scheme) ||
        !md_valid_forward_euler_substeps(
            config->md_forward_euler_substeps) ||
        !md_valid_recovery_mode(config->md_recovery_mode) ||
        (config->md_integration_scheme !=
             NL_MD_INTEGRATION_FORWARD_EULER &&
         config->md_forward_euler_substeps !=
             NL_MD_FORWARD_EULER_SUBSTEPS_FIXED) ||
        (config->md_integration_scheme == NL_MD_INTEGRATION_FORWARD_EULER &&
         config->md_recovery_mode != NL_MD_RECOVERY_STRICT))
        return -1;

    const int local_stats_interval =
        nonlinear_gpu_workload_interval(myID);
    int minimum_stats_interval = 0;
    int maximum_stats_interval = 0;
    if (MPI_Allreduce(&local_stats_interval, &minimum_stats_interval, 1,
                      MPI_INT, MPI_MIN, comm_solver) != MPI_SUCCESS ||
        MPI_Allreduce(&local_stats_interval, &maximum_stats_interval, 1,
                      MPI_INT, MPI_MAX, comm_solver) != MPI_SUCCESS)
        return -1;
    if (minimum_stats_interval < 0) return -1;
    if (minimum_stats_interval != maximum_stats_interval) {
        if (myID == 0) {
            fprintf(stderr,
                    "HERCULES_NONLINEAR_STATS_INTERVAL must have the same "
                    "value on every solver rank (range=%d..%d)\n",
                    minimum_stats_interval, maximum_stats_interval);
        }
        return -1;
    }
    g_workload.interval = config->material_model == MANZARI_DAFALIAS
        ? maximum_stats_interval : 0;
    if (myID == 0 && maximum_stats_interval > 0 &&
        config->material_model != MANZARI_DAFALIAS) {
        fprintf(stderr,
                "HERCULES_NONLINEAR_STATS_INTERVAL is currently available "
                "only for Manzari-Dafalias; diagnostics are disabled\n");
    }

    nonlinear_gpu_element_t *elements = NULL;
    if (count > 0) {
        elements = (nonlinear_gpu_element_t *)malloc(
            (size_t)count*sizeof(nonlinear_gpu_element_t));
        if (elements == NULL) return -1;
        for (int32_t i = 0; i < count; i++) {
            const int32_t eindex = mapping[i];
            edata_t *edata = (edata_t *)mesh->elemTable[eindex].data;
            elements[i].eindex = eindex;
            elements[i].h = (double)edata->edgesize;
        }
    }

    g_nl.count = count;
    g_config = *config;
    if (nonlinear_alloc_copy(myID, &g_nl.elements, elements, count,
                             "element metadata") != 0) {
        free(elements);
        return -1;
    }
    free(elements);

    NL_ALLOC_COPY(constants, constants);
    NL_ALLOC_COPY(stresses, stresses);
    NL_ALLOC_COPY(strains, strains);
    const int allocate_previous_state_shadows =
        config->material_model != MANZARI_DAFALIAS &&
        config->material_model != PJ_LIQ_V8;
    if (allocate_previous_state_shadows) {
        NL_ALLOC_COPY(strains1, strains1);
        NL_ALLOC_COPY(pstrains1, pstrains1);
    }
    NL_ALLOC_COPY(pstrains2, pstrains2);
    if (allocate_previous_state_shadows)
        NL_ALLOC_COPY(alphastress1, alphastress1);
    NL_ALLOC_COPY(alphastress2, alphastress2);
    if (allocate_previous_state_shadows)
        NL_ALLOC_COPY(ep1, ep1);
    NL_ALLOC_COPY(ep2, ep2);
    NL_ALLOC_COPY(psi_n, psi_n);
    NL_ALLOC_COPY(load_unload_n, LoUnlo_n);
    NL_ALLOC_COPY(sv_n, Sv_n);
    NL_ALLOC_COPY(sv_max, Sv_max);
    NL_ALLOC_COPY(kappa, kappa);
    NL_ALLOC_COPY(sref, Sref);
    if (config->material_model == MANZARI_DAFALIAS) {
        NL_ALLOC_COPY(fabric, fabric);
    }
    NL_ALLOC_COPY(gamma1d, gamma1d);
    NL_ALLOC_COPY(tao1d, tao1d);
    NL_ALLOC_COPY(ggmax, GGmax);

    if (config->material_model == PRESSURE_DEPEND_MULTIYIELD_03 && count > 0) {
        int32_t stride = 0;
        for (int32_t i = 0; i < count; i++) {
            if (host_state->constants[i].pdmy03.surface_count > stride)
                stride = host_state->constants[i].pdmy03.surface_count;
        }
        if (stride < 2 || stride > PDMY03_MAX_SURFACES) return -1;
        g_nl.pdmy03_surface_stride = stride;
        const size_t point_count = (size_t)count*NL_GPU_QP;
        cudaError_t pdmy_error = cudaMalloc(
            (void **)&g_nl.pdmy03_states,
            point_count*sizeof(*g_nl.pdmy03_states));
        if (pdmy_error == cudaSuccess)
            pdmy_error = cudaMemset(g_nl.pdmy03_states, 0,
                point_count*sizeof(*g_nl.pdmy03_states));
        if (pdmy_error == cudaSuccess)
            pdmy_error = cudaMalloc(
                (void **)&g_nl.pdmy03_centers,
                point_count*(size_t)stride*sizeof(*g_nl.pdmy03_centers));
        if (pdmy_error == cudaSuccess)
            pdmy_error = cudaMemset(g_nl.pdmy03_centers, 0,
                point_count*(size_t)stride*sizeof(*g_nl.pdmy03_centers));
        if (nonlinear_cuda_failure(myID, pdmy_error,
                "PDMY03 history allocation") != 0) return -1;
    }

    if (config->material_model == PJ_LIQ_V8 && count > 0) {
        const size_t point_count = (size_t)count*NL_GPU_QP;
        const pj8_parameters_t parameters = pj8_reference_parameters(1000.0);
        cudaError_t pj8_error = cudaMalloc(
            &g_nl.pj8_parameters, sizeof(parameters));
        if (pj8_error == cudaSuccess)
            pj8_error = cudaMemcpy(g_nl.pj8_parameters, &parameters,
                sizeof(parameters), cudaMemcpyHostToDevice);
        if (pj8_error == cudaSuccess)
            pj8_error = cudaMalloc(&g_nl.pj8_states,
                point_count*sizeof(pj8_state_t));
        if (pj8_error == cudaSuccess)
            pj8_error = cudaMemset(g_nl.pj8_states, 0,
                point_count*sizeof(pj8_state_t));
        if (nonlinear_cuda_failure(myID, pj8_error,
                "PJ-Liq V8 history allocation") != 0) return -1;
    }

    if (nonlinear_cuda_failure(myID,
            cudaMalloc((void **)&g_nl.error, sizeof(nonlinear_gpu_error_t)),
            "error buffer allocation") != 0 ||
        nonlinear_cuda_failure(myID,
            cudaMemset(g_nl.error, 0, sizeof(nonlinear_gpu_error_t)),
            "error buffer initialization") != 0) return -1;
    if (g_workload.interval > 0 && count > 0) {
        if (nonlinear_cuda_failure(myID,
                cudaMalloc((void **)&g_workload.device,
                           sizeof(*g_workload.device)),
                "workload-statistics buffer allocation") != 0 ||
            nonlinear_cuda_failure(myID,
                cudaMemset(g_workload.device, 0,
                           sizeof(*g_workload.device)),
                "workload-statistics buffer initialization") != 0)
            return -1;
    }

    g_mesh_host = mesh;
    g_active = 1;
    if (myID == 0) {
        fprintf(stderr,
                "Nonlinear backend: GPU (default backend, double-precision "
                "constitutive state, %d local elements on rank 0)\n", count);
        if (!allocate_previous_state_shadows) {
            const double saved_mib = (double)((size_t)count*
                (3*sizeof(qptensors_t)+sizeof(qpvectors_t)))/(1024.0*1024.0);
            fprintf(stderr,
                    "Nonlinear GPU specialized state: omitted unused legacy "
                    "previous-state shadows (%.1f MiB on rank 0)\n",
                    saved_mib);
        }
        if (config->material_model == MANZARI_DAFALIAS) {
            if (config->md_integration_scheme ==
                    NL_MD_INTEGRATION_FORWARD_EULER) {
                if (config->md_forward_euler_substeps ==
                        NL_MD_FORWARD_EULER_SUBSTEPS_DEMAND) {
                    fprintf(stderr,
                            "Manzari-Dafalias integration: Forward Euler "
                            "(substeps=demand, max_substeps=%d, "
                            "yield_drift=reported, geostatic=%s)\n",
                            config->no_substeps,
                            config->approximate_geostatic_state == YES
                                ? "approximate K0 effective stress"
                                : "solver stress history");
                } else {
                    fprintf(stderr,
                            "Manzari-Dafalias integration: Forward Euler "
                            "(substeps=fixed, fixed_substeps=%d over the "
                            "plastic remainder, yield_drift=reported, "
                            "geostatic=%s)\n",
                            config->no_substeps,
                            config->approximate_geostatic_state == YES
                                ? "approximate K0 effective stress"
                                : "solver stress history");
                }
            } else {
                const char *recovery = "strict";
                if (config->md_recovery_mode == NL_MD_RECOVERY_FORCE_ACCEPT)
                    recovery = "force_accept";
                else if (config->md_recovery_mode ==
                         NL_MD_RECOVERY_UNLIMITED_SUBSTEPS)
                    recovery = "unlimited_substeps";
                else if (config->md_recovery_mode == NL_MD_RECOVERY_OPENSEES)
                    recovery = "opensees";
                if (md_unlimited_substeps_enabled(config->md_recovery_mode)) {
                    fprintf(stderr,
                        "Manzari-Dafalias integration: Modified Euler "
                        "(local_error=%g, substeps=unlimited, "
                        "no_substeps=%d correction cap, "
                        "yield_residual=%g Pa, recovery=%s, "
                        "geostatic=%s)\n",
                        config->backbone_error_tolerance,
                        config->no_substeps, config->error_tolerance,
                        recovery,
                        config->approximate_geostatic_state == YES
                            ? "approximate K0 effective stress"
                            : "solver stress history");
                } else {
                    fprintf(stderr,
                        "Manzari-Dafalias integration: Modified Euler "
                        "(local_error=%g, max_substeps=%d, "
                        "yield_residual=%g Pa, recovery=%s, "
                        "geostatic=%s)\n",
                        config->backbone_error_tolerance,
                        config->no_substeps, config->error_tolerance,
                        recovery,
                        config->approximate_geostatic_state == YES
                            ? "approximate K0 effective stress"
                            : "solver stress history");
                }
                if (md_force_accept_enabled(config->md_recovery_mode)) {
                    fprintf(stderr,
                        "Manzari-Dafalias warning: exhausted adaptive "
                        "recovery may force-accept a finite admissible "
                        "state above backbone_errTol\n");
                }
            }
            if (g_workload.interval > 0) {
                fprintf(stderr,
                        "Manzari-Dafalias GPU workload statistics: every "
                        "%d time step%s (opt-in device scan and global MPI "
                        "reduction)\n",
                        g_workload.interval,
                        g_workload.interval == 1 ? "" : "s");
            }
        } else if (config->material_model == PRESSURE_DEPEND_MULTIYIELD_03) {
            const double state_mib = (double)((size_t)count*NL_GPU_QP*
                (sizeof(pdmy03_state_t)+(size_t)g_nl.pdmy03_surface_stride*
                 sizeof(tensor_t)))/(1024.0*1024.0);
            fprintf(stderr,
                    "PressureDependMultiYield03 integration: standard "
                    "automatic nested surfaces (max_substeps=%d, "
                    "surface_stride=%d, local_history=%.1f MiB, "
                    "geostatic=%s)\n",
                    config->no_substeps, g_nl.pdmy03_surface_stride,
                    state_mib,
                    config->approximate_geostatic_state == YES
                        ? "approximate K0 effective stress"
                        : "solver stress history");
        } else if (config->material_model == PJ_LIQ_V8) {
            const double state_mib = (double)((size_t)count*NL_GPU_QP*
                sizeof(pj8_state_t))/(1024.0*1024.0);
            fprintf(stderr,
                    "PJ-Liq V8 integration: frozen calibration %s "
                    "(fixed_substeps=%d, local_history=%.1f MiB, "
                    "effective-skeleton stress with BrickUP pore pressure "
                    "and Biot total-stress coupling, geostatic=%s)\n",
                    PJ8_PARAMETER_SHA256, config->no_substeps, state_mib,
                    config->approximate_geostatic_state == YES
                        ? "approximate K0 effective stress"
                        : "solver stress history");
        }
        fflush(stderr);
    }
    return 0;
}

#undef NL_ALLOC_COPY

void nonlinear_gpu_enqueue(mysolver_t *solver, int32_t step,
                           double delta_t, double delta_t_squared)
{
    if (!g_active) return;
    const int sample_workload = g_workload.interval > 0 &&
        (step % g_workload.interval) == 0;
    if (sample_workload) {
        g_workload.pending = 1;
        g_workload.step = step;
        memset(&g_workload.host, 0, sizeof(g_workload.host));
    }
    if (g_nl.count == 0) return;

    cudaStream_t stream = solver->streams[CUDA_STREAM_MAIN];
    cudaError_t error = cudaMemsetAsync(
        g_nl.error, 0, sizeof(nonlinear_gpu_error_t), stream);
    if (sample_workload && error == cudaSuccess) {
        error = cudaMemsetAsync(
            g_nl.load_unload_n, 0,
            (size_t)g_nl.count*sizeof(*g_nl.load_unload_n), stream);
    }
    if (sample_workload && error == cudaSuccess) {
        error = cudaMemsetAsync(g_workload.device, 0,
                                sizeof(*g_workload.device), stream);
    }
    if (error != cudaSuccess) {
        fprintf(stderr, "Nonlinear GPU diagnostic setup failed: %s\n",
                cudaGetErrorString(error));
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
    }

    const int block_size = NL_GPU_UPDATE_BLOCK_SIZE;
    const int grid_size = (g_nl.count + block_size - 1)/block_size;
    if (g_config.material_model == MANZARI_DAFALIAS) {
        nonlinear_gpu_launch_manzari_dafalias(
            g_nl, g_config, solver->gpuDataDevice, grid_size, block_size,
            step, delta_t, delta_t_squared, stream);
    } else if (g_config.material_model == PRESSURE_DEPEND_MULTIYIELD_03) {
        nonlinear_gpu_launch_pdmy03(
            g_nl, g_config, solver->gpuDataDevice, grid_size, block_size,
            step, delta_t, delta_t_squared, stream);
    } else if (g_config.material_model == PJ_LIQ_V8) {
        nonlinear_gpu_launch_pj_liq_v8(
            g_nl, g_config, solver->gpuDataDevice, grid_size, block_size,
            step, delta_t, delta_t_squared, stream);
    } else {
        nonlinear_gpu_launch_classic(
            g_nl, g_config, solver->gpuDataDevice, grid_size, block_size,
            step, delta_t, delta_t_squared, stream);
    }
    error = cudaGetLastError();
    if (error != cudaSuccess) {
        fprintf(stderr, "Nonlinear GPU kernel launch failed: %s\n",
                cudaGetErrorString(error));
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
    }

    if (sample_workload) {
        const size_t point_count = (size_t)g_nl.count*NL_GPU_QP;
        const int stats_block_size = 256;
        size_t stats_grid_size =
            (point_count+(size_t)stats_block_size-1)/stats_block_size;
        if (stats_grid_size > 1024) stats_grid_size = 1024;
        kernel_nonlinear_workload_stats<<<
            (int)stats_grid_size, stats_block_size, 0, stream>>>(
                g_nl.load_unload_n, point_count, g_workload.device);
        error = cudaGetLastError();
        if (error != cudaSuccess) {
            fprintf(stderr,
                    "Nonlinear GPU workload-statistics kernel launch "
                    "failed: %s\n", cudaGetErrorString(error));
            MPI_Abort(MPI_COMM_WORLD, ERROR);
            exit(1);
        }
    }
}

static void nonlinear_gpu_release_station_output(void)
{
    if (g_station.indices != NULL) cudaFree(g_station.indices);
    if (g_station.device != NULL) cudaFree(g_station.device);
    if (g_station.host != NULL) cudaFreeHost(g_station.host);
    memset(&g_station, 0, sizeof(g_station));
}

int nonlinear_gpu_configure_station_output(
    int32_t myID, int32_t station_count,
    const int32_t *station_nonlinear_indices)
{
    if (!g_active || station_count < 0 ||
        (station_count > 0 && station_nonlinear_indices == NULL)) return -1;

    nonlinear_gpu_release_station_output();
    if (station_count == 0) return 0;

    for (int32_t i = 0; i < station_count; i++) {
        if (station_nonlinear_indices[i] < 0 ||
            station_nonlinear_indices[i] >= g_nl.count) {
            fprintf(stderr,
                    "Thread %d: nonlinear GPU station %d has invalid "
                    "nonlinear index %d (count=%d)\n",
                    myID, i, station_nonlinear_indices[i], g_nl.count);
            return -1;
        }
    }

    g_station.count = station_count;
    cudaError_t error = cudaMalloc(
        (void **)&g_station.indices,
        (size_t)station_count*sizeof(*g_station.indices));
    if (nonlinear_cuda_failure(myID, error,
                               "station index allocation") != 0) goto fail;
    error = cudaMemcpy(g_station.indices, station_nonlinear_indices,
                       (size_t)station_count*sizeof(*g_station.indices),
                       cudaMemcpyHostToDevice);
    if (nonlinear_cuda_failure(myID, error,
                               "station index upload") != 0) goto fail;
    error = cudaMalloc(
        (void **)&g_station.device,
        (size_t)station_count*sizeof(*g_station.device));
    if (nonlinear_cuda_failure(myID, error,
                               "station state allocation") != 0) goto fail;
    error = cudaMallocHost(
        (void **)&g_station.host,
        (size_t)station_count*sizeof(*g_station.host));
    if (nonlinear_cuda_failure(myID, error,
                               "station host buffer allocation") != 0) {
        goto fail;
    }
    return 0;

fail:
    nonlinear_gpu_release_station_output();
    return -1;
}

void nonlinear_gpu_download_station_state(mysolver_t *solver,
                                           nlsolver_t *host_state)
{
    if (!g_active || g_station.count == 0 || host_state == NULL) return;

    cudaStream_t stream = solver->streams[CUDA_STREAM_MAIN];
    const int block_size = 128;
    const int grid_size = (g_station.count + block_size - 1)/block_size;
    kernel_nonlinear_gather_station_state<<<grid_size, block_size, 0, stream>>>(
        g_nl, g_station.count, g_station.indices, g_station.device);
    cudaError_t error = cudaGetLastError();
    if (error == cudaSuccess) {
        error = cudaMemcpyAsync(
            g_station.host, g_station.device,
            (size_t)g_station.count*sizeof(*g_station.host),
            cudaMemcpyDeviceToHost, stream);
    }
    if (error == cudaSuccess) error = cudaStreamSynchronize(stream);
    if (error != cudaSuccess) {
        fprintf(stderr, "Nonlinear GPU station-state gather failed: %s\n",
                cudaGetErrorString(error));
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
    }

    for (int32_t i = 0; i < g_station.count; i++) {
        const nonlinear_gpu_station_state_t *record = g_station.host + i;
        const int32_t nl_eindex = record->nonlinear_index;
        host_state->stresses[nl_eindex].qp[4] = record->stress;
        host_state->strains[nl_eindex].qp[4] = record->strain;
        host_state->constants[nl_eindex].fs[4] = record->fs;
        host_state->kappa[nl_eindex].qv[4] = record->kappa;
        host_state->gamma1d[nl_eindex].qv[4] = record->gamma1d;
        host_state->tao1d[nl_eindex].qv[4] = record->tao1d;
        host_state->GGmax[nl_eindex].qv[4] = record->ggmax;
        host_state->constants[nl_eindex].sigmaZ_st = record->sigma_z_static;
        host_state->constants[nl_eindex].phi = record->phi;
    }
}

#define NL_DOWNLOAD(field, host_field)                                        \
    cudaMemcpy(host_state->host_field, g_nl.field,                            \
               (size_t)g_nl.count*sizeof(*host_state->host_field),            \
               cudaMemcpyDeviceToHost)

void nonlinear_gpu_download_state(mysolver_t *solver, nlsolver_t *host_state)
{
    if (!g_active || g_nl.count == 0 || host_state == NULL) return;
    cudaStreamSynchronize(solver->streams[CUDA_STREAM_MAIN]);
    cudaError_t error = cudaSuccess;
    error = NL_DOWNLOAD(constants, constants);
    if (error == cudaSuccess) error = NL_DOWNLOAD(stresses, stresses);
    if (error == cudaSuccess) error = NL_DOWNLOAD(strains, strains);
    if (error == cudaSuccess && g_nl.strains1 != NULL)
        error = NL_DOWNLOAD(strains1, strains1);
    if (error == cudaSuccess && g_nl.pstrains1 != NULL)
        error = NL_DOWNLOAD(pstrains1, pstrains1);
    if (error == cudaSuccess) error = NL_DOWNLOAD(pstrains2, pstrains2);
    if (error == cudaSuccess && g_nl.alphastress1 != NULL)
        error = NL_DOWNLOAD(alphastress1, alphastress1);
    if (error == cudaSuccess) error = NL_DOWNLOAD(alphastress2, alphastress2);
    if (error == cudaSuccess && g_nl.ep1 != NULL)
        error = NL_DOWNLOAD(ep1, ep1);
    if (error == cudaSuccess) error = NL_DOWNLOAD(ep2, ep2);
    if (error == cudaSuccess) error = NL_DOWNLOAD(psi_n, psi_n);
    if (error == cudaSuccess) error = NL_DOWNLOAD(load_unload_n, LoUnlo_n);
    if (error == cudaSuccess) error = NL_DOWNLOAD(sv_n, Sv_n);
    if (error == cudaSuccess) error = NL_DOWNLOAD(sv_max, Sv_max);
    if (error == cudaSuccess) error = NL_DOWNLOAD(kappa, kappa);
    if (error == cudaSuccess) error = NL_DOWNLOAD(sref, Sref);
    if (error == cudaSuccess && g_config.material_model == MANZARI_DAFALIAS)
        error = NL_DOWNLOAD(fabric, fabric);
    if (error == cudaSuccess) error = NL_DOWNLOAD(gamma1d, gamma1d);
    if (error == cudaSuccess) error = NL_DOWNLOAD(tao1d, tao1d);
    if (error == cudaSuccess) error = NL_DOWNLOAD(ggmax, GGmax);
    if (error != cudaSuccess) {
        fprintf(stderr, "Nonlinear GPU state download failed: %s\n",
                cudaGetErrorString(error));
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
    }
}

#undef NL_DOWNLOAD

static void nonlinear_gpu_report_workload(int32_t myID)
{
    if (!g_workload.pending) return;

    if (g_nl.count > 0) {
        const cudaError_t error = cudaMemcpy(
            &g_workload.host, g_workload.device,
            sizeof(g_workload.host), cudaMemcpyDeviceToHost);
        if (error != cudaSuccess) {
            fprintf(stderr,
                    "Thread %d: nonlinear GPU workload-statistics download "
                    "failed: %s\n", myID, cudaGetErrorString(error));
            MPI_Abort(MPI_COMM_WORLD, ERROR);
            exit(1);
        }
    }

    const unsigned long long local_sums[4] = {
        g_workload.host.zero_or_inactive_points,
        g_workload.host.one_step_points,
        g_workload.host.multi_step_points,
        g_workload.host.accepted_substeps
    };
    unsigned long long global_sums[4] = {};
    const unsigned long long local_maximum =
        g_workload.host.maximum_substeps;
    unsigned long long global_maximum = 0;
    const int sum_status = MPI_Reduce(
        local_sums, global_sums, 4, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0,
        comm_solver);
    const int max_status = MPI_Reduce(
        &local_maximum, &global_maximum, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX,
        0, comm_solver);
    if (sum_status != MPI_SUCCESS || max_status != MPI_SUCCESS) {
        fprintf(stderr,
                "Thread %d: nonlinear GPU workload-statistics MPI "
                "reduction failed\n", myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
    }

    if (myID == 0) {
        const unsigned long long updated_points =
            global_sums[1]+global_sums[2];
        const unsigned long long point_slots =
            global_sums[0]+updated_points;
        const double multi_fraction = updated_points > 0
            ? 100.0*(double)global_sums[2]/(double)updated_points : 0.0;
        const double mean_substeps = updated_points > 0
            ? (double)global_sums[3]/(double)updated_points : 0.0;
        fprintf(stderr,
                "Nonlinear GPU workload: step=%d point_slots=%llu "
                "zero_or_inactive=%llu updated=%llu one_step=%llu "
                "multi_step=%llu multi_step_fraction=%.2f%% "
                "accepted_substeps=%llu mean_substeps=%.3f "
                "maximum_substeps=%llu\n",
                g_workload.step, point_slots, global_sums[0],
                updated_points, global_sums[1], global_sums[2],
                multi_fraction, global_sums[3], mean_substeps,
                global_maximum);
        fflush(stderr);
    }
    g_workload.pending = 0;
}

void nonlinear_gpu_check_error(int32_t myID)
{
    if (!g_active) return;
    if (g_nl.count == 0) {
        nonlinear_gpu_report_workload(myID);
        return;
    }
    const cudaError_t error = cudaMemcpy(&g_error_host, g_nl.error,
        sizeof(g_error_host), cudaMemcpyDeviceToHost);
    if (error != cudaSuccess) {
        fprintf(stderr, "Thread %d: nonlinear GPU error check failed: %s\n",
                myID, cudaGetErrorString(error));
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
    }
    if (g_error_host.code != 0) {
        fprintf(stderr,
                "Thread %d: nonlinear GPU constitutive failure "
                "(code=%d element=%d qp=%d residual=%12.8E value=%12.8E)\n",
                myID, g_error_host.code, g_error_host.eindex,
                g_error_host.qp, g_error_host.residual, g_error_host.value);
        if (g_error_host.code == 3) {
            const char *reason = "unknown";
            switch (g_error_host.md_reason) {
            case NL_MD_FAILURE_INITIAL_PRESSURE:
                reason = "non-positive initial effective mean stress"; break;
            case NL_MD_FAILURE_INVALID_CONTROLS:
                reason = "invalid integration controls or void ratio"; break;
            case NL_MD_FAILURE_ELASTIC_MODULI:
                reason = "invalid pressure-dependent elastic moduli"; break;
            case NL_MD_FAILURE_YIELD_INTERSECTION:
                reason = "elastic/yield intersection did not converge"; break;
            case NL_MD_FAILURE_ATTEMPT_LIMIT:
                reason = "maximum total constitutive attempts reached"; break;
            case NL_MD_FAILURE_FIRST_STAGE:
                reason = g_config.md_integration_scheme ==
                        NL_MD_INTEGRATION_FORWARD_EULER
                    ? "Forward-Euler stage is invalid"
                    : "first Modified-Euler stage is invalid";
                break;
            case NL_MD_FAILURE_TENSILE_PREDICTOR:
                reason = "predictor reached the minimum effective pressure";
                break;
            case NL_MD_FAILURE_SECOND_STAGE:
                reason = "second Modified-Euler stage is invalid"; break;
            case NL_MD_FAILURE_NONFINITE_LOCAL_ERROR:
                reason = "non-finite local integration error"; break;
            case NL_MD_FAILURE_LOCAL_ERROR_MIN_STEP:
                reason = "local-error target failed at minimum substep";
                break;
            case NL_MD_FAILURE_SUBSTEP_LIMIT:
                reason = "maximum accepted constitutive substeps reached";
                break;
            case NL_MD_FAILURE_STRESS_CORRECTION_MIN_STEP:
                reason = "yield correction failed at minimum substep"; break;
            case NL_MD_FAILURE_FINAL_STATE:
                reason = "final constitutive state is inadmissible"; break;
            default:
                break;
            }
            fprintf(stderr,
                    "  Manzari-Dafalias reason: %s (reason=%d)\n"
                    "  attempts=%d accepted_substeps=%d "
                    "correction_retries=%d max_correction_iterations=%d "
                    "segment_projections=%d force_accepts=%d\n"
                    "  progress_fraction=%12.8E substep_fraction=%12.8E "
                    "local_error=%12.8E "
                    "yield_residual=%12.8E effective_mean_stress=%12.8E Pa "
                    "strain_increment_norm=%12.8E\n",
                    reason, g_error_host.md_reason,
                    g_error_host.md_attempts,
                    g_error_host.md_accepted_substeps,
                    g_error_host.md_correction_retries,
                    g_error_host.md_correction_iterations,
                    g_error_host.md_projection_fallbacks,
                    g_error_host.md_force_accepts,
                    g_error_host.md_progress_fraction,
                    g_error_host.md_substep_fraction, g_error_host.value,
                    g_error_host.residual,
                    g_error_host.md_effective_mean_stress,
                    g_error_host.md_strain_increment_norm);
            fprintf(stderr,
                    "  %s stage%s: "
                    "dgamma=(%12.8E, %12.8E) "
                    "raw_dgamma=(%12.8E, %12.8E) "
                    "denominator=(%12.8E, %12.8E)\n"
                    "  stage h=(%12.8E, %12.8E) "
                    "hardening_gap=(%12.8E, %12.8E)\n",
                    g_config.md_integration_scheme ==
                            NL_MD_INTEGRATION_FORWARD_EULER
                        ? "Forward-Euler" : "Modified-Euler",
                    g_config.md_integration_scheme ==
                            NL_MD_INTEGRATION_FORWARD_EULER
                        ? "" : "s",
                    g_error_host.md_first_dgamma,
                    g_error_host.md_second_dgamma,
                    g_error_host.md_first_raw_dgamma,
                    g_error_host.md_second_raw_dgamma,
                    g_error_host.md_first_denominator,
                    g_error_host.md_second_denominator,
                    g_error_host.md_first_h,
                    g_error_host.md_second_h,
                    g_error_host.md_first_hardening_gap,
                    g_error_host.md_second_hardening_gap);
            fprintf(stderr,
                    "  pre-step effective stress (compression positive): "
                    "xx=%12.8E yy=%12.8E zz=%12.8E xy=%12.8E "
                    "yz=%12.8E xz=%12.8E Pa\n"
                    "  Hercules strain increment: "
                    "xx=%12.8E yy=%12.8E zz=%12.8E xy=%12.8E "
                    "yz=%12.8E xz=%12.8E\n",
                    g_error_host.md_stress.xx,
                    g_error_host.md_stress.yy,
                    g_error_host.md_stress.zz,
                    g_error_host.md_stress.xy,
                    g_error_host.md_stress.yz,
                    g_error_host.md_stress.xz,
                    g_error_host.md_strain_increment.xx,
                    g_error_host.md_strain_increment.yy,
                    g_error_host.md_strain_increment.zz,
                    g_error_host.md_strain_increment.xy,
                    g_error_host.md_strain_increment.yz,
                    g_error_host.md_strain_increment.xz);
            if (g_mesh_host != NULL && g_error_host.eindex >= 0 &&
                g_error_host.eindex < g_mesh_host->lenum) {
                const elem_t *element =
                    &g_mesh_host->elemTable[g_error_host.eindex];
                const node_t *node =
                    &g_mesh_host->nodeTable[element->lnid[0]];
                const edata_t *edata = (const edata_t *)element->data;
                const double h = edata != NULL ? edata->edgesize : 0.0;
                const double x = node->x*g_mesh_host->ticksize;
                const double y = node->y*g_mesh_host->ticksize;
                const double z = node->z*g_mesh_host->ticksize;
                fprintf(stderr,
                        "  element origin=(%.6f, %.6f, %.6f) m "
                        "center=(%.6f, %.6f, %.6f) m h=%.6f m\n",
                        x, y, z, x+0.5*h, y+0.5*h, z+0.5*h, h);
            }
        } else if (g_error_host.code == 4) {
            const char *reason = "unknown";
            switch (g_error_host.md_reason) {
            case PDMY03_FAILURE_CONSTANTS:
                reason = "invalid constants or substep control"; break;
            case PDMY03_FAILURE_INITIAL_GEOMETRY:
                reason = "initial nested-surface geometry is invalid"; break;
            case PDMY03_FAILURE_NONFINITE:
                reason = "non-finite constitutive state"; break;
            case PDMY03_FAILURE_SURFACE_GEOMETRY:
                reason = "nested-surface geometry update failed"; break;
            }
            fprintf(stderr,
                    "  PDMY03 reason: %s (reason=%d)\n"
                    "  requested_substeps=%d accepted_substeps=%d "
                    "crossed_surfaces=%d substeps_clamped=%d\n"
                    "  normalized_yield_residual=%12.8E "
                    "effective_mean_stress=%12.8E Pa "
                    "strain_increment_norm=%12.8E\n",
                    reason, g_error_host.md_reason,
                    g_error_host.md_attempts,
                    g_error_host.md_accepted_substeps,
                    g_error_host.md_projection_fallbacks,
                    g_error_host.md_force_accepts,
                    g_error_host.residual,
                    g_error_host.md_effective_mean_stress,
                    g_error_host.md_strain_increment_norm);
        } else if (g_error_host.code == 5) {
            fprintf(stderr,
                    "  PJ-Liq V8 reason: invalid/non-finite explicit state "
                    "or initialization; compatibility_residual=%12.8E, "
                    "effective_mean_stress=%12.8E Pa\n",
                    g_error_host.residual, g_error_host.value);
        }
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
    }
    nonlinear_gpu_report_workload(myID);
}
#endif

#ifndef NL_GPU_PRODUCTION_ONLY
int nonlinear_gpu_material_point_history(
    const nlconstants_t *constants,
    const nonlinear_gpu_config_t *config,
    const nonlinear_material_point_initial_state_t *initial_state,
    const tensor_t *strain_history,
    int32_t steps,
    nonlinear_material_point_result_t *results)
{
    if (constants == NULL || config == NULL || initial_state == NULL ||
        strain_history == NULL ||
        results == NULL || steps <= 0) {
        return -1;
    }

    tensor_t *strain_device = NULL;
    nonlinear_material_point_result_t *result_device = NULL;
    cudaError_t error = cudaMalloc((void **)&strain_device,
                                   (size_t)steps*sizeof(*strain_device));
    if (error == cudaSuccess) {
        error = cudaMalloc((void **)&result_device,
                           (size_t)steps*sizeof(*result_device));
    }
    if (error == cudaSuccess) {
        error = cudaMemcpy(strain_device, strain_history,
                           (size_t)steps*sizeof(*strain_device),
                           cudaMemcpyHostToDevice);
    }
    if (error == cudaSuccess) {
        kernel_nonlinear_material_point_history<<<1, 1>>>(
            *constants, *config, *initial_state, strain_device, steps,
            result_device);
        error = cudaGetLastError();
    }
    if (error == cudaSuccess) {
        error = cudaDeviceSynchronize();
    }
    if (error == cudaSuccess) {
        error = cudaMemcpy(results, result_device,
                           (size_t)steps*sizeof(*result_device),
                           cudaMemcpyDeviceToHost);
    }

    if (result_device != NULL) cudaFree(result_device);
    if (strain_device != NULL) cudaFree(strain_device);
    if (error != cudaSuccess) {
        fprintf(stderr, "Nonlinear GPU material-point test failed: %s\n",
                cudaGetErrorString(error));
        return -1;
    }
    return 0;
}

int nonlinear_manzari_dafalias_host_history(
    const nlconstants_t *constants,
    const nonlinear_gpu_config_t *config,
    const nonlinear_material_point_initial_state_t *initial_state,
    const tensor_t *strain_history,
    int32_t steps,
    nonlinear_material_point_result_t *results)
{
    if (constants == NULL || config == NULL || initial_state == NULL ||
        strain_history == NULL ||
        results == NULL || steps <= 0 ||
        config->material_model != MANZARI_DAFALIAS) return -1;

    const int seeded = initial_state->provided != 0;
    tensor_t strain_prev = seeded ? initial_state->strain : md_zero();
    tensor_t stress = seeded ? initial_state->stress : md_zero();
    tensor_t plastic_strain = seeded ? initial_state->plastic_strain : md_zero();
    tensor_t alpha = seeded ? initial_state->alpha : md_zero();
    tensor_t alpha_in = seeded ? initial_state->alpha_in : md_zero();
    tensor_t fabric = seeded ? initial_state->fabric : md_zero();
    double initialized = seeded ? 1.0 : 0.0;

    for (int32_t step = 0; step < steps; step++) {
        double void_ratio = constants->md.e_init;
        double dgamma = 0.0;
        double yield_error = 0.0;
        double integration_error = 0.0;
        int32_t accepted_substeps = 0;
        nl_md_diagnostics_t diagnostics = {};
        const int ok = nl_update_manzari_dafalias(
            constants, config, strain_history[step], strain_prev, md_zero(),
            &stress, &plastic_strain, &alpha, &alpha_in, &fabric,
            &initialized, &void_ratio, &dgamma, &yield_error,
            &integration_error, &accepted_substeps, &diagnostics);

        nonlinear_material_point_result_t result = {};
        result.stress = stress;
        result.plastic_strain = plastic_strain;
        result.backstress = alpha;
        result.ep_bar = void_ratio;
        result.fs = yield_error;
        result.psi = dgamma;
        result.load_unload = (double)accepted_substeps;
        result.stress_reference = alpha_in;
        result.model_error = yield_error;
        result.fabric = fabric;
        result.alpha_in = alpha_in;
        result.void_ratio = void_ratio;
        result.plastic_multiplier = dgamma;
        result.integration_error = integration_error;
        result.accepted_substeps = accepted_substeps;
        result.status = ok ? 0 : 1;
        results[step] = result;
        strain_prev = strain_history[step];
    }
    return 0;
}

int nonlinear_pdmy03_host_history(
    const nlconstants_t *constants,
    const nonlinear_gpu_config_t *config,
    const nonlinear_material_point_initial_state_t *initial_state,
    const tensor_t *strain_history,
    int32_t steps,
    nonlinear_material_point_result_t *results)
{
    if (constants == NULL || config == NULL || initial_state == NULL ||
        strain_history == NULL || results == NULL || steps <= 0 ||
        config->material_model != PRESSURE_DEPEND_MULTIYIELD_03 ||
        constants->pdmy03.surface_count < 2 ||
        constants->pdmy03.surface_count > PDMY03_MAX_SURFACES) return -1;

    const int seeded = initial_state->provided != 0;
    tensor_t strain_previous = seeded ? initial_state->strain : pdmy03_zero();
    tensor_t stress = seeded ? initial_state->stress : pdmy03_zero();
    tensor_t plastic_strain = seeded ? initial_state->plastic_strain
                                     : pdmy03_zero();
    pdmy03_state_t state = {};
    tensor_t centers[PDMY03_MAX_SURFACES] = {};

    for (int32_t step = 0; step < steps; step++) {
        pdmy03_diagnostics_t diagnostics = {};
        const int ok = pdmy03_update(
            &constants->pdmy03, config->no_substeps, strain_history[step],
            strain_previous, pdmy03_zero(), &stress, &plastic_strain,
            &state, centers, &diagnostics);
        nonlinear_material_point_result_t result = {};
        result.stress = stress;
        result.plastic_strain = plastic_strain;
        result.fs = diagnostics.yield_residual;
        result.model_error = diagnostics.yield_residual;
        result.accepted_substeps = diagnostics.accepted_substeps;
        result.pdmy03_active_surface = state.active_surface;
        result.pdmy03_requested_substeps = diagnostics.requested_substeps;
        result.pdmy03_crossed_surfaces = diagnostics.crossed_surfaces;
        result.pdmy03_substeps_clamped = diagnostics.substeps_clamped;
        result.pdmy03_pressure_d = state.pressure_d;
        double center_norm = 0.0;
        for (int i = 0; i < constants->pdmy03.surface_count; i++)
            center_norm += pdmy03_ddot(centers[i], centers[i]);
        result.pdmy03_center_norm = sqrt(center_norm);
        result.status = ok ? 0 : 1;
        results[step] = result;
        strain_previous = strain_history[step];
    }
    return 0;
}

int nonlinear_pj8_host_history(
    const nlconstants_t *constants,
    const nonlinear_gpu_config_t *config,
    const nonlinear_material_point_initial_state_t *initial_state,
    const tensor_t *strain_history,
    int32_t steps,
    nonlinear_material_point_result_t *results)
{
    if (constants == NULL || config == NULL || initial_state == NULL ||
        strain_history == NULL || results == NULL || steps <= 0 ||
        config->material_model != PJ_LIQ_V8 ||
        !initial_state->provided || config->no_substeps < 1) return -1;

    const pj8_parameters_t parameters = pj8_reference_parameters(1000.0);
    tensor_t strain_previous = initial_state->strain;
    pj8_state_t state = {};
    tensor_t reference = initial_state->stress;
    reference.xy = reference.yz = reference.xz = 0.0;
    if (!pj8_initialize(&parameters, initial_state->stress,
                        constants->pj8_void_ratio, &state) ||
        !pj8_begin_dynamic_phase(&parameters, reference, &state)) return -1;

    for (int32_t step = 0; step < steps; step++) {
        pj8_update_info_t info = {};
        tensor_t stress = state.stress;
        const int ok = pj8_update(&parameters,
            pj8_sub(strain_history[step], strain_previous),
            config->no_substeps, &state, &stress, &info);
        nonlinear_material_point_result_t result = {};
        result.stress = stress;
        result.backstress = state.alpha;
        result.ep_bar = state.lambda_total;
        result.fs = fabs(pj8_compatibility_residual(&state));
        result.psi = state.D;
        result.load_unload = (double)info.accepted_substeps;
        result.kappa = pj8_pressure(state.stress);
        result.stress_reference = state.alpha01;
        result.gamma1d = state.eps_v_irreversible;
        result.tao1d = state.eps_v_reversible;
        result.ggmax = state.amplitude_factor;
        result.model_error = result.fs;
        result.fabric = state.fabric;
        result.alpha_in = state.alpha0;
        result.void_ratio = state.void_ratio;
        result.plastic_multiplier = state.lambda_total;
        result.accepted_substeps = info.accepted_substeps;
        pj8_state_values(&state, result.pj8_state_values);
        result.pj8_reversal_registered = info.reversal_registered;
        result.status = ok ? 0 : 1;
        results[step] = result;
        strain_previous = strain_history[step];
    }
    return 0;
}
#endif

#ifdef NL_GPU_RUNTIME_ONLY
enum {
    NL_CHECKPOINT_BRICKUP_PRESSURE = 1u,
    NL_CHECKPOINT_PJ8_STATE = 2u
};

typedef struct nonlinear_checkpoint_header_t {
    char magic[8];
    uint32_t version;
    uint32_t flags;
    uint32_t rank_count;
    uint32_t pj8_state_bytes;
    uint64_t data_bytes;
    char pj8_parameter_sha256[65];
    char reserved[7];
} nonlinear_checkpoint_header_t;

static const char g_nonlinear_checkpoint_magic[8] = {
    'H', 'N', 'L', 'V', '8', 'C', 'P', '\0'
};

static size_t nonlinear_checkpoint_rank_bytes(
    uint32_t flags, uint64_t nodes, uint64_t points)
{
    size_t bytes = 0;
    if ((flags & NL_CHECKPOINT_BRICKUP_PRESSURE) != 0)
        bytes += (size_t)(2*nodes)*sizeof(solver_float);
    if ((flags & NL_CHECKPOINT_PJ8_STATE) != 0) {
        const uint64_t elements = points/NL_GPU_QP;
        bytes += (size_t)(2*elements)*sizeof(qptensors_t);
        bytes += (size_t)points*sizeof(pj8_state_t);
    }
    return bytes;
}

static int nonlinear_checkpoint_write_device(
    FILE *fp, const void *device, size_t bytes, void *scratch,
    size_t scratch_bytes)
{
    const unsigned char *source = (const unsigned char *)device;
    while (bytes > 0) {
        const size_t chunk = bytes < scratch_bytes ? bytes : scratch_bytes;
        if (cudaMemcpy(scratch, source, chunk, cudaMemcpyDeviceToHost) !=
                cudaSuccess || fwrite(scratch, 1, chunk, fp) != chunk)
            return -1;
        source += chunk;
        bytes -= chunk;
    }
    return 0;
}

static int nonlinear_checkpoint_read_device(
    FILE *fp, void *device, size_t bytes, void *scratch,
    size_t scratch_bytes)
{
    unsigned char *target = (unsigned char *)device;
    while (bytes > 0) {
        const size_t chunk = bytes < scratch_bytes ? bytes : scratch_bytes;
        if (fread(scratch, 1, chunk, fp) != chunk ||
            cudaMemcpy(target, scratch, chunk, cudaMemcpyHostToDevice) !=
                cudaSuccess)
            return -1;
        target += chunk;
        bytes -= chunk;
    }
    return 0;
}

static uint32_t nonlinear_checkpoint_local_flags(mysolver_t *solver)
{
    uint32_t flags = 0;
    if (solver != NULL && solver->pore_pressure_tm1 != NULL &&
        solver->pore_pressure_tm2 != NULL && solver->gpuData != NULL &&
        solver->gpuData->porePressureTm1Device != NULL &&
        solver->gpuData->porePressureTm2Device != NULL)
        flags |= NL_CHECKPOINT_BRICKUP_PRESSURE;
    /* The field schema must agree on ranks that happen to own zero nonlinear
     * elements; their point count is zero and no device bytes are copied. */
    if (g_active && g_config.material_model == PJ_LIQ_V8)
        flags |= NL_CHECKPOINT_PJ8_STATE;
    return flags;
}

int nonlinear_gpu_checkpoint_write(
    const char *filename, off_t extension_offset, int32_t myID,
    int32_t group_size, int32_t node_count, mysolver_t *solver,
    MPI_Comm comm_solver)
{
    const uint32_t local_flags = nonlinear_checkpoint_local_flags(solver);
    uint32_t minimum_flags = 0, maximum_flags = 0;
    MPI_Allreduce(&local_flags, &minimum_flags, 1, MPI_UNSIGNED, MPI_MIN,
                  comm_solver);
    MPI_Allreduce(&local_flags, &maximum_flags, 1, MPI_UNSIGNED, MPI_MAX,
                  comm_solver);
    if (minimum_flags != maximum_flags) {
        if (myID == 0)
            fprintf(stderr, "Checkpoint extension fields differ by rank\n");
        return -1;
    }
    if (local_flags == 0) return 0;

    const uint64_t local_nodes =
        (local_flags & NL_CHECKPOINT_BRICKUP_PRESSURE) != 0
            ? (uint64_t)node_count : 0;
    const uint64_t local_points =
        (local_flags & NL_CHECKPOINT_PJ8_STATE) != 0
            ? (uint64_t)g_nl.count*NL_GPU_QP : 0;
    uint64_t *counts = (uint64_t *)malloc(
        (size_t)2*group_size*sizeof(uint64_t));
    int local_error = counts == NULL;
    int global_error = 0;
    MPI_Allreduce(&local_error, &global_error, 1, MPI_INT, MPI_MAX,
                  comm_solver);
    if (global_error) {
        free(counts);
        return -1;
    }
    MPI_Allgather(&local_nodes, 1, MPI_UNSIGNED_LONG_LONG,
                  counts, 1, MPI_UNSIGNED_LONG_LONG, comm_solver);
    MPI_Allgather(&local_points, 1, MPI_UNSIGNED_LONG_LONG,
                  counts+group_size, 1, MPI_UNSIGNED_LONG_LONG,
                  comm_solver);

    nonlinear_checkpoint_header_t header = {};
    memcpy(header.magic, g_nonlinear_checkpoint_magic, sizeof(header.magic));
    header.version = 1;
    header.flags = local_flags;
    header.rank_count = (uint32_t)group_size;
    header.pj8_state_bytes = (uint32_t)sizeof(pj8_state_t);
    if ((local_flags & NL_CHECKPOINT_PJ8_STATE) != 0)
        memcpy(header.pj8_parameter_sha256, PJ8_PARAMETER_SHA256,
               sizeof(PJ8_PARAMETER_SHA256));
    for (int rank = 0; rank < group_size; rank++)
        header.data_bytes += nonlinear_checkpoint_rank_bytes(
            local_flags, counts[rank], counts[group_size+rank]);

    if (myID == 0) {
        FILE *fp = fopen(filename, "rb+");
        if (fp == NULL || fseeko(fp, extension_offset, SEEK_SET) != 0 ||
            fwrite(&header, sizeof(header), 1, fp) != 1 ||
            fwrite(counts, sizeof(uint64_t), (size_t)2*group_size, fp) !=
                (size_t)2*group_size) {
            local_error = 1;
        }
        if (fp != NULL) fclose(fp);
    }
    MPI_Bcast(&local_error, 1, MPI_INT, 0, comm_solver);
    if (local_error) {
        free(counts);
        return -1;
    }
    MPI_Barrier(comm_solver);

    const size_t scratch_bytes = 64u*1024u*1024u;
    void *scratch = malloc(scratch_bytes);
    local_error = scratch == NULL;
    MPI_Allreduce(&local_error, &global_error, 1, MPI_INT, MPI_MAX,
                  comm_solver);
    if (global_error) {
        free(scratch);
        free(counts);
        return -1;
    }

    off_t record_offset = extension_offset+(off_t)sizeof(header)+
        (off_t)2*group_size*sizeof(uint64_t);
    for (int rank = 0; rank < myID; rank++)
        record_offset += (off_t)nonlinear_checkpoint_rank_bytes(
            local_flags, counts[rank], counts[group_size+rank]);
    for (int group_start = 0; group_start < group_size;
         group_start += 512) {
        if (myID >= group_start && myID < group_start+512) {
            FILE *fp = fopen(filename, "rb+");
            if (fp == NULL || fseeko(fp, record_offset, SEEK_SET) != 0) {
                local_error = 1;
            } else {
                if ((local_flags & NL_CHECKPOINT_BRICKUP_PRESSURE) != 0 &&
                    (nonlinear_checkpoint_write_device(
                         fp, solver->gpuData->porePressureTm1Device,
                         (size_t)local_nodes*sizeof(solver_float), scratch,
                         scratch_bytes) != 0 ||
                     nonlinear_checkpoint_write_device(
                         fp, solver->gpuData->porePressureTm2Device,
                         (size_t)local_nodes*sizeof(solver_float), scratch,
                         scratch_bytes) != 0))
                    local_error = 1;
                if (!local_error &&
                    (local_flags & NL_CHECKPOINT_PJ8_STATE) != 0) {
                    const size_t element_bytes =
                        (size_t)g_nl.count*sizeof(qptensors_t);
                    if (nonlinear_checkpoint_write_device(
                            fp, g_nl.strains, element_bytes, scratch,
                            scratch_bytes) != 0 ||
                        nonlinear_checkpoint_write_device(
                            fp, g_nl.stresses, element_bytes, scratch,
                            scratch_bytes) != 0 ||
                        nonlinear_checkpoint_write_device(
                            fp, g_nl.pj8_states,
                            (size_t)local_points*sizeof(pj8_state_t), scratch,
                            scratch_bytes) != 0)
                        local_error = 1;
                }
            }
            if (fp != NULL) fclose(fp);
        }
        MPI_Barrier(comm_solver);
    }
    MPI_Allreduce(&local_error, &global_error, 1, MPI_INT, MPI_MAX,
                  comm_solver);
    free(scratch);
    free(counts);
    if (!global_error && myID == 0) {
        fprintf(stderr,
                "Checkpoint extension: saved%s%s state (%llu bytes)\n",
                (local_flags & NL_CHECKPOINT_BRICKUP_PRESSURE) != 0
                    ? " BrickUP pressure" : "",
                (local_flags & NL_CHECKPOINT_PJ8_STATE) != 0
                    ? " PJ-Liq V8" : "",
                (unsigned long long)header.data_bytes);
    }
    return global_error ? -1 : 0;
}

int nonlinear_gpu_checkpoint_read(
    const char *filename, off_t extension_offset, int32_t myID,
    int32_t group_size, int32_t node_count, mysolver_t *solver,
    MPI_Comm comm_solver)
{
    const uint32_t local_flags = nonlinear_checkpoint_local_flags(solver);
    uint32_t minimum_flags = 0, maximum_flags = 0;
    MPI_Allreduce(&local_flags, &minimum_flags, 1, MPI_UNSIGNED, MPI_MIN,
                  comm_solver);
    MPI_Allreduce(&local_flags, &maximum_flags, 1, MPI_UNSIGNED, MPI_MAX,
                  comm_solver);
    if (minimum_flags != maximum_flags) return -1;
    if (local_flags == 0) return 0;

    nonlinear_checkpoint_header_t header = {};
    int local_error = 0;
    if (myID == 0) {
        FILE *fp = fopen(filename, "rb");
        if (fp == NULL || fseeko(fp, extension_offset, SEEK_SET) != 0 ||
            fread(&header, sizeof(header), 1, fp) != 1)
            local_error = 1;
        if (fp != NULL) fclose(fp);
    }
    MPI_Bcast(&local_error, 1, MPI_INT, 0, comm_solver);
    MPI_Bcast(&header, sizeof(header), MPI_BYTE, 0, comm_solver);
    if (local_error ||
        memcmp(header.magic, g_nonlinear_checkpoint_magic,
               sizeof(header.magic)) != 0 ||
        header.version != 1 || header.flags != local_flags ||
        header.rank_count != (uint32_t)group_size ||
        ((local_flags & NL_CHECKPOINT_PJ8_STATE) != 0 &&
         (header.pj8_state_bytes != sizeof(pj8_state_t) ||
          strcmp(header.pj8_parameter_sha256, PJ8_PARAMETER_SHA256) != 0))) {
        if (myID == 0)
            fprintf(stderr,
                    "Checkpoint lacks a compatible BrickUP/PJ-Liq V8 "
                    "extension; a coupled restart cannot reset pressure or "
                    "constitutive history\n");
        return -1;
    }

    uint64_t *counts = (uint64_t *)malloc(
        (size_t)2*group_size*sizeof(uint64_t));
    local_error = counts == NULL;
    int global_error = 0;
    MPI_Allreduce(&local_error, &global_error, 1, MPI_INT, MPI_MAX,
                  comm_solver);
    if (global_error) {
        free(counts);
        return -1;
    }
    if (myID == 0) {
        FILE *fp = fopen(filename, "rb");
        if (fp == NULL ||
            fseeko(fp, extension_offset+(off_t)sizeof(header), SEEK_SET) != 0 ||
            fread(counts, sizeof(uint64_t), (size_t)2*group_size, fp) !=
                (size_t)2*group_size)
            local_error = 1;
        if (fp != NULL) fclose(fp);
    }
    MPI_Bcast(&local_error, 1, MPI_INT, 0, comm_solver);
    MPI_Bcast(counts, 2*group_size, MPI_UNSIGNED_LONG_LONG, 0, comm_solver);
    const uint64_t local_nodes =
        (local_flags & NL_CHECKPOINT_BRICKUP_PRESSURE) != 0
            ? (uint64_t)node_count : 0;
    const uint64_t local_points =
        (local_flags & NL_CHECKPOINT_PJ8_STATE) != 0
            ? (uint64_t)g_nl.count*NL_GPU_QP : 0;
    uint64_t expected_data_bytes = 0;
    for (int rank = 0; rank < group_size; rank++)
        expected_data_bytes += nonlinear_checkpoint_rank_bytes(
            local_flags, counts[rank], counts[group_size+rank]);
    if (local_error || header.data_bytes != expected_data_bytes ||
        counts[myID] != local_nodes ||
        counts[group_size+myID] != local_points) {
        free(counts);
        return -1;
    }

    const size_t scratch_bytes = 64u*1024u*1024u;
    void *scratch = malloc(scratch_bytes);
    local_error = scratch == NULL;
    MPI_Allreduce(&local_error, &global_error, 1, MPI_INT, MPI_MAX,
                  comm_solver);
    if (global_error) {
        free(scratch);
        free(counts);
        return -1;
    }
    off_t record_offset = extension_offset+(off_t)sizeof(header)+
        (off_t)2*group_size*sizeof(uint64_t);
    for (int rank = 0; rank < myID; rank++)
        record_offset += (off_t)nonlinear_checkpoint_rank_bytes(
            local_flags, counts[rank], counts[group_size+rank]);
    for (int group_start = 0; group_start < group_size;
         group_start += 512) {
        if (myID >= group_start && myID < group_start+512) {
            FILE *fp = fopen(filename, "rb");
            if (fp == NULL || fseeko(fp, record_offset, SEEK_SET) != 0) {
                local_error = 1;
            } else {
                if ((local_flags & NL_CHECKPOINT_BRICKUP_PRESSURE) != 0 &&
                    (nonlinear_checkpoint_read_device(
                         fp, solver->gpuData->porePressureTm1Device,
                         (size_t)local_nodes*sizeof(solver_float), scratch,
                         scratch_bytes) != 0 ||
                     nonlinear_checkpoint_read_device(
                         fp, solver->gpuData->porePressureTm2Device,
                         (size_t)local_nodes*sizeof(solver_float), scratch,
                         scratch_bytes) != 0))
                    local_error = 1;
                if (!local_error &&
                    (local_flags & NL_CHECKPOINT_PJ8_STATE) != 0) {
                    const size_t element_bytes =
                        (size_t)g_nl.count*sizeof(qptensors_t);
                    if (nonlinear_checkpoint_read_device(
                            fp, g_nl.strains, element_bytes, scratch,
                            scratch_bytes) != 0 ||
                        nonlinear_checkpoint_read_device(
                            fp, g_nl.stresses, element_bytes, scratch,
                            scratch_bytes) != 0 ||
                        nonlinear_checkpoint_read_device(
                            fp, g_nl.pj8_states,
                            (size_t)local_points*sizeof(pj8_state_t), scratch,
                            scratch_bytes) != 0)
                        local_error = 1;
                }
            }
            if (fp != NULL) fclose(fp);
            if (!local_error &&
                (local_flags & NL_CHECKPOINT_BRICKUP_PRESSURE) != 0 &&
                (cudaMemcpy(solver->pore_pressure_tm1,
                            solver->gpuData->porePressureTm1Device,
                            (size_t)local_nodes*sizeof(solver_float),
                            cudaMemcpyDeviceToHost) != cudaSuccess ||
                 cudaMemcpy(solver->pore_pressure_tm2,
                            solver->gpuData->porePressureTm2Device,
                            (size_t)local_nodes*sizeof(solver_float),
                            cudaMemcpyDeviceToHost) != cudaSuccess))
                local_error = 1;
        }
        MPI_Barrier(comm_solver);
    }
    MPI_Allreduce(&local_error, &global_error, 1, MPI_INT, MPI_MAX,
                  comm_solver);
    free(scratch);
    free(counts);
    if (!global_error && myID == 0)
        fprintf(stderr,
                "Checkpoint extension: restored BrickUP/PJ-Liq V8 "
                "coupled state\n");
    return global_error ? -1 : 0;
}

void nonlinear_gpu_finalize(void)
{
    nonlinear_gpu_release_station_output();
    if (g_workload.device != NULL) cudaFree(g_workload.device);
    memset(&g_workload, 0, sizeof(g_workload));
#define NL_FREE(field) do { if (g_nl.field != NULL) cudaFree(g_nl.field); } while (0)
    NL_FREE(elements); NL_FREE(constants);
    NL_FREE(stresses); NL_FREE(strains); NL_FREE(strains1);
    NL_FREE(pstrains1); NL_FREE(pstrains2);
    NL_FREE(alphastress1); NL_FREE(alphastress2);
    NL_FREE(ep1); NL_FREE(ep2); NL_FREE(psi_n); NL_FREE(load_unload_n);
    NL_FREE(sv_n); NL_FREE(sv_max); NL_FREE(kappa); NL_FREE(sref);
    NL_FREE(fabric);
    NL_FREE(pdmy03_states); NL_FREE(pdmy03_centers);
    NL_FREE(pj8_parameters); NL_FREE(pj8_states);
    NL_FREE(gamma1d); NL_FREE(tao1d); NL_FREE(ggmax); NL_FREE(error);
#undef NL_FREE
    memset(&g_nl, 0, sizeof(g_nl));
    g_mesh_host = NULL;
    g_active = 0;
}
#endif
