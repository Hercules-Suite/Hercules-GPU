/* -*- C++ -*- */
/*
 * Device-resident nonlinear constitutive update and force assembly.
 *
 * The legacy nonlinear.cc implementation remains the CPU/GSL reference.  This
 * translation unit intentionally contains no GSL calls and performs no dynamic
 * allocation from device code.
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

#define NL_GPU_QP 8
#define NL_GPU_KAPPA_MIN 1.0e-10
#define NL_GPU_PI 3.141592653589793238462643383279502884

#if defined(NL_GPU_PRODUCTION_ONLY) && defined(NL_GPU_MATERIAL_POINT_ONLY)
#error "nonlinear GPU compilation modes are mutually exclusive"
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

typedef struct nonlinear_gpu_state_t {
    int32_t count;
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
    nonlinear_gpu_error_t *error;
} nonlinear_gpu_state_t;

static nonlinear_gpu_state_t g_nl = {};
static nonlinear_gpu_station_buffers_t g_station = {};
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
    double denominator;
    double h;
    double hardening_gap;
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
    out.denominator = denominator;
    out.h = dep.h;
    out.hardening_gap = dep.hardening_gap;

    if (fabs(denominator) <= DBL_MIN) return out;
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
    tensor_t *alpha, tensor_t fabric, tensor_t alpha_in,
    double void_ratio, tensor_t direction_hint, tensor_t *plastic_strain,
    double K, double G, double tolerance, int maximum_iterations,
    double *residual, int32_t *iterations_used, int32_t *projection_used)
{
    *iterations_used = 0;
    *projection_used = 0;
    const tensor_t proposed_stress = *stress;
    const tensor_t proposed_alpha = *alpha;
    const tensor_t proposed_plastic_strain = *plastic_strain;
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
            m, *stress, *alpha, fabric, void_ratio, alpha_in);
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

NL_HD static int nl_update_manzari_dafalias(
    const nlconstants_t *c, const nonlinear_gpu_config_t *cfg,
    tensor_t strain, tensor_t strain_prev, tensor_t static_stress,
    tensor_t *dynamic_stress, tensor_t *plastic_strain, tensor_t *alpha,
    tensor_t *alpha_in, tensor_t *fabric, double *initialized,
    double *void_ratio, double *plastic_multiplier, double *yield_error,
    double *integration_error, int32_t *accepted_substeps,
    nl_md_diagnostics_t *diagnostics)
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
    nl_md_diagnostics_t empty_diagnostics = {};
    *diagnostics = empty_diagnostics;
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
    const tensor_t elastic_trial = md_elastic_increment(de, K, G);
    if (md_ddot(md_sub(alpha_work, alpha_in_work), elastic_trial) < 0.0)
        alpha_in_work = alpha_work;

    const tensor_t trial = md_add(stress, elastic_trial);
    const double trial_p = md_trace(trial)/3.0+1.0e-2*m->P_atm;
    const double trial_f = md_yield(m, trial, alpha_work);
    if (trial_p >= 1.0e-2*m->P_atm && trial_f <= cfg->error_tolerance) {
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
    /* OpenSees stops reducing the normalized constitutive fraction at 1e-6
     * and then accepts a special first-order state even when its embedded
     * error estimate is too large. Hercules keeps the configured error target
     * strict. Permit up to three additional decades when the imposed strain
     * increment is large enough to represent them accurately, while retaining
     * the 1e-6 floor for already tiny global increments. */
    const double strain_roundoff =
        64.0*DBL_EPSILON*(1.0+md_norm(strain_prev));
    const double dT_min = fmin(
        1.0e-6,
        fmax(1.0e-9, strain_roundoff/fmax(md_norm(de), DBL_MIN)));
    /* no_substeps caps accepted work. Rejected error-control trials are
     * bounded separately so one difficult point cannot spin forever. */
    int attempts = 0;
    while (T < 1.0-DBL_EPSILON) {
        if (++attempts > 4*cfg->no_substeps) {
            diagnostics->failure_reason = NL_MD_FAILURE_ATTEMPT_LIMIT;
            diagnostics->attempts = attempts;
            return 0;
        }
        diagnostics->attempts = attempts;
        if (*accepted_substeps >= cfg->no_substeps) {
            diagnostics->failure_reason = NL_MD_FAILURE_SUBSTEP_LIMIT;
            return 0;
        }
        dT = fmin(dT, 1.0-T);
        diagnostics->substep_fraction = dT;
        diagnostics->progress_fraction = T;
        diagnostics->effective_mean_stress = md_trace(stress)/3.0;
        const tensor_t de_step = md_scale(de, dT);
        const tensor_t strain_at_end = md_add(
            md_scale(strain_prev, -1.0), md_scale(de, T+dT));
        const double e_step = m->e_init-
            (1.0+m->e_init)*md_trace(strain_at_end);
        const nl_md_stage_t first = md_stage(
            m, stress, alpha_work, fabric_work, alpha_in_work,
            e_step, de_step, K, G);
        if (!first.ok) {
            diagnostics->failure_reason = NL_MD_FAILURE_FIRST_STAGE;
            return 0;
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
        const double predictor_p = md_trace(predictor_stress)/3.0+
                                   1.0e-2*m->P_atm;
        if (predictor_p < 1.0e-2*m->P_atm) {
            diagnostics->effective_mean_stress =
                md_trace(predictor_stress)/3.0;
            if (dT <= dT_min) {
                diagnostics->failure_reason = NL_MD_FAILURE_TENSILE_PREDICTOR;
                return 0;
            }
            dT = fmax(0.1*dT, dT_min);
            continue;
        }
        const double predictor_f = md_yield(
            m, predictor_stress, alpha_work);
        if (!isfinite(predictor_f)) {
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
        if (first.dgamma == 0.0 &&
            predictor_f <= cfg->error_tolerance) {
            stress = predictor_stress;
            diagnostics->local_error = 0.0;
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
            diagnostics->failure_reason = NL_MD_FAILURE_SECOND_STAGE;
            return 0;
        }
        diagnostics->second_dgamma = second.dgamma;
        diagnostics->second_raw_dgamma = second.raw_dgamma;
        diagnostics->second_denominator = second.denominator;
        diagnostics->second_h = second.h;
        diagnostics->second_hardening_gap = second.hardening_gap;
        const double stress_norm = md_norm(stress);
        const double error_denominator = stress_norm < 0.5
            ? 1.0 : 2.0*stress_norm;
        const double step_error =
            md_norm(md_sub(second.dsigma, first.dsigma))/error_denominator;
        diagnostics->local_error = step_error;
        if (!isfinite(step_error)) {
            diagnostics->failure_reason = NL_MD_FAILURE_NONFINITE_LOCAL_ERROR;
            return 0;
        }
        if (step_error > cfg->backbone_error_tolerance) {
            if (dT <= dT_min) {
                *integration_error = step_error;
                diagnostics->failure_reason =
                    NL_MD_FAILURE_LOCAL_ERROR_MIN_STEP;
                return 0;
            }
            const double q = fmax(0.8*sqrt(
                cfg->backbone_error_tolerance/step_error), 0.1);
            dT = fmax(q*dT, dT_min);
            continue;
        }

        /* Build the entire candidate in local temporaries. A failed yield
         * correction must not modify committed history; reject it, reduce the
         * constitutive fraction, and retry from the same accepted state. */
        tensor_t candidate_stress = md_add(stress, md_scale(
            md_add(first.dsigma, second.dsigma), 0.5));
        tensor_t candidate_alpha = md_add(alpha_work, md_scale(
            md_add(first.dalpha, second.dalpha), 0.5));
        tensor_t candidate_fabric = md_add(fabric_work, md_scale(
            md_add(first.dfabric, second.dfabric), 0.5));
        tensor_t candidate_pstrain = md_add(pstrain, md_scale(
            md_add(first.dpstrain, second.dpstrain), 0.5));
        const double candidate_multiplier = *plastic_multiplier+
            0.5*(first.dgamma+second.dgamma);

        double correction_residual = 0.0;
        int32_t correction_iterations = 0;
        int32_t projection_used = 0;
        if (!md_stress_correction(
                m, stress, &candidate_stress, &candidate_alpha,
                candidate_fabric, alpha_in_work, e_step,
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
        diagnostics->projection_fallbacks += projection_used;

        stress = candidate_stress;
        alpha_work = candidate_alpha;
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
            *sigma = nl_spectral_recompose(corrected, vec);
            *pstrain = nl_sub(strain,
                              nl_elastic_strain(*sigma, c->mu, bulk));
            *epbar = 0.0;
            *fs = corrected[0];
            return 1;
        }
        if (shear_tension > 0.0 && eig[0] <= 0.0) {
            *sigma = sigma_trial;
            *epbar = epbar_prev;
            *fs = eig[0];
            return 1;
        }
    }
    double yield = eig[0] - eig[2] + (eig[0] + eig[2])*sin(c->phi) -
                   2.0*(c->c + c->h*epbar_prev)*cos(c->phi);
    if (yield < 0.0) {
        *sigma = sigma_trial;
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
    *sigma = nl_spectral_recompose(corrected, vec);
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
    tensor_t relative_dev = nl_sub(nl_dev(stress_predictor), backstress_prev);
    const double j2 = nl_j2(relative_dev);
    const double i1 = nl_i1(stress_predictor);
    double hardening;
    if (model == VONMISES_EP) hardening = c->c * 2.0/sqrt(3.0);
    else if (model == DRUCKERPRAGER)
        hardening = c->gamma*(c->c+c->h*epbar_prev);
    else if (model == VONMISES_FA) hardening = c->Sstrain0*c->mu;
    else if (model == MOHR_COULOMB) {
        double eig[3], vec[3][3];
        nl_jacobi_eigen(stress_predictor, eig, vec);
        hardening = 2.0*(c->c+c->h*epbar_prev)*cos(c->phi);
        const double pred = eig[0]-eig[2]+(eig[0]+eig[2])*sin(c->phi)-hardening;
        if (pred < 1.0e-3) {
            *pstrain = pstrain_prev; *sigma = stress_predictor;
            *epbar = epbar_prev; *fs = pred; return 1;
        }
        return nl_update_mohr_coulomb(c, strain, stress_predictor, pstrain,
                                      sigma, epbar_prev, epbar, fs,
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
    tensor_t updated_dev = nl_dev(*sigma);
    double updated_j2 = nl_j2(updated_dev);
    *fs = ((model == DRUCKERPRAGER) ? c->alpha*nl_i1(*sigma) : 0.0) +
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
        *sigma = nl_iso(skk/3.0);
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

    const int seeded = config.material_model == MANZARI_DAFALIAS &&
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

    for (int32_t step = 0; step < steps; step++) {
        tensor_t plastic_next = plastic_strain;
        tensor_t backstress_next = backstress;
        double ep_bar_next = ep_bar;
        double fs = 0.0;
        double model_error = 0.0;
        nl_md_diagnostics_t md_diagnostics = {};

        const int ok = nl_update_material(
            &constants, &config, strain_history[step], strain_prev, nl_zero(),
            plastic_strain, backstress, ep_bar,
            &plastic_next, &backstress_next, &stress, &ep_bar_next, &fs,
            &psi, &load_unload, &tao, &tao_max, &kappa,
            &stress_reference, &gamma1d, &tao1d, &ggmax, &model_error,
            &fabric, &md_initialized, &md_integration_error, &md_substeps,
            &md_diagnostics);

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

#ifndef NL_GPU_MATERIAL_POINT_ONLY
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

__device__ static tensor_t nl_md_static_stress(
    const nlconstants_t *c, const nonlinear_gpu_config_t *cfg,
    double h, int qp)
{
    if (cfg->approximate_geostatic_state != YES || c->sigmaZ_st <= 0.0)
        return nl_zero();
    const double lz = c->isTopoNonlin ? 0.0 :
        ((qp & 4) ? 0.577350269189 : -0.577350269189);
    const double vertical = nl_max(
        c->sigmaZ_st+9.8*c->md.density*h*lz*0.5, 0.0);
    const double k0 = 1.0-sin(c->phi);
    tensor_t stress = {-k0*vertical, -k0*vertical, -vertical,
                       0.0, 0.0, 0.0};
    return stress;
}

__global__ static void kernel_nonlinear_update_force(
    nonlinear_gpu_state_t state, nonlinear_gpu_config_t cfg,
    gpu_data_t *gpu, int32_t step, double delta_t, double delta_t_squared)
{
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

    if (step == 0 && MATERIAL_MODEL_IS_BACKBONE(cfg.material_model)) {
        for (int qp = 0; qp < NL_GPU_QP; qp++) {
            state.kappa[nl_eindex].qv[qp] = 100.0;
            state.gamma1d[nl_eindex].qv[qp] = 1.0e-10;
            state.ggmax[nl_eindex].qv[qp] = 1.0;
        }
    }

    if (any_displacement) {
        const int ngp = c->isTopoNonlin ? 5 : 8;
        for (int qp = 0; qp < ngp; qp++) {
            if (c->isTopoNonlin && c->tetraVol[qp] == 0.0) continue;

            const tensor_t strain_prev = state.strains[nl_eindex].qp[qp];
            state.strains1[nl_eindex].qp[qp] = strain_prev;
            const tensor_t strain = c->isTopoNonlin
                ? nl_tetra_strain(u, c->topoPart, qp, element.h)
                : nl_hex_strain(u, qp, element.h);
            state.strains[nl_eindex].qp[qp] = strain;
            state.pstrains1[nl_eindex].qp[qp] =
                state.pstrains2[nl_eindex].qp[qp];
            state.alphastress1[nl_eindex].qp[qp] =
                state.alphastress2[nl_eindex].qp[qp];
            state.ep1[nl_eindex].qv[qp] = state.ep2[nl_eindex].qv[qp];

            if (cfg.material_model == LINEAR ||
                step <= cfg.geostatic_final_step) {
                state.stresses[nl_eindex].qp[qp] =
                    nl_stress(strain, c->mu, c->lambda);
                continue;
            }

            double model_error = 0.0;
            int32_t md_substeps = 0;
            nl_md_diagnostics_t md_diagnostics = {};
            const tensor_t static_stress =
                (cfg.material_model == MANZARI_DAFALIAS)
                ? nl_md_static_stress(c, &cfg, element.h, qp) : nl_zero();
            const int ok = nl_update_material(
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
                &state.ggmax[nl_eindex].qv[qp], &model_error,
                (cfg.material_model == MANZARI_DAFALIAS)
                    ? &state.fabric[nl_eindex].qp[qp] : NULL,
                &state.sv_max[nl_eindex].qv[qp],
                &state.kappa[nl_eindex].qv[qp], &md_substeps,
                &md_diagnostics);

            if (MATERIAL_MODEL_IS_BACKBONE(cfg.material_model)) {
                c->fs[qp] = model_error;
                if (!ok || !isfinite(model_error) ||
                    model_error > cfg.error_tolerance) {
                    nl_record_error(state.error, 2, eindex, qp,
                                    model_error,
                                    state.stresses[nl_eindex].qp[qp].xx);
                }
            } else if (cfg.material_model == MANZARI_DAFALIAS) {
                c->fs[qp] = model_error;
                if (!ok || !isfinite(model_error) ||
                    model_error > cfg.error_tolerance) {
                    const tensor_t md_stress = md_scale(md_add(
                        state.stresses[nl_eindex].qp[qp], static_stress),
                        -1.0);
                    nl_record_md_error(
                        state.error, eindex, qp,
                        fmax(model_error,
                             md_diagnostics.correction_residual),
                        fmax(state.kappa[nl_eindex].qv[qp],
                             md_diagnostics.local_error),
                        md_substeps,
                        &md_diagnostics, md_stress,
                        md_sub(strain, strain_prev));
                }
            } else if (!ok || !isfinite(c->fs[qp])) {
                nl_record_error(state.error, 1, eindex, qp, c->fs[qp],
                                state.stresses[nl_eindex].qp[qp].xx);
            }
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
        count < 0) return -1;

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
    NL_ALLOC_COPY(strains1, strains1);
    NL_ALLOC_COPY(pstrains1, pstrains1);
    NL_ALLOC_COPY(pstrains2, pstrains2);
    NL_ALLOC_COPY(alphastress1, alphastress1);
    NL_ALLOC_COPY(alphastress2, alphastress2);
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

    if (nonlinear_cuda_failure(myID,
            cudaMalloc((void **)&g_nl.error, sizeof(nonlinear_gpu_error_t)),
            "error buffer allocation") != 0 ||
        nonlinear_cuda_failure(myID,
            cudaMemset(g_nl.error, 0, sizeof(nonlinear_gpu_error_t)),
            "error buffer initialization") != 0) return -1;

    g_mesh_host = mesh;
    g_active = 1;
    if (myID == 0) {
        fprintf(stderr,
                "Nonlinear backend: GPU (default backend, double-precision "
                "constitutive state, %d local elements on rank 0)\n", count);
        if (config->material_model == MANZARI_DAFALIAS) {
            fprintf(stderr,
                    "Manzari-Dafalias integration: Modified Euler "
                    "(local_error=%g, max_substeps=%d, "
                    "yield_residual=%g Pa, "
                    "geostatic=%s)\n",
                    config->backbone_error_tolerance, config->no_substeps,
                    config->error_tolerance,
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
    if (!g_active || g_nl.count == 0) return;
    cudaStream_t stream = solver->streams[CUDA_STREAM_MAIN];
    cudaMemsetAsync(g_nl.error, 0, sizeof(nonlinear_gpu_error_t), stream);
    const int block_size = 128;
    const int grid_size = (g_nl.count + block_size - 1)/block_size;
#ifdef __CUDACC__
    kernel_nonlinear_update_force<<<grid_size, block_size, 0, stream>>>(
        g_nl, g_config, solver->gpuDataDevice, step, delta_t,
        delta_t_squared);
#else
    (void)grid_size;
#endif
    const cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess) {
        fprintf(stderr, "Nonlinear GPU kernel launch failed: %s\n",
                cudaGetErrorString(error));
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
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
    if (error == cudaSuccess) error = NL_DOWNLOAD(strains1, strains1);
    if (error == cudaSuccess) error = NL_DOWNLOAD(pstrains1, pstrains1);
    if (error == cudaSuccess) error = NL_DOWNLOAD(pstrains2, pstrains2);
    if (error == cudaSuccess) error = NL_DOWNLOAD(alphastress1, alphastress1);
    if (error == cudaSuccess) error = NL_DOWNLOAD(alphastress2, alphastress2);
    if (error == cudaSuccess) error = NL_DOWNLOAD(ep1, ep1);
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

void nonlinear_gpu_check_error(int32_t myID)
{
    if (!g_active || g_nl.count == 0) return;
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
                reason = "first Modified-Euler stage is invalid"; break;
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
                    "segment_projections=%d\n"
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
                    g_error_host.md_progress_fraction,
                    g_error_host.md_substep_fraction, g_error_host.value,
                    g_error_host.residual,
                    g_error_host.md_effective_mean_stress,
                    g_error_host.md_strain_increment_norm);
            fprintf(stderr,
                    "  Modified-Euler stages: "
                    "dgamma=(%12.8E, %12.8E) "
                    "raw_dgamma=(%12.8E, %12.8E) "
                    "denominator=(%12.8E, %12.8E)\n"
                    "  stage h=(%12.8E, %12.8E) "
                    "hardening_gap=(%12.8E, %12.8E)\n",
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
        }
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
    }
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
#endif

#ifndef NL_GPU_MATERIAL_POINT_ONLY
void nonlinear_gpu_finalize(void)
{
    nonlinear_gpu_release_station_output();
#define NL_FREE(field) do { if (g_nl.field != NULL) cudaFree(g_nl.field); } while (0)
    NL_FREE(elements); NL_FREE(constants);
    NL_FREE(stresses); NL_FREE(strains); NL_FREE(strains1);
    NL_FREE(pstrains1); NL_FREE(pstrains2);
    NL_FREE(alphastress1); NL_FREE(alphastress2);
    NL_FREE(ep1); NL_FREE(ep2); NL_FREE(psi_n); NL_FREE(load_unload_n);
    NL_FREE(sv_n); NL_FREE(sv_max); NL_FREE(kappa); NL_FREE(sref);
    NL_FREE(fabric);
    NL_FREE(gamma1d); NL_FREE(tao1d); NL_FREE(ggmax); NL_FREE(error);
#undef NL_FREE
    memset(&g_nl, 0, sizeof(g_nl));
    g_mesh_host = NULL;
    g_active = 0;
}
#endif
