/* -*- C++ -*- */
#ifndef HERCULES_NONLINEAR_PJ_LIQ_V8_H
#define HERCULES_NONLINEAR_PJ_LIQ_V8_H

/* Native, allocation-free PJ-Liq V8 constitutive kernel.
 *
 * Source of truth: pj_model_v8_port/model.py and V8_FORMULATION.md.
 * Tensor order is Hercules's [xx, yy, zz, xy, yz, xz], with physical
 * (not engineering) shear components. Stress and strain are tension-positive.
 *
 * Coupled BrickUP ownership is deliberately outside this allocation-free
 * kernel.  The pressure reconstructed from eps_v_confining is effective
 * skeleton mean stress p', including the response to V8 plastic volumetric
 * strain; it is never written to the BrickUP fluid-pressure field.  In a
 * coupled element Hercules independently solves fluid mass balance
 *
 *     d(zeta) = alpha*d(eps_v) + d(p_f)/M,
 *
 * and equilibrium uses sigma_total=sigma_effective-alpha*p_f*I.  Therefore
 * the standalone diagnostic 1-p'/p'_anchor must not be injected as p_f.
 */

#include <float.h>
#include <math.h>
#include <stdint.h>

#ifdef PJ8_STANDALONE
typedef struct tensor_t {
    double xx, yy, zz, xy, yz, xz;
} tensor_t;
#else
#include "nonlinear.h"
#endif

#if defined(__CUDACC__)
#define PJ8_HD __host__ __device__
#else
#define PJ8_HD
#endif

#define PJ8_PARAMETER_COUNT 116
#define PJ8_CONSTITUTIVE_PARAMETER_COUNT 112
#define PJ8_STATE_VALUE_COUNT 93
#define PJ8_PARAMETER_SHA256 \
    "9585f0c155c9444885c5115e7753a1a5d97c783e2c685938a49a65435c9e8f83"

typedef struct pj8_parameters_t {
    double E_ref, nu, h, m, M, kd, zeta, beta0, overshoot_strain;
    double p_ref, p_min, n_G, q_H, denominator_floor_ratio;
    int32_t state_enabled;
    double e_max, e_min, p_atm, Q, R, n_b, n_d, relative_state_limit;
    int32_t fabric_enabled;
    double c_z, z_max;
    int32_t contraction_gate_enabled;
    double C_D, C_in_ratio;
    int32_t adaptive;
    double integration_rtol, integration_atol;
    int32_t max_refinement_depth, localize_events;
    int32_t compatibility_enabled, reversible_enabled, irreversible_enabled;
    double reversible_ratio, reversible_release_rate, irreversible_decay;
    int32_t amplitude_scaling_enabled;
    double cyclic_amplitude_reference, cyclic_amplitude_exponent;
    double low_amplitude_knee_ratio, low_amplitude_exponent;
    double amplitude_factor_minimum, amplitude_factor_maximum;
    int32_t state_contraction_enabled;
    double reference_relative_state, dense_state_contraction_exponent;
    double loose_state_contraction_exponent, state_factor_minimum;
    double state_factor_maximum;
    int32_t cyclic_flow_correction_enabled;
    double cyclic_shear_modulus_reduction, cyclic_hardening_boost;
    double cyclic_flow_pressure_exponent;
    int32_t cyclic_flow_minimum_reversals;
    int32_t state_dependent_knee_enabled;
    double loose_knee_state_exponent, dense_knee_state_exponent;
    double effective_knee_minimum;
    int32_t static_bias_enabled;
    double bias_hardening_intercept, bias_intercept_exponent;
    double bias_crossing_decay, bias_hardening_scale, bias_margin_exponent;
    double bias_amplitude_ratio, bias_pressure_exponent;
    double bias_reference_pressure, bias_confinement_exponent;
    int32_t bias_minimum_reversals;
    double bias_contraction_scale, bias_contraction_exponent;
    double bias_reversible_scale, bias_reversible_exponent;
    int32_t bias_reversible_volume_enabled;
    double bias_reversible_volume_amplitude;
    double bias_reversible_volume_reference_bias;
    double bias_reversible_volume_bias_exponent;
    double bias_reversible_volume_pressure_exponent;
    double bias_reversible_volume_buildup_reversals;
    double bias_reversible_mean_scale;
    double bias_reversible_mean_transition_pressure;
    double bias_reversible_mean_buildup_reversals;
    int32_t bias_ratchet_enabled;
    double bias_ratchet_rate, bias_ratchet_limit;
    double bias_ratchet_amplitude_onset, bias_ratchet_amplitude_full;
    double bias_ratchet_ratio_full, bias_ratchet_ratio_cutoff;
    double bias_ratchet_reference_bias, bias_ratchet_bias_exponent;
    double bias_ratchet_pressure_exponent;
    int32_t state_shakedown_enabled;
    double state_shakedown_scale, state_shakedown_state_sensitivity;
    double state_shakedown_reference_void_ratio;
    double state_shakedown_anchor_width, state_shakedown_bias_reference;
    double state_shakedown_bias_exponent, state_shakedown_pressure_exponent;
    int32_t state_shakedown_minimum_reversals;
    double state_shakedown_factor_maximum;
    double state_shakedown_compliance_exponent;
    double state_shakedown_shear_multiplier_minimum;
    double state_shakedown_dense_hardening_scale;
    double state_shakedown_dense_hardening_amplitude_onset;
    double state_shakedown_dense_hardening_amplitude_decay;
    double state_shakedown_dense_compliance_decay;
    int32_t objective_reversal_enabled;
    double reversal_direction_cosine, reversal_strain_deadband;
    double reversal_stress_deadband_ratio;
} pj8_parameters_t;

typedef struct pj8_state_t {
    tensor_t stress, alpha, alpha0, alpha01, n, fabric;
    double D, beta, lambda_total, ep_eq_since_reversal, void_ratio;
    int64_t reversals, pressure_floor_hits, denominator_floor_hits;
    int64_t beta_fallbacks;
    double eps_v_total, eps_v_confining, eps_v_irreversible;
    double eps_v_reversible, pressure_anchor, D_ir, D_re;
    tensor_t last_reversal_deviator;
    double cyclic_amplitude, amplitude_factor;
    int64_t amplitude_reversals;
    double initial_relative_state, state_contraction_factor;
    double effective_knee_ratio;
    tensor_t geostatic_deviator, static_bias_tensor, cyclic_direction;
    double static_bias_index;
    int32_t cyclic_phase_active;
    double bias_ratchet_strain, physical_eps_v_total;
    double bias_reversible_volume;
    tensor_t last_host_deviatoric_strain_direction;
    int32_t initialized;
} pj8_state_t;

typedef struct pj8_update_info_t {
    int32_t accepted_substeps;
    int32_t reversal_registered;
} pj8_update_info_t;

PJ8_HD static inline double pj8_min(double a, double b) { return a < b ? a : b; }
PJ8_HD static inline double pj8_max(double a, double b) { return a > b ? a : b; }
PJ8_HD static inline double pj8_clip(double x, double a, double b)
{ return pj8_min(pj8_max(x, a), b); }
PJ8_HD static inline tensor_t pj8_zero(void)
{ tensor_t a = {0,0,0,0,0,0}; return a; }
PJ8_HD static inline tensor_t pj8_add(tensor_t a, tensor_t b)
{ tensor_t c={a.xx+b.xx,a.yy+b.yy,a.zz+b.zz,a.xy+b.xy,a.yz+b.yz,a.xz+b.xz}; return c; }
PJ8_HD static inline tensor_t pj8_sub(tensor_t a, tensor_t b)
{ tensor_t c={a.xx-b.xx,a.yy-b.yy,a.zz-b.zz,a.xy-b.xy,a.yz-b.yz,a.xz-b.xz}; return c; }
PJ8_HD static inline tensor_t pj8_scale(tensor_t a, double x)
{ tensor_t c={a.xx*x,a.yy*x,a.zz*x,a.xy*x,a.yz*x,a.xz*x}; return c; }
PJ8_HD static inline tensor_t pj8_iso(double x)
{ tensor_t a={x,x,x,0,0,0}; return a; }
PJ8_HD static inline double pj8_trace(tensor_t a) { return a.xx+a.yy+a.zz; }
PJ8_HD static inline double pj8_ddot(tensor_t a, tensor_t b)
{ return a.xx*b.xx+a.yy*b.yy+a.zz*b.zz+2.0*(a.xy*b.xy+a.yz*b.yz+a.xz*b.xz); }
PJ8_HD static inline double pj8_norm(tensor_t a)
{ return sqrt(pj8_max(pj8_ddot(a,a),0.0)); }
PJ8_HD static inline tensor_t pj8_dev(tensor_t a)
{ return pj8_sub(a,pj8_iso(pj8_trace(a)/3.0)); }
PJ8_HD static inline double pj8_pressure(tensor_t stress)
{ return -pj8_trace(stress)/3.0; }
PJ8_HD static inline double pj8_q(tensor_t stress)
{ tensor_t s=pj8_dev(stress); return sqrt(pj8_max(1.5*pj8_ddot(s,s),0.0)); }
PJ8_HD static inline int pj8_finite_tensor(tensor_t a)
{ return isfinite(a.xx)&&isfinite(a.yy)&&isfinite(a.zz)&&isfinite(a.xy)&&isfinite(a.yz)&&isfinite(a.xz); }

/* stress_scale converts the frozen kPa calibration to the caller's stress
 * unit. Use 1 for the golden oracle and 1000 for Hercules SI production. */
PJ8_HD static inline pj8_parameters_t pj8_reference_parameters(double stress_scale)
{
    pj8_parameters_t p = {};
    p.E_ref=127339.75550887753*stress_scale; p.nu=0.30;
    p.h=122.44207260468994; p.m=0.945; p.M=1.25; p.kd=1.125;
    p.zeta=0.025; p.beta0=1.0e3; p.overshoot_strain=5.0e-5;
    p.p_ref=101.3*stress_scale; p.p_min=1.0e-3*stress_scale;
    p.n_G=0.65; p.q_H=0.35; p.denominator_floor_ratio=0.02;
    p.state_enabled=1; p.e_max=0.78; p.e_min=0.51;
    p.p_atm=101.3*stress_scale; p.Q=10.0; p.R=1.5;
    p.n_b=0.05; p.n_d=0.02; p.relative_state_limit=0.75;
    p.fabric_enabled=1; p.c_z=600.0; p.z_max=4.0;
    p.contraction_gate_enabled=1; p.C_D=0.10; p.C_in_ratio=0.01;
    p.adaptive=1; p.integration_rtol=3.0e-3; p.integration_atol=1.0e-8;
    p.max_refinement_depth=6; p.localize_events=1;
    p.compatibility_enabled=1; p.reversible_enabled=1;
    p.irreversible_enabled=1; p.reversible_ratio=7.10;
    p.reversible_release_rate=300.0; p.irreversible_decay=0.0;
    p.amplitude_scaling_enabled=1; p.cyclic_amplitude_reference=0.409;
    p.cyclic_amplitude_exponent=3.3; p.low_amplitude_knee_ratio=0.85;
    p.low_amplitude_exponent=6.0; p.amplitude_factor_minimum=0.05;
    p.amplitude_factor_maximum=12.0; p.state_contraction_enabled=1;
    p.reference_relative_state=0.4405633323;
    p.dense_state_contraction_exponent=17.0;
    p.loose_state_contraction_exponent=6.5; p.state_factor_minimum=0.005;
    p.state_factor_maximum=5.0; p.cyclic_flow_correction_enabled=1;
    p.cyclic_shear_modulus_reduction=0.85; p.cyclic_hardening_boost=0.0;
    p.cyclic_flow_pressure_exponent=0.5; p.cyclic_flow_minimum_reversals=1;
    p.state_dependent_knee_enabled=1; p.loose_knee_state_exponent=3.7;
    p.dense_knee_state_exponent=0.0; p.effective_knee_minimum=0.20;
    p.static_bias_enabled=1; p.bias_hardening_intercept=50.0;
    p.bias_intercept_exponent=2.3; p.bias_crossing_decay=3.0;
    p.bias_hardening_scale=335.0; p.bias_margin_exponent=1.5;
    p.bias_amplitude_ratio=1.0; p.bias_pressure_exponent=0.5;
    p.bias_reference_pressure=26.266666666666666*stress_scale;
    p.bias_confinement_exponent=0.85; p.bias_minimum_reversals=1;
    p.bias_contraction_scale=0.0; p.bias_contraction_exponent=1.0;
    p.bias_reversible_scale=0.0; p.bias_reversible_exponent=1.0;
    p.bias_reversible_volume_enabled=1;
    p.bias_reversible_volume_amplitude=7.6e-5;
    p.bias_reversible_volume_reference_bias=0.646;
    p.bias_reversible_volume_bias_exponent=2.0;
    p.bias_reversible_volume_pressure_exponent=0.22;
    p.bias_reversible_volume_buildup_reversals=6.0;
    p.bias_reversible_mean_scale=9.0e-5;
    p.bias_reversible_mean_transition_pressure=40.0*stress_scale;
    p.bias_reversible_mean_buildup_reversals=12.0;
    p.bias_ratchet_enabled=1; p.bias_ratchet_rate=2.0;
    p.bias_ratchet_limit=0.0215; p.bias_ratchet_amplitude_onset=0.30;
    p.bias_ratchet_amplitude_full=0.43; p.bias_ratchet_ratio_full=0.70;
    p.bias_ratchet_ratio_cutoff=0.80; p.bias_ratchet_reference_bias=0.646;
    p.bias_ratchet_bias_exponent=0.62; p.bias_ratchet_pressure_exponent=1.0;
    p.state_shakedown_enabled=1; p.state_shakedown_scale=5.0;
    p.state_shakedown_state_sensitivity=2.0;
    p.state_shakedown_reference_void_ratio=0.601;
    p.state_shakedown_anchor_width=0.02;
    p.state_shakedown_bias_reference=0.646;
    p.state_shakedown_bias_exponent=0.5;
    p.state_shakedown_pressure_exponent=0.0;
    p.state_shakedown_minimum_reversals=1;
    p.state_shakedown_factor_maximum=50.0;
    p.state_shakedown_compliance_exponent=1.5;
    p.state_shakedown_shear_multiplier_minimum=0.03;
    p.state_shakedown_dense_hardening_scale=1.0;
    p.state_shakedown_dense_hardening_amplitude_onset=0.82;
    p.state_shakedown_dense_hardening_amplitude_decay=3.5;
    p.state_shakedown_dense_compliance_decay=200.0;
    p.objective_reversal_enabled=1; p.reversal_direction_cosine=-0.20;
    p.reversal_strain_deadband=1.0e-12;
    p.reversal_stress_deadband_ratio=1.0e-4;
    return p;
}

/* PJ-Liq V8 freezes e_max and e_min, so relative density is the clearer
 * external state descriptor. Keep void ratio as an internal constitutive
 * variable and perform this affine conversion at the input boundary. */
PJ8_HD static inline double pj8_void_ratio_from_relative_density(
    const pj8_parameters_t *p, double relative_density)
{
    return p->e_max-relative_density*(p->e_max-p->e_min);
}

PJ8_HD static inline double pj8_relative_state(
    const pj8_parameters_t *p, double pressure, double void_ratio)
{
    if (!p->state_enabled) return 0.0;
    const double dr=pj8_clip((p->e_max-void_ratio)/(p->e_max-p->e_min),0.0,1.0);
    const double arg=pj8_max(100.0*pj8_max(pressure,p->p_min)/p->p_atm,1.0e-12);
    const double critical=pj8_clip(p->R/pj8_max(p->Q-log(arg),0.10),0.0,1.0);
    return pj8_clip(dr-critical,-p->relative_state_limit,p->relative_state_limit);
}

PJ8_HD static inline void pj8_surfaces(const pj8_parameters_t *p,
    double pressure,double void_ratio,double *mb,double *md,double *xi)
{
    *xi=pj8_relative_state(p,pressure,void_ratio);
    *mb=p->M*exp(p->n_b*(*xi)); *md=p->kd*exp(-p->n_d*(*xi));
}

PJ8_HD static inline void pj8_moduli(const pj8_parameters_t *p,double pressure,
    double *shear,double *bulk)
{
    const double ratio=pj8_max(pressure,p->p_min)/p->p_ref;
    const double young=p->E_ref*pow(ratio,p->n_G);
    *shear=young/(2.0*(1.0+p->nu));
    *bulk=young/(3.0*(1.0-2.0*p->nu));
}

PJ8_HD static inline double pj8_state_factor(const pj8_parameters_t *p,double xi)
{
    if (!p->state_contraction_enabled) return 1.0;
    const double d=xi-p->reference_relative_state;
    const double exponent=d>=0.0?p->dense_state_contraction_exponent:p->loose_state_contraction_exponent;
    return pj8_clip(exp(-exponent*d),p->state_factor_minimum,p->state_factor_maximum);
}

PJ8_HD static inline double pj8_effective_knee(const pj8_parameters_t *p,double xi)
{
    if (!p->state_dependent_knee_enabled) return p->low_amplitude_knee_ratio;
    const double d=xi-p->reference_relative_state;
    const double exponent=d>=0.0?p->dense_knee_state_exponent:p->loose_knee_state_exponent;
    return pj8_clip(p->low_amplitude_knee_ratio*exp(exponent*d),p->effective_knee_minimum,1.0);
}

PJ8_HD static inline double pj8_amplitude_factor(const pj8_parameters_t *p,
    double amplitude,double knee)
{
    if (!p->amplitude_scaling_enabled) return 1.0;
    const double ratio=pj8_max(amplitude/p->cyclic_amplitude_reference,1.0e-12);
    const double f=pow(ratio,p->cyclic_amplitude_exponent)*
        pow(pj8_min(1.0,ratio/knee),p->low_amplitude_exponent);
    return pj8_clip(f,p->amplitude_factor_minimum,p->amplitude_factor_maximum);
}

PJ8_HD static inline double pj8_projected_bias(const pj8_state_t *s)
{
    if (!s->cyclic_phase_active) return 0.0;
    const double norm=pj8_norm(s->cyclic_direction);
    return norm<=1.0e-12?0.0:fabs(pj8_ddot(s->static_bias_tensor,s->cyclic_direction))/norm;
}

PJ8_HD static inline double pj8_initial_relative_density(
    const pj8_parameters_t *p,const pj8_state_t *s)
{
    const double arg=pj8_max(100.0*pj8_max(s->pressure_anchor,p->p_min)/p->p_atm,1.0e-12);
    return pj8_clip(s->initial_relative_state+p->R/pj8_max(p->Q-log(arg),0.10),0.0,1.0);
}

PJ8_HD static inline double pj8_reference_relative_density(const pj8_parameters_t *p)
{ return (p->e_max-p->state_shakedown_reference_void_ratio)/(p->e_max-p->e_min); }

PJ8_HD static inline double pj8_raw_shakedown(const pj8_parameters_t *p,
    double pressure,const pj8_state_t *s)
{
    if (!p->state_shakedown_enabled || !s->cyclic_phase_active ||
        s->amplitude_reversals<p->state_shakedown_minimum_reversals) return 1.0;
    const double bias=pj8_projected_bias(s);
    if (bias<=1.0e-14) return 1.0;
    const double distance=pj8_initial_relative_density(p,s)-pj8_reference_relative_density(p);
    const double state_term=exp(p->state_shakedown_state_sensitivity*distance);
    const double anchor_release=1.0-exp(-pow(distance/p->state_shakedown_anchor_width,2.0));
    const double bias_term=pow(bias/p->state_shakedown_bias_reference,
                               p->state_shakedown_bias_exponent);
    const double pressure_term=pow(pj8_clip(pressure/pj8_max(s->pressure_anchor,p->p_min),0.0,1.0),
                                   p->state_shakedown_pressure_exponent);
    return pj8_clip(1.0+p->state_shakedown_scale*anchor_release*state_term*
                    bias_term*pressure_term,1.0,p->state_shakedown_factor_maximum);
}

PJ8_HD static inline double pj8_dense_state_weight(const pj8_parameters_t *p,
    const pj8_state_t *s)
{
    const double d=pj8_max(pj8_initial_relative_density(p,s)-pj8_reference_relative_density(p),0.0);
    return 1.0-exp(-pow(d/p->state_shakedown_anchor_width,2.0));
}

PJ8_HD static inline double pj8_shakedown(const pj8_parameters_t *p,
    double pressure,const pj8_state_t *s)
{
    const double raw=pj8_raw_shakedown(p,pressure,s);
    const double excess=pj8_max(s->cyclic_amplitude-
        p->state_shakedown_dense_hardening_amplitude_onset,0.0);
    const double dense_scale=p->state_shakedown_dense_hardening_scale*
        exp(-p->state_shakedown_dense_hardening_amplitude_decay*excess);
    const double branch=1.0+pj8_dense_state_weight(p,s)*(dense_scale-1.0);
    return pj8_clip(1.0+(raw-1.0)*branch,1.0,p->state_shakedown_factor_maximum);
}

PJ8_HD static inline double pj8_shakedown_compliance(const pj8_parameters_t *p,
    double pressure,const pj8_state_t *s)
{
    const double raw=pj8_raw_shakedown(p,pressure,s);
    const double d=pj8_max(pj8_initial_relative_density(p,s)-pj8_reference_relative_density(p),0.0);
    const double attenuation=exp(-p->state_shakedown_dense_compliance_decay*d*d);
    return 1.0+(raw-1.0)*attenuation;
}

PJ8_HD static inline void pj8_cyclic_flow(const pj8_parameters_t *p,
    double pressure,const pj8_state_t *s,double *shear_factor,double *hardening_factor)
{
    if (!p->cyclic_flow_correction_enabled ||
        s->amplitude_reversals<p->cyclic_flow_minimum_reversals) {
        *shear_factor=1.0; *hardening_factor=1.0; return;
    }
    const double ratio=pj8_clip(pressure/pj8_max(s->pressure_anchor,p->p_min),0.0,1.0);
    const double activity=pow(ratio,p->cyclic_flow_pressure_exponent);
    *shear_factor=1.0-p->cyclic_shear_modulus_reduction*activity;
    *hardening_factor=1.0+p->cyclic_hardening_boost*activity;
}

PJ8_HD static inline double pj8_bias_hardening(const pj8_parameters_t *p,
    const pj8_state_t *s)
{
    if (!p->static_bias_enabled || !s->cyclic_phase_active ||
        s->amplitude_reversals<p->bias_minimum_reversals) return 0.0;
    const double bias=pj8_projected_bias(s);
    if (bias<=1.0e-14) return 0.0;
    const double amplitude=p->bias_amplitude_ratio*s->cyclic_amplitude;
    const double margin=pj8_max(bias-amplitude,0.0);
    const double crossing=pj8_max(amplitude/bias-1.0,0.0);
    return p->bias_hardening_intercept*pow(bias,p->bias_intercept_exponent)*
        exp(-p->bias_crossing_decay*crossing)+
        p->bias_hardening_scale*pow(margin,p->bias_margin_exponent);
}

PJ8_HD static inline void pj8_moduli_for_state(const pj8_parameters_t *p,
    double pressure,const pj8_state_t *s,double *shear,double *bulk)
{
    double sf,hf; pj8_moduli(p,pressure,shear,bulk);
    pj8_cyclic_flow(p,pressure,s,&sf,&hf); (void)hf;
    *shear *= sf;
    const double compliance=pj8_shakedown_compliance(p,pressure,s);
    *shear *= pj8_clip(pow(compliance,-p->state_shakedown_compliance_exponent),
                       p->state_shakedown_shear_multiplier_minimum,1.0);
}

PJ8_HD static inline double pj8_hardening_for_state(const pj8_parameters_t *p,
    double pressure,const pj8_state_t *s)
{
    const double ratio=pj8_max(pressure,p->p_min)/p->p_ref;
    double value=p->h*pow(ratio,-p->q_H);
    double sf,hf; pj8_cyclic_flow(p,pressure,s,&sf,&hf); (void)sf;
    value *= hf;
    const double boost=pj8_bias_hardening(p,s);
    const double pressure_ratio=pj8_clip(pressure/pj8_max(s->pressure_anchor,p->p_min),0.0,1.0);
    const double confinement=pow(pj8_clip(p->bias_reference_pressure/
        pj8_max(s->pressure_anchor,p->p_min),0.25,4.0),p->bias_confinement_exponent);
    value *= 1.0+boost*confinement*pow(pressure_ratio,p->bias_pressure_exponent);
    return value*pj8_shakedown(p,pressure,s);
}

PJ8_HD static inline double pj8_bulk_reference(const pj8_parameters_t *p)
{ double g,k; pj8_moduli(p,p->p_ref,&g,&k); return k; }

PJ8_HD static inline double pj8_confining_strain(const pj8_parameters_t *p,
    double pressure,double anchor)
{
    pressure=pj8_max(pressure,p->p_min); anchor=pj8_max(anchor,p->p_min);
    const double bulk=pj8_bulk_reference(p);
    if (fabs(p->n_G-1.0)<1.0e-12) return p->p_ref/bulk*log(anchor/pressure);
    const double factor=pow(p->p_ref,p->n_G)/(bulk*(1.0-p->n_G));
    return factor*(pow(anchor,1.0-p->n_G)-pow(pressure,1.0-p->n_G));
}

PJ8_HD static inline double pj8_pressure_from_confining(
    const pj8_parameters_t *p,double eps_vc,double anchor,int *at_floor)
{
    const double threshold=pj8_confining_strain(p,p->p_min,anchor);
    if (eps_vc>=threshold) { *at_floor=1; return p->p_min; }
    *at_floor=0;
    const double bulk=pj8_bulk_reference(p);
    double pressure;
    if (fabs(p->n_G-1.0)<1.0e-12)
        pressure=anchor*exp(-bulk*eps_vc/p->p_ref);
    else {
        double base=pow(anchor,1.0-p->n_G)-
            (1.0-p->n_G)*bulk*eps_vc/pow(p->p_ref,p->n_G);
        base=pj8_max(base,pow(p->p_min,1.0-p->n_G));
        pressure=pow(base,1.0/(1.0-p->n_G));
    }
    return pj8_max(pressure,p->p_min);
}

PJ8_HD static inline int pj8_initialize(const pj8_parameters_t *p,
    tensor_t stress,double void_ratio,pj8_state_t *state)
{
    const double pressure=pj8_pressure(stress);
    if (!(pressure>0.0)||!isfinite(pressure)||!isfinite(void_ratio)) return 0;
    *state=pj8_state_t{};
    const tensor_t dev=pj8_dev(stress);
    double mb,md,xi; pj8_surfaces(p,pressure,void_ratio,&mb,&md,&xi); (void)mb;
    state->stress=stress; state->alpha=pj8_scale(dev,1.0/pressure);
    state->D_ir=p->irreversible_enabled?p->zeta*sqrt(2.0/3.0)*md:0.0;
    state->D=state->D_ir; state->beta=p->beta0;
    state->void_ratio=void_ratio; state->pressure_anchor=pressure;
    state->last_reversal_deviator=dev; state->amplitude_factor=1.0;
    state->initial_relative_state=xi; state->state_contraction_factor=pj8_state_factor(p,xi);
    state->effective_knee_ratio=pj8_effective_knee(p,xi);
    state->geostatic_deviator=dev; state->initialized=1;
    return 1;
}

PJ8_HD static inline int pj8_begin_dynamic_phase(const pj8_parameters_t *p,
    tensor_t reference_stress,pj8_state_t *state)
{
    if (!state->initialized) return 0;
    const tensor_t current=pj8_dev(state->stress);
    const tensor_t reference=pj8_dev(reference_stress);
    const tensor_t bias=pj8_sub(current,reference);
    const double anchor=pj8_max(state->pressure_anchor,p->p_min);
    state->static_bias_tensor=pj8_scale(bias,1.0/anchor);
    state->static_bias_index=pj8_norm(bias)/anchor;
    state->cyclic_direction=pj8_zero(); state->cyclic_phase_active=1;
    state->bias_ratchet_strain=0.0; state->bias_reversible_volume=0.0;
    state->eps_v_total=state->physical_eps_v_total;
    state->last_host_deviatoric_strain_direction=pj8_zero();
    return 1;
}

PJ8_HD static inline void pj8_dilatancy(const pj8_parameters_t *p,
    tensor_t alpha,tensor_t normal,double beta,tensor_t fabric,
    double pressure,double void_ratio,double eps_v_ir,double eps_v_re,
    double *D_ir,double *D_re)
{
    double mb,md,xi; pj8_surfaces(p,pressure,void_ratio,&mb,&md,&xi); (void)xi;
    const double eta=pj8_ddot(alpha,normal);
    const double eta_d=sqrt(2.0/3.0)*md;
    *D_ir=0.0;
    if (p->irreversible_enabled) {
        const double distance=pj8_max(eta_d-eta,0.0);
        const double amplitude=p->zeta*(1.0+pj8_max(pj8_ddot(fabric,normal),0.0));
        double gate=1.0;
        if (p->contraction_gate_enabled && p->C_D>0.0) {
            const double reversal=(sqrt(2.0/3.0)*mb-eta)/pj8_max(beta,1.0e-12);
            const double c_in=p->C_in_ratio*sqrt(2.0/3.0)*mb;
            gate=(reversal+c_in)*(reversal+c_in)/(reversal*reversal+p->C_D);
        }
        *D_ir=amplitude*gate*distance*
            exp(-p->irreversible_decay*pj8_max(-eps_v_ir,0.0));
    }
    *D_re=0.0;
    if (p->reversible_enabled) {
        if (eta>eta_d) *D_re=-p->zeta*p->reversible_ratio*(eta-eta_d);
        else if (eps_v_re>0.0)
            *D_re=pj8_min(p->reversible_release_rate*eps_v_re,pj8_max(eta_d,1.0e-12));
    }
}

PJ8_HD static inline double pj8_irreversible_factor(const pj8_parameters_t *p,
    const pj8_state_t *s)
{
    double factor=s->amplitude_factor*s->state_contraction_factor;
    const double bias=pj8_projected_bias(s);
    factor *= 1.0+p->bias_contraction_scale*pow(bias,p->bias_contraction_exponent);
    return factor;
}

PJ8_HD static inline void pj8_register_reversal(const pj8_parameters_t *p,
    const pj8_state_t *old,tensor_t *alpha0,tensor_t *alpha01,
    double *ep_eq,double *beta,int64_t *reversals)
{
    const double ratio=(*ep_eq)/p->overshoot_strain;
    const double modifier=ratio>1.0?0.0:1.0-ratio;
    const tensor_t previous=*alpha0;
    *alpha0=pj8_norm(previous)==0.0?old->alpha:
        pj8_add(pj8_scale(*alpha01,modifier),pj8_scale(old->alpha,1.0-modifier));
    *alpha01=previous; *ep_eq=0.0; *beta=p->beta0;
    *reversals=old->reversals+1;
}

PJ8_HD static inline pj8_state_t pj8_backbone_forward_euler(
    const pj8_parameters_t *p,const pj8_state_t *old,tensor_t deps,
    int force_reversal,int allow_legacy_reversal)
{
    pj8_state_t state=*old;
    const tensor_t s_t=pj8_dev(old->stress);
    const double p_t=pj8_max(pj8_pressure(old->stress),p->p_min);
    double shear,bulk; pj8_moduli_for_state(p,p_t,old,&shear,&bulk);
    const double h_eff=pj8_hardening_for_state(p,p_t,old);
    const double threshold=pj8_confining_strain(p,p->p_min,old->pressure_anchor);
    const int floor_active=p->compatibility_enabled && old->eps_v_confining>=threshold;
    const double coupling_bulk=floor_active?0.0:bulk;

    const double deps_v=pj8_trace(deps);
    const tensor_t deps_d=pj8_sub(deps,pj8_iso(deps_v/3.0));
    const tensor_t s_trial=pj8_add(s_t,pj8_scale(deps_d,2.0*shear));
    double p_trial;
    if (p->compatibility_enabled) {
        int at_floor=0;
        p_trial=pj8_pressure_from_confining(p,old->eps_v_confining+deps_v,
                                            old->pressure_anchor,&at_floor);
    } else p_trial=pj8_max(p_t-bulk*deps_v,p->p_min);
    const tensor_t alpha_trial=pj8_scale(s_trial,1.0/pj8_max(p_trial,p->p_min));
    const tensor_t delta_alpha_trial=pj8_sub(alpha_trial,old->alpha);

    tensor_t alpha0=old->alpha0,alpha01=old->alpha01;
    double beta_t=old->beta,ep_eq=old->ep_eq_since_reversal;
    int64_t reversals=old->reversals;
    const double legacy=pj8_ddot(pj8_sub(alpha_trial,alpha0),delta_alpha_trial);
    const int reversal=force_reversal || (allow_legacy_reversal && legacy<0.0);
    if (reversal) pj8_register_reversal(p,old,&alpha0,&alpha01,&ep_eq,&beta_t,&reversals);

    const double alpha_n_t=pj8_ddot(old->alpha,old->n);
    const double hardening=p_t*h_eff*pow(pj8_max(beta_t,1.0e-12),p->m);
    double denominator=2.0*shear+(2.0/3.0)*hardening-
        coupling_bulk*old->D*alpha_n_t;
    const double floor=p->denominator_floor_ratio*2.0*shear;
    int64_t denominator_hits=old->denominator_floor_hits;
    if (denominator<floor) { denominator=floor; denominator_hits++; }
    const double numerator=2.0*shear*pj8_ddot(deps_d,old->n)+
        coupling_bulk*deps_v*alpha_n_t;
    const double d_lambda=pj8_max(0.0,numerator/denominator);
    tensor_t s=pj8_add(s_t,pj8_scale(pj8_sub(deps_d,pj8_scale(old->n,d_lambda)),
                                     2.0*shear));

    const double d_ir=pj8_max(old->D_ir,0.0);
    double d_re=old->D_re;
    if (d_lambda>0.0 && d_re>0.0)
        d_re=pj8_min(d_re,old->eps_v_reversible/d_lambda);
    const double eps_v_ir=old->eps_v_irreversible-d_lambda*d_ir;
    const double eps_v_re=pj8_max(old->eps_v_reversible-d_lambda*d_re,0.0);
    const double eps_v_total=old->eps_v_total+deps_v;
    const double eps_v_confining=eps_v_total-eps_v_ir-eps_v_re;
    double pressure; int at_floor=0;
    if (p->compatibility_enabled)
        pressure=pj8_pressure_from_confining(p,eps_v_confining,
                                             old->pressure_anchor,&at_floor);
    else {
        pressure=pj8_max(p_t-bulk*(deps_v+d_lambda*(d_ir+d_re)),p->p_min);
        at_floor=pressure<=p->p_min;
    }
    tensor_t stress=pj8_sub(s,pj8_iso(pressure));
    ep_eq += d_lambda;
    const double void_ratio=old->void_ratio+(1.0+old->void_ratio)*deps_v;
    tensor_t fabric=old->fabric;
    if (p->fabric_enabled && d_lambda>0.0) {
        const tensor_t term=pj8_add(pj8_scale(old->n,sqrt(2.0/3.0)*p->z_max),fabric);
        fabric=pj8_add(fabric,pj8_scale(term,-p->c_z*d_lambda*pj8_max(-d_re,0.0)));
    }

    double mb,md,xi; pj8_surfaces(p,pressure,void_ratio,&mb,&md,&xi); (void)md; (void)xi;
    tensor_t alpha=pj8_scale(s,1.0/pressure);
    const double alpha_limit=sqrt(2.0/3.0)*mb;
    const double alpha_norm=pj8_norm(alpha);
    if (alpha_norm>alpha_limit) {
        alpha=pj8_scale(alpha,alpha_limit/alpha_norm);
        s=pj8_scale(alpha,pressure); stress=pj8_sub(s,pj8_iso(pressure));
    }
    const tensor_t offset=pj8_sub(alpha,alpha0);
    const double qa=pj8_ddot(offset,offset);
    const double qb=2.0*pj8_ddot(offset,alpha);
    const double qc=pj8_ddot(alpha,alpha)-(2.0/3.0)*mb*mb;
    const double disc=qb*qb-4.0*qa*qc;
    int64_t beta_fallbacks=old->beta_fallbacks;
    double beta;
    if (qa<=1.0e-14 || disc<0.0) { beta=p->beta0; beta_fallbacks++; }
    else {
        beta=(-qb+sqrt(pj8_max(disc,0.0)))/(2.0*qa);
        if (!isfinite(beta)) { beta=p->beta0; beta_fallbacks++; }
    }
    beta=pj8_max(beta,1.0e-6);
    const tensor_t mapped=pj8_add(alpha,pj8_scale(offset,beta));
    const double mapped_norm=pj8_norm(mapped);
    const tensor_t normal=mapped_norm>1.0e-12?pj8_scale(mapped,1.0/mapped_norm):pj8_zero();
    double D_ir,D_re;
    pj8_dilatancy(p,alpha,normal,beta,fabric,pressure,void_ratio,
                  eps_v_ir,eps_v_re,&D_ir,&D_re);
    D_ir *= pj8_irreversible_factor(p,old);

    state.stress=stress; state.alpha=alpha; state.alpha0=alpha0;
    state.alpha01=alpha01; state.n=normal; state.fabric=fabric;
    state.D_ir=D_ir; state.D_re=D_re; state.D=D_ir+D_re; state.beta=beta;
    state.lambda_total=old->lambda_total+d_lambda;
    state.ep_eq_since_reversal=ep_eq; state.void_ratio=void_ratio;
    state.reversals=reversals;
    state.pressure_floor_hits=old->pressure_floor_hits+(at_floor?1:0);
    state.denominator_floor_hits=denominator_hits;
    state.beta_fallbacks=beta_fallbacks;
    state.eps_v_total=eps_v_total; state.eps_v_confining=eps_v_confining;
    state.eps_v_irreversible=eps_v_ir; state.eps_v_reversible=eps_v_re;
    if (state.reversals>old->reversals) {
        const tensor_t reversal_dev=pj8_dev(old->stress);
        const double excursion=pj8_norm(pj8_sub(reversal_dev,old->last_reversal_deviator));
        const double divisor=old->amplitude_reversals==0?1.0:2.0;
        const double amplitude=excursion/(divisor*pj8_max(old->pressure_anchor,p->p_min));
        if (amplitude>1.0e-12) {
            state.cyclic_amplitude=amplitude;
            state.amplitude_factor=pj8_amplitude_factor(p,amplitude,old->effective_knee_ratio);
        }
        state.last_reversal_deviator=reversal_dev;
        state.amplitude_reversals=old->amplitude_reversals+1;
    }
    return state;
}

PJ8_HD static inline double pj8_bias_reversible_factor(
    const pj8_parameters_t *p,const pj8_state_t *s)
{
    if (!p->static_bias_enabled || !s->cyclic_phase_active) return 1.0;
    return 1.0+p->bias_reversible_scale*
        pow(pj8_projected_bias(s),p->bias_reversible_exponent);
}

PJ8_HD static inline double pj8_bias_reversible_volume_target(
    const pj8_parameters_t *p,const pj8_state_t *s)
{
    if (!p->static_bias_enabled || !p->bias_reversible_volume_enabled ||
        !s->cyclic_phase_active || s->amplitude_reversals<1) return 0.0;
    const double direction_norm=pj8_norm(s->cyclic_direction);
    if (direction_norm<=1.0e-14 || s->cyclic_amplitude<=1.0e-14) return 0.0;
    const double bias=pj8_projected_bias(s);
    if (bias<=1.0e-14) return 0.0;
    const tensor_t phase_direction=pj8_scale(s->cyclic_direction,1.0/direction_norm);
    const tensor_t current=pj8_dev(s->stress);
    const tensor_t static_dev=pj8_add(s->geostatic_deviator,
        pj8_scale(s->static_bias_tensor,s->pressure_anchor));
    const double dynamic_projection=pj8_ddot(pj8_sub(current,static_dev),phase_direction);
    const double phase=pj8_clip(dynamic_projection/
        (s->pressure_anchor*s->cyclic_amplitude),-1.0,1.0);
    const double bias_factor=pow(p->bias_reversible_volume_reference_bias/bias,
                                 p->bias_reversible_volume_bias_exponent);
    const double pressure_factor=pow(pj8_max(s->pressure_anchor,p->p_min)/
        p->bias_reference_pressure,p->bias_reversible_volume_pressure_exponent);
    const double oscillation_buildup=1.0-exp(-(double)s->amplitude_reversals/
        p->bias_reversible_volume_buildup_reversals);
    const double oscillation=-p->bias_reversible_volume_amplitude*bias_factor*
        pressure_factor*oscillation_buildup*phase;
    const double mean_buildup=1.0-exp(-(double)s->amplitude_reversals/
        p->bias_reversible_mean_buildup_reversals);
    const double mean_shift=p->bias_reversible_mean_scale*
        (p->bias_reversible_mean_transition_pressure-s->pressure_anchor)/
        p->bias_reference_pressure*(bias/p->bias_reversible_volume_reference_bias)*
        mean_buildup;
    return oscillation+mean_shift;
}

PJ8_HD static inline double pj8_smoothstep(double value)
{ value=pj8_clip(value,0.0,1.0); return value*value*(3.0-2.0*value); }

PJ8_HD static inline double pj8_bias_ratchet_activity(
    const pj8_parameters_t *p,const pj8_state_t *s)
{
    if (!p->static_bias_enabled || !p->bias_ratchet_enabled ||
        !s->cyclic_phase_active || s->amplitude_reversals<p->bias_minimum_reversals)
        return 0.0;
    const double bias=pj8_projected_bias(s);
    if (bias<=1.0e-14) return 0.0;
    const double amplitude_position=(s->cyclic_amplitude-p->bias_ratchet_amplitude_onset)/
        (p->bias_ratchet_amplitude_full-p->bias_ratchet_amplitude_onset);
    const double ratio_position=(s->cyclic_amplitude/bias-p->bias_ratchet_ratio_full)/
        (p->bias_ratchet_ratio_cutoff-p->bias_ratchet_ratio_full);
    return pj8_smoothstep(amplitude_position)*(1.0-pj8_smoothstep(ratio_position));
}

PJ8_HD static inline double pj8_bias_ratchet_capacity(
    const pj8_parameters_t *p,const pj8_state_t *s)
{
    const double bias=pj8_projected_bias(s);
    if (bias<=1.0e-14) return 0.0;
    const double factor=pow(p->bias_ratchet_reference_bias/bias,
                            p->bias_ratchet_bias_exponent);
    /* The repeated factor is intentional and preserves the frozen V8 oracle. */
    return p->bias_ratchet_limit*pj8_bias_ratchet_activity(p,s)*factor*factor;
}

PJ8_HD static inline double pj8_bias_ratchet_rate_factor(
    const pj8_parameters_t *p,const pj8_state_t *s)
{
    const double capacity=pj8_bias_ratchet_capacity(p,s);
    if (capacity<=1.0e-14 || s->bias_ratchet_strain>=capacity) return 0.0;
    const double stress_level=pow(pj8_max(s->pressure_anchor/
        p->bias_reference_pressure-1.0,0.0),p->bias_ratchet_pressure_exponent);
    const double bias=pj8_projected_bias(s);
    const double bias_factor=pow(p->bias_ratchet_reference_bias/
        pj8_max(bias,1.0e-14),p->bias_ratchet_bias_exponent);
    const double saturation=pj8_max(1.0-s->bias_ratchet_strain/capacity,0.0);
    return p->bias_ratchet_rate*stress_level*bias_factor*saturation;
}

PJ8_HD static inline double pj8_bias_ratchet_increment(
    const pj8_parameters_t *p,const pj8_state_t *old,const pj8_state_t *trial,
    tensor_t *direction)
{
    *direction=pj8_zero();
    const double rate=pj8_bias_ratchet_rate_factor(p,old);
    const double dl=pj8_max(trial->lambda_total-old->lambda_total,0.0);
    const double norm=pj8_norm(old->cyclic_direction);
    const double projection=pj8_ddot(old->static_bias_tensor,old->cyclic_direction);
    if (rate<=0.0 || dl<=0.0 || norm<=1.0e-12) return 0.0;
    const double sign=projection>0.0?1.0:(projection<0.0?-1.0:0.0);
    *direction=pj8_scale(old->cyclic_direction,sign/norm);
    const double capacity=pj8_bias_ratchet_capacity(p,old);
    return pj8_min(rate*dl,pj8_max(capacity-old->bias_ratchet_strain,0.0));
}

PJ8_HD static inline pj8_state_t pj8_without_reversible_bias(
    const pj8_parameters_t *p,const pj8_state_t *state)
{
    if (fabs(state->bias_reversible_volume)<=1.0e-16) return *state;
    pj8_state_t mechanical=*state;
    mechanical.eps_v_total=state->physical_eps_v_total;
    mechanical.eps_v_confining=mechanical.eps_v_total-
        mechanical.eps_v_irreversible-mechanical.eps_v_reversible;
    const tensor_t dev=pj8_dev(state->stress);
    int at_floor=0;
    const double pressure=pj8_pressure_from_confining(p,mechanical.eps_v_confining,
        mechanical.pressure_anchor,&at_floor); (void)at_floor;
    mechanical.alpha=pj8_scale(dev,1.0/pressure);
    mechanical.stress=pj8_sub(dev,pj8_iso(pressure));
    return mechanical;
}

PJ8_HD static inline pj8_state_t pj8_forward_euler(const pj8_parameters_t *p,
    const pj8_state_t *old,tensor_t deps,int force_reversal,int allow_legacy)
{
    const double physical_deps_v=pj8_trace(deps);
    const pj8_state_t mechanical=pj8_without_reversible_bias(p,old);
    const pj8_state_t trial=pj8_backbone_forward_euler(
        p,&mechanical,deps,force_reversal,allow_legacy);
    tensor_t ratchet_direction;
    const double ratchet=pj8_bias_ratchet_increment(p,old,&trial,&ratchet_direction);
    pj8_state_t provisional=trial;
    if (ratchet>0.0) {
        const tensor_t effective_deps=pj8_sub(deps,pj8_scale(ratchet_direction,ratchet));
        provisional=pj8_backbone_forward_euler(
            p,&mechanical,effective_deps,force_reversal,allow_legacy);
    }
    if (provisional.amplitude_reversals>old->amplitude_reversals) {
        const tensor_t reversal_dev=pj8_dev(old->stress);
        const tensor_t excursion=pj8_sub(reversal_dev,old->last_reversal_deviator);
        const double norm=pj8_norm(excursion);
        if (norm>1.0e-12) provisional.cyclic_direction=pj8_scale(excursion,1.0/norm);
    }
    const double target=pj8_bias_reversible_volume_target(p,&provisional);
    pj8_state_t state=provisional;
    state.bias_ratchet_strain=old->bias_ratchet_strain+ratchet;
    state.physical_eps_v_total=old->physical_eps_v_total+physical_deps_v;
    state.bias_reversible_volume=target;
    const double effective_total=state.physical_eps_v_total+target;
    const int changed=fabs(effective_total-state.eps_v_total)>1.0e-16;
    state.eps_v_total=effective_total;
    state.eps_v_confining=state.eps_v_total-state.eps_v_irreversible-
        state.eps_v_reversible;
    if (changed) {
        const tensor_t dev=pj8_dev(state.stress); int at_floor=0;
        const double pressure=pj8_pressure_from_confining(p,state.eps_v_confining,
            state.pressure_anchor,&at_floor);
        state.alpha=pj8_scale(dev,1.0/pressure);
        state.stress=pj8_sub(dev,pj8_iso(pressure));
        state.pressure_floor_hits += at_floor?1:0;
    }
    state.void_ratio=old->void_ratio+(1.0+old->void_ratio)*physical_deps_v;
    state.D_re *= pj8_bias_reversible_factor(p,&state);
    state.D=state.D_ir+state.D_re;
    return state;
}

PJ8_HD static inline tensor_t pj8_host_direction(const pj8_parameters_t *p,
    tensor_t deps,int *valid)
{
    const tensor_t dev=pj8_dev(deps); const double norm=pj8_norm(dev);
    *valid=norm>p->reversal_strain_deadband;
    return *valid?pj8_scale(dev,1.0/norm):pj8_zero();
}

PJ8_HD static inline int pj8_host_reversal(const pj8_parameters_t *p,
    const pj8_state_t *s,tensor_t direction,int direction_valid)
{
    if (!direction_valid || !s->cyclic_phase_active) return 0;
    const double previous_norm=pj8_norm(s->last_host_deviatoric_strain_direction);
    if (previous_norm<=1.0e-14) return 0;
    const double cosine=pj8_ddot(direction,s->last_host_deviatoric_strain_direction)/previous_norm;
    if (cosine>p->reversal_direction_cosine) return 0;
    const double excursion=pj8_norm(pj8_sub(pj8_dev(s->stress),s->last_reversal_deviator));
    return excursion>=p->reversal_stress_deadband_ratio*
        pj8_max(s->pressure_anchor,p->p_min);
}

PJ8_HD static inline int pj8_update(const pj8_parameters_t *p,tensor_t deps,
    int32_t fixed_substeps,pj8_state_t *state,tensor_t *stress_new,
    pj8_update_info_t *info)
{
    if (!state || !state->initialized || fixed_substeps<1 ||
        !pj8_finite_tensor(deps)) return 0;
    /* A rate-independent zero increment must not advance compatibility,
     * reversal, or diagnostic history merely because the host took a step. */
    if (pj8_ddot(deps,deps) == 0.0) {
        if (stress_new) *stress_new=state->stress;
        if (info) { info->accepted_substeps=0;
                    info->reversal_registered=0; }
        return 1;
    }
    int valid=0;
    const tensor_t direction=p->objective_reversal_enabled?
        pj8_host_direction(p,deps,&valid):pj8_zero();
    const int reversal=p->objective_reversal_enabled?
        pj8_host_reversal(p,state,direction,valid):0;
    pj8_state_t current=*state;
    const tensor_t sub=pj8_scale(deps,1.0/(double)fixed_substeps);
    for (int32_t i=0;i<fixed_substeps;i++)
        current=pj8_forward_euler(p,&current,sub,
            p->objective_reversal_enabled && reversal && i==0,
            !p->objective_reversal_enabled);
    if (p->objective_reversal_enabled && valid)
        current.last_host_deviatoric_strain_direction=direction;
    if (!pj8_finite_tensor(current.stress) ||
        !isfinite(current.lambda_total) || !isfinite(current.void_ratio)) return 0;
    *state=current; if (stress_new) *stress_new=current.stress;
    if (info) { info->accepted_substeps=fixed_substeps;
                info->reversal_registered=reversal; }
    return 1;
}

PJ8_HD static inline double pj8_compatibility_residual(const pj8_state_t *s)
{ return s->eps_v_total-s->eps_v_confining-s->eps_v_irreversible-s->eps_v_reversible; }

PJ8_HD static inline int pj8_state_values(
    const pj8_state_t *state,double values[PJ8_STATE_VALUE_COUNT])
{
    int i=0;
#define PJ8_STATE_TENSOR(t) do { \
    values[i++]=(t).xx; values[i++]=(t).yy; values[i++]=(t).zz; \
    values[i++]=(t).xy; values[i++]=(t).yz; values[i++]=(t).xz; \
} while (0)
    PJ8_STATE_TENSOR(state->stress); PJ8_STATE_TENSOR(state->alpha);
    PJ8_STATE_TENSOR(state->alpha0); PJ8_STATE_TENSOR(state->alpha01);
    PJ8_STATE_TENSOR(state->n); PJ8_STATE_TENSOR(state->fabric);
    values[i++]=state->D; values[i++]=state->beta;
    values[i++]=state->lambda_total; values[i++]=state->ep_eq_since_reversal;
    values[i++]=state->void_ratio; values[i++]=(double)state->reversals;
    values[i++]=(double)state->pressure_floor_hits;
    values[i++]=(double)state->denominator_floor_hits;
    values[i++]=(double)state->beta_fallbacks;
    values[i++]=state->eps_v_total; values[i++]=state->eps_v_confining;
    values[i++]=state->eps_v_irreversible; values[i++]=state->eps_v_reversible;
    values[i++]=state->pressure_anchor; values[i++]=state->D_ir;
    values[i++]=state->D_re; PJ8_STATE_TENSOR(state->last_reversal_deviator);
    values[i++]=state->cyclic_amplitude; values[i++]=state->amplitude_factor;
    values[i++]=(double)state->amplitude_reversals;
    values[i++]=state->initial_relative_state;
    values[i++]=state->state_contraction_factor;
    values[i++]=state->effective_knee_ratio;
    PJ8_STATE_TENSOR(state->geostatic_deviator);
    PJ8_STATE_TENSOR(state->static_bias_tensor);
    PJ8_STATE_TENSOR(state->cyclic_direction);
    values[i++]=state->static_bias_index;
    values[i++]=(double)state->cyclic_phase_active;
    values[i++]=state->bias_ratchet_strain;
    values[i++]=state->physical_eps_v_total;
    values[i++]=state->bias_reversible_volume;
    PJ8_STATE_TENSOR(state->last_host_deviatoric_strain_direction);
#undef PJ8_STATE_TENSOR
    return i;
}

#undef PJ8_HD
#endif /* HERCULES_NONLINEAR_PJ_LIQ_V8_H */
