/* -*- C++ -*- */
#ifndef HERCULES_NONLINEAR_RIVA_SAND_INTERMEDIATE_BIAS_RESEARCH_H
#define HERCULES_NONLINEAR_RIVA_SAND_INTERMEDIATE_BIAS_RESEARCH_H

/* Native, allocation-free translation of the private intermediate-bias
 * RIVA-Sand research successor at OpenSees_RIVA commit 637d207e867b863d,
 * with the completed-half-cycle plastic gate from commit 5e4ebd30cb5f43b6.
 *
 * This is deliberately a separate material.  It does not change the frozen
 * production RIVASand/RIVASandCustom kernel in nonlinear_riva_sand.h.  The
 * Python handoff is an offline oracle only; neither this header nor Hercules
 * calls Python at runtime.
 *
 * Tensors use physical components [xx yy zz xy yz xz], stress is tension
 * positive, and p'=-trace(sigma')/3.  The returned stress is effective
 * skeleton stress.  BrickUP continues to own its fluid-pressure DOF and total
 * stress assembly.
 */

#include "nonlinear_riva_sand.h"

#if defined(__CUDACC__)
#define RIVA_IB_HD __host__ __device__
#else
#define RIVA_IB_HD
#endif

#define RIVA_IB_PARAMETER_COUNT 249
#define RIVA_IB_LOGICAL_STATE_COUNT 63
#define RIVA_IB_STATE_VALUE_COUNT 138
#define RIVA_IB_KERNEL_REVISION 3u
#define RIVA_IB_PARAMETER_SHA256 \
    "3c17e962e64d3a3b29d4797399380dfc85669a5b9fd20b187af359f57aa5b56e"
#define RIVA_IB_SOURCE_COMMIT \
    "5e4ebd30cb5f43b677c61306c7f1fbbdcb24905c"

/* The first 116 calibrated values are held in base.  p_residual is a
 * Hercules-only production research control and remains zero here, so it is
 * not counted among the handoff's 249 values. */
typedef struct riva_ib_parameters_t {
    riva_parameters_t base;

    int32_t phase_transformation_enabled;
    double phase_cyclic_shear_modulus_reduction;
    int32_t branch_compliance_enabled;
    double branch_compliance_bias_reference;
    double branch_compliance_bias_exponent;
    double branch_compliance_minimum;
    double phase_compliance_peak;
    double phase_compliance_bell_gain;
    double phase_memory_hardening_reduction;
    double phase_memory_reference_volume;
    double phase_bias_hardening_intercept;
    double phase_bias_hardening_exponent;
    double phase_compliance_shape;
    double phase_compliance_location;
    double phase_compliance_half_width;
    double branch_directional_balance;
    double branch_balance_bias_exponent;
    double branch_balance_bias_cap;
    double phase_ratio;
    double phase_width;
    double phase_contraction_rate;
    double phase_reversible_scale;
    double phase_reversible_mean_scale;
    double phase_reversible_relaxation_strain;
    double phase_reversible_relaxation_bias_exponent;
    double phase_potential_anchor_fraction;
    double phase_pressure_exponent;
    double phase_wave_pressure_exponent;
    double phase_mean_pressure_exponent;
    double phase_amplitude_onset_ratio;
    double phase_amplitude_full_ratio;
    double phase_volume_density_onset;
    double phase_volume_density_full;
    double phase_volume_replacement_density_full;
    double phase_intermediate_wave_multiplier;
    double phase_intermediate_mean_multiplier;
    double phase_intermediate_relaxation_ratio;
    double phase_bias_exponent;
    double phase_wave_bias_exponent;
    double phase_mean_bias_exponent;

    int32_t phase_accumulation_control_enabled;
    double phase_accumulation_density_onset;
    double phase_accumulation_density_peak;
    double phase_accumulation_density_cutoff;
    double phase_accumulation_reference_amplitude;
    double phase_accumulation_memory_onset;
    double phase_accumulation_memory_width;
    double phase_accumulation_amplitude_exponent;
    double phase_accumulation_gain_full_amplitude;
    double phase_accumulation_high_amplitude_gain_ratio;
    double phase_accumulation_hardening_gain;

    int32_t loose_stabilization_enabled;
    double loose_stabilization_density_full;
    double loose_stabilization_density_cutoff;
    double loose_stabilization_bias_onset;
    double loose_stabilization_hardening_multiplier;
    double loose_stabilization_contraction_scale;
    int32_t loose_phase_enabled;
    double loose_phase_activity_scale;
    double loose_phase_wave_activity_scale;
    double loose_phase_mean_activity_scale;
    double loose_phase_replacement_fraction;

    int32_t unbiased_phase_enabled;
    double unbiased_phase_density_onset;
    double unbiased_phase_density_full;
    double unbiased_phase_bias_cutoff;
    double unbiased_phase_activity_scale;
    double unbiased_phase_potential_center;

    int32_t loose_shear_flow_enabled;
    double loose_shear_cyclic_modulus_reduction;
    double loose_shear_modulus_scale;
    double loose_shear_early_hardening_multiplier;
    double loose_shear_late_hardening_multiplier;
    double loose_shear_hardening_memory_onset;
    double loose_shear_hardening_memory_width;
    double loose_shear_branch_compliance_gain;
    double loose_shear_branch_compliance_onset;
    double loose_shear_branch_compliance_full;
    double loose_shear_ratchet_rate;
    double loose_shear_ratchet_capacity;
    double loose_shear_ratchet_pressure_exponent;

    int32_t mapping_backstress_enabled;
    double mapping_bias_onset;
    double mapping_bias_full;
    double mapping_density_full;
    double mapping_density_cutoff;
    double mapping_backstress_rate;
    double mapping_backstress_capacity_fraction;
    double mapping_center_limit_ratio;
    double mapping_core_radius_ratio;
    double mapping_ray_flow_weight;
    double mapping_fabric_dilation_rate;
    double mapping_fabric_recovery_rate;
    double mapping_fabric_saturation;
    double mapping_fabric_flow_weight;
    double mapping_fabric_dilatancy_weight;
    double mapping_directional_ratchet_weight;
    double mapping_high_bias_contraction_onset;
    double mapping_high_bias_contraction_full;
    double mapping_high_bias_contraction_gain;
    double mapping_memory_shear_minimum_ratio;
    double mapping_memory_shear_activation;
    int32_t mapping_corrector_iterations;
    double mapping_corrector_relaxation;
    double mapping_outer_tolerance;

    int32_t intermediate_bias_flow_enabled;
    double intermediate_bias_density_onset;
    double intermediate_bias_density_peak;
    double intermediate_bias_density_cutoff;
    double intermediate_bias_onset;
    double intermediate_bias_full;
    double intermediate_bias_cutoff_onset;
    double intermediate_bias_cutoff;
    double intermediate_compliance_peak;
    double intermediate_compliance_bell_gain;
    double intermediate_directional_balance;
    double intermediate_balance_bias_exponent;
    double intermediate_compliance_minimum;
    double intermediate_high_bias_onset;
    double intermediate_high_bias_full;
    double intermediate_high_bias_cutoff_onset;
    double intermediate_high_bias_cutoff;
    double intermediate_high_bias_compliance_peak;
    double intermediate_high_bias_compliance_bell_gain;
    double intermediate_high_bias_directional_balance;
    double intermediate_high_bias_amplitude_onset;
    double intermediate_high_bias_amplitude_full;
    double intermediate_first_reversal_amplitude_ratio;
    double intermediate_high_bias_phase_anchor_fraction;
    double intermediate_high_bias_phase_relaxation_multiplier;
    double intermediate_high_bias_phase_activation_reversals;
    double intermediate_high_bias_pre_reversal_phase_scale;

    /* Last completed half-cycle plastic multiplier needed for full activity. */
    double bias_reversible_volume_ep_ref;

    /* Derived once for the frozen material record; not a calibrated input. */
    double reference_relative_density_value;
} riva_ib_parameters_t;

/* base contains the first 38 logical handoff fields plus two Hercules lifecycle
 * flags.  The fields below complete the exact 63-field restart contract. */
typedef struct riva_ib_state_t {
    riva_state_t base;
    double phase_irreversible_volume;
    double phase_reversible_volume;
    double phase_potential_anchor;
    double phase_accumulation_lambda_anchor;
    double phase_accumulation_hardening_state;
    tensor_t unbiased_phase_direction;
    double loose_shear_lambda_anchor;
    double loose_shear_hardening_state;
    double loose_shear_gate_value;
    tensor_t mapping_anchor;
    tensor_t mapping_backstress;
    tensor_t mapping_directional_fabric;
    double mapping_gate_value;
    double mapping_capacity;
    double mapping_kinematic_denominator;
    double mapping_shear_modulus_ratio;
    double mapping_phase_contraction_scale;
    double mapping_outer_residual;
    int64_t mapping_stress_corrections;
    int64_t mapping_corrector_passes;
    int64_t mapping_monotone_caps;
    double initial_relative_density_value;
    double intermediate_low_gate_value;
    double intermediate_high_gate_base;
    /* Committed activity before the most recent reversal reset. */
    double ep_half_last;
} riva_ib_state_t;

/* The Python oracle uses NumPy's component-wise tensor division.  Keep that
 * operation distinct from reciprocal multiplication: the two are not
 * bit-identical and the difference is measurable after long low-confinement
 * histories. */
RIVA_IB_HD static inline tensor_t riva_ib_div(tensor_t value,double divisor)
{
    tensor_t result={value.xx/divisor,value.yy/divisor,value.zz/divisor,
                     value.xy/divisor,value.yz/divisor,value.xz/divisor};
    return result;
}

/* NumPy reduces the full row-major 3x3 product, including each symmetric
 * off-diagonal entry twice.  For nine binary64 entries its pairwise reducer
 * combines the first eight values as two four-value trees, then adds zz.
 * Preserve that operation tree for the research-oracle translation; the
 * mathematically equivalent compressed contraction has different rounding at
 * low-confinement zero crossings. */
RIVA_IB_HD static inline double riva_ib_add_rn(double a,double b)
{
#if defined(__CUDA_ARCH__)
    return __dadd_rn(a,b);
#else
    volatile double result=a+b;
    return result;
#endif
}

RIVA_IB_HD static inline double riva_ib_mul_rn(double a,double b)
{
#if defined(__CUDA_ARCH__)
    return __dmul_rn(a,b);
#else
    volatile double result=a*b;
    return result;
#endif
}

RIVA_IB_HD static inline double riva_ib_sub_rn(double a,double b)
{
#if defined(__CUDA_ARCH__)
    return __dsub_rn(a,b);
#else
    volatile double result=a-b;
    return result;
#endif
}

RIVA_IB_HD static inline double riva_ib_div_rn(double a,double b)
{
#if defined(__CUDA_ARCH__)
    return __ddiv_rn(a,b);
#else
    volatile double result=a/b;
    return result;
#endif
}

/* NumPy 2.0's ARM64 f64 tanh ufunc uses its FMA3/SVML interval polynomial.
 * The frozen oracle was generated by that path.  Preserve the [0,0.1875) and
 * [0.5,6) intervals exercised by the dense phase coordinates and unbiased
 * potential; system libm is oracle-equivalent in the remaining calibrated
 * ranges.  Coefficients and Horner order are from NumPy v2.0.2
 * loops_hyperbolic.dispatch.c.src. */
RIVA_IB_HD static inline double riva_ib_numpy_tanh(double x)
{
    const double ax=fabs(x);
    if (ax>=0.1875 && ax<0.5) return tanh(x);
    if (ax>=6.0) return tanh(x);
    double y,r;
    if (ax<0.1875) {
        y=ax;
        r=fma(0x1.5e67ab76a26e7p-12,y,-0x1.5d7e76dc56871p-10);
        r=fma(r,y,-0x1.a1306713a4f3ap-13);
        r=fma(r,y,0x1.e3be689423841p-9);
        r=fma(r,y,-0x1.dd99a221ed573p-16);
        r=fma(r,y,-0x1.22404577aa9ddp-7);
        r=fma(r,y,-0x1.8c4c1fd7852fep-21);
        r=fma(r,y,0x1.664f94e6ac14ep-6);
        r=fma(r,y,-0x1.51ca7f096011fp-28);
        r=fma(r,y,-0x1.ba1ba1990520bp-5);
        r=fma(r,y,-0x1.1ea19ddddb3b4p-37);
        r=fma(r,y,0x1.1111111112ab5p-3);
        r=fma(r,y,-0x1.6863ee44ed636p-49);
        r=fma(r,y,-0x1.5555555555555p-2);
        r=fma(r,y,-0x1.0b3ea3fdfaa19p-64);
        r=fma(r,y,0x1.0000000000000p+0);
        r=fma(r,y,0x0.0p+0);
    } else if (ax<0.75) {
        y=riva_ib_sub_rn(ax,0x1.4000000000000p-1);
        r=fma(0x1.69d8374520edap-4,y,0x1.246332a2fcba5p-5);
        r=fma(r,y,-0x1.4949d60113d63p-8);
        r=fma(r,y,-0x1.8fd89fe05e0d1p-10);
        r=fma(r,y,-0x1.045b9076cc487p-9);
        r=fma(r,y,0x1.5f87d903aaac8p-11);
        r=fma(r,y,0x1.5b29bb02cf69bp-8);
        r=fma(r,y,-0x1.129a092de747ap-7);
        r=fma(r,y,-0x1.31fe77c9c60afp-8);
        r=fma(r,y,0x1.e7f8380184b45p-6);
        r=fma(r,y,-0x1.87d27ccff4291p-6);
        r=fma(r,y,-0x1.c3c021789a786p-5);
        r=fma(r,y,0x1.1a686f6ab2533p-3);
        r=fma(r,y,-0x1.2426c751e48a2p-6);
        r=fma(r,y,-0x1.893b59c35c882p-2);
        r=fma(r,y,0x1.6284c3374f815p-1);
        r=fma(r,y,0x1.1bf47eabb8f95p-1);
    } else if (ax<1.0) {
        y=riva_ib_sub_rn(ax,0x1.c000000000000p-1);
        r=fma(-0x1.ded519f981716p-4,y,-0x1.29d851a896fcdp-4);
        r=fma(r,y,0x1.c9fd6200d0adep-8);
        r=fma(r,y,0x1.3f7af01d5af7ap-8);
        r=fma(r,y,0x1.085ee7e8ac170p-13);
        r=fma(r,y,-0x1.e104671036300p-10);
        r=fma(r,y,0x1.07df0f9f90c17p-9);
        r=fma(r,y,0x1.0c85b4d538746p-9);
        r=fma(r,y,-0x1.4a6046865ec7dp-7);
        r=fma(r,y,0x1.69543e7c420d4p-7);
        r=fma(r,y,0x1.b2ca62572b098p-7);
        r=fma(r,y,-0x1.e2196b7326859p-5);
        r=fma(r,y,0x1.f203c316ce730p-5);
        r=fma(r,y,0x1.4f152b2bad124p-4);
        r=fma(r,y,-0x1.6ba7cb7576538p-2);
        r=fma(r,y,0x1.02500a09f8d6ep-1);
        r=fma(r,y,0x1.686650b8c2015p-1);
    } else if (ax<1.5) {
        y=riva_ib_sub_rn(ax,0x1.4000000000000p+0);
        r=fma(-0x1.02d288b5b3371p-16,y,0x1.9065ae369b212p-18);
        r=fma(r,y,0x1.2cd40e0ad0a9fp-15);
        r=fma(r,y,-0x1.e40bdead17e6bp-14);
        r=fma(r,y,0x1.3524622610430p-13);
        r=fma(r,y,0x1.9bc98ddf0f340p-14);
        r=fma(r,y,-0x1.b852a6e0758d5p-11);
        r=fma(r,y,0x1.be9392199ec18p-10);
        r=fma(r,y,-0x1.ca3f1f2b9192bp-11);
        r=fma(r,y,-0x1.326bd4914222ap-8);
        r=fma(r,y,0x1.f1cf6c7f5b00ap-7);
        r=fma(r,y,-0x1.3a7a011ff8c2ap-6);
        r=fma(r,y,-0x1.9c7a02788557cp-7);
        r=fma(r,y,0x1.bba40cbef72bep-4);
        r=fma(r,y,-0x1.e7291743d7555p-3);
        r=fma(r,y,0x1.1f25131e3a8c0p-2);
        r=fma(r,y,0x1.b2523bb6b2deep-1);
    } else if (ax<2.0) {
        y=riva_ib_sub_rn(ax,0x1.c000000000000p+0);
        r=fma(0x1.290981209c1a6p-20,y,-0x1.8e1ba4c98a030p-20);
        r=fma(r,y,-0x1.58ab8e019f311p-23);
        r=fma(r,y,0x1.224cd6c4513e5p-17);
        r=fma(r,y,-0x1.f12a6626911b4p-16);
        r=fma(r,y,0x1.d4304bc9246e8p-15);
        r=fma(r,y,-0x1.078c63d1b8445p-15);
        r=fma(r,y,-0x1.a0c68a4489f10p-13);
        r=fma(r,y,0x1.c77dee0afd227p-11);
        r=fma(r,y,-0x1.fc15b0a9d98fap-10);
        r=fma(r,y,0x1.0379811e43dd5p-9);
        r=fma(r,y,0x1.e4709c7e8430ep-9);
        r=fma(r,y,-0x1.8157e26e0d541p-6);
        r=fma(r,y,0x1.01ba038be6a3dp-4);
        r=fma(r,y,-0x1.b6d85a01efb80p-4);
        r=fma(r,y,0x1.d22ca1c24a139p-4);
        r=fma(r,y,0x1.e1fbf97e33527p-1);
    } else if (ax<3.0) {
        y=riva_ib_sub_rn(ax,0x1.4000000000000p+1);
        r=fma(-0x1.67e924bf5ff6ep-26,y,0x1.ffd0766ad4016p-25);
        r=fma(r,y,-0x1.92fa6323b7cf8p-24);
        r=fma(r,y,-0x1.4b645e68eeaa3p-29);
        r=fma(r,y,0x1.b9008bca408afp-21);
        r=fma(r,y,-0x1.13c415f7b9d41p-18);
        r=fma(r,y,0x1.c12eadd55be7ap-17);
        r=fma(r,y,-0x1.0462601dc2faap-15);
        r=fma(r,y,0x1.4055bce68597ap-15);
        r=fma(r,y,0x1.4cffcfa69fbb6p-14);
        r=fma(r,y,-0x1.793826f78537ep-11);
        r=fma(r,y,0x1.7682afa611151p-9);
        r=fma(r,y,-0x1.07b55c1c7d278p-7);
        r=fma(r,y,0x1.16df44871efc8p-6);
        r=fma(r,y,-0x1.addae58c7141ap-6);
        r=fma(r,y,0x1.b3afe1fba5c76p-6);
        r=fma(r,y,0x1.f9258260a71c2p-1);
    } else if (ax<4.0) {
        y=riva_ib_sub_rn(ax,0x1.c000000000000p+1);
        r=fma(0x1.3f7f7de6b0eb6p-33,y,-0x1.c63c29f505f5bp-31);
        r=fma(r,y,0x1.df04d67876402p-29);
        r=fma(r,y,-0x1.abfebfb72bc83p-27);
        r=fma(r,y,0x1.34df71865f620p-25);
        r=fma(r,y,-0x1.22b8d9720cdb0p-24);
        r=fma(r,y,-0x1.fa600f593181bp-25);
        r=fma(r,y,0x1.7b6a219dea9f4p-20);
        r=fma(r,y,-0x1.2bf0cb4a71647p-17);
        r=fma(r,y,0x1.57e48e5b79d10p-15);
        r=fma(r,y,-0x1.405695e36240fp-13);
        r=fma(r,y,0x1.ef2ee77717cbfp-12);
        r=fma(r,y,-0x1.3a18d5843190fp-10);
        r=fma(r,y,0x1.3c6869dfc8870p-9);
        r=fma(r,y,-0x1.dc59376c7aa19p-9);
        r=fma(r,y,0x1.dd37d19b22b21p-9);
        r=fma(r,y,0x1.ff112c63a9077p-1);
    } else {
        y=riva_ib_sub_rn(ax,0x1.4000000000000p+2);
        r=fma(0x1.9ed18bae3ebbcp-41,y,-0x1.fab216b9e0e49p-40);
        r=fma(r,y,-0x1.5c72be95e4d2cp-38);
        r=fma(r,y,0x1.51c38f8695ed3p-34);
        r=fma(r,y,-0x1.5bb1bcf83ca73p-31);
        r=fma(r,y,0x1.22666d739bec0p-28);
        r=fma(r,y,-0x1.a3c935dce3f7dp-26);
        r=fma(r,y,0x1.0cbcc8d4c5c8ap-23);
        r=fma(r,y,-0x1.31eaafe73efd5p-21);
        r=fma(r,y,0x1.33b66d7d77264p-19);
        r=fma(r,y,-0x1.0e08de39ce756p-17);
        r=fma(r,y,0x1.95a4482f180b7p-16);
        r=fma(r,y,-0x1.fb6bbc89b1a5bp-15);
        r=fma(r,y,0x1.fb9aef915d828p-14);
        r=fma(r,y,-0x1.7cc5e74677410p-13);
        r=fma(r,y,0x1.7ccec13a9ef96p-13);
        r=fma(r,y,0x1.fff419668df11p-1);
    }
    return x<0.0?-r:r;
}

RIVA_IB_HD static inline tensor_t riva_ib_add_rn_tensor(tensor_t a,tensor_t b)
{
    tensor_t c={riva_ib_add_rn(a.xx,b.xx),riva_ib_add_rn(a.yy,b.yy),
        riva_ib_add_rn(a.zz,b.zz),riva_ib_add_rn(a.xy,b.xy),
        riva_ib_add_rn(a.yz,b.yz),riva_ib_add_rn(a.xz,b.xz)};
    return c;
}

RIVA_IB_HD static inline tensor_t riva_ib_sub_rn_tensor(tensor_t a,tensor_t b)
{
    tensor_t c={riva_ib_sub_rn(a.xx,b.xx),riva_ib_sub_rn(a.yy,b.yy),
        riva_ib_sub_rn(a.zz,b.zz),riva_ib_sub_rn(a.xy,b.xy),
        riva_ib_sub_rn(a.yz,b.yz),riva_ib_sub_rn(a.xz,b.xz)};
    return c;
}

RIVA_IB_HD static inline tensor_t riva_ib_scale_rn_tensor(tensor_t a,double x)
{
    tensor_t c={riva_ib_mul_rn(a.xx,x),riva_ib_mul_rn(a.yy,x),
        riva_ib_mul_rn(a.zz,x),riva_ib_mul_rn(a.xy,x),
        riva_ib_mul_rn(a.yz,x),riva_ib_mul_rn(a.xz,x)};
    return c;
}

RIVA_IB_HD static inline double riva_ib_trace_rn(tensor_t a)
{ return riva_ib_add_rn(riva_ib_add_rn(a.xx,a.yy),a.zz); }

RIVA_IB_HD static inline tensor_t riva_ib_dev_rn(tensor_t a)
{
    const double mean=riva_ib_div_rn(riva_ib_trace_rn(a),3.0);
    tensor_t isotropic={mean,mean,mean,0.0,0.0,0.0};
    return riva_ib_sub_rn_tensor(a,isotropic);
}

RIVA_IB_HD static inline double riva_ib_numpy_ddot(tensor_t a,tensor_t b)
{
    const double v0=riva_ib_mul_rn(a.xx,b.xx);
    const double v1=riva_ib_mul_rn(a.xy,b.xy);
    const double v2=riva_ib_mul_rn(a.xz,b.xz);
    const double v3=riva_ib_mul_rn(a.xy,b.xy);
    const double v4=riva_ib_mul_rn(a.yy,b.yy);
    const double v5=riva_ib_mul_rn(a.yz,b.yz);
    const double v6=riva_ib_mul_rn(a.xz,b.xz);
    const double v7=riva_ib_mul_rn(a.yz,b.yz);
    const double v8=riva_ib_mul_rn(a.zz,b.zz);
    const double first=riva_ib_add_rn(
        riva_ib_add_rn(v0,v1),riva_ib_add_rn(v2,v3));
    const double second=riva_ib_add_rn(
        riva_ib_add_rn(v4,v5),riva_ib_add_rn(v6,v7));
    return riva_ib_add_rn(riva_ib_add_rn(first,second),v8);
}

RIVA_IB_HD static inline double riva_ib_numpy_norm(tensor_t value)
{ return sqrt(riva_max(riva_ib_numpy_ddot(value,value),0.0)); }

/* Redirect contractions written below while leaving production RIVA-Sand's
 * already-defined helpers and ABI untouched. */
#define riva_ddot riva_ib_numpy_ddot
#define riva_norm riva_ib_numpy_norm
#define riva_add riva_ib_add_rn_tensor
#define riva_sub riva_ib_sub_rn_tensor
#define riva_scale riva_ib_scale_rn_tensor
#define riva_trace riva_ib_trace_rn
#define riva_dev riva_ib_dev_rn

RIVA_IB_HD static inline riva_ib_parameters_t
riva_ib_reference_parameters(double stress_scale)
{
    riva_ib_parameters_t p = {};
    p.base=riva_reference_parameters(stress_scale);
    p.base.p_residual=0.0;
    p.phase_transformation_enabled=1;
    p.phase_cyclic_shear_modulus_reduction=0.813;
    p.branch_compliance_enabled=1;
    p.branch_compliance_bias_reference=0.5384061785684372;
    p.branch_compliance_bias_exponent=1.0;
    p.branch_compliance_minimum=0.10;
    p.phase_compliance_peak=6.0; p.phase_compliance_bell_gain=3.5;
    p.phase_memory_hardening_reduction=2.8;
    p.phase_memory_reference_volume=5.0e-5;
    p.phase_bias_hardening_intercept=87.291302;
    p.phase_bias_hardening_exponent=3.2;
    p.phase_compliance_shape=1.55; p.phase_compliance_location=0.442;
    p.phase_compliance_half_width=0.787;
    p.branch_directional_balance=0.050;
    p.branch_balance_bias_exponent=4.0; p.branch_balance_bias_cap=3.9;
    p.phase_ratio=0.62; p.phase_width=0.50;
    p.phase_contraction_rate=4.0e-4;
    p.phase_reversible_scale=-3.0e-3;
    p.phase_reversible_mean_scale=6.5e-4;
    p.phase_reversible_relaxation_strain=0.0048;
    p.phase_reversible_relaxation_bias_exponent=1.6;
    p.phase_potential_anchor_fraction=1.0;
    p.phase_pressure_exponent=1.6; p.phase_wave_pressure_exponent=1.0;
    p.phase_mean_pressure_exponent=0.0;
    p.phase_amplitude_onset_ratio=0.55; p.phase_amplitude_full_ratio=0.90;
    p.phase_volume_density_onset=0.53; p.phase_volume_density_full=0.90;
    p.phase_volume_replacement_density_full=0.66;
    p.phase_intermediate_wave_multiplier=1.1333333333333333;
    p.phase_intermediate_mean_multiplier=10.76923076923077;
    p.phase_intermediate_relaxation_ratio=0.10416666666666667;
    p.phase_bias_exponent=2.0; p.phase_wave_bias_exponent=-0.25;
    p.phase_mean_bias_exponent=12.0;
    p.phase_accumulation_control_enabled=1;
    p.phase_accumulation_density_onset=0.53;
    p.phase_accumulation_density_peak=0.663;
    p.phase_accumulation_density_cutoff=0.82;
    p.phase_accumulation_reference_amplitude=0.43;
    p.phase_accumulation_memory_onset=0.030;
    p.phase_accumulation_memory_width=0.006;
    p.phase_accumulation_amplitude_exponent=1.0;
    p.phase_accumulation_gain_full_amplitude=0.65;
    p.phase_accumulation_high_amplitude_gain_ratio=0.6666666666666666;
    p.phase_accumulation_hardening_gain=3.0;
    p.loose_stabilization_enabled=1;
    p.loose_stabilization_density_full=0.50;
    p.loose_stabilization_density_cutoff=0.56;
    p.loose_stabilization_bias_onset=0.13;
    p.loose_stabilization_hardening_multiplier=1.40;
    p.loose_stabilization_contraction_scale=0.50;
    p.loose_phase_enabled=1; p.loose_phase_activity_scale=0.20;
    p.loose_phase_wave_activity_scale=0.45;
    p.loose_phase_mean_activity_scale=1.40;
    p.loose_phase_replacement_fraction=1.0;
    p.unbiased_phase_enabled=1; p.unbiased_phase_density_onset=0.82;
    p.unbiased_phase_density_full=0.89; p.unbiased_phase_bias_cutoff=0.08;
    p.unbiased_phase_activity_scale=0.70;
    p.unbiased_phase_potential_center=0.80;
    p.loose_shear_flow_enabled=1;
    p.loose_shear_cyclic_modulus_reduction=0.85;
    p.loose_shear_modulus_scale=0.75;
    p.loose_shear_early_hardening_multiplier=1.40;
    p.loose_shear_late_hardening_multiplier=3.50;
    p.loose_shear_hardening_memory_onset=0.035;
    p.loose_shear_hardening_memory_width=0.010;
    p.loose_shear_branch_compliance_gain=1.0;
    p.loose_shear_branch_compliance_onset=0.60;
    p.loose_shear_branch_compliance_full=0.95;
    p.loose_shear_ratchet_rate=4.0; p.loose_shear_ratchet_capacity=0.0172;
    p.loose_shear_ratchet_pressure_exponent=1.0;
    p.mapping_backstress_enabled=1; p.mapping_bias_onset=0.34;
    p.mapping_bias_full=0.54; p.mapping_density_full=0.55;
    p.mapping_density_cutoff=0.64; p.mapping_backstress_rate=11000.0;
    p.mapping_backstress_capacity_fraction=0.90;
    p.mapping_center_limit_ratio=0.97; p.mapping_core_radius_ratio=0.020;
    p.mapping_ray_flow_weight=0.16; p.mapping_fabric_dilation_rate=120.0;
    p.mapping_fabric_recovery_rate=30.0; p.mapping_fabric_saturation=0.75;
    p.mapping_fabric_flow_weight=0.10;
    p.mapping_fabric_dilatancy_weight=0.20;
    p.mapping_directional_ratchet_weight=0.020;
    p.mapping_high_bias_contraction_onset=0.62;
    p.mapping_high_bias_contraction_full=0.82;
    p.mapping_high_bias_contraction_gain=16.0;
    p.mapping_memory_shear_minimum_ratio=1.0;
    p.mapping_memory_shear_activation=3.0;
    p.mapping_corrector_iterations=5; p.mapping_corrector_relaxation=1.0;
    p.mapping_outer_tolerance=1.0e-10;
    p.intermediate_bias_flow_enabled=1;
    p.intermediate_bias_density_onset=0.60;
    p.intermediate_bias_density_peak=0.663;
    p.intermediate_bias_density_cutoff=0.76;
    p.intermediate_bias_onset=0.30; p.intermediate_bias_full=0.50;
    p.intermediate_bias_cutoff_onset=0.60;
    p.intermediate_bias_cutoff=0.66;
    p.intermediate_compliance_peak=1.0;
    p.intermediate_compliance_bell_gain=0.50;
    p.intermediate_directional_balance=0.0;
    p.intermediate_balance_bias_exponent=1.0;
    p.intermediate_compliance_minimum=0.20;
    p.intermediate_high_bias_onset=0.62;
    p.intermediate_high_bias_full=0.67;
    p.intermediate_high_bias_cutoff_onset=0.72;
    p.intermediate_high_bias_cutoff=0.78;
    p.intermediate_high_bias_compliance_peak=3.0;
    p.intermediate_high_bias_compliance_bell_gain=1.50;
    p.intermediate_high_bias_directional_balance=0.12;
    p.intermediate_high_bias_amplitude_onset=0.40;
    p.intermediate_high_bias_amplitude_full=0.45;
    p.intermediate_first_reversal_amplitude_ratio=0.40;
    p.intermediate_high_bias_phase_anchor_fraction=0.90;
    p.intermediate_high_bias_phase_relaxation_multiplier=0.25;
    p.intermediate_high_bias_phase_activation_reversals=6.0;
    p.intermediate_high_bias_pre_reversal_phase_scale=1.50;
    p.bias_reversible_volume_ep_ref=1.5e-5;
    p.reference_relative_density_value=
        (p.base.e_max-p.base.state_shakedown_reference_void_ratio)/
        (p.base.e_max-p.base.e_min);
    return p;
}

RIVA_IB_HD static inline double riva_ib_initial_density(
    const riva_ib_state_t *s)
{ return s->initial_relative_density_value; }

RIVA_IB_HD static inline double riva_ib_projected_bias(
    const riva_ib_state_t *s)
{ return riva_projected_bias(&s->base); }

RIVA_IB_HD static inline double riva_ib_density_weight(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s)
{ const double d=riva_max(riva_ib_initial_density(s)-
      p->reference_relative_density_value,0.0);
  return 1.0-exp(-pow(d/p->base.state_shakedown_anchor_width,2.0)); }

RIVA_IB_HD static inline double riva_ib_phase_amplitude_activity(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s)
{ const double ratio=s->base.cyclic_amplitude/p->base.cyclic_amplitude_reference;
  return riva_smoothstep((ratio-p->phase_amplitude_onset_ratio)/
      (p->phase_amplitude_full_ratio-p->phase_amplitude_onset_ratio)); }

RIVA_IB_HD static inline double riva_ib_loose_stabilization_gate_raw(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s)
{
    if (!p->loose_stabilization_enabled) return 0.0;
    const double dg=1.0-riva_smoothstep((riva_ib_initial_density(s)-
        p->loose_stabilization_density_full)/
        (p->loose_stabilization_density_cutoff-
         p->loose_stabilization_density_full));
    const double bg=riva_smoothstep(s->base.static_bias_index/
        p->loose_stabilization_bias_onset);
    return dg*bg;
}

RIVA_IB_HD static inline double riva_ib_loose_gate(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s)
{
    if (!p->loose_phase_enabled) return 0.0;
    return s->base.cyclic_phase_active?s->loose_shear_gate_value:
        riva_ib_loose_stabilization_gate_raw(p,s);
}

RIVA_IB_HD static inline double riva_ib_unbiased_gate(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s)
{
    if (!p->unbiased_phase_enabled) return 0.0;
    const double dg=riva_smoothstep((riva_ib_initial_density(s)-
        p->unbiased_phase_density_onset)/
        (p->unbiased_phase_density_full-p->unbiased_phase_density_onset));
    const double bg=1.0-riva_smoothstep(s->base.static_bias_index/
        p->unbiased_phase_bias_cutoff);
    return dg*bg;
}

RIVA_IB_HD static inline double riva_ib_mapping_gate_raw(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s)
{
    if (!p->mapping_backstress_enabled || !s->base.cyclic_phase_active) return 0.0;
    const double bg=riva_smoothstep((s->base.static_bias_index-p->mapping_bias_onset)/
        (p->mapping_bias_full-p->mapping_bias_onset));
    const double dg=1.0-riva_smoothstep((riva_ib_initial_density(s)-
        p->mapping_density_full)/(p->mapping_density_cutoff-p->mapping_density_full));
    return bg*dg;
}

RIVA_IB_HD static inline double riva_ib_mapping_gate(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s)
{ return s->base.cyclic_phase_active?s->mapping_gate_value:
    riva_ib_mapping_gate_raw(p,s); }

RIVA_IB_HD static inline void riva_ib_intermediate_gates_raw(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s,double *low,double *high)
{
    *low=0.0; *high=0.0;
    if (!p->intermediate_bias_flow_enabled) return;
    const double d=riva_ib_initial_density(s);
    double dg;
    if (d<=p->intermediate_bias_density_peak)
        dg=riva_smoothstep((d-p->intermediate_bias_density_onset)/
            (p->intermediate_bias_density_peak-p->intermediate_bias_density_onset));
    else dg=1.0-riva_smoothstep((d-p->intermediate_bias_density_peak)/
            (p->intermediate_bias_density_cutoff-p->intermediate_bias_density_peak));
    double lg=riva_smoothstep((s->base.static_bias_index-p->intermediate_bias_onset)/
        (p->intermediate_bias_full-p->intermediate_bias_onset));
    lg*=1.0-riva_smoothstep((s->base.static_bias_index-
        p->intermediate_bias_cutoff_onset)/
        (p->intermediate_bias_cutoff-p->intermediate_bias_cutoff_onset));
    double hg=riva_smoothstep((s->base.static_bias_index-
        p->intermediate_high_bias_onset)/
        (p->intermediate_high_bias_full-p->intermediate_high_bias_onset));
    hg*=1.0-riva_smoothstep((s->base.static_bias_index-
        p->intermediate_high_bias_cutoff_onset)/
        (p->intermediate_high_bias_cutoff-p->intermediate_high_bias_cutoff_onset));
    *low=dg*lg; *high=dg*hg;
}

RIVA_IB_HD static inline double riva_ib_intermediate_low_gate(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s)
{ return p->intermediate_bias_flow_enabled && s->base.cyclic_phase_active?
    s->intermediate_low_gate_value:0.0; }

RIVA_IB_HD static inline double riva_ib_intermediate_high_gate(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s)
{
    if (!p->intermediate_bias_flow_enabled || !s->base.cyclic_phase_active ||
        s->intermediate_high_gate_base==0.0) return 0.0;
    double amplitude=s->base.cyclic_amplitude;
    if (s->base.amplitude_reversals==1)
        amplitude*=p->intermediate_first_reversal_amplitude_ratio;
    return s->intermediate_high_gate_base*riva_smoothstep(
        (amplitude-p->intermediate_high_bias_amplitude_onset)/
        (p->intermediate_high_bias_amplitude_full-
         p->intermediate_high_bias_amplitude_onset));
}

RIVA_IB_HD static inline double riva_ib_branch_progress(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s)
{
    if (!p->phase_transformation_enabled || !p->branch_compliance_enabled ||
        !s->base.cyclic_phase_active || s->base.amplitude_reversals<1 ||
        s->base.cyclic_amplitude<=1.0e-14) return 0.0;
    const double excursion=riva_norm(riva_sub(riva_dev(s->base.stress),
        s->base.last_reversal_deviator));
    const double span=2.0*riva_max(s->base.pressure_anchor,p->base.p_min)*
        s->base.cyclic_amplitude;
    return riva_clip(excursion/riva_max(span,1.0e-14),0.0,1.0);
}

RIVA_IB_HD static inline double riva_ib_transformation_progress(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s)
{
    if (!p->phase_transformation_enabled) return 0.0;
    const double x=riva_ib_branch_progress(p,s);
    const double lo=p->phase_compliance_location-p->phase_compliance_half_width;
    const double hi=p->phase_compliance_location+p->phase_compliance_half_width;
    if (x<=lo) return 0.0;
    if (x>=hi) return 1.0;
    return pow(riva_smoothstep((x-lo)/riva_max(hi-lo,1.0e-12)),
               p->phase_compliance_shape);
}

RIVA_IB_HD static inline double riva_ib_phase_flow_bias_activity(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s,double exponent)
{
    const double bias=riva_ib_projected_bias(s);
    const double ref=p->branch_compliance_bias_reference;
    if (bias<=1.0e-14) return 0.0;
    const double gate=riva_smoothstep(bias/(0.25*ref));
    return riva_min(gate*pow(ref/bias,exponent),2.0);
}

RIVA_IB_HD static inline double riva_ib_phase_volume_density_weight(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s)
{ return riva_smoothstep((riva_ib_initial_density(s)-p->phase_volume_density_onset)/
    (p->phase_volume_density_full-p->phase_volume_density_onset)); }

RIVA_IB_HD static inline double riva_ib_phase_dense_bias_gate(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s)
{ const double onset=0.5*p->branch_compliance_bias_reference;
  return riva_ib_density_weight(p,s)*riva_smoothstep(
      riva_ib_projected_bias(s)/riva_max(onset,1.0e-14)); }

RIVA_IB_HD static inline double riva_ib_phase_volume_gate_parent(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s)
{ const double onset=0.5*p->branch_compliance_bias_reference;
  return riva_ib_phase_volume_density_weight(p,s)*riva_smoothstep(
      riva_ib_projected_bias(s)/riva_max(onset,1.0e-14)); }

RIVA_IB_HD static inline double riva_ib_phase_volume_replacement_gate_parent(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s)
{
    const double dg=riva_smoothstep((riva_ib_initial_density(s)-
        p->phase_volume_density_onset)/(p->phase_volume_replacement_density_full-
        p->phase_volume_density_onset));
    const double bg=riva_smoothstep(riva_ib_projected_bias(s)/
        riva_max(0.5*p->branch_compliance_bias_reference,1.0e-14));
    return dg*bg;
}

RIVA_IB_HD static inline double riva_ib_phase_volume_replacement_gate(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s)
{
    return riva_max(riva_ib_phase_volume_replacement_gate_parent(p,s),
        riva_max(riva_ib_unbiased_gate(p,s),
            riva_ib_loose_gate(p,s)*p->loose_phase_replacement_fraction));
}

RIVA_IB_HD static inline double riva_ib_phase_volume_gate(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s)
{
    const double loose=riva_ib_loose_gate(p,s);
    const double loose_activation=riva_max(riva_max(p->loose_phase_activity_scale,
        p->loose_phase_wave_activity_scale),riva_max(p->loose_phase_mean_activity_scale,
        p->loose_phase_replacement_fraction));
    double gate=riva_max(riva_ib_phase_volume_gate_parent(p,s),
        riva_max(riva_ib_unbiased_gate(p,s),loose*loose_activation));
    if (s->base.cyclic_phase_active && s->base.amplitude_reversals<1)
        gate=riva_max(gate,s->intermediate_high_gate_base*
            riva_ib_phase_volume_density_weight(p,s));
    return gate;
}

RIVA_IB_HD static inline double riva_ib_signed_phase_potential(
    const riva_ib_parameters_t *p,const riva_material_parameters_t *m,
    const riva_ib_state_t *s)
{
    tensor_t direction;
    if (riva_ib_unbiased_gate(p,s)>1.0e-14) {
        const double norm=riva_norm(s->unbiased_phase_direction);
        if (norm<=1.0e-14) return 0.0;
        direction=riva_ib_div(s->unbiased_phase_direction,norm);
        const double eta=riva_ddot(s->base.alpha,direction);
        double mb,md,xi; riva_surfaces(&p->base,m,riva_pressure(s->base.stress),
            s->base.void_ratio,&mb,&md,&xi); (void)mb; (void)xi;
        const double eta_pt=riva_ib_mul_rn(riva_ib_mul_rn(
            p->phase_ratio,sqrt(2.0/3.0)),md);
        const double smooth=riva_ib_sub_rn(sqrt(riva_ib_add_rn(
            riva_ib_mul_rn(eta,eta),1.0e-12)),1.0e-6);
        return riva_ib_sub_rn(riva_ib_numpy_tanh(riva_ib_div_rn(smooth,
            riva_max(riva_ib_mul_rn(p->phase_width,eta_pt),1.0e-12))),
            p->unbiased_phase_potential_center);
    }
    const double norm=riva_norm(s->base.static_bias_tensor);
    if (norm<=1.0e-14) return 0.0;
    direction=riva_ib_div(s->base.static_bias_tensor,norm);
    const double eta=riva_ddot(s->base.alpha,direction);
    double mb,md,xi; riva_surfaces(&p->base,m,riva_pressure(s->base.stress),
        s->base.void_ratio,&mb,&md,&xi); (void)mb; (void)xi;
    const double eta_pt=riva_ib_mul_rn(riva_ib_mul_rn(
        p->phase_ratio,sqrt(2.0/3.0)),md);
    return riva_ib_numpy_tanh(riva_ib_div_rn(eta,
        riva_max(riva_ib_mul_rn(p->phase_width,eta_pt),1.0e-12)));
}

RIVA_IB_HD static inline void riva_ib_phase_coordinates(
    const riva_ib_parameters_t *p,const riva_material_parameters_t *m,
    const riva_ib_state_t *s,double *loading_eta,double *eta_pt,double *signed_phase)
{
    tensor_t direction;
    const int unbiased=riva_ib_unbiased_gate(p,s)>1.0e-14;
    if (unbiased) direction=s->unbiased_phase_direction;
    else direction=s->base.cyclic_direction;
    const double norm=riva_norm(direction);
    if (norm<=1.0e-14) { *loading_eta=0.0; *eta_pt=0.0; *signed_phase=-1.0; return; }
    direction=riva_ib_div(direction,norm);
    direction=riva_scale(direction,unbiased?1.0:-1.0);
    *loading_eta=riva_ddot(s->base.alpha,direction);
    if (unbiased) *loading_eta=fabs(*loading_eta);
    double mb,md,xi; riva_surfaces(&p->base,m,riva_pressure(s->base.stress),
        s->base.void_ratio,&mb,&md,&xi); (void)mb; (void)xi;
    *eta_pt=riva_ib_mul_rn(riva_ib_mul_rn(
        p->phase_ratio,sqrt(2.0/3.0)),md);
    *signed_phase=riva_ib_numpy_tanh(riva_ib_div_rn(
        riva_ib_sub_rn(*loading_eta,*eta_pt),
        riva_max(riva_ib_mul_rn(p->phase_width,*eta_pt),1.0e-12)));
}

RIVA_IB_HD static inline double riva_ib_transformation_zone(
    const riva_ib_parameters_t *p,const riva_material_parameters_t *m,
    const riva_ib_state_t *s)
{ if (!p->phase_transformation_enabled) return 0.0;
  const double x=riva_ib_signed_phase_potential(p,m,s);
  return riva_max(riva_ib_sub_rn(1.0,riva_ib_mul_rn(x,x)),0.0); }

RIVA_IB_HD static inline double riva_ib_intermediate_phase_activation(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s)
{
    const double gate=s->intermediate_high_gate_base;
    if (gate==0.0) return 1.0;
    const double ramp=riva_smoothstep((double)s->base.amplitude_reversals/
        p->intermediate_high_bias_phase_activation_reversals);
    return 1.0+gate*(ramp-1.0);
}

RIVA_IB_HD static inline double riva_ib_phase_activity_parent(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s,
    double bias_exponent,double pressure_exponent)
{
    double biased=0.0;
    if (s->base.cyclic_phase_active && s->base.amplitude_reversals>=1) {
        const double pr=pow(riva_clip(riva_pressure(s->base.stress)/
            riva_max(s->base.pressure_anchor,p->base.p_min),0.0,1.0),pressure_exponent);
        biased=riva_ib_phase_volume_density_weight(p,s)*
            riva_ib_phase_flow_bias_activity(p,s,bias_exponent)*
            riva_ib_phase_amplitude_activity(p,s)*pr;
    }
    const double ug=riva_ib_unbiased_gate(p,s),lg=riva_ib_loose_gate(p,s);
    if (riva_max(ug,lg)<=1.0e-14 || !s->base.cyclic_phase_active) return biased;
    const double pr=pow(riva_clip(riva_pressure(s->base.stress)/
        riva_max(s->base.pressure_anchor,p->base.p_min),0.0,1.0),pressure_exponent);
    const double unbiased=p->unbiased_phase_activity_scale*ug*
        riva_ib_phase_amplitude_activity(p,s)*pr;
    double loose_scale=p->loose_phase_activity_scale;
    if (bias_exponent==p->phase_wave_bias_exponent)
        loose_scale=p->loose_phase_wave_activity_scale;
    else if (bias_exponent==p->phase_mean_bias_exponent)
        loose_scale=p->loose_phase_mean_activity_scale;
    return biased+unbiased+loose_scale*lg*riva_ib_phase_amplitude_activity(p,s)*pr;
}

RIVA_IB_HD static inline double riva_ib_pre_reversal_excursion(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s)
{
    const tensor_t static_dev=riva_add(s->base.geostatic_deviator,
        riva_scale(s->base.static_bias_tensor,s->base.pressure_anchor));
    return riva_norm(riva_sub(riva_dev(s->base.stress),static_dev))/
        riva_max(s->base.pressure_anchor,p->base.p_min);
}

RIVA_IB_HD static inline double riva_ib_phase_activity(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s,
    double bias_exponent,double pressure_exponent)
{
    double activity=riva_ib_phase_activity_parent(p,s,bias_exponent,pressure_exponent);
    if (activity<=0.0 && s->base.cyclic_phase_active &&
        s->intermediate_high_gate_base>1.0e-14 && s->base.amplitude_reversals<1) {
        const double ratio=riva_ib_pre_reversal_excursion(p,s)/
            p->base.cyclic_amplitude_reference;
        const double aa=riva_smoothstep((ratio-p->phase_amplitude_onset_ratio)/
            (p->phase_amplitude_full_ratio-p->phase_amplitude_onset_ratio));
        const double pr=pow(riva_clip(riva_pressure(s->base.stress)/
            riva_max(s->base.pressure_anchor,p->base.p_min),0.0,1.0),pressure_exponent);
        const double bias=riva_max(s->base.static_bias_index,0.0);
        const double ref=p->branch_compliance_bias_reference;
        const double bg=riva_smoothstep(bias/riva_max(0.25*ref,1.0e-14));
        const double ba=riva_min(bg*pow(ref/riva_max(bias,1.0e-14),bias_exponent),2.0);
        return s->intermediate_high_gate_base*
            p->intermediate_high_bias_pre_reversal_phase_scale*
            riva_ib_phase_volume_density_weight(p,s)*ba*aa*pr;
    }
    return activity*riva_ib_intermediate_phase_activation(p,s);
}

RIVA_IB_HD static inline double riva_ib_branch_compliance_multiplier(
    const riva_ib_parameters_t *p,const riva_material_parameters_t *m,
    const riva_ib_state_t *s)
{
    if (!p->phase_transformation_enabled) return 1.0;
    const double density=riva_ib_density_weight(p,s);
    const double bias=riva_ib_projected_bias(s);
    const double ba=bias<=1.0e-14?0.0:riva_min(pow(bias/
        p->branch_compliance_bias_reference,p->branch_compliance_bias_exponent),2.0);
    const double transition=riva_ib_transformation_progress(p,s);
    if (density<=1.0e-14 || ba<=1.0e-14 || transition<=1.0e-14) return 1.0;
    const double dn=riva_norm(s->base.cyclic_direction);
    const double projection=riva_ddot(s->base.static_bias_tensor,
        s->base.cyclic_direction)/riva_max(dn,1.0e-14);
    const double sign=projection>0.0?1.0:(projection<0.0?-1.0:0.0);
    const double signed_bias=sign*riva_min(pow(fabs(projection)/
        p->branch_compliance_bias_reference,p->branch_balance_bias_exponent),
        p->branch_balance_bias_cap);
    const double compliance=density*ba*(p->phase_compliance_peak*transition*
        exp(p->branch_directional_balance*signed_bias)+
        p->phase_compliance_bell_gain*riva_ib_transformation_zone(p,m,s));
    return riva_clip(1.0/(1.0+compliance),p->branch_compliance_minimum,1.0);
}

RIVA_IB_HD static inline double riva_ib_intermediate_branch_multiplier(
    const riva_ib_parameters_t *p,const riva_material_parameters_t *m,
    const riva_ib_state_t *s)
{
    const double low=riva_ib_intermediate_low_gate(p,s);
    const double high=riva_ib_intermediate_high_gate(p,s);
    if (riva_max(low,high)<=1.0e-14) return 1.0;
    const double transition=riva_ib_transformation_progress(p,s);
    if (transition<=1.0e-14) return 1.0;
    const double dn=riva_norm(s->base.cyclic_direction);
    const double projection=riva_ddot(s->base.static_bias_tensor,
        s->base.cyclic_direction)/riva_max(dn,1.0e-14);
    const double sign=projection>0.0?1.0:(projection<0.0?-1.0:0.0);
    const double sb=sign*pow(fabs(projection)/
        riva_max(p->branch_compliance_bias_reference,1.0e-14),
        p->intermediate_balance_bias_exponent);
    const double zone=riva_ib_transformation_zone(p,m,s);
    const double compliance=low*(p->intermediate_compliance_peak*transition*
        exp(p->intermediate_directional_balance*sb)+
        p->intermediate_compliance_bell_gain*zone)+
        high*(p->intermediate_high_bias_compliance_peak*transition*
        exp(p->intermediate_high_bias_directional_balance*sb)+
        p->intermediate_high_bias_compliance_bell_gain*zone);
    return riva_clip(1.0/(1.0+compliance),p->intermediate_compliance_minimum,1.0);
}

RIVA_IB_HD static inline void riva_ib_cyclic_flow(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s,double pressure,
    double *shear,double *hardening)
{
    const riva_parameters_t *b=&p->base;
    if (!b->cyclic_flow_correction_enabled ||
        s->base.amplitude_reversals<b->cyclic_flow_minimum_reversals) {
        *shear=1.0; *hardening=1.0; return;
    }
    const double ratio=riva_clip(riva_ib_div_rn(pressure,riva_max(
        s->base.pressure_anchor,b->p_min)),0.0,1.0);
    const double activity=pow(ratio,b->cyclic_flow_pressure_exponent);
    *shear=riva_ib_sub_rn(1.0,riva_ib_mul_rn(
        b->cyclic_shear_modulus_reduction,activity));
    *hardening=riva_ib_add_rn(1.0,riva_ib_mul_rn(
        b->cyclic_hardening_boost,activity));
    if (p->phase_transformation_enabled) {
        const double phase=riva_ib_sub_rn(1.0,riva_ib_mul_rn(
            p->phase_cyclic_shear_modulus_reduction,activity));
        const double gate=riva_ib_phase_dense_bias_gate(p,s);
        *shear=riva_ib_add_rn(*shear,riva_ib_mul_rn(
            gate,riva_ib_sub_rn(phase,*shear)));
    }
    if (p->loose_shear_flow_enabled) {
        const double gate=riva_ib_loose_gate(p,s);
        if (gate>1.0e-14) {
            const double loose=riva_ib_sub_rn(1.0,riva_ib_mul_rn(
                p->loose_shear_cyclic_modulus_reduction,activity));
            *shear=riva_ib_add_rn(*shear,riva_ib_mul_rn(
                gate,riva_ib_sub_rn(loose,*shear)));
        }
    }
}

RIVA_IB_HD static inline double riva_ib_raw_shakedown(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s,double pressure)
{
    const riva_parameters_t *b=&p->base;
    if (!b->state_shakedown_enabled || !s->base.cyclic_phase_active ||
        s->base.amplitude_reversals<b->state_shakedown_minimum_reversals) return 1.0;
    const double bias=riva_ib_projected_bias(s); if (bias<=1.0e-14) return 1.0;
    const double distance=riva_ib_initial_density(s)-p->reference_relative_density_value;
    const double state_term=exp(b->state_shakedown_state_sensitivity*distance);
    const double release=1.0-exp(-pow(distance/b->state_shakedown_anchor_width,2.0));
    const double bias_term=pow(bias/b->state_shakedown_bias_reference,
        b->state_shakedown_bias_exponent);
    const double pressure_term=pow(riva_clip(pressure/riva_max(
        s->base.pressure_anchor,b->p_min),0.0,1.0),b->state_shakedown_pressure_exponent);
    return riva_clip(1.0+b->state_shakedown_scale*release*state_term*bias_term*
        pressure_term,1.0,b->state_shakedown_factor_maximum);
}

RIVA_IB_HD static inline double riva_ib_shakedown(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s,double pressure)
{
    const riva_parameters_t *b=&p->base;
    const double raw=riva_ib_raw_shakedown(p,s,pressure);
    const double excess=riva_max(s->base.cyclic_amplitude-
        b->state_shakedown_dense_hardening_amplitude_onset,0.0);
    const double dense=b->state_shakedown_dense_hardening_scale*
        exp(-b->state_shakedown_dense_hardening_amplitude_decay*excess);
    const double branch=1.0+riva_ib_density_weight(p,s)*(dense-1.0);
    return riva_clip(1.0+(raw-1.0)*branch,1.0,b->state_shakedown_factor_maximum);
}

RIVA_IB_HD static inline double riva_ib_shakedown_shear(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s,double pressure)
{
    const riva_parameters_t *b=&p->base;
    const double raw=riva_ib_raw_shakedown(p,s,pressure);
    const double distance=riva_max(riva_ib_initial_density(s)-
        p->reference_relative_density_value,0.0);
    const double attenuation=exp(-b->state_shakedown_dense_compliance_decay*
        distance*distance);
    const double compliance=1.0+(raw-1.0)*attenuation;
    return riva_clip(pow(compliance,-b->state_shakedown_compliance_exponent),
        b->state_shakedown_shear_multiplier_minimum,1.0);
}

RIVA_IB_HD static inline double riva_ib_phase_bias_hardening(
    const riva_ib_parameters_t *p,const riva_material_parameters_t *m,
    const riva_ib_state_t *s)
{
    const riva_parameters_t *b=&p->base;
    const double frozen=riva_bias_hardening(b,&s->base);
    if (!p->phase_transformation_enabled || !b->static_bias_enabled ||
        !s->base.cyclic_phase_active ||
        s->base.amplitude_reversals<b->bias_minimum_reversals) return frozen;
    const double bias=riva_ib_projected_bias(s);
    const double amplitude=b->bias_amplitude_ratio*s->base.cyclic_amplitude;
    const double margin=riva_max(bias-amplitude,0.0);
    const double crossing=riva_max(amplitude/riva_max(bias,1.0e-14)-1.0,0.0);
    const double phase=p->phase_bias_hardening_intercept*
        pow(bias,p->phase_bias_hardening_exponent)*exp(-b->bias_crossing_decay*crossing)+
        b->bias_hardening_scale*pow(margin,b->bias_margin_exponent);
    const double gate=riva_ib_phase_dense_bias_gate(p,s);
    double boost=frozen+gate*(phase-frozen);
    const double high=riva_smoothstep((bias/p->branch_compliance_bias_reference-1.0)/0.50);
    const double memory=1.0-exp(-riva_max(-s->phase_irreversible_volume,0.0)/
        p->phase_memory_reference_volume);
    const double log_mult=-p->phase_memory_hardening_reduction*
        riva_ib_density_weight(p,s)*riva_ib_transformation_zone(p,m,s)*high*memory;
    return boost*exp(log_mult);
}

RIVA_IB_HD static inline double riva_ib_phase_accumulation_memory(
    const riva_ib_state_t *s)
{ return riva_max(s->base.lambda_total-s->phase_accumulation_lambda_anchor,0.0); }

RIVA_IB_HD static inline double riva_ib_phase_accumulation_density_gate(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s)
{
    const double d=riva_ib_initial_density(s);
    if (d<=p->phase_accumulation_density_peak)
        return riva_smoothstep((d-p->phase_accumulation_density_onset)/
            (p->phase_accumulation_density_peak-p->phase_accumulation_density_onset));
    return 1.0-riva_smoothstep((d-p->phase_accumulation_density_peak)/
        (p->phase_accumulation_density_cutoff-p->phase_accumulation_density_peak));
}

RIVA_IB_HD static inline double riva_ib_phase_accumulation_target_activity(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s)
{
    if (!p->phase_transformation_enabled || !p->phase_accumulation_control_enabled ||
        !s->base.cyclic_phase_active) return 0.0;
    const double amplitude=riva_max(s->base.cyclic_amplitude,
        p->phase_accumulation_reference_amplitude);
    const double onset=p->phase_accumulation_memory_onset*pow(
        p->phase_accumulation_reference_amplitude/amplitude,
        p->phase_accumulation_amplitude_exponent);
    const double transition=riva_smoothstep((riva_ib_phase_accumulation_memory(s)-onset)/
        p->phase_accumulation_memory_width);
    const double biasgate=riva_smoothstep(riva_ib_projected_bias(s)/
        riva_max(0.5*p->branch_compliance_bias_reference,1.0e-14));
    return riva_ib_phase_accumulation_density_gate(p,s)*biasgate*transition;
}

RIVA_IB_HD static inline double riva_ib_phase_accumulation_target_hardening(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s)
{
    const double position=(s->base.cyclic_amplitude-
        p->phase_accumulation_reference_amplitude)/
        (p->phase_accumulation_gain_full_amplitude-
         p->phase_accumulation_reference_amplitude);
    const double gain=1.0+riva_smoothstep(position)*
        (p->phase_accumulation_high_amplitude_gain_ratio-1.0);
    return gain*riva_ib_phase_accumulation_target_activity(p,s);
}

RIVA_IB_HD static inline double riva_ib_loose_shear_memory(const riva_ib_state_t *s)
{ return riva_max(s->base.lambda_total-s->loose_shear_lambda_anchor,0.0); }

RIVA_IB_HD static inline double riva_ib_loose_target_hardening(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s)
{
    if (!p->loose_shear_flow_enabled || !s->base.cyclic_phase_active) return 0.0;
    return riva_ib_loose_gate(p,s)*riva_smoothstep((riva_ib_loose_shear_memory(s)-
        p->loose_shear_hardening_memory_onset)/p->loose_shear_hardening_memory_width);
}

RIVA_IB_HD static inline double riva_ib_loose_hardening_multiplier(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s)
{
    const double activity=riva_max(s->loose_shear_hardening_state,0.0);
    return p->loose_shear_early_hardening_multiplier+activity*
        (p->loose_shear_late_hardening_multiplier-
         p->loose_shear_early_hardening_multiplier);
}

RIVA_IB_HD static inline void riva_ib_moduli_for_state(
    const riva_ib_parameters_t *p,const riva_material_parameters_t *m,
    const riva_ib_state_t *s,double pressure,double *shear,double *bulk)
{
    riva_moduli(&p->base,m,pressure,shear,bulk);
    double sf,hf; riva_ib_cyclic_flow(p,s,pressure,&sf,&hf); (void)hf;
    *shear=riva_ib_mul_rn(*shear,sf);
    *shear=riva_ib_mul_rn(*shear,riva_ib_shakedown_shear(p,s,pressure));
    *shear=riva_ib_mul_rn(*shear,
        riva_ib_branch_compliance_multiplier(p,m,s));
    if (p->loose_shear_flow_enabled && s->base.cyclic_phase_active) {
        const double gate=riva_ib_loose_gate(p,s);
        if (gate>1.0e-14) {
            const double scale=1.0+gate*(p->loose_shear_modulus_scale-1.0);
            const double transition=riva_smoothstep((riva_ib_branch_progress(p,s)-
                p->loose_shear_branch_compliance_onset)/
                (p->loose_shear_branch_compliance_full-
                 p->loose_shear_branch_compliance_onset));
            *shear=riva_ib_mul_rn(*shear,riva_ib_div_rn(scale,
                riva_ib_add_rn(1.0,riva_ib_mul_rn(riva_ib_mul_rn(
                    gate,p->loose_shear_branch_compliance_gain),transition))));
        }
    }
    *shear=riva_ib_mul_rn(*shear,
        riva_ib_intermediate_branch_multiplier(p,m,s));
}

RIVA_IB_HD static inline double riva_ib_hardening_for_state(
    const riva_ib_parameters_t *p,const riva_material_parameters_t *m,
    const riva_ib_state_t *s,double pressure)
{
    const riva_parameters_t *b=&p->base;
    double hardening=m->h*pow(riva_max(pressure,b->p_min)/b->p_ref,-b->q_H);
    double sf,hf; riva_ib_cyclic_flow(p,s,pressure,&sf,&hf); (void)sf;
    hardening*=hf;
    const double boost=riva_ib_phase_bias_hardening(p,m,s);
    const double pr=riva_clip(pressure/riva_max(s->base.pressure_anchor,b->p_min),0.0,1.0);
    const double confinement=pow(riva_clip(b->bias_reference_pressure/
        riva_max(s->base.pressure_anchor,b->p_min),0.25,4.0),b->bias_confinement_exponent);
    hardening*=1.0+boost*confinement*pow(pr,b->bias_pressure_exponent);
    hardening*=riva_ib_shakedown(p,s,pressure);
    hardening*=1.0+p->phase_accumulation_hardening_gain*
        riva_max(s->phase_accumulation_hardening_state,0.0);
    const double loose_gate=riva_ib_loose_gate(p,s);
    hardening*=1.0+loose_gate*(p->loose_stabilization_hardening_multiplier-1.0);
    if (p->loose_shear_flow_enabled && loose_gate>1.0e-14 &&
        s->base.cyclic_phase_active) {
        const double parent=1.0+loose_gate*(p->loose_stabilization_hardening_multiplier-1.0);
        const double target=1.0+loose_gate*(riva_ib_loose_hardening_multiplier(p,s)-1.0);
        hardening=hardening/riva_max(parent,1.0e-14)*target;
    }
    return hardening;
}

RIVA_IB_HD static inline double riva_ib_irreversible_factor(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s)
{
    const riva_parameters_t *b=&p->base;
    double factor=s->base.amplitude_factor*s->base.state_contraction_factor;
    factor*=1.0+b->bias_contraction_scale*pow(riva_ib_projected_bias(s),
        b->bias_contraction_exponent);
    const double gate=riva_ib_loose_stabilization_gate_raw(p,s);
    return factor*(1.0-gate*(1.0-p->loose_stabilization_contraction_scale));
}

RIVA_IB_HD static inline int riva_ib_initialize_material(
    const riva_ib_parameters_t *p,const riva_material_parameters_t *m,
    tensor_t stress,double void_ratio,riva_ib_state_t *state)
{
    if (!p || !state || !riva_initialize_material(&p->base,m,stress,void_ratio,
                                                   &state->base)) return 0;
    state->base.alpha=riva_ib_div(riva_dev(stress),riva_pressure(stress));
    /* riva_initialize_material zeroes the base.  Explicitly clear the
     * extension too so reused integration-point storage is deterministic. */
    state->phase_irreversible_volume=0.0;
    state->phase_reversible_volume=0.0;
    state->phase_potential_anchor=0.0;
    state->phase_accumulation_lambda_anchor=0.0;
    state->phase_accumulation_hardening_state=0.0;
    state->unbiased_phase_direction=riva_zero();
    state->loose_shear_lambda_anchor=0.0;
    state->loose_shear_hardening_state=0.0;
    state->loose_shear_gate_value=0.0;
    state->mapping_anchor=riva_zero();
    state->mapping_backstress=riva_zero();
    state->mapping_directional_fabric=riva_zero();
    state->mapping_gate_value=0.0;
    state->mapping_capacity=0.0;
    state->mapping_kinematic_denominator=0.0;
    state->mapping_shear_modulus_ratio=1.0;
    state->mapping_phase_contraction_scale=1.0;
    state->mapping_outer_residual=0.0;
    state->mapping_stress_corrections=0;
    state->mapping_corrector_passes=0;
    state->mapping_monotone_caps=0;
    state->initial_relative_density_value=riva_initial_relative_density(
        &p->base,m,&state->base);
    state->intermediate_low_gate_value=0.0;
    state->intermediate_high_gate_base=0.0;
    state->ep_half_last=0.0;
    return 1;
}

RIVA_IB_HD static inline int riva_ib_initialize(
    const riva_ib_parameters_t *p,tensor_t stress,double void_ratio,
    riva_ib_state_t *state)
{
    const riva_material_parameters_t m=riva_reference_material_parameters(&p->base);
    return riva_ib_initialize_material(p,&m,stress,void_ratio,state);
}

RIVA_IB_HD static inline int riva_ib_admit_geostatic_state(
    const riva_ib_parameters_t *p,const riva_material_parameters_t *m,
    riva_ib_state_t *state,int32_t *admitted_out)
{ return riva_admit_geostatic_state(&p->base,m,&state->base,admitted_out); }

RIVA_IB_HD static inline int riva_ib_begin_dynamic_phase(
    const riva_ib_parameters_t *p,const riva_material_parameters_t *m,
    const tensor_t *reference_stress,riva_ib_state_t *s)
{
    if (!p || !m || !s || !s->base.initialized) return 0;
    /* Restarting an already-active point must never re-anchor its memories. */
    if (s->base.cyclic_phase_active) return 1;
    const tensor_t current=riva_dev(s->base.stress);
    const tensor_t reference=reference_stress?riva_dev(*reference_stress):
        s->base.geostatic_deviator;
    const tensor_t bias=riva_sub(current,reference);
    const double anchor=riva_max(s->base.pressure_anchor,p->base.p_min);
    s->base.static_bias_tensor=riva_ib_div(bias,anchor);
    s->base.static_bias_index=riva_norm(bias)/anchor;
    s->base.cyclic_direction=riva_zero();
    s->base.cyclic_phase_active=1;
    s->base.bias_ratchet_strain=0.0;
    s->base.bias_reversible_volume=0.0;
    s->base.eps_v_total=s->base.physical_eps_v_total;
    s->base.last_host_deviatoric_strain_direction=riva_zero();
    s->ep_half_last=0.0;

    s->phase_irreversible_volume=0.0;
    s->phase_reversible_volume=0.0;
    s->phase_potential_anchor=p->phase_potential_anchor_fraction*
        riva_ib_signed_phase_potential(p,m,s);
    s->phase_accumulation_lambda_anchor=s->base.lambda_total;
    s->phase_accumulation_hardening_state=0.0;
    s->unbiased_phase_direction=riva_zero();
    s->loose_shear_lambda_anchor=s->base.lambda_total;
    s->loose_shear_hardening_state=0.0;
    s->loose_shear_gate_value=riva_ib_loose_stabilization_gate_raw(p,s);
    s->mapping_anchor=s->base.alpha;
    s->mapping_backstress=riva_zero();
    s->mapping_directional_fabric=riva_zero();
    s->mapping_gate_value=riva_ib_mapping_gate_raw(p,s);
    s->mapping_capacity=0.0;
    s->mapping_kinematic_denominator=0.0;
    s->mapping_shear_modulus_ratio=1.0;
    s->mapping_phase_contraction_scale=1.0;
    s->mapping_outer_residual=0.0;
    s->mapping_stress_corrections=0;
    s->mapping_corrector_passes=0;
    s->mapping_monotone_caps=0;
    riva_ib_intermediate_gates_raw(p,s,&s->intermediate_low_gate_value,
                                    &s->intermediate_high_gate_base);
    s->phase_potential_anchor*=1.0+s->intermediate_high_gate_base*
        (p->intermediate_high_bias_phase_anchor_fraction-1.0);
    return 1;
}

RIVA_IB_HD static inline void riva_ib_dilatancy(
    const riva_parameters_t *p,const riva_material_parameters_t *m,
    tensor_t alpha,tensor_t normal,double beta,tensor_t fabric,
    double pressure,double void_ratio,double eps_v_ir,double eps_v_re,
    double *D_ir,double *D_re)
{
    double mb,md,xi;
    riva_surfaces(p,m,pressure,void_ratio,&mb,&md,&xi); (void)xi;
    const double eta=riva_ddot(alpha,normal);
    const double sqrt23=sqrt(2.0/3.0);
    const double eta_d=riva_ib_mul_rn(sqrt23,md);
    *D_ir=0.0;
    if (p->irreversible_enabled) {
        const double distance=riva_max(riva_ib_sub_rn(eta_d,eta),0.0);
        const double amplitude=riva_ib_mul_rn(m->zeta,riva_ib_add_rn(
            1.0,riva_max(riva_ddot(fabric,normal),0.0)));
        double gate=1.0;
        if (p->contraction_gate_enabled && p->C_D>0.0) {
            const double reversal=riva_ib_div_rn(riva_ib_sub_rn(
                riva_ib_mul_rn(sqrt23,mb),eta),
                riva_max(beta,1.0e-12));
            const double c_in=riva_ib_mul_rn(riva_ib_mul_rn(
                p->C_in_ratio,sqrt23),mb);
            const double shifted=riva_ib_add_rn(reversal,c_in);
            gate=riva_ib_div_rn(riva_ib_mul_rn(shifted,shifted),
                riva_ib_add_rn(riva_ib_mul_rn(reversal,reversal),p->C_D));
        }
        const double decay=exp(riva_ib_mul_rn(
            -p->irreversible_decay,riva_max(-eps_v_ir,0.0)));
        *D_ir=riva_ib_mul_rn(riva_ib_mul_rn(
            riva_ib_mul_rn(amplitude,gate),distance),decay);
    }
    *D_re=0.0;
    if (p->reversible_enabled) {
        if (eta>eta_d)
            *D_re=riva_ib_mul_rn(riva_ib_mul_rn(
                -m->zeta,p->reversible_ratio),riva_ib_sub_rn(eta,eta_d));
        else if (eps_v_re>0.0)
            *D_re=riva_min(riva_ib_mul_rn(
                               p->reversible_release_rate,eps_v_re),
                           riva_max(eta_d,1.0e-12));
    }
}

RIVA_IB_HD static inline void riva_ib_register_reversal(
    const riva_parameters_t *p,const riva_state_t *old,tensor_t *alpha0,
    tensor_t *alpha01,double *ep_eq,double *beta,int64_t *reversals)
{
    const double ratio=(*ep_eq)/p->overshoot_strain;
    const double modifier=ratio>1.0?0.0:1.0-ratio;
    const tensor_t previous=*alpha0;
    *alpha0=riva_norm(previous)==0.0?old->alpha:
        riva_add(riva_scale(*alpha01,modifier),
                 riva_scale(old->alpha,1.0-modifier));
    *alpha01=previous;
    *ep_eq=0.0;
    *beta=p->beta0;
    *reversals=old->reversals+1;
}

RIVA_IB_HD static inline tensor_t riva_ib_host_direction(
    const riva_parameters_t *p,tensor_t deps,int *valid)
{
    const tensor_t dev=riva_dev(deps);
    const double norm=riva_norm(dev);
    *valid=norm>p->reversal_strain_deadband;
    return *valid?riva_ib_div(dev,norm):riva_zero();
}

RIVA_IB_HD static inline int riva_ib_host_reversal(
    const riva_parameters_t *p,const riva_state_t *s,tensor_t direction,
    int direction_valid)
{
    if (!direction_valid || !s->cyclic_phase_active) return 0;
    const double previous_norm=riva_norm(
        s->last_host_deviatoric_strain_direction);
    if (previous_norm<=1.0e-14) return 0;
    const double cosine=riva_ddot(
        direction,s->last_host_deviatoric_strain_direction)/previous_norm;
    if (cosine>p->reversal_direction_cosine) return 0;
    const double excursion=riva_norm(riva_sub(
        riva_dev(s->stress),s->last_reversal_deviator));
    return excursion>=p->reversal_stress_deadband_ratio*
        riva_max(s->pressure_anchor,p->p_min);
}

RIVA_IB_HD static inline riva_ib_state_t riva_ib_base_backbone(
    const riva_ib_parameters_t *p,const riva_material_parameters_t *m,
    const riva_ib_state_t *old,tensor_t deps,int force_reversal,
    int allow_legacy_reversal)
{
    const riva_parameters_t *b=&p->base;
    riva_ib_state_t state=*old;
    const tensor_t s_t=riva_dev(old->base.stress);
    const double p_t=riva_max(riva_pressure(old->base.stress),b->p_min);
    double shear,bulk; riva_ib_moduli_for_state(p,m,old,p_t,&shear,&bulk);
    const double h_eff=riva_ib_hardening_for_state(p,m,old,p_t);
    const double threshold=riva_confining_strain(b,m,b->p_min,
        old->base.pressure_anchor);
    const int floor_active=b->compatibility_enabled &&
        old->base.eps_v_confining>=threshold;
    const double coupling_bulk=floor_active?0.0:bulk;
    const double deps_v=riva_trace(deps);
    const tensor_t deps_d=riva_dev(deps);
    const tensor_t s_trial=riva_add(s_t,riva_scale(
        deps_d,riva_ib_mul_rn(2.0,shear)));
    double p_trial;
    if (b->compatibility_enabled) {
        int hit=0; p_trial=riva_pressure_from_confining(b,m,
            old->base.eps_v_confining+deps_v,old->base.pressure_anchor,&hit);
    } else p_trial=riva_max(p_t-bulk*deps_v,b->p_min);
    const tensor_t alpha_trial=riva_ib_div(s_trial,riva_max(p_trial,b->p_min));
    const tensor_t delta_alpha=riva_sub(alpha_trial,old->base.alpha);
    tensor_t alpha0=old->base.alpha0,alpha01=old->base.alpha01;
    double beta_t=old->base.beta,ep_eq=old->base.ep_eq_since_reversal;
    int64_t reversals=old->base.reversals;
    const double legacy=riva_ddot(riva_sub(alpha_trial,alpha0),delta_alpha);
    if (force_reversal || (allow_legacy_reversal && legacy<0.0))
        riva_ib_register_reversal(b,&old->base,&alpha0,&alpha01,&ep_eq,
                                  &beta_t,&reversals);

    const double alpha_n=riva_ddot(old->base.alpha,old->base.n);
    const double hardening=riva_ib_mul_rn(riva_ib_mul_rn(p_t,h_eff),
        pow(riva_max(beta_t,1.0e-12),m->m));
    double denominator=riva_ib_sub_rn(riva_ib_add_rn(
        riva_ib_mul_rn(2.0,shear),riva_ib_mul_rn(2.0/3.0,hardening)),
        riva_ib_mul_rn(riva_ib_mul_rn(coupling_bulk,old->base.D),alpha_n));
    const double denom_floor=riva_ib_mul_rn(riva_ib_mul_rn(
        b->denominator_floor_ratio,2.0),shear);
    int64_t denom_hits=old->base.denominator_floor_hits;
    if (denominator<denom_floor) { denominator=denom_floor; denom_hits++; }
    const double numerator=riva_ib_add_rn(riva_ib_mul_rn(riva_ib_mul_rn(
        2.0,shear),riva_ddot(deps_d,old->base.n)),riva_ib_mul_rn(
        riva_ib_mul_rn(coupling_bulk,deps_v),alpha_n));
    const double dl=riva_max(0.0,riva_ib_div_rn(numerator,denominator));
    tensor_t s=riva_add(s_t,riva_scale(riva_sub(deps_d,
        riva_scale(old->base.n,dl)),riva_ib_mul_rn(2.0,shear)));
    const double d_ir=riva_max(old->base.D_ir,0.0);
    double d_re=old->base.D_re;
    if (dl>0.0 && d_re>0.0)
        d_re=riva_min(d_re,old->base.eps_v_reversible/dl);
    const double eps_ir=riva_ib_sub_rn(old->base.eps_v_irreversible,
        riva_ib_mul_rn(dl,d_ir));
    const double eps_re=riva_max(riva_ib_sub_rn(
        old->base.eps_v_reversible,riva_ib_mul_rn(dl,d_re)),0.0);
    const double eps_total=riva_ib_add_rn(old->base.eps_v_total,deps_v);
    const double eps_confining=riva_ib_sub_rn(
        riva_ib_sub_rn(eps_total,eps_ir),eps_re);
    double pressure; int at_floor=0;
    if (b->compatibility_enabled)
        pressure=riva_pressure_from_confining(b,m,eps_confining,
            old->base.pressure_anchor,&at_floor);
    else {
        pressure=riva_max(p_t-bulk*(deps_v+dl*(d_ir+d_re)),b->p_min);
        at_floor=pressure<=b->p_min;
    }
    ep_eq=riva_ib_add_rn(ep_eq,dl);
    const double void_ratio=riva_ib_add_rn(old->base.void_ratio,
        riva_ib_mul_rn(riva_ib_add_rn(1.0,old->base.void_ratio),deps_v));
    tensor_t fabric=old->base.fabric;
    if (b->fabric_enabled && dl>0.0) {
        const tensor_t term=riva_add(riva_scale(old->base.n,
            riva_ib_mul_rn(sqrt(2.0/3.0),b->z_max)),fabric);
        fabric=riva_add(fabric,riva_scale(term,riva_ib_mul_rn(
            riva_ib_mul_rn(-b->c_z,dl),riva_max(-d_re,0.0))));
    }
    double mb,md,xi; riva_surfaces(b,m,pressure,void_ratio,&mb,&md,&xi);
    (void)md; (void)xi;
    tensor_t alpha=riva_ib_div(s,pressure);
    double alpha_limit=riva_ib_mul_rn(sqrt(2.0/3.0),mb);
    if (old->base.geostatic_admitted)
        alpha_limit=riva_max(alpha_limit,riva_norm(old->base.alpha));
    const double alpha_norm=riva_norm(alpha);
    if (alpha_norm>alpha_limit) {
        alpha=riva_scale(alpha,alpha_limit/alpha_norm);
        s=riva_scale(alpha,pressure);
    }
    const tensor_t offset=riva_sub(alpha,alpha0);
    const double qa=riva_ddot(offset,offset);
    const double qb=riva_ib_mul_rn(2.0,riva_ddot(offset,alpha));
    const double qc=riva_ib_sub_rn(riva_ddot(alpha,alpha),
        riva_ib_mul_rn(2.0/3.0,riva_ib_mul_rn(mb,mb)));
    const double disc=riva_ib_sub_rn(riva_ib_mul_rn(qb,qb),
        riva_ib_mul_rn(riva_ib_mul_rn(4.0,qa),qc));
    int64_t fallbacks=old->base.beta_fallbacks; double beta;
    if (qa<=1.0e-14 || disc<0.0) { beta=b->beta0; fallbacks++; }
    else { beta=riva_ib_div_rn(riva_ib_add_rn(
               -qb,sqrt(riva_max(disc,0.0))),riva_ib_mul_rn(2.0,qa));
           if (!isfinite(beta)) { beta=b->beta0; fallbacks++; } }
    beta=riva_max(beta,1.0e-6);
    const tensor_t mapped=riva_add(alpha,riva_scale(offset,beta));
    const double mapped_norm=riva_norm(mapped);
    const tensor_t normal=mapped_norm>1.0e-12?
        riva_ib_div(mapped,mapped_norm):riva_zero();
    double D_ir,D_re; riva_ib_dilatancy(b,m,alpha,normal,beta,fabric,pressure,
        void_ratio,eps_ir,eps_re,&D_ir,&D_re);
    D_ir=riva_ib_mul_rn(D_ir,riva_ib_irreversible_factor(p,old));

    state.base.stress=riva_sub(s,riva_iso(pressure));
    state.base.alpha=alpha; state.base.alpha0=alpha0; state.base.alpha01=alpha01;
    state.base.n=normal; state.base.fabric=fabric;
    state.base.D_ir=D_ir; state.base.D_re=D_re;
    state.base.D=riva_ib_add_rn(D_ir,D_re);
    state.base.beta=beta; state.base.lambda_total=riva_ib_add_rn(
        old->base.lambda_total,dl);
    state.base.ep_eq_since_reversal=ep_eq; state.base.void_ratio=void_ratio;
    state.base.reversals=reversals;
    state.base.pressure_floor_hits=old->base.pressure_floor_hits+(at_floor?1:0);
    state.base.denominator_floor_hits=denom_hits;
    state.base.beta_fallbacks=fallbacks;
    state.base.eps_v_total=eps_total; state.base.eps_v_confining=eps_confining;
    state.base.eps_v_irreversible=eps_ir; state.base.eps_v_reversible=eps_re;
    if (reversals>old->base.reversals) {
        /* Capture the completed branch before its plastic counter is reset. */
        state.ep_half_last=old->base.ep_eq_since_reversal;
        const tensor_t reversal_dev=riva_dev(old->base.stress);
        const double excursion=riva_norm(riva_sub(reversal_dev,
            old->base.last_reversal_deviator));
        const double divisor=old->base.amplitude_reversals==0?1.0:2.0;
        const double amplitude=excursion/(divisor*riva_max(
            old->base.pressure_anchor,b->p_min));
        if (amplitude>1.0e-12) {
            state.base.cyclic_amplitude=amplitude;
            state.base.amplitude_factor=riva_amplitude_factor(b,amplitude,
                old->base.effective_knee_ratio);
        }
        state.base.last_reversal_deviator=reversal_dev;
        state.base.amplitude_reversals=old->base.amplitude_reversals+1;
    }
    return state;
}

RIVA_IB_HD static inline tensor_t riva_ib_unit(tensor_t value)
{ const double n=riva_norm(value); return n<=1.0e-14?riva_zero():riva_ib_div(value,n); }

RIVA_IB_HD static inline double riva_ib_mapping_capacity(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s,double mb)
{
    const double radius=sqrt(2.0/3.0)*mb;
    const double clearance=riva_max(p->mapping_center_limit_ratio*radius-
        riva_norm(s->mapping_anchor),0.0);
    return p->mapping_backstress_capacity_fraction*clearance;
}

RIVA_IB_HD static inline void riva_ib_mapping_intersection(
    const riva_ib_parameters_t *p,tensor_t alpha,tensor_t center,double mb,
    tensor_t fallback,double *beta,tensor_t *normal,tensor_t *ray_out,int *failed)
{
    const double radius=sqrt(2.0/3.0)*mb;
    const double core=p->mapping_core_radius_ratio*radius;
    const tensor_t raw=riva_sub(alpha,center);
    const tensor_t branch=riva_ib_unit(fallback);
    tensor_t ray=raw;
    if (riva_norm(branch)>1.0e-14) {
        const double along=riva_ddot(ray,branch);
        if (along<core) ray=riva_add(ray,riva_scale(branch,core-along));
    }
    if (riva_norm(ray)<=1.0e-14) {
        if (riva_norm(branch)<=1.0e-14) {
            *beta=p->base.beta0; *normal=riva_zero(); *ray_out=riva_zero();
            *failed=1; return;
        }
        ray=riva_scale(branch,core);
    }
    const double qa=riva_ddot(ray,ray);
    const double qb=2.0*riva_ddot(center,ray);
    const double qc=riva_ddot(center,center)-radius*radius;
    const double disc=qb*qb-4.0*qa*qc;
    if (qa<=1.0e-16 || disc<0.0) {
        *beta=p->base.beta0; *normal=riva_ib_unit(ray); *ray_out=*normal;
        *failed=1; return;
    }
    const double rt=sqrt(riva_max(disc,0.0));
    const double r1=(-qb+rt)/(2.0*qa),r2=(-qb-rt)/(2.0*qa);
    double scale=-1.0;
    if (r1>=0.0 && isfinite(r1)) scale=r1;
    if (r2>=0.0 && isfinite(r2)) scale=riva_max(scale,r2);
    if (scale<0.0) {
        *beta=p->base.beta0; *normal=riva_ib_unit(ray); *ray_out=*normal;
        *failed=1; return;
    }
    const tensor_t mapped=riva_add(center,riva_scale(ray,scale));
    *beta=riva_max(riva_norm(riva_sub(mapped,alpha))/
        riva_max(riva_norm(raw),core),1.0e-6);
    *normal=riva_ib_unit(mapped); *ray_out=riva_ib_unit(ray); *failed=0;
}

RIVA_IB_HD static inline tensor_t riva_ib_mapping_flow(
    const riva_ib_parameters_t *p,tensor_t normal,tensor_t ray,tensor_t fabric,
    const tensor_t *bias_direction)
{
    const tensor_t transverse=riva_sub(fabric,riva_scale(normal,
        riva_ddot(fabric,normal)));
    tensor_t flow=riva_add(riva_add(riva_scale(normal,1.0-p->mapping_ray_flow_weight),
        riva_scale(ray,p->mapping_ray_flow_weight)),
        riva_scale(transverse,p->mapping_fabric_flow_weight));
    if (bias_direction) flow=riva_add(flow,riva_scale(riva_ib_unit(*bias_direction),
        p->mapping_directional_ratchet_weight));
    if (riva_norm(flow)<=1.0e-14) return normal;
    flow=riva_ib_unit(flow);
    if (riva_ddot(flow,normal)<=0.05)
        flow=riva_ib_unit(riva_add(normal,riva_scale(transverse,0.05)));
    return flow;
}

RIVA_IB_HD static inline tensor_t riva_ib_backstress_update(
    const riva_ib_parameters_t *p,tensor_t old,tensor_t flow,double capacity,
    double dl,tensor_t *rate)
{
    const tensor_t target=riva_scale(flow,capacity);
    if (capacity<=1.0e-14 || dl<=0.0) {
        *rate=riva_scale(riva_sub(target,old),p->mapping_backstress_rate);
        return old;
    }
    const double decay=exp(-p->mapping_backstress_rate*dl);
    const tensor_t end=riva_add(target,riva_scale(riva_sub(old,target),decay));
    const tensor_t midpoint=riva_scale(riva_add(old,end),0.5);
    *rate=riva_scale(riva_sub(target,midpoint),p->mapping_backstress_rate);
    return end;
}

RIVA_IB_HD static inline tensor_t riva_ib_fabric_update(
    const riva_ib_parameters_t *p,tensor_t old,tensor_t normal,double dilatancy,
    double dl)
{
    const double dilation=riva_max(-dilatancy,0.0);
    const double contraction=riva_max(dilatancy,0.0);
    const double rate=p->mapping_fabric_dilation_rate*dilation+
        p->mapping_fabric_recovery_rate*contraction;
    if (rate<=1.0e-14 || dl<=0.0) return old;
    const tensor_t target=riva_scale(normal,-p->mapping_fabric_dilation_rate*
        dilation*p->mapping_fabric_saturation/rate);
    return riva_add(target,riva_scale(riva_sub(old,target),exp(-rate*dl)));
}

RIVA_IB_HD static inline double riva_ib_mapping_shear_ratio(
    const riva_ib_parameters_t *p,tensor_t backstress,double capacity)
{
    const double ratio=riva_norm(backstress)/riva_max(capacity,1.0e-14);
    const double activity=1.0-exp(-p->mapping_memory_shear_activation*ratio*ratio);
    return 1.0-(1.0-p->mapping_memory_shear_minimum_ratio)*activity;
}

RIVA_IB_HD static inline double riva_ib_directional_phase_scale(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s,tensor_t backstress,
    double capacity)
{
    const double high=riva_smoothstep((s->base.static_bias_index-
        p->mapping_high_bias_contraction_onset)/
        (p->mapping_high_bias_contraction_full-
         p->mapping_high_bias_contraction_onset));
    const tensor_t direction=riva_ib_unit(s->base.static_bias_tensor);
    const double memory=fabs(riva_ddot(backstress,direction))/
        riva_max(capacity,1.0e-14);
    return exp(-p->mapping_high_bias_contraction_gain*riva_ib_mapping_gate(p,s)*
        high*memory);
}

RIVA_IB_HD static inline riva_ib_state_t riva_ib_mapping_backbone(
    const riva_ib_parameters_t *p,const riva_material_parameters_t *m,
    const riva_ib_state_t *old,tensor_t deps,int force_reversal,
    int allow_legacy_reversal)
{
    const double gate=riva_ib_mapping_gate(p,old);
    if (gate<=1.0e-14 || !old->base.cyclic_phase_active)
        return riva_ib_base_backbone(p,m,old,deps,force_reversal,
                                     allow_legacy_reversal);
    const riva_parameters_t *b=&p->base;
    riva_ib_state_t state=*old;
    const double deps_v=riva_trace(deps);
    const tensor_t deps_d=riva_dev(deps);
    const tensor_t s_old=riva_dev(old->base.stress);
    const double p_old=riva_max(riva_pressure(old->base.stress),b->p_min);
    double shear_max,bulk; riva_moduli(b,m,p_old,&shear_max,&bulk);
    double mb_old,md_old,xi_old; riva_surfaces(b,m,p_old,old->base.void_ratio,
        &mb_old,&md_old,&xi_old); (void)md_old; (void)xi_old;
    const double capacity=riva_ib_mapping_capacity(p,old,mb_old);
    const double shear_ratio=riva_ib_mapping_shear_ratio(p,
        old->mapping_backstress,capacity);
    const double shear=shear_max*shear_ratio;
    const double threshold=riva_confining_strain(b,m,b->p_min,
        old->base.pressure_anchor);
    const int floor_active=b->compatibility_enabled &&
        old->base.eps_v_confining>=threshold;
    const double coupling_bulk=floor_active?0.0:bulk;
    tensor_t alpha0=old->base.alpha0,alpha01=old->base.alpha01;
    double ep_eq=old->base.ep_eq_since_reversal;
    int64_t reversals=old->base.reversals;
    double p_elastic;
    if (b->compatibility_enabled) { int hit=0;
        p_elastic=riva_pressure_from_confining(b,m,old->base.eps_v_confining+
            deps_v,old->base.pressure_anchor,&hit);
    } else p_elastic=riva_max(p_old-bulk*deps_v,b->p_min);
    const tensor_t alpha_trial=riva_ib_div(riva_add(s_old,
        riva_scale(deps_d,2.0*shear)),p_elastic);
    const double legacy=riva_ddot(riva_sub(alpha_trial,alpha0),
        riva_sub(alpha_trial,old->base.alpha));
    if (force_reversal || (allow_legacy_reversal && legacy<0.0)) {
        double beta_unused=old->base.beta;
        riva_ib_register_reversal(b,&old->base,&alpha0,&alpha01,&ep_eq,
            &beta_unused,&reversals);
    }
    const tensor_t loading=riva_ib_unit(deps_d);
    const tensor_t center_anchor=riva_add(riva_scale(alpha0,1.0-gate),
        riva_scale(old->mapping_anchor,gate));
    tensor_t backstress_end=old->mapping_backstress;
    tensor_t fabric_end=old->mapping_directional_fabric;
    double dl=0.0,d_ir_mid=riva_max(old->base.D_ir,0.0),d_re_mid=old->base.D_re;
    tensor_t flow=old->base.n;
    if (riva_norm(flow)<=1.0e-14) flow=loading;
    tensor_t normal=flow;
    double beta=riva_max(old->base.beta,1.0e-6);
    double kinematic=0.0;
    int64_t denom_hits=old->base.denominator_floor_hits;
    int64_t mapping_fallbacks=0,monotone_caps=old->mapping_monotone_caps;
    int at_floor=0;
    for (int32_t iteration=0;iteration<p->mapping_corrector_iterations;iteration++) {
        const tensor_t backstress_mid=riva_scale(riva_add(old->mapping_backstress,
            backstress_end),0.5);
        const tensor_t fabric_mid=riva_scale(riva_add(
            old->mapping_directional_fabric,fabric_end),0.5);
        const double total_d=d_ir_mid+d_re_mid;
        const double eps_mid=old->base.eps_v_confining+
            0.5*(deps_v+dl*total_d);
        double p_mid;
        if (b->compatibility_enabled) { int hit=0;
            p_mid=riva_pressure_from_confining(b,m,eps_mid,
                old->base.pressure_anchor,&hit); at_floor|=hit;
        } else { p_mid=riva_max(p_old-0.5*bulk*(deps_v+dl*total_d),b->p_min);
                 at_floor|=p_mid<=b->p_min; }
        const tensor_t s_mid=riva_add(s_old,riva_scale(riva_sub(deps_d,
            riva_scale(flow,dl)),shear));
        const tensor_t alpha_mid=riva_ib_div(s_mid,p_mid);
        const double void_mid=old->base.void_ratio+
            0.5*(1.0+old->base.void_ratio)*deps_v;
        double mb_mid,md_mid,xi_mid; riva_surfaces(b,m,p_mid,void_mid,
            &mb_mid,&md_mid,&xi_mid); (void)md_mid; (void)xi_mid;
        const tensor_t center_mid=riva_add(center_anchor,
            riva_scale(backstress_mid,gate));
        tensor_t ray; int failed=0;
        riva_ib_mapping_intersection(p,alpha_mid,center_mid,mb_mid,loading,
            &beta,&normal,&ray,&failed); mapping_fallbacks+=failed;
        flow=riva_ib_mapping_flow(p,normal,ray,fabric_mid,
                                  &old->base.static_bias_tensor);
        const tensor_t dil_fabric=riva_add(old->base.fabric,
            riva_scale(fabric_mid,gate*p->mapping_fabric_dilatancy_weight));
        riva_ib_dilatancy(b,m,alpha_mid,flow,beta,dil_fabric,p_mid,void_mid,
            old->base.eps_v_irreversible,old->base.eps_v_reversible,
            &d_ir_mid,&d_re_mid);
        d_ir_mid*=riva_ib_irreversible_factor(p,old);
        d_ir_mid*=riva_ib_directional_phase_scale(p,old,backstress_mid,capacity);
        if (dl>0.0 && d_re_mid>0.0)
            d_re_mid=riva_min(d_re_mid,old->base.eps_v_reversible/dl);
        tensor_t backstress_rate;
        const tensor_t candidate_backstress=riva_ib_backstress_update(p,
            old->mapping_backstress,flow,capacity,dl,&backstress_rate);
        (void)candidate_backstress;
        const double normal_flow=riva_max(riva_ddot(normal,flow),0.05);
        const double alpha_normal=riva_ddot(alpha_mid,normal);
        const double h=m->h*pow(riva_max(p_mid,b->p_min)/b->p_ref,-b->q_H);
        const double hardening=p_mid*h*pow(riva_max(beta,1.0e-12),m->m);
        kinematic=gate*p_mid*riva_ddot(normal,backstress_rate);
        double denominator=2.0*shear*normal_flow+(2.0/3.0)*hardening+
            kinematic-coupling_bulk*(d_ir_mid+d_re_mid)*alpha_normal;
        const double floor=b->denominator_floor_ratio*2.0*shear;
        if (denominator<floor) { denominator=floor; denom_hits++; }
        const double numerator=2.0*shear*riva_ddot(deps_d,normal)+
            coupling_bulk*deps_v*alpha_normal;
        double candidate=riva_max(0.0,numerator/denominator);
        candidate=(1.0-p->mapping_corrector_relaxation)*dl+
            p->mapping_corrector_relaxation*candidate;
        const double projection=riva_ddot(deps_d,flow);
        if (projection>1.0e-16) {
            const double limit=riva_ddot(deps_d,deps_d)/projection;
            if (candidate>limit) { candidate=limit; monotone_caps++; }
        }
        dl=candidate;
        backstress_end=riva_ib_backstress_update(p,old->mapping_backstress,
            flow,capacity,dl,&backstress_rate);
        fabric_end=riva_ib_fabric_update(p,old->mapping_directional_fabric,
            normal,d_ir_mid+d_re_mid,dl);
    }
    const tensor_t s_trial_end=riva_add(s_old,riva_scale(deps_d,2.0*shear));
    tensor_t s_end=riva_sub(s_trial_end,riva_scale(flow,2.0*shear*dl));
    double eps_ir=old->base.eps_v_irreversible-dl*riva_max(d_ir_mid,0.0);
    if (dl>0.0 && d_re_mid>0.0)
        d_re_mid=riva_min(d_re_mid,old->base.eps_v_reversible/dl);
    double eps_re=riva_max(old->base.eps_v_reversible-dl*d_re_mid,0.0);
    const double eps_total=old->base.eps_v_total+deps_v;
    double eps_confining=eps_total-eps_ir-eps_re;
    double pressure;
    if (b->compatibility_enabled) { int hit=0;
        pressure=riva_pressure_from_confining(b,m,eps_confining,
            old->base.pressure_anchor,&hit); at_floor|=hit;
    } else { pressure=riva_max(p_old-bulk*(deps_v+dl*(d_ir_mid+d_re_mid)),
                               b->p_min); at_floor|=pressure<=b->p_min; }
    const double void_ratio=old->base.void_ratio+
        (1.0+old->base.void_ratio)*deps_v;
    double mb,md,xi; riva_surfaces(b,m,pressure,void_ratio,&mb,&md,&xi);
    (void)md; (void)xi;
    double radius=pressure*sqrt(2.0/3.0)*mb;
    tensor_t alpha=riva_ib_div(s_end,pressure);
    double residual=riva_max(riva_norm(alpha)-sqrt(2.0/3.0)*mb,0.0);
    int corrections=0;
    if (residual>p->mapping_outer_tolerance) {
        const tensor_t plastic=riva_scale(flow,2.0*shear);
        const double qa=riva_ddot(plastic,plastic);
        const double qb=-2.0*riva_ddot(s_trial_end,plastic);
        const double qc=riva_ddot(s_trial_end,s_trial_end)-radius*radius;
        const double disc=qb*qb-4.0*qa*qc;
        double corrected=-1.0;
        if (qa>1.0e-16 && disc>=0.0) {
            const double rt=sqrt(riva_max(disc,0.0));
            const double r1=(-qb-rt)/(2.0*qa),r2=(-qb+rt)/(2.0*qa);
            if (r1>=0.0 && isfinite(r1)) corrected=r1;
            if (r2>=0.0 && isfinite(r2) &&
                (corrected<0.0 || fabs(r2-dl)<fabs(corrected-dl))) corrected=r2;
        }
        if (corrected>=0.0) {
            const double projection=riva_ddot(deps_d,flow);
            if (projection>1.0e-16) corrected=riva_min(corrected,
                riva_ddot(deps_d,deps_d)/projection);
            if (fabs(corrected-dl)>1.0e-14) {
                dl=corrected; corrections=1; tensor_t rate;
                backstress_end=riva_ib_backstress_update(p,old->mapping_backstress,
                    flow,capacity,dl,&rate);
                fabric_end=riva_ib_fabric_update(p,old->mapping_directional_fabric,
                    normal,d_ir_mid+d_re_mid,dl);
                s_end=riva_sub(s_trial_end,riva_scale(flow,2.0*shear*dl));
                eps_ir=old->base.eps_v_irreversible-dl*riva_max(d_ir_mid,0.0);
                if (dl>0.0 && d_re_mid>0.0)
                    d_re_mid=riva_min(d_re_mid,old->base.eps_v_reversible/dl);
                eps_re=riva_max(old->base.eps_v_reversible-dl*d_re_mid,0.0);
                eps_confining=eps_total-eps_ir-eps_re;
                if (b->compatibility_enabled) { int hit=0;
                    pressure=riva_pressure_from_confining(b,m,eps_confining,
                        old->base.pressure_anchor,&hit); at_floor|=hit;
                } else pressure=riva_max(p_old-bulk*(deps_v+dl*(d_ir_mid+d_re_mid)),
                                         b->p_min);
                alpha=riva_ib_div(s_end,pressure);
                riva_surfaces(b,m,pressure,void_ratio,&mb,&md,&xi);
                residual=riva_max(riva_norm(alpha)-sqrt(2.0/3.0)*mb,0.0);
            }
        }
    }
    const tensor_t center_end=riva_add(center_anchor,riva_scale(backstress_end,gate));
    tensor_t normal_end,ray_end; int failed=0; double beta_end;
    riva_ib_mapping_intersection(p,alpha,center_end,mb,flow,&beta_end,
        &normal_end,&ray_end,&failed); mapping_fallbacks+=failed;
    const tensor_t flow_end=riva_ib_mapping_flow(p,normal_end,ray_end,fabric_end,
        &old->base.static_bias_tensor);
    const tensor_t dil_fabric_end=riva_add(old->base.fabric,
        riva_scale(fabric_end,gate*p->mapping_fabric_dilatancy_weight));
    double d_ir_end,d_re_end; riva_ib_dilatancy(b,m,alpha,flow_end,beta_end,
        dil_fabric_end,pressure,void_ratio,eps_ir,eps_re,&d_ir_end,&d_re_end);
    d_ir_end*=riva_ib_irreversible_factor(p,old);
    d_ir_end*=riva_ib_directional_phase_scale(p,old,backstress_end,capacity);
    ep_eq=riva_ib_add_rn(ep_eq,dl);
    state.base.stress=riva_sub(s_end,riva_iso(pressure));
    state.base.alpha=alpha; state.base.alpha0=alpha0; state.base.alpha01=alpha01;
    state.base.n=flow_end; state.base.fabric=old->base.fabric;
    state.base.D_ir=d_ir_end; state.base.D_re=d_re_end; state.base.D=d_ir_end+d_re_end;
    state.base.beta=beta_end; state.base.lambda_total=riva_ib_add_rn(
        old->base.lambda_total,dl);
    state.base.ep_eq_since_reversal=ep_eq; state.base.void_ratio=void_ratio;
    state.base.reversals=reversals;
    state.base.pressure_floor_hits=old->base.pressure_floor_hits+(at_floor?1:0);
    state.base.denominator_floor_hits=denom_hits;
    state.base.beta_fallbacks=old->base.beta_fallbacks+mapping_fallbacks;
    state.base.eps_v_total=eps_total; state.base.eps_v_confining=eps_confining;
    state.base.eps_v_irreversible=eps_ir; state.base.eps_v_reversible=eps_re;
    state.mapping_backstress=backstress_end;
    state.mapping_directional_fabric=fabric_end;
    state.mapping_capacity=capacity; state.mapping_kinematic_denominator=kinematic;
    state.mapping_shear_modulus_ratio=shear_ratio;
    state.mapping_outer_residual=residual;
    state.mapping_stress_corrections=old->mapping_stress_corrections+corrections;
    state.mapping_corrector_passes=p->mapping_corrector_iterations;
    state.mapping_monotone_caps=monotone_caps;
    if (reversals>old->base.reversals) {
        /* The directional mapping backbone registers its own reversals. */
        state.ep_half_last=old->base.ep_eq_since_reversal;
        const tensor_t reversal_dev=riva_dev(old->base.stress);
        const double excursion=riva_norm(riva_sub(reversal_dev,
            old->base.last_reversal_deviator));
        const double divisor=old->base.amplitude_reversals==0?1.0:2.0;
        const double amplitude=excursion/(divisor*riva_max(
            old->base.pressure_anchor,b->p_min));
        if (amplitude>1.0e-12) {
            state.base.cyclic_amplitude=amplitude;
            state.base.amplitude_factor=riva_amplitude_factor(b,amplitude,
                old->base.effective_knee_ratio);
        }
        state.base.last_reversal_deviator=reversal_dev;
        state.base.amplitude_reversals=old->base.amplitude_reversals+1;
    }
    return state;
}

RIVA_IB_HD static inline double riva_ib_ratchet_increment(
    const riva_ib_parameters_t *p,const riva_ib_state_t *old,
    const riva_ib_state_t *trial,tensor_t *direction)
{
    const riva_parameters_t *b=&p->base;
    *direction=riva_zero();
    const double loose_gate=riva_ib_loose_gate(p,old);
    double capacity,rate;
    if (p->loose_shear_flow_enabled && loose_gate>1.0e-14) {
        const double parent_activity=riva_bias_ratchet_activity(b,&old->base);
        const double activity=parent_activity+loose_gate*
            (riva_ib_phase_amplitude_activity(p,old)-parent_activity);
        const double bias=riva_ib_projected_bias(old);
        const double parent_capacity=bias<=1.0e-14?0.0:
            b->bias_ratchet_limit*activity*pow(b->bias_ratchet_reference_bias/bias,
                b->bias_ratchet_bias_exponent);
        const double loose_capacity=p->loose_shear_ratchet_capacity*activity;
        capacity=parent_capacity+loose_gate*(loose_capacity-parent_capacity);
        if (capacity<=1.0e-14 || old->base.bias_ratchet_strain>=capacity) return 0.0;
        const double pressure=riva_pressure(old->base.stress);
        const double pr=pow(riva_clip(pressure/riva_max(old->base.pressure_anchor,
            b->p_min),0.0,1.0),p->loose_shear_ratchet_pressure_exponent);
        const double saturation=riva_max(1.0-old->base.bias_ratchet_strain/capacity,0.0);
        const double loose_rate=p->loose_shear_ratchet_rate*activity*pr*saturation;
        const double stress_level=pow(riva_max(old->base.pressure_anchor/
            b->bias_reference_pressure-1.0,0.0),b->bias_ratchet_pressure_exponent);
        const double bias_factor=pow(b->bias_ratchet_reference_bias/
            riva_max(bias,1.0e-14),b->bias_ratchet_bias_exponent);
        const double parent_rate=b->bias_ratchet_rate*stress_level*bias_factor*saturation;
        rate=parent_rate+loose_gate*(loose_rate-parent_rate);
    } else {
        const double bias=riva_ib_projected_bias(old);
        if (bias<=1.0e-14) return 0.0;
        capacity=b->bias_ratchet_limit*riva_bias_ratchet_activity(b,&old->base)*
            pow(b->bias_ratchet_reference_bias/bias,b->bias_ratchet_bias_exponent);
        if (capacity<=1.0e-14 || old->base.bias_ratchet_strain>=capacity) return 0.0;
        const double stress_level=pow(riva_max(old->base.pressure_anchor/
            b->bias_reference_pressure-1.0,0.0),b->bias_ratchet_pressure_exponent);
        const double bias_factor=pow(b->bias_ratchet_reference_bias/
            riva_max(bias,1.0e-14),b->bias_ratchet_bias_exponent);
        const double saturation=riva_max(1.0-old->base.bias_ratchet_strain/capacity,0.0);
        rate=b->bias_ratchet_rate*stress_level*bias_factor*saturation;
    }
    const double dl=riva_max(trial->base.lambda_total-old->base.lambda_total,0.0);
    const double norm=riva_norm(old->base.cyclic_direction);
    const double projection=riva_ddot(old->base.static_bias_tensor,
        old->base.cyclic_direction);
    if (dl<=0.0 || norm<=1.0e-12 || rate<=0.0) return 0.0;
    const double sign=projection>0.0?1.0:(projection<0.0?-1.0:0.0);
    *direction=riva_scale(old->base.cyclic_direction,sign/norm);
    const double increment=riva_min(rate*dl,
        riva_max(capacity-old->base.bias_ratchet_strain,0.0));
    return (1.0-riva_ib_mapping_gate(p,old))*increment;
}

RIVA_IB_HD static inline double riva_ib_bias_volume_target(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s,int inside_substeps)
{
    double target=riva_bias_reversible_volume_target(&p->base,&s->base);
    const double plastic_gate=riva_smoothstep(
        s->ep_half_last/p->bias_reversible_volume_ep_ref);
    target=riva_ib_mul_rn(target,plastic_gate);
    if (!p->phase_transformation_enabled) return target;
    target*=1.0-riva_ib_phase_volume_replacement_gate(p,s);
    if (inside_substeps) target+=s->phase_reversible_volume;
    return target;
}

RIVA_IB_HD static inline riva_ib_state_t riva_ib_mechanical_state(
    const riva_ib_parameters_t *p,const riva_material_parameters_t *m,
    const riva_ib_state_t *s)
{
    if (fabs(s->base.bias_reversible_volume)<=1.0e-16) return *s;
    riva_ib_state_t out=*s;
    out.base.eps_v_total=s->base.physical_eps_v_total;
    out.base.eps_v_confining=out.base.eps_v_total-out.base.eps_v_irreversible-
        out.base.eps_v_reversible;
    tensor_t dev=riva_dev(s->base.stress); int hit=0;
    const double pressure=riva_pressure_from_confining(&p->base,m,
        out.base.eps_v_confining,out.base.pressure_anchor,&hit); (void)hit;
    if (s->base.geostatic_admitted)
        dev=riva_regularize_deviator(&p->base,m,pressure,out.base.void_ratio,
            s->base.alpha,1,dev,&out.base.alpha);
    else out.base.alpha=riva_ib_div(dev,pressure);
    out.base.stress=riva_sub(dev,riva_iso(pressure));
    return out;
}

RIVA_IB_HD static inline riva_ib_state_t riva_ib_backbone(
    const riva_ib_parameters_t *p,const riva_material_parameters_t *m,
    const riva_ib_state_t *old,tensor_t deps,int force_reversal,
    int allow_legacy_reversal)
{
    return riva_ib_mapping_gate(p,old)>1.0e-14 && old->base.cyclic_phase_active?
        riva_ib_mapping_backbone(p,m,old,deps,force_reversal,allow_legacy_reversal):
        riva_ib_base_backbone(p,m,old,deps,force_reversal,allow_legacy_reversal);
}

RIVA_IB_HD static inline riva_ib_state_t riva_ib_forward_euler(
    const riva_ib_parameters_t *p,const riva_material_parameters_t *m,
    const riva_ib_state_t *old,tensor_t deps,int force_reversal,
    int allow_legacy_reversal,int inside_substeps)
{
    const double physical_deps_v=riva_trace(deps);
    const riva_ib_state_t mechanical=riva_ib_mechanical_state(p,m,old);
    const riva_ib_state_t trial=riva_ib_backbone(p,m,&mechanical,deps,
        force_reversal,allow_legacy_reversal);
    tensor_t ratchet_direction;
    const double ratchet=riva_ib_ratchet_increment(p,old,&trial,&ratchet_direction);
    riva_ib_state_t provisional=trial;
    if (ratchet>0.0) {
        const tensor_t effective=riva_sub(deps,riva_scale(ratchet_direction,ratchet));
        provisional=riva_ib_backbone(p,m,&mechanical,effective,force_reversal,
                                     allow_legacy_reversal);
    }
    if (provisional.base.amplitude_reversals>old->base.amplitude_reversals) {
        const tensor_t excursion=riva_sub(riva_dev(old->base.stress),
            old->base.last_reversal_deviator);
        const double norm=riva_norm(excursion);
        if (norm>1.0e-12)
            provisional.base.cyclic_direction=riva_ib_div(excursion,norm);
    }
    const double target=riva_ib_bias_volume_target(p,&provisional,inside_substeps);
    riva_ib_state_t state=provisional;
    state.base.bias_ratchet_strain=old->base.bias_ratchet_strain+ratchet;
    state.base.physical_eps_v_total=old->base.physical_eps_v_total+physical_deps_v;
    state.base.bias_reversible_volume=target;
    const double effective_total=state.base.physical_eps_v_total+target;
    const int changed=fabs(effective_total-state.base.eps_v_total)>1.0e-16;
    state.base.eps_v_total=effective_total;
    state.base.eps_v_confining=state.base.eps_v_total-
        state.base.eps_v_irreversible-state.base.eps_v_reversible;
    if (changed) {
        tensor_t dev=riva_dev(state.base.stress); int hit=0;
        const double pressure=riva_pressure_from_confining(&p->base,m,
            state.base.eps_v_confining,state.base.pressure_anchor,&hit);
        if (state.base.geostatic_admitted)
            dev=riva_regularize_deviator(&p->base,m,pressure,state.base.void_ratio,
                state.base.alpha,1,dev,&state.base.alpha);
        else state.base.alpha=riva_ib_div(dev,pressure);
        state.base.stress=riva_sub(dev,riva_iso(pressure));
        state.base.pressure_floor_hits+=hit?1:0;
    }
    state.base.void_ratio=old->base.void_ratio+
        (1.0+old->base.void_ratio)*physical_deps_v;
    state.base.D_re*=riva_bias_reversible_factor(&p->base,&state.base);
    state.base.D=state.base.D_ir+state.base.D_re;
    return state;
}

RIVA_IB_HD static inline void riva_ib_rebuild_pressure(
    const riva_ib_parameters_t *p,const riva_material_parameters_t *m,
    riva_ib_state_t *s)
{
    s->base.eps_v_total=riva_ib_add_rn(s->base.physical_eps_v_total,
        s->base.bias_reversible_volume);
    s->base.eps_v_confining=riva_ib_sub_rn(riva_ib_sub_rn(
        s->base.eps_v_total,s->base.eps_v_irreversible),
        s->base.eps_v_reversible);
    tensor_t dev=riva_dev(s->base.stress); int hit=0;
    const double pressure=riva_pressure_from_confining(&p->base,m,
        s->base.eps_v_confining,s->base.pressure_anchor,&hit);
    if (s->base.geostatic_admitted)
        dev=riva_regularize_deviator(&p->base,m,pressure,s->base.void_ratio,
            s->base.alpha,1,dev,&s->base.alpha);
    else s->base.alpha=riva_ib_div(dev,pressure);
    s->base.stress=riva_sub(dev,riva_iso(pressure));
    s->base.pressure_floor_hits+=hit?1:0;
}

RIVA_IB_HD static inline double riva_ib_phase_dense_transition(
    const riva_ib_parameters_t *p,const riva_ib_state_t *s)
{ return riva_smoothstep((riva_ib_initial_density(s)-
    p->reference_relative_density_value)/riva_max(p->phase_volume_density_full-
    p->reference_relative_density_value,1.0e-14)); }

RIVA_IB_HD static inline riva_ib_state_t riva_ib_apply_host_phase_volume(
    const riva_ib_parameters_t *p,const riva_material_parameters_t *m,
    const riva_ib_state_t *old,const riva_ib_state_t *input,tensor_t deps)
{
    riva_ib_state_t state=*input;
    const double path=riva_norm(riva_dev(deps));
    const double old_activity=riva_ib_phase_activity(p,old,p->phase_bias_exponent,
        p->phase_pressure_exponent);
    const double new_activity=riva_ib_phase_activity(p,&state,p->phase_bias_exponent,
        p->phase_pressure_exponent);
    const double activity=riva_ib_mul_rn(
        0.5,riva_ib_add_rn(old_activity,new_activity));
    const double old_wave=riva_ib_phase_activity(p,old,p->phase_wave_bias_exponent,
        p->phase_wave_pressure_exponent);
    const double new_wave=riva_ib_phase_activity(p,&state,p->phase_wave_bias_exponent,
        p->phase_wave_pressure_exponent);
    const double old_mean=riva_ib_phase_activity(p,old,p->phase_mean_bias_exponent,
        p->phase_mean_pressure_exponent);
    const double new_mean=riva_ib_phase_activity(p,&state,p->phase_mean_bias_exponent,
        p->phase_mean_pressure_exponent);
    double loading,eta_pt,signed_phase;
    riva_ib_phase_coordinates(p,m,&state,&loading,&eta_pt,&signed_phase);
    (void)loading; (void)eta_pt;
    const double contractile=riva_ib_mul_rn(
        0.5,riva_ib_sub_rn(1.0,signed_phase));
    double ir_increment=riva_ib_mul_rn(-p->phase_contraction_rate,activity);
    ir_increment=riva_ib_mul_rn(
        ir_increment,riva_ib_transformation_zone(p,m,&state));
    ir_increment=riva_ib_mul_rn(ir_increment,contractile);
    ir_increment=riva_ib_mul_rn(ir_increment,path);
    state.phase_irreversible_volume=riva_ib_add_rn(
        old->phase_irreversible_volume,ir_increment);
    state.base.eps_v_irreversible=riva_ib_add_rn(
        state.base.eps_v_irreversible,ir_increment);
    const double old_potential=riva_ib_signed_phase_potential(p,m,old);
    const double new_potential=riva_ib_signed_phase_potential(p,m,&state);
    const double old_memory=riva_ib_sub_rn(1.0,exp(riva_ib_div_rn(
        -riva_max(-old->phase_irreversible_volume,0.0),
        p->phase_memory_reference_volume)));
    const double new_memory=riva_ib_sub_rn(1.0,exp(riva_ib_div_rn(
        -riva_max(-state.phase_irreversible_volume,0.0),
        p->phase_memory_reference_volume)));
    const double transition=riva_ib_mul_rn(0.5,riva_ib_add_rn(
        riva_ib_phase_dense_transition(p,old),
        riva_ib_phase_dense_transition(p,&state)));
    const double wave_mult=riva_ib_add_rn(
        p->phase_intermediate_wave_multiplier,riva_ib_mul_rn(transition,
        riva_ib_sub_rn(1.0,p->phase_intermediate_wave_multiplier)));
    const double mean_mult=riva_ib_add_rn(
        p->phase_intermediate_mean_multiplier,riva_ib_mul_rn(transition,
        riva_ib_sub_rn(1.0,p->phase_intermediate_mean_multiplier)));
    double old_wave_target=riva_ib_mul_rn(p->phase_reversible_scale,wave_mult);
    old_wave_target=riva_ib_mul_rn(old_wave_target,old_wave);
    old_wave_target=riva_ib_mul_rn(old_wave_target,riva_ib_sub_rn(
        old_potential,old->phase_potential_anchor));
    double old_mean_target=riva_ib_mul_rn(
        p->phase_reversible_mean_scale,mean_mult);
    old_mean_target=riva_ib_mul_rn(old_mean_target,old_mean);
    old_mean_target=riva_ib_mul_rn(old_mean_target,old_memory);
    const double old_target=riva_ib_add_rn(old_wave_target,old_mean_target);
    double new_wave_target=riva_ib_mul_rn(p->phase_reversible_scale,wave_mult);
    new_wave_target=riva_ib_mul_rn(new_wave_target,new_wave);
    new_wave_target=riva_ib_mul_rn(new_wave_target,riva_ib_sub_rn(
        new_potential,state.phase_potential_anchor));
    double new_mean_target=riva_ib_mul_rn(
        p->phase_reversible_mean_scale,mean_mult);
    new_mean_target=riva_ib_mul_rn(new_mean_target,new_mean);
    new_mean_target=riva_ib_mul_rn(new_mean_target,new_memory);
    const double new_target=riva_ib_add_rn(new_wave_target,new_mean_target);
    const double target=riva_ib_mul_rn(
        0.5,riva_ib_add_rn(old_target,new_target));
    const double bias_ratio=riva_max(riva_ib_div_rn(
        riva_ib_projected_bias(&state),p->branch_compliance_bias_reference),1.0);
    const double relaxation_blend=riva_ib_add_rn(
        p->phase_intermediate_relaxation_ratio,riva_ib_mul_rn(transition,
        riva_ib_sub_rn(1.0,p->phase_intermediate_relaxation_ratio)));
    double relaxation_strain=riva_ib_mul_rn(
        p->phase_reversible_relaxation_strain,relaxation_blend);
    relaxation_strain=riva_ib_mul_rn(relaxation_strain,
        pow(bias_ratio,p->phase_reversible_relaxation_bias_exponent));
    const double high=riva_ib_intermediate_high_gate(p,&state);
    const double multiplier=riva_max(riva_ib_add_rn(1.0,riva_ib_mul_rn(high,
        riva_ib_sub_rn(
            p->intermediate_high_bias_phase_relaxation_multiplier,1.0))),0.0);
    const double exponent=riva_ib_div_rn(
        riva_ib_mul_rn(-multiplier,path),relaxation_strain);
    const double relaxation=riva_ib_sub_rn(1.0,exp(exponent));
    state.phase_reversible_volume=riva_ib_add_rn(
        old->phase_reversible_volume,riva_ib_mul_rn(relaxation,
        riva_ib_sub_rn(target,old->phase_reversible_volume)));
    state.base.bias_reversible_volume=riva_ib_add_rn(riva_ib_sub_rn(
        state.base.bias_reversible_volume,old->phase_reversible_volume),
        state.phase_reversible_volume);
    riva_ib_rebuild_pressure(p,m,&state);
    if (riva_ib_mapping_gate(p,old)>1.0e-14) {
        const double scale=0.5*(riva_ib_directional_phase_scale(p,old,
            old->mapping_backstress,old->mapping_capacity)+
            riva_ib_directional_phase_scale(p,&state,state.mapping_backstress,
                state.mapping_capacity));
        const double increment=state.phase_irreversible_volume-
            old->phase_irreversible_volume;
        const double correction=(scale-1.0)*increment;
        state.phase_irreversible_volume=old->phase_irreversible_volume+scale*increment;
        state.base.eps_v_irreversible+=correction;
        state.mapping_phase_contraction_scale=scale;
        riva_ib_rebuild_pressure(p,m,&state);
    }
    return state;
}

RIVA_IB_HD static inline riva_ib_state_t riva_ib_host_outer_correction(
    const riva_ib_parameters_t *p,const riva_material_parameters_t *m,
    const riva_ib_state_t *input)
{
    if (riva_ib_mapping_gate(p,input)<=1.0e-14) return *input;
    riva_ib_state_t s=*input; int correction_count=0;
    for (int pass=0;pass<2;pass++) {
        tensor_t dev=riva_dev(s.base.stress);
        const double pressure=riva_max(riva_pressure(s.base.stress),p->base.p_min);
        double mb,md,xi; riva_surfaces(&p->base,m,pressure,s.base.void_ratio,
            &mb,&md,&xi); (void)md; (void)xi;
        const double radius=pressure*sqrt(2.0/3.0)*mb;
        const double norm=riva_norm(dev);
        const double residual=riva_max(norm/pressure-sqrt(2.0/3.0)*mb,0.0);
        s.mapping_outer_residual=residual;
        if (residual<=p->mapping_outer_tolerance) break;
        const tensor_t direction=riva_ib_unit(dev);
        double shear,bulk; riva_moduli(&p->base,m,pressure,&shear,&bulk); (void)bulk;
        const double dl=riva_max((norm-radius)/riva_max(2.0*shear,1.0e-14),0.0);
        if (dl<=0.0) break;
        tensor_t rate;
        s.mapping_backstress=riva_ib_backstress_update(p,s.mapping_backstress,
            direction,s.mapping_capacity,dl,&rate);
        s.mapping_directional_fabric=riva_ib_fabric_update(p,
            s.mapping_directional_fabric,direction,s.base.D,dl);
        dev=riva_scale(direction,radius);
        s.base.stress=riva_sub(dev,riva_iso(pressure));
        s.base.alpha=riva_ib_div(dev,pressure);
        s.base.lambda_total+=dl; s.base.ep_eq_since_reversal+=dl;
        const double gate=riva_ib_mapping_gate(p,&s);
        const tensor_t center=riva_add(riva_scale(s.base.alpha0,1.0-gate),
            riva_scale(riva_add(s.mapping_anchor,s.mapping_backstress),gate));
        tensor_t normal,ray; double beta; int failed=0;
        riva_ib_mapping_intersection(p,s.base.alpha,center,mb,direction,
            &beta,&normal,&ray,&failed);
        s.base.beta=beta; s.base.n=riva_ib_mapping_flow(p,normal,ray,
            s.mapping_directional_fabric,&s.base.static_bias_tensor);
        s.base.beta_fallbacks+=failed; correction_count++;
    }
    s.mapping_stress_corrections+=correction_count;
    const double pressure=riva_max(riva_pressure(s.base.stress),p->base.p_min);
    double mb,md,xi; riva_surfaces(&p->base,m,pressure,s.base.void_ratio,
        &mb,&md,&xi); (void)md; (void)xi;
    s.mapping_outer_residual=riva_max(riva_norm(riva_dev(s.base.stress))/pressure-
        sqrt(2.0/3.0)*mb,0.0);
    return s;
}

RIVA_IB_HD static inline int riva_ib_update_material(
    const riva_ib_parameters_t *p,const riva_material_parameters_t *m,
    tensor_t deps,int32_t nsub,riva_ib_state_t *state,tensor_t *stress_new,
    riva_update_info_t *info)
{
    if (!p || !m || !state || !state->base.initialized || nsub<1 ||
        !riva_material_parameters_valid(&p->base,m) || !riva_finite_tensor(deps))
        return 0;
    const int objective=p->base.objective_reversal_enabled;
    int valid=0;
    const tensor_t direction=objective?riva_ib_host_direction(&p->base,deps,&valid):
        riva_zero();
    const int reversal=objective?riva_ib_host_reversal(&p->base,&state->base,
        direction,valid):0;
    const riva_ib_state_t initial=*state;
    riva_ib_state_t current=initial;
    const int phase_active=p->phase_transformation_enabled &&
        (riva_ib_phase_dense_bias_gate(p,&initial)>1.0e-14 ||
         riva_ib_phase_volume_gate(p,&initial)>1.0e-14);
    const tensor_t sub=riva_ib_div(deps,(double)nsub);
    for (int32_t i=0;i<nsub;i++)
        current=riva_ib_forward_euler(p,m,&current,sub,
            objective && reversal && i==0,!objective,phase_active);
    if (objective && valid)
        current.base.last_host_deviatoric_strain_direction=direction;
    if (phase_active)
        current=riva_ib_apply_host_phase_volume(p,m,&initial,&current,deps);
    current.phase_accumulation_hardening_state=
        riva_ib_phase_accumulation_target_hardening(p,&current);
    if (riva_ib_unbiased_gate(p,&current)>1.0e-14 &&
        riva_norm(current.unbiased_phase_direction)<=1.0e-14 &&
        riva_norm(current.base.cyclic_direction)>1.0e-14)
        current.unbiased_phase_direction=riva_ib_unit(current.base.cyclic_direction);
    current.loose_shear_hardening_state=riva_ib_loose_target_hardening(p,&current);
    if (riva_ib_mapping_gate(p,&current)>1.0e-14) {
        current=riva_ib_host_outer_correction(p,m,&current);
        current.loose_shear_hardening_state=riva_ib_loose_target_hardening(p,&current);
    }
    if (!riva_finite_tensor(current.base.stress) ||
        !isfinite(current.base.lambda_total) || !isfinite(current.base.void_ratio))
        return 0;
    *state=current; if (stress_new) *stress_new=current.base.stress;
    if (info) { info->accepted_substeps=nsub; info->reversal_registered=reversal; }
    return 1;
}

RIVA_IB_HD static inline int riva_ib_update(
    const riva_ib_parameters_t *p,tensor_t deps,int32_t nsub,
    riva_ib_state_t *state,tensor_t *stress_new,riva_update_info_t *info)
{
    const riva_material_parameters_t m=riva_reference_material_parameters(&p->base);
    return riva_ib_update_material(p,&m,deps,nsub,state,stress_new,info);
}

RIVA_IB_HD static inline int riva_ib_state_values(
    const riva_ib_state_t *s,double values[RIVA_IB_STATE_VALUE_COUNT])
{
    int i=riva_state_values(&s->base,values);
#define RIVA_IB_STATE_TENSOR(t) do { \
    values[i++]=(t).xx; values[i++]=(t).yy; values[i++]=(t).zz; \
    values[i++]=(t).xy; values[i++]=(t).yz; values[i++]=(t).xz; \
} while (0)
    values[i++]=s->phase_irreversible_volume;
    values[i++]=s->phase_reversible_volume;
    values[i++]=s->phase_potential_anchor;
    values[i++]=s->phase_accumulation_lambda_anchor;
    values[i++]=s->phase_accumulation_hardening_state;
    RIVA_IB_STATE_TENSOR(s->unbiased_phase_direction);
    values[i++]=s->loose_shear_lambda_anchor;
    values[i++]=s->loose_shear_hardening_state;
    values[i++]=s->loose_shear_gate_value;
    RIVA_IB_STATE_TENSOR(s->mapping_anchor);
    RIVA_IB_STATE_TENSOR(s->mapping_backstress);
    RIVA_IB_STATE_TENSOR(s->mapping_directional_fabric);
    values[i++]=s->mapping_gate_value;
    values[i++]=s->mapping_capacity;
    values[i++]=s->mapping_kinematic_denominator;
    values[i++]=s->mapping_shear_modulus_ratio;
    values[i++]=s->mapping_phase_contraction_scale;
    values[i++]=s->mapping_outer_residual;
    values[i++]=(double)s->mapping_stress_corrections;
    values[i++]=(double)s->mapping_corrector_passes;
    values[i++]=(double)s->mapping_monotone_caps;
    values[i++]=s->initial_relative_density_value;
    values[i++]=s->intermediate_low_gate_value;
    values[i++]=s->intermediate_high_gate_base;
    values[i++]=s->ep_half_last;
#undef RIVA_IB_STATE_TENSOR
    return i;
}

#undef riva_dev
#undef riva_trace
#undef riva_scale
#undef riva_sub
#undef riva_add
#undef riva_norm
#undef riva_ddot
#undef RIVA_IB_HD
#endif /* HERCULES_NONLINEAR_RIVA_SAND_INTERMEDIATE_BIAS_RESEARCH_H */
