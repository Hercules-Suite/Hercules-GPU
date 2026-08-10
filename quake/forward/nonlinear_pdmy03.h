/* Shared host/CUDA PressureDependMultiYield03 constitutive integration. */
#ifndef HERCULES_NONLINEAR_PDMY03_H
#define HERCULES_NONLINEAR_PDMY03_H

#include <math.h>
#include <stdint.h>

#include "nonlinear.h"

#ifdef __CUDACC__
#define PDMY03_HD __host__ __device__
#else
#define PDMY03_HD
#endif

#define PDMY03_SMALL 1.0e-12
#define PDMY03_SURFACE_LIMIT 1.0e30

typedef enum pdmy03_failure_t {
    PDMY03_SUCCESS = 0,
    PDMY03_FAILURE_CONSTANTS = 1,
    PDMY03_FAILURE_INITIAL_PRESSURE = 2,
    PDMY03_FAILURE_NONFINITE = 3,
    PDMY03_FAILURE_SURFACE_GEOMETRY = 4
} pdmy03_failure_t;

typedef struct pdmy03_diagnostics_t {
    int32_t reason;
    int32_t requested_substeps;
    int32_t accepted_substeps;
    int32_t crossed_surfaces;
    int32_t substeps_clamped;
    double yield_residual;
    double mean_stress;
    double strain_increment_norm;
} pdmy03_diagnostics_t;

PDMY03_HD static inline double pdmy03_max(double a, double b)
{ return a > b ? a : b; }

PDMY03_HD static inline double pdmy03_min(double a, double b)
{ return a < b ? a : b; }

PDMY03_HD static inline tensor_t pdmy03_zero(void)
{
    tensor_t a = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    return a;
}

PDMY03_HD static inline tensor_t pdmy03_add(tensor_t a, tensor_t b)
{
    tensor_t r = {a.xx+b.xx, a.yy+b.yy, a.zz+b.zz,
                  a.xy+b.xy, a.yz+b.yz, a.xz+b.xz};
    return r;
}

PDMY03_HD static inline tensor_t pdmy03_sub(tensor_t a, tensor_t b)
{
    tensor_t r = {a.xx-b.xx, a.yy-b.yy, a.zz-b.zz,
                  a.xy-b.xy, a.yz-b.yz, a.xz-b.xz};
    return r;
}

PDMY03_HD static inline tensor_t pdmy03_scale(tensor_t a, double s)
{
    tensor_t r = {a.xx*s, a.yy*s, a.zz*s, a.xy*s, a.yz*s, a.xz*s};
    return r;
}

PDMY03_HD static inline double pdmy03_mean(tensor_t a)
{ return (a.xx+a.yy+a.zz)/3.0; }

PDMY03_HD static inline tensor_t pdmy03_deviator(tensor_t a)
{
    const double p = pdmy03_mean(a);
    a.xx -= p; a.yy -= p; a.zz -= p;
    return a;
}

PDMY03_HD static inline tensor_t pdmy03_from_dev_mean(tensor_t d, double p)
{
    d = pdmy03_deviator(d);
    d.xx += p; d.yy += p; d.zz += p;
    return d;
}

PDMY03_HD static inline double pdmy03_ddot(tensor_t a, tensor_t b)
{
    return a.xx*b.xx+a.yy*b.yy+a.zz*b.zz+
           2.0*(a.xy*b.xy+a.yz*b.yz+a.xz*b.xz);
}

PDMY03_HD static inline double pdmy03_norm(tensor_t a)
{ return sqrt(pdmy03_max(pdmy03_ddot(a, a), 0.0)); }

PDMY03_HD static inline int pdmy03_finite_tensor(tensor_t a)
{
    return isfinite(a.xx) && isfinite(a.yy) && isfinite(a.zz) &&
           isfinite(a.xy) && isfinite(a.yz) && isfinite(a.xz);
}

PDMY03_HD static inline double pdmy03_octahedral_shear(tensor_t a,
                                                       int engineering)
{
    const double value = pdmy03_norm(pdmy03_deviator(a))/sqrt(3.0);
    return engineering ? 2.0*value : value;
}

PDMY03_HD static inline double pdmy03_surface_size(
    const pdmy03_constants_t *c, int surface)
{
    return c->failure_ratio*(double)surface/(double)c->surface_count;
}

PDMY03_HD static inline double pdmy03_yield(
    const pdmy03_constants_t *c, tensor_t stress, const tensor_t *centers,
    int surface)
{
    const double cone = pdmy03_mean(stress)-c->residual_pressure;
    const tensor_t shifted = pdmy03_sub(
        pdmy03_deviator(stress), pdmy03_scale(centers[surface-1], cone));
    const double radius = pdmy03_surface_size(c, surface)*cone;
    return 1.5*pdmy03_ddot(shifted, shifted)-radius*radius;
}

/* Prepare the standard automatically generated nested surfaces. Input angles
 * are in degrees and reference pressure is a positive compression magnitude. */
PDMY03_HD static inline int pdmy03_prepare_constants(pdmy03_constants_t *c)
{
    const double pi = 3.1415926535897932384626433832795;
    if (c == NULL || c->surface_count < 2 ||
        c->surface_count > PDMY03_MAX_SURFACES ||
        !(c->ref_shear_modulus > 0.0) || !(c->ref_bulk_modulus > 0.0) ||
        !(c->peak_shear_strain > 0.0) || !(c->ref_pressure > 0.0) ||
        !(c->atmospheric_pressure > 0.0) || c->pressure_dependence < 0.0 ||
        c->loading_type < 0 || c->loading_type > 2) return 0;

    const double sin_phi = sin(c->friction_angle*pi/180.0);
    const double sin_pt = sin(c->phase_transform_angle*pi/180.0);
    if (!(sin_phi > 0.0) || sin_phi >= 1.0 || sin_pt < 0.0 || sin_pt >= 1.0)
        return 0;
    if (c->loading_type == 1) {
        c->failure_ratio = 6.0*sin_phi/(3.0+sin_phi);
        c->phase_transform_ratio = 6.0*sin_pt/(3.0+sin_pt);
    } else if (c->loading_type == 2) {
        c->failure_ratio = sqrt(3.0)*sin_phi;
        c->phase_transform_ratio = sqrt(3.0)*sin_pt;
    } else {
        c->failure_ratio = 6.0*sin_phi/(3.0-sin_phi);
        c->phase_transform_ratio = 6.0*sin_pt/(3.0-sin_pt);
    }
    if (!(c->failure_ratio > 0.0)) return 0;

    c->residual_pressure = 2.0*c->cohesion/c->failure_ratio;
    c->residual_pressure = pdmy03_max(
        c->residual_pressure, 1.0e-4*c->atmospheric_pressure);
    const double cone_height = c->ref_pressure+c->residual_pressure;
    const double peak_shear = sqrt(2.0)*cone_height*c->failure_ratio/3.0;
    const double denominator =
        c->ref_shear_modulus*c->peak_shear_strain-peak_shear;
    if (!(denominator > 0.0)) return 0;
    const double ref_strain = c->peak_shear_strain*peak_shear/denominator;
    c->reference_strain = ref_strain;
    c->reference_peak_shear = peak_shear;
    const double stress_increment = peak_shear/(double)c->surface_count;
    c->phase_transform_strain = 0.0;

    for (int i = 1; i <= c->surface_count; i++) {
        const double stress1 = i*stress_increment;
        const double stress2 = stress1+stress_increment;
        const double ratio1 = 3.0*stress1/(sqrt(2.0)*cone_height);
        const double ratio2 = 3.0*stress2/(sqrt(2.0)*cone_height);
        const double strain1 = stress1*ref_strain/
            (c->ref_shear_modulus*ref_strain-stress1);
        const double strain2 = stress2*ref_strain/
            (c->ref_shear_modulus*ref_strain-stress2);
        if (ratio1 <= c->phase_transform_ratio &&
            ratio2 >= c->phase_transform_ratio && ratio2 > ratio1) {
            const double fraction = (ratio2-c->phase_transform_ratio)/
                                    (ratio2-ratio1);
            c->phase_transform_strain =
                strain2-fraction*(strain2-strain1);
        }
    }
    return isfinite(c->phase_transform_strain);
}

/* Automatic-surface moduli are inexpensive to reconstruct and are therefore
 * not embedded as a forty-double array in every nonlinear element. */
PDMY03_HD static inline double pdmy03_surface_modulus(
    const pdmy03_constants_t *c, int surface)
{
    if (surface >= c->surface_count) return 0.0;
    const double stress_increment =
        c->reference_peak_shear/(double)c->surface_count;
    const double stress1 = surface*stress_increment;
    const double stress2 = stress1+stress_increment;
    const double strain1 = stress1*c->reference_strain/
        (c->ref_shear_modulus*c->reference_strain-stress1);
    const double strain2 = stress2*c->reference_strain/
        (c->ref_shear_modulus*c->reference_strain-stress2);
    const double ep_modulus = 2.0*(stress2-stress1)/(strain2-strain1);
    double plastic_modulus = PDMY03_SURFACE_LIMIT;
    if (2.0*c->ref_shear_modulus-ep_modulus > 0.0) {
        plastic_modulus = 2.0*c->ref_shear_modulus*ep_modulus/
            (2.0*c->ref_shear_modulus-ep_modulus);
    }
    if (plastic_modulus < 0.0) plastic_modulus = 0.0;
    return pdmy03_min(plastic_modulus, PDMY03_SURFACE_LIMIT);
}

PDMY03_HD static inline double pdmy03_modulus_factor(
    const pdmy03_constants_t *c, tensor_t stress)
{
    double scale = (pdmy03_mean(stress)-c->residual_pressure)/
                   (-c->ref_pressure-c->residual_pressure);
    scale = pdmy03_max(scale, 1.0e-10);
    scale = pow(scale, c->pressure_dependence);
    return pdmy03_max(scale, 1.0e-10);
}

PDMY03_HD static inline tensor_t pdmy03_elastic_predict(
    const pdmy03_constants_t *c, tensor_t stress, tensor_t strain_increment,
    double *factor_out)
{
    const double factor = pdmy03_modulus_factor(c, stress);
    tensor_t dev = pdmy03_add(pdmy03_deviator(stress),
        pdmy03_scale(pdmy03_deviator(strain_increment),
                     2.0*c->ref_shear_modulus*factor));
    double mean = pdmy03_mean(stress)+3.0*c->ref_bulk_modulus*factor*
                                      pdmy03_mean(strain_increment);
    if (mean > 0.0) mean = 0.0;
    if (factor_out != NULL) *factor_out = factor;
    return pdmy03_from_dev_mean(dev, mean);
}

PDMY03_HD static inline void pdmy03_scale_to_surface(
    const pdmy03_constants_t *c, tensor_t *stress, tensor_t *centers,
    int surface)
{
    for (int iteration = 0; iteration < 20; iteration++) {
        const double difference = pdmy03_yield(c, *stress, centers, surface);
        const double cone = pdmy03_mean(*stress)-c->residual_pressure;
        const double radius = -pdmy03_surface_size(c, surface)*cone;
        if (surface == c->surface_count) {
            if (fabs(difference) <= PDMY03_SMALL*pdmy03_max(radius*radius, 1.0))
                return;
            tensor_t dev = pdmy03_deviator(*stress);
            const double length = sqrt(pdmy03_max(difference+radius*radius, 0.0));
            if (length > PDMY03_SMALL)
                *stress = pdmy03_from_dev_mean(
                    pdmy03_scale(dev, radius/length), pdmy03_mean(*stress));
            return;
        }
        if (difference >= 0.0) return;
        const double shifted_length =
            sqrt(pdmy03_max(radius*radius+difference, 0.0));
        if (!(shifted_length > PDMY03_SMALL)) return;
        tensor_t dev = pdmy03_deviator(*stress);
        tensor_t shifted = pdmy03_sub(dev,
            pdmy03_scale(centers[surface-1], cone));
        double coefficient = (radius-shifted_length)/shifted_length;
        coefficient = pdmy03_max(coefficient, 1.0e-13);
        dev = pdmy03_add(dev, pdmy03_scale(shifted, coefficient));
        *stress = pdmy03_from_dev_mean(dev, pdmy03_mean(*stress));
    }
}

PDMY03_HD static inline int pdmy03_initialize_state(
    const pdmy03_constants_t *c, pdmy03_state_t *state, tensor_t *centers,
    tensor_t *stress)
{
    /* OpenSees elast2Plast() removes tensile mean stress but retains the
     * deviator.  An exactly zero state is a valid elastic/plastic-stage
     * transition and must not become a constitutive failure. */
    if (pdmy03_mean(*stress) > 0.0)
        *stress = pdmy03_from_dev_mean(pdmy03_deviator(*stress), 0.0);
    const pdmy03_state_t empty_state = {};
    *state = empty_state;
    state->initialized = 1;
    state->on_ppz = -1;
    for (int i = 0; i < c->surface_count; i++) centers[i] = pdmy03_zero();

    int active = 0;
    while (active < c->surface_count &&
           pdmy03_yield(c, *stress, centers, active+1) > 0.0) active++;
    if (active == c->surface_count &&
        pdmy03_yield(c, *stress, centers, active) > 0.0)
        pdmy03_scale_to_surface(c, stress, centers, active);
    state->active_surface = active;
    /* OpenSees captures stress0 after elast2Plast(), including its tension
     * clamp and any projection back to the failure surface. */
    state->initial_stress = *stress;
    state->maximum_pressure = pdmy03_mean(*stress);
    if (active == 0) return 1;

    const double signed_cone = pdmy03_mean(*stress)-c->residual_pressure;
    const double cone_height = -signed_cone;
    const tensor_t dev = pdmy03_deviator(*stress);
    const double stress_ratio_length = sqrt(1.5*pdmy03_ddot(dev, dev));
    if (!(cone_height > 0.0) || !(stress_ratio_length > PDMY03_SMALL)) return 0;
    /* Translate all active inner surfaces, but never translate the failure
     * surface away from the hydrostatic axis. */
    for (int surface = 1;
         surface <= active && surface < c->surface_count; surface++) {
        const double coefficient = 1.0-
            pdmy03_surface_size(c, surface)*cone_height/stress_ratio_length;
        centers[surface-1] = pdmy03_scale(dev, coefficient/signed_cone);
    }
    return 1;
}

PDMY03_HD static inline tensor_t pdmy03_surface_normal(
    const pdmy03_constants_t *c, tensor_t stress, const tensor_t *centers,
    int active)
{
    const double cone = pdmy03_mean(stress)-c->residual_pressure;
    const tensor_t center = centers[active-1];
    const tensor_t dev = pdmy03_deviator(stress);
    const double size = pdmy03_surface_size(c, active);
    const double volume = cone*(pdmy03_ddot(center, center)-
                           2.0*size*size/3.0)-pdmy03_ddot(dev, center);
    tensor_t normal = pdmy03_from_dev_mean(
        pdmy03_scale(pdmy03_sub(dev, pdmy03_scale(center, cone)), 3.0),
        volume);
    const double length = pdmy03_norm(normal);
    return length > PDMY03_SMALL ? pdmy03_scale(normal, 1.0/length)
                                : pdmy03_zero();
}

PDMY03_HD static inline tensor_t pdmy03_contact_stress(
    const pdmy03_constants_t *c, tensor_t stress, const tensor_t *centers,
    int active)
{
    const double cone = pdmy03_mean(stress)-c->residual_pressure;
    const tensor_t center = centers[active-1];
    tensor_t shifted = pdmy03_sub(pdmy03_deviator(stress),
                                  pdmy03_scale(center, cone));
    const double length = sqrt(1.5*pdmy03_ddot(shifted, shifted));
    if (!(length > PDMY03_SMALL)) return stress;
    tensor_t dev = pdmy03_add(
        pdmy03_scale(shifted,
            pdmy03_surface_size(c, active)*(-cone)/length),
        pdmy03_scale(center, cone));
    return pdmy03_from_dev_mean(dev, pdmy03_mean(stress));
}

PDMY03_HD static inline double pdmy03_deviator_ratio(
    const pdmy03_constants_t *c, tensor_t stress)
{
    return sqrt(1.5*pdmy03_ddot(pdmy03_deviator(stress),
                                pdmy03_deviator(stress)))/
           pdmy03_max(fabs(pdmy03_mean(stress))+c->residual_pressure,
                      PDMY03_SMALL);
}

PDMY03_HD static inline int pdmy03_is_reversal(
    const pdmy03_constants_t *c, tensor_t normal_stress, tensor_t trial,
    tensor_t committed_stress,
    const tensor_t *centers, int active)
{
    if (active == 0) return 0;
    const tensor_t normal = pdmy03_surface_normal(
        c, normal_stress, centers, active);
    return pdmy03_ddot(pdmy03_sub(trial, committed_stress), normal) < 0.0;
}

PDMY03_HD static inline void pdmy03_update_inner_surfaces(
    const pdmy03_constants_t *c, tensor_t stress, tensor_t *centers,
    int active)
{
    if (active <= 1) return;
    const double cone = pdmy03_mean(stress)-c->residual_pressure;
    if (fabs(cone) <= PDMY03_SMALL) return;
    const tensor_t dev = pdmy03_deviator(stress);
    const tensor_t center = centers[active-1];
    const double outer_size = pdmy03_surface_size(c, active);
    for (int surface = 1; surface < active; surface++) {
        tensor_t value = pdmy03_sub(pdmy03_scale(center, cone), dev);
        value = pdmy03_scale(value,
            pdmy03_surface_size(c, surface)/outer_size);
        centers[surface-1] = pdmy03_scale(pdmy03_add(value, dev), 1.0/cone);
    }
}

PDMY03_HD static inline double pdmy03_contraction_factor(
    const pdmy03_constants_t *c, const pdmy03_state_t *state,
    tensor_t current_stress)
{
    const double p0 = fabs(state->initial_stress.xx+state->initial_stress.yy+
                           state->initial_stress.zz)/3.0+
                      c->residual_pressure;
    const double csr0 = sqrt(state->initial_stress.xy*state->initial_stress.xy+
                             state->initial_stress.yz*state->initial_stress.yz+
                             state->initial_stress.xz*state->initial_stress.xz)/
                          pdmy03_max(p0, PDMY03_SMALL);
    const double csr = sqrt(current_stress.xy*current_stress.xy+
                            current_stress.yz*current_stress.yz+
                            current_stress.xz*current_stress.xz)/
                         pdmy03_max(p0, PDMY03_SMALL);
    return (1.0+pow(fabs(csr-csr0)*c->contraction[3], 3.0))*
           pow(1.0+c->contraction[4]*csr0, 2.0);
}

PDMY03_HD static inline int pdmy03_critical_state(
    const pdmy03_constants_t *c, tensor_t contact, tensor_t updated_stress,
    tensor_t current_strain, tensor_t trial_strain)
{
    const double e_initial = 0.6;
    const double cs1 = 0.9, cs2 = 0.02, cs3 = 0.7;
    const double trial_volume = 3.0*pdmy03_mean(trial_strain);
    const double current_volume = 3.0*pdmy03_mean(current_strain);
    const double e_trial = e_initial+trial_volume*(1.0+e_initial);
    const double e_current = e_initial+current_volume*(1.0+e_initial);
    const double ecr_trial = cs1-cs2*pow(
        fabs(-pdmy03_mean(contact)/c->atmospheric_pressure), cs3);
    const double ecr_current = cs1-cs2*pow(
        fabs(-pdmy03_mean(updated_stress)/c->atmospheric_pressure), cs3);
    if (e_current < ecr_current && e_trial < ecr_trial) return 0;
    if (e_current > ecr_current && e_trial > ecr_trial) return 0;
    return 1;
}

PDMY03_HD static inline void pdmy03_ppz_translate(
    const pdmy03_constants_t *c, pdmy03_state_t *state,
    tensor_t current_stress, tensor_t trial_strain, tensor_t full_increment,
    tensor_t committed_pivot_rate)
{
    if (c->liquefaction[0] == 0.0) return;
    const double denominator = state->maximum_pressure-c->residual_pressure;
    state->damage = 0.0;
    if (fabs(denominator) > PDMY03_SMALL) {
        const double ratio = (state->maximum_pressure-
                              pdmy03_mean(current_stress))/denominator;
        if (ratio > 0.0) state->damage = pow(ratio, 0.25);
    }
    const double reversal = pdmy03_ddot(pdmy03_deviator(full_increment),
                                         pdmy03_deviator(committed_pivot_rate));
    if (reversal < 0.0) {
        const double movement = pdmy03_octahedral_shear(
            pdmy03_sub(pdmy03_deviator(trial_strain),
                       pdmy03_deviator(state->ppz_pivot)), 1);
        const double candidate = state->damage*c->liquefaction[1]*movement;
        if (candidate > state->cumulative_translation)
            state->cumulative_translation = candidate;
    }
}

PDMY03_HD static inline void pdmy03_update_ppz(
    const pdmy03_constants_t *c, pdmy03_state_t *state,
    tensor_t contact_stress, tensor_t current_stress,
    tensor_t current_strain, tensor_t trial_strain,
    tensor_t full_increment, tensor_t sub_increment,
    tensor_t committed_pivot_rate)
{
    const double reversal = pdmy03_ddot(pdmy03_deviator(full_increment),
                                         pdmy03_deviator(committed_pivot_rate));
    if (state->on_ppz < 1) {
        state->damage = 0.0;
        const double denominator = state->maximum_pressure-c->residual_pressure;
        if (fabs(denominator) > PDMY03_SMALL) {
            const double ratio = (state->maximum_pressure-
                                  pdmy03_mean(current_stress))/denominator;
            if (ratio > 0.0) state->damage = pow(ratio, 0.25);
        }
    }
    if (c->liquefaction[0] == 0.0) {
        if (state->on_ppz == 2) {
            state->ppz_pivot = trial_strain;
            state->cumulative_dilation +=
                pdmy03_octahedral_shear(sub_increment, 1);
        } else {
            state->on_ppz = 2;
            state->ppz_pivot = trial_strain;
            state->pivot_strain_rate = pdmy03_deviator(full_increment);
            if (reversal < 0.0) state->cumulative_dilation = 0.0;
        }
        return;
    }
    if (state->on_ppz == 2) {
        state->ppz_pivot = trial_strain;
        const double shear = pdmy03_octahedral_shear(sub_increment, 1);
        state->cumulative_dilation += shear;
        state->maximum_cumulative_dilation +=
            state->damage*c->liquefaction[0]*shear;
        return;
    }
    if ((state->on_ppz == -1 || state->on_ppz == 0) && reversal < 0.0) {
        state->opposite_pre_ppz_strain = state->pre_ppz_strain;
        double ratio = (-pdmy03_mean(contact_stress)+c->residual_pressure)/
                       (c->ref_pressure+c->residual_pressure);
        ratio = pow(pdmy03_max(ratio, 0.0), 1.0-c->pressure_dependence);
        state->pre_ppz_strain = ratio*c->phase_transform_strain;
        if (state->opposite_pre_ppz_strain == 0.0)
            state->opposite_pre_ppz_strain = state->pre_ppz_strain;
    }
    state->ppz_size = 0.5*(state->cumulative_translation+
                           state->maximum_cumulative_dilation);
    if (state->on_ppz == 0 || (state->on_ppz == 1 && reversal < 0.0)) {
        tensor_t direction = pdmy03_sub(state->ppz_pivot, state->ppz_center);
        const double distance = pdmy03_octahedral_shear(direction, 1);
        const double coefficient = distance > 0.0 ?
            (state->ppz_size-state->cumulative_translation)/distance : 0.0;
        state->ppz_center = pdmy03_sub(
            state->ppz_pivot, pdmy03_scale(direction, coefficient));
    }
    const double distance = pdmy03_octahedral_shear(
        pdmy03_sub(trial_strain, state->ppz_center), 1);
    if (distance > state->ppz_size) {
        state->cumulative_dilation = 0.0;
        state->on_ppz = 2;
        state->ppz_pivot = trial_strain;
        state->pivot_strain_rate = pdmy03_deviator(full_increment);
        state->cumulative_translation = 0.0;
    } else {
        if (state->on_ppz == 0 || state->on_ppz == 1)
            pdmy03_ppz_translate(c, state, current_stress, trial_strain,
                                 full_increment, committed_pivot_rate);
        if (state->on_ppz == -1 || state->on_ppz == 0) state->on_ppz = 1;
    }
    (void)current_strain;
}

PDMY03_HD static inline double pdmy03_plastic_potential(
    const pdmy03_constants_t *c, pdmy03_state_t *state,
    tensor_t contact, tensor_t trial_stress, tensor_t surface_normal,
    tensor_t current_stress,
    tensor_t updated_stress, tensor_t current_strain, tensor_t trial_strain,
    tensor_t full_increment, tensor_t sub_increment,
    tensor_t committed_pivot_rate)
{
    double factor_pt = pdmy03_deviator_ratio(c, contact)/
                       pdmy03_max(c->phase_transform_ratio, PDMY03_SMALL);
    const double current_ratio = pdmy03_deviator_ratio(c, updated_stress);
    const double trial_ratio = pdmy03_deviator_ratio(c, trial_stress);
    const double shear_loading = pdmy03_ddot(pdmy03_deviator(updated_stress),
                                             pdmy03_deviator(trial_stress));
    double potential;
    if (factor_pt >= 1.0 && trial_ratio >= current_ratio &&
        shear_loading >= 0.0) {
        pdmy03_update_ppz(c, state, contact, current_stress, current_strain,
                         trial_strain, full_increment, sub_increment,
                         committed_pivot_rate);
        if (state->on_ppz == 1) potential = 0.0;
        else {
            factor_pt -= 1.0;
            const double pressure_factor = pow(
                (fabs(pdmy03_mean(contact))+c->residual_pressure)/
                 c->atmospheric_pressure, -c->dilation[2]);
            potential = pressure_factor*factor_pt*factor_pt*
                (c->dilation[0]+pow(state->cumulative_dilation,
                                     c->dilation[1]));
            potential = pdmy03_min(fabs(potential), 5.0e4);
        }
    } else {
        double angle = 1.0;
        if (current_ratio != 0.0) {
            tensor_t direction = pdmy03_sub(
                pdmy03_scale(pdmy03_deviator(trial_stress),
                    1.0/(fabs(pdmy03_mean(trial_stress))+
                         c->residual_pressure)),
                pdmy03_scale(pdmy03_deviator(updated_stress),
                    1.0/(fabs(pdmy03_mean(updated_stress))+
                         c->residual_pressure)));
            const double denominator = pdmy03_norm(direction)*
                                       pdmy03_norm(pdmy03_deviator(updated_stress));
            if (denominator > PDMY03_SMALL)
                angle = pdmy03_ddot(pdmy03_deviator(updated_stress), direction)/
                        denominator;
        }
        factor_pt = factor_pt*angle-1.0;
        double rule = pow((fabs(pdmy03_mean(contact))+c->residual_pressure)/
                          c->atmospheric_pressure, c->contraction[2]);
        rule = pdmy03_max(rule, 0.1);
        potential = -factor_pt*factor_pt*
            (c->contraction[0]+state->maximum_cumulative_dilation*
             c->contraction[1])*rule*
            pdmy03_contraction_factor(c, state, current_stress);
        if (potential > 0.0) potential = -potential;
        if (state->on_ppz > 0) state->on_ppz = 0;
        if (state->on_ppz != -1)
            pdmy03_ppz_translate(c, state, current_stress, trial_strain,
                                 full_increment, committed_pivot_rate);
    }
    if (pdmy03_critical_state(c, contact, updated_stress,
                              current_strain, trial_strain)) potential = 0.0;
    (void)surface_normal;
    return potential;
}

PDMY03_HD static inline double pdmy03_positive_quadratic_root(
    double a, double b, double c, int context)
{
    if (c == 0.0) return 0.0;
    if (fabs(a) <= PDMY03_SMALL)
        return fabs(b) > PDMY03_SMALL ? -c/b : 0.0;
    const double discriminant = pdmy03_max(b*b-4.0*a*c, 0.0);
    const double root1 = (-b+sqrt(discriminant))/(2.0*a);
    const double root2 = (-b-sqrt(discriminant))/(2.0*a);
    double low = PDMY03_SURFACE_LIMIT, high = -PDMY03_SURFACE_LIMIT;
    if (root1 >= 0.0) { low = root1; high = root1; }
    if (root2 >= 0.0) { low = pdmy03_min(low, root2); high = pdmy03_max(high, root2); }
    if (high < 0.0) return pdmy03_max(root1, root2);
    /* OpenSees secondOrderEqn() always returns the smaller nonnegative root;
     * its final integer argument selects diagnostics, not the root. */
    (void)context;
    return low;
}

PDMY03_HD static inline int pdmy03_update_active_surface(
    const pdmy03_constants_t *c, tensor_t stress, tensor_t *centers,
    int active)
{
    if (active >= c->surface_count) return 1;
    const double cone = pdmy03_mean(stress)-c->residual_pressure;
    tensor_t center = centers[active-1];
    const tensor_t outer_center = centers[active];
    const double size = pdmy03_surface_size(c, active);
    const double outer_size = pdmy03_surface_size(c, active+1);
    tensor_t t1 = pdmy03_sub(pdmy03_deviator(stress),
                             pdmy03_scale(center, cone));
    tensor_t t2 = pdmy03_scale(pdmy03_sub(center, outer_center), cone);
    const double a = pdmy03_ddot(t1, t1);
    const double b = 2.0*pdmy03_ddot(t1, t2);
    const double cc = pdmy03_ddot(t2, t2)-
        2.0*outer_size*outer_size*cone*cone/3.0;
    double x = pdmy03_positive_quadratic_root(a, b, cc, 0);
    if (fabs(x-1.0) < 1.0e-10) x = 1.0;
    if (x < 1.0) return 1;

    tensor_t direction = pdmy03_add(pdmy03_scale(t1, x),
                                    pdmy03_scale(center, cone));
    direction = pdmy03_scale(direction, 1.0-size/outer_size);
    direction = pdmy03_sub(direction,
        pdmy03_scale(pdmy03_sub(center,
            pdmy03_scale(outer_center, size/outer_size)), cone));
    direction = pdmy03_deviator(direction);
    if (pdmy03_norm(direction) < PDMY03_SMALL) return 1;

    const double aa = cone*cone*pdmy03_ddot(direction, direction);
    double bb = 2.0*cone*pdmy03_ddot(t1, direction);
    if (fabs(bb) < PDMY03_SMALL) bb = 0.0;
    const double ccc = pdmy03_ddot(t1, t1)-
        2.0*size*size*cone*cone/3.0;
    if (fabs(ccc) < PDMY03_SMALL ||
        fabs(ccc)/pdmy03_max(pdmy03_ddot(t1, t1), PDMY03_SMALL) <
        PDMY03_SMALL) return 1;
    if (bb > 1.0e-8 || ccc < -1.0e-8) return 0;
    x = pdmy03_positive_quadratic_root(aa, bb, ccc, 1);
    centers[active-1] = pdmy03_sub(center, pdmy03_scale(direction, x));
    return pdmy03_finite_tensor(centers[active-1]);
}

PDMY03_HD static inline int pdmy03_correct_stress(
    const pdmy03_constants_t *c, pdmy03_state_t *state,
    tensor_t *centers, tensor_t *trial_stress, tensor_t updated_stress,
    tensor_t current_stress, tensor_t current_strain, tensor_t trial_strain,
    tensor_t full_increment, tensor_t sub_increment,
    tensor_t committed_pivot_rate, double modulus_factor,
    pdmy03_diagnostics_t *diagnostics)
{
    int crossed = 0;
    while (state->active_surface <= c->surface_count) {
        const int active = state->active_surface;
        tensor_t contact = pdmy03_contact_stress(c, *trial_stress,
                                                 centers, active);
        tensor_t normal = pdmy03_surface_normal(c, contact, centers, active);
        double potential = pdmy03_plastic_potential(
            c, state, contact, *trial_stress, normal, current_stress,
            updated_stress,
            current_strain, trial_strain, full_increment, sub_increment,
            committed_pivot_rate);
        const double normal_dev2 = pdmy03_ddot(pdmy03_deviator(normal),
                                               pdmy03_deviator(normal));
        const double normal_mean = pdmy03_mean(normal);
        const double modulus = pdmy03_surface_modulus(c, active);
        const double term1 = 2.0*c->ref_shear_modulus*modulus_factor*normal_dev2;
        const double term2 = 9.0*c->ref_bulk_modulus*modulus_factor*
                             normal_mean*potential;
        double denominator = term1+term2+modulus*modulus_factor;
        const double limit = ((active == c->surface_count)
            ? pdmy03_surface_modulus(c, active-1) : modulus)*
            modulus_factor/2.0;
        if (denominator < limit &&
            fabs(9.0*c->ref_bulk_modulus*modulus_factor*normal_mean) >
            PDMY03_SMALL) {
            potential = (term2+limit-denominator)/
                (9.0*c->ref_bulk_modulus*modulus_factor*normal_mean);
            denominator = limit;
        }
        tensor_t difference = pdmy03_sub(pdmy03_deviator(*trial_stress),
                                         pdmy03_deviator(contact));
        double loading = denominator > PDMY03_SMALL ?
            pdmy03_ddot(normal, difference)/denominator : 0.0;
        if (loading < 0.0) loading = 0.0;
        const double inner_modulus = active > 1 ?
            pdmy03_surface_modulus(c, active-1) : 0.0;
        if (crossed && inner_modulus > 0.0)
            loading *= (inner_modulus-modulus)/inner_modulus;

        const double old_mean = pdmy03_mean(*trial_stress);
        double mean = old_mean-potential*3.0*c->ref_bulk_modulus*
                                  modulus_factor*loading;
        tensor_t dev = pdmy03_deviator(*trial_stress);
        double shear_coefficient =
            -2.0*c->ref_shear_modulus*modulus_factor*loading;
        if (mean > 0.0 && mean != old_mean) {
            shear_coefficient *= old_mean/(old_mean-mean);
            mean = 0.0;
        } else if (mean > 0.0) mean = 0.0;
        dev = pdmy03_add(dev,
            pdmy03_scale(pdmy03_deviator(normal), shear_coefficient));
        *trial_stress = pdmy03_from_dev_mean(dev, mean);
        pdmy03_scale_to_surface(c, trial_stress, centers, active);

        if (active >= c->surface_count ||
            pdmy03_yield(c, *trial_stress, centers, active+1) <= 0.0) break;
        state->active_surface++;
        crossed = 1;
        diagnostics->crossed_surfaces++;
    }
    return pdmy03_finite_tensor(*trial_stress);
}

PDMY03_HD static inline int pdmy03_requested_substeps(
    const pdmy03_constants_t *c, const pdmy03_state_t *state,
    tensor_t current_stress, tensor_t increment)
{
    if (pdmy03_norm(increment) == 0.0) return 0;
    const double factor = pdmy03_modulus_factor(c, current_stress);
    double ep_modulus = 2.0*c->ref_shear_modulus*factor;
    if (state->active_surface > 0) {
        const double plastic =
            pdmy03_surface_modulus(c, state->active_surface)*factor;
        ep_modulus = 2.0*c->ref_shear_modulus*factor*plastic/
                     pdmy03_max(2.0*c->ref_shear_modulus*factor+plastic,
                                PDMY03_SMALL);
    }
    const tensor_t stress_increment =
        pdmy03_scale(pdmy03_deviator(increment), ep_modulus);
    const double cone = -(pdmy03_mean(current_stress)-c->residual_pressure);
    const double single = c->failure_ratio*cone/
                          (double)c->surface_count;
    const double total = 3.0*pdmy03_octahedral_shear(stress_increment, 0)/
                         sqrt(2.0);
    int count = single > PDMY03_SMALL ? (int)(total/single)+1 : 1;
    if (count > c->surface_count) count = c->surface_count;
    const int shear_count = (int)(pdmy03_octahedral_shear(increment, 1)/1.0e-5);
    const int volume_count = (int)(pdmy03_mean(increment)/1.0e-5);
    if (shear_count > count) count = shear_count;
    if (volume_count > count) count = volume_count;
    return count > 0 ? count : 1;
}

/* One committed, explicit PDMY03 update. Stress uses the Hercules/OpenSees
 * tension-positive convention; compression therefore has negative mean. */
PDMY03_HD static inline int pdmy03_update(
    const pdmy03_constants_t *c, int maximum_substeps,
    tensor_t strain, tensor_t previous_strain, tensor_t static_stress,
    tensor_t *dynamic_stress, tensor_t *plastic_strain,
    pdmy03_state_t *state, tensor_t *centers,
    pdmy03_diagnostics_t *diagnostics)
{
    const pdmy03_diagnostics_t empty_diagnostics = {};
    *diagnostics = empty_diagnostics;
    const tensor_t increment = pdmy03_sub(strain, previous_strain);
    diagnostics->strain_increment_norm = pdmy03_norm(increment);
    tensor_t current = pdmy03_add(*dynamic_stress, static_stress);
    if (diagnostics->strain_increment_norm == 0.0) {
        diagnostics->mean_stress = pdmy03_mean(current);
        return 1;
    }
    if (!state->initialized) {
        if (!pdmy03_initialize_state(c, state, centers, &current)) {
            diagnostics->reason = PDMY03_FAILURE_INITIAL_PRESSURE;
            diagnostics->mean_stress = pdmy03_mean(current);
            return 0;
        }
        *dynamic_stress = pdmy03_sub(current, static_stress);
    }

    double factor = 1.0;
    tensor_t trial = pdmy03_elastic_predict(c, current, increment, &factor);
    if (state->active_surface > 0 &&
        pdmy03_is_reversal(c, current, trial, current, centers,
                           state->active_surface)) {
        pdmy03_update_inner_surfaces(c, current, centers,
                                     state->active_surface);
        state->active_surface = 0;
    }
    if (state->active_surface == 0 &&
        pdmy03_yield(c, trial, centers, 1) <= 0.0) {
        *dynamic_stress = pdmy03_sub(trial, static_stress);
        diagnostics->mean_stress = pdmy03_mean(trial);
        diagnostics->requested_substeps = 1;
        diagnostics->accepted_substeps = 1;
        state->maximum_pressure = pdmy03_min(state->maximum_pressure,
                                             pdmy03_mean(trial));
        return pdmy03_finite_tensor(*dynamic_stress);
    }

    const int requested = pdmy03_requested_substeps(c, state, current,
                                                     increment);
    const int used = requested > maximum_substeps ? maximum_substeps : requested;
    if (used <= 0) {
        diagnostics->reason = PDMY03_FAILURE_CONSTANTS;
        return 0;
    }
    diagnostics->requested_substeps = requested;
    diagnostics->substeps_clamped = requested > used;
    const tensor_t sub_increment = pdmy03_scale(increment, 1.0/(double)used);
    const tensor_t committed_pivot_rate = state->pivot_strain_rate;
    tensor_t plastic = *plastic_strain;
    tensor_t updated = current;

    for (int substep = 0; substep < used; substep++) {
        const tensor_t base = substep == 0 ? current : trial;
        updated = base;
        trial = pdmy03_elastic_predict(c, base, sub_increment, &factor);
        if (state->active_surface > 0 &&
            pdmy03_is_reversal(c, base, trial, current, centers,
                               state->active_surface)) {
            pdmy03_update_inner_surfaces(c, current, centers,
                                         state->active_surface);
            state->active_surface = 0;
        }
        if (state->active_surface == 0 &&
            pdmy03_yield(c, trial, centers, 1) <= 0.0) {
            diagnostics->accepted_substeps++;
            continue;
        }
        if (state->active_surface == 0) state->active_surface = 1;
        const tensor_t trial_strain = pdmy03_add(
            previous_strain, pdmy03_scale(sub_increment, substep+1.0));
        if (!pdmy03_correct_stress(
                c, state, centers, &trial, updated, current,
                previous_strain, trial_strain, increment, sub_increment,
                committed_pivot_rate, factor, diagnostics) ||
            !pdmy03_update_active_surface(c, trial, centers,
                                          state->active_surface)) {
            diagnostics->reason = PDMY03_FAILURE_SURFACE_GEOMETRY;
            diagnostics->mean_stress = pdmy03_mean(trial);
            return 0;
        }
        tensor_t elastic_increment = pdmy03_scale(
            pdmy03_deviator(pdmy03_sub(trial, updated)),
            1.0/(2.0*c->ref_shear_modulus*factor));
        elastic_increment = pdmy03_from_dev_mean(
            elastic_increment,
            (pdmy03_mean(trial)-pdmy03_mean(updated))/
            (3.0*c->ref_bulk_modulus*factor));
        plastic = pdmy03_add(
            plastic, pdmy03_sub(sub_increment, elastic_increment));
        const double bulk = c->ref_bulk_modulus*factor;
        if (bulk > PDMY03_SMALL) {
            state->pressure_d += 3.0*pdmy03_mean(sub_increment)-
                (pdmy03_mean(trial)-pdmy03_mean(updated))/bulk;
            if (state->pressure_d < 0.0) state->pressure_d = 0.0;
        }
        diagnostics->accepted_substeps++;
    }

    state->maximum_pressure = pdmy03_min(state->maximum_pressure,
                                         pdmy03_mean(trial));
    *dynamic_stress = pdmy03_sub(trial, static_stress);
    *plastic_strain = plastic;
    const double active_yield = state->active_surface > 0 ?
        fabs(pdmy03_yield(c, trial, centers, state->active_surface)) : 0.0;
    const double radius = state->active_surface > 0 ?
        fabs(pdmy03_surface_size(c, state->active_surface)*
             (pdmy03_mean(trial)-c->residual_pressure)) : 1.0;
    diagnostics->yield_residual = active_yield/
                                  pdmy03_max(radius*radius, 1.0);
    diagnostics->mean_stress = pdmy03_mean(trial);
    if (!pdmy03_finite_tensor(*dynamic_stress) ||
        !isfinite(diagnostics->yield_residual)) {
        diagnostics->reason = PDMY03_FAILURE_NONFINITE;
        return 0;
    }

    return 1;
}

#undef PDMY03_HD

#endif /* HERCULES_NONLINEAR_PDMY03_H */
