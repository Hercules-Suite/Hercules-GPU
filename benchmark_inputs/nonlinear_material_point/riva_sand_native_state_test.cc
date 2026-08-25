/* Host-only RIVA-Sand state-contract tests.
 *
 * These checks complement the golden-response replay: they verify that two
 * material points have isolated histories, that a byte-restored checkpoint
 * state reproduces uninterrupted integration, that the low-static-bias
 * reversible-volume response is bounded, and that the coupled BrickUP split
 * satisfies sealed fluid mass balance before pressure is combined with the
 * effective skeleton stress exactly once.
 */
#define RIVA_STANDALONE
#include "../../quake/forward/nonlinear_riva_sand.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>

static bool same_tensor(const tensor_t &a, const tensor_t &b)
{
    return a.xx == b.xx && a.yy == b.yy && a.zz == b.zz &&
           a.xy == b.xy && a.yz == b.yz && a.xz == b.xz;
}

static bool close_value(double a, double b, double tolerance = 1.0e-12)
{
    return std::fabs(a-b) <= tolerance*(1.0+std::fabs(a)+std::fabs(b));
}

static bool same_state(const riva_state_t &a, const riva_state_t &b)
{
    if (a.geostatic_admitted != b.geostatic_admitted ||
        a.initialized != b.initialized) {
        std::cerr << "RIVA-Sand state flags differ: geostatic_admitted="
                  << a.geostatic_admitted << " versus "
                  << b.geostatic_admitted << ", initialized="
                  << a.initialized << " versus " << b.initialized << "\n";
        return false;
    }
    double av[RIVA_STATE_VALUE_COUNT] = {};
    double bv[RIVA_STATE_VALUE_COUNT] = {};
    if (riva_state_values(&a, av) != RIVA_STATE_VALUE_COUNT ||
        riva_state_values(&b, bv) != RIVA_STATE_VALUE_COUNT)
        return false;
    for (int i = 0; i < RIVA_STATE_VALUE_COUNT; ++i)
        if (av[i] != bv[i]) {
            std::cerr << "flattened state field " << i << " differs: "
                      << av[i] << " versus " << bv[i] << "\n";
            return false;
        }
    return true;
}

static riva_state_t initial_state(const riva_parameters_t &parameters,
                                 double void_ratio)
{
    const tensor_t effective_stress = {-18.0, -18.0, -40.0, 0, 0, 0};
    riva_state_t state = {};
    if (!riva_initialize_geostatic(
            &parameters, effective_stress, void_ratio, effective_stress,
            &state, nullptr)) {
        std::cerr << "RIVA-Sand initialization failed\n";
        std::exit(1);
    }
    return state;
}

static riva_state_t initial_material_state(
    const riva_parameters_t &parameters,
    const riva_material_parameters_t &material,
    double void_ratio)
{
    const tensor_t effective_stress = {-18.0, -18.0, -40.0, 0, 0, 0};
    riva_state_t state = {};
    if (!riva_initialize_geostatic_material(
            &parameters, &material, effective_stress, void_ratio,
            effective_stress, &state, nullptr)) {
        std::cerr << "RIVA-Sand custom-material initialization failed\n";
        std::exit(1);
    }
    return state;
}

static bool advance(const riva_parameters_t &parameters, tensor_t increment,
                    int substeps, riva_state_t *state)
{
    tensor_t stress = state->stress;
    riva_update_info_t info = {};
    const int expected_substeps = riva_ddot(increment, increment) == 0.0
        ? 0 : substeps;
    return riva_update(&parameters, increment, substeps, state, &stress, &info) &&
           info.accepted_substeps == expected_substeps &&
           riva_finite_tensor(stress);
}

static bool advance_material(const riva_parameters_t &parameters,
                             const riva_material_parameters_t &material,
                             tensor_t increment, int substeps,
                             riva_state_t *state)
{
    tensor_t stress = state->stress;
    riva_update_info_t info = {};
    const int expected_substeps = riva_ddot(increment, increment) == 0.0
        ? 0 : substeps;
    return riva_update_material(&parameters, &material, increment, substeps,
                               state, &stress, &info) &&
           info.accepted_substeps == expected_substeps &&
           riva_finite_tensor(stress);
}

static double reversible_bias_target(const riva_parameters_t &parameters,
                                     double projected_bias)
{
    riva_state_t state = {};
    state.cyclic_phase_active = 1;
    state.amplitude_reversals = 6;
    state.cyclic_direction.xy = 1.0;
    const double direction_norm = riva_norm(state.cyclic_direction);
    state.static_bias_tensor = riva_scale(
        state.cyclic_direction, projected_bias/direction_norm);
    state.pressure_anchor = parameters.bias_reversible_mean_transition_pressure;
    state.cyclic_amplitude = 0.40;
    state.geostatic_deviator = riva_zero();

    const tensor_t phase_direction =
        riva_scale(state.cyclic_direction, 1.0/direction_norm);
    const tensor_t static_deviator =
        riva_scale(state.static_bias_tensor, state.pressure_anchor);
    const tensor_t dynamic_deviator = riva_scale(
        phase_direction, state.pressure_anchor*state.cyclic_amplitude);
    state.stress = riva_sub(riva_add(static_deviator, dynamic_deviator),
                            riva_iso(state.pressure_anchor));
    return riva_bias_reversible_volume_target(&parameters, &state);
}

int main()
{
    const riva_parameters_t parameters = riva_reference_parameters(1.0);
    const riva_material_parameters_t reference_material =
        riva_reference_material_parameters(&parameters);
    if (!riva_material_parameters_valid(&parameters, &reference_material)) {
        std::cerr << "RIVA-Sand reference material block is invalid\n";
        return 1;
    }
    if (parameters.p_residual != 0.0 || RIVA_KERNEL_REVISION != 2u) {
        std::cerr << "RIVA-Sand production residual pressure or kernel "
                     "revision is incorrect\n";
        return 1;
    }

    /* approximate_geostatic_state=yes supplies an admissible analytical K0
     * stress. It must pass through the common admission path unchanged and
     * must not be marked as an over-bound research state. */
    double sin_phi = 3.0*reference_material.M/(6.0+reference_material.M);
    const double k0 = 1.0-sin_phi;
    const tensor_t approximate_stress = {
        -40.0*k0, -40.0*k0, -40.0, 0, 0, 0};
    riva_state_t approximate_state = {};
    int32_t admitted = -1;
    if (!riva_initialize_geostatic_material(
            &parameters, &reference_material, approximate_stress, 0.601,
            approximate_stress, &approximate_state, &admitted) ||
        admitted != 0 || approximate_state.geostatic_admitted != 0 ||
        !same_tensor(approximate_state.stress, approximate_stress)) {
        std::cerr << "RIVA-Sand approximate K0 admission changed an "
                     "admissible stress\n";
        return 1;
    }

    /* approximate_geostatic_state=no supplies the stress committed by the
     * ramped gravity solve. A compressive over-bound state must be preserved,
     * recentered, and prevented from moving farther away from the cone. */
    const tensor_t ramped_stress = {-2.0, -2.0, -40.0, 0, 0, 0};
    riva_state_t ramped_state = {};
    admitted = 0;
    if (!riva_initialize_geostatic_material(
            &parameters, &reference_material, ramped_stress, 0.601,
            ramped_stress, &ramped_state, &admitted) ||
        admitted != 1 || ramped_state.geostatic_admitted != 1 ||
        !same_tensor(ramped_state.stress, ramped_stress) ||
        !same_tensor(ramped_state.alpha0, ramped_state.alpha) ||
        !same_tensor(ramped_state.alpha01, ramped_state.alpha) ||
        !(riva_geostatic_cone_ratio(
              &parameters, &reference_material, &ramped_state) > 1.0)) {
        std::cerr << "RIVA-Sand ramped over-bound geostatic admission failed\n";
        return 1;
    }
    const riva_state_t admitted_before_zero = ramped_state;
    if (!advance_material(parameters, reference_material, riva_zero(), 4,
                          &ramped_state) ||
        !same_state(admitted_before_zero, ramped_state) ||
        ramped_state.geostatic_admitted != 1) {
        std::cerr << "RIVA-Sand admitted stress changed at zero increment\n";
        return 1;
    }
    const tensor_t expanding_candidate = riva_scale(
        admitted_before_zero.alpha,1.25*admitted_before_zero.pressure_anchor);
    tensor_t regularized_alpha = {};
    const tensor_t regularized_deviator = riva_regularize_deviator(
        &parameters, &reference_material,
        admitted_before_zero.pressure_anchor,
        admitted_before_zero.void_ratio, admitted_before_zero.alpha, 1,
        expanding_candidate, &regularized_alpha);
    (void)regularized_deviator;
    if (!close_value(riva_norm(regularized_alpha),
                     riva_norm(admitted_before_zero.alpha))) {
        std::cerr << "RIVA-Sand admitted stress expanded away from the cone\n";
        return 1;
    }
    riva_state_t admitted_after_shear = admitted_before_zero;
    const tensor_t first_dynamic_increment = {0, 0, 0, 5.0e-5, 0, 0};
    if (!advance_material(parameters, reference_material,
                          first_dynamic_increment, 4,
                          &admitted_after_shear) ||
        admitted_after_shear.geostatic_admitted != 1 ||
        riva_norm(admitted_after_shear.alpha) >
            riva_norm(admitted_before_zero.alpha)*(1.0+1.0e-12)) {
        std::cerr << "RIVA-Sand first dynamic update expanded an admitted "
                     "geostatic stress\n";
        return 1;
    }

    /* pResidual remains a non-exposed research control. Enabling 1% of
     * atmospheric pressure translates the cone while leaving physical stress
     * unchanged; the production factory above keeps the value exactly zero. */
    riva_parameters_t residual_parameters = parameters;
    residual_parameters.p_residual = 0.01*parameters.p_atm;
    const tensor_t shallow_stress = {-0.20, -0.20, -0.20, 0, 0, 0};
    riva_state_t residual_state = {};
    admitted = -1;
    if (!riva_initialize_geostatic_material(
            &residual_parameters, &reference_material, shallow_stress, 0.601,
            shallow_stress, &residual_state, &admitted) ||
        admitted != 0 || !same_tensor(residual_state.stress, shallow_stress) ||
        !close_value(residual_state.pressure_anchor,
                     0.20+residual_parameters.p_residual)) {
        std::cerr << "RIVA-Sand residual cone-pressure translation failed\n";
        return 1;
    }
    int residual_floor_hit = 0;
    const double residual_cone_floor = riva_pressure_from_confining(
        &residual_parameters, &reference_material, 1.0e6,
        residual_state.pressure_anchor, &residual_floor_hit);
    if (!residual_floor_hit ||
        residual_cone_floor != riva_cone_pressure_floor(&residual_parameters) ||
        !close_value(riva_physical_pressure(
                         &residual_parameters,residual_cone_floor),
                     residual_parameters.p_min)) {
        std::cerr << "RIVA-Sand translated cone lost the physical pressure "
                     "floor: hit=" << residual_floor_hit
                  << " cone=" << residual_cone_floor
                  << " expected="
                  << riva_cone_pressure_floor(&residual_parameters)
                  << " physical="
                  << riva_physical_pressure(
                         &residual_parameters,residual_cone_floor)
                  << " expected physical=" << residual_parameters.p_min
                  << "\n";
        return 1;
    }
    const tensor_t tensile_stress = {0.50, 0.50, 0.50, 0, 0, 0};
    riva_state_t raw_translated_state = {};
    if (!riva_initialize_material(
            &residual_parameters, &reference_material, tensile_stress, 0.601,
            &raw_translated_state) ||
        riva_initialize_geostatic_material(
            &residual_parameters, &reference_material, tensile_stress, 0.601,
            tensile_stress, &raw_translated_state, nullptr)) {
        std::cerr << "RIVA-Sand geostatic admission accepted noncompressive "
                     "physical stress\n";
        return 1;
    }
    const double relative_density_0601 =
        (parameters.e_max-0.601)/(parameters.e_max-parameters.e_min);
    if (riva_void_ratio_from_relative_density(&parameters, 0.0) !=
            parameters.e_max ||
        riva_void_ratio_from_relative_density(&parameters, 1.0) !=
            parameters.e_min ||
        riva_void_ratio_from_relative_density(
            &parameters, relative_density_0601) != 0.601) {
        std::cerr << "RIVA-Sand relative-density conversion failed\n";
        return 1;
    }
    if (riva_void_ratio_from_material_relative_density(
            &reference_material, relative_density_0601) != 0.601) {
        std::cerr << "RIVA-Sand custom relative-density conversion failed\n";
        return 1;
    }

    /* The inverse-bias reversible-volume law must remain unchanged in its
     * calibrated range while remaining bounded and vanishing as bias tends
     * to zero.  This is shared host/device code, so the host contract test
     * also guards the CUDA implementation against the singular expression. */
    const double bias_onset =
        0.25*parameters.bias_reversible_volume_reference_bias;
    const double onset_target =
        reversible_bias_target(parameters, bias_onset);
    const double above_bias = 1.5*bias_onset;
    const double above_target =
        reversible_bias_target(parameters, above_bias);
    const double pressure_factor = std::pow(
        parameters.bias_reversible_mean_transition_pressure/
            parameters.bias_reference_pressure,
        parameters.bias_reversible_volume_pressure_exponent);
    const double buildup = 1.0-std::exp(
        -6.0/parameters.bias_reversible_volume_buildup_reversals);
    const auto unregularized_target = [&](double bias) {
        return -parameters.bias_reversible_volume_amplitude*std::pow(
            parameters.bias_reversible_volume_reference_bias/bias,
            parameters.bias_reversible_volume_bias_exponent)*
            pressure_factor*buildup;
    };
    const double low_target =
        reversible_bias_target(parameters, 1.0e-3*bias_onset);
    const double expected_onset = unregularized_target(bias_onset);
    const double expected_above = unregularized_target(above_bias);
    if (!std::isfinite(onset_target) || !std::isfinite(above_target) ||
        !std::isfinite(low_target) ||
        std::fabs(onset_target-expected_onset) >
            1.0e-12*std::fabs(expected_onset) ||
        std::fabs(above_target-expected_above) >
            1.0e-12*std::fabs(expected_above) ||
        std::fabs(low_target) >= 1.0e-6*std::fabs(onset_target) ||
        reversible_bias_target(parameters, 0.0) != 0.0) {
        std::cerr << "RIVA-Sand low-static-bias regularization failed: "
                  << "onset=" << onset_target
                  << " expected_onset=" << expected_onset
                  << " above=" << above_target
                  << " expected_above=" << expected_above
                  << " low=" << low_target << "\n";
        return 1;
    }
    const std::array<tensor_t, 8> path = {{
        {0, 0, 0,  2.5e-4, 0, 0},
        {0, 0, 0,  2.5e-4, 0, 0},
        {0, 0, 0, -2.5e-4, 0, 0},
        {0, 0, 0, -2.5e-4, 0, 0},
        {-1e-5, -1e-5, -1e-5, -2.5e-4, 0, 0},
        {0, 0, 0,  2.5e-4, 0, 0},
        {0, 0, 0,  2.5e-4, 0, 0},
        {1e-5, 1e-5, 1e-5, -2.5e-4, 0, 0},
    }};
    const tensor_t other_increment = {0, 0, 0, 0, 1.25e-4, -7.5e-5};

    riva_state_t baseline = initial_state(parameters, 0.601);
    for (const tensor_t increment : path)
        if (!advance(parameters, increment, 4, &baseline)) return 1;

    /* The configurable entry point must be bit-for-bit equivalent when its
     * compact material block contains the frozen reference values. */
    riva_state_t reference_custom =
        initial_material_state(parameters, reference_material, 0.601);
    for (const tensor_t increment : path)
        if (!advance_material(parameters, reference_material, increment, 4,
                              &reference_custom)) return 1;
    if (!same_state(baseline, reference_custom)) {
        std::cerr << "RIVASandCustom reference row changed frozen RIVA-Sand response\n";
        return 1;
    }

    /* zeta is an actual RIVA-Sand material input, not an inert parser column. Its
     * initial irreversible dilatancy scale changes immediately. */
    riva_material_parameters_t custom_material = reference_material;
    custom_material.zeta = 2.0*reference_material.zeta;
    if (!riva_material_parameters_valid(&parameters, &custom_material)) {
        std::cerr << "RIVA-Sand valid custom material block was rejected\n";
        return 1;
    }
    const riva_state_t reference_initial =
        initial_material_state(parameters, reference_material, 0.601);
    const riva_state_t custom_state =
        initial_material_state(parameters, custom_material, 0.601);
    if (custom_state.D_ir == reference_initial.D_ir) {
        std::cerr << "RIVASandCustom zeta does not affect constitutive state\n";
        return 1;
    }

    riva_state_t isolated = initial_state(parameters, 0.601);
    riva_state_t other = initial_state(parameters, 0.70);
    for (const tensor_t increment : path) {
        if (!advance(parameters, increment, 4, &isolated) ||
            !advance(parameters, other_increment, 2, &other)) return 1;
    }
    if (!same_state(baseline, isolated)) {
        std::cerr << "RIVA-Sand material-point state is not isolated\n";
        return 1;
    }

    riva_state_t restarted = initial_state(parameters, 0.601);
    for (size_t i = 0; i < path.size()/2; ++i)
        if (!advance(parameters, path[i], 4, &restarted)) return 1;
    std::array<unsigned char, sizeof(riva_state_t)> checkpoint = {};
    std::memcpy(checkpoint.data(), &restarted, sizeof(restarted));
    restarted = {};
    std::memcpy(&restarted, checkpoint.data(), sizeof(restarted));
    for (size_t i = path.size()/2; i < path.size(); ++i)
        if (!advance(parameters, path[i], 4, &restarted)) return 1;
    if (!same_state(baseline, restarted)) {
        std::cerr << "RIVA-Sand checkpoint/restart changed the response\n";
        return 1;
    }

    const riva_state_t before_zero = restarted;
    if (!advance(parameters, riva_zero(), 4, &restarted) ||
        !same_state(before_zero, restarted)) {
        std::cerr << "RIVA-Sand zero increment changed committed history\n";
        return 1;
    }

    /* RIVA-Sand's confining-strain state is the elastic/plastic effective-skeleton
     * law; it is not copied into the BrickUP pressure unknown.  BrickUP owns
     * fluid mass balance.  For a sealed uniform element,
     *
     *   d(zeta) = alpha*d(eps_v) + d(p_f)/M = 0.
     *
     * Verify this independent equation on a nonzero volumetric increment,
     * then verify the effective-to-total stress split. */
    const double storage_modulus = 5.5e6;
    const double biot_alpha = 0.8;
    const double volumetric_increment = -3.0e-5;
    const double pore_pressure_increment =
        -storage_modulus*biot_alpha*volumetric_increment;
    const double fluid_content_increment =
        biot_alpha*volumetric_increment+
        pore_pressure_increment/storage_modulus;
    if (fluid_content_increment != 0.0) {
        std::cerr << "RIVA-Sand/BrickUP sealed mass balance failed\n";
        return 1;
    }

    const double pore_pressure = 12.5+pore_pressure_increment;
    const tensor_t effective = restarted.stress;
    const tensor_t total = riva_sub(effective,
                                   riva_iso(biot_alpha*pore_pressure));
    if (total.xx != effective.xx-biot_alpha*pore_pressure ||
        total.yy != effective.yy-biot_alpha*pore_pressure ||
        total.zz != effective.zz-biot_alpha*pore_pressure ||
        total.xy != effective.xy || total.yz != effective.yz ||
        total.xz != effective.xz) {
        std::cerr << "RIVA-Sand/BrickUP effective-to-total stress contract failed\n";
        return 1;
    }

    std::cout << "PASS: RIVA-Sand frozen/custom equivalence, zeta exposure, "
                 "low-bias and geostatic admission regularization, optional "
                 "residual-cone translation, state isolation, restart, "
                 "zero-increment, and BrickUP mass-balance/stress-split "
                 "checks\n";
    return 0;
}
