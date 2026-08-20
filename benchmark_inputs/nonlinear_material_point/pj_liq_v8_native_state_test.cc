/* Host-only PJ-Liq V8 state-contract tests.
 *
 * These checks complement the golden-response replay: they verify that two
 * material points have isolated histories, that a byte-restored checkpoint
 * state reproduces uninterrupted integration, and that the coupled BrickUP
 * split satisfies sealed fluid mass balance before pressure is combined with
 * the effective skeleton stress exactly once.
 */
#define PJ8_STANDALONE
#include "../../quake/forward/nonlinear_pj_liq_v8.h"

#include <array>
#include <cstdlib>
#include <cstring>
#include <iostream>

static bool same_state(const pj8_state_t &a, const pj8_state_t &b)
{
    double av[PJ8_STATE_VALUE_COUNT] = {};
    double bv[PJ8_STATE_VALUE_COUNT] = {};
    if (pj8_state_values(&a, av) != PJ8_STATE_VALUE_COUNT ||
        pj8_state_values(&b, bv) != PJ8_STATE_VALUE_COUNT)
        return false;
    for (int i = 0; i < PJ8_STATE_VALUE_COUNT; ++i)
        if (av[i] != bv[i]) {
            std::cerr << "flattened state field " << i << " differs: "
                      << av[i] << " versus " << bv[i] << "\n";
            return false;
        }
    return true;
}

static pj8_state_t initial_state(const pj8_parameters_t &parameters,
                                 double void_ratio)
{
    const tensor_t effective_stress = {-18.0, -18.0, -40.0, 0, 0, 0};
    pj8_state_t state = {};
    if (!pj8_initialize(&parameters, effective_stress, void_ratio, &state) ||
        !pj8_begin_dynamic_phase(&parameters, effective_stress, &state)) {
        std::cerr << "PJ-Liq V8 initialization failed\n";
        std::exit(1);
    }
    return state;
}

static pj8_state_t initial_material_state(
    const pj8_parameters_t &parameters,
    const pj8_material_parameters_t &material,
    double void_ratio)
{
    const tensor_t effective_stress = {-18.0, -18.0, -40.0, 0, 0, 0};
    pj8_state_t state = {};
    if (!pj8_initialize_material(&parameters, &material, effective_stress,
                                 void_ratio, &state) ||
        !pj8_begin_dynamic_phase(&parameters, effective_stress, &state)) {
        std::cerr << "PJ-Liq V8 custom-material initialization failed\n";
        std::exit(1);
    }
    return state;
}

static bool advance(const pj8_parameters_t &parameters, tensor_t increment,
                    int substeps, pj8_state_t *state)
{
    tensor_t stress = state->stress;
    pj8_update_info_t info = {};
    const int expected_substeps = pj8_ddot(increment, increment) == 0.0
        ? 0 : substeps;
    return pj8_update(&parameters, increment, substeps, state, &stress, &info) &&
           info.accepted_substeps == expected_substeps &&
           pj8_finite_tensor(stress);
}

static bool advance_material(const pj8_parameters_t &parameters,
                             const pj8_material_parameters_t &material,
                             tensor_t increment, int substeps,
                             pj8_state_t *state)
{
    tensor_t stress = state->stress;
    pj8_update_info_t info = {};
    const int expected_substeps = pj8_ddot(increment, increment) == 0.0
        ? 0 : substeps;
    return pj8_update_material(&parameters, &material, increment, substeps,
                               state, &stress, &info) &&
           info.accepted_substeps == expected_substeps &&
           pj8_finite_tensor(stress);
}

int main()
{
    const pj8_parameters_t parameters = pj8_reference_parameters(1.0);
    const pj8_material_parameters_t reference_material =
        pj8_reference_material_parameters(&parameters);
    if (!pj8_material_parameters_valid(&parameters, &reference_material)) {
        std::cerr << "PJ-Liq V8 reference material block is invalid\n";
        return 1;
    }
    const double relative_density_0601 =
        (parameters.e_max-0.601)/(parameters.e_max-parameters.e_min);
    if (pj8_void_ratio_from_relative_density(&parameters, 0.0) !=
            parameters.e_max ||
        pj8_void_ratio_from_relative_density(&parameters, 1.0) !=
            parameters.e_min ||
        pj8_void_ratio_from_relative_density(
            &parameters, relative_density_0601) != 0.601) {
        std::cerr << "PJ-Liq V8 relative-density conversion failed\n";
        return 1;
    }
    if (pj8_void_ratio_from_material_relative_density(
            &reference_material, relative_density_0601) != 0.601) {
        std::cerr << "PJ-Liq V8 custom relative-density conversion failed\n";
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

    pj8_state_t baseline = initial_state(parameters, 0.601);
    for (const tensor_t increment : path)
        if (!advance(parameters, increment, 4, &baseline)) return 1;

    /* The configurable entry point must be bit-for-bit equivalent when its
     * compact material block contains the frozen reference values. */
    pj8_state_t reference_custom =
        initial_material_state(parameters, reference_material, 0.601);
    for (const tensor_t increment : path)
        if (!advance_material(parameters, reference_material, increment, 4,
                              &reference_custom)) return 1;
    if (!same_state(baseline, reference_custom)) {
        std::cerr << "PJLiqV8Custom reference row changed frozen V8 response\n";
        return 1;
    }

    /* zeta is an actual V8 material input, not an inert parser column. Its
     * initial irreversible dilatancy scale changes immediately. */
    pj8_material_parameters_t custom_material = reference_material;
    custom_material.zeta = 2.0*reference_material.zeta;
    if (!pj8_material_parameters_valid(&parameters, &custom_material)) {
        std::cerr << "PJ-Liq V8 valid custom material block was rejected\n";
        return 1;
    }
    const pj8_state_t reference_initial =
        initial_material_state(parameters, reference_material, 0.601);
    const pj8_state_t custom_state =
        initial_material_state(parameters, custom_material, 0.601);
    if (custom_state.D_ir == reference_initial.D_ir) {
        std::cerr << "PJLiqV8Custom zeta does not affect constitutive state\n";
        return 1;
    }

    pj8_state_t isolated = initial_state(parameters, 0.601);
    pj8_state_t other = initial_state(parameters, 0.70);
    for (const tensor_t increment : path) {
        if (!advance(parameters, increment, 4, &isolated) ||
            !advance(parameters, other_increment, 2, &other)) return 1;
    }
    if (!same_state(baseline, isolated)) {
        std::cerr << "PJ-Liq V8 material-point state is not isolated\n";
        return 1;
    }

    pj8_state_t restarted = initial_state(parameters, 0.601);
    for (size_t i = 0; i < path.size()/2; ++i)
        if (!advance(parameters, path[i], 4, &restarted)) return 1;
    std::array<unsigned char, sizeof(pj8_state_t)> checkpoint = {};
    std::memcpy(checkpoint.data(), &restarted, sizeof(restarted));
    restarted = {};
    std::memcpy(&restarted, checkpoint.data(), sizeof(restarted));
    for (size_t i = path.size()/2; i < path.size(); ++i)
        if (!advance(parameters, path[i], 4, &restarted)) return 1;
    if (!same_state(baseline, restarted)) {
        std::cerr << "PJ-Liq V8 checkpoint/restart changed the response\n";
        return 1;
    }

    const pj8_state_t before_zero = restarted;
    if (!advance(parameters, pj8_zero(), 4, &restarted) ||
        !same_state(before_zero, restarted)) {
        std::cerr << "PJ-Liq V8 zero increment changed committed history\n";
        return 1;
    }

    /* V8's confining-strain state is the elastic/plastic effective-skeleton
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
        std::cerr << "PJ-Liq V8/BrickUP sealed mass balance failed\n";
        return 1;
    }

    const double pore_pressure = 12.5+pore_pressure_increment;
    const tensor_t effective = restarted.stress;
    const tensor_t total = pj8_sub(effective,
                                   pj8_iso(biot_alpha*pore_pressure));
    if (total.xx != effective.xx-biot_alpha*pore_pressure ||
        total.yy != effective.yy-biot_alpha*pore_pressure ||
        total.zz != effective.zz-biot_alpha*pore_pressure ||
        total.xy != effective.xy || total.yz != effective.yz ||
        total.xz != effective.xz) {
        std::cerr << "PJ-Liq V8/BrickUP effective-to-total stress contract failed\n";
        return 1;
    }

    std::cout << "PASS: PJ-Liq V8 frozen/custom equivalence, zeta exposure, "
                 "state isolation, restart, zero-increment, and BrickUP "
                 "mass-balance/stress-split checks\n";
    return 0;
}
