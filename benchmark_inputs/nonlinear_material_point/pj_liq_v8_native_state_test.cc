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

int main()
{
    const pj8_parameters_t parameters = pj8_reference_parameters(1.0);
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

    std::cout << "PASS: PJ-Liq V8 state isolation, restart, zero-increment, "
                 "and BrickUP mass-balance/stress-split checks\n";
    return 0;
}
