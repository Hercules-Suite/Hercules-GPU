/* Host state-contract checks for RIVASand02.
 *
 * The test uses the frozen low-bias input path to exercise byte restart at a
 * static preload checkpoint, immediately after one-time activation, and next
 * to the first host reversal.  It also verifies point-state isolation, the
 * 63-field/138-value flatten contract, activation idempotence after restart,
 * and BrickUP ownership of fluid pressure and total-stress assembly.
 *
 * Usage:
 *   riva_sand_intermediate_bias_research_state_test INPUT_PATH_DIR
 */
#define RIVA_STANDALONE
#include "../../quake/forward/nonlinear_riva_sand_intermediate_bias_research.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using Row = std::unordered_map<std::string, std::string>;

struct Increment {
    std::string phase;
    tensor_t deps;
    int32_t nsub;
    bool activate_before;
};

static std::vector<std::string> split_csv(const std::string &line)
{
    std::vector<std::string> values;
    std::stringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ',')) values.push_back(value);
    if (!line.empty() && line.back() == ',') values.emplace_back();
    if (!values.empty() && !values.back().empty() && values.back().back() == '\r')
        values.back().pop_back();
    return values;
}

static tensor_t read_tensor(const Row &row, const std::string &prefix)
{
    const auto value = [&](const char *component) {
        return std::stod(row.at(prefix + "_" + component));
    };
    return {value("xx"), value("yy"), value("zz"), value("xy"),
            value("yz"), value("xz")};
}

static std::vector<Increment> read_path(const std::string &path)
{
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open " + path);
    std::string line;
    if (!std::getline(input, line)) throw std::runtime_error("empty " + path);
    const std::vector<std::string> header = split_csv(line);
    std::vector<Increment> increments;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != header.size())
            throw std::runtime_error("invalid CSV width in " + path);
        Row row;
        for (size_t i = 0; i < header.size(); ++i) row[header[i]] = values[i];
        increments.push_back({row.at("phase"), read_tensor(row, "deps"),
                              std::stoi(row.at("nsub")),
                              std::stoi(row.at("activate_before")) != 0});
    }
    return increments;
}

static bool same_tensor(const tensor_t &a, const tensor_t &b)
{
    return a.xx == b.xx && a.yy == b.yy && a.zz == b.zz &&
           a.xy == b.xy && a.yz == b.yz && a.xz == b.xz;
}

static bool same_state(const riva_ib_state_t &a, const riva_ib_state_t &b,
                       const char *context)
{
    if (a.base.initialized != b.base.initialized ||
        a.base.geostatic_admitted != b.base.geostatic_admitted) {
        std::cerr << context << ": lifecycle flags differ\n";
        return false;
    }
    double av[RIVA_IB_STATE_VALUE_COUNT] = {};
    double bv[RIVA_IB_STATE_VALUE_COUNT] = {};
    const int acount = riva_ib_state_values(&a, av);
    const int bcount = riva_ib_state_values(&b, bv);
    if (acount != RIVA_IB_STATE_VALUE_COUNT ||
        bcount != RIVA_IB_STATE_VALUE_COUNT) {
        std::cerr << context << ": flattened state count differs\n";
        return false;
    }
    for (int i = 0; i < RIVA_IB_STATE_VALUE_COUNT; ++i) {
        if (av[i] != bv[i]) {
            std::cerr << context << ": flattened field " << i << " differs: "
                      << av[i] << " versus " << bv[i] << '\n';
            return false;
        }
    }
    return true;
}

static bool close_value(double a, double b, double tolerance = 1.0e-14)
{
    return std::fabs(a - b) <=
           tolerance * (1.0 + std::fabs(a) + std::fabs(b));
}

static riva_ib_state_t initialize(const riva_ib_parameters_t &parameters,
                                  const riva_material_parameters_t &material,
                                  double void_ratio = 0.601)
{
    const tensor_t initial_stress = {-19.4, -19.4, -40.0, 0.0, 0.0, 0.0};
    riva_ib_state_t state = {};
    if (!riva_ib_initialize_material(&parameters, &material, initial_stress,
                                     void_ratio, &state))
        throw std::runtime_error("native research initialization failed");
    return state;
}

static bool advance(const riva_ib_parameters_t &parameters,
                    const riva_material_parameters_t &material,
                    const Increment &increment, riva_ib_state_t *state,
                    riva_update_info_t *info)
{
    if (increment.activate_before &&
        !riva_ib_begin_dynamic_phase(&parameters, &material, nullptr, state))
        return false;
    tensor_t stress = state->base.stress;
    return riva_ib_update_material(&parameters, &material, increment.deps,
                                   increment.nsub, state, &stress, info) &&
           info->accepted_substeps == increment.nsub &&
           same_tensor(stress, state->base.stress) &&
           riva_finite_tensor(stress);
}

static riva_ib_state_t run_history(
    const riva_ib_parameters_t &parameters,
    const riva_material_parameters_t &material,
    const std::vector<Increment> &increments, int restart_after,
    bool verify_active_restart)
{
    riva_ib_state_t state = initialize(parameters, material);
    for (size_t i = 0; i < increments.size(); ++i) {
        riva_update_info_t info = {};
        if (!advance(parameters, material, increments[i], &state, &info))
            throw std::runtime_error("update failed at input row " +
                                     std::to_string(i + 1));
        if (static_cast<int>(i) == restart_after) {
            std::array<unsigned char, sizeof(riva_ib_state_t)> bytes = {};
            std::memcpy(bytes.data(), &state, sizeof(state));
            const riva_ib_state_t before = state;
            state = {};
            std::memcpy(&state, bytes.data(), sizeof(state));
            if (!same_state(before, state, "byte restart"))
                throw std::runtime_error("byte restart did not preserve state");
            if (verify_active_restart) {
                if (!state.base.cyclic_phase_active)
                    throw std::runtime_error("post-activation restart is inactive");
                const riva_ib_state_t active_before = state;
                if (!riva_ib_begin_dynamic_phase(&parameters, &material, nullptr,
                                                 &state) ||
                    !same_state(active_before, state,
                                "idempotent restart activation"))
                    throw std::runtime_error(
                        "active restart re-anchored dynamic memories");
            }
        }
    }
    return state;
}

static riva_ib_state_t biased_state(const riva_ib_parameters_t &parameters,
                                    const riva_material_parameters_t &material,
                                    double void_ratio, double shear)
{
    const tensor_t reference = {-19.4, -19.4, -40.0, 0, 0, 0};
    tensor_t stress = reference;
    stress.xz = shear;
    riva_ib_state_t state = {};
    if (!riva_ib_initialize_material(&parameters, &material, stress, void_ratio,
                                     &state) ||
        !riva_ib_begin_dynamic_phase(&parameters, &material, &reference, &state))
        throw std::runtime_error("biased gate fixture failed to initialize");
    if (state.ep_half_last != 0.0 || state.base.amplitude_reversals != 0)
        throw std::runtime_error("activation fabricated completed cyclic activity");
    return state;
}

static void check_plastic_gate(const riva_ib_parameters_t &parameters,
                               const riva_material_parameters_t &material)
{
    if (parameters.bias_reversible_volume_ep_ref != 1.5e-5)
        throw std::runtime_error("plastic activity reference changed");
    riva_ib_parameters_t isolated = parameters;
    isolated.phase_transformation_enabled = 0;
    riva_ib_state_t state = biased_state(isolated, material, 0.601, 10.0);
    state.base.cyclic_direction = {0, 0, 0, 0, 0, 1};
    state.base.amplitude_reversals = 6;
    state.base.cyclic_amplitude = 0.2;
    const double raw = riva_bias_reversible_volume_target(&isolated.base,
                                                         &state.base);
    if (std::fabs(raw) <= 1.0e-12)
        throw std::runtime_error("gate test requires a nonzero inherited target");
    const double fractions[] = {0.0, 0.25, 0.5, 1.0, 2.0};
    const double factors[] = {0.0, 0.15625, 0.5, 1.0, 1.0};
    for (size_t i = 0; i < 5; ++i) {
        state.ep_half_last = fractions[i]*parameters.bias_reversible_volume_ep_ref;
        if (!close_value(riva_ib_bias_volume_target(&isolated, &state, 0),
                         raw*factors[i], 1.0e-18))
            throw std::runtime_error("plastic gate endpoints or transition failed");
    }
    // The independent PT volume remains present when inherited activity is zero.
    state.ep_half_last = 0.0;
    state.phase_reversible_volume = 2.5e-5;
    if (riva_ib_bias_volume_target(&parameters, &state, 1) !=
        state.phase_reversible_volume)
        throw std::runtime_error("gate suppressed independent PT volume");

    // Reused point storage must not retain an activity scalar from a prior soil.
    state.ep_half_last = 1.0;
    if (!riva_ib_initialize_material(&parameters, &material, state.base.stress,
                                     0.601, &state) || state.ep_half_last != 0.0)
        throw std::runtime_error("reinitialization retained plastic activity");

    // Exercise both real backbones, checking the capture before reset and its
    // preservation across every subsequent constitutive substep in that host step.
    for (bool mapping : {false, true}) {
        for (int nsub : {1, 4, 16}) {
            for (double amplitude : {1.0e-5, 0.01}) {
                riva_ib_parameters_t p = parameters;
                if (!mapping) {
                    p.phase_transformation_enabled = 0;
                    p.mapping_backstress_enabled = 0;
                    p.intermediate_bias_flow_enabled = 0;
                }
                riva_ib_state_t s = biased_state(
                    p, material, mapping ? 0.662 : 0.601, mapping ? 15.0 : 10.0);
                if (mapping && riva_ib_mapping_gate(&p, &s) <= 1.0e-14)
                    throw std::runtime_error("mapping test did not activate mapping");
                int reversals = 0;
                double previous = 0.0;
                for (int step = 1; step <= 64; ++step) {
                    const double gamma = amplitude*std::sin(
                        2.0*3.14159265358979323846*step/32.0);
                    const riva_ib_state_t before = s;
                    const Increment increment = {
                        "dynamic", {0, 0, 0, 0, 0, 0.5*(gamma-previous)},
                        nsub, false};
                    riva_update_info_t info = {};
                    if (!advance(p, material, increment, &s, &info))
                        throw std::runtime_error("gate cycle integration failed");
                    const double expected = info.reversal_registered
                        ? before.base.ep_eq_since_reversal : before.ep_half_last;
                    if (s.ep_half_last != expected || s.ep_half_last < 0.0 ||
                        !std::isfinite(s.ep_half_last))
                        throw std::runtime_error("half-cycle capture/reset order failed");
                    if (info.reversal_registered) ++reversals;
                    previous = gamma;
                }
                if (reversals != 4)
                    throw std::runtime_error("gate fixture lost a host reversal");
                if (amplitude == 0.01 &&
                    s.ep_half_last <= p.bias_reversible_volume_ep_ref)
                    throw std::runtime_error("plastic cycle did not saturate gate");
                if (!mapping && amplitude == 1.0e-5 &&
                    (s.ep_half_last >= p.bias_reversible_volume_ep_ref ||
                     std::fabs(s.base.bias_reversible_volume) >= 1.0e-10))
                    throw std::runtime_error("gate did not suppress weak microcycles");
                double values[RIVA_IB_STATE_VALUE_COUNT] = {};
                if (riva_ib_state_values(&s, values) != RIVA_IB_STATE_VALUE_COUNT ||
                    values[RIVA_IB_STATE_VALUE_COUNT-1] != s.ep_half_last)
                    throw std::runtime_error("activity missing from flattened restart");
                const riva_ib_state_t active = s;
                if (!riva_ib_begin_dynamic_phase(&p, &material, nullptr, &s) ||
                    !same_state(active, s, "active plastic-memory restart"))
                    throw std::runtime_error("active restart reset plastic memory");
            }
        }
    }
    std::cout << "PASS: plastic gate endpoints, PT isolation, microcycles, both "
                 "backbones, and capture-once ordering with 1/4/16 substeps\n";
}

int main(int argc, char **argv)
{
    if (argc != 2) {
        std::cerr << "usage: " << argv[0] << " INPUT_PATH_DIR\n";
        return 2;
    }
    try {
        if (RIVA_IB_PARAMETER_COUNT != 249 ||
            RIVA_IB_LOGICAL_STATE_COUNT != 63 ||
            RIVA_IB_STATE_VALUE_COUNT != 138 || RIVA_IB_KERNEL_REVISION != 4u) {
            throw std::runtime_error("compile-time parameter/state contract changed");
        }
        const riva_ib_parameters_t parameters = riva_ib_reference_parameters(1.0);
        const riva_material_parameters_t material =
            riva_reference_material_parameters(&parameters.base);
        check_plastic_gate(parameters, material);
        riva_ib_state_t initial = initialize(parameters, material);
        double flattened[RIVA_IB_STATE_VALUE_COUNT] = {};
        if (riva_ib_state_values(&initial, flattened) !=
            RIVA_IB_STATE_VALUE_COUNT)
            throw std::runtime_error("runtime state flatten count changed");

        const std::string directory(argv[1]);
        const std::vector<Increment> increments = read_path(
            directory + "/3484_intermediate_bias025.csv");
        if (increments.size() != 484)
            throw std::runtime_error("unexpected low-bias path length");

        int activation_index = -1;
        int first_reversal_index = -1;
        riva_ib_state_t baseline = initialize(parameters, material);
        for (size_t i = 0; i < increments.size(); ++i) {
            if (increments[i].activate_before) activation_index = static_cast<int>(i);
            riva_update_info_t info = {};
            if (!advance(parameters, material, increments[i], &baseline, &info))
                throw std::runtime_error("baseline update failed");
            if (first_reversal_index < 0 && info.reversal_registered)
                first_reversal_index = static_cast<int>(i);
        }
        if (activation_index != 100 || first_reversal_index <= activation_index)
            throw std::runtime_error("activation/reversal fixture changed");

        const riva_ib_state_t static_restart = run_history(
            parameters, material, increments, 49, false);
        if (!same_state(baseline, static_restart, "static-preload restart"))
            return 1;

        const riva_ib_state_t active_restart = run_history(
            parameters, material, increments, activation_index, true);
        if (!same_state(baseline, active_restart, "post-activation restart"))
            return 1;

        const riva_ib_state_t reversal_restart = run_history(
            parameters, material, increments, first_reversal_index - 1, true);
        if (!same_state(baseline, reversal_restart, "reversal-adjacent restart"))
            return 1;

        const riva_ib_state_t completed_restart = run_history(
            parameters, material, increments, first_reversal_index, true);
        if (!same_state(baseline, completed_restart, "completed-half-cycle restart"))
            return 1;

        /* Interleave a different point history.  The primary point must remain
         * bit-identical to its uninterrupted baseline. */
        riva_ib_state_t isolated = initialize(parameters, material);
        riva_ib_state_t other = initialize(parameters, material, 0.662);
        const Increment other_increment = {
            "dynamic", {0.0, 0.0, 0.0, 0.0, 1.25e-5, -7.5e-6}, 2, false};
        bool other_activated = false;
        for (const Increment &increment : increments) {
            riva_update_info_t primary_info = {};
            if (!advance(parameters, material, increment, &isolated,
                         &primary_info))
                throw std::runtime_error("isolated primary update failed");
            if (increment.activate_before && !other_activated) {
                if (!riva_ib_begin_dynamic_phase(&parameters, &material, nullptr,
                                                 &other))
                    throw std::runtime_error("other point activation failed");
                other_activated = true;
            }
            riva_update_info_t other_info = {};
            if (!advance(parameters, material, other_increment, &other,
                         &other_info))
                throw std::runtime_error("other point update failed");
        }
        if (!same_state(baseline, isolated, "integration-point isolation"))
            return 1;
        if (same_tensor(isolated.base.stress, other.base.stress) &&
            isolated.base.lambda_total == other.base.lambda_total &&
            isolated.base.void_ratio == other.base.void_ratio)
            throw std::runtime_error("different point histories aliased");

        /* An affine displacement field produces one homogeneous strain at all
         * eight Hex8 quadrature points.  Independent point storage must then
         * reproduce the same material-point result at every point. */
        std::array<riva_ib_state_t, 8> homogeneous_points;
        homogeneous_points.fill(active_restart);
        const tensor_t affine_increment = {
            1.0e-6, -2.0e-6, 0.5e-6, 0.75e-6, -0.25e-6, 3.0e-6};
        riva_ib_state_t homogeneous_reference = active_restart;
        riva_update_info_t homogeneous_reference_info = {};
        tensor_t homogeneous_stress = homogeneous_reference.base.stress;
        if (!riva_ib_update_material(
                &parameters, &material, affine_increment, 4,
                &homogeneous_reference, &homogeneous_stress,
                &homogeneous_reference_info))
            throw std::runtime_error("homogeneous reference update failed");
        for (size_t point = 0; point < homogeneous_points.size(); ++point) {
            riva_update_info_t point_info = {};
            tensor_t point_stress = homogeneous_points[point].base.stress;
            if (!riva_ib_update_material(
                    &parameters, &material, affine_increment, 4,
                    &homogeneous_points[point], &point_stress, &point_info) ||
                !same_state(homogeneous_reference, homogeneous_points[point],
                            "homogeneous Hex8 point") ||
                !same_tensor(homogeneous_stress, point_stress))
                throw std::runtime_error(
                    "homogeneous Hex8 point differs from material replay");
        }

        /* The field adapter suppresses a zero strain increment even when a
         * rigid translation makes nodal displacement nonzero. */
        const tensor_t rigid_increment = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        for (size_t point = 0; point < homogeneous_points.size(); ++point) {
            const riva_ib_state_t before_rigid = homogeneous_points[point];
            if (riva_ddot(rigid_increment, rigid_increment) > 0.0) {
                riva_update_info_t rigid_info = {};
                tensor_t rigid_stress = homogeneous_points[point].base.stress;
                if (!riva_ib_update_material(
                        &parameters, &material, rigid_increment, 4,
                        &homogeneous_points[point], &rigid_stress,
                        &rigid_info))
                    throw std::runtime_error("rigid-motion update failed");
            }
            if (!same_state(before_rigid, homogeneous_points[point],
                            "rigid Hex8 motion"))
                throw std::runtime_error(
                    "rigid Hex8 motion changed constitutive history");
        }

        /* BrickUP owns pore pressure, fluid mass balance, and total stress.
         * The constitutive call receives no pore-pressure variable and returns
         * effective skeleton stress only. */
        double pore_pressure = 12.5;
        const double pore_pressure_before = pore_pressure;
        riva_ib_state_t skeleton = active_restart;
        const Increment skeleton_increment = {
            "dynamic", {0.0, 0.0, -3.0e-6, 0.0, 0.0, 2.0e-5}, 4, false};
        riva_update_info_t skeleton_info = {};
        tensor_t effective_stress = skeleton.base.stress;
        if (!riva_ib_update_material(&parameters, &material,
                                     skeleton_increment.deps,
                                     skeleton_increment.nsub, &skeleton,
                                     &effective_stress, &skeleton_info) ||
            pore_pressure != pore_pressure_before ||
            !same_tensor(effective_stress, skeleton.base.stress))
            throw std::runtime_error("material did not preserve BrickUP ownership");

        const double biot_alpha = 0.8;
        const double storage_modulus = 5.5e6;
        const double volumetric_increment = -3.0e-5;
        const double pore_pressure_increment =
            -storage_modulus * biot_alpha * volumetric_increment;
        const double fluid_content_increment =
            biot_alpha * volumetric_increment +
            pore_pressure_increment / storage_modulus;
        if (!close_value(fluid_content_increment, 0.0))
            throw std::runtime_error("BrickUP sealed mass balance failed");
        pore_pressure += pore_pressure_increment;
        const tensor_t total_stress = riva_sub(
            effective_stress, riva_iso(biot_alpha * pore_pressure));
        if (!close_value(total_stress.xx,
                         effective_stress.xx - biot_alpha * pore_pressure) ||
            !close_value(total_stress.yy,
                         effective_stress.yy - biot_alpha * pore_pressure) ||
            !close_value(total_stress.zz,
                         effective_stress.zz - biot_alpha * pore_pressure) ||
            total_stress.xy != effective_stress.xy ||
            total_stress.yz != effective_stress.yz ||
            total_stress.xz != effective_stress.xz)
            throw std::runtime_error("effective-to-total stress assembled incorrectly");

        std::cout
            << "PASS: RIVASand02 249-parameter/63-state/"
               "138-value flatten contract, static/post-activation/reversal restart, "
               "activation idempotence, point isolation, homogeneous/rigid "
               "Hex8 behavior, and BrickUP effective-stress ownership\n";
    } catch (const std::exception &error) {
        std::cerr << "RIVASand02 state test failed: "
                  << error.what() << '\n';
        return 1;
    }
    return 0;
}
