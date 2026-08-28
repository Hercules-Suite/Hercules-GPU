/* Standalone prescribed-strain replay for the intermediate-bias RIVA-Sand
 * research successor.  The executable consumes only frozen input_paths CSVs
 * and writes the complete 235-column oracle schema.  It never reads golden
 * response files and has no Python runtime dependency.
 *
 * Usage:
 *   riva_sand_intermediate_bias_research_native_driver INPUT_DIR OUTPUT_DIR
 */
#define RIVA_STANDALONE
#include "../../quake/forward/nonlinear_riva_sand_intermediate_bias_research.h"

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;
using Row = std::unordered_map<std::string, std::string>;

struct CaseDefinition {
    const char *name;
    double vertical_effective_stress;
    double void_ratio;
};

static constexpr double kK0 = 0.485;
static const CaseDefinition kCases[] = {
    {"3484_zero_bias", 40.0, 0.601},
    {"3484_intermediate_bias025", 40.0, 0.601},
    {"3484_intermediate_bias030", 100.0, 0.601},
    {"4666_loose_bias015", 40.0, 0.643},
    {"4666_dense_zero_bias", 40.0, 0.538},
    {"4666_loose_high_bias", 40.0, 0.662},
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

static void put_text(Row &row, const std::string &name, const std::string &value)
{
    row[name] = value;
}

template <class Value>
static void put(Row &row, const std::string &name, Value value)
{
    std::ostringstream text;
    text << std::setprecision(17) << value;
    row[name] = text.str();
}

static void put_tensor(Row &row, const std::string &prefix, tensor_t value)
{
    put(row, prefix + "_xx", value.xx);
    put(row, prefix + "_yy", value.yy);
    put(row, prefix + "_zz", value.zz);
    put(row, prefix + "_xy", value.xy);
    put(row, prefix + "_yz", value.yz);
    put(row, prefix + "_xz", value.xz);
}

static void append_tensor_columns(std::vector<std::string> &columns,
                                  const std::string &prefix)
{
    for (const char *component : {"xx", "yy", "zz", "xy", "yz", "xz"})
        columns.push_back(prefix + "_" + component);
}

static void append_state_columns(std::vector<std::string> &columns)
{
    for (const char *name : {"stress", "alpha", "alpha0", "alpha01", "n",
                             "fabric"})
        append_tensor_columns(columns, std::string("state_") + name);
    for (const char *name : {
             "D", "beta", "lambda_total", "ep_eq_since_reversal", "void_ratio",
             "reversals", "pressure_floor_hits", "denominator_floor_hits",
             "beta_fallbacks", "eps_v_total", "eps_v_confining",
             "eps_v_irreversible", "eps_v_reversible", "pressure_anchor",
             "D_ir", "D_re"})
        columns.push_back(std::string("state_") + name);
    append_tensor_columns(columns, "state_last_reversal_deviator");
    for (const char *name : {"cyclic_amplitude", "amplitude_factor",
                             "amplitude_reversals", "initial_relative_state",
                             "state_contraction_factor", "effective_knee_ratio"})
        columns.push_back(std::string("state_") + name);
    for (const char *name : {"geostatic_deviator", "static_bias_tensor",
                             "cyclic_direction"})
        append_tensor_columns(columns, std::string("state_") + name);
    for (const char *name : {"static_bias_index", "cyclic_phase_active",
                             "bias_ratchet_strain", "physical_eps_v_total",
                             "bias_reversible_volume"})
        columns.push_back(std::string("state_") + name);
    append_tensor_columns(columns,
                          "state_last_host_deviatoric_strain_direction");
    for (const char *name : {"phase_irreversible_volume",
                             "phase_reversible_volume", "phase_potential_anchor",
                             "phase_accumulation_lambda_anchor",
                             "phase_accumulation_hardening_state"})
        columns.push_back(std::string("state_") + name);
    append_tensor_columns(columns, "state_unbiased_phase_direction");
    for (const char *name : {"loose_shear_lambda_anchor",
                             "loose_shear_hardening_state",
                             "loose_shear_gate_value"})
        columns.push_back(std::string("state_") + name);
    for (const char *name : {"mapping_anchor", "mapping_backstress",
                             "mapping_directional_fabric"})
        append_tensor_columns(columns, std::string("state_") + name);
    for (const char *name : {
             "mapping_gate_value", "mapping_capacity",
             "mapping_kinematic_denominator", "mapping_shear_modulus_ratio",
             "mapping_phase_contraction_scale", "mapping_outer_residual",
             "mapping_stress_corrections", "mapping_corrector_passes",
             "mapping_monotone_caps", "initial_relative_density_value",
             "intermediate_low_gate_value", "intermediate_high_gate_base"})
        columns.push_back(std::string("state_") + name);
}

static const std::vector<std::string> &diagnostic_names()
{
    static const std::vector<std::string> names = {
        "cyclic_amplitude", "amplitude_factor", "amplitude_reversals",
        "initial_relative_state", "initial_relative_density",
        "state_contraction_factor", "effective_knee_ratio",
        "shear_modulus_factor", "cyclic_hardening_factor",
        "static_bias_index", "projected_static_bias", "bias_hardening_boost",
        "bias_reversible_volume", "bias_ratchet_activity",
        "bias_ratchet_capacity", "bias_ratchet_strain",
        "state_shakedown_factor", "state_shakedown_compliance_factor",
        "state_shakedown_shear_multiplier", "lambda_total",
        "ep_eq_since_reversal", "branch_progress",
        "branch_compliance_multiplier", "branch_bias_sign",
        "phase_dense_bias_gate", "phase_volume_density_weight",
        "phase_volume_gate", "phase_volume_replacement_gate",
        "phase_volume_dense_transition", "transformation_progress",
        "transformation_zone", "branch_stress_phase", "phase_loading_eta",
        "phase_transformation_eta", "phase_signed_state", "phase_activity",
        "phase_wave_activity", "phase_mean_activity", "phase_potential",
        "phase_irreversible_volume", "phase_reversible_volume",
        "phase_accumulation_density_gate", "phase_accumulation_bias_gate",
        "phase_accumulation_memory", "phase_accumulation_onset",
        "phase_accumulation_activity", "phase_accumulation_target_activity",
        "phase_accumulation_amplitude_gain",
        "phase_accumulation_target_hardening_state",
        "phase_accumulation_hardening_multiplier", "loose_stabilization_gate",
        "loose_phase_gate", "unbiased_phase_gate",
        "unbiased_phase_direction_norm", "loose_shear_memory",
        "loose_shear_hardening_state", "loose_shear_gate_value",
        "loose_shear_hardening_multiplier",
        "loose_shear_branch_compliance_multiplier", "mapping_gate",
        "mapping_backstress_norm", "mapping_fabric_norm", "mapping_capacity",
        "mapping_kinematic_denominator", "mapping_shear_modulus_ratio",
        "mapping_phase_contraction_scale", "mapping_outer_residual",
        "mapping_stress_corrections", "mapping_corrector_passes",
        "mapping_monotone_caps", "intermediate_bias_flow_gate",
        "intermediate_high_bias_flow_gate", "intermediate_branch_multiplier",
        "intermediate_stable_amplitude", "intermediate_phase_potential_anchor",
        "intermediate_high_bias_phase_activation",
        "intermediate_pre_reversal_cyclic_excursion",
    };
    return names;
}

static std::vector<std::string> output_header()
{
    std::vector<std::string> columns = {
        "step", "phase", "activate_before", "nsub"};
    append_tensor_columns(columns, "deps");
    append_tensor_columns(columns, "stress");
    for (const char *name : {"mean_effective_pressure", "q", "ru_vertical",
                             "accepted_substeps", "reversal_registered"})
        columns.push_back(name);
    append_state_columns(columns);
    for (const std::string &name : diagnostic_names())
        columns.push_back("diagnostic_" + name);
    if (columns.size() != 235)
        throw std::runtime_error("internal output schema is not 235 columns");
    return columns;
}

static void put_state(Row &row, const riva_ib_state_t &state)
{
    const riva_state_t &s = state.base;
    put_tensor(row, "state_stress", s.stress);
    put_tensor(row, "state_alpha", s.alpha);
    put_tensor(row, "state_alpha0", s.alpha0);
    put_tensor(row, "state_alpha01", s.alpha01);
    put_tensor(row, "state_n", s.n);
    put_tensor(row, "state_fabric", s.fabric);
    put(row, "state_D", s.D);
    put(row, "state_beta", s.beta);
    put(row, "state_lambda_total", s.lambda_total);
    put(row, "state_ep_eq_since_reversal", s.ep_eq_since_reversal);
    put(row, "state_void_ratio", s.void_ratio);
    put(row, "state_reversals", s.reversals);
    put(row, "state_pressure_floor_hits", s.pressure_floor_hits);
    put(row, "state_denominator_floor_hits", s.denominator_floor_hits);
    put(row, "state_beta_fallbacks", s.beta_fallbacks);
    put(row, "state_eps_v_total", s.eps_v_total);
    put(row, "state_eps_v_confining", s.eps_v_confining);
    put(row, "state_eps_v_irreversible", s.eps_v_irreversible);
    put(row, "state_eps_v_reversible", s.eps_v_reversible);
    put(row, "state_pressure_anchor", s.pressure_anchor);
    put(row, "state_D_ir", s.D_ir);
    put(row, "state_D_re", s.D_re);
    put_tensor(row, "state_last_reversal_deviator", s.last_reversal_deviator);
    put(row, "state_cyclic_amplitude", s.cyclic_amplitude);
    put(row, "state_amplitude_factor", s.amplitude_factor);
    put(row, "state_amplitude_reversals", s.amplitude_reversals);
    put(row, "state_initial_relative_state", s.initial_relative_state);
    put(row, "state_state_contraction_factor", s.state_contraction_factor);
    put(row, "state_effective_knee_ratio", s.effective_knee_ratio);
    put_tensor(row, "state_geostatic_deviator", s.geostatic_deviator);
    put_tensor(row, "state_static_bias_tensor", s.static_bias_tensor);
    put_tensor(row, "state_cyclic_direction", s.cyclic_direction);
    put(row, "state_static_bias_index", s.static_bias_index);
    put(row, "state_cyclic_phase_active", s.cyclic_phase_active);
    put(row, "state_bias_ratchet_strain", s.bias_ratchet_strain);
    put(row, "state_physical_eps_v_total", s.physical_eps_v_total);
    put(row, "state_bias_reversible_volume", s.bias_reversible_volume);
    put_tensor(row, "state_last_host_deviatoric_strain_direction",
               s.last_host_deviatoric_strain_direction);
    put(row, "state_phase_irreversible_volume", state.phase_irreversible_volume);
    put(row, "state_phase_reversible_volume", state.phase_reversible_volume);
    put(row, "state_phase_potential_anchor", state.phase_potential_anchor);
    put(row, "state_phase_accumulation_lambda_anchor",
        state.phase_accumulation_lambda_anchor);
    put(row, "state_phase_accumulation_hardening_state",
        state.phase_accumulation_hardening_state);
    put_tensor(row, "state_unbiased_phase_direction",
               state.unbiased_phase_direction);
    put(row, "state_loose_shear_lambda_anchor", state.loose_shear_lambda_anchor);
    put(row, "state_loose_shear_hardening_state",
        state.loose_shear_hardening_state);
    put(row, "state_loose_shear_gate_value", state.loose_shear_gate_value);
    put_tensor(row, "state_mapping_anchor", state.mapping_anchor);
    put_tensor(row, "state_mapping_backstress", state.mapping_backstress);
    put_tensor(row, "state_mapping_directional_fabric",
               state.mapping_directional_fabric);
    put(row, "state_mapping_gate_value", state.mapping_gate_value);
    put(row, "state_mapping_capacity", state.mapping_capacity);
    put(row, "state_mapping_kinematic_denominator",
        state.mapping_kinematic_denominator);
    put(row, "state_mapping_shear_modulus_ratio",
        state.mapping_shear_modulus_ratio);
    put(row, "state_mapping_phase_contraction_scale",
        state.mapping_phase_contraction_scale);
    put(row, "state_mapping_outer_residual", state.mapping_outer_residual);
    put(row, "state_mapping_stress_corrections",
        state.mapping_stress_corrections);
    put(row, "state_mapping_corrector_passes", state.mapping_corrector_passes);
    put(row, "state_mapping_monotone_caps", state.mapping_monotone_caps);
    put(row, "state_initial_relative_density_value",
        state.initial_relative_density_value);
    put(row, "state_intermediate_low_gate_value",
        state.intermediate_low_gate_value);
    put(row, "state_intermediate_high_gate_base",
        state.intermediate_high_gate_base);
}

static double base_ratchet_activity(const riva_ib_parameters_t &p,
                                    const riva_ib_state_t &s)
{
    const riva_parameters_t &b = p.base;
    if (!b.static_bias_enabled || !b.bias_ratchet_enabled ||
        !s.base.cyclic_phase_active ||
        s.base.amplitude_reversals < b.bias_minimum_reversals)
        return 0.0;
    const double bias = riva_ib_projected_bias(&s);
    if (bias <= 1.0e-14) return 0.0;
    const double amplitude_position =
        (s.base.cyclic_amplitude - b.bias_ratchet_amplitude_onset) /
        (b.bias_ratchet_amplitude_full - b.bias_ratchet_amplitude_onset);
    const double ratio_position =
        (s.base.cyclic_amplitude / bias - b.bias_ratchet_ratio_full) /
        (b.bias_ratchet_ratio_cutoff - b.bias_ratchet_ratio_full);
    return riva_smoothstep(amplitude_position) *
           (1.0 - riva_smoothstep(ratio_position));
}

static double ratchet_activity(const riva_ib_parameters_t &p,
                               const riva_ib_state_t &s)
{
    const double parent = base_ratchet_activity(p, s);
    if (!p.loose_shear_flow_enabled) return parent;
    const double gate = riva_ib_loose_gate(&p, &s);
    const double loose = riva_ib_phase_amplitude_activity(&p, &s);
    return parent + gate * (loose - parent);
}

static double ratchet_capacity(const riva_ib_parameters_t &p,
                               const riva_ib_state_t &s)
{
    const double activity = ratchet_activity(p, s);
    const double bias = riva_ib_projected_bias(&s);
    double parent = 0.0;
    if (bias > 1.0e-14)
        parent = p.base.bias_ratchet_limit * activity *
                 std::pow(p.base.bias_ratchet_reference_bias / bias,
                          p.base.bias_ratchet_bias_exponent);
    if (!p.loose_shear_flow_enabled) return parent;
    const double gate = riva_ib_loose_gate(&p, &s);
    const double loose = p.loose_shear_ratchet_capacity * activity;
    return parent + gate * (loose - parent);
}

static double shakedown_compliance(const riva_ib_parameters_t &p,
                                   const riva_ib_state_t &s, double pressure)
{
    const double raw = riva_ib_raw_shakedown(&p, &s, pressure);
    const double distance = riva_max(
        riva_ib_initial_density(&s) - p.reference_relative_density_value, 0.0);
    const double attenuation = std::exp(
        -p.base.state_shakedown_dense_compliance_decay * distance * distance);
    return 1.0 + (raw - 1.0) * attenuation;
}

static double branch_stress_phase(const riva_ib_state_t &s)
{
    const double norm = riva_norm(s.base.cyclic_direction);
    if (norm <= 1.0e-14 || s.base.cyclic_amplitude <= 1.0e-14) return 1.0;
    const tensor_t direction = riva_scale(s.base.cyclic_direction, 1.0 / norm);
    const tensor_t static_deviator = riva_add(
        s.base.geostatic_deviator,
        riva_scale(s.base.static_bias_tensor, s.base.pressure_anchor));
    const double projection = riva_ddot(
        riva_sub(riva_dev(s.base.stress), static_deviator), direction);
    return riva_clip(projection /
                         (s.base.pressure_anchor * s.base.cyclic_amplitude),
                     -1.0, 1.0);
}

static double signum(double value)
{
    return value > 0.0 ? 1.0 : (value < 0.0 ? -1.0 : 0.0);
}

static void put_diagnostics(Row &row, const riva_ib_parameters_t &p,
                            const riva_material_parameters_t &m,
                            const riva_ib_state_t &state)
{
    const riva_state_t &s = state.base;
    const double pressure = riva_pressure(s.stress);
    double shear_factor = 1.0, hardening_factor = 1.0;
    riva_ib_cyclic_flow(&p, &state, pressure, &shear_factor, &hardening_factor);
    const double projected_bias = riva_ib_projected_bias(&state);
    const double compliance = shakedown_compliance(p, state, pressure);
    double loading_eta = 0.0, eta_pt = 0.0, signed_phase = -1.0;
    riva_ib_phase_coordinates(&p, &m, &state, &loading_eta, &eta_pt,
                              &signed_phase);
    const double accumulation_bias_gate = riva_smoothstep(
        projected_bias /
        riva_max(0.5 * p.branch_compliance_bias_reference, 1.0e-14));
    const double accumulation_memory = riva_ib_phase_accumulation_memory(&state);
    const double accumulation_amplitude = riva_max(
        s.cyclic_amplitude, p.phase_accumulation_reference_amplitude);
    const double accumulation_onset = p.phase_accumulation_memory_onset *
        std::pow(p.phase_accumulation_reference_amplitude /
                     accumulation_amplitude,
                 p.phase_accumulation_amplitude_exponent);
    const double accumulation_gain = 1.0 +
        riva_smoothstep(
            (s.cyclic_amplitude - p.phase_accumulation_reference_amplitude) /
            (p.phase_accumulation_gain_full_amplitude -
             p.phase_accumulation_reference_amplitude)) *
            (p.phase_accumulation_high_amplitude_gain_ratio - 1.0);
    const double loose_transition = riva_smoothstep(
        (riva_ib_branch_progress(&p, &state) -
         p.loose_shear_branch_compliance_onset) /
        (p.loose_shear_branch_compliance_full -
         p.loose_shear_branch_compliance_onset));
    const double loose_gate = riva_ib_loose_gate(&p, &state);
    const double stable_amplitude =
        s.cyclic_amplitude *
        (s.amplitude_reversals == 1
             ? p.intermediate_first_reversal_amplitude_ratio
             : 1.0);

    const auto diagnostic = [&](const char *name, double value) {
        put(row, std::string("diagnostic_") + name, value);
    };
    diagnostic("cyclic_amplitude", s.cyclic_amplitude);
    diagnostic("amplitude_factor", s.amplitude_factor);
    diagnostic("amplitude_reversals", static_cast<double>(s.amplitude_reversals));
    diagnostic("initial_relative_state", s.initial_relative_state);
    diagnostic("initial_relative_density", riva_ib_initial_density(&state));
    diagnostic("state_contraction_factor", s.state_contraction_factor);
    diagnostic("effective_knee_ratio", s.effective_knee_ratio);
    diagnostic("shear_modulus_factor", shear_factor);
    diagnostic("cyclic_hardening_factor", hardening_factor);
    diagnostic("static_bias_index", s.static_bias_index);
    diagnostic("projected_static_bias", projected_bias);
    diagnostic("bias_hardening_boost",
               riva_ib_phase_bias_hardening(&p, &m, &state));
    diagnostic("bias_reversible_volume", s.bias_reversible_volume);
    diagnostic("bias_ratchet_activity", ratchet_activity(p, state));
    diagnostic("bias_ratchet_capacity", ratchet_capacity(p, state));
    diagnostic("bias_ratchet_strain", s.bias_ratchet_strain);
    diagnostic("state_shakedown_factor",
               riva_ib_shakedown(&p, &state, pressure));
    diagnostic("state_shakedown_compliance_factor", compliance);
    diagnostic("state_shakedown_shear_multiplier",
               riva_ib_shakedown_shear(&p, &state, pressure));
    diagnostic("lambda_total", s.lambda_total);
    diagnostic("ep_eq_since_reversal", s.ep_eq_since_reversal);
    diagnostic("branch_progress", riva_ib_branch_progress(&p, &state));
    diagnostic("branch_compliance_multiplier",
               riva_ib_branch_compliance_multiplier(&p, &m, &state));
    diagnostic("branch_bias_sign",
               signum(riva_ddot(s.static_bias_tensor, s.cyclic_direction)));
    diagnostic("phase_dense_bias_gate",
               riva_ib_phase_dense_bias_gate(&p, &state));
    diagnostic("phase_volume_density_weight",
               riva_ib_phase_volume_density_weight(&p, &state));
    diagnostic("phase_volume_gate", riva_ib_phase_volume_gate(&p, &state));
    diagnostic("phase_volume_replacement_gate",
               riva_ib_phase_volume_replacement_gate(&p, &state));
    diagnostic("phase_volume_dense_transition",
               riva_ib_phase_dense_transition(&p, &state));
    diagnostic("transformation_progress",
               riva_ib_transformation_progress(&p, &state));
    diagnostic("transformation_zone",
               riva_ib_transformation_zone(&p, &m, &state));
    diagnostic("branch_stress_phase", branch_stress_phase(state));
    diagnostic("phase_loading_eta", loading_eta);
    diagnostic("phase_transformation_eta", eta_pt);
    diagnostic("phase_signed_state", signed_phase);
    diagnostic("phase_activity",
               riva_ib_phase_activity(&p, &state, p.phase_bias_exponent,
                                      p.phase_pressure_exponent));
    diagnostic("phase_wave_activity",
               riva_ib_phase_activity(&p, &state, p.phase_wave_bias_exponent,
                                      p.phase_wave_pressure_exponent));
    diagnostic("phase_mean_activity",
               riva_ib_phase_activity(&p, &state, p.phase_mean_bias_exponent,
                                      p.phase_mean_pressure_exponent));
    diagnostic("phase_potential",
               riva_ib_signed_phase_potential(&p, &m, &state));
    diagnostic("phase_irreversible_volume", state.phase_irreversible_volume);
    diagnostic("phase_reversible_volume", state.phase_reversible_volume);
    diagnostic("phase_accumulation_density_gate",
               riva_ib_phase_accumulation_density_gate(&p, &state));
    diagnostic("phase_accumulation_bias_gate", accumulation_bias_gate);
    diagnostic("phase_accumulation_memory", accumulation_memory);
    diagnostic("phase_accumulation_onset", accumulation_onset);
    diagnostic("phase_accumulation_activity",
               riva_max(state.phase_accumulation_hardening_state, 0.0));
    diagnostic("phase_accumulation_target_activity",
               riva_ib_phase_accumulation_target_activity(&p, &state));
    diagnostic("phase_accumulation_amplitude_gain", accumulation_gain);
    diagnostic("phase_accumulation_target_hardening_state",
               riva_ib_phase_accumulation_target_hardening(&p, &state));
    diagnostic("phase_accumulation_hardening_multiplier",
               1.0 + p.phase_accumulation_hardening_gain *
                         riva_max(state.phase_accumulation_hardening_state, 0.0));
    diagnostic("loose_stabilization_gate",
               riva_ib_loose_stabilization_gate_raw(&p, &state));
    diagnostic("loose_phase_gate", loose_gate);
    diagnostic("unbiased_phase_gate", riva_ib_unbiased_gate(&p, &state));
    diagnostic("unbiased_phase_direction_norm",
               riva_norm(state.unbiased_phase_direction));
    diagnostic("loose_shear_memory", riva_ib_loose_shear_memory(&state));
    diagnostic("loose_shear_hardening_state", state.loose_shear_hardening_state);
    diagnostic("loose_shear_gate_value", state.loose_shear_gate_value);
    diagnostic("loose_shear_hardening_multiplier",
               riva_ib_loose_hardening_multiplier(&p, &state));
    diagnostic("loose_shear_branch_compliance_multiplier",
               1.0 /
                   (1.0 + loose_gate * p.loose_shear_branch_compliance_gain *
                              loose_transition));
    diagnostic("mapping_gate", riva_ib_mapping_gate(&p, &state));
    diagnostic("mapping_backstress_norm", riva_norm(state.mapping_backstress));
    diagnostic("mapping_fabric_norm",
               riva_norm(state.mapping_directional_fabric));
    diagnostic("mapping_capacity", state.mapping_capacity);
    diagnostic("mapping_kinematic_denominator",
               state.mapping_kinematic_denominator);
    diagnostic("mapping_shear_modulus_ratio", state.mapping_shear_modulus_ratio);
    diagnostic("mapping_phase_contraction_scale",
               state.mapping_phase_contraction_scale);
    diagnostic("mapping_outer_residual", state.mapping_outer_residual);
    diagnostic("mapping_stress_corrections",
               static_cast<double>(state.mapping_stress_corrections));
    diagnostic("mapping_corrector_passes",
               static_cast<double>(state.mapping_corrector_passes));
    diagnostic("mapping_monotone_caps",
               static_cast<double>(state.mapping_monotone_caps));
    diagnostic("intermediate_bias_flow_gate",
               riva_ib_intermediate_low_gate(&p, &state));
    diagnostic("intermediate_high_bias_flow_gate",
               riva_ib_intermediate_high_gate(&p, &state));
    diagnostic("intermediate_branch_multiplier",
               riva_ib_intermediate_branch_multiplier(&p, &m, &state));
    diagnostic("intermediate_stable_amplitude", stable_amplitude);
    diagnostic("intermediate_phase_potential_anchor",
               state.phase_potential_anchor);
    diagnostic("intermediate_high_bias_phase_activation",
               riva_ib_intermediate_phase_activation(&p, &state));
    diagnostic("intermediate_pre_reversal_cyclic_excursion",
               riva_ib_pre_reversal_excursion(&p, &state));
}

static Row response_row(const riva_ib_parameters_t &parameters,
                        const riva_material_parameters_t &material,
                        const riva_ib_state_t &state, int step,
                        const std::string &phase, int activate_before, int nsub,
                        tensor_t deps, const riva_update_info_t &info,
                        double vertical_effective_stress)
{
    Row row;
    put(row, "step", step);
    put_text(row, "phase", phase);
    put(row, "activate_before", activate_before);
    put(row, "nsub", nsub);
    put_tensor(row, "deps", deps);
    put_tensor(row, "stress", state.base.stress);
    put(row, "mean_effective_pressure", riva_pressure(state.base.stress));
    put(row, "q", riva_q(state.base.stress));
    put(row, "ru_vertical",
        1.0 - (-state.base.stress.zz) / vertical_effective_stress);
    put(row, "accepted_substeps", info.accepted_substeps);
    put(row, "reversal_registered", info.reversal_registered);
    put_state(row, state);
    put_diagnostics(row, parameters, material, state);
    return row;
}

static std::vector<Row> read_input_path(const fs::path &path)
{
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open " + path.string());
    std::string line;
    if (!std::getline(input, line))
        throw std::runtime_error("empty input path " + path.string());
    const std::vector<std::string> header = split_csv(line);
    std::vector<Row> rows;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != header.size())
            throw std::runtime_error("invalid CSV width in " + path.string());
        Row row;
        for (size_t i = 0; i < header.size(); ++i) row[header[i]] = values[i];
        rows.push_back(std::move(row));
    }
    if (rows.empty()) throw std::runtime_error("no increments in " + path.string());
    return rows;
}

static void write_row(std::ofstream &output,
                      const std::vector<std::string> &header, const Row &row)
{
    for (size_t i = 0; i < header.size(); ++i) {
        const auto found = row.find(header[i]);
        if (found == row.end())
            throw std::runtime_error("missing output column " + header[i]);
        output << (i ? "," : "") << found->second;
    }
    output << '\n';
}

static void replay_case(const CaseDefinition &definition,
                        const fs::path &input_directory,
                        const fs::path &output_directory,
                        const std::vector<std::string> &header)
{
    const fs::path input_path = input_directory / (std::string(definition.name) + ".csv");
    const fs::path output_path = output_directory / (std::string(definition.name) + ".csv");
    const std::vector<Row> increments = read_input_path(input_path);
    const riva_ib_parameters_t parameters = riva_ib_reference_parameters(1.0);
    const riva_material_parameters_t material =
        riva_reference_material_parameters(&parameters.base);
    const tensor_t initial_stress = {
        -kK0 * definition.vertical_effective_stress,
        -kK0 * definition.vertical_effective_stress,
        -definition.vertical_effective_stress, 0.0, 0.0, 0.0};
    riva_ib_state_t state = {};
    if (!riva_ib_initialize_material(&parameters, &material, initial_stress,
                                     definition.void_ratio, &state))
        throw std::runtime_error(std::string("initialization failed for ") +
                                 definition.name);

    std::ofstream output(output_path);
    if (!output) throw std::runtime_error("cannot create " + output_path.string());
    for (size_t i = 0; i < header.size(); ++i)
        output << (i ? "," : "") << header[i];
    output << '\n';

    const riva_update_info_t initialized_info = {};
    write_row(output, header,
              response_row(parameters, material, state, 0, "initialized", 0, 0,
                           riva_zero(), initialized_info,
                           definition.vertical_effective_stress));

    for (const Row &increment : increments) {
        const int step = std::stoi(increment.at("step"));
        const int activate_before = std::stoi(increment.at("activate_before"));
        const int nsub = std::stoi(increment.at("nsub"));
        const tensor_t deps = read_tensor(increment, "deps");
        if (activate_before &&
            !riva_ib_begin_dynamic_phase(&parameters, &material, nullptr, &state))
            throw std::runtime_error(std::string("activation failed for ") +
                                     definition.name + " step " +
                                     std::to_string(step));
        riva_update_info_t info = {};
        if (!riva_ib_update_material(&parameters, &material, deps, nsub, &state,
                                     nullptr, &info))
            throw std::runtime_error(std::string("update failed for ") +
                                     definition.name + " step " +
                                     std::to_string(step));
        write_row(output, header,
                  response_row(parameters, material, state, step,
                               increment.at("phase"), activate_before, nsub, deps,
                               info, definition.vertical_effective_stress));
    }
    std::cout << definition.name << ": " << increments.size() + 1
              << " states -> " << output_path << '\n';
}

int main(int argc, char **argv)
{
    if (argc != 3) {
        std::cerr << "usage: " << argv[0] << " INPUT_PATH_DIR OUTPUT_DIR\n";
        return 2;
    }
    try {
        const fs::path input_directory(argv[1]);
        const fs::path output_directory(argv[2]);
        if (!fs::is_directory(input_directory))
            throw std::runtime_error("input path is not a directory: " +
                                     input_directory.string());
        fs::create_directories(output_directory);
        const std::vector<std::string> header = output_header();
        for (const CaseDefinition &definition : kCases)
            replay_case(definition, input_directory, output_directory, header);
    } catch (const std::exception &error) {
        std::cerr << "RIVASandIntermediateBiasResearch replay failed: "
                  << error.what() << '\n';
        return 1;
    }
    return 0;
}
