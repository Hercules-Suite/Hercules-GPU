/* Native host-side RIVA-Sand golden-data driver.
 *
 * This executable reads the prescribed deps_* columns and row-zero initial
 * state from one oracle CSV, then writes a same-schema native result. It does
 * not call Python or consume any expected constitutive response after row 0.
 */
#define RIVA_STANDALONE
#include "../../quake/forward/nonlinear_riva_sand.h"

#include <cerrno>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using Row = std::unordered_map<std::string, std::string>;

static std::vector<std::string> split_csv(const std::string &line)
{
    std::vector<std::string> values;
    std::stringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ',')) values.push_back(value);
    if (!line.empty() && line.back() == ',') values.emplace_back();
    return values;
}

static tensor_t read_tensor(const Row &row, const std::string &prefix)
{
    auto value = [&](const char *component) {
        return std::stod(row.at(prefix + "_" + component));
    };
    return {value("xx"), value("yy"), value("zz"), value("xy"),
            value("yz"), value("xz")};
}

static void put_tensor(Row &row, const std::string &prefix, tensor_t value)
{
    auto put = [&](const char *component, double scalar) {
        std::ostringstream text;
        text << std::setprecision(17) << scalar;
        row[prefix + "_" + component] = text.str();
    };
    put("xx", value.xx); put("yy", value.yy); put("zz", value.zz);
    put("xy", value.xy); put("yz", value.yz); put("xz", value.xz);
}

template <class Value>
static void put(Row &row, const std::string &name, Value value)
{
    std::ostringstream text;
    text << std::setprecision(17) << value;
    row[name] = text.str();
}

static void put_state(Row &row, const riva_state_t &state)
{
    put_tensor(row, "stress", state.stress);
    put_tensor(row, "alpha", state.alpha);
    put_tensor(row, "alpha0", state.alpha0);
    put_tensor(row, "alpha01", state.alpha01);
    put_tensor(row, "n", state.n);
    put_tensor(row, "fabric", state.fabric);
    put(row, "D", state.D); put(row, "beta", state.beta);
    put(row, "lambda_total", state.lambda_total);
    put(row, "ep_eq_since_reversal", state.ep_eq_since_reversal);
    put(row, "void_ratio", state.void_ratio);
    put(row, "reversals", state.reversals);
    put(row, "pressure_floor_hits", state.pressure_floor_hits);
    put(row, "denominator_floor_hits", state.denominator_floor_hits);
    put(row, "beta_fallbacks", state.beta_fallbacks);
    put(row, "eps_v_total", state.eps_v_total);
    put(row, "eps_v_confining", state.eps_v_confining);
    put(row, "eps_v_irreversible", state.eps_v_irreversible);
    put(row, "eps_v_reversible", state.eps_v_reversible);
    put(row, "pressure_anchor", state.pressure_anchor);
    put(row, "D_ir", state.D_ir); put(row, "D_re", state.D_re);
    put_tensor(row, "last_reversal_deviator", state.last_reversal_deviator);
    put(row, "cyclic_amplitude", state.cyclic_amplitude);
    put(row, "amplitude_factor", state.amplitude_factor);
    put(row, "amplitude_reversals", state.amplitude_reversals);
    put(row, "initial_relative_state", state.initial_relative_state);
    put(row, "state_contraction_factor", state.state_contraction_factor);
    put(row, "effective_knee_ratio", state.effective_knee_ratio);
    put_tensor(row, "geostatic_deviator", state.geostatic_deviator);
    put_tensor(row, "static_bias_tensor", state.static_bias_tensor);
    put_tensor(row, "cyclic_direction", state.cyclic_direction);
    put(row, "static_bias_index", state.static_bias_index);
    put(row, "cyclic_phase_active", state.cyclic_phase_active);
    put(row, "bias_ratchet_strain", state.bias_ratchet_strain);
    put(row, "physical_eps_v_total", state.physical_eps_v_total);
    put(row, "bias_reversible_volume", state.bias_reversible_volume);
    put_tensor(row, "last_host_deviatoric_strain_direction",
               state.last_host_deviatoric_strain_direction);
}

int main(int argc, char **argv)
{
    if (argc != 3) {
        std::cerr << "usage: " << argv[0] << " golden-input.csv native-output.csv\n";
        return 2;
    }
    std::ifstream input(argv[1]);
    if (!input) {
        std::cerr << "cannot open " << argv[1] << '\n';
        return 2;
    }
    std::string line;
    if (!std::getline(input, line)) return 2;
    const std::vector<std::string> header = split_csv(line);
    std::vector<Row> rows;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != header.size()) {
            std::cerr << "invalid CSV row width in " << argv[1] << '\n';
            return 2;
        }
        Row row;
        for (size_t i = 0; i < header.size(); ++i) row[header[i]] = values[i];
        rows.push_back(std::move(row));
    }
    if (rows.empty()) return 2;

    const riva_parameters_t parameters = riva_reference_parameters(1.0);
    riva_state_t state = {};
    const tensor_t initial_stress = read_tensor(rows.front(), "stress");
    tensor_t reference_stress = initial_stress;
    reference_stress.xy = reference_stress.yz = reference_stress.xz = 0.0;
    if (!riva_initialize(&parameters, initial_stress,
                        std::stod(rows.front().at("void_ratio")), &state) ||
        !riva_begin_dynamic_phase(&parameters, reference_stress, &state)) {
        std::cerr << "RIVA-Sand initialization failed\n";
        return 1;
    }
    const double p0 = riva_pressure(initial_stress);
    const double sv0 = -initial_stress.zz;

    for (size_t index = 0; index < rows.size(); ++index) {
        riva_update_info_t info = {};
        if (index > 0) {
            const tensor_t deps = read_tensor(rows[index], "deps");
            if (!riva_update(&parameters, deps,
                            std::stoi(rows[index].at("nsub")), &state,
                            nullptr, &info)) {
                std::cerr << "RIVA-Sand update failed at row " << index + 2
                          << '\n';
                return 1;
            }
        }
        put(rows[index], "accepted_substeps", info.accepted_substeps);
        put(rows[index], "reversal_registered", info.reversal_registered);
        put(rows[index], "mean_effective_pressure", riva_pressure(state.stress));
        put(rows[index], "q", riva_q(state.stress));
        put(rows[index], "ru_vertical", 1.0 - (-state.stress.zz) / sv0);
        put(rows[index], "ru_mean", 1.0 - riva_pressure(state.stress) / p0);
        put(rows[index], "compatibility_residual",
            riva_compatibility_residual(&state));
        put_state(rows[index], state);
    }

    std::ofstream output(argv[2]);
    if (!output) {
        std::cerr << "cannot create " << argv[2] << '\n';
        return 2;
    }
    for (size_t i = 0; i < header.size(); ++i)
        output << (i ? "," : "") << header[i];
    output << '\n';
    for (const Row &row : rows) {
        for (size_t i = 0; i < header.size(); ++i)
            output << (i ? "," : "") << row.at(header[i]);
        output << '\n';
    }
    std::cout << "RIVA-Sand native history: " << rows.size()
              << " states -> " << argv[2] << '\n';
    return 0;
}
