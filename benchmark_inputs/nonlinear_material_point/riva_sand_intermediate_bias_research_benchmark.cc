/* Native research-vs-production RIVA-Sand host benchmark.
 *
 * Both kernels replay the same frozen 3484_intermediate_bias025 input rows,
 * analytical K0 initial state, activation row, and row-specific nsub values.
 * Parsing and allocation occur outside the timed region.  Python is not used.
 *
 * Usage:
 *   riva_sand_intermediate_bias_research_benchmark INPUT_PATH_DIR [REPETITIONS]
 */
#define RIVA_STANDALONE
#include "../../quake/forward/nonlinear_riva_sand_intermediate_bias_research.h"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using Row = std::unordered_map<std::string, std::string>;

struct Increment {
    tensor_t deps;
    int32_t nsub;
    bool activate_before;
};

struct TimingResult {
    double seconds;
    double checksum;
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
        increments.push_back({read_tensor(row, "deps"),
                              std::stoi(row.at("nsub")),
                              std::stoi(row.at("activate_before")) != 0});
    }
    if (increments.size() != 484)
        throw std::runtime_error("expected 484 low-bias input rows");
    return increments;
}

static const tensor_t kInitialStress = {-19.4, -19.4, -40.0, 0.0, 0.0, 0.0};
static constexpr double kVoidRatio = 0.601;

static double production_checksum(const riva_state_t &state)
{
    double values[RIVA_STATE_VALUE_COUNT] = {};
    const int count = riva_state_values(&state, values);
    if (count != RIVA_STATE_VALUE_COUNT)
        throw std::runtime_error("production flatten count changed");
    double checksum = 0.0;
    for (int i = 0; i < count; ++i)
        checksum += static_cast<double>(i + 1) * values[i];
    return checksum;
}

static double research_checksum(const riva_ib_state_t &state)
{
    double values[RIVA_IB_STATE_VALUE_COUNT] = {};
    const int count = riva_ib_state_values(&state, values);
    if (count != RIVA_IB_STATE_VALUE_COUNT)
        throw std::runtime_error("research flatten count changed");
    double checksum = 0.0;
    for (int i = 0; i < count; ++i)
        checksum += static_cast<double>(i + 1) * values[i];
    return checksum;
}

static TimingResult benchmark_production(const std::vector<Increment> &increments,
                                         int repetitions)
{
    const riva_parameters_t parameters = riva_reference_parameters(1.0);
    volatile double checksum_sink = 0.0;
    const auto start = std::chrono::steady_clock::now();
    for (int repetition = 0; repetition < repetitions; ++repetition) {
        riva_state_t state = {};
        if (!riva_initialize(&parameters, kInitialStress, kVoidRatio, &state))
            throw std::runtime_error("production initialization failed");
        for (const Increment &increment : increments) {
            if (increment.activate_before &&
                !riva_begin_dynamic_phase(&parameters, kInitialStress, &state))
                throw std::runtime_error("production activation failed");
            riva_update_info_t info = {};
            if (!riva_update(&parameters, increment.deps, increment.nsub, &state,
                             nullptr, &info))
                throw std::runtime_error("production update failed");
            const int32_t expected_substeps =
                riva_ddot(increment.deps, increment.deps) == 0.0
                    ? 0
                    : increment.nsub;
            if (info.accepted_substeps != expected_substeps)
                throw std::runtime_error("production substep contract changed");
        }
        checksum_sink = checksum_sink + production_checksum(state);
    }
    const auto stop = std::chrono::steady_clock::now();
    return {std::chrono::duration<double>(stop - start).count(), checksum_sink};
}

static TimingResult benchmark_research(const std::vector<Increment> &increments,
                                       int repetitions)
{
    const riva_ib_parameters_t parameters = riva_ib_reference_parameters(1.0,0);
    const riva_material_parameters_t material =
        riva_reference_material_parameters(&parameters.base);
    volatile double checksum_sink = 0.0;
    const auto start = std::chrono::steady_clock::now();
    for (int repetition = 0; repetition < repetitions; ++repetition) {
        riva_ib_state_t state = {};
        if (!riva_ib_initialize_material(&parameters, &material, kInitialStress,
                                         kVoidRatio, &state))
            throw std::runtime_error("research initialization failed");
        for (const Increment &increment : increments) {
            if (increment.activate_before &&
                !riva_ib_begin_dynamic_phase(&parameters, &material,
                                             &kInitialStress, &state))
                throw std::runtime_error("research activation failed");
            riva_update_info_t info = {};
            if (!riva_ib_update_material(&parameters, &material, increment.deps,
                                         increment.nsub, &state, nullptr, &info,1) ||
                info.accepted_substeps != increment.nsub)
                throw std::runtime_error("research update failed");
        }
        checksum_sink = checksum_sink + research_checksum(state);
    }
    const auto stop = std::chrono::steady_clock::now();
    return {std::chrono::duration<double>(stop - start).count(), checksum_sink};
}

int main(int argc, char **argv)
{
    if (argc < 2 || argc > 3) {
        std::cerr << "usage: " << argv[0]
                  << " INPUT_PATH_DIR [REPETITIONS]\n";
        return 2;
    }
    try {
        const int repetitions = argc == 3 ? std::stoi(argv[2]) : 2000;
        if (repetitions < 1)
            throw std::runtime_error("repetitions must be positive");
        const std::string path = std::string(argv[1]) +
            "/3484_intermediate_bias025.csv";
        const std::vector<Increment> increments = read_path(path);

        /* Warm both instruction/data paths before measuring. */
        (void)benchmark_production(increments, 1);
        (void)benchmark_research(increments, 1);

        const TimingResult production =
            benchmark_production(increments, repetitions);
        const TimingResult research = benchmark_research(increments, repetitions);
        const double total_rows =
            static_cast<double>(increments.size()) * repetitions;
        const double production_ns_per_row =
            production.seconds * 1.0e9 / total_rows;
        const double research_ns_per_row = research.seconds * 1.0e9 / total_rows;

        std::cout << std::setprecision(17)
                  << "compiler=" << __VERSION__ << '\n'
                  << "case=3484_intermediate_bias025\n"
                  << "rows_per_repetition=" << increments.size() << '\n'
                  << "repetitions=" << repetitions << '\n'
                  << "total_rows=" << static_cast<uint64_t>(total_rows) << '\n'
                  << "production_seconds=" << production.seconds << '\n'
                  << "production_ns_per_row=" << production_ns_per_row << '\n'
                  << "production_checksum=" << production.checksum << '\n'
                  << "research_seconds=" << research.seconds << '\n'
                  << "research_ns_per_row=" << research_ns_per_row << '\n'
                  << "research_checksum=" << research.checksum << '\n'
                  << "research_over_production="
                  << research_ns_per_row / production_ns_per_row << '\n'
                  << "sizeof_riva_state_t=" << sizeof(riva_state_t) << '\n'
                  << "sizeof_riva_ib_state_t=" << sizeof(riva_ib_state_t) << '\n'
                  << "sizeof_riva_parameters_t=" << sizeof(riva_parameters_t)
                  << '\n'
                  << "sizeof_riva_ib_parameters_t="
                  << sizeof(riva_ib_parameters_t) << '\n';
    } catch (const std::exception &error) {
        std::cerr << "RIVA-Sand native benchmark failed: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
