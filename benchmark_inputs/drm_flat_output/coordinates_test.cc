#include "../../quake/forward/drm_output_coordinates.h"
#include <cassert>
#include <iostream>

// A non-C locale must not introduce decimal commas or digit grouping.
struct CommaNumbers : std::numpunct<char> {
    char do_decimal_point() const override { return ','; }
    char do_thousands_sep() const override { return '.'; }
    std::string do_grouping() const override { return "\3"; }
};

int main(int argc, char** argv) {
    assert(argc == 2);
    const std::string path = argv[1];
    const drm_output::BoxSpec specs[] = {
        {{1000, 2000, 4, 6, 3, 1, 37}},
        {{10, 20, 0.3, 0.5, 0.4, 0.1, 0}}};
    const std::locale original = std::locale();
    std::locale::global(std::locale(original, new CommaNumbers));
    // The second box is shorter, so reusing the path also checks truncation.
    for (const auto& spec : specs) {
        auto geometry = drm_output::box_faces(spec);
        for (auto& p : geometry.points)
            p[2] += 25.125 + 0.0125 * p[0] - 0.0025 * p[1]; // Synthetic terrain/shift.
        drm_output::write_coordinates(path, geometry.points);
        std::ifstream input(path);
        std::string line;
        size_t row = 0;
        bool columns = false, count = false;
        while (std::getline(input, line)) {
            assert(!line.empty());
            if (line[0] == '#') {
                if (line == "# X_m Y_m Z_m") columns = true;
                if (line == "# point_count = " + std::to_string(geometry.points.size())) count = true;
                continue;
            }
            std::istringstream values(line);
            values.imbue(std::locale::classic());
            drm_output::Point point;
            std::string extra;
            assert(values >> point[0] >> point[1] >> point[2]);
            assert(!(values >> extra));
            assert(row < geometry.points.size());
            assert(point == geometry.points[row++]); // Exact float64 round trip and face order.
        }
        assert(!input.bad() && columns && count && row == geometry.points.size());
    }
    std::locale::global(original);
    // A regular file cannot serve as a parent directory.
    bool failed = false;
    try { drm_output::write_coordinates(path + "/drmcoordinates.0", {}); }
    catch (const std::runtime_error&) { failed = true; }
    assert(failed);
#ifdef __linux__
    // Disk-full errors may surface only on flush/close.
    failed = false;
    try { drm_output::write_coordinates("/dev/full", drm_output::box_faces(specs[0]).points); }
    catch (const std::runtime_error&) { failed = true; }
    assert(failed);
#endif
    std::cout << "DRM ASCII coordinates: exact order, duplicates, rotation, shifted depths, precision, locale and truncation passed\n";
}
