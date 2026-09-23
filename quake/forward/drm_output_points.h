/* Geometry for the independent, single-file DRM displacement recorder.
 * Kept independent of MPI/CUDA so input and ordering can be tested locally. */
#ifndef DRM_OUTPUT_POINTS_H
#define DRM_OUTPUT_POINTS_H

#include <array>
#include <algorithm>
#include <cmath>
#include <climits>
#include <istream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace drm_output {
typedef std::array<double, 3> Point;
// Coordinates are broadcast as triples; MPI counts and local values use int.
static const size_t max_points = INT_MAX / 3;

inline std::vector<Point> read_points(std::istream& input) {
    std::vector<Point> points;
    std::string line;
    size_t line_number = 0;
    bool first_record = true;
    while (std::getline(input, line)) {
        ++line_number;
        line = line.substr(0, line.find('#'));
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream row(line);
        std::string a, b, c, extra;
        if (!(row >> a)) continue;
        if (!(row >> b >> c) || (row >> extra))
            throw std::runtime_error("Expected three columns at points line " +
                                     std::to_string(line_number));
        if (first_record && a == "latitude" && b == "longitude" && c == "depth") {
            first_record = false;
            continue;
        }
        first_record = false;
        Point point;
        const std::string tokens[3] = {a, b, c};
        for (int axis = 0; axis < 3; ++axis) {
            size_t used = 0;
            try { point[axis] = std::stod(tokens[axis], &used); }
            catch (const std::exception&) { used = 0; }
            if (used != tokens[axis].size() || !std::isfinite(point[axis]))
                throw std::runtime_error("Invalid coordinate at points line " +
                                         std::to_string(line_number));
        }
        if (std::abs(point[0]) > 90 || std::abs(point[1]) > 180 || point[2] < 0)
            throw std::runtime_error("Invalid latitude/longitude/depth at points line " +
                                     std::to_string(line_number));
        if (points.size() == max_points)
            throw std::runtime_error("Too many DRM points for MPI int counts");
        points.push_back(point);
    }
    if (input.bad()) throw std::runtime_error("Error reading DRM points");
    if (points.empty()) throw std::runtime_error("DRM points file is empty");
    return points;
}

inline int edge_count(double length, double spacing) {
    const double count = length / spacing;
    const double rounded = std::round(count);
    if (!std::isfinite(count) || rounded < 1 || rounded > INT_MAX - 3 ||
        std::abs(count - rounded) > 1e-8 +
            16 * std::numeric_limits<double>::epsilon() * std::abs(count))
        throw std::runtime_error("DRM widths and depth must be positive integer multiples of drm_edgesize");
    return static_cast<int>(rounded);
}

inline size_t surface_count(int nx, int ny, int nz) {
    // Four walls through nz-1 plus the entire bottom. No top lid.
    const long double count = (2.0L * nx + 2.0L * ny) * nz +
                              (static_cast<long double>(nx) + 1) * (ny + 1);
    if (count > max_points)
        throw std::runtime_error("DRM box exceeds MPI point-count limits");
    return static_cast<size_t>(count);
}

inline void append_surface(std::vector<Point>& points, double x0, double y0,
                           double spacing, int nx, int ny, int nz) {
    // X outer, Y middle, Z fastest; shared wall edges appear only once.
    for (int ix = 0; ix <= nx; ++ix)
        for (int iy = 0; iy <= ny; ++iy) {
            const bool wall = ix == 0 || ix == nx || iy == 0 || iy == ny;
            for (int iz = wall ? 0 : nz; iz <= nz; ++iz)
                points.push_back(Point{{x0 + ix * spacing, y0 + iy * spacing,
                                        iz * spacing}});
        }
}

inline std::vector<Point> box_points(const double boundary[5], double spacing,
                                     double offset_x, double offset_y) {
    if (!std::isfinite(spacing) || spacing <= 0 ||
        !std::isfinite(offset_x) || !std::isfinite(offset_y))
        throw std::runtime_error("DRM spacing must be positive and offsets finite");
    for (int i = 0; i < 5; ++i)
        if (!std::isfinite(boundary[i]))
            throw std::runtime_error("DRM boundary must be finite");
    const int nx = edge_count(boundary[2] - boundary[0], spacing);
    const int ny = edge_count(boundary[3] - boundary[1], spacing);
    const int nz = edge_count(boundary[4], spacing);
    const size_t inner = surface_count(nx, ny, nz);
    const size_t outer = surface_count(nx + 2, ny + 2, nz + 1);
    if (inner > max_points - outer)
        throw std::runtime_error("DRM box exceeds MPI point-count limits");
    std::vector<Point> points;
    points.reserve(inner + outer);
    const double x0 = boundary[0] + offset_x;
    const double y0 = boundary[1] + offset_y;
    append_surface(points, x0, y0, spacing, nx, ny, nz);
    append_surface(points, x0 - spacing, y0 - spacing, spacing,
                   nx + 2, ny + 2, nz + 1);
    return points;
}
} // namespace drm_output
#endif
