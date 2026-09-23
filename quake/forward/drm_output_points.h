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

// One row: top-center X/Y, exterior X/Y lengths, exterior depth, spacing,
// clockwise strike from north (degrees). Box-local X/Y rotate about the center.
typedef std::array<double, 7> BoxSpec;
struct Face {
    const char* name;
    int rows, columns;
    size_t offset;
};
struct BoxGeometry {
    std::vector<Point> points;
    std::array<Face, 10> faces;
};

inline std::vector<BoxSpec> read_boxes(std::istream& input, const char* key, int count) {
    if (count <= 0) throw std::runtime_error("DRM box count must be positive");
    std::vector<BoxSpec> boxes;
    std::string line;
    bool found = false;
    while (std::getline(input, line)) {
        line = line.substr(0, line.find('#'));
        if (!found) {
            std::replace(line.begin(), line.end(), '=', ' ');
            std::istringstream heading(line);
            std::string name, extra;
            if (!(heading >> name) || name != key) continue;
            if (heading >> extra)
                throw std::runtime_error("Put output_drm_boxes rows on the following lines");
            found = true;
            continue;
        }
        std::istringstream row(line);
        std::string token;
        if (!(row >> token)) continue;
        if (boxes.size() == static_cast<size_t>(count)) {
            // A following parameter is allowed, but another numeric row is not.
            if (token.find_first_of("+-0123456789.") == 0)
                throw std::runtime_error("More output_drm_boxes rows than number_output_drm_boxes");
            break;
        }
        BoxSpec box;
        for (int column = 0; column < 7; ++column) {
            if (column && !(row >> token))
                throw std::runtime_error("Each output_drm_boxes row needs exactly seven values");
            size_t used = 0;
            try { box[column] = std::stod(token, &used); }
            catch (const std::exception&) { used = 0; }
            if (used != token.size() || !std::isfinite(box[column]))
                throw std::runtime_error("Invalid value in output_drm_boxes row " + std::to_string(boxes.size()));
        }
        if (row >> token)
            throw std::runtime_error("Each output_drm_boxes row needs exactly seven values");
        boxes.push_back(box);
    }
    if (input.bad() || boxes.size() != static_cast<size_t>(count))
        throw std::runtime_error("Missing or incomplete output_drm_boxes table");
    return boxes;
}

inline BoxGeometry box_faces(const BoxSpec& box) {
    for (double value : box)
        if (!std::isfinite(value)) throw std::runtime_error("DRM box values must be finite");
    const double h = box[5];
    if (h <= 0) throw std::runtime_error("DRM box spacing must be positive");
    const int nx = edge_count(box[2], h), ny = edge_count(box[3], h), nz = edge_count(box[4], h);
    if (nx < 3 || ny < 3 || nz < 2)
        throw std::runtime_error("DRM exterior lengths must exceed 2*spacing and depth must exceed spacing");
    const char* names[10] = {"exterior_front", "exterior_back", "exterior_left",
        "exterior_right", "exterior_bottom", "interior_front", "interior_back",
        "interior_left", "interior_right", "interior_bottom"};
    BoxGeometry geometry;
    size_t total = 0;
    for (int side = 0; side < 2; ++side) {
        const int x = nx - 2 * side, y = ny - 2 * side, z = nz - side;
        for (int face = 0; face < 5; ++face) {
            const int rows = face < 2 ? y + 1 : x + 1;
            const int columns = face < 4 ? z + 1 : y + 1;
            const long double size = static_cast<long double>(rows) * columns;
            if (size > max_points - total)
                throw std::runtime_error("DRM box faces exceed MPI point-count limits");
            const int index = side * 5 + face;
            geometry.faces[index] = Face{names[index], rows, columns, total};
            total += static_cast<size_t>(size);
        }
    }
    geometry.points.reserve(total);
    const double angle = std::remainder(box[6], 360.0) * std::acos(-1.0) / 180.0;
    double c = std::cos(angle), s = std::sin(angle);
    // Cardinal strikes should not put a boundary point just outside the
    // domain because sin(pi) or cos(pi/2) is a tiny nonzero value.
    if (std::abs(c) < 4 * std::numeric_limits<double>::epsilon()) c = 0;
    if (std::abs(s) < 4 * std::numeric_limits<double>::epsilon()) s = 0;
    for (int side = 0; side < 2; ++side) {
        const double hx = box[2] / 2 - side * h, hy = box[3] / 2 - side * h;
        const double depth = box[4] - side * h;
        for (int face = 0; face < 5; ++face) {
            const Face& grid = geometry.faces[side * 5 + face];
            for (int row = 0; row < grid.rows; ++row)
                for (int column = 0; column < grid.columns; ++column) {
                    double x, y, z;
                    if (face < 2) { // front=-local X, back=+local X; [Y,Z]
                        x = face == 0 ? -hx : hx;
                        y = row == grid.rows - 1 ? hy : -hy + row * h;
                        z = column == grid.columns - 1 ? depth : column * h;
                    } else if (face < 4) { // left=-local Y, right=+local Y; [X,Z]
                        x = row == grid.rows - 1 ? hx : -hx + row * h;
                        y = face == 2 ? -hy : hy;
                        z = column == grid.columns - 1 ? depth : column * h;
                    } else { // bottom; [X,Y]
                        x = row == grid.rows - 1 ? hx : -hx + row * h;
                        y = column == grid.columns - 1 ? hy : -hy + column * h;
                        z = depth;
                    }
                    geometry.points.push_back(Point{{box[0] + c*x - s*y,
                                                     box[1] + s*x + c*y, z}});
                }
        }
    }
    return geometry;
}
} // namespace drm_output
#endif
