/* ASCII coordinates in the same sample order as a DRM displacement snapshot.
 * Independent of MPI/CUDA; the caller supplies the final solver coordinates. */
#ifndef DRM_OUTPUT_COORDINATES_H
#define DRM_OUTPUT_COORDINATES_H

#include "drm_output_points.h"
#include <fstream>
#include <iomanip>
#include <locale>

namespace drm_output {
inline void write_coordinates(const std::string& path, const std::vector<Point>& points) {
    std::ofstream output(path.c_str(), std::ios::out | std::ios::trunc);
    if (!output) throw std::runtime_error("Cannot open DRM coordinate file: " + path);
    output.imbue(std::locale::classic());
    output << "# Hercules DRM sample coordinates\n"
           << "# Solver X north, Y east, Z positive down; metres.\n"
           << "# Includes rotation and building/topography depth shifts.\n"
           << "# Data row i matches binary sample i (zero-based); shared face edges/corners repeat.\n"
           << "# point_count = " << points.size() << "\n"
           << "# X_m Y_m Z_m\n"
           << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (const Point& point : points)
        output << point[0] << ' ' << point[1] << ' ' << point[2] << '\n';
    output.close();
    if (!output) throw std::runtime_error("Cannot write DRM coordinate file: " + path);
}
} // namespace drm_output
#endif
