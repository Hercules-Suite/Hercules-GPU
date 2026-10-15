/* ASCII coordinates/materials in the sample order of a DRM displacement snapshot.
 * Independent of MPI/CUDA; the caller supplies final coordinates and properties. */
#ifndef DRM_OUTPUT_COORDINATES_H
#define DRM_OUTPUT_COORDINATES_H

#include "drm_output_points.h"
#include <fstream>
#include <iomanip>
#include <locale>

namespace drm_output {
// Initial containing-element properties: Vs, Vp, density, Qs, Qp.
typedef std::array<double, 5> Material;

inline void write_coordinates(const std::string& path, const std::vector<Point>& points,
                              const std::vector<Material>& materials) {
    if (materials.size() != points.size())
        throw std::runtime_error("DRM coordinate/material row counts differ");
    std::ofstream output(path.c_str(), std::ios::out | std::ios::trunc);
    if (!output) throw std::runtime_error("Cannot open DRM coordinate file: " + path);
    output.imbue(std::locale::classic());
    output << "# Hercules DRM sample coordinates\n"
           << "# Solver-local X/Y, Z positive down; metres; see coordinate_frame.txt for PROJ mapping.\n"
           << "# Includes rotation and building/topography depth shifts.\n"
           << "# Initial containing-element Vs/Vp (m/s), density (kg/m^3), dimensionless Qs/Qp; no spatial averaging across elements.\n"
           << "# Qs/Qp: selected nominal BKT attenuation after branch overrides/cutoffs; Qp uses final Vs/Vp, not input Qp.\n"
           << "# inf means zero selected material loss; nan means no applicable scalar Q (e.g. Rayleigh/mass).\n"
           << "# Q values do not describe effective Rayleigh, plastic/hysteretic or poroelastic attenuation.\n"
           << "# Data row i matches binary sample i (zero-based); shared face edges/corners repeat.\n"
           << "# point_count = " << points.size() << "\n"
           << "# X_m Y_m Z_m Vs_m_s Vp_m_s density_kg_m3 Qs Qp\n"
           << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (size_t i = 0; i < points.size(); ++i) {
        const Point& point = points[i];
        output << point[0] << ' ' << point[1] << ' ' << point[2];
        for (double value : materials[i]) output << ' ' << value;
        output << '\n';
    }
    output.close();
    if (!output) throw std::runtime_error("Cannot write DRM coordinate file: " + path);
}
} // namespace drm_output
#endif
