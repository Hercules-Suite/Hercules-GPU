/* Nominal P-wave Q implied by the selected shear/bulk attenuation branches.
 * Output metadata only: this does not fit or modify damping coefficients. */
#ifndef DRM_OUTPUT_QUALITY_H
#define DRM_OUTPUT_QUALITY_H
#include <cmath>
#include <limits>

// Keep this distinct from the legacy global drm_output() function in drm.h.
namespace drm_output_quality {
// Disabled branches use positive infinity (zero loss), not the solver's
// finite high-Q sentinel. Use the same final velocities exported with Qp.
inline double selected_qp(double vs, double vp, double qs, double qk) {
    const double unavailable = std::numeric_limits<double>::quiet_NaN();
    if (!std::isfinite(vs) || vs <= 0 || !std::isfinite(vp) || vp <= 0 ||
        std::isnan(qs) || qs <= 0 || std::isnan(qk) || qk <= 0)
        return unavailable;
    const double ratio = vs / vp;
    const double L = (4.0 / 3.0) * ratio * ratio;
    if (L >= 1.0) return unavailable;
    const double inverse_qp = L / qs + (1.0 - L) / qk;
    return inverse_qp == 0.0 ? std::numeric_limits<double>::infinity() : 1.0 / inverse_qp;
}
} // namespace drm_output_quality
#endif
