#ifndef HERCULES_ACOUSTIC_OPEN_PORE_QUADRATURE_H
#define HERCULES_ACOUSTIC_OPEN_PORE_QUADRATURE_H

#ifdef __CUDACC__
#define OPEN_PORE_HOST_DEVICE __host__ __device__
#else
#define OPEN_PORE_HOST_DEVICE
#endif

/* Four equal-weight Gauss points on a rectangular interface fragment.
 * Each side supplies its own Hex8 natural-coordinate center and half widths;
 * the tangential axes and q index identify the same physical point on both.
 * The fragment's area/compliance is already in the caller's scale, so each
 * point contributes scale/4. Share this rule between setup capacity and the
 * time-dependent pressure target to preserve constant-pressure consistency. */
OPEN_PORE_HOST_DEVICE static inline void acoustic_open_pore_shape_values(
    const double center[3], const double half_width[3],
    const int tangent_axes[2], int q, double shape[8]) {
  const double gauss = 0.57735026918962576451;
  double point[3] = {center[0], center[1], center[2]};
  for (int t = 0; t < 2; ++t) {
    const int axis = tangent_axes[t];
    point[axis] += half_width[axis] * ((q & (1 << t)) ? gauss : -gauss);
  }
  for (int node = 0; node < 8; ++node) {
    const double sx = (node & 1) ? 1.0 : -1.0;
    const double sy = (node & 2) ? 1.0 : -1.0;
    const double sz = (node & 4) ? 1.0 : -1.0;
    shape[node] = 0.125 * (1.0 + sx * point[0]) *
        (1.0 + sy * point[1]) * (1.0 + sz * point[2]);
  }
}

#undef OPEN_PORE_HOST_DEVICE
#endif
