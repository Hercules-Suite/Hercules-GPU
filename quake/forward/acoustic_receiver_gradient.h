#ifndef HERCULES_ACOUSTIC_RECEIVER_GRADIENT_H
#define HERCULES_ACOUSTIC_RECEIVER_GRADIENT_H

/* Setup-only observation algebra. A Hex8 derivative in direction axis is
 * constant along that axis. Recover its value at the receiver by interpolating
 * the derivatives at two adjacent cell centers (one-sided extrapolation at a
 * material/domain boundary). Exact for quadratic scalar fields, including
 * mixed terms, on unequal Cartesian cells. This does not change the FEM field.
 */
namespace acoustic_receiver {
inline void gradient_weights(int axis, double h, const double local[3],
                             double weights[8]) {
  for (int c = 0; c < 8; ++c) {
    double w = (c & (1 << axis)) ? 1.0 : -1.0;
    for (int d = 0; d < 3; ++d)
      if (d != axis) w *= 1.0 + ((c & (1 << d)) ? 1.0 : -1.0) * local[d];
    weights[c] = w / (4.0 * h);
  }
}

inline double neighbor_fraction(double receiver, double center,
                                double neighbor_center) {
  return (receiver - center) / (neighbor_center - center);
}
}
#endif
