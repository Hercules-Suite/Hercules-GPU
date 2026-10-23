/* Setup helpers for the scalar acoustic first-order radiation boundary.
 * q is pressure or density-scaled velocity potential:
 *   d_n q + q_t / c = 0,  C_ij = integral_face N_i N_j / (rho*c).
 * Row-sum lumping gives h^2/(4*rho*c) at each of a Hex8 face's four nodes.
 * Node bits are X=1,Y=2,Z=4; face bits are -X,+X,-Y,+Y,-Z,+Z.
 * The top (-Z) is always excluded: its existing free-surface rule wins.
 */
#ifndef HERCULES_ACOUSTIC_DASHPOT_H
#define HERCULES_ACOUSTIC_DASHPOT_H

#include <stdint.h>

template <typename Tick>
static inline uint8_t acoustic_dashpot_boundary_faces(
    const Tick lower[3], const Tick upper[3],
    const Tick domain_lower[3], const Tick domain_upper[3]) {
  uint8_t faces = 0;
  for (int axis = 0; axis < 3; ++axis) {
    if (axis < 2 && lower[axis] == domain_lower[axis])
      faces |= (uint8_t)(1u << (2 * axis));
    if (upper[axis] == domain_upper[axis])
      faces |= (uint8_t)(1u << (2 * axis + 1));
  }
  return faces;
}

/* Return dt*C per node. Contributions from meeting outer faces add;
 * internal material interfaces and MPI partition boundaries never enter
 * the geometric domain-face mask. No allocation or runtime state. */
static inline void acoustic_dashpot_lumped_dt_c(
    uint8_t faces, double size, double density, double speed, double dt,
    double nodal_dt_c[8]) {
  const double face_dt_c = dt * size * size / (4.0 * density * speed);
  faces &= (uint8_t)~(1u << 4); /* Never damp the free-surface normal. */
  for (int node = 0; node < 8; ++node) {
    unsigned count = 0;
    for (int axis = 0; axis < 3; ++axis) {
      const unsigned side = (node >> axis) & 1;
      count += (faces >> (2 * axis + side)) & 1;
    }
    nodal_dt_c[node] = face_dt_c * count;
  }
}

#endif
