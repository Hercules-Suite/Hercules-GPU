/* Full low-frequency Biot equations, using w = porosity * (U_f - u).
 * SI units. Material rho is the saturated mixture density, not grain density.
 * These small algebraic operations are shared by initialization and CUDA.
 */
#ifndef HERCULES_POROELASTIC_BIOT_H
#define HERCULES_POROELASTIC_BIOT_H

#include <math.h>

#if defined(__CUDACC__)
#define BIOT_HD __host__ __device__
#else
#define BIOT_HD
#endif

/* Native VT geometry: five P1 tetrahedra, with the same retained volume
 * fractions and partition as the solid operator. No remapping to Hex8. */
struct biot_cut_t { double fraction[5]; int partition; };
struct biot_boundary_t {
  int node, modes;
  double basis[3][3], drag[3]; /* orthonormal allowed modes, diagonal Darcy drag */
};

BIOT_HD static inline double biot_tet_geometry(int partition, int t,
    int nodes[4], double grad[4][3]) {
  const int table[2][5][4] = {
    {{0,2,1,4},{3,1,2,7},{6,4,7,2},{5,7,4,1},{2,4,7,1}},
    {{0,3,1,5},{0,2,3,6},{4,5,6,0},{6,5,7,3},{0,6,3,5}}};
  double J[3][3];
  for (int a=0; a<4; ++a) nodes[a]=table[partition-1][t][a];
  for (int d=0; d<3; ++d) for (int a=0; a<3; ++a)
    J[d][a]=((nodes[a+1]>>d)&1)-((nodes[0]>>d)&1);
  const double det=J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])
      -J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])
      +J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);
  for (int a=0; a<3; ++a) for (int d=0; d<3; ++d)
    grad[a+1][d]=(J[(d+1)%3][(a+1)%3]*J[(d+2)%3][(a+2)%3]
                    -J[(d+1)%3][(a+2)%3]*J[(d+2)%3][(a+1)%3])/det;
  for (int d=0; d<3; ++d) grad[0][d]=-grad[1][d]-grad[2][d]-grad[3][d];
  return fabs(det)/6;
}

/* Unit-cube point, P1 weights embedded in the eight cube vertices. */
BIOT_HD static inline int biot_tet_weights(int partition, const double x[3],
    double weights[8]) {
  for (int j=0; j<8; ++j) weights[j]=0;
  for (int t=0; t<5; ++t) {
    int nodes[4]; double grad[4][3], N[4]={1,0,0,0};
    biot_tet_geometry(partition,t,nodes,grad);
    for (int a=1; a<4; ++a) {
      for (int d=0; d<3; ++d) N[a]+=grad[a][d]*(x[d]-((nodes[0]>>d)&1));
      N[0]-=N[a];
    }
    if (N[0]>=-1e-10 && N[1]>=-1e-10 && N[2]>=-1e-10 && N[3]>=-1e-10) {
      for (int a=0; a<4; ++a) weights[nodes[a]]=N[a];
      return t;
    }
  }
  return -1;
}

/* Intersect a terrain tangent plane with impermeable coordinate axes, then
 * diagonalize the projected anisotropic drag. Setup only; at most two modes. */
static inline void biot_boundary_modes(const double normal[3],
    const double drag[3], biot_boundary_t* out) {
  double constraints[3][3]={}, basis[2][3]={}; int nc=0, nb=0;
  for (int d=0; d<3; ++d) if (drag[d]<0) constraints[nc++][d]=1;
  double n[3]; double norm=0;
  for (int d=0; d<3; ++d) { n[d]=drag[d]<0?0:normal[d]; norm+=n[d]*n[d]; }
  if (norm>1e-24) {
    for (int d=0; d<3; ++d) constraints[nc][d]=n[d]/sqrt(norm);
    ++nc;
  }
  for (int axis=0; axis<3 && nb<3-nc; ++axis) {
    double v[3]={}; v[axis]=1;
    for (int a=0; a<nc; ++a) {
      double dot=0; for (int d=0; d<3; ++d) dot+=v[d]*constraints[a][d];
      for (int d=0; d<3; ++d) v[d]-=dot*constraints[a][d];
    }
    for (int a=0; a<nb; ++a) {
      double dot=0; for (int d=0; d<3; ++d) dot+=v[d]*basis[a][d];
      for (int d=0; d<3; ++d) v[d]-=dot*basis[a][d];
    }
    double length=0; for (int d=0; d<3; ++d) length+=v[d]*v[d];
    if (length>1e-24) {
      for (int d=0; d<3; ++d) basis[nb][d]=v[d]/sqrt(length);
      ++nb;
    }
  }
  if (nb==2) {
    double a=0,b=0,c=0;
    for (int d=0; d<3; ++d) {
      const double D=fmax(0.0,drag[d]);
      a+=D*basis[0][d]*basis[0][d]; b+=D*basis[0][d]*basis[1][d];
      c+=D*basis[1][d]*basis[1][d];
    }
    const double theta=0.5*atan2(2*b,a-c), co=cos(theta), si=sin(theta);
    for (int d=0; d<3; ++d) {
      const double u=basis[0][d],v=basis[1][d];
      basis[0][d]=co*u+si*v; basis[1][d]=-si*u+co*v;
    }
  }
  out->modes=nb;
  for (int a=0; a<nb; ++a) {
    out->drag[a]=0;
    for (int d=0; d<3; ++d) {
      out->basis[a][d]=basis[a][d];
      out->drag[a]+=fmax(0.0,drag[d])*basis[a][d]*basis[a][d];
    }
  }
}

/* Eigenvalues of H x = c^2 R x. The product formula avoids cancellation
 * in the slow-wave eigenvalue when the two speeds are very different. */
BIOT_HD static inline void biot_compressional_speeds_squared(
    double rho, double rho_f, double inertia, double drained_p_modulus,
    double alpha, double storage, double* fast, double* slow) {
  const double det_r = rho * inertia - rho_f * rho_f;
  const double h00 = drained_p_modulus + alpha * alpha * storage;
  const double h01 = alpha * storage;
  const double trace = (inertia * h00 + rho * storage - 2 * rho_f * h01) / det_r;
  const double product = drained_p_modulus * storage / det_r;
  const double disc = fmax(0.0, trace * trace - 4 * product);
  *fast = 0.5 * (trace + sqrt(disc));
  *slow = product / *fast;
}

/* Advance staggered velocities, represented as displacement increments.
 * Fu and Fw are dt^2 times the assembled forces at integer time n.
 * Integrate constant-force Darcy relaxation exactly over one velocity step;
 * r=0 is the inviscid limit. This removes the drag relaxation CFL restriction
 * without removing the wave CFL restriction. ms*mw-mf^2 must be positive.
 */
BIOT_HD static inline void biot_increment(
    double ms, double mf, double mw, double drag_dt,
    double previous_du, double previous_dw, double Fu, double Fw,
    double* du, double* dw) {
  const double ratio = mf / ms;
  const double effective_mass = mw - mf * ratio;
  const double r = drag_dt / effective_mass;
  const double decay = exp(-r);
  const double phi = r == 0.0 ? 1.0 : -expm1(-r) / r;
  *dw = decay * previous_dw + phi * (Fw - ratio * Fu) / effective_mass;
  *du = previous_du + Fu / ms - ratio * (*dw - previous_dw);
}

#undef BIOT_HD
#endif
