#ifndef HERCULES_CPML_H
#define HERCULES_CPML_H

#include <stddef.h>

/* Second-order CFS-PML, Xie et al. (2014), eqs. (8)-(10).
 * S_j(s)=kappa_j+d_j/(s+alpha_j), J=S_x*S_y*S_z.
 * Inertia: J(s)*(s*s+a*s). Stiffness: C_ijkl*J/(S_j*S_l).
 * Products are realized as cascades, NOT residues divided by differences
 * between poles. Equal/zero alpha and edge/corner products are well defined.
 * Each first-order section uses the trapezoidal (bilinear) recurrence.
 * Coefficients are prepared on the host; history stays on the device.
 */
#ifdef __CUDACC__
#define CPML_HD __host__ __device__
#else
#define CPML_HD
#endif

typedef struct {
    double direct, decay, feed;
} cpml_filter_t;

typedef struct {
    double mass, damping, stiffness;
    double residue[3];
    cpml_filter_t pole[3];
    int order;
} cpml_inertia_t;

typedef struct {
    cpml_filter_t forward[3], inverse[3];
    cpml_inertia_t inertia;
    double nodal_mass_dt2;
    double mu, lambda; /* dt^2*h times moduli; acoustic mu=dt^2*h/rho */
    double b_over_dt;
} cpml_element_t;

typedef struct {
    cpml_inertia_t inertia;
    double nodal_mass_dt2;
} cpml_surface_t;

typedef struct {
    cpml_filter_t forward[3];
} cpml_interface_t;

#ifdef __cplusplus

/* c+g/(s+p), transposed direct form with one persistent scalar.
 * y[n]=direct*x[n]+memory; memory'=decay*memory+feed*x[n]. */
static inline cpml_filter_t cpml_filter(double c, double g, double p,
                                         double dt) {
    cpml_filter_t f;
    const double beta = dt / (2.0 + p*dt);
    f.decay = (2.0-p*dt)/(2.0+p*dt);
    f.direct = c + g*beta;
    f.feed = g*beta*(1.0+f.decay);
    return f;
}

/* All coefficient/state arrays are small fixed-size objects. No allocation
 * or transcendental function is used in the device routines below. */
template <typename T>
CPML_HD static inline double cpml_apply(const cpml_filter_t& f, double x,
                                        T* state) {
    const double old = (double)*state;
    const double y = f.direct*x + old;
    *state = (T)(f.decay*old + f.feed*x);
    return y;
}

static inline void cpml_poly_factor(double p[6], int degree, double root) {
    for (int k=degree+1; k>0; --k) p[k] = p[k-1] + root*p[k];
    p[0] *= root;
}

/* Divide J(s)*(s^2+a*s) into M*s^2+C*s+K and a proper rational
 * remainder. Expand the latter in the nested basis
 * 1/(s+p0), 1/((s+p0)(s+p1)), ... . Zero-damping factors cancel BEFORE
 * division, so inactive axes introduce neither poles nor roundoff drift. */
static inline cpml_inertia_t cpml_inertia(const double d[3],
    const double alpha[3], const double kappa[3], int axes,
    double material_a, double dt) {
    cpml_inertia_t out = {};
    double numerator[6] = {1,0,0,0,0,0};
    double denominator[6] = {1,0,0,0,0,0};
    double poles[3] = {0,0,0};
    double scale = 1.0;
    int order = 0;
    for (int j=0; j<axes; ++j) {
        scale *= kappa[j];
        if (d[j] == 0.0) continue;
        poles[order] = alpha[j];
        cpml_poly_factor(numerator, order, alpha[j]+d[j]/kappa[j]);
        cpml_poly_factor(denominator, order, alpha[j]);
        ++order;
    }
    double remainder[6] = {0,0,0,0,0,0};
    for (int k=0; k<=order; ++k) {
        remainder[k+2] += scale*numerator[k];
        remainder[k+1] += scale*material_a*numerator[k];
    }
    double quotient[3] = {0,0,0};
    for (int k=2; k>=0; --k) {
        quotient[k] = remainder[order+k];
        for (int j=0; j<=order; ++j)
            remainder[j+k] -= quotient[k]*denominator[j];
    }
    out.mass = quotient[2];
    out.damping = quotient[1];
    out.stiffness = quotient[0];
    out.order = order;
    for (int j=0; j<order; ++j) {
        const int degree = order-j-1;
        const double r = remainder[degree];
        out.residue[j] = r;
        double tail[6] = {1,0,0,0,0,0};
        for (int k=j+1; k<order; ++k)
            cpml_poly_factor(tail, k-j-1, poles[k]);
        for (int k=0; k<=degree; ++k) remainder[k] -= r*tail[k];
        out.pole[j] = cpml_filter(0.0, 1.0, poles[j], dt);
    }
    return out;
}

static inline cpml_element_t cpml_element(const double d[3],
    const double alpha[3], const double kappa[3], double dt,
    double material_a, double b_over_dt, double mass, double mu,
    double lambda) {
    cpml_element_t c = {};
    for (int j=0; j<3; ++j) {
        c.forward[j] = cpml_filter(kappa[j], d[j], alpha[j], dt);
        c.inverse[j] = cpml_filter(1.0/kappa[j], -d[j]/(kappa[j]*kappa[j]),
                                  alpha[j]+d[j]/kappa[j], dt);
    }
    c.inertia = cpml_inertia(d, alpha, kappa, 3, material_a, dt);
    c.nodal_mass_dt2 = dt*dt*mass;
    c.mu = mu; c.lambda = lambda; c.b_over_dt = b_over_dt;
    return c;
}

static inline cpml_interface_t cpml_interface(const double d[3],
    const double alpha[3], const double kappa[3], double dt) {
    cpml_interface_t c;
    for (int j=0; j<3; ++j)
        c.forward[j] = cpml_filter(kappa[j], d[j], alpha[j], dt);
    return c;
}

/* The cofactor J/S_j multiplies normal component j in each interface weak
 * form. Retain all components for sloping normals and all four quadrature
 * points; a PML interface must not silently revert to centroid integration. */
template <typename T>
CPML_HD static inline void cpml_interface_apply(const cpml_interface_t& c,
    double v[3], T* history, size_t stride) {
    for (int j=0; j<3; ++j) {
        v[j] = cpml_apply(c.forward[(j+1)%3], v[j], history+2*j*stride);
        v[j] = cpml_apply(c.forward[(j+2)%3], v[j], history+(2*j+1)*stride);
    }
}

/* State channels are SoA: state[channel*count+element]. Native integration
 * uses every Hex8 trilinear mode; reduced integration keeps only the linear
 * modes. The latter is rank deficient and is a diagnostic option only. */
CPML_HD static inline int cpml_stiffness_states(int components, int native) {
    const int moments = native ? 4 : 1;
    return 9*components*moments + (components == 3 ? (native ? 24 : 12) : 0);
}
CPML_HD static inline int cpml_element_states(int components, int native) {
    return cpml_stiffness_states(components, native) + 24*components;
}
CPML_HD static inline double cpml_sign(int node, int mode) {
    int bits = ((~node)&mode);
    bits = (bits&1) + ((bits>>1)&1) + ((bits>>2)&1);
    return (bits&1) ? -1.0 : 1.0;
}

template <typename T>
CPML_HD static inline double cpml_diagonal(const cpml_element_t& c, int axis,
    double x, T* history, size_t stride) {
    /* J/S_axis^2 = S_other1*S_other2/S_axis. */
    x = cpml_apply(c.forward[(axis+1)%3], x, history);
    x = cpml_apply(c.forward[(axis+2)%3], x, history+stride);
    return cpml_apply(c.inverse[axis], x, history+2*stride);
}

template <typename T>
CPML_HD static inline double cpml_inertia_force(const cpml_inertia_t& c,
    double u, T* history, size_t stride) {
    double value = c.stiffness*u;
    for (int j=0; j<c.order; ++j) {
        u = cpml_apply(c.pole[j], u, history+j*stride);
        value += c.residue[j]*u;
    }
    return -value;
}

/* Exact spatial integration for element-constant profiles on a cube,
 * expressed in orthogonal sign modes instead of eight Gauss-point histories.
 * Input/output arrays have a fixed 3-component stride; acoustic uses comp 0.
 * Returns dt^2 times stiffness force; inertia is assembled separately using
 * the same nodal lumping as the leading mass, never a centroid average. */
template <typename T>
CPML_HD static inline void cpml_stiffness(const cpml_element_t& c,
    int components, int native, const double u[8][3], double force[8][3],
    T* history, size_t stride) {
    double modal[8][3] = {};
    double fm[8][3] = {};
    for (int m=1; m<8; ++m)
        for (int n=0; n<8; ++n)
            for (int i=0; i<components; ++i)
                modal[m][i] += 0.125*cpml_sign(n,m)*u[n][i];
    const int moments = native ? 4 : 1;
    for (int j=0; j<3; ++j) {
        const int k=(j+1)%3, l=(j+2)%3;
        for (int t=0; t<moments; ++t) {
            const int m=(1<<j) | ((t&1) ? (1<<k) : 0) |
                                  ((t&2) ? (1<<l) : 0);
            const double weight=0.5/((t&1) ? 3.0 : 1.0)/((t&2) ? 3.0 : 1.0);
            for (int i=0; i<components; ++i) {
                const int channel=3*((j*components+i)*moments+t);
                const double v=cpml_diagonal(c,j,modal[m][i],history+channel*stride,stride);
                const double modulus=(components==3 && i==j) ? c.lambda+2*c.mu : c.mu;
                fm[m][i] += weight*modulus*v;
            }
        }
    }
    if (components==3) {
        const int offset=27*moments, cross_moments=native ? 2 : 1;
        for (int pair=0; pair<3; ++pair) {
            const int j=pair, k=(pair+1)%3, l=(pair+2)%3;
            for (int t=0; t<cross_moments; ++t) {
                const int common=t ? (1<<l) : 0;
                const double weight=t ? (1.0/6.0) : 0.5;
                for (int dir=0; dir<2; ++dir) {
                    const int a=dir ? k : j, b=dir ? j : k;
                    const int in=(1<<b)|common, out=(1<<a)|common;
                    const int channel=offset+2*((pair*cross_moments+t)*2+dir);
                    const double normal=cpml_apply(c.forward[l],modal[in][b],history+channel*stride);
                    const double shear=cpml_apply(c.forward[l],modal[in][a],history+(channel+1)*stride);
                    fm[out][a] += weight*c.lambda*normal;
                    fm[out][b] += weight*c.mu*shear;
                }
            }
        }
    }
    for (int n=0; n<8; ++n)
        for (int i=0; i<components; ++i) {
            force[n][i]=0.0;
            for (int m=1; m<8; ++m)
                force[n][i] -= cpml_sign(n,m)*fm[m][i];
        }
}

#endif /* __cplusplus */

#undef CPML_HD
#endif
