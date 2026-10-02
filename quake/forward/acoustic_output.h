#ifndef HERCULES_ACOUSTIC_OUTPUT_H
#define HERCULES_ACOUSTIC_OUTPUT_H

/* Output-only inverse coordinate stretch. Coefficients are prepared once
 * from the containing element's solver CPML profile, including M-PML. The
 * memory belongs to this observation point, not to a centroid/GPU stencil.
 * This is the same gradient-side ADE used by kernelPMLAcousticStep:
 *   psi[n] = b*psi[n-1] + a*grad(p[n])
 *   grad_stretched(p[n]) = grad(p[n])/kappa + psi[n].
 * Outside CPML, active=0 is an exact identity. */
typedef struct acoustic_output_cpml_t {
    int active;
    double a[3], b[3], inv_kappa[3], memory[3];
} acoustic_output_cpml_t;

static inline void acoustic_output_stretch_gradient(
    acoustic_output_cpml_t* state, double grad[3]) {
    if (!state->active) return;
    for (int axis = 0; axis < 3; ++axis) {
        state->memory[axis] = state->b[axis] * state->memory[axis] +
                              state->a[axis] * grad[axis];
        grad[axis] = state->inv_kappa[axis] * grad[axis] + state->memory[axis];
    }
}

#endif
