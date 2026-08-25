# RIVA-Sand geostatic admission specification

This file is the equation-level source of truth for the Hercules-specific
transition from an equilibrated elastic geostatic stress to RIVA-Sand. The
frozen cyclic equations and calibrated constants remain specified by
`rivasand_port/model.py` and `RIVASAND_FORMULATION.md`; this extension neither
recalibrates them nor changes an ordinarily admissible initial state.

Stress is tension-positive and mean effective pressure is

```text
p_phys = -tr(sigma')/3.
```

The optional research cone translation is

```text
p_cone  = p_phys + p_residual
p_floor = p_min  + p_residual,
```

where production fixes `p_residual=0`. Admission always checks the physical
pressure and rejects `p_phys <= p_min`; a positive translated cone pressure
must never make a tensile physical state acceptable.

For the initialized void ratio, evaluate the ordinary RIVA-Sand bounding
surface at `p_cone`:

```text
xi_r    = relative_state(p_cone,e)
M_b     = M exp(n_b xi_r)
b       = sqrt(2/3) M_b
alpha   = dev(sigma')/p_cone
rho_geo = ||alpha||/b.
```

If `rho_geo <= 1`, initialization is unchanged and
`geostatic_admitted=0`. If `rho_geo > 1`, Hercules preserves the equilibrated
stress exactly, sets

```text
alpha0  = alpha
alpha01 = alpha
geostatic_admitted = 1,
```

and reports the admission. It does not project the initial stress onto the
cone.

During subsequent explicit updates, the ordinary radial bound is `b_new`.
For an admitted state the normalized trial deviator is instead limited by

```text
b_limit = max(b_new, ||alpha_old||).
```

If `||alpha_trial|| > b_limit`, it is radially scaled to `b_limit`; otherwise
it is unchanged. Consequently an over-bound normalized state can remain at its
current radius or move inward, but it cannot move farther from the current
bounding cone. This is a closed-form projection and adds no local iteration or
substep.

`geostatic_admitted` is committed constitutive state because it selects this
future regularization rule. Hercules therefore compares the flag explicitly
between host and CUDA histories and checkpoints it through the full
`riva_state_t`. The frozen 93-value numerical oracle remains unchanged; the
flag is a separate comparison field rather than a 94th oracle value.

At the first post-geostatic step, production initializes every active
RIVA-Sand integration point and reports the global initialized count, admitted
count, and maximum admitted `rho_geo`. No automatic upper cutoff is imposed because a
physically defensible cutoff is calibration- and initialization-dependent;
large ratios must be investigated rather than silently accepted as validated.
