# Acoustic-Elastic Support Notes

This is the first GPU-Hercules acoustic foundation. It is intentionally
opt-in and uses the existing mesh payload format.

## Material Marker

Set `include_acoustic = yes` in the input file. With that enabled, any material
sample with:

- `Vp > 0`
- `rho > 0`
- `Vs <= 1.0e-6`

is treated as an acoustic material. This is meant for water/ocean regions from a
1D profile, CVM, or 3D velocity model. When `include_acoustic = no`, `Vs <= 0`
is treated as an input error.

## Octree Rule

Elastic elements still use the existing rule:

```text
h <= Vs / (simulation_node_per_wavelength * simulation_wave_max_freq_hz)
```

Acoustic elements use:

```text
h <= Vp / (simulation_node_per_wavelength * simulation_wave_max_freq_hz)
```

This prevents a water layer with `Vs = 0` from forcing infinite refinement.

## Solver Constants

Acoustic elements are initialized with:

```text
mu     = 0
lambda = rho * Vp^2
```

The old shear-velocity Rayleigh damping rule is disabled for acoustic elements.
BKT-family damping coefficients are also zeroed for acoustic elements.

## Current Limitation

This first pass does not yet implement a physically exact inviscid
acoustic-elastic interface. The current mesh still shares vector displacement
DOFs across material boundaries. That is useful as a runnable foundation, but a
full water-solid coupling should duplicate/interface-couple only the normal
component or introduce pressure DOFs in acoustic elements.

The next implementation step should therefore be interface-aware meshing and
coupling: identify octants cut by the water-solid boundary, refine them to the
requested interface tolerance, and build either normal-only constraints or a
pressure-displacement coupling on the interface faces.
