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

For non-conforming water-solid boundaries, set `maximum_octant_level` in the
input file. When acoustic support is enabled, mixed acoustic/elastic cells are
refined up to that level. At the final level, cells that still contain both
water and solid samples are kept elastic; this avoids removing shear stiffness
from solid material at the seafloor.

## Solver Constants

Acoustic elements are initialized with:

```text
mu     = 0
lambda = rho * Vp^2
```

The old shear-velocity Rayleigh damping rule is disabled for acoustic elements.
BKT-family damping coefficients are also zeroed for acoustic elements.

## Imported Meshes

When `mesh_input_mode = 1`, Hercules reuses the material payload stored in
`mesh_etree_output_file` once those properties are available. In that mode,
`3D_velocity_model_directory` is not re-sampled for every element during the
solve. For acoustic-water runs, regenerate `mesh.e` with `include_acoustic =
yes` after changing the water model or after switching to this branch.

Imported acoustic elements are sanitized after loading: `Vs` is kept at zero,
BKT-family damping coefficients are cleared, and those elements retain their
undamped pressure stiffness in the BKT/BKT2 force path.

## Current Limitation

This first pass does not yet implement a physically exact inviscid
acoustic-elastic interface. The current mesh still shares vector displacement
DOFs across material boundaries. That is useful as a runnable foundation, but a
full water-solid coupling should duplicate/interface-couple only the normal
component or introduce pressure DOFs in acoustic elements.

Station or plane output inside acoustic water should therefore not be
interpreted as physical three-component water motion; pressure output and a
proper acoustic interface remain future work.

The next implementation step should therefore be interface-aware meshing and
coupling: identify octants cut by the water-solid boundary, refine them to the
requested interface tolerance, and build either normal-only constraints or a
pressure-displacement coupling on the interface faces.
