# RIVA-Sand G0 input migration and host checks — 2026-09-08

All four names (`RIVASand`, `RIVASandCustom`, `RIVASand02`, and
`RIVASandCustom02`) now use the same complete material-table format:

```text
Vs_limit Dr G0 M kd h m zeta e_max e_min Q R n_G
```

The standard names no longer accept Dr-only rows. Custom names remain
accepted as full-input aliases within their respective generation. Model
identifiers and the distinction between original and 02 equations remain.
The material-point driver takes the same twelve material values without
Vs_limit. Vs_limit selects/interpolates rows; it does not derive stiffness.

G0 is positive and dimensionless. Reference `483.48301127222084` retains
the original elastic scale exactly. Other values scale elastic shear and bulk
stiffness together at fixed Poisson ratio. The pressure exponent and existing
pressure guards are unchanged. The conversion is cached once per material;
it adds no per-substep arithmetic and adds one double to the compact material
record. All names now allocate one compact material record per nonlinear
element, rather than sharing a single record for the standard names.

## Restart migration

The current native checkpoint markers are 3 for original RIVA-Sand and 5 for
both 02 names. These adapter/checkpoint markers are distinct from the unchanged
equation revisions. Pre-G0 checkpoints are rejected: rerun initialization.
Retain the identical complete material inputs when restarting a new checkpoint,
since material properties are reconstructed from the input rather than stored
in the constitutive history. State layouts and reference fingerprints are unchanged.

## Completed local checks

- All five original golden histories and all six revision-3 successor
  histories passed native replay. Reference files were not changed, and
  their existing tolerances were retained.
- Successor replay compared all 236 output columns, including the full state.
- Native state, geostatic admission, synthetic BrickUP coupling, plastic
  activity, restart, and material-isolation checks passed.
- Custom-density regression retained agreement with 150 Python diagnostics.
- Custom02 parameter activity, reference equivalence, isolation, and restart
  tests passed at 1, 4, and 16 fixed substeps.
- `riva_sand_g0_test.cc` tested G0=350/reference/700 with Pa/kPa scaling,
  analytical shear/bulk moduli, exact reference preservation, rejection of
  invalid stiffness values, and exact cyclic equivalence to an independent
  elastic-scale override in both kernels.
- Documentation frozen-lockfile install, lint, TypeScript check, build, and
  all 30 rendered/source-wiring tests passed.

The successor suite now includes the G0 test. Run:

```sh
bash benchmark_inputs/nonlinear_material_point/run_riva_intermediate_bias_research_suite.sh \
  /path/to/rivasand_intermediate_bias_flow_port /path/to/c++
```

This record covers native host kernels and static parser/CUDA wiring review.
No full CUDA/MPI executable or mesh run was performed locally. Compile on the
target NVIDIA/MPI system and run the migrated material-point fixtures before
field analysis; cross-device behavior and full checkpoint I/O remain target-
system acceptance checks. Non-reference G0 values are not a new calibration,
and these tests do not establish improved LEAP acceleration agreement.
