# Compact RIVA-Sand inputs and velocity-derived stiffness

The standard and Custom input interfaces are distinct again. This supersedes
the full-row-for-all input contract in `RIVA_G0_INPUT_VALIDATION.md`.
OpenSees inputs are unaffected by this Hercules change.

## Standard names: Dr only, with optional G0

For `RIVASand` or `RIVASand02`, this example uses mesh-derived stiffness:

```text
material_model = RIVASand02
material_properties_count = 2
material_properties_list =
# Vs_limit Dr
  150 0.663
  450 0.663
```

Append a positive, dimensionless G0 to **every** row for an explicit override:

```text
material_properties_list =
# Vs_limit Dr G0
  150 0.663 483.48301127222084
  450 0.663 483.48301127222084
```

This reference G0 retains the former reference stiffness exactly. Standard
names use their generation's reference physical parameters for all other
fields: M=1.25, kd=1.125, h=122.44207260468994, m=0.945, zeta=0.025,
e_max=0.78, e_min=0.51, Q=10, R=1.5 and n_G=0.65. These are Ottawa F65
reference values, not a site-specific calibration.

## Custom names: full material input

`RIVASandCustom` and `RIVASandCustom02` require thirteen columns:

```text
material_model = RIVASandCustom02
material_properties_count = 2
material_properties_list =
# Vs_limit Dr G0 M kd h m zeta e_max e_min Q R n_G
  150 0.663 483.48301127222084 1.25 1.125 122.44207260468994 0.945 0.025 0.78 0.51 10 1.5 0.65
  450 0.663 483.48301127222084 1.25 1.125 122.44207260468994 0.945 0.025 0.78 0.51 10 1.5 0.65
```

Set **G0=0 in every row** to request automatic mesh stiffness while retaining
the other custom parameters. Zero is a selector, never a zero elastic modulus.
The same selector is accepted in the compact three-column format.

Use one numeric row per line, equal row widths, and increasing Vs_limit.
Blank lines and `#` comments are allowed. Do not mix zero and positive G0
rows: interpolating between an automatic selector and a physical coefficient
would be ambiguous. Positive G0 and other continuous inputs are linearly
interpolated by actual element Vs. All fields must be finite; Dr and n_G are
in [0,1], M/kd/h/Q positive, m/zeta/R nonnegative, and e_max>e_min>=0 with
the reference anchor e=0.601 inside the interval.

## What automatic stiffness matches

- Target initial elastic shear modulus: **mesh bulk density × mesh Vs²**,
  in Pa for kg/m³ and m/s. Vs_limit is not the target velocity; buoyant
  density is not used.
- Each integration point first admits its actual compressive effective
  geostatic stress without changing it. Approximate=yes uses the approximate
  field; ramped no uses the numerically solved stress at dynamic activation.
- The calibration evaluates the complete initialized modulus of the selected
  kernel, including pressure guards and initial bias-dependent stiffness
  factors, and adjusts only its dimensional elastic scale.
- That scale is then fixed during shaking. It is not fitted again as pressure,
  strain or cyclic memory evolves. Pressure-dependent degradation remains.
- Bulk stiffness follows the unchanged nu=0.30; matching Vs does **not**
  simultaneously match mesh Vp or an algorithmic elastoplastic tangent.
- Physical mean effective stress must still pass the existing admission
  checks. No extra surcharge, pressure shift or silent clipping is introduced.
  A high Vs at very small confinement can imply a very large reference G0;
  review that velocity/stress profile instead of assuming it is calibrated.

Automatic matching is an initialization option, not a new cyclic calibration.
It can change both shear and pore-pressure histories relative to fixed
reference G0. It does not establish improved strong-motion accelerations.

## Restart and runtime

One extra double per element stores the target shear modulus (zero for explicit
G0), and one per integration point caches the fitted elastic scale. The latter
is saved/restored with constitutive history. Updates load the cache; fitting
occurs only at first dynamic initialization, not per substep. There are no new
local iterations. Actual CUDA register/runtime impact still requires measurement.

On the branch-reversal research branch, checkpoint markers are **4 for original
RIVA-Sand and 7 for 02**. The earlier velocity-stiffness release used marker 6
for 02; the [changed branch reference](RIVA_BRANCH_REVERSAL_VALIDATION.md) requires
fresh initialization. Older
checkpoints are incompatible even with explicit G0: rerun initialization.
Keep the same model selection, mesh and complete material inputs on restart.
The saved elastic-scale layout and calibration hashes are unchanged; the six
02 reversal-reference values now store accumulated branch strain.

## Verification and target-system commands

The new host regression `riva_sand_velocity_test.cc` covers compact/Custom
parsing and rejection, automatic versus initialized elastic stiffness,
Pa/kPa scaling, loose/dense/bias/low-pressure cases, per-point isolation,
and fixed-cache restart at 1/4/16 substeps. It is part of:

```sh
bash benchmark_inputs/nonlinear_material_point/run_riva_intermediate_bias_research_suite.sh \
  /path/to/rivasand_intermediate_bias_flow_port /path/to/c++
```

The material-point subcommand has no mesh. Automatic mode therefore requires
`material_point_riva_density` (kg/m³) and `material_point_riva_vs` (m/s), plus
the existing initial effective-stress arrays. Standard material-point constants
are `Dr [G0]`; Custom constants are the twelve values without Vs_limit.
Two synthetic automatic fixtures are included in the full CPU/CUDA suite:

```sh
mpirun -np 1 /path/to/psolve --nonlinear-material-point riva_sand_velocity_synthetic.in
mpirun -np 1 /path/to/psolve --nonlinear-material-point riva_sand02_velocity_synthetic.in
```

Run from this benchmark directory on the NVIDIA/MPI target system. The
existing explicit-G0 fixtures retain the reference histories. Local host
checks and source review do not replace CUDA/MPI compilation, both geostatic
mesh workflows, full checkpoint I/O, or subsequent field qualification.

Local verification on 2026-09-08 passed: all five original reference histories
(31,742 full-output values at the established 1e-6 tolerances), all six
successor histories including all 236 output columns (established 1e-8/1e-10
tolerances), state/admission/plastic-gate/BrickUP host contracts, 150 custom-
density diagnostics, and the new velocity/input/restart tests. All eight
current RIVA-Sand material-point fixture rows parse successfully. Documentation
install, lint, TypeScript check, production build and all 30 regression tests
passed. No golden histories or frozen calibration records were edited.
