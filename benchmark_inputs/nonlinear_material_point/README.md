# Nonlinear material-point CPU/GPU test

This test drives one constitutive point with a prescribed strain history. It
does not build a mesh, read a source, or advance the wave solver. Its purpose is
to isolate CPU/CUDA constitutive parity before running an element or production
benchmark.

Build `quake/forward/psolve` normally on a CUDA system, then run the complete
suite with exactly one MPI rank:

```bash
cd benchmark_inputs/nonlinear_material_point
./run_suite.sh ../../quake/forward/psolve
```

To isolate one model, invoke the runner directly, for example:

```bash
mpirun -np 1 ../../quake/forward/psolve \
  --nonlinear-material-point vonmises_gqh_synthetic.in
```

The command runs the complete history through the production CPU
`material_update` routine and the CUDA material-update routine, writes
`<material_point_output_prefix>.cpu.csv` and `.gpu.csv`, and returns nonzero on
a status, finite-value, or tolerance mismatch.

`material_point_strain_history` contains tensor strains in this order:

```text
epsilon_xx epsilon_yy epsilon_zz epsilon_xy epsilon_yz epsilon_xz
```

Off-diagonal values are tensorial shear strains. For an engineering shear
history `gamma_xy`, provide `epsilon_xy = gamma_xy / 2`.

`material_point_constants` contains 24 values in this order:

```text
lambda mu alpha beta gamma cohesion phi dilatancy hardening
yield_strain psi0 m
theta1 theta2 theta3 theta4 theta5
beta_MKZ s_MKZ phi_MKZ gamma_ref_MKZ
alpha_RO eta_RO phi_RO
```

Stress-like values use Pa; angles use radians. Unused model-specific entries
must be present and may be zero. The supplied fixtures are synthetic numerical
regression cases, not production soil calibrations:

- `vonmises_ep_synthetic.in` is the shortest smoke test for yielding and
  plastic-state advancement.
- `vonmises_fa_quartic_synthetic.in` exercises the FA quartic corrector.
- `vonmises_fam_quartic_reversal_synthetic.in` exercises both the FAM quartic
  corrector and its load-reversal recalibration.
- `vonmises_bae_legacy_id2_synthetic.in` uses the constants currently assigned
  to legacy Istanbul material ID 2. Its CPU history enters the GSL-backed
  Borja--Amies root solver, while its GPU history uses the CUDA root solver.
- `vonmises_bah_synthetic.in` exercises the hyperbolic Borja--Amies hardening
  law.
- `vonmises_gqh_synthetic.in` exercises the dedicated CPU GQH backbone
  inversion and the corresponding CUDA inversion.
- `vonmises_mkz_synthetic.in` exercises the MKZ backbone inversion.
- `vonmises_ro_synthetic.in` exercises the Ramberg--Osgood hardening law.

The suite intentionally preserves the current CPU dispatch ordering as its
compatibility baseline. First establish CPU/CUDA agreement for all fixtures.
Treat any later change that sends BAE, BAH, GQH, MKZ, or Ramberg--Osgood ahead
of the generic CPU elastic-predictor return as a separate physics change: save
the pre-change CSV histories, update both backends consistently, and review
the changed unloading/reversal histories before accepting a new baseline.

For scientific validation, replace the constants and strain history with
values derived from the calibrated case being studied and retain its
geostatic, tolerance, substep, damping, and material-property provenance.

The material-point runner currently exercises rate-independent constitutive
updates, matching the implemented CUDA path. Do not use it to claim validation
of the unfinished rate-dependent branch.

The default comparison tolerances are `rtol=1e-6` and `atol=1e-9`; the fixture
may override them with `material_point_compare_rtol` and
`material_point_compare_atol`. Physical state fields are compared directly.
For the bounding-surface families, the CPU/GSL and CUDA root residuals are
instead checked independently against `error_tolerance`, because different
root algorithms may stop at different residuals within the same accepted
bound. A single point is a correctness test only and must not be used to
measure GPU speedup.
