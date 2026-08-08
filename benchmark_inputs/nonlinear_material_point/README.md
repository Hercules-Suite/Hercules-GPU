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

For `material_model=ManzariDafalias`, the same key instead contains 17 values:

```text
G0 nu e_init Mc c lambda_c e0 ksi P_atm m h0 ch nb A0 nd z_max cz
```

This follows the standard Manzari--Dafalias parameter order except that the
material-point test does not need density. The production 17-column
`material_properties_list` starts with `Vs_limit`, omits `nu`, and obtains both
Poisson's ratio and density from the element's `Vp`, `Vs`, and `rho`.

- `vonmises_ep_synthetic.in` is the shortest smoke test for yielding and
  plastic-state advancement.
- `vonmises_fa_quartic_synthetic.in` exercises the FA quartic corrector.
- `vonmises_fam_quartic_reversal_synthetic.in` exercises both the FAM quartic
  corrector and its load-reversal recalibration.
- `vonmises_bae_legacy_id2_synthetic.in` uses the constants currently assigned
  to legacy Istanbul material ID 2. Its CPU history enters the GSL-backed
  Borja--Amies root solver, while its GPU history uses the CUDA root solver. It
  distinguishes unloading through zero total strain from a repeated, true
  zero strain increment.
- `vonmises_bae_virgin_onset_synthetic.in` starts with a purely volumetric
  increment, accumulates sub-millipascal virgin deviatoric stress, crosses the
  numerical onset floor, and applies another purely volumetric increment after
  cyclic history exists. It guards the production wavefront case in which the
  scalar kappa equation otherwise has no finite or numerically resolvable root.
- `vonmises_bah_synthetic.in` exercises the hyperbolic Borja--Amies hardening
  law.
- `vonmises_gqh_synthetic.in` exercises the dedicated CPU GQH backbone
  inversion and the corresponding CUDA inversion.
- `vonmises_mkz_synthetic.in` exercises the MKZ backbone inversion.
- `vonmises_ro_synthetic.in` exercises the Ramberg--Osgood hardening law.
- `manzari_dafalias_modified_euler_synthetic.in` establishes isotropic
  confinement, applies a cyclic tensor-shear path, and checks the same
  allocation-free Modified-Euler implementation on host and CUDA. Its CSV adds
  alpha, alpha-in, fabric, void ratio, plastic multiplier, absolute yield
  residual, local integration error, and accepted-substep diagnostics.
- `manzari_dafalias_seeded_onset_synthetic.in` starts from a committed
  isotropically consolidated state and crosses first yield in the next strain
  increment. It verifies the elastic intersection followed by Modified-Euler
  integration of only the plastic remainder.

Manzari--Dafalias reuses the production nonlinear integration controls rather
than adding scheme-specific input keys: `backbone_errTol` is the dimensionless
local Modified-Euler error target, `no_substeps` caps accepted adaptive
substeps (with rejected trials also bounded), and `error_tolerance` is the
maximum absolute yield-surface residual in Pa after a plastic stress
correction. The elastic predictor retains the one-sided admissibility check
`f_trial <= error_tolerance`, so a genuinely elastic state may remain inside
the surface. An increment that starts strictly inside but ends outside is split
at its elastic yield intersection before integrating the plastic remainder,
matching the OpenSees transition. Elastic stages preserve alpha, fabric, and
plastic strain; only
plastic flow evolves those history variables. After an accepted
plastic substep, the stress-correction loop
permits at most `min(no_substeps, 50)` iterations. If it cannot reduce the
absolute residual to the tolerance, Hercules rejects the constitutive update
and the production GPU error check aborts with Manzari--Dafalias failure code
3 rather than continuing from an off-surface state. The model
is currently GPU-only in production; the host invocation here is a regression
driver for the shared constitutive arithmetic, not a second production
backend. This host/CUDA parity check is not independent physical validation;
compare the recorded stress and state history with OpenSees or another trusted
implementation before using a calibration.

The independent OpenSees characterization uses the same 21-step path:

```bash
OpenSees opensees_manzari_dafalias_synthetic.tcl opensees_md.csv
python3 compare_opensees_manzari_dafalias.py \
  --hercules-csv manzari_dafalias_modified_euler_synthetic.cpu.csv \
  --opensees-csv opensees_md.csv
```

The OpenSees driver uses kPa, the native calibration scale of its
low-pressure dilatancy regularization; the comparator converts its stresses to
Pa before comparing stress, alpha, fabric, alpha-in, and void-ratio histories.
The comparison reports differences but does not impose a pass threshold,
because it is an independent model characterization rather than a host/device
identity check.

## Canonical Manzari--Dafalias comparisons

Two longer OpenSees drivers characterize the model beyond the short cyclic
CPU/CUDA fixture:

- `opensees_manzari_dafalias_drained_triaxial.tcl` first consolidates a
  stabilized single-point unit brick isotropically to 100 kPa and then applies
  15% monotonic axial
  compression while keeping both lateral face loads at 100 kPa. It exercises
  critical-state strength and the contraction-to-dilation response. Failed
  outer increments use recursive bisection; each attempt is limited to the 20
  Newton iterations used by the published OpenSees triaxial example. For this
  uniform affine path, single-point integration avoids evaluating eight
  identical material states in a full-integration brick.
- `opensees_manzari_dafalias_undrained_cyclic_shear.tcl` consolidates
  isotropically and then applies five constant-volume shear cycles with 0.5%
  engineering-shear amplitude. It calculates the excess pore pressure needed
  to enforce undrained conditions as `delta_u = p'_0 - p'` and reports
  `r_u = delta_u/p'_0`, effective-stress reduction, and stress reversals.

The drained problem is stress-controlled in OpenSees. The conversion tool
places the lateral strain that OpenSees solved into a Hercules material-point
input. Hercules therefore sees the identical strain history; departure of its
radial stress from the prescribed confinement measures disagreement in the
model's contraction/dilation tendency without mixing in a second outer
stress-control algorithm.

Both canonical OpenSees drivers use a forced-elastic consolidation stage and
then enable the elastoplastic material. The conversion tool seeds Hercules
with the recorded stress, strain, alpha, alpha-in, and fabric at that stage
boundary and replays only the common elastoplastic history. Replaying the
elastic consolidation strains from zero through Hercules' pressure-dependent
elastoplastic stage would create a different initial effective stress and is
not a valid constitutive comparison. The first generated history row is the
seeded state itself, so it also verifies zero-increment state preservation.

The consolidation boundary is theoretically isotropic, but an OpenSees finite
element solve may leave roundoff-scale shear stress or directional internal
state. The conversion tool verifies that stress isotropy and zero alpha,
alpha-in, and fabric hold within strict cleanup tolerances, then seeds the
exact mean normal stress with zero directional state. This prevents
roundoff-only symmetry breaking from selecting different host/device adaptive
branches while retaining the recorded post-consolidation strain path.

Run either standalone driver on a machine with OpenSees:

```bash
OpenSees opensees_manzari_dafalias_drained_triaxial.tcl \
  opensees_md_drained_triaxial.csv
OpenSees opensees_manzari_dafalias_undrained_cyclic_shear.tcl \
  opensees_md_undrained_cyclic_shear.csv
```

Copy the CSVs to the Hercules machine, generate the matching material-point
inputs, and run them through `psolve`:

```bash
python3 generate_hercules_manzari_dafalias_path.py \
  --opensees-csv opensees_md_drained_triaxial.csv \
  --output manzari_md_drained_triaxial.in \
  --output-prefix manzari_md_drained_triaxial
mpirun -np 1 ../../quake/forward/psolve \
  --nonlinear-material-point manzari_md_drained_triaxial.in

python3 generate_hercules_manzari_dafalias_path.py \
  --opensees-csv opensees_md_undrained_cyclic_shear.csv \
  --output manzari_md_undrained_cyclic_shear.in \
  --output-prefix manzari_md_undrained_cyclic_shear
mpirun -np 1 ../../quake/forward/psolve \
  --nonlinear-material-point manzari_md_undrained_cyclic_shear.in
```

Then create invariant, state, confinement, and pore-pressure comparisons:

```bash
python3 compare_opensees_manzari_dafalias_canonical.py \
  --mode drained \
  --hercules-csv manzari_md_drained_triaxial.cpu.csv \
  --opensees-csv opensees_md_drained_triaxial.csv

python3 compare_opensees_manzari_dafalias_canonical.py \
  --mode undrained \
  --hercules-csv manzari_md_undrained_cyclic_shear.cpu.csv \
  --opensees-csv opensees_md_undrained_cyclic_shear.csv
```

The comparator is report-only by default. After documenting calibration- and
path-specific acceptance values, CI or a benchmark script may add
`--max-stress-nrmse`, `--max-radial-confinement-error-kpa`, or
`--max-ru-difference`.

Generated canonical histories use `material_point_compare_rtol=1e-6` and
`material_point_compare_atol=1e-8` for CPU/CUDA state parity over hundreds of
adaptive constitutive increments. Each tensor component uses the largest
resolved component of the same six-component tensor as its relative scale;
this avoids ill-conditioned relative errors for physically inactive
components without loosening scalar or dimensionless-state checks. CPU and GPU
yield residuals are checked independently against the absolute
`error_tolerance`. Modified-Euler local-error estimates and accepted-substep
counts are likewise checked independently against `backbone_errTol` and
`no_substeps`; host and device arithmetic need not select identical valid
adaptive substep histories. Per-step plastic multipliers are nonnegative
integration diagnostics; the evolved plastic-strain tensor remains part of
the direct state comparison.

The constant-volume test is an undrained *constitutive* characterization. Its
pore pressure is inferred from effective mean stress under the stated constant
total-pressure condition. It does not exercise Hercules U-P pressure degrees
of freedom, storage, permeability/Darcy flow, boundary drainage, or coupling.
A separate BrickUP production benchmark is required before claiming that the
complete liquefaction pathway is validated.

BAE, BAH, GQH, MKZ, and Ramberg--Osgood are dispatched directly to their
incremental bounding/backbone update. The model uses the strain increment and
stored stress/history to identify loading and reversal; it does not use the
generic total-strain yield predictor applied to conventional yield-surface
models. A true zero strain increment preserves stress and history. This follows
the update boundary used by OpenSees
[`MultiaxialCyclicPlasticity`](https://github.com/OpenSees/OpenSees/blob/master/SRC/material/nD/cyclicSoil/MultiaxialCyclicPlasticity.cpp),
which owns its zero-increment and loading/reversal decisions inside the cyclic
integrator.

Because these von Mises backbones are pressure-insensitive, a purely
volumetric increment updates only hydrostatic stress and preserves all cyclic
history; it does not enter the scalar kappa solver. At untouched virgin state
(`kappa=100`), deviatoric trial stress below `1.0e-3 Pa` uses the elastic
asymptote and accumulates incrementally until the backbone root is resolvable.
That numerical-onset rule is not used after cyclic history has been
established, so unloading and reloading remain model-controlled.

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

## OpenSees BAE comparison

`compare_opensees_bae.py` runs the dedicated
`vonmises_bae_opensees_synthetic.in` history through Hercules and an OpenSees
one-brick affine simple-shear test. It maps the Hercules BAE constants to
OpenSees `MultiaxialCyclicPlasticity` as follows:

```text
K = lambda + 2*mu/3    G = mu        Su = cohesion
Ho = 0                 h = psi0*G    m = m
beta = 1
```

`beta=1` is the closest match to the Hercules end-state hardening rule. Run a
second comparison with `--beta 0.5` to measure sensitivity to the integration
parameter normally recommended by OpenSees. OpenSees consumes engineering
shear, so the driver converts Hercules tensorial `epsilon_xy` to
`gamma_xy=2*epsilon_xy`.

`opensees_bae_synthetic.tcl` is a standalone version of the one-brick test. It
contains the mapped constants and complete 20-step engineering-shear history,
so Hercules and OpenSees do not need to be installed on the same machine. The
17-step `vonmises_bae_legacy_id2_synthetic` history from `run_suite.sh` is not
interchangeable with this dedicated 20-step cross-code fixture.

The OpenSees driver applies that history as a `Path` time series multiplying a
prescribed top-face displacement. It checks the actual control-node
displacement before recording each stress. This is essential at the repeated
zero point: advancing with `LoadControl 0` would release displacement control
and re-equilibrate the brick instead of applying a true zero strain increment.
Version-2 CSVs include `driver_version` and `domain_time`; the comparison tool
rejects older CSVs, which must be regenerated. Penalty enforcement can leave
roundoff-scale differences between requested and actual displacement. The
comparison accepts absolute shear-path error up to `1e-12`, reports the maximum
observed error, and retains the actual OpenSees strain in `comparison.csv`.

First generate or reuse the Hercules CPU history, then provide an OpenSees
executable that includes `MultiaxialCyclicPlasticity`:

```bash
python3 compare_opensees_bae.py \
  --fixture vonmises_bae_opensees_synthetic.in \
  --psolve ../../quake/forward/psolve \
  --opensees /path/to/OpenSees
```

To compare an existing Hercules CSV without rerunning `psolve`:

```bash
python3 compare_opensees_bae.py \
  --fixture vonmises_bae_opensees_synthetic.in \
  --hercules-csv /path/to/vonmises_bae_opensees_synthetic.cpu.csv \
  --opensees /path/to/OpenSees
```

For a split Horizon/Stampede3 run, first produce the matching Hercules history
on Horizon:

```bash
mpirun -np 1 ../../quake/forward/psolve \
  --nonlinear-material-point vonmises_bae_opensees_synthetic.in
```

Run the standalone case from the same directory on Stampede3. The first
command uses the closest `beta=1` mapping; the second is the OpenSees
`beta=0.5` sensitivity run:

```bash
OpenSees opensees_bae_synthetic.tcl
OpenSees opensees_bae_synthetic.tcl 0.5 opensees_bae_beta05.csv
```

Copy the resulting CSV back to Horizon and compare it without launching
OpenSees there:

```bash
python3 compare_opensees_bae.py \
  --hercules-csv vonmises_bae_opensees_synthetic.cpu.csv \
  --opensees-csv opensees_bae_beta_1.csv \
  --beta 1.0
```

The script writes the generated Tcl driver, raw OpenSees history,
component-by-component stress comparison, solver logs, and `summary.json` under
`opensees_bae_comparison_output/`. It reports normalized shear-stress RMSE and
the relative difference in integrated shear work. The default is report-only;
use `--max-shear-nrmse` or `--max-work-relative-error` only after choosing and
documenting case-appropriate acceptance thresholds.

This is a characterization comparison, not an equivalence proof. OpenSees uses
a coupled Newton update, strain splitting at the bounding surface, and resets
its trial mapping modulus and kappa on a committed near-zero deviatoric
increment. Hercules uses a scalar bracketed root and preserves its complete
state on a true zero increment. The supplied history includes a repeated zero
point followed by reloading so that difference is visible rather than hidden.

Run the parser/driver-generation check without either solver:

```bash
python3 compare_opensees_bae.py --selftest
```
