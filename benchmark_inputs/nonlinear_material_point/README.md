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
  Borja--Amies root solver, while its GPU history uses the CUDA root solver. It
  distinguishes unloading through zero total strain from a repeated, true
  zero strain increment.
- `vonmises_bah_synthetic.in` exercises the hyperbolic Borja--Amies hardening
  law.
- `vonmises_gqh_synthetic.in` exercises the dedicated CPU GQH backbone
  inversion and the corresponding CUDA inversion.
- `vonmises_mkz_synthetic.in` exercises the MKZ backbone inversion.
- `vonmises_ro_synthetic.in` exercises the Ramberg--Osgood hardening law.

BAE, BAH, GQH, MKZ, and Ramberg--Osgood are dispatched directly to their
incremental bounding/backbone update. The model uses the strain increment and
stored stress/history to identify loading and reversal; it does not use the
generic total-strain yield predictor applied to conventional yield-surface
models. A true zero strain increment preserves stress and history. This follows
the update boundary used by OpenSees
[`MultiaxialCyclicPlasticity`](https://github.com/OpenSees/OpenSees/blob/master/SRC/material/nD/cyclicSoil/MultiaxialCyclicPlasticity.cpp),
which owns its zero-increment and loading/reversal decisions inside the cyclic
integrator.

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
