# RIVASand02 reversal and BiasVolume selectors

Both `RIVASand02` and `RIVASandCustom02` accept these optional scalars in
`parameters.in` and `psolve --nonlinear-material-point` inputs. Omission selects:

```ini
material_model = RIVASand02
rivasand02_reversaltype = 3
rivasand02_biasvolume = 1
```

They apply to every RS02 integration point and are not material-row columns.
Invalid values, duplicate entries, and explicit use with other models are
rejected. The retired `reversalType` and `rivasand02_reversalType` keys remain
rejected. Original `RIVASand` and `RIVASandCustom` behavior is unchanged.

| Reversal type | Reference |
| --- | --- |
| 1 | Previous valid host increment |
| 2 | Accumulated loading-branch strain |
| 3 (default) | Full-host elastic stress-ratio predictor |

Type 3 ports the unguarded rule from OpenSees_RIVA
`355d8f04295ab3e69f9f8fb98516b1eb8f817446`:
`(alpha_trial-alpha0):(alpha_trial-alpha)<0`, evaluated from the committed
mechanical state. All rules check once per host increment before material
substeps. Type 3 applies at most one event, in the first substep, and retains
the strain-based initialization rule before cyclic activation. The optional
OpenSees `reversalGuard` and `reversalLatch` are not imported.

| BiasVolume mode | Behavior |
| --- | --- |
| 0 | Original bias-driven reversible volume |
| 1 (default) | Disable bias-driven reversible volume; retain separate phase-transformation volume |
| 2 | Original bias volume with the upstream gated correction to its contractive mean contribution |

Mode 2 uses the same fixed-bias, reference-pressure, effective-pressure and
plastic-activity gates as the upstream revision. Its committed activity is
held fixed within substeps and updated once after each host increment. It
adds one double per point (`field_bias_mean_activity`); modes 0/1 leave it
zero. No per-step allocation, extra exchange, or device-to-host transfer is
introduced. Both selectors are broadcast during setup, reported at startup,
and used by the shared field and CPU/GPU material-point constitutive core.

The 257-value reference fingerprint is
`16c8567fd2370c562899b9ec00cc2fc313be137f239cde5d0e22cc1cd3eee800`.
This identifies the upstream reference record, with execution-mode overrides
recorded separately. User material rows are unchanged. Mode 2 is an optional
research choice; the port alone does not establish stability or convergence.

## Restart behavior

Kernel revision 5 uses 64 logical fields and 139 flattened values. The first
138 retain their order; `field_bias_mean_activity` is appended. All previous
RS02 checkpoint markers through 9 are rejected. Restart requires both
selectors, model, mesh, and material inputs unchanged. The marker is
`10 + 3*BiasVolume + (reversaltype-1)`:

| BiasVolume | Reversal 1 | Reversal 2 | Reversal 3 |
| --- | --- | --- | --- |
| 0 | 10 | 11 | 12 |
| 1 | 13 | 14 | 15 (default) |
| 2 | 16 | 17 | 18 |

Rerun old cases from initialization. To reproduce an earlier response rule,
explicitly select mode 0 and its reversal type. Original-generation marker 4
is unchanged. New failure snapshots use `RIVASAND02_FAILURE_V2`, include both
selectors, and record 139 values; they are diagnostics, not checkpoints.

## Verification

On 2026-10-09, all **189 host comparisons** against the exact OpenSees revision
above matched all 139 fields bit for bit, with finite states and clean
ASan/UBSan results. The matrix covers three reversal rules, three BiasVolume
modes, six existing strain paths plus a high-pressure low-bias variant, and
three substep settings. Mode 2 changes the original 138 response fields in
the high-pressure low-bias variant; that active correction also matches upstream.
The selector test verifies omitted defaults, strict parsing, phase-volume
retention in mode 1, same-mode restart continuation, and rejection of all
cross-selector and old-layout checkpoints. The historical native suite also
passes with its old controls pinned explicitly.

Horizon CUDA compilation passed. Job 53060 passed 23 CPU/GPU material-point
comparisons: standard/Custom02 reference rows at three resolutions in each
mode, three biased cases exercising the correction, and two inputs omitting
both selectors. Every complete state comparison passed at rtol=1e-8 and
atol=1e-10. Standard/Custom02 histories and omitted/explicit defaults were
byte-identical on each backend.

The same job completed the full 30-s Horizon BrickUP/CPML case with both
selectors omitted (type 3/mode 1), dt=0.005 s and ten material substeps. The
launcher and Slurm job exited successfully; 6000 station samples extend to
29.995 s. Displacement, velocity, acceleration, pressure, strain and stress
were finite. The existing undefined surface-station Ru column remains NaN.
The original case, earlier comparison outputs, and executables were preserved.
This validates completion of this case, not timestep/mesh convergence or
experimental agreement. Native MPI checkpoint file I/O was not tested. Comparison helpers and upstream snapshots stay
local-only under `.local-diagnostics/rs02-biasvolume-20261009/`.

The earlier type-3-only change passed 54 mode-0 host comparisons of the shared
138 values and six GPU comparisons. Its field reproduction below retains
that scope and does not validate the new mode-1 default.

### Earlier mode-0 field reproduction on 2026-10-09 (35c539c)

The maintainer's coupled BrickUP/CPML case used `RIVASand02`, automatic
stiffness from mesh density and Vs, Dr=0.66296, a 0.005-s global timestep,
and ten constitutive substeps. Both the original executable (job 53057)
and the diagnostic build (job 53056) failed with type 1 near 8.3--8.6 s.
The diagnostic build captured rank 8, element 697232, point 0 at step 1729:
the largest incoming strain component was approximately 7.23e32 and the
stored void ratio was already 1.47e61. The cached E_ref was finite and
positive. Replaying that exact state and increment with Hercules and
OpenSees `355d8f042` produced identical rejection traces, with void-ratio
overflow in substep 8. The near-1-Pa printed pressure was therefore not a
complete description of the already unphysical input state.

The setup scan reported a CPML critical dt of 0.015625 s and an undrained
minimum h/Vp of 0.013224 s. The failed timestep was below both values;
these initial linear estimates do not establish stability of the evolving
nonlinear coupled state.

Changing only the reversal selector to type 3 in the diagnostic build
completed the full 30-s case in job 53056, with exit status zero and
6000 station samples through 29.995 s. Station displacement, velocity,
acceleration, pressure, strain, and stress were finite; the surface-station
Ru column remained NaN. The original case, executable, and outputs were
preserved. A separate 12-s preliminary run completed, but its shortened
source history did not reproduce the full-duration failure and is not
used as evidence of a fix.

Keeping type 1 and changing only the global timestep to 0.0025 s also
completed the full 30-s case (job 53058, exit status zero), with 12000
station samples through 29.9975 s. Its primary station fields were finite,
with the same undefined surface Ru column. The constitutive substep count
remained ten. This demonstrates timestep sensitivity of the failure; it
does not isolate the first instability mechanism or establish convergence.

This establishes completion of this particular case with type 3, not
experimental validation or timestep/mesh convergence. Type 1 was the
compatibility default at that revision, before the defaults changed to type 3
and BiasVolume mode 1. No user material coefficients or pore-pressure coupling
law were changed in that earlier comparison. The new failure snapshot helps distinguish an invalid
incoming coupled state from a constitutive port discrepancy.

### CUDA fixture utility

`run_riva_branch_reversal_cuda.py` requires Python 3 and `--output DIR`.
Without `--psolve`, it only generates six fixtures; executing them requires
an NVIDIA/MPI build and optionally `--launcher 'mpirun -np 1'`.
`--reversal-type` accepts 1/2/3 (default 3); `--bias-volume` accepts 0/1/2
(default 1). Existing fixture inputs, logs, CSVs and `summary.json` in that
output directory are overwritten. Use a separate directory for each pair.

```sh
python3 benchmark_inputs/nonlinear_material_point/run_riva_branch_reversal_cuda.py \
  --psolve /path/to/psolve --launcher 'mpirun -np 1' \
  --reversal-type 3 --bias-volume 1 --output /path/to/results/type3-mode1
```

The runner checks complete host/GPU state parity, finite outputs, exact
zero-increment preservation, expected reversal counts, selector diagnostics,
and standard/Custom02 equality. It does not test MPI checkpoint file I/O,
field stability, or convergence. Historical native oracle helpers explicitly
pin BiasVolume mode 0; their default type-1 replay remains a historical test,
not the production default.

### Historical type-1/type-2 checks

Host and documentation checks passed on 2026-09-22. The
[verification record](RIVA_REVERSAL_SELECTOR_RESULTS.json) includes source
hashes and the twelve unchanged before/after history hashes. Documentation
lint, TypeScript, production build and all 33 parser/rendering tests pass.

The focused host test checks omitted/explicit defaults, invalid and duplicate
inputs, rejection of the retired key, both reversal modes, rejected updates,
restart continuation and cross-type restart rejection. It uses the same strict scalar reader as both
native input paths. The existing branch-reference test explicitly selects
type 2 and retains its mixed-path refinement, pure-shear, rotation, unit,
trial-isolation and checkpoint checks.

The full native suite and byte-identical comparisons below were run for
`36af020`, before the input-key rename. The constitutive kernel and native
replay driver are unchanged; the focused host test, generated-input audit
and documentation checks were repeated after the rename.

The native suite replays the six frozen input histories with default type 1
and explicit type 2. The original frozen oracle checks type 1; the separate
branch-reference oracle checks type 2. The old reference data is preserved.
Direct comparisons with executables built from unchanged `9975658` check
that all six complete CSV histories are byte-identical for each mode: the
old offline previous-increment path for type 1 and the old default for type 2.

Run the host selector test with a C++17 compiler:

```sh
c++ -std=c++17 -O2 riva_sand02_reversal_selector_test.cc -o /tmp/riva-selector-test
/tmp/riva-selector-test
```

The existing host suite accepts the unchanged private handoff path and the
optional private branch-reference verifier. These private references are
not copied into the repository.

CUDA/MPI execution and field restart verification require the target system.
The fixture generator can check both standard and Custom02 names, complete
CPU/GPU states, zero increments and expected smooth-path reversal counts:

```sh
python run_riva_branch_reversal_cuda.py --reversal-type 1 \
  --output /tmp/riva-type1 --psolve /path/to/psolve --launcher 'mpirun -np 1'
python run_riva_branch_reversal_cuda.py --reversal-type 2 \
  --output /tmp/riva-type2 --psolve /path/to/psolve --launcher 'mpirun -np 1'
python run_riva_branch_reversal_cuda.py --reversal-type 3 \
  --output /tmp/riva-type3 --psolve /path/to/psolve --launcher 'mpirun -np 1'
```

Omit `--psolve` to generate inputs without claiming GPU execution. These are
verification paths, not newly calibrated field inputs. The selector preserves
each rule's prior numerical behavior; it does not establish timestep or mesh
independence or validate the free-field calibration in Hercules.
