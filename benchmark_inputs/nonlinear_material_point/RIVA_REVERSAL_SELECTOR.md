# RIVASand02 reversal selector

Both `RIVASand02` and `RIVASandCustom02` accept this optional scalar in
`parameters.in` and in `psolve --nonlinear-material-point` inputs:

```ini
material_model = RIVASand02
rivasand02_reversaltype = 1
```

Omitting `rivasand02_reversaltype` also selects **1**. It is a separate setting, not an
additional column in `material_properties_list` or `material_point_constants`.
The selector applies to all RS02 integration points in the simulation.

| Value | Reversal reference | Default | Checkpoint marker |
|---|---|---|---|
| 1 | Previous valid host increment | Yes | 8 |
| 2 | Accumulated loading-branch strain | No | 7 |
| 3 | Full-host elastic stress-ratio predictor | No | 9 |

Use `rivasand02_reversaltype = 2` to preserve the behavior of Hercules `main` at
`9975658`, including the previously merged loading-branch correction.
Type 3 selects the unguarded host stress-ratio predictor ported from
OpenSees_RIVA `355d8f04295ab3e69f9f8fb98516b1eb8f817446`. It evaluates
`(alpha_trial-alpha0):(alpha_trial-alpha)<0` once from the full host increment
and the committed mechanical state. It forces at most one reversal in the
first substep, with the strain-based initialization rule before cyclic activation.
Only integers 1, 2 and 3 are accepted. Fractional values, duplicate entries and explicit use with another material model are rejected.
The former `reversalType` and `rivasand02_reversalType` keys are rejected;
rename either to `rivasand02_reversaltype` in existing inputs, preserving its value.
Original `RIVASand` and `RIVASandCustom` behavior is unchanged.

All rules check reversal once per host increment before the fixed material
substeps. Type-1/type-2 behavior and the 138-value state layout are unchanged. Type 3
uses the literal OpenSees rule without the optional `reversalGuard`. The selector does not add
OpenSees `reversalLatch` or `BiasVolume`. The 249-value reference calibration
hash is unchanged; `rivasand02_reversaltype` is an exposed execution choice rather than
a recalibrated material coefficient.

The root input value is broadcast to all solver ranks and carried in the
existing GPU configuration. The field kernel and CPU/GPU material-point
adapters use the same selected update. No new per-point allocation, device
transfer or synchronization is introduced during time stepping. Startup and
material-point output identify the selected rule and checkpoint marker.

## Restart behavior

The type-1 and type-2 references occupy the same six history slots but have different
meanings. Cross-type restart is rejected before state is loaded. Type 2 keeps
marker 7 and can continue an existing marker-7 history with explicit
`rivasand02_reversaltype = 2`, subject to the usual unchanged model, mesh and material
inputs. The default type 1 writes marker 8; type 3 writes marker 9 and cannot
continue type-1/type-2 checkpoints. Marker-6 and earlier RS02
checkpoints remain unsupported; rerun their initialization. Original-generation
marker 4 is unchanged. A change of the default does not silently convert a
saved branch reference into a previous-increment direction.

## Verification

On 2026-10-09, the three rules were compared with the exact OpenSees revision
above on six prescribed strain histories, using each history's original
substep count and overrides of 10 and 50. All 54 host comparisons produced
identical values for the 138 shared state fields. AddressSanitizer and
UndefinedBehaviorSanitizer reported no errors. The Horizon CUDA build also passed. Job 53056 passed six CPU/GPU
material-point comparisons for standard and Custom02 reference rows at
128, 512 and 2048 steps per cycle, with type 3 and four material substeps.
All registered four reversals over two cycles, and complete state histories
agreed under the runner's relative/absolute tolerances of 1e-8/1e-10.
These material-point checks do not resolve the separate coupled-field failure.
The local comparison harness and upstream snapshot remain ignored in
`.local-diagnostics/rivasand02-horizon-20261009/`.

### Horizon field reproduction on 2026-10-09

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
experimental validation or timestep/mesh convergence. Type 1 remains the
compatibility default; no material coefficients or pore-pressure coupling
law were changed. The new failure snapshot helps distinguish an invalid
incoming coupled state from a constitutive port discrepancy.

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
