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

Use `rivasand02_reversaltype = 2` to preserve the behavior of Hercules `main` at
`9975658`, including the previously merged loading-branch correction.
Only the integers 1 and 2 are accepted. Type 3, fractional values, duplicate
entries and explicit use with another material model are rejected.
The former `reversalType` and `rivasand02_reversalType` keys are rejected;
rename either to `rivasand02_reversaltype` in existing inputs, preserving its value.
Original `RIVASand` and `RIVASandCustom` behavior is unchanged.

Both rules check reversal once per host increment before the fixed material
substeps. The existing equations, coefficient values, deadbands, integration
order and 138-value state layout are unchanged. The selector does not add
OpenSees `reversalLatch` or `BiasVolume`. The 249-value reference calibration
hash is unchanged; `rivasand02_reversaltype` is an exposed execution choice rather than
a recalibrated material coefficient.

The root input value is broadcast to all solver ranks and carried in the
existing GPU configuration. The field kernel and CPU/GPU material-point
adapters use the same selected update. No new per-point allocation, device
transfer or synchronization is introduced during time stepping. Startup and
material-point output identify the selected rule and checkpoint marker.

## Restart behavior

The two references occupy the same six history slots but have different
meanings. Cross-type restart is rejected before state is loaded. Type 2 keeps
marker 7 and can continue an existing marker-7 history with explicit
`rivasand02_reversaltype = 2`, subject to the usual unchanged model, mesh and material
inputs. The new default type 1 writes marker 8. Marker-6 and earlier RS02
checkpoints remain unsupported; rerun their initialization. Original-generation
marker 4 is unchanged. A change of the default does not silently convert a
saved branch reference into a previous-increment direction.

## Verification

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
```

Omit `--psolve` to generate inputs without claiming GPU execution. These are
verification paths, not newly calibrated field inputs. The selector preserves
each rule's prior numerical behavior; it does not establish timestep or mesh
independence or validate the free-field calibration in Hercules.
