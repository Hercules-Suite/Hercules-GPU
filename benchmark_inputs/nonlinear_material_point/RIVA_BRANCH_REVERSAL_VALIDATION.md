# Hercules 02 loading-branch reversal research

**Current interface:** the [rivasand02_reversaltype selector](RIVA_REVERSAL_SELECTOR.md)
defaults to type 1. Set `rivasand02_reversaltype = 2` to reproduce this historical
branch-reference study or continue its marker-7 checkpoints. The validation
record below describes the original type-2 port.

On `research/rivasand02-branch-reversal`, both `RIVASand02` and
`RIVASandCustom02` select the **kernel-revision-4 loading-branch research
successor**. This branch starts from Hercules `04a5cfa` and ports the reversal
reference correction from OpenSees_RIVA `e97c12fbc`. The numerical kernel
is shared by the field CUDA update and the host/CUDA material-point adapter.
Existing names and material-row formats are retained on this research branch.
This is not production acceptance or a new parameter calibration.

## Change and checkpoint migration

The prior detector saved the direction of the previous host increment. Smooth
combined shear and vertical strain can rotate that direction gradually enough
to miss reversals when the host step is refined. The new reference stores
accumulated deviatoric strain along the loading branch and resets when a
reversal is registered. Detection still occurs once per host increment before
fixed material substeps; existing activity and stress deadbands remain.
Failed trial updates do not replace the committed point state.

Only the 02 research event reference changes. Physical inputs, calibration
hash, pressure-dependent stiffness, plastic-activity gate, and the original
`RIVASand`/`RIVASandCustom` kernel are unchanged. Hercules's fluid-pressure
and damping ownership are unchanged. No SSPbrickUP patch is copied into
Hercules. No OpenSees `BiasVolume` or `reversalLatch` option is imported;
the Hercules revision number does not imply those OpenSees features.

Both 02 names now write and require **checkpoint marker 7**. All older markers,
including 6, are rejected. Rerun those cases from initialization even with
explicit G0. Original-generation checkpoints retain marker 4. The flattened
state remains 138 values, but the six historical
`last_host_deviatoric_strain_direction` slots now hold accumulated branch
strain. The old unit direction cannot reconstruct that history. The saved
per-point elastic scale and fluid state remain in the checkpoint payload.

The actual checkpoint reader uses the revision/size/calibration check exercised
by the new host tests. New-format byte restart is checked before/after reversals,
during a loading branch, and at a cycle boundary for both 02 selectors. Restored
and uninterrupted trajectories are identical. Old/future markers, wrong state
widths, and wrong calibration hashes are rejected. This is not verification
of GPU checkpoint file I/O or multi-rank MPI restart.

## Local evidence

- A separately identified private Python successor replays all six existing
  inputs with the correction. All **236 columns** match the current native
  driver at the existing `rtol=1e-8`, `atol=1e-10`; stress/state are finite.
  The revision-3 package and its golden histories are not rewritten.
- The offline replay driver's explicit `--legacy-reversal-reference` option
  reproduces the historical oracle in core and full modes. It is not a field
  input or acceptance evidence for the corrected mode.
- The corrected mixed path retains four reversals at 128, 512, 2,048 and 8,192
  steps per cycle. The previous reference gives four at 128 and zero at finer
  resolutions. Maximum tensor-stress differences from the 8,192-step reference
  decrease under refinement (four material substeps):

| Steps per cycle | Maximum difference (kPa) |
| --- | --- |
| 128 | 1.19938 |
| 512 | 0.229917 |
| 2,048 | 0.0479737 |

- Pure shear matches the previous history exactly except for the six redefined
  reference values. Additional 1/4/16-substep checks, rotation, Pa/kPa scaling,
  discarded trials, invalid-trial rejection, and checkpoint continuation pass.
- On dense and loose prescribed paths shared with the corrected OpenSees
  material, all 4,096 updates have identical reversal counts and stress
  differences below **2.1e-13 kPa**. These are the tested zero-static-bias paths;
  the two hosts retain different option sets.
- Standard02/Custom02 equality, input activity, custom-density bounds,
  automatic stiffness, cached elastic scale, stage/admission, and BrickUP
  ownership host tests pass. Original RIVA-Sand passes its five full-output
  histories (31,742 values at established 1e-6 tolerances) and state tests.
- Documentation build, all 31 rendered/source checks, TypeScript and lint pass.
  See [the compact results](RIVA_BRANCH_REVERSAL_RESULTS.json).

## Reproduction and target validation

Run the host suite with the unchanged private revision-3 handoff and the
separately identified corrected verifier:

```sh
bash benchmark_inputs/nonlinear_material_point/run_riva_intermediate_bias_research_suite.sh \
  /path/to/rivasand_intermediate_bias_flow_port /path/to/c++ \
  /path/to/rivasand_hercules_branch_reversal_research/verify.py
```

This includes branch checks, historical parity, corrected native replays, and
the full 236-column corrected-oracle comparison. A verifier failure makes the
suite fail. The third argument is optional for environments that do not have
the private corrected oracle; when omitted, the suite explicitly reports that
corrected-default oracle acceptance was skipped.
The corrected Python oracle is separately retained in the private
`rivasand_hercules_branch_reversal_research` workspace; the JSON record identifies
its source hash and six-case comparisons. It does not become a runtime dependency.

Generate six CUDA inputs without claiming execution:

```sh
python3 benchmark_inputs/nonlinear_material_point/run_riva_branch_reversal_cuda.py \
  --output standalone/branch_cuda
```

On an NVIDIA/MPI host, rebuild psolve from this branch, then run:

```sh
python3 benchmark_inputs/nonlinear_material_point/run_riva_branch_reversal_cuda.py \
  --psolve /path/to/psolve --launcher 'mpirun -np 1' \
  --output standalone/branch_cuda
bash benchmark_inputs/nonlinear_material_point/run_suite.sh /path/to/psolve
```

The runner requires four reversals at every resolution, finite complete outputs,
successful full-state CPU/GPU comparison inside psolve, and exact Standard02/
Custom02 output equality separately on each backend. Fixture rows also pass
their event-count checks in the local host kernel.

**Pending:** CUDA compilation/execution, device checkpoint file I/O, multi-rank
MPI restart, both geostatic field workflows, coupled wave/BrickUP time-step and
mesh checks, and experimental validation. The development Mac has no CUDA/MPI
toolchain. The OpenSees column time step is not a Hercules stability qualification.
