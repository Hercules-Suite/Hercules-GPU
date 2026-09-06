# RIVA-Sand intermediate-bias research validation

This record applies to `RIVASand02` and `RIVASandCustom02`, translated
from OpenSees_RIVA commit
`637d207e867b863dee7b2bf14d9aeb75fdc32e7b`, with the plastic-activity gate
from `5e4ebd30cb5f43b677c61306c7f1fbbdcb24905c`. It does not qualify or modify
production `RIVASand` or `RIVASandCustom`.

Native contract:

- 249 frozen calibrated inputs; SHA-256
  `3c17e962e64d3a3b29d4797399380dfc85669a5b9fd20b187af359f57aa5b56e`.
- 63 logical point-state fields, flattened to 138 schema values.
- Kernel revision 3; revision-2 research checkpoints are rejected by the
  existing state-size, kernel-revision, and parameter-hash checks.
- 1,112 bytes per native point state and 1,984 bytes in the shared immutable
  parameter block on the arm64 GCC 15 host used below.
- Effective skeleton stress only. BrickUP retains pore-pressure, storage,
  drainage, and total-stress assembly ownership.

## Local validation status

- [x] Private package manifest, schemas, oracle tests, and all six golden
  fixtures verified.
- [x] Separate parser/factory identity and `Vs_limit Dr` input boundary.
- [x] Native host/CUDA-compatible kernel with fixed substeps, host reversal
  order, ratchet reevaluation, compatible-volume rebuild, mapping correction,
  and activation caches.
- [x] Native 236-column replay driver reading `input_paths` only, including
  the new `state_ep_half_last` column.
- [x] Exact byte restart at static preload, immediately post-activation, and
  reversal-adjacent checkpoints.
- [x] Activation idempotence, point-state isolation, eight-point homogeneous
  update, rigid-motion suppression, and BrickUP ownership host checks.
- [x] Production frozen/custom host-state regression remains passing.
- [x] Gate endpoints and transition, independent PT volume, weak microcycles,
  and completed-half-cycle capture in both real backbones pass with 1, 4,
  and 16 fixed substeps. Capture occurs exactly once per host reversal.
- [x] The new scalar survives byte restart after a completed half-cycle and
  idempotent activation; reused point initialization clears it.
- [x] Documentation lint, TypeScript checks, build, and all 29 rendered/source
  tests pass after the 02 interface update.
- [x] All six cases pass core parity at `rtol=1e-8`, `atol=1e-10`.
- [x] All six cases pass full 236-column parity at the same tolerance. The
  native research kernel preserves NumPy's pairwise tensor-contraction tree,
  binary64 operation boundaries, and the frozen NumPy 2.0.2 FMA `tanh`
  intervals needed by the dense unbiased phase path. No tolerance was
  loosened and no calibration value was changed.
- [x] The standard Hercules material-point runner now has a dedicated
  intermediate-bias research fixture. It replays the same native core on host
  and CUDA, compares all 138 state values and lifecycle selectors, checks an
  exact zero increment, and invokes the uniform sealed BrickUP pressure path.
- [ ] Target-host CUDA compilation and execution of that standard fixture,
  device checkpoint round trip,
  multi-rank MPI restart, approximate/ramped geostatic wave cases, and complete
  BrickUP motions require an NVIDIA/MPI Hercules host. This local machine has
  no `nvcc`, CUDA headers, `mpirun`, or `mpicxx`.

Run the standard host/CUDA fixture on the NVIDIA target with:

```bash
mpirun -np 1 quake/forward/psolve --nonlinear-material-point \
  benchmark_inputs/nonlinear_material_point/riva_sand_intermediate_bias_research_synthetic.in
```

Run the complete local suite with:

```bash
benchmark_inputs/nonlinear_material_point/run_riva_intermediate_bias_research_suite.sh \
  /Users/wzhang/software/PJ_model/rivasand_intermediate_bias_flow_port \
  /opt/homebrew/bin/g++-15
```

The suite runs the state checks and native benchmark even after a parity
failure so all validation evidence is retained in one invocation.

## Plastic-activity update (2026-09-06)

The new scalar stores the plastic multiplier accumulated in the last
completed half-cycle. Both the base and directional-mapping backbones capture
it before their reversal counters reset; later substeps retain that captured
value. Its bounded activity factor acts only on the inherited biased-volume
target, not on the separate phase-transformation volume. Initialization clears
the scalar, while reactivating an already-active restored point preserves it.
No fictitious cyclic reversals are introduced at activation.

Use the private revision-3 handoff package for oracle verification. No golden
responses or tolerances were changed to make the native implementation pass.
The new state is appended after the previous 137 numerical values, and the
host/CUDA result buffers now include it. Full GPU/MPI checkpoint execution
still requires the target system; the local byte-restart test is not a
substitute for that qualification.

The gate originally retained the `RIVASandIntermediateBiasResearch` name.
The subsequent interface change below retires that name. Neither change
imports the later OpenSees `BiasVolume` modes or OpenSees-specific
`reversalLatch`. Original `RIVASand`/`RIVASandCustom` are unchanged.

## 02 interface rename and configurable adapter (2026-09-06)

- `RIVASand02`: renamed default interface with `Vs_limit Dr` rows.
- `RIVASandCustom02`: the same revision-3 kernel with `Vs_limit Dr M kd h m
  zeta e_max e_min Q R n_G` rows. Custom material blocks are interpolated and
  stored per element on both host and device; default 02 retains one block
  per rank. Material-point inputs omit only the Vs selector.
- `RIVASandIntermediateBiasResearch` is no longer accepted by either parser.
  Historical source/fixture filenames remain unchanged for traceability.
- No constitutive code or state-layout changes. Default model ID remains
  15, new Custom02 is 16, and revision-3 reference hash/state remain unchanged.
  Resume with the same model selection and physical inputs. This does not
  make revision-2 checkpoints compatible.
- Native exact default/custom reference equality, all ten independently
  perturbed physical inputs, per-material isolation, void-ratio bounds, and
  custom byte restart pass with 1/4/16 substeps. This tests input activity,
  not calibration of the perturbed parameter sets.
- All six successor histories still pass full 236-column parity at the
  established `rtol=1e-8`, `atol=1e-10`. The original five histories retain
  their pre-existing host tolerance; no golden data is regenerated.
- The standard target-host suite adds `riva_sand_custom02_equivalence.in`
  and requires byte-identical default/custom full CSVs independently on CPU
  and CUDA. Source-wiring tests cover both parser names, column layouts,
  per-element GPU parameter indexing/allocation, launch selection, geostatic
  reporting, and checkpoint state-family selection.
- Full NVIDIA/MPI compilation and execution are still pending on the target
  machine. Host tests and source inspection do not substitute for those runs.

### Incremental host cost

Three alternating before/after samples used GCC 15.1.0 `-O2`, 1,000
repetitions of the same 484-increment `3484_intermediate_bias025` path, and
the same local arm64 machine. Median research-kernel time was 2.983728 s
before the gate and 2.975631 s after it (ratio 0.9973). This difference is
within ordinary timing noise and does not establish a speedup or a GPU cost.
The gate adds one 8-byte state scalar and one 8-byte shared parameter;
it introduces no additional constitutive iterations or substeps. The
research kernel remains approximately twice the cost of the original model
in this host microbenchmark; that is not overhead introduced by this gate.

## Historical revision-2 native host benchmark

The benchmark replays the 484-row `3484_intermediate_bias025` path with the
same initial state, activation row, row-specific fixed substeps, compiler, and
process for both kernels. A 2,000-repetition GCC 15 `-O2` sample measured:

| Kernel | Time | Time per row |
| --- | ---: | ---: |
| Production RIVA-Sand | 2.727994 s | 2,818.18 ns |
| Intermediate-bias research | 5.788280 s | 5,979.63 ns |

The research/production ratio was 2.122 for that local host sample. The strict
oracle arithmetic adds host-side cost to a research kernel that also carries a
larger state and parameter block. This laptop benchmark does not pin cores or
control system load. It is a native constitutive microbenchmark, not the
required GPU/MPI production benchmark; target-system register pressure,
occupancy, transfer, and end-to-end wave runtime remain unmeasured.
