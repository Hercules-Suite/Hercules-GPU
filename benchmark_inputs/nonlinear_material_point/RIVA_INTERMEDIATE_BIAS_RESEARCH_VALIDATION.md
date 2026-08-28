# RIVA-Sand intermediate-bias research validation

This record applies only to `RIVASandIntermediateBiasResearch`, translated
from OpenSees_RIVA commit
`637d207e867b863dee7b2bf14d9aeb75fdc32e7b`. It does not qualify or modify
production `RIVASand` or `RIVASandCustom`.

Native contract:

- 248 frozen calibrated inputs; SHA-256
  `85bd23d467e3e5f29a3b44da9903606a9dee95d5f16d9f25f16e1b1611b34397`.
- 62 logical point-state fields, flattened to 137 schema values.
- 1,104 bytes per native point state and 1,976 bytes in the shared immutable
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
- [x] Native 235-column replay driver reading `input_paths` only.
- [x] Exact byte restart at static preload, immediately post-activation, and
  reversal-adjacent checkpoints.
- [x] Activation idempotence, point-state isolation, eight-point homogeneous
  update, rigid-motion suppression, and BrickUP ownership host checks.
- [x] Production frozen/custom host-state regression remains passing.
- [x] Documentation build and all 28 rendered-source tests pass.
- [x] All six cases pass core parity at `rtol=1e-8`, `atol=1e-10`.
- [x] All six cases pass full 235-column parity at the same tolerance. The
  native research kernel preserves NumPy's pairwise tensor-contraction tree,
  binary64 operation boundaries, and the frozen NumPy 2.0.2 FMA `tanh`
  intervals needed by the dense unbiased phase path. No tolerance was
  loosened and no calibration value was changed.
- [ ] CUDA compile, host/device parity, device checkpoint round trip,
  multi-rank MPI restart, approximate/ramped geostatic wave cases, and complete
  BrickUP motions require an NVIDIA/MPI Hercules host. This local machine has
  no `nvcc`, CUDA headers, `mpirun`, or `mpicxx`.

Run the complete local suite with:

```bash
benchmark_inputs/nonlinear_material_point/run_riva_intermediate_bias_research_suite.sh \
  /Users/wzhang/software/PJ_model/rivasand_intermediate_bias_flow_port \
  /opt/homebrew/bin/g++-15
```

The suite runs the state checks and native benchmark even after a parity
failure so all validation evidence is retained in one invocation.

## Native host benchmark

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
