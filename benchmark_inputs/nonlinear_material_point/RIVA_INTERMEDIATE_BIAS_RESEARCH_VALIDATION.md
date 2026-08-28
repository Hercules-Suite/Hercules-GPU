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
- [ ] Six-case core parity at `rtol=1e-8`, `atol=1e-10`: five pass. The
  `3484_zero_bias` trace first exceeds tolerance late in the low-confinement
  history (`state_alpha` at output row 99; stress at row 131). The largest
  reported stress difference is approximately `1.11e-7` in the fixture's kPa
  units. Events, substep counts, and gates remain aligned.
- [ ] Six-case full parity at the same tolerance: four pass. In addition to
  the zero-bias drift, `4666_dense_zero_bias` has a near-zero normal-component
  difference of approximately `1.73e-10` at output row 291. No tolerance was
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

The suite deliberately returns nonzero until both core and full six-case
parity pass. It still runs the state checks and native benchmark after a parity
failure so the remaining evidence is retained.

## Native host benchmark

The benchmark replays the 484-row `3484_intermediate_bias025` path with the
same initial state, activation row, row-specific fixed substeps, compiler, and
process for both kernels. A 2,000-repetition GCC 15 `-O2` sample measured:

| Kernel | Time | Time per row |
| --- | ---: | ---: |
| Production RIVA-Sand | 2.830349 s | 2,923.91 ns |
| Intermediate-bias research | 4.678666 s | 4,833.33 ns |

The research/production ratio was 1.653 for that local host sample. A repeated
2,000-repetition run measured 2.844910 s and 4.358742 s (1.532), respectively;
the spread is reported rather than hidden because this laptop benchmark does
not pin cores or control system load. This is a native constitutive
microbenchmark, not the required GPU/MPI production benchmark; target-system
register pressure, occupancy, transfer, and end-to-end wave runtime remain
unmeasured.
