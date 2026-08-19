# PJ-Liq V8 Hercules integration status

This repository-local copy records the status of the checklist supplied with
the frozen V8 port package. The authoritative constitutive parameter SHA-256 is
`9585f0c155c9444885c5115e7753a1a5d97c783e2c685938a49a65435c9e8f83`.

## Completed locally

- [x] Registered `PJLiqV8` and documented its `Vs_limit relative_density`
  production table; the input boundary converts dimensionless `D_r` in
  `[0,1]` to internal void ratio using frozen `e_max=0.78`, `e_min=0.51`.
- [x] Mapped `[xx, yy, zz, xy, yz, xz]`, tensor shear, tension-positive stress,
  and the factor-two symmetric double contraction explicitly.
- [x] Translated all 112 frozen parameters, four host-event controls, and all
  38 logical state fields (93 flattened scalar values) to allocation-free
  binary64 host/device code.
- [x] Initialized from equilibrated effective stress and activated dynamics
  once with the zero-static-shear reference stress.
- [x] Implemented fixed Forward-Euler subdivision with one host-level reversal
  decision and default `no_substeps=1`.
- [x] Replayed all five golden histories with a native C++ driver. Core fields
  pass at `rtol=1e-6, atol=1e-8`; all 93 values pass at
  `rtol=1e-6, atol=1e-6`. Event flags and integer counters remain exact. Four
  histories also pass the package's tighter defaults; the final near-collapse
  part of the zero-bias history is sensitive to platform `libm` arithmetic.
- [x] Verified per-point state isolation, byte-copy restart reproduction, and
  exact zero-increment preservation in a host-only contract test.
- [x] Added a separate CUDA production object and per-integration-point state.
- [x] Added versioned checkpoint fields for prior strain/stress, complete V8
  history, and both BrickUP pressure time levels. Legacy displacement-only
  checkpoints are rejected for PJ/BrickUP restarts instead of resetting state.
- [x] Kept Python out of the production executable.
- [x] Implemented the BrickUP effective-stress split. V8 supplies effective
  skeleton stress; BrickUP independently solves
  `d(zeta)=alpha*d(eps_v)+d(p_f)/M` and equilibrium applies
  `sigma=sigma'-alpha*p_f*I`. V8's standalone effective-pressure-loss
  diagnostic is never injected into `p_f`.
- [x] Added a prescribed-strain PJ/BrickUP regression that compares all 93
  host/CUDA state values and invokes the production uniform pressure kernels.

## Maintainer/Horizon verification still required

- [ ] Build the full CUDA/MPI solver and run `run_suite.sh` (including the PJ
  fixture) on the supported NVIDIA toolchain.
- [ ] Run a full homogeneous one-element solver case to verify element strain
  mapping, rigid motion, force assembly, and non-shared integration-point
  history outside standalone material-point mode.
- [ ] Restart that solver case at an intermediate step and compare it with the
  uninterrupted history.
- [ ] Run zero-bias, statically biased, and dense wave cases; inspect pressure,
  denominator-floor, beta-fallback, and nonfinite diagnostics.
- [ ] Compare global timestep refinement and `no_substeps=1,2,4`; retain one
  unless subdivision materially changes the response.
- [ ] Run a stress-controlled coupled BrickUP element characterization with
  storage, permeability and boundary conditions representative of production.
  The current sealed uniform fixture validates field ownership and arithmetic,
  not the full boundary-value liquefaction response.
- [ ] Benchmark memory and update time against MD and PDMY03 under identical
  mesh, timestep, hardware, compiler, and output settings.
- [ ] Review those results before describing PJ-Liq V8 as production-validated.
