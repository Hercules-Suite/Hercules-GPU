# RIVA-Sand Hercules integration status

This repository-local copy records the status of the checklist supplied with
the frozen RIVA-Sand port package. The authoritative constitutive parameter SHA-256 is
`9585f0c155c9444885c5115e7753a1a5d97c783e2c685938a49a65435c9e8f83`.

## Completed locally

- [x] Registered `RIVASand` and documented its `Vs_limit relative_density`
  production table; the input boundary converts dimensionless `D_r` in
  `[0,1]` to internal void ratio using frozen `e_max=0.78`, `e_min=0.51`.
- [x] Registered `RIVASandCustom` as the same RIVA-Sand kernel with model-specific,
  compact per-element material storage and the 12-column production table
  `Vs_limit Dr M kd h m zeta e_max e_min Q R n_G`. The original two-column
  `RIVASand` interface and frozen response remain unchanged; `h_p0` is not
  included because it belongs to the separate pre-name V8.1 experiment. Frozen RIVA-Sand
  shares one material block per rank; unrelated nonlinear families no longer
  reserve RIVA-Sand material or void-ratio fields in `nlconstants_t`.
- [x] Mapped `[xx, yy, zz, xy, yz, xz]`, tensor shear, tension-positive stress,
  and the factor-two symmetric double contraction explicitly.
- [x] Translated all 112 frozen parameters, four host-event controls, and all
  39 logical state fields (93 frozen flattened scalar values plus the separate
  committed geostatic-admission selector) to allocation-free
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
  exact zero-increment preservation in a host-only contract test. The same
  test verifies frozen/custom reference equivalence and confirms that exposed
  `zeta` changes the constitutive state.
- [x] Added a separate CUDA production object and per-integration-point state.
- [x] Added versioned checkpoint fields for prior strain/stress, complete RIVA-Sand
  history, and both BrickUP pressure time levels. Legacy displacement-only
  checkpoints are rejected for RIVA-Sand/BrickUP restarts instead of resetting state.
  The extension also records `RIVA_KERNEL_REVISION` separately from the frozen
  parameter SHA, so a pre-revision RIVA-Sand checkpoint cannot silently resume
  under changed constitutive equations; BrickUP-only checkpoint compatibility
  is unchanged.
- [x] Kept Python out of the production executable.
- [x] Implemented the BrickUP effective-stress split. RIVA-Sand supplies effective
  skeleton stress; BrickUP independently solves
  `d(zeta)=alpha*d(eps_v)+d(p_f)/M` and equilibrium applies
  `sigma=sigma'-alpha*p_f*I`. RIVA-Sand's standalone effective-pressure-loss
  diagnostic is never injected into `p_f`.
- [x] Added a prescribed-strain RIVA-Sand/BrickUP regression that compares all 93
  frozen numerical host/CUDA state values, separately compares the committed
  `geostatic_admitted` flag, and invokes the production uniform pressure kernels.
- [x] Bounded the inverse-square reversible-volume law below one quarter of its
  reference projected bias, synchronized the Python oracle and formulation,
  and added a low-but-nonzero static-bias host/CUDA material-point fixture. The
  gate tends continuously to zero and is exactly one throughout the retained
  calibrated range.
- [x] Added a common compressive geostatic initializer for both approximate-K0
  and ramped-gravity effective stresses. Admissible states remain unchanged;
  over-bound states preserve their equilibrated stress, recenter the mapping
  origins, and use a closed-form non-expansive radial limit with no added
  iterations or substeps. Noncompressive physical states are rejected.
- [x] Specified the Hercules-only equations in
  `RIVA_SAND_GEOSTATIC_ADMISSION.md`. Production initializes all active points
  at the first post-geostatic step and reports the global initialized/admitted
  counts plus maximum admitted cone ratio instead of silently discarding the
  admission result.
- [x] Added an internal residual-cone-pressure control with production default
  zero. The standalone contract test enables `0.01*p_atm` only to verify the
  translated cone and the separate physical-pressure admission check.
- [x] Added an over-bound ramped-geostatic host/CUDA material-point fixture;
  the existing K0/BrickUP fixture covers the admissible approximate-geostatic
  path through the same initializer. The committed admission flag is compared
  independently of the frozen 93-value numerical oracle.
- [x] Added a custom-interface material-point fixture whose reference-valued
  CPU and CUDA histories must exactly reproduce the frozen RIVA-Sand fixture.

## Maintainer/Horizon verification still required

- [ ] Build the full CUDA/MPI solver and run `run_suite.sh` (including the RIVA-Sand
  fixture) on the supported NVIDIA toolchain.
- [ ] Run a full homogeneous one-element solver case to verify element strain
  mapping, rigid motion, force assembly, and non-shared integration-point
  history outside standalone material-point mode.
- [ ] Restart that solver case at an intermediate step and compare it with the
  uninterrupted history.
- [ ] Run zero-bias, statically biased, and dense wave cases; inspect pressure,
  denominator-floor, beta-fallback, and nonfinite diagnostics.
- [ ] Run one full `approximate_geostatic_state=yes` RIVA-Sand/BrickUP case and
  one `approximate_geostatic_state=no` ramped-gravity case on the supported GPU
  solver; inspect the printed admission counts and maximum admitted cone ratio, then
  confirm first-dynamic-step stress continuity and finite station histories.
- [ ] Characterize optional `p_residual=0.01*p_atm` at shallow confinement
  before exposing or enabling it in production input. The default remains
  exactly zero.
- [ ] Characterize the reversible-volume branch at low cyclic amplitude and low
  effective confinement. The contained low-bias gate removes the singularity
  but does not validate or recalibrate those two separate regimes.
- [ ] Compare global timestep refinement and `no_substeps=1,2,4`; retain one
  unless subdivision materially changes the response.
- [ ] Run a stress-controlled coupled BrickUP element characterization with
  storage, permeability and boundary conditions representative of production.
  The current sealed uniform fixture validates field ownership and arithmetic,
  not the full boundary-value liquefaction response.
- [ ] Benchmark memory and update time against MD and PDMY03 under identical
  mesh, timestep, hardware, compiler, and output settings.
- [ ] Review those results before describing RIVA-Sand as production-validated.
