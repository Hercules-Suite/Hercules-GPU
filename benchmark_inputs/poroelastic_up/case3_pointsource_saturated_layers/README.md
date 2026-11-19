# Case 3 — Point source under saturated soil layers (realistic target)

**What this validates:** the realistic configuration — a point earthquake source in
elastic bedrock radiating up through shallow **saturated poroelastic** layers to
surface receivers — gives the physically correct ground motion. The analytical
handle is the **undrained equivalence**: with `permeability = 0`, each saturated
poroelastic layer is *exactly* an elastic layer with the undrained `Vp`
(`Vp_u = sqrt(Vp_drained^2 + alpha^2*M/rho)`; `Vs`, `rho` unchanged). So the
poroelastic run is checked against an elastic twin built from those undrained moduli.

**Model:** 2×2×1.5 km box; 3-layer profile — saturated soil A (0–150 m), saturated
soil B (150–400 m), elastic bedrock (≥400 m). `poroelastic_depth_max_m = 400` makes
only the top two layers poroelastic. Point double-couple at 800 m depth (in the
bedrock). Surface receivers at 0–800 m east offset.

Undrained-equivalent `Vp` (used in `elastic_twin/`):
`layer A 1804 m/s`, `layer B 1936 m/s`.

## Run

```bash
mpirun -np <N> <psolve> parameters.in                                   # poroelastic
cd elastic_twin && mpirun -np <N> <psolve> parameters.in && cd ..       # equivalent elastic
python3 post/case3_compare.py
```

## Pass criteria
- **Primary:** the cross-correlation **lag** of poro vs the undrained twin (body-wave window,
  `post/case3_compare.py`) is **sub-ms** on reliable (high-correlation) components. The lag is
  amplitude-independent and sub-sample, so it isn't fooled by the ~10 % peak-amplitude difference
  that biases a relative-threshold first-break pick. Kernel analysis predicts ~0.1 ms; a correct
  coupling gives a sub-ms lag (not the +3..11 ms a threshold picker reports). Run `--selftest` to
  verify the estimator (sign, sub-sample accuracy, amplitude-independence).
- The body-wave NRMS is kept for context; the coda diverges from staggered-scheme dispersion
  (expected), so NRMS over the full trace is not the right pass metric.
- Pore-pressure traces are reported as a diagnostic (poroelastic run only).

## Gotchas (important)
- **`dt` is sized for the fast undrained `Vp`** (≥1936 m/s in the soil, 2500 in bedrock),
  not the dry `Vp`. The solver's auto-CFL uses the dry `Vp` and would under-size it.
- **Keep the source below `poroelastic_depth_max_m`** (here 800 m vs 400 m) so it sits
  in genuinely elastic material — a source inside a poroelastic element breaks the clean
  elastic-equivalent reference.
- **`permeability = 0`** is required for the equivalence to be exact; any drainage adds
  diffusion absent from the elastic twin.
- **`poroelastic_surface_pressure_bc = undrained`** so no spurious `p=0` plane is imposed
  that the elastic twin cannot reproduce.
- **`poroelastic_undrained_fast_path = yes`** (default) folds `alpha^2 M` into the native
  elastic volumetric stiffness and leaves pore pressure as a diagnostic; this is the path
  intended to match the undrained-equivalent elastic twin.
- Active sealed poroelastic cells now use the undrained normal dashpot
  impedance in both fast and explicit pressure paths, with unchanged shear
  impedance. This case meets the zero-permeability, undrained-surface and
  no-open-pore eligibility conditions. An exact undrained elastic twin therefore
  has matching mechanical dashpots as well as fast-path stiffness, up to
  material rounding. Historical results at `2fd4f81` used drained porous
  coefficients and had different outer reflections. Local dashpots still
  reflect oblique waves and have no Darcy absorber; twin agreement is a
  shared-code check, not independent boundary validation. In the measured
  fine-grid follow-up, the full-record fast-path difference falls from 1.407%
  to 0.119%; the remaining difference is not declared solved. See the
  [measured follow-up](../WOLF_2022_VALIDATION.md#sealed-dashpot-correction-follow-up).
- The optional CPS reflectivity cross-check needs the equivalent-elastic model
  and a domain/record for which artificial outer-boundary effects are excluded.

This case is a template — adjust the profile, source, offsets, and domain size to your
actual site; keep the gotchas above.
