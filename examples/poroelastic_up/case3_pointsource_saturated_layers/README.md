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
- Per-offset NRMS of surface velocity (`vx, vy, vz`) between the poroelastic run and
  the elastic twin **< 3 %** at every offset. (Exact at `perm=0` up to discretization.)
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
- The poro-vs-twin check is insensitive to absorbing-boundary reflections (shared by both
  runs). Scale the domain/record up only if you add the **optional** CPS reflectivity
  cross-check (build the equivalent-elastic `.mod` under `benchmark_inputs/`).

This case is a template — adjust the profile, source, offsets, and domain size to your
actual site; keep the gotchas above.
