# Poroelastic BrickUP (u-p) validation suite

Validation cases for the reduced-Biot `u-p` (BrickUP-style) poroelastic element,
against **analytical references** (no second code needed). Hercules' u-p formulation
neglects fluid inertia / the slow P-wave, treats the input `Vp/Vs` as the **drained
skeleton**, and uses `M = K_f/porosity`, `alpha=1`. The key fact that makes exact
analytical references possible:

> In the undrained limit (`permeability = 0`), a poroelastic element is *exactly* an
> elastic element with undrained moduli: `Vp_u = sqrt(Vp_drained^2 + alpha^2*M/rho)`,
> `Vs` and `rho` unchanged.

| Case | Folder | Reference | Status |
|------|--------|-----------|--------|
| 1 | `case1_undrained_pwave_speed/` | Gassmann `Vp_u` + undrained & drained twins | active — P-speed slope + which-twin verdict |
| 2 | `case2_sh_transfer_function/` | Kausel/Roesset 1D SH transfer function | **DEFERRED** — base-force thin column is killed by absorbing side boundaries; needs the plane-wave/DRM injector in a wider domain (see its README) |
| 3 | `case3_pointsource_saturated_layers/` | undrained-equivalent elastic twin | active — per-offset surface-motion NRMS < 3 % |

Each case is self-contained (`parameters.in`, `sourcefiles/`, an `elastic_twin/` or
`elastic/` reference run, a `post/*.py` analysis script, and a `README.md`).
`homogeneous_column_base_force/` is the original smoke-test column the cases are built from.

## Cross-cutting conventions (every case)
- `poroelastic_permeability_* = 0`, `poroelastic_biot_alpha = 1`,
  `poroelastic_surface_pressure_bc = undrained` (so the elastic equivalence is exact).
- **Size `simulation_delta_time_sec` from the undrained `Vp`**, not the dry `Vp` — the
  solver's auto-CFL only sees the dry `Vp` and would under-size the step.
- Analysis scripts require `numpy` (run on the host with the station output).
- **Octree domain rule:** divide (`region_length_north_m`, `region_length_east_m`,
  `region_depth_deep_m`) by their integer GCD; the largest reduced value must be a power
  of 2 (octor.c embeds the domain in a `2^31`-tick cube along that axis). E.g. `1×1×64` ✓
  (`60` ✗ → element sizes don't tile the 1 m lateral); `2000×2000×1500` → `4×4×3` ✓.
- **Vp/Vs cap & the elastic twin:** `mu_and_lambda` clamps `Vp` to `the_threshold_Vp_over_Vs · Vs`
  (default 3.0). The undrained-equivalent **elastic twin** uses the *undrained* `Vp` with the
  *low* skeleton `Vs`, so its `Vp/Vs` is high (e.g. saturated soft soil `1804/300 = 6`) and gets
  clamped unless you raise `the_threshold_Vp_over_Vs` above it — otherwise the twin silently runs
  at `3·Vs` and is not the undrained-equivalent. The poroelastic run does **not** need this (its
  drained `Vp/Vs` is low, and the Biot coupling adds `α²M` as a force that bypasses the cap) — so
  the coupling can represent an undrained `Vp/Vs` the elastic path cannot. Always set the twin's
  threshold above the saturated layers' undrained `Vp/Vs`.

## Run pattern
```bash
cd caseN_xxx && mpirun -np <N> <path-to-psolve> parameters.in
cd elastic_twin && mpirun -np <N> <path-to-psolve> parameters.in && cd ..   # reference
python3 post/<script>.py
```
