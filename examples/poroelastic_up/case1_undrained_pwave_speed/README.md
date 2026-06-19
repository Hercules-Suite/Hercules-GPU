# Case 1 — Undrained P-wave speed (Gassmann limit)

**What this validates:** that the BrickUP u-p coupling produces the correct
*undrained* compressional stiffness. In the undrained limit (`permeability = 0`)
a poroelastic element is exactly an elastic element with the undrained modulus,
so the propagating P-wave speed must be

```
Vp_u = sqrt(Vp_drained^2 + alpha^2 * M / rho),   M = K_f / porosity
     = sqrt(400^2 + 1.0^2 * (2.2e9/0.4) / 2000) = 1705.9 m/s
```

The dry skeleton `Vp` is only 400 m/s, so the wave you actually see is ~4.3×
faster — this case also confirms you sized the time step for `Vp_u` (the solver's
auto-CFL uses the *dry* `Vp` and would under-size it).

**Setup:** homogeneous 1×1×32 m saturated column, `perm = 0`,
`poroelastic_surface_pressure_bc = undrained`. A **vertical** compliant-base
force (`sourcefiles/source.in`, reusing the 10 Hz Ricker velocity pulse) launches
a P pulse up from the base plane at z = 32 m. Nine receivers run down the column.

## Run

```bash
# poroelastic run
mpirun -np <N> <path-to-psolve> parameters.in
# undrained-equivalent elastic twin (Vp = Vp_u, poroelasticity off)
cd elastic_twin && mpirun -np <N> <path-to-psolve> parameters.in && cd ..
# analyze
python3 post/case1_pwave_speed.py
```

## Pass criteria
- Measured P-wave speed (first-arrival slope fit) within **2 %** of 1705.9 m/s.
- NRMS of `vz` between the poroelastic run and the elastic twin **< 2 %** at every station.
- Vertical compression produces excess pore pressure of consistent sign (diagnostic).

## Notes
- `simulation_delta_time_sec = 5e-5` gives Courant ≈ 0.17 at `Vp_u` (safe). Do **not**
  raise it toward the dry-`Vp` CFL.
- The elastic twin reuses the same source via `source_directory = ../sourcefiles`; run it
  from inside `elastic_twin/`.
