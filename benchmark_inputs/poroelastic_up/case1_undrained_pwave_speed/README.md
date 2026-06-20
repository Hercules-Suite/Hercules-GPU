# Case 1 — Undrained P-wave speed (Gassmann limit)

**What this validates:** that the BrickUP u-p coupling produces the correct
*undrained* compressional stiffness. At `permeability = 0` a poroelastic element
is exactly an elastic element with the undrained modulus, so the propagating
P-wave speed must be

```
Vp_u = sqrt(Vp_drained^2 + alpha^2 * M / rho),   M = K_f / porosity
     = sqrt(400^2 + (2.2e9/0.4)/2000) = 1705.9 m/s
```

**Mesh:** `simulation_wave_max_freq_hz = 25` with `Vs_min = 200` and 8 nodes/λ
gives **1 m elements** (`200/(25·8) = 1.0`). The source is the `fc = 10 Hz`
Ricker base-velocity pulse (energy within the 25 Hz mesh band).

**Test regime matters.** The column must be several *undrained* wavelengths long
or the wave crosses it quasi-statically and the speed can't be measured. At
`fc = 10 Hz`, `lambda_u = Vp_u/fc ≈ 171 m`, so the **1024 m** column ≈ **6 λ_u**.
`dt` is sized for `Vp_u` (not the dry `Vp`).

**Three runs bracket the answer:**
- `./` (poroelastic, perm=0)
- `elastic_twin/` — Vp = `Vp_u` = 1705.872 (undrained-equivalent)
- `drained_twin/` — Vp = 400 (dry skeleton; coupling absent)

## Run

```bash
mpirun -np <N> <psolve> parameters.in                                # poroelastic
cd elastic_twin && mpirun -np <N> <psolve> parameters.in && cd ..     # undrained ref
cd drained_twin && mpirun -np <N> <psolve> parameters.in && cd ..     # drained ref
python3 post/case1_pwave_speed.py
```

## Pass criteria
- Measured P-wave speed (first-arrival slope fit) within **3 %** of 1705.9 m/s.
- Mean NRMS(vz) of poroelastic vs the **undrained** twin **< 5 %**.
- The script's verdict reports which medium the poroelastic run matches:
  - **undrained twin** → coupling delivers `alpha^2*M` correctly (PASS);
  - **drained twin** → coupling not adding the undrained stiffening (FAIL — solver issue);
  - **neither** → partial coupling / dispersion / setup — investigate.

## Notes / debugging
- If the speed comes out near 400 (or the run matches the drained twin), the
  undrained coupling isn't engaging — check `include_poroelastic`, the element
  classification (monitor log poroelastic element/node counts), and `perm=0`.
- If the measured speed is nonsensical and *decreasing with distance*, the column
  is too short / source too low-frequency for the wave (quasi-static regime) —
  lengthen the column or raise the source frequency.
- Eyeball a trace (e.g. `vz` at z=0) for poro vs each twin to see directly whether
  it's a speed difference, amplitude difference, or phase shift.
