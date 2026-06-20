# Case 2 — SH site-response transfer function  ⚠️ DEFERRED

> **STATUS: deferred — does not work as built.** The base-force-in-a-thin-column
> approach below is killed by the absorbing (dashpot) **side** boundaries, which
> damp the laterally-uniform shear motion so the column never resonates. This was
> confirmed with the diagnostic `post/case2_sh_tf.py`: the **elastic** run also
> shows no amplification (`surf/within RMS ≈ 1`, no peak near `f0`), so it is a
> 1D-column test limitation, **not** a poroelastic-solver problem.
>
> **Correct fix (TODO):** drive a vertically-incident SV plane wave via the
> `drm_planewaves` injector (`include_incident_planewaves = yes`,
> `type_of_wave = SV`, `planewave_Z_angle = 0`) inside a **wider** domain so the
> DRM box fits with margins and the receivers sit in the central column; the side
> dashpots then only see the residual/scattered field. A wider domain alone (still
> base-forced) will not help a thin column. The `post/case2_sh_tf.py` diagnostic
> (reports the TF for both runs + a verdict) can be reused once that is built.
>
> Shear fluid-independence is meanwhile exercised indirectly by Case 3 (poroelastic
> vs undrained-equivalent elastic twin includes the shear response).

The original design (kept for reference / the eventual plane-wave rebuild):

**What this validates:** (1) the poroelastic implementation does **not** corrupt the
shear/site response — shear is fluid-independent, so the poroelastic and elastic
runs must give the same SH spectral ratio (the firm correctness check); and (2) the
surface site amplification matches the analytical 1D SH transfer function for a soft
layer over an elastic halfspace.

**Profile:** soft soil (0–20 m: `Vs=150`, `rho=1900`) over bedrock (`Vs=450`,
`rho=2200`), 1×1×64 m column. Fundamental site frequency `f0 = Vs/(4H) = 1.875 Hz`;
impedance ratio `alpha_z = (rho1*Vs1)/(rho2*Vs2) = 0.288`; peak amplification
`~1/alpha_z ~ 3.5`.

**Input:** a horizontal compliant-base force (`sourcefiles/source.in`) drives a
vertically-propagating shear wave from the base; the STF is a broadband (fc=4 Hz)
Ricker velocity pulse (`sourcefiles/baseforcefunction.in`). The incident-plane-wave
injector (`include_incident_planewaves=yes`, `type_of_wave=SV`, `planewave_Z_angle=0`)
is the documented alternative, but the base force is simpler/robust in a thin column.

## Run

```bash
cd poroelastic && mpirun -np <N> <psolve> parameters.in && cd ..
cd elastic     && mpirun -np <N> <psolve> parameters.in && cd ..
python3 post/case2_sh_tf.py
```

## Pass criteria
- **Resonance frequency**: the surface TF peak is within **5 %** of `f0 = 1.875 Hz`
  (and overtones near 5.6, 9.4 Hz). This is the firm physics check.
- **Fluid independence**: poroelastic vs elastic surface `vx` NRMS **< 3 %**.
- **Peak amplitude** ≈ `1/alpha_z ≈ 3.5` (approximate, within ~20 %). The absolute
  value depends on the dashpot incident/outcrop convention (a factor ~2) and on
  damping (here `Q≈1e6`, ~none), so treat amplitude as a sanity check, not a tight gate.

## Notes
- `simulation_delta_time_sec = 1e-4` is sized against the undrained `Vp` of the
  poroelastic run (~1.8 km/s), not the dry `Vp`.
- Stations: surface (z=0) and within bedrock (z=25 m). The TF is referenced to the
  known input velocity history.
