# Topography example: boundary-reaching gentle tilt (flatten-PML test)

Purpose-built to **actually exercise the C-PML flatten-PML re-route**
(`cpml_pml_flattened_element`), which the earlier Gaussian-hill example
**cannot** test: that hill decays to ~0 over the PML band, so no PML
element is topography-intersected and the re-route never fires
(`flattened-pml-skipped ≈ 0`).

This example is a **gentle planar tilt that is non-zero and sloping
across the entire domain, including all four lateral PML bands**, so for
`computation_method=vt` essentially every near-surface lateral-PML
element is topography-deformed → the flatten-PML predicate returns 1 for
all of them → **`flattened-pml-skipped` is large and uniform**. That is
the functional test.

## Files

| File | Purpose |
|---|---|
| `gen_topography_tilt.py` | reproducible generator (parameterized) |
| `topography.in` | generated elevations (the only file Hercules needs) |

Regenerate / reparameterize:
```
python3 gen_topography_tilt.py                     # defaults below
python3 gen_topography_tilt.py --sx 0.02 --sy 0    # steeper, W-E flat
python3 gen_topography_tilt.py --bump-h 400        # + interior bump too
```

## The surface

Pure plane `elev(x,y) = z0 + sx·x + sy·y`, defaults `sx=sy=0.0125 m/m`,
`z0=0`, domain 40000×40000 m → SW corner 0 m, NE corner **1000 m**,
slope ≈ **1.01°** (deliberately gentle — we test the re-route
*mechanism*, not steep-boundary pathology). Format/ordering as the other
example (line 1 = spacing `100`; first point = SW corner; inner x
south→north, outer y west→east). Verified: SW=0, centre=500, NE (last
line)=1000; lateral-PML-band points are 17.5 … ~1000 m (clearly
non-flat in the PML, unlike the Gaussian).

## Validation recipe — isolating the flatten-PML re-route

**Critical:** use `type_of_damping = none`, NOT `bkt2`. VT + BKT is a
*separate, pre-existing* catastrophic instability (NaN ~9 s,
ABC-independent — see the spun-off VT+BKT root-cause task); running
`vt`+`bkt2` here would blow up from *that* and mask the flatten-PML
behavior. `type_of_damping=none` + `vt` + `cpml` is known stable for
flat PML (run `…topo_v4`), so any failure with this tilted topography is
attributable to the flatten-PML re-route specifically.

Run three cases (same 40 km elastic-halfspace / cps parameters; ~30–60 s
is plenty), all `type_of_damping=none`, `abc_type=cpml`:

| # | `computation_method` | topography | Tests |
|---|---|---|---|
| **A** | `vt` | this tilt | **the flatten-PML re-route** (predicate fires for all lateral-PML elements) |
| **B** | `FEM` | this tilt | FEM path with boundary topo (flatten-PML is a strict no-op for FEM by construction) — comparison wavefield |
| **C** | (`vt` or `FEM`) | none (flat) | reference baseline |

`topography_directory = <repo>/benchmark_inputs/topography_tilt_boundary`

### What to check

1. **Diagnostic engaged:** Run A startup log
   `Topography force path: gpu (... flattened-pml-skipped=N)` with **N
   large** (≈ all near-surface lateral-PML elements). Contrast: the
   Gaussian example gives N≈0. A large N here *proves the re-route is
   being exercised* — the precondition the Gaussian lacked.
2. **Stability (primary pass/fail for the re-route):** Run A finite over
   the full duration — no NaN/blow-up/growth. If the re-route were
   buggy (double-count, missing native force, count≠mapping), A would
   crash or diverge. (B/C are sanity baselines.)
3. **Transition-band reflection (the documented bounded tradeoff):**
   flatten-PML makes the lateral-PML elements regular flat cubes while
   the interior just inside is VT-tilted → a small geometric step at the
   PML inner edge. Quantify the spurious reflection: late-coda energy in
   A vs the flat reference C; expect it small/bounded (same class as the
   other discretization residuals). This *measures* the tradeoff we
   accepted for the SPECFEM-style regular PML.
4. **A vs B:** both have boundary topography; differences reflect the
   VT-flatten treatment vs FEM’s own boundary-topo handling.

This complements the *by-construction* byte-identity safety of the
predicate (already established: returns 0 unless VT+topo+cpml+PML+
non-air, so all flat/FEM/non-CPML runs — incl. Stage 1–4 validation —
are bit-identical). What it adds is the **functional** correctness of
the re-route when it actually engages.

Note: VT = Restrepo & Bielak 2014, *IJNME* 100:504–533 (fictitious-
domain 5-tet cut-cell; regular cubic octree preserved). FEM = whole-cube
air/void staircase. The flatten-PML re-route applies only to VT
(`get_topo_meth()==VT`); it is a strict no-op for FEM.
