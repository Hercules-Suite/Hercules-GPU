# Topography example: Gaussian hill (40 km × 40 km domain)

A ready-to-run Hercules topography example, geometry-matched to the
elastic-halfspace / `cps` benchmark (`region_length_east_m =
region_length_north_m = 40000`), so a topography run can be compared
directly against the validated flat-surface runs, and so it exercises
the VT + C-PML path (including the SPECFEM-style flatten-PML fix).

## Files

| File | Purpose |
|---|---|
| `gen_topography_gaussian_hill.py` | reproducible generator (single source of truth; parameterized) |
| `topography.in` | the generated elevation file (this is the **only** file Hercules needs in `topography_directory`) |

Regenerate / reparameterize:

```
python3 gen_topography_gaussian_hill.py            # defaults below
python3 gen_topography_gaussian_hill.py --height 1500 --sigma 2000   # steeper
```

## `topography.in` format (as implemented here)

- **Line 1**: grid spacing `d` in metres, the **same** in x and y.
- **Lines 2 …**: one elevation (metres) per grid point.
- The **first** point is the **SW (bottom-left) corner** `(x=0, y=0)`.
- Ordering (Hercules convention `x = south→north`, `y = west→east`):
  finish **all** points sharing a `y` before advancing `y` — i.e.
  inner loop = x (south→north), outer loop = y (west→east).

This example: `d = 100 m`, domain 40000 × 40000 m → **401 × 401 =
160801** elevation lines (+1 spacing line = 160802 total). Verified:
line 1 = `100`; line 2 = SW corner ≈ 0; line 80402 = centre = 800.0000;
last line = NE corner ≈ 0.

## Feature: Gaussian hill (interior, gentle by design)

```
elev(x,y) = H · exp( −((x−x0)² + (y−y0)²) / (2·σ²) )
```

Defaults: `H = 800 m`, `σ = 3000 m`, centre `(20000, 20000)`. Max slope
≈ `H/(σ·√e)` ≈ **9.2°** — gentle/moderate, well-resolved by the ~78 m
surface mesh. The hill is **centred and decays to ~0 well before the
lateral boundaries** (edge elevations are 0.0000). This is deliberate:

- Relief at the absorbing boundary is ~0 → the **well-behaved regime**;
- topography↔PML overlap is minimal → the C-PML **flatten-PML** path is
  exercised cleanly and the startup diagnostic
  `flattened-pml-skipped=N` should report a small/zero `N`.

For a flatten-PML **stress test**, widen `σ` / raise `H` / move the
centre so the hill reaches a lateral boundary, and watch
`flattened-pml-skipped` rise.

## Enabling it in the parameter file

Add to your `parameters.in` (values matched to the cps/elastic-halfspace
benchmark; only the topography block differs from a flat run):

```
include_topography          = yes
computation_method          = vt          # VT (fictitious-domain, Restrepo & Bielak 2014); recommended -- exercises VT + C-PML + flatten-PML.  Use FEM for the cube/air staircase baseline.
topography_directory        = <path>/benchmark_inputs/topography_gaussian_hill
```

Everything else (region size, source, stations, `abc_type`,
`type_of_damping`, …) stays as in the flat-surface elastic-halfspace
run so the only difference is the free-surface geometry.

## What to check at runtime / for validation

- Startup log: `Topography force path: gpu (... flattened-pml-skipped=N)`.
  With the default interior hill, `N` ≈ 0 (topography is ~flat in the
  PML band). With `abc_type=cpml` + `computation_method=vt`, any
  VT-deformed PML cube is re-routed to the regular cube path by the
  flatten-PML fix (`cpml_pml_flattened_element`); FEM is a strict no-op.
- Compare against the matching flat-surface run: away from the hill the
  wavefield should agree; over/around the hill, expect topographic
  scattering/amplification (the physical effect being modelled).
- This example is also the recommended basis for the still-pending
  **topo + C-PML validation** (stability + late-coda reflection from the
  topo→flat PML transition band, same gate style as the Stage 1–4 CPML
  validation).

Notes: VT (Virtual Topography) = Restrepo & Bielak, *IJNME* 2014,
100:504–533 (DOI 10.1002/nme.4756): a fictitious-domain scheme — regular
cubic octree preserved, topography-intersected cubes split into 5
non-conforming tetrahedra (empty tets discarded), order-2 convergence
preserved. FEM = the whole-cube air/void "staircase" baseline.
