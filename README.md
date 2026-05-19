For instructions on how to build, install and execute Hercules, please read
the documentation: https://clifford-yen.github.io/HerculesDoc/index.html 

Some general description of the core components in this directory are given as follows.

| File           | Description |
| -------------- | ----------- |
| `octor`        | Hercules mesh generator/partitioner/load balancer module
| `quake`        | Hercules solver module
| `user.mk`      | User-specified compiler and linker settings. This should be the only file we need to change.
| `systemdef.mk` | generic system-specific compile and linker settings. Tailored by entries in user.mk
| `common.mk`    | system-independent settings

Runtime performance toggles (no source edits needed):

Defaults listed below are either exported by `run_phaseA_phaseB.sh` or selected inside the solver/octor code when the corresponding envs are unset.

- `HERCULES_ADJUST_ASSIGN_OMP=1|0`
	Enables or disables OpenMP parallelization in assignment adjust core. Default is `1` (enabled).
- `HERCULES_ADJUST_DIST_OMP=1|0`
	Enables or disables OpenMP parallelization in parent-centric distribution adjust core. Default is `1` (enabled).
- `HERCULES_SPARSE_CLEAR_MODE=auto|range|nodes|full` (or `2|0|1|3`)
	Selects sparse host force reset mode in GPU solver path.
	Default is `full`/`3`.
	`auto`/`2` prefers `full` when topography touches a large fraction of local nodes, otherwise falls back to `range` or `nodes`.
	`range`/`0` uses exact contiguous ranges built from the topography node set and forces range `memset`.
	`nodes`/`1` can be faster for highly fragmented 3D model cases.
	`full`/`3` clears the full local host force array with one contiguous `memset`.
- `HERCULES_FULL_CLEAR_OMP=1|0`
	Enables or disables OpenMP parallelization for full-array host force clear when `HERCULES_SPARSE_CLEAR_MODE=full`. Default is `1`
	(enabled).
- `HERCULES_STATION_IO_OMP=1|0`
	Enables or disables OpenMP parallelization for station displacement interpolation and per-station file writes. Default is `1` (enabled).
- `HERCULES_STATION_ROUTE_MODE=interval|all`
	Controls station ownership lookup during setup. Default is `interval`, which uses the octree partition interval to skip stations that cannot belong to the current rank before calling the local octree search. `all` restores the older behavior where every rank searches every station.
- `HERCULES_PLANE_IO_OMP=1|0`
	Enables or disables OpenMP parallelization for per-strip plane
	interpolation before MPI send/write. Default is `1` (enabled).
- `HERCULES_CONTRIB_PACK_OMP=1|0`
	Enables or disables OpenMP parallelization while packing contribution messages for MPI communication. Default is `1` (enabled).
- `HERCULES_CONTRIB_UNPACK_OMP=1|0`
	Enables or disables OpenMP parallelization while unpacking contribution messages when a message's destination mapping has no duplicates. Default is `1` (enabled).
- `HERCULES_CONTRIB_RECV_MODE=waitsome|waitall`
	Selects how contribution messages are received and unpacked. Default is `waitsome`, which unpacks completed receives as they arrive so MPI receive wait and local unpack work can overlap. `waitall` restores the older deterministic receive behavior. Setting `HERCULES_REPRO_MODE=1` or `2` also forces deterministic `waitall`.
- `HERCULES_SPARSE_UPLOAD_MODE=auto|dense|sparse|range` (or `0|1|2|3`)
	Selects host-to-device force upload mode for topography sparse upload path. Default is `dense`/`1`.
	Use `auto`/`0` to enable the heuristic path selection explicitly.
	`dense`/`1` forces one full-array H2D copy and avoids sparse host gather, often better for large models where `Sparse force upload build` dominates.
	`sparse`/`2` forces sparse gather + scatter path.
	`range`/`3` copies cached merged contiguous spans derived from topography node ranges plus source nodes; this is an experimental middle ground between `dense` and `sparse`.
- `HERCULES_TOPO_FORCE_MODE=cpu|gpu` (or `0|1`)
	Selects CPU or GPU-native topography effective-force assembly. Default is `gpu`/`1`.
	`gpu`/`1` computes eligible VT topography effective-force elements directly into the device force buffer before stiffness/damping. It skips nonlinear-topography and topo-BKT elements the same way the CPU effective-force path does.
	Use `cpu`/`0` to force the legacy host topography effective-force path.
	For topography runs with `computation_method=fem`, there is no per-step effective-force topography assembly; the solver now treats that as no host-side topography force so the GPU source path can still skip the dense host-to-device force upload.
	The run log prints a one-line `Force path:` summary showing whether source force is on CPU/GPU, whether host force upload is required, and which topography effective-force path was selected.
- `HERCULES_DISP_SYNC_MODE=auto|full|boundary` (or `0|1`)
	Selects how updated displacements are synchronized between GPU and host. Default is `auto`.
	`full`/`0` keeps the legacy full-array copy each step.
	`boundary`/`1` keeps displacement history resident on the GPU and copies only the boundary/adjustment nodes needed by MPI displacement exchange, then patches those host-updated nodes back to the GPU for the next step.
	For correctness, `boundary` currently falls back to `full` when station/plane/checkpoint/DRM/nonlinear/fixed-base/geostatic/incident-plane-wave paths require full host displacement arrays. The run log prints `Displacement sync path:` and, when active, the boundary node count.

Advanced/diagnostic toggles:

- `HERCULES_REPRO_MODE=0|1|2`
	`0` disables repro mode. `1` disables remaining OpenMP-sensitive
	deterministic communication paths and uses deterministic MPI unpack ordering. `2` is stricter and also disables additional CPU OpenMP paths guarded by repro mode; `run_phaseA_phaseB.sh` does not auto-set this env.
- `HERCULES_ALLOW_MULTI_RANK_PER_GPU=0|1`
	Allows more than one local MPI rank to share the same GPU. Default is unset/`0`, which aborts when local MPI ranks exceed visible GPUs. This guard is intentional; set it only for diagnostic oversubscription runs.
- `HERCULES_MESH_IMPORT_VERIFY=0|1`
	Controls the extra convergence-check refinement pass when `mesh_input_mode=1`.
	Default is `0`. This applies only when `HERCULES_MESH_IMPORT_TREE_MODE=refine`.
- `HERCULES_MESH_IMPORT_TREE_MODE=cursor|refine`
	Controls how Phase B rebuilds the octree from an existing `mesh.e`.
	Default is `cursor`, which streams saved mesh leaves from the etree and
	imports them directly into the local octor tree. Set to `refine` to restore
	the older recursive `octor_refinetree()` search path.
- `OCTOR_IMPORT_PATH_CACHE=0|1`
	Controls a small parent-path cache used by cursor mesh import.
	Default is `1`, which reuses the shared root-to-leaf prefix between
	consecutive etree leaves. Set to `0`/`off`/`false`/`no` to disable it for
	comparison runs.
- `OCTOR_EXTRACT_PROFILE=0|1|2`
	Controls `octor_extractmesh()` profiling:
	`0` disables the Extractmesh profile.
	`1` prints the coarse stage profile: `setcom`, `vertex`, `indirect`,
	`node`, `elem`, `dnode`, `cleanup`, and `total`.
	`2` prints the coarse profile plus detailed `vertex`, `indirect`, and
	`node` sub-stage profiles. Mesh-input runs set this to `1` by default
	unless you define it yourself. Use `2` only for diagnostic runs because it
	adds extra timing/reduction work.
- `OCTOR_NODE_SORT_MODE=radix|qsort`
	Controls the node-table Morton/Z-order sort inside `octor_extractmesh()`.
	Default is `radix`, a linear-time 8-bucket radix sort over the 30 Morton
	bit planes. Set to `qsort` or `legacy` to restore the original C library
	`qsort()` path.
- `OCTOR_INIT_TASKS_PER_RANK=<integer>`
	Controls the number of coarse bootstrap octree tasks per MPI rank before
	parallel refinement. Default is the original octor value, `10`.
	Higher values are experimental; `1024` improved initial refinement balance
	in one run but caused an unsafe partitioning tree shape.

Build note:

- `USE_OPENMP` now defaults to `1` in `quake/forward/Makefile`.
	Use `USE_OPENMP=0` only if you want to disable OpenMP explicitly.
	Examples: `make` (OpenMP on by default) or `make USE_OPENMP=0`.

Two-phase script defaults:

- `run_phaseA_phaseB.sh` now auto-sets the following when they are unset:
	`OMP_PROC_BIND=close`, `OMP_PLACES=cores`,
	`HERCULES_ADJUST_ASSIGN_OMP=1`, `HERCULES_ADJUST_DIST_OMP=1`,
	`HERCULES_SPARSE_CLEAR_MODE=full`, `HERCULES_FULL_CLEAR_OMP=1`,
	`HERCULES_STATION_IO_OMP=1`, `HERCULES_STATION_ROUTE_MODE=interval`,
	`HERCULES_PLANE_IO_OMP=1`,
	`HERCULES_CONTRIB_PACK_OMP=1`, `HERCULES_CONTRIB_UNPACK_OMP=1`,
	`HERCULES_CONTRIB_RECV_MODE=waitsome`,
	`HERCULES_SPARSE_UPLOAD_MODE=dense`.
	These are the current recommended runtime defaults.
- OpenMP thread defaults:
	Phase A defaults to `PHASE_A_OMP_NUM_THREADS=1` when unset.
	Phase B uses `PHASE_B_OMP_NUM_THREADS` when set; otherwise it uses
	`OMP_NUM_THREADS` when set; otherwise it auto-derives a per-rank value
	for multi-rank-per-node launches.
- Phase B mapping defaults:
	`PHASE_B_MAP_MODE=auto` is the default and is the conservative recommended setting.
	For multi-rank-per-node GPU runs it uses node-level Open MPI placement such as
	`--map-by ppr:<ranks_per_node>:node:PE=<threads> --bind-to core --rank-by slot`.
	Available explicit modes are:
	`auto` (default conservative behavior),
	`per-node` (force one rank per node),
	`single-node` (single-node slot binding),
	`numa` (opt-in NUMA-aware placement for nodes where ranks/node divides evenly across local NUMA domains).
	The `numa` mode is experimental and should be used only when you want to test locality-sensitive placement explicitly.
- You can still override any of these by exporting your own value before running the script.

Two-phase job run (Phase A then Phase B):

- Script: `run_phaseA_phaseB.sh`
- Positional arguments: `PHASE_A_NP PHASE_B_NP [BASE_IN]`

Example (cluster case directory + external GPU binary):

```bash
cd /scratch/08458/wz_tacc/share/hercules/Hercules_Example/Test_case

export PHASE_A_OMP_NUM_THREADS=1
export PHASE_B_OMP_NUM_THREADS=64
# Optional: override script defaults when needed
# export OMP_PROC_BIND=spread
# export OMP_PLACES=threads
# export HERCULES_ADJUST_ASSIGN_OMP=0
# export HERCULES_ADJUST_DIST_OMP=0
# export HERCULES_SPARSE_CLEAR_MODE=auto
# export HERCULES_FULL_CLEAR_OMP=0
# export HERCULES_STATION_IO_OMP=0
# export HERCULES_STATION_ROUTE_MODE=all
# export HERCULES_PLANE_IO_OMP=0
# export HERCULES_SPARSE_UPLOAD_MODE=dense
# export HERCULES_SPARSE_UPLOAD_MODE=range
# export HERCULES_TOPO_FORCE_MODE=cpu
# Advanced/debug overrides:
# export HERCULES_CONTRIB_PACK_OMP=0
# export HERCULES_CONTRIB_UNPACK_OMP=0
# export HERCULES_CONTRIB_RECV_MODE=waitall
# export HERCULES_REPRO_MODE=1
# export HERCULES_MESH_IMPORT_VERIFY=1
# export PHASE_B_MAP_MODE=numa

PSOLVE_BIN=<path to psolve> \
<path to run_phaseA_phaseB.sh> \
288 4 \
inputfiles/parameter.in
```

Notes:

- With `BASE_IN`, the script auto-generates temporary Phase A/Phase B inputs by overriding `output_mesh`, `mesh_only_mode`, and `mesh_input_mode`.
- `BASE_IN` should define `mesh_etree_output_file`; both generated phase inputs are rewritten to use the same resolved mesh path. Relative mesh paths are resolved from the run directory, except when only the input-file-relative path already exists.
- Create the mesh output directory before running Phase A. For example, if the input uses `mesh_etree_output_file = outputfiles/mesh.e`, run `mkdir -p outputfiles` from the case directory first. The two-phase script prints `Mesh etree: ...` when it finds and rewrites this setting.
- Phase B mapping is controlled by `PHASE_B_MAP_MODE`.
	`auto` is the default and is recommended unless you are intentionally testing an alternative placement.
	For example, on NVL72 with `4` GPU ranks per node and `PHASE_B_OMP_NUM_THREADS=35`,
	`auto` typically uses a node-level mapping like `ppr:4:node:PE=35`,
	while `PHASE_B_MAP_MODE=numa` tries `ppr:2:numa:PE=35` when the local node exposes two NUMA domains.

Mesh mode input flags (input file):

- `output_mesh`
	`1` writes the generated mesh etree (`mesh_etree_output_file`) to disk.
	`0` skips mesh etree output.
- `mesh_only_mode`
	`1` runs meshing/setup and exits before time-stepping solver iterations.
	`0` runs the full solver after setup.
- `mesh_input_mode`
	`1` builds mesh from an existing mesh etree (import path).
	`0` builds mesh from the velocity/topography models (normal meshing path).

Typical two-phase settings:

- Phase A (generate mesh database):
	`output_mesh=1`, `mesh_only_mode=1`, `mesh_input_mode=0`
- Phase B (import mesh and solve):
	`output_mesh=0`, `mesh_only_mode=0`, `mesh_input_mode=1`

Acoustic-elastic input flags (input file):

- `include_acoustic`
	`yes` enables acoustic material handling for water/fluid regions. Elements
	whose material properties satisfy `Vp > 0`, `rho > 0`, and
	`Vs <= 1.0e-6` are treated as acoustic. Acoustic elements have `Vs` forced
	to zero and their material damping coefficients are cleared during mesh
	setup/import. If this is `no`, the solver aborts when it encounters
	nonpositive `Vs` in the material model.
- `acoustic_interface_coupling`
	Optional. `yes` enables the current normal-projected coupling treatment on
	acoustic element faces adjacent to elastic material. Default is `yes`, but
	it is only active when `include_acoustic=yes`.
- `acoustic_interface_max_level`
	Optional. When set to a nonnegative octree level, mixed acoustic/elastic
	octants are refined up to that level during meshing. Default is `-1`, which
	disables this extra interface-refinement limit.
- `acoustic_mesh_refinement_factor`
	Optional. Extra refinement multiplier for fully acoustic elements during
	fresh mesh generation. The acoustic mesh-size target becomes
	`Vp / (simulation_node_per_wavelength * simulation_wave_max_freq_hz *
	acoustic_mesh_refinement_factor)`, while elastic elements still use the
	usual shear-wave rule. Default is `1.0`. Values must be `>= 1.0`. This does
	not replace `acoustic_interface_max_level`: that parameter refines mixed
	water/solid interface octants, while this parameter refines the acoustic
	water volume. When using a two-phase workflow, set this before generating
	the mesh etree in phase A; importing an already-generated mesh preserves the
	existing mesh resolution.
- `acoustic_pressure_formulation`
	Optional. `yes` enables the pressure-DOF acoustic solver path. Acoustic
	elements assemble scalar pressure mass/stiffness instead of displacement
	mass/stiffness, while elastic elements continue to use displacement DOFs.
	This requires `include_acoustic=yes`. On a flat model surface, acoustic
	pressure nodes on the top external domain face are constrained to `p=0`
	for a pressure-release water surface. When `acoustic_interface_coupling=yes`,
	the pressure formulation builds an explicit acoustic-solid interface table.
	In octree meshes, nonconforming interface faces are split into area-weighted
	fragments with fragment centers and normals, and CUDA kernels apply the
	normal traction/acceleration coupling only on those fragments. Each fragment
	interpolates pressure on the acoustic element and displacement acceleration
	on the adjacent elastic element, so refined elastic neighbors can couple to
	coarser acoustic faces. Fragment subdivision follows
	`acoustic_interface_max_level` and finer local neighbors.
	When the neighboring element is local, the table uses that element's
	mesh-level acoustic/elastic classification, so mixed elements kept elastic
	still couple to adjacent acoustic elements. Subdivision has a per-face cap of
	`32 x 32` fragments to avoid accidental table blow-up.
- `acoustic_interface_diagnostics_file`
	Optional pressure-interface diagnostic output path. When this path is set
	and `acoustic_pressure_formulation=yes`, the solver writes one row per local
	acoustic-solid interface fragment after mesh setup. The file records the
	acoustic element id, face id, fragment area, fragment center in local and
	model coordinates, outward normal, element level/size, and adjacent elastic
	element local index/id and local fragment coordinates when the adjacent
	element is local. This file can be very large, so it is disabled when unset.
- `print_station_pressures`
	Optional. `yes` writes scalar pressure station histories when
	`acoustic_pressure_formulation=yes`. Pressure is appended as a
	`Pressure(Pa)` column in the normal `station.N` files after any requested
	displacement, velocity, and acceleration columns. For stations whose
	containing element is acoustic, the motion columns are reconstructed from
	the pressure gradient; for elastic stations they continue to use the
	elastic displacement field. The pressure column is `nan` for stations whose
	containing element is elastic, because scalar acoustic pressure is not
	defined there.
- `reconstruct_acoustic_plane_displacements`
	Optional. `yes` reconstructs displacement at acoustic output-plane points
	from the pressure gradient and writes those reconstructed values into the
	standard `planedisplacements.N` files. Elastic plane points continue to use
	the elastic displacement field, so mixed water/solid planes remain a single
	three-component displacement file. This requires `include_acoustic=yes` and
	`acoustic_pressure_formulation=yes`. Default is `no`, which preserves the
	legacy plane-output behavior.
- `basin_material_type`
	Optional. Controls material assigned inside `basin_velocity_model` geometry.
	Default is `elastic`, which preserves the historical hard-coded basin
	material. Set to `acoustic_water` to assign water properties inside the
	basin geometry: `Vs=0`, `Qs=0`, `Qp=0`, and configurable `Vp/rho`. This
	requires `include_acoustic=yes`.
- `basin_water_vp`, `basin_water_rho`
	Optional acoustic-water basin properties used only when
	`basin_material_type=acoustic_water`. Defaults are `1500.0 m/s` and
	`1000.0 kg/m3`.
- `station_mediums_output_file`
	Optional station diagnostic output path. The solver writes this one-time
	metadata file during station setup. It records each station's owning medium
	(`acoustic` or `elastic`), whether pressure is defined there, owner rank,
	element id/level, local coordinates, and material properties. If unset, it
	defaults to `output_stations_directory/station_mediums.txt`.
- `acoustic_formulation`
	Optional. **`velocity_potential` (default)** or `pressure`. Selects
	the acoustic–elastic interface coupling variable. `velocity_potential`
	is the stable, cross-validated default; `pressure` is the legacy path,
	kept as an explicit opt-in and bit-identical to earlier validated
	pressure runs. Only affects simulations with acoustic elements /
	acoustic–elastic interfaces — pure-elastic runs are unaffected and
	bit-identical either way. The default removes the added-mass / CFL
	instability the `pressure` path showed for a thin stiff water layer
	over an elastic halfspace (displacement validated against CPS at
	corr ~0.99–0.998). Requires `include_acoustic=yes` and
	`acoustic_pressure_formulation=yes`. Set `= pressure` to restore the
	exact legacy behavior.
- `acoustic_quality_factor`
	Optional. Default unset / `<= 0` ⇒ acoustic is exactly lossless and
	bit-identical to the validated baseline. A positive quality factor
	`Q` adds a small, high-frequency-selective damping to acoustic
	(water) elements that attenuates the synthetic lossless-cavity grid
	ring while barely touching the physical signal. Applies to interior
	acoustic elements only (skipped in the CPML zone); works for both
	`pressure` and `velocity_potential`. **This is a numerical
	regularizer for a low-order-discretisation artifact, not physical
	ocean attenuation** (and is non-SPECFEM) — use it only to clean the
	synthetic no-attenuation benchmark; runs with a realistically
	attenuating Earth do not need it. Empirical sweet spot
	`acoustic_quality_factor = 100` (sensible band 100–200): lower
	over-perturbs, higher under-damps.

Acoustic-pressure formulation — correctness fix & known limitation:

- **Interior stiffness factor-of-2 correction.** An earlier acoustic
	Laplacian coefficient was a factor of two too small, making the
	discrete acoustic speed ~0.707·Vp everywhere in acoustic media. It is
	now correct. This **changes all acoustic-pressure results** versus
	pre-fix builds (in-water vertical-displacement correlation vs the
	flat-water CPS analytic improved from ~0.83 to ~0.997). Elastic-only
	/ non-acoustic runs are unaffected.
- **Lossless in-water pressure-coda ring (known limitation).** With no
	attenuation, the in-water *pressure* coda carries a non-decaying
	grid-scale ring absent from the exact solution. It is
	displacement-invisible (the validated 0.997 displacement is
	unaffected) and is a fundamental property of the lossless pressure
	formulation, not a code bug. With any realistic ocean Q it damps
	physically; for a clean lossless-benchmark comparison apply a
	post-processing low-pass to the output pressure, or opt into the
	in-model `acoustic_quality_factor` regularizer (see that flag).

Absorbing boundary input flags (input file):

Hercules supports two absorbing boundary conditions, selected at run
time. The default is the legacy first-order Lysmer dashpot. The newer
option is an unsplit Convolutional Perfectly Matched Layer (C-PML,
Komatitsch & Martin 2007 — the same second-order displacement form
SPECFEM3D's FE-PML uses). The C-PML is per-element, wraps the five
absorbing faces (`x_min`, `x_max`, `y_min`, `y_max`, `z_max`) and leaves
`z_min` as the free surface, and works for both elastic and
acoustic-pressure elements. Validated against CPS analytical
seismograms it is uniformly more accurate than the Lysmer dashpot at
every offset/component (see the summary below).

- `abc_type`
	Optional. `dashpot` (default) = legacy Lysmer dashpot. `cpml` =
	unsplit C-PML. If unset, behavior is unchanged from earlier versions.
- `cpml_thickness_m`
	Target PML thickness in metres. When `> 0`, used directly. When `0`
	or unset (default), the solver uses `cpml_layers ×` the element size
	at the absorbing boundary and logs the realized thickness at startup.
- `cpml_layers`
	Optional integer. Default `10`. Used only when `cpml_thickness_m` is
	`0` — "this many element layers deep" at the boundary, without needing
	to know the local mesh resolution.
- `cpml_reflection_target`
	Optional. Target theoretical reflection coefficient. Default `1.0e-3`.
	Smaller = stronger absorption but tighter stability margins.
- `cpml_poly_order`
	Optional integer in `[1, 4]`. Default `2` (the standard
	K&M/SPECFEM3D choice). Controls how steeply the damping ramps up
	across the PML.
- `cpml_kappa_max`
	Optional. Default `1.0` — no extra coordinate stretch, and the
	validated setting. Values `> 1` improve grazing/evanescent absorption
	and are now consistently implemented, but still need re-validation;
	treat `2–3` as promising-but-unconfirmed and prefer a thicker PML
	(`cpml_layers` ~30–50) for the deep-corner vertical residual. Values
	above `3` are untested. `kappa = 1` keeps the PML↔interior interface
	exact.
- `cpml_alpha_max_hz`
	Optional. When `> 0`, sets the CFS frequency directly (Hz). When `< 0`
	or unset (default), it is derived from `simulation_wave_max_freq_hz`.
	Set it to your source peak frequency if that differs markedly from
	`simulation_wave_max_freq_hz` (e.g. a Ricker pulse well below the
	mesh resolution frequency).
- `cpml_cfl_safety`
	Optional. Default `0.9`. The critical time step is multiplied by this
	factor when `abc_type=cpml`, for margin against the corner-element
	instability. `1.0` disables the reduction; `< 0.5` is rarely needed.
- `cpml_mpml_cross_fraction`
	Optional. Default `0.1`. Strength of the M-PML cross-axis coupling
	(Meza-Fajardo & Papageorgiou 2008) that stabilises long elastic runs.
	`0` disables it (plain unsplit CPML, known to drift over long runs).
	Recommended range `[0.05, 0.20]`; valid `[0, 0.5]`. Lower it if
	receivers very close to the PML inner edge show extra coda; raise it
	for stability margin in pathological cases.
- `cpml_alpha_floor_fraction`
	Optional. Default `0.1`. The CFS profile otherwise reaches zero at
	the outer PML edge, which destabilises corner elements; flooring it
	at 10% of the maximum keeps a non-zero decay rate everywhere. Set to
	`0` to disable (only safe with M-PML cross-coupling on). Valid
	`[0, 1]`.
- `cpml_pml_material_damping`
	Optional. `none` (default) | `bkt` | `rayleigh`. Whether the elastic
	PML zone is physically attenuating or lossless-elastic. `none` skips
	all material damping on PML elements — correct **only** for a
	non-attenuating interior (`type_of_damping=none`); it is the
	bit-identical historical baseline there. **Important:**
	`type_of_damping=bkt*` with `abc_type=cpml` under `none` is
	catastrophically unstable (diverges to NaN by t≈40 s). Resolution
	policy with an attenuating interior (BKT-family or Rayleigh):
	(1) **unset** → auto-promoted to the matching viscoelastic PML
	(`bkt`/`rayleigh`), safe by default; (2) **explicit `none`** →
	honored, with a loud KNOWN-UNSTABLE warning (a deliberate escape
	hatch for attribution/regression); (3) **explicit `bkt`/`rayleigh`**
	→ used as set (a mismatched value is a hard error). `bkt` is fully
	validated (stable, correct, corr ~0.99 vs the matching-attenuation
	CPS analytic); `rayleigh` is implemented and pending a validation
	run. Note: a related correctness fix ported stiffness-proportional
	Rayleigh damping to the GPU, so **existing `type_of_damping=rayleigh`
	GPU results change** (they now match the CPU path); every
	non-Rayleigh run stays byte-identical.
- `cpml_clamp_outer`
	Optional. `yes` or `no` (default `no`). `yes` pins displacement to
	zero on the outer PML faces (never `z_min`, the free surface),
	converting round-trip absorption into single-pass. Useful for very
	thin PMLs or very long runs that leak round-trip reflections; default
	`no` keeps the classic free outer face.
- `cpml_stiffness_mode`
	Optional. `native` (default) or `reduced`. How the elastic PML
	element stiffness is formed.
	* `native` — full native element stiffness plus a CPML correction
		that is exactly zero where there is no damping, so the
		PML↔interior operator is continuous and not rank-deficient.
		Recommended.
	* `reduced` — standalone reduced CPML stiffness only (legacy).
		Differs from `native` by ~0.2–0.4% on the homogeneous CPS
		benchmark; kept for comparison/debugging. `native` is preferred
		for coarse/heterogeneous meshes and long runs.
- `cpml_stiffness_mode_acoustic`
	Optional. `native` or `reduced`. Independent override for **acoustic**
	PML elements; if unset it inherits `cpml_stiffness_mode` (so existing
	configs are unchanged). Lets you exercise the elastic and acoustic
	PML stiffness paths independently. The active modes are printed at
	startup.
- `topography_bkt_gpu_stabilize`
	Optional, tri-state: `auto` (default), `yes`, or `no`. Fixes a
	**pre-existing** GPU instability where `computation_method=vt`
	(Virtual Topography) together with any `type_of_damping=bkt*`
	diverges to NaN (ABC-independent; unrelated to the C-PML work).
	**`auto` (default) is correct-by-default:** it enables the fix
	exactly for `vt`+`bkt*` and is a strict no-op, bit-identical to the
	validated baseline, for every other configuration (FEM,
	flat/no-topography, non-VT, non-BKT, all CPML paths). `yes` forces it
	on (the VT guard still applies); `no` is an explicit escape hatch
	that restores the pre-fix behaviour (the VT+BKT NaN) for A/B
	comparison. The only behaviour change vs. the old default is that
	`vt`+`bkt*` — previously an unconditional NaN — now produces the
	correct viscoelastic solution, matching the validated CPU operator.
	Flat-PML topography cells stay on the cube-stiffness + CPML path
	exactly as before.

Notes on C-PML behavior:

- Each PML element is automatically classified as elastic or acoustic
	(the same test the acoustic-pressure path uses), so heterogeneous
	models and acoustic/elastic coupling need no extra configuration.
- Material damping on PML elements follows `cpml_pml_material_damping`
	(see that flag). With the default `none` all material damping is
	skipped on PML elements — correct only for a non-attenuating
	interior; `bkt`/`rayleigh` carry the attenuation into the PML
	consistently. Non-PML elements always use material damping as usual;
	dashpot mode and no-attenuation runs are bit-identical.
- The startup log prints how many elastic/acoustic PML elements were
	found and the realized PML thickness — use it to confirm the PML zone
	is non-empty and split between materials as expected.
- **Topography + C-PML — taper your relief to flat before the PML
	(recommended), do not rely on flatten-PML alone.** The C-PML buffer
	is by construction a flat axis-aligned layer (as in SPECFEM3D); real
	relief is not carried into it. For `computation_method=vt` a
	flatten-PML re-route forces the lateral-PML elements flat so
	VT+topography+C-PML stays stable, and the startup log reports how
	many elements it skipped. This is a **stability backstop, not a
	substitute for tapering**: relief still tens-to-hundreds of metres
	high at the PML inner edge gets flattened there into an artificial
	step (a "cliff" for steep boundary relief) that scatters spuriously.
	For real models, taper `topography.in` smoothly down to flat inside
	an interior buffer that ends flat — with margin — before the PML
	inner edge (flat collar ≥ PML thickness, plus a smooth ramp). The
	helper `benchmark_inputs/taper_topography.py` does exactly this for
	any `topography.in`; the two example generators also take
	`--edge-taper`. (FEM is a strict no-op for flatten-PML; non-CPML /
	non-VT / flat runs are byte-identical. The
	`benchmark_inputs/topography_tilt_boundary` example is a deliberate
	worst-case test, not a model-setup template.)

Using `benchmark_inputs/taper_topography.py` (the recommended
preprocessing step for any real topography + C-PML run):

- **What it is.** A standalone, dependency-free Python 3 script (stdlib
	only — runs anywhere, including a cluster login node). It reads a
	`topography.in`, smoothly tapers the relief to a flat reference near
	the boundary, and writes a new `topography.in` that is flat in (and
	well before) the PML band while leaving the deep interior untouched.
	Same file format in and out; no solver change.
- **Basic use** (defaults target the 40 km × 40 km benchmark domain;
	spacing is read from the file's line 1):

	```
	cd benchmark_inputs
	python3 taper_topography.py path/to/topography.in --out topography_tapered.in
	```

	Then point Hercules at the tapered file — Hercules expects the
	filename `topography.in` in `topography_directory`, e.g.:

	```
	mkdir -p model/topo_tapered
	python3 taper_topography.py raw/topography.in \
	    --out model/topo_tapered/topography.in --force
	# parameters.in:  topography_directory = .../model/topo_tapered
	```

- **Domain / PML options** (tune to your model, not the 40 km default):
	* `--lx` / `--ly` — domain extent S–N (x) and W–E (y) in metres
		(defaults 40000 / 40000). The script infers the grid from these +
		the file's spacing and errors out if the point count does not
		match, so a wrong `--lx/--ly` is caught immediately.
	* `--buffer` — total taper band width inward from each lateral
		boundary (m); default 6000.
	* `--flat-margin` — width of the guaranteed-dead-flat collar adjacent
		to every boundary (default 1500). **Set this ≥ your PML thickness
		with margin** (the startup log line `CPML: … (thickness X m …)`
		gives X; the example domain's ~10-layer PML is ≈ 780 m). The
		script warns if it is below `--pml-thickness`.
	* `--pml-thickness` — PML band thickness (m), used only for that
		safety warning (default 780).
	* `--z-ref` — the flat reference: a number in metres (default `0.0`),
		or `auto` to use the mean elevation just inside the taper band.
		Use `auto` for real DEMs sitting at a large offset so the taper
		ramps to a representative interior level instead of an arbitrary
		0.
	* `--force` — allow overwriting the output (or in-place when `--out`
		equals the input).
- **What it prints** — the resolved grid, the buffer/collar/ramp widths,
	`z_ref`, how many points were modified vs left untouched, the max
	elevation change, and a perimeter-flatness self-check that must print
	`OK`. If it is not `OK`, do not use the file.
- **`--help`** lists every option with its default.
- **Alternative for the two bundled examples.** Regenerate either
	example surface already tapered with `python3 gen_topography_tilt.py
	--edge-taper` (or `gen_topography_gaussian_hill.py --edge-taper`),
	with matching `--taper-buffer` / `--taper-flat-margin` /
	`--taper-zref`. Output is identical to running the standalone utility
	on the un-tapered surface; the flag is off by default so the
	committed example files stay reproducible.

Stability fixes baked in:

These are automatic and tuned by the fractions above; you normally do
not touch them.

- **Alpha floor** (`cpml_alpha_floor_fraction`): keeps a non-zero decay
	rate everywhere in the PML so corner elements can't lose all damping
	at once. No-op at the inner edge.
- **M-PML cross-axis coupling** (`cpml_mpml_cross_fraction`): mixes a
	little of each axis's damping into the other two to stabilise long
	elastic runs, at the cost of slightly more absorption leaking onto
	auxiliary axes.
- **Centroid-paired CFS terms**: the grazing-incidence CFS terms are
	evaluated consistently at the element centroid so they cancel exactly
	where there is no damping (otherwise PML corners grow exponentially).
- **Inner-edge half-cell sampling**: the profile is sampled at the
	element centroid, so the first PML layer carries a small, bounded,
	well-characterised discrete reflection — shrink it with more
	`cpml_layers` or a higher `cpml_poly_order` if needed. (An alternative
	inner-face-shift refinement from another branch was evaluated and not
	adopted here; the centroid baseline is already validation-grade.)

CPS analytical validation (elastic half-space, dip-slip point source):
in the late coda — where absorbing-boundary quality dominates — the
C-PML is **2–5× more accurate than the Lysmer dashpot** at near-PML
receivers and best of every implementation tested, while matching CPS
as closely as Lysmer in the ABC-independent direct-arrival window (it
does not contaminate the physical wavefield). Far-field receivers
(beyond ~1 PML thickness) are already accurate with any method. **For
benchmarking PML quality, keep receivers at least one PML thickness
away from the inner edge** — a receiver inside the PML sees the
modified strain field, not a physical seismogram.

Recommended (validated) parameters:

```
abc_type = cpml
cpml_clamp_outer = yes      # single-pass absorption
cpml_layers = 30            # raise to 30–50 to close the deep-corner Z residual
# cpml_kappa_max:           leave at 1.0 (validated); 2–3 promising but unconfirmed
# cpml_stiffness_mode:      leave at native (recommended)
# cpml_alpha_max_hz:        leave default unless the source peak freq differs
# cpml_mpml_cross_fraction / cpml_alpha_floor_fraction: leave at 0.1
```

Notes:

- The validated configuration is `cpml_kappa_max = 1`; the deep-corner
	vertical residual is reduced most reliably by a thicker PML
	(`cpml_layers` 30–50), not by `kappa`.
- `cpml_clamp_outer = yes` is recommended for long runs / thinner PMLs.
- `cpml_layers = 10` already beats the dashpot 2–5× on the horizontal
	components; 30–50 helps the deep-corner vertical and low-frequency
	absorption.
- For double-precision correctness checks build with `solver_float =
	double`; production single-precision runs are verified back-to-back.


