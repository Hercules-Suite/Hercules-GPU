For instructions on how to build, install, configure, and run Hercules, see the
[Hercules GPU Documentation](https://hercules-gpu-docs.liling-dpim.chatgpt.site).

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

Phase launcher defaults:

- `run_phaseA_phaseB.sh` now auto-sets the following when they are unset:
	`OMP_PROC_BIND=close`, `OMP_PLACES=cores`,
	`HERCULES_ADJUST_ASSIGN_OMP=1`, `HERCULES_ADJUST_DIST_OMP=1`,
	`HERCULES_SPARSE_CLEAR_MODE=full`, `HERCULES_FULL_CLEAR_OMP=1`,
	`HERCULES_STATION_IO_OMP=1`, `HERCULES_STATION_ROUTE_MODE=interval`,
	`HERCULES_PLANE_IO_OMP=1`,
	`HERCULES_CONTRIB_PACK_OMP=1`, `HERCULES_CONTRIB_UNPACK_OMP=1`,
	`HERCULES_CONTRIB_RECV_MODE=waitsome`,
	`HERCULES_SPARSE_UPLOAD_MODE=dense`,
	`HERCULES_PHASE_MODE=unified`.
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

Phase launcher (unified by default):

- Script: `run_phaseA_phaseB.sh`
- Positional arguments: `PHASE_A_NP PHASE_B_NP [BASE_IN]`

Example (cluster case directory + external GPU binary):

```bash
cd /scratch/08458/wz_tacc/share/hercules/Hercules_Example/Test_case

export PHASE_A_OMP_NUM_THREADS=1
export PHASE_B_OMP_NUM_THREADS=64
# To request separate Phase-A and Phase-B launches instead:
# export HERCULES_PHASE_MODE=two-phase
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

- With `BASE_IN`, the script auto-generates temporary Phase A, Phase B, and
	unified inputs. Two-phase inputs override all mesh-mode flags; unified mode
	preserves the base input's `output_mesh` behavior while setting
	`mesh_only_mode=0` and `mesh_input_mode=0`.
- `BASE_IN` must define `mesh_etree_output_file` for two-phase runs and unified
	runs with `output_mesh=1`; it is optional for unified `output_mesh=0`.
	Relative paths are resolved from the run directory, except when only the
	input-file-relative path already exists.
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
	Optional. **`velocity_potential` (default)** or `pressure`. Selects the
	acoustic-elastic interface coupling variable. `velocity_potential` is
	now the default because it is the cross-validated, stable coupling;
	`pressure` is the legacy path, retained as an explicit opt-in
	(`acoustic_formulation = pressure`) and bit-identical to all earlier
	validated pressure runs. This parameter only affects simulations with
	acoustic elements / acoustic-elastic interfaces — **pure-elastic runs
	are unaffected and bit-identical** regardless of the setting.
	`pressure` (legacy): the fluid feels the solid normal **acceleration**
	(a 2nd time-difference, `~omega^2`) and the solid feels the fluid
	pressure. `velocity_potential` (default) reinterprets the acoustic
	field as a velocity potential `phi` (the bulk solver — Laplacian,
	lumped mass, central difference — is reused unchanged): the interface
	exchanges the solid normal **velocity** `(u^{n+1}-u^n)/dt . n` (a 1st
	time-difference, `~omega`) and the solid feels `p = -dphi/dt`; the
	pressure-release free surface becomes `phi^{n+1}=phi^n`; station
	pressure value/gradient output is recovered as `-dphi/dt`. The milder
	first-order-in-time coupling removes the explicit added-mass / CFL
	instability and the 0 km symmetry-axis secular growth that the
	`pressure` path exhibits for a thin stiff water layer over an elastic
	halfspace, with displacement validated against CPS at corr
	~0.99–0.998. Requires `include_acoustic=yes` and
	`acoustic_pressure_formulation=yes`. (Independently reproduces the
	`acoustic-pressure-pml-gpu-v3` branch's `velocity_potential` result —
	used as a cross-validation.) To restore the exact legacy behavior of
	pre-default-flip runs, set `acoustic_formulation = pressure`.
- `acoustic_interface_filter`
	Optional. `no` (default) or `yes`. Enables a per-face 1-pole IIR /
	exponential-moving-average low-pass on the acoustic-elastic interface
	source term, intended to attenuate the grid-scale (`~c/4h ~ 4.6 Hz` in
	the flat-water benchmark) content that pumps the lossless water-column
	eigenmode while passing the physical sub-1 Hz band. Default `no` is
	bit-identical to earlier runs.
	**NOT RECOMMENDED — found destabilizing.** Empirically (independent
	test, default `alpha=0.15`/`strength=1.0`), enabling this on top of
	`acoustic_formulation=velocity_potential` makes the run physical for
	~25 s and then diverge exponentially (e-fold ~1.5 s; `|p|` -> 1e16).
	The filter inserts a stateful IIR pole onto the coupling source
	*inside* the two-way fluid-solid feedback loop, adding a destabilizing
	lag to the (delicate, explicit, loosely-coupled but otherwise stable)
	`velocity_potential` scheme — the same class of in-loop-extra-dynamics
	instability seen with the chi predictor/lagged-coupling experiments.
	The residual lossless-fluid pressure ring is better treated as the
	documented fundamental limitation (it damps physically with any
	realistic Q) or, for a clean lossless-benchmark *comparison*, by a
	post-processing low-pass on the OUTPUT pressure (outside the solver
	loop, cannot affect stability) — not by this in-loop filter. Kept
	in-tree (default off, code is a faithful EMA) for the record / future
	passive-form experiments.
- `acoustic_interface_filter_alpha`
	Optional. Default `0.15`. EMA coefficient in `filtered = old +
	alpha*(source - old)`; smaller `alpha` = stronger low-pass (lower
	corner frequency). Only used when `acoustic_interface_filter=yes`.
- `acoustic_interface_filter_strength`
	Optional. Default `1.0`. Blend toward the low-passed value:
	`source += strength*(filtered - source)`. `1.0` fully replaces the raw
	source with the filtered one; `0.0` is a no-op. Only used when
	`acoustic_interface_filter=yes`.
- `acoustic_quality_factor`
	Optional. Default unset / `<= 0` ⇒ acoustic is **exactly lossless and
	bit-identical** to the validated baseline. When set to a positive
	quality factor `Q`, acoustic (water) elements get a small
	**stiffness-proportional (Kelvin–Voigt)** damping added directly in
	the acoustic force kernel: `f -= coeff · beta · K · (pⁿ − pⁿ⁻¹)`,
	with the same stiffness operator `K` the acoustic force already uses.
	`beta` is derived self-containedly from `Q`: the leapfrog modal EOM
	for this term is `p'' + (beta·ω²·Δt)·p' + ω²·p = 0`, i.e. damping
	ratio `ζ(ω) = beta·ω·Δt/2` and `Q(ω) = 1/(beta·ω·Δt)`; pinning `Q`
	at the simulation reference frequency `theFreq` gives
	**`beta = 1 / (Q · 2π·theFreq · Δt)`**. It deliberately does **not**
	depend on `theBBase`: `compute_setab` zeroes `theBBase` for every
	non-`RAYLEIGH`/`MASS` run (`NONE` and all `BKT*`) — exactly the
	no-attenuation lossless benchmark this knob targets — so an earlier
	`theBBase`-scaled `beta` silently no-op'd there (caught when the
	`Q=500` run reproduced the no-damping baseline within GPU float
	noise; fixed). It is also deliberately **not** the elastic `zeta`
	route: acoustic elements never traverse the elastic `c3/c4` `b·K`
	path (verified — `kernel.cu` does not reference `c3/c4`), and the
	shared `zeta` would add only the wrong-sign mass-proportional `a·M`.
	The damping ratio ∝ ω (Q ∝ 1/ω, high-frequency
	selective): it attenuates the lossless-cavity `~c/4h` (~4.6 Hz) grid
	ring while barely touching the sub-source-band physical signal.
	Applies to **interior** acoustic elements only — it is skipped on
	CPML PML elements so the PML zone stays pure CPML (consistent with
	the material-damping skip in commit `5d91b09`). Works for both
	`pressure` and `velocity_potential`. **Framing (important):** this is
	a *numerical regularizer for a low-order-discretisation artifact, not
	physical ocean attenuation* — real seawater Q at the seismic band is
	~1e5 (PREM ocean Qκ ≈ 5.8e4), which would change nothing observable —
	and it is **non-SPECFEM**: SPECFEM3D treats the acoustic/fluid domain
	as lossless (no SLS, no artificial viscosity) and avoids this ring
	via high-order spectral elements, not fluid damping. Use it only to
	clean the synthetic no-attenuation lossless benchmark; any run with a
	realistically attenuating solid Earth does not need it. Validate
	before relying on it: long-run stable, displacement still corr ~0.997
	vs CPS (damping small enough not to perturb the validated physics),
	and quantified ring reduction.
	**Empirically-determined sweet spot: `acoustic_quality_factor = 100`**
	for the flat-water no-attenuation lossless benchmark. A Q sweep
	(500 → 100 → 50, the v23/v24/v25 runs) located the knee at Q≈100:
	at the diagnostic far water stations it removes ~50–82 % of the
	>2 Hz ring energy and lifts the pressure correlation with the exact
	ring-free CPS to ~0.86 (from ~0.68 lossless), while Z-displacement
	corr vs CPS stays bit-for-bit-identical at 0.9965 (mean perturbation
	≤ 0.5 %) and the run is stable. Q=50 doubles the perturbation cost
	for ≲1 % additional ring removal (the remaining high-frequency
	content is physical reverberation that CPS also has, and correctly
	does not respond to more damping); Q=500 under-damps. Sensible band
	Q ∈ [100, 200].

Acoustic-pressure formulation — correctness fix & known limitation:

- **Interior stiffness factor-of-2 correction.** The consistent
	trilinear-hex acoustic Laplacian coefficient was `dt^2*h/(32*rho)`,
	a factor of two too small, which made the discrete acoustic speed
	`Vp/sqrt(2)` (~0.707 Vp) everywhere in acoustic media. It is now
	`dt^2*h/(16*rho)` (derived two independent ways; the lumped `1/kappa`
	mass is exact, so the stiffness prefactor alone sets the speed). This
	**changes all acoustic-pressure results** versus pre-fix builds:
	against the flat-water CPS analytic reference the in-water vertical
	displacement correlation went from ~0.83 to ~0.997 and a ~1.75x
	amplitude overshoot was eliminated. Elastic-only / dashpot /
	non-acoustic runs are unaffected (the change is gated by
	`element_is_acoustic`).
- **Lossless in-water pressure-coda ring (known limitation).** In the
	no-attenuation flat-water benchmark, the in-water *pressure* coda
	carries a non-decaying `~c/4h` (~4.6 Hz) ring absent from the exact
	CPS solution. Root cause (diagnosed in depth): the interface coupling
	injects the solid's discrete acceleration, which is intrinsically
	grid-amplified (`M^-1 K u`), into a zero-loss resonator. It is
	*displacement-invisible* (the validated 0.997 displacement is
	unaffected) and is reproduced identically by two independent
	implementations, i.e. it is a fundamental property of the lossless
	pressure formulation, not a code bug. `acoustic_formulation=
	velocity_potential` removes the associated instability/growth and
	weakens the pump (`~omega` vs `~omega^2`) but does not eliminate the
	ring in a strictly lossless fluid. The in-loop `acoustic_interface_
	filter` was tried as a targeted suppressor and **rejected — it
	destabilizes `velocity_potential`** (an IIR lag inside the FSI
	feedback loop; see that flag's entry). The residual ring is therefore
	left as a documented fundamental limitation: with any realistic Q
	(attenuating ocean) it damps physically, and for a clean lossless-
	benchmark comparison a post-processing low-pass on the OUTPUT pressure
	(outside the solver loop) is the safe route. An opt-in in-model
	alternative also exists — `acoustic_quality_factor` adds a small
	stiffness-proportional (Kelvin–Voigt) damping that attenuates the
	ring at source (high-frequency selective; default off / lossless;
	see that flag's entry) — but it is a non-SPECFEM numerical
	regularizer, so the post-processing route remains the conservative
	default. SPECFEM-spirit
	conservative-mass and
	displacement-potential alternatives were tried and rejected (the
	latter is explicit-CFL-unstable here); see the git history /
	`visualizations/` for the full v3→v19 diagnostic chain.

Absorbing boundary input flags (input file):

Hercules supports two absorbing boundary conditions, selected at run
time. The default is the legacy first-order Lysmer dashpot, folded into
the lumped mass matrix at setup. The newer option is an unsplit
Convolutional Perfectly Matched Layer (C-PML) in ADE form (Komatitsch &
Martin, 2007).

The C-PML is the **complete two-sided second-order displacement form**
(Matzen 2011; Xie et al. 2014 — the formulation SPECFEM3D's FE-PML
uses). Both occurrences of the spatial derivative in the elastic
momentum equation are coordinate-stretched:

- **Gradient/strain side:** a per-element `psi` ADE memory modifies the
	strain fed to the constitutive law (the stretched stress).
- **Divergence side:** multiplying the momentum equation by `s_j` turns
	the stress-divergence stretch into nodal terms — a mass-proportional
	volumetric damping `rho*d*du/dt` (injected into the mass arrays,
	mirroring the Lysmer discretisation) plus a Complex-Frequency-Shifted
	(CFS) restoring `+rho*d*alpha*u` and displacement-convolution memory
	`-rho*d*alpha^2*zeta`. The CFS pair is evaluated from the element-
	centroid displacement so it cancels exactly at DC and only acts at
	grazing incidence. (A naive "convolve the stress" divergence memory —
	a first-order velocity-stress idea — was tried and rejected: it
	corrupts even the ABC-independent direct arrival under 1-point
	reduced integration. See the git history for the diagnosis.)

Two robustness fixes are layered on: an `alpha` floor and M-PML
cross-axis coupling (Meza-Fajardo & Papageorgiou, 2008). The C-PML
implementation is per-element (each element appears in exactly one
elastic or acoustic mapper), wraps the five absorbing faces (`x_min`,
`x_max`, `y_min`, `y_max`, `z_max`), and leaves `z_min` as the free
surface. The acoustic-pressure PML path carries the same two-sided
treatment (a scalar mirror of the validated elastic form): a gradient-
side `psi` on the pressure Laplacian plus the divergence-side
matched-`kappa` inertia, mass-proportional damping, and centroid-paired
CFS terms. Validated against CPS analytical seismograms the elastic
C-PML is uniformly more accurate than the Lysmer dashpot at every
offset/component (see the table below).

- `abc_type`
	Optional. `dashpot` (default) selects the legacy Lysmer dashpot path
	used historically by Hercules. `cpml` selects the unsplit C-PML path.
	If unset, behavior is unchanged from earlier versions.
- `cpml_thickness_m`
	Target PML thickness in metres. When `> 0`, this value is used directly.
	When `0` or unset (the default), the solver falls back to `cpml_layers ×
	(max element edgesize among elements on an absorbing face)`, which is
	logged at startup so you can confirm the realized thickness. The
	boundary max — rather than the global element max — is used because
	octree meshes refine by local Vs, so the global max is usually a
	coarse interior element that would over-thicken the PML and push its
	inner edge deep into the model; the boundary max makes `cpml_layers`
	mean "this many element layers at the boundary". The face mask is
	computed independently per element from this thickness, so every rank
	must agree on the value; the fallback uses an MPI max-reduce over the
	boundary elements to keep ranks consistent.
- `cpml_layers`
	Optional integer. Default `10`. Used only when `cpml_thickness_m` is `0`.
	Provides a coarse "this many element layers deep" specification without
	requiring the user to know the local mesh resolution at the boundary.
- `cpml_reflection_target`
	Optional. Target theoretical reflection coefficient `R` used in the
	`d_max` formula `d_max = -(p+1) Vp ln(R) / (2 L)`. Default `1.0e-3`.
	Smaller values give stronger absorption but can compress the spatial
	profile and degrade stability margins.
- `cpml_poly_order`
	Optional integer, must be in `[1, 4]`. Default `2`. The polynomial order
	`p` in the spatial profile `d_j(xi) ~ (xi/L)^p`. `p = 2` is the standard
	K&M/SPECFEM3D choice.
- `cpml_kappa_max`
	Optional. Default `1.0` (no real-coordinate stretch). The `kappa`
	stretch (`kappa(xi) = 1 + (kappa_max-1)(xi/L)^p`, 1 at the inner
	edge) improves evanescent/grazing absorption. **It is now
	consistently implemented:** the gradient side stretches strain by
	`1/kappa` *and* the element inertia is scaled by the matched
	`kappa_x*kappa_y*kappa_z` (the leading term of the second-order
	displacement CFS-CPML, Matzen 2011). Without that matched inertia
	the effective impedance would be `rho*c0/sqrt(kappa) != rho*c0`,
	i.e. a spurious inner-edge reflection for `kappa>1`; the matched
	pair restores `rho*c0` (no reflection) and *raises* the local CFL
	margin (heavier elements). `kappa=1` at the inner edge keeps the
	PML↔interior interface exact.
	Multi-axis (edge/corner): `kappa` is **formal-multiaxial** — the
	normalized profile is cross-mixed with the transverse axes (the same
	M-PML mechanism as `cpml_mpml_cross_fraction`, and what the v3
	"formal multiaxial PML" / SPECFEM-class codes ship) and `kappa` is
	derived from the mixed profile, so a corner element's `kappa`
	reflects all active damping directions, consistent with the
	cross-coupled `d` and the matched `kappa`-inertia. Note: the full
	`Pi = s_x s_y s_z` double/triple-delta cross-convolution form is a
	theoretical construct **not implemented by any production seismic
	code (including the v3 reference and SPECFEM3D)**; M-PML
	cross-mixing is the accepted practical multi-axis treatment and is
	what is implemented here.
	History/caveat: an earlier build applied only the `1/kappa` strain
	stretch without the matched inertia, so `kappa_max>1` was
	*inconsistent* (this is why v15 at `kappa=3` only marginally helped
	and perturbed the direct window). With the matched-inertia +
	formal-multiaxial-`kappa` fixes `kappa_max` in `2-3` is expected to
	deliver the grazing benefit, but **still needs re-validation** before
	a firm production recommendation; until then the safest setting is
	`1.0`, and the deep-corner vertical residual is most reliably reduced
	by a **thicker PML** (`cpml_layers` ~30-50, which is what actually
	drove the v12 Z-16 km improvement). Values above `3` remain untested.
- `cpml_alpha_max_hz`
	Optional. When `> 0`, sets `alpha_max = 2 * pi * cpml_alpha_max_hz` in
	rad/s. When `< 0` or unset (default), the solver uses
	`alpha_max = pi * simulation_wave_max_freq_hz`. The standard recipe
	from the PML literature is `alpha_max = pi * f_dominant_source`. Use
	this knob if your source has a peak frequency distinctly different
	from `simulation_wave_max_freq_hz` (e.g. a Ricker pulse at `fc=2.5 Hz`
	in a 1 Hz-resolved mesh) — set `cpml_alpha_max_hz = fc`.
- `cpml_cfl_safety`
	Optional. Default `0.9`. The critical time-step computed from `min h/Vp`
	is multiplied by this factor when `abc_type=cpml` to give a margin
	against the corner-element multi-axis ADE recurrence. Set to `1.0` to
	disable the reduction; values `< 0.5` are likely unnecessary.
- `cpml_mpml_cross_fraction`
	Optional. Default `0.1`. Cross-axis coupling fraction `p` for the
	M-PML stabilisation (Meza-Fajardo & Papageorgiou 2008):
	`d_eff_j = d_j + p * (d_{k1} + d_{k2})`. Setting `p = 0` disables M-PML
	(reverting to plain unsplit CPML, which is known to drift for long
	elastic simulations). The literature-recommended range is `[0.05,
	0.20]`. Lowering `p` reduces the "modification bleed" onto auxiliary
	axes — useful if receivers very close to the PML inner edge show
	extra coda content; raising `p` improves stability margin in
	pathological cases. Valid range `[0, 0.5]`.
- `cpml_alpha_floor_fraction`
	Optional. Default `0.1`. Fraction of `alpha_max` that the per-axis
	`alpha` profile is clamped to from below. The K&M profile
	`alpha(xi) = alpha_max * (1 - xi/L)` reaches zero at the outer edge,
	which destabilises corner elements; flooring at 10% of `alpha_max`
	keeps a non-zero decay rate everywhere in the PML. Set to `0` to
	disable the floor (only safe with explicit M-PML cross-coupling).
	Valid range `[0, 1]`.
- `cpml_pml_material_damping`
	Optional. `none` (default) | `bkt` | `rayleigh`. Selects whether the
	elastic PML zone is **physically attenuating** (viscoelastic
	CFS-CPML) or lossless-elastic. `none` is the historical behavior:
	all material damping skipped on PML elements (bit-identical to the
	pre-feature baseline; this is what every prior `type_of_damping=none`
	CPML validation used). **Important — a previously-hidden critical
	bug:** `type_of_damping=bkt*` **with** `abc_type=cpml` is
	*catastrophically unstable* under `none` — the field diverges to NaN
	by t≈40 s. It was never caught because all earlier CPML validation
	used `type_of_damping=none`, which never exercises the BKT+CPML path.
	`bkt` is the validated fix: the BKT SLS relaxation is carried into
	the PML *stretch-consistently* (native+delta; the relaxed stress
	reduces exactly to the elastic stress at the inner edge, `d=0`, so
	the boundary stays reflectionless). Validated: stable to 240 s, no
	growth, correct physical static offset, and corr **0.994 (Z) /
	0.992 (R)** vs the matching-attenuation CPS analytic. **Resolution
	policy:** with `abc_type=cpml` and an attenuating interior
	(`type_of_damping` BKT-family or `rayleigh`) —
	(1) **unset** → auto-promoted to the matching viscoelastic PML
	(`bkt`/`rayleigh`) with a loud startup log: safe by default, no user
	action; (2) **explicit `none`** → *honored* (not overridden) with a
	loud KNOWN-UNSTABLE warning — a deliberate escape hatch to run the
	raw NaN-prone baseline for attribution/regression (this is how the
	`bkt` instability was proven via v3/v5; it likewise re-enables the
	`rayleigh` attribution control); (3) **explicit `bkt`/`rayleigh`** →
	used as set, with a mismatched value (e.g. `=rayleigh` over a BKT
	interior) still a hard error. So the knob is meaningful again: unset
	= safe auto, explicit `none` = "give me the raw unstable baseline,
	I know what I'm doing", explicit `bkt`/`rayleigh` = force that path.
	`rayleigh` (stiffness-proportional `b*K`, same native+delta path,
	single effective displacement `u_eff=(1+b/Δt)uⁿ−(b/Δt)uⁿ⁻¹`,
	reusing the relaxed-ψ memory) is **implemented (Stage 3)** and
	pending the user's `type_of_damping=rayleigh` validation run; `bkt`
	is fully validated.
	**Behavior change — Rayleigh `b·K` ported to the GPU:** Stage-3 recon
	found stiffness-proportional `b·K` was *never applied on the main
	(non-topography) GPU stiffness path* — `kernelStiffnessCalcLocal`
	used `tm1` only, so `type_of_damping=rayleigh` GPU runs silently lost
	`b·K` (only mass-proportional `a·M` was active). That is now fixed
	(`u_eff` fed to the same operator, matching the CPU conventional and
	GPU-topo paths). Consequence: **existing `type_of_damping=rayleigh`
	GPU results change** — they gain the previously-missing `b·K` and now
	match the CPU path. This is a correctness fix, independent of CPML.
	Every non-Rayleigh run is byte-identical (`c3=0` ⇒ `u_eff≡tm1`), so
	all BKT/none CPML validation (incl. Stage 2) is unaffected.
- `cpml_clamp_outer`
	Optional. `yes` or `no`, default `no`. When `yes`, forces `u = 0` at
	every node sitting on a domain-extent face that is part of the PML
	zone (`x_min`, `x_max`, `y_min`, `y_max`, `z_max`). `z_min` is the
	free surface and is never clamped. This converts the PML problem
	from a round-trip absorption (wave goes out, reflects off the outer
	free surface, comes back through the PML) into a single-pass
	absorption (wave goes out, gets killed at the outer face). Useful
	when very thin PMLs or very long simulations leak round-trip
	reflections. Default `no` preserves the K&M-style free outer face,
	which is what most pre-2010 CPML papers assume.
- `cpml_stiffness_mode`
	Optional. `native` (default) or `reduced`. Selects how the elastic
	PML element stiffness is formed.
	* `native` — the PML element gets the **full native effective
		stiffness** (the same `aTransposeU/firstVector/au` operator as
		interior elements) **plus** a CPML correction `-(sigma_pml -
		sigma_raw)` evaluated at the element centre in the native
		(`mu = 9*c1`, `lambda = 9*c2`) normalization. The correction is
		exactly zero where `d = 0`, so the PML↔interior operator is
		continuous (patch-test exact) and PML elements retain the native
		higher (hourglass-resisting) modes. Recommended.
	* `reduced` — the PML element gets a **standalone 1-point reduced
		CPML stiffness** only (no native higher modes; the pre-`a1c0fb8`
		behaviour). Kept for comparison/debugging. Validated to differ
		from `native` by only ~0.2–0.4% on the homogeneous CPS benchmark
		(the operator mismatch is an O((h/λ)²) higher-order effect), but
		`native` is preferred for coarser/heterogeneous meshes and very
		long runs because it removes the inner-edge operator discontinuity
		and is not rank-deficient. Both modes share the same coefficient
		pipeline and the same divergence-side damping + CFS terms.
- `cpml_stiffness_mode_acoustic`
	Optional. `native` or `reduced`. Independent override of the PML
	stiffness mode for **acoustic** PML elements, decoupled from
	`cpml_stiffness_mode` (which then governs only the elastic PML). If
	unset it **inherits the resolved `cpml_stiffness_mode` value**, so
	existing configs are bit-unchanged. This split exists because the
	elastic and acoustic PML operators have different maturity: it lets
	you run, e.g., elastic `reduced` while exercising the acoustic
	`native+correction` path (or vice versa) without confounding the two.
	The startup log prints the elastic and acoustic PML stiffness modes
	on separate lines so the active configuration is verifiable.
- `topography_bkt_gpu_stabilize`
	Optional, tri-state: `auto` (default), `yes`, or `no`.
	**`auto` (default) = correct-by-default:** the fix auto-enables
	exactly when `computation_method=vt` **and** `type_of_damping=bkt*`
	(the only configuration it applies to), and is a **strict no-op,
	bit-identical** to the validated baseline for every other
	configuration (FEM topography, flat/no-topography, non-VT, non-BKT,
	and every CPML/flat-PML path — the builder hard-guards
	`computation_method=vt` and never builds a mapper or launches the
	new kernel otherwise). `yes` forces it on (VT-guard still applies).
	**`no` is an explicit escape hatch** that forces it fully off,
	restoring the pre-fix behaviour — i.e. the VT+BKT NaN — for A/B
	comparison or bisection; with `no`, *every* configuration including
	non-VT/non-BKT is bit-identical to the pre-fix baseline. Validated
	configs are unaffected by the default flip: the only behaviour change
	vs. the old default is that `vt`+`bkt*` (previously an unconditional
	NaN, no valid result to preserve) now produces the correct
	visco-elastic solution. (No status logging is emitted — the
	stabilization runs silently.) This mirrors the in-tree precedent
	`a36b27e` (auto-enable viscoelastic PML for BKT*+CPML once validated,
	with an explicit opt-out).
	**Pre-existing VT + BKT GPU instability.** With
	`computation_method=vt` (Virtual Topography, Restrepo & Bielak 2014)
	**and** any `type_of_damping=bkt*`, the GPU solver diverges to NaN
	(25-station hill case: blow-up at the hill peak ~t≈2.7 s, NaN by
	≈8.8 s). It is **ABC-independent** (dashpot and cpml alike) and
	**pre-existing** (independent of the C-PML work). Root cause: the GPU
	port never implemented the CPU cut-cell visco-elastic operator. On
	the CPU, VT cut cells flagged `topoBkt` are handled by
	`TetraForcesBKT` — the 5-tet operator with every term scaled by the
	per-tet solid fill fraction. On the GPU the BKT damping kernels have
	only the full-cube operator (plus the PML skip) and no `topoBkt`/
	5-tet branch, while `kernelTopoForceEffective` excludes `topoBkt`
	cells. The cut cell therefore got a **full-cube** force on its
	**fill-fraction-scaled (tiny) lumped mass**; the explicit
	central-difference update divides full-cube force by tiny mass →
	CFL collapse → blow-up.
	When active (default `auto` resolves to on for `vt`+`bkt*`, or
	explicit `yes`), the missing operator is supplied: a faithful CUDA
	port of the CPU `TetraForcesBKT` (`topoTetraForcesBKT` /
	`kernelTopoBktForceEffective` in `psolve.cu`, per-tet algebra copied
	verbatim from the CPU source, both `cube_part` variants, all five
	tets) computes the fill-scaled 5-tet **visco-elastic** force for the
	`topoBkt` cut cells. Its shear/kappa inputs are the **same**
	`damping_vector_shear/kappa` the BKT conv kernel already writes (the
	cut cells stay in the full damping mapper for the conv kernel, so the
	SLS convolution memory is advanced exactly as on the CPU); the
	`mu = rho·Vs²`, `kappa = rho·(Vp²−4/3 Vs²)` moduli are computed
	host-side with the exact CPU expression (bit-identical). The cut
	cells are excluded only from the full-cube force kernel
	(`kernelDampingCalcLocal` is launched with a `topoBkt`-removed
	mapper), so each is processed by exactly one internal-force path.
	GPU `VT`+`BKT` then matches the validated CPU `TetraForcesBKT`
	reference to the usual ported-kernel float tolerance, with Q retained
	in the cut cells.
	**CPML interaction (this branch).** Flat-PML topo cells
	(`cpml_pml_flattened_element`, the SPECFEM-style PML∩topography
	overlap) are **excluded** from the new kernel and left on the
	cube-stiffness + CPML path exactly as on the v2 baseline; the new
	kernel additionally early-exits on any PML-flagged element (mirrors
	the BKT conv/local kernels' PML skip; a strict no-op under dashpot
	where `pmlFlagArrayDevice` is NULL). No `kernel.cu` / `psolve.h` /
	gpu_data_t changes; `solver_init_damping_gpu_mapper` and
	`solver_init_topography_gpu_mapper` are byte-identical to the v2
	baseline; the gate-off path launches the damping kernels with the
	original mappers — so every existing CPML/non-VT/non-BKT
	configuration is bit-identical to the validated v2 baseline by
	construction.
	The stabilization runs silently (no per-rank or mode status
	logging); fatal allocation/CUDA errors still abort with a message.
	Validated (`vt-bkt-cpml`, Vista): `computation_method=vt`,
	`type_of_damping=bkt2` converges under both `abc_type=dashpot` and
	`abc_type=cpml`; the explicit escape hatch `=no` still NaNs
	(pre-fix baseline preserved); FEM+bkt2+cpml on this branch matches
	the latest v2 baseline to within run-to-run non-determinism
	(~1e-10); and the GPU port matches the CPU `TetraForcesBKT` reference
	(`HERCULES_TOPO_FORCE_MODE=cpu`) to the same ~1e-10 / ~1e-9-relative
	floor — i.e. numerically equivalent to the validated CPU operator,
	not merely stable.

Notes on C-PML behavior:

- Per-element classification (acoustic vs elastic) uses the same
	`element_is_acoustic(edata)` test the existing acoustic-pressure path
	uses, so heterogeneous models and acoustic/elastic coupling work without
	extra configuration. Each PML element appears in exactly one of two
	per-type mappers (elastic or acoustic).
- **Material damping on PML elements — `cpml_pml_material_damping`
	dependent.** With the default `cpml_pml_material_damping=none`, all
	material damping is skipped on PML elements when `abc_type=cpml`
	(BKT/Q convolution memory, stiffness-proportional Rayleigh `b*K`,
	the Lysmer dashpot, **and** mass-proportional Rayleigh `a*M` — the
	last was previously applied to PML while the others were skipped, an
	inconsistency that double-counted decay and perturbed the
	divergence-side/CFS DC balance; now uniformly skipped). This was
	believed to be the correct treatment, but it is **only valid when
	the interior is non-attenuating** (`type_of_damping=none`, which is
	what all prior CPML validation used). With an attenuating BKT
	interior the elastic-PML/viscoelastic-interior mismatch makes the
	`none` path **catastrophically unstable (NaN by t≈40 s)** — see the
	`cpml_pml_material_damping` flag entry. The fix is to carry the
	material damping *into* the PML stretch-consistently (viscoelastic
	CFS-CPML), which is **auto-enabled for BKT+CPML**. So: `none` skips
	everything (correct only for a lossless interior, bit-identical
	baseline there); `bkt`/`rayleigh` carry BKT/Rayleigh relaxation
	through the CPML stretch (native+delta, reflectionless at `d=0`).
	Non-PML elements always use material damping as usual. Dashpot mode,
	`BOUNDARY`-undefined builds, and no-attenuation runs are
	bit-identical. Acoustic elements have `zeta=0`, so this is an
	elastic-PML concern.
- The startup log prints
	`CPML: N elastic + M acoustic PML elements globally (thickness X m, alpha_max Y rad/s)`,
	and the derived thickness when `cpml_thickness_m` was zero. Use this to
	confirm the PML zone is non-empty and split between materials the way you
	expect.
- **Topography + C-PML — taper your relief to flat before the PML
	(recommended workflow), do not rely on flatten-PML alone.** The C-PML
	buffer is, by construction, a regular axis-aligned flat-topped layer
	(this is also how SPECFEM3D builds its FE-PML): real free-surface
	relief is not carried into the PML. For `computation_method=vt`
	Hercules has a *flatten-PML re-route* (`cpml_pml_flattened_element`)
	that forces the lateral-PML elements flat so VT+topography+C-PML stays
	stable; the startup log reports how many elements it skipped as
	`Topography force path: gpu (... flattened-pml-skipped=N)`. **This is
	a stability backstop, not a substitute for tapering.** If your
	topography is still tens-to-hundreds of metres high *at the PML inner
	edge*, flattening it there manufactures an artificial step (a
	several-hundred-metre "cliff" for steep boundary relief) that is a
	strong spurious scatterer; the residual is only small/bounded when the
	boundary relief is already near-flat. **Recommended setup for real
	models:** smoothly taper your `topography.in` down to a flat reference
	inside an interior buffer that ends flat — with margin — *before* the
	PML inner edge (flat collar ≥ PML thickness; a multi-km smooth ramp).
	Then the PML and the flatten-PML re-route see no relief and there is
	no cliff. The helper `benchmark_inputs/taper_topography.py` applies
	exactly this 2-D raised-cosine (Tukey) edge taper to any
	`topography.in` (incl. real DEM-derived); the two example generators
	also take `--edge-taper`. (FEM is a strict no-op for flatten-PML;
	non-CPML / non-VT / flat / air runs are byte-identical — the re-route
	never fires. The `benchmark_inputs/topography_tilt_boundary` example
	is a *deliberate boundary-reaching worst-case* to functionally
	exercise the re-route, **not** a model-setup template.)

Using `benchmark_inputs/taper_topography.py` (the recommended
preprocessing step for any real topography + C-PML run):

- **What it is.** A standalone, dependency-free Python 3 script (stdlib
	only — runs anywhere Hercules builds, including the cluster login
	node). It reads a `topography.in`, multiplies the *deviation from a
	flat reference* by a 2-D raised-cosine (Tukey) edge window, and writes
	a new `topography.in` that is **bit-flat in (and well before) the PML
	band** while leaving the deep interior untouched. It does **not**
	change the file format or the solver — `topography.in` in, the same
	format out.
- **Basic use** (defaults target the 40 km × 40 km benchmark domain;
	spacing is read from the file's line 1):

	```
	cd benchmark_inputs
	python3 taper_topography.py path/to/topography.in --out topography_tapered.in
	```

	Then point Hercules at the tapered file — set
	`topography_directory` to the directory containing
	`topography_tapered.in` renamed to `topography.in` (Hercules expects
	exactly that filename), e.g.:

	```
	mkdir -p model/topo_tapered
	python3 taper_topography.py raw/topography.in \
	    --out model/topo_tapered/topography.in --force
	# parameters.in:  topography_directory = .../model/topo_tapered
	```

- **Domain / PML options** (tune to your model, not the 40 km default):
	* `--lx` / `--ly` — domain extent S–N (x) and W–E (y) in metres
		(defaults 40000 / 40000). The script infers the grid from these +
		the file's spacing and **errors out** if the point count does not
		match, so a wrong `--lx/--ly` is caught immediately.
	* `--buffer` — total taper band width inward from each lateral
		boundary (m), = flat collar + ramp (default 6000).
	* `--flat-margin` — width of the **guaranteed-dead-flat collar**
		adjacent to every boundary (default 1500). **Set this ≥ your PML
		thickness with margin** (the startup log line
		`CPML: … (thickness X m …)` gives X; the example domain's ~10-layer
		PML is ≈ 780 m). The script warns if `--flat-margin` is below
		`--pml-thickness`.
	* `--pml-thickness` — PML band thickness (m) used only for that
		safety warning (default 780).
	* `--z-ref` — the flat reference: a number in metres (default `0.0`),
		or `auto` to use the mean elevation of the first interior ring just
		inside the taper band. Use `auto` for real DEMs that sit at a large
		offset (e.g. a plateau at +1500 m) so the taper ramps to a
		representative interior level instead of an arbitrary 0, which
		minimises the tapered amplitude and the ramp's own gradient.
	* `--force` — allow overwriting the output (or in-place when `--out`
		equals the input).
- **What it prints** — the resolved grid, the buffer/collar/ramp widths,
	`z_ref` (and its source), how many points were modified vs left
	untouched, the max elevation change, and a **perimeter-flatness
	self-check** (`max |z_boundary − z_ref|` must be `0` → prints `OK`).
	If that check is not `OK`, do not use the file.
- **`--help`** lists every option with its default
	(`python3 taper_topography.py --help`).
- **Alternative for the two bundled examples.** Instead of running the
	utility you can regenerate either example surface already tapered:
	`python3 gen_topography_tilt.py --edge-taper` (or
	`gen_topography_gaussian_hill.py --edge-taper`), with matching
	`--taper-buffer` / `--taper-flat-margin` / `--taper-zref` options. The
	`--edge-taper` output is bit-identical to running the standalone
	utility on the un-tapered surface (same window math); the flag is
	**off by default** so the committed example `topography.in` files stay
	byte-reproducible.

Stability fixes baked in (tunable via the fractions above):

- **Alpha floor** (`cpml_alpha_floor_fraction`, default `0.1`). The K&M
	linear profile `alpha(xi) = alpha_max * (1 - xi/L)` falls to zero at
	the outer PML edge, which lets corner elements lose all ADE decay
	simultaneously on multiple axes and triggers a fast exponential
	blow-up (this was observed in early testing). The solver clamps
	`alpha` from below at `cpml_alpha_floor_fraction * alpha_max` so
	corner elements always retain some decay rate. No-op at the inner
	edge.
- **M-PML cross-axis coupling** (`cpml_mpml_cross_fraction`, default
	`0.1`). Plain unsplit CPML has a known long-time creeping
	instability for the elastic wave equation even with the alpha
	floor. Following Meza-Fajardo & Papageorgiou (2008), each axis's
	damping rate `d_j` is augmented at coefficient-build time with a `p`
	contribution from the other two axes:

	    d_eff_j = d_j + p * (d_{k1} + d_{k2})

	`d_eff_j` is what enters the `b` and `a` ADE coefficients (the
	constitutive operator `1/kappa_j` is unchanged). The cross-coupling
	restores inter-axis symmetry in the corner and stabilises long
	elastic runs at the cost of slightly more absorption "leaking" onto
	auxiliary axes inside the PML elements.
- **Centroid-paired CFS terms.** The divergence-side CFS restoring
	`+rho*d*alpha*u` and memory `-rho*d*alpha^2*zeta` must cancel at DC.
	Applying the restoring per-node (via the mass arrays) while advancing
	`zeta` from the element-centroid displacement breaks that discrete
	cancellation and grows exponentially at PML corners. Both halves are
	therefore evaluated together from the same element-centroid
	displacement in the PML kernel, which makes the per-element net CFS
	force an exact discrete null at DC. No-op where `d = 0`.
- **Inner-edge half-cell sampling (known bounded residual; centroid kept
	by design).** PML membership and the `(d, kappa, alpha)` profile are
	evaluated at the element **centroid**, so the first PML layer samples
	`xi = h/2` and lumps a small non-zero damping onto the nodes that lie
	exactly on the regular/PML interface, where continuum CPML wants
	`d = 0`. This is the standard, well-characterised cell-centred
	discrete-PML reflection: bounded, second-order, and controlled by the
	**existing knobs** — more `cpml_layers` (smaller `h/L`) or higher
	`cpml_poly_order` shrink it with no code change. The
	`acoustic-pressure-pml-gpu-v3` branch's commit `4a75279` ("shift CPML
	profile samples to the PML inner face", `xi -> xi - h/2`) was
	cross-validated against the CPS analytic solution on the 25-station
	homogeneous elastic half-space: it is **stable** (it stays
	element-level, so it does *not* break the centroid-paired DC balance
	above) and **systematically improves the vertical component** (lower
	L2 and late-coda misfit at 20/20 meaningful stations, ≈ −4 % L2 /
	−7 % late-misfit), but as committed it also renormalises the profile
	to a thinner `active_L = L - h` and **mildly regresses the radial
	component at the grazing far offsets** (12/16 km) by reducing total
	grazing absorption. Net it is a bounded second-order refinement, not a
	correctness fix; the centroid baseline is already validation-grade
	(Z corr 0.992, R 0.990 vs CPS). **Decision: not adopted on this
	branch.** Revisit only for long-duration, very-few-layer, or
	tight-reflection-budget studies — and then as the *compensated*
	face-shift (inner-edge `xi=0` fix **without** the `L-h`
	renormalisation, or `+1–2` layers) with full elastic + acoustic +
	corner-stability re-validation.

CPS analytical validation (elastic half-space, dip-slip point source,
2 s source, 60 s simulation, depth 1500 m). Metric: per-component RMS
misfit vs the CPS analytic, normalised by the CPS RMS — `0` = perfect,
lower = better. The **late coda (25–40 s)** isolates absorbing-boundary
quality (the direct arrival is gone; only boundary reflections remain).
Configuration that produced the table: `abc_type=cpml`,
`cpml_clamp_outer=yes`, a thick PML (`cpml_layers ≈ 50`), and
`cpml_kappa_max=3` in the *earlier inconsistent* form (strain-only
`1/kappa`, no matched inertia). The matched `kappa`-inertia has since
been added; the numbers below are therefore a conservative reference —
re-validation of the consistent `kappa` build (and of the
native+correction stiffness) is pending. The thick PML is the main
driver of the deep-corner Z accuracy here.

Late-coda RMS misfit, near-PML receivers (where ABC quality matters):

| receiver        | Lysmer | CFS-CPML (k=3) | other-branch CPML |
|-----------------|--------|----------------|-------------------|
| R east, 12 km   | 0.222  | **0.107**      | 0.121             |
| R east, 16 km   | 0.649  | **0.153**      | 0.400             |
| Z up,   12 km   | 0.561  | **0.171**      | 0.886             |
| Z up,   16 km   | 2.195  | **1.110**      | 2.735             |

The CFS-CPML is **2–5× more accurate than the Lysmer dashpot** in the
coda at the near-PML receivers, and best of every implementation tested.
At far-field receivers (0–8 km, > 1 PML thickness from the inner edge)
all methods are already accurate and differences are negligible.

In the **direct-arrival window (10–22 s), which is ABC-independent**, the
CFS-CPML matches CPS as closely as Lysmer (e.g. R-east 16 km: 0.151 vs
Lysmer 0.154) — i.e. it does not contaminate the physical wavefield, a
property the rejected stress-convolution variant did not have. Over the
**full record (0–40 s)** it is likewise the most accurate at every
near-PML receiver (e.g. R-east 16 km: 0.170 vs Lysmer 0.366).

**For benchmarking PML quality, keep receivers at least one PML
thickness away from the inner edge** (a receiver embedded in the PML
sees the modified strain field, not a physical seismogram).

Recommended (validated) parameters — this is the configuration behind
the table above:

```
abc_type = cpml
cpml_clamp_outer = yes      # single-pass absorption (Dirichlet outer face)
cpml_layers = 30            # thicker PML closes the deep-corner Z residual
# cpml_kappa_max: 1.0 is the validated setting; 2-3 (now consistently
#   implemented via the matched kappa-inertia) is promising but pending
#   re-validation -- see the cpml_kappa_max entry above
# cpml_stiffness_mode: leave at default 'native' (recommended)
# cpml_alpha_max_hz: leave at default unless source f_dominant != f_max
# cpml_mpml_cross_fraction / cpml_alpha_floor_fraction: leave at 0.1
```

Notes:

- The **validated** configuration is `cpml_kappa_max = 1` (this is what
	the table above and the v9–v14 runs used). `kappa_max > 1` is now
	*consistently* implemented (matched `1/kappa` strain + `kappa`
	inertia) but has not yet been re-validated in the consistent form;
	treat `2-3` as promising-but-unconfirmed.
- The deep-corner vertical (Z, 16 km) residual is **thickness-driven**:
	v12 reached Z-16 km ≈ 1.1 because it used `cpml_layers ≈ 50`, not
	because of `kappa`. For that residual prefer `cpml_layers` 30–50.
- `cpml_clamp_outer = yes` converts round-trip to single-pass
	absorption and is recommended for long runs / thinner PMLs.
- `cpml_layers = 10` already beats the Lysmer dashpot by 2–5× on the
	horizontal components; raise it (30–50) for the deep-corner Z and
	for stronger low-frequency Rayleigh absorption at lateral faces.
- For double-precision correctness checks build with `solver_float =
	double`; production single-precision runs are verified back-to-back.

Typical two-phase settings:

- Phase A (generate mesh database):
	`output_mesh=1`, `mesh_only_mode=1`, `mesh_input_mode=0`
- Phase B (import mesh and solve):
	`output_mesh=0`, `mesh_only_mode=0`, `mesh_input_mode=1`
