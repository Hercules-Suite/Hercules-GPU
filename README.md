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

These are also the current default values used by `run_phaseA_phaseB.sh`
when the corresponding envs are unset.

- `HERCULES_ADJUST_ASSIGN_OMP=1|0`
	Enables or disables OpenMP parallelization in assignment adjust core.
	Default is `1` (enabled).
- `HERCULES_ADJUST_DIST_OMP=1|0`
	Enables or disables OpenMP parallelization in parent-centric distribution
	adjust core.
	Default is `1` (enabled).

- `HERCULES_SPARSE_CLEAR_MODE=auto|range|nodes|full` (or `2|0|1|3`)
	Selects sparse host force reset mode in GPU solver path.
	Default is `full`/`3`.
	`auto`/`2` prefers `full` when topography touches a large fraction of
	local nodes, otherwise falls back to `range` or `nodes`.
	`range`/`0` uses exact contiguous ranges built from the topography node
	set and forces range `memset`.
	`nodes`/`1` can be faster for highly fragmented 3D model cases.
	`full`/`3` clears the full local host force array with one contiguous
	`memset`.
- `HERCULES_SPARSE_CLEAR_OMP=1|0`
	Enables or disables OpenMP parallelization only for sparse host force
	clearing. Default is `0` (disabled) because enabling this path
	reintroduced incorrect plane output in testing. `run_phaseA_phaseB.sh`
	does not auto-set this env.
- `HERCULES_SPARSE_HOST_OMP=1|0`
	Enables or disables OpenMP parallelization in sparse host force staging
	and gather/build loops. Default is `0` (disabled) because this path
	remains the strongest candidate for the high-thread correctness issue
	outside full repro mode.
- `HERCULES_STATION_IO_OMP=1|0`
	Enables or disables OpenMP parallelization for station displacement
	interpolation and per-station file writes. Default is `1` (enabled).
- `HERCULES_PLANE_IO_OMP=1|0`
	Enables or disables OpenMP parallelization for per-strip plane
	interpolation before MPI send/write. Default is `1` (enabled).

Advanced/diagnostic toggles:

- `HERCULES_SHARING_OMP=1|0`
	Enables or disables OpenMP parallelization in shared displacement
	communication pack/unpack paths. Solver default is `0` (disabled)
	because these paths can still make results depend on `OMP_NUM_THREADS`
	even when plane/station/adjust OpenMP toggles are off;
	`run_phaseA_phaseB.sh` does not auto-set this env.
- `HERCULES_CONTRIB_PACK_OMP=1|0`
	Enables or disables OpenMP parallelization while packing contribution
	messages for MPI communication. Solver default is `0` (disabled)
	unless explicitly enabled; `run_phaseA_phaseB.sh` does not auto-set
	this env.
- `HERCULES_REPRO_MODE=0|1|2`
	`0` disables repro mode. `1` disables remaining OpenMP-sensitive
	deterministic communication paths and uses deterministic MPI unpack
	ordering. `2` is stricter and also disables additional CPU OpenMP
	paths guarded by repro mode; `run_phaseA_phaseB.sh` does not auto-set
	this env.

Build note:

- `USE_OPENMP` now defaults to `1` in `quake/forward/Makefile`.
	Use `USE_OPENMP=0` only if you want to disable OpenMP explicitly.
	Examples: `make` (OpenMP on by default) or `make USE_OPENMP=0`.

Two-phase script defaults:

- `run_phaseA_phaseB.sh` now auto-sets the following when they are unset:
	`OMP_PROC_BIND=close`, `OMP_PLACES=cores`,
	`HERCULES_ADJUST_ASSIGN_OMP=1`, `HERCULES_ADJUST_DIST_OMP=1`,
	`HERCULES_SPARSE_CLEAR_MODE=full`, `HERCULES_SPARSE_HOST_OMP=0`,
	`HERCULES_STATION_IO_OMP=1`, `HERCULES_PLANE_IO_OMP=1`,
	These are the current recommended runtime defaults.
- You can still override any of these by exporting your own value before
	running the script.

Two-phase job run (Phase A then Phase B):

- Script: `run_phaseA_phaseB.sh`
- Positional arguments: `PHASE_A_NP PHASE_B_NP [BASE_IN]`

Example (cluster case directory + external GPU binary):

```bash
cd /scratch/08458/wz_tacc/share/hercules/Hercules_Example/Test_case

export OMP_NUM_THREADS=64
# Optional: override script defaults when needed
# export OMP_PROC_BIND=spread
# export OMP_PLACES=threads
# export HERCULES_ADJUST_ASSIGN_OMP=0
# export HERCULES_ADJUST_DIST_OMP=0
# export HERCULES_SPARSE_CLEAR_MODE=auto
# export HERCULES_SPARSE_HOST_OMP=0
# export HERCULES_STATION_IO_OMP=0
# export HERCULES_PLANE_IO_OMP=0
# Advanced/debug overrides:
# export HERCULES_SHARING_OMP=0
# export HERCULES_CONTRIB_PACK_OMP=0
# export HERCULES_REPRO_MODE=1

PSOLVE_BIN=/work/08458/wz_tacc/vista/hercules_GPU_new/quake/forward/psolve \
/work/08458/wz_tacc/vista/hercules_GPU_new/run_phaseA_phaseB.sh \
288 4 \
/scratch/08458/wz_tacc/share/hercules/Hercules_Example/Test_case/inputfiles/Test_merge.in
```

Notes:

- With `BASE_IN`, the script auto-generates temporary Phase A/Phase B inputs
	by overriding `output_mesh`, `mesh_only_mode`, and `mesh_input_mode`.
- Phase B enforces one-rank-per-node via OpenMPI mapping inside the script.

Mesh mode input flags (`numerical.in`):

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
