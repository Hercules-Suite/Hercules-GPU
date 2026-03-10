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

- `HERCULES_ADJUST_ASSIGN_OMP=1|0`
	Enables or disables OpenMP parallelization in assignment adjust core.
- `HERCULES_ADJUST_DIST_OMP=1|0`
	Enables or disables OpenMP parallelization in parent-centric distribution
	adjust core.
- `HERCULES_SPARSE_CLEAR_MODE=auto|range|nodes` (or `2|0|1`)
	Selects sparse host force reset mode in GPU solver path.
	`auto`/`2` is default and prefers `nodes` unless range layout is compact.
	`range`/`0` forces range memset; `nodes`/`1` can be faster for highly
	fragmented 3D model cases.
- `HERCULES_STATION_IO_OMP=1|0`
	Enables or disables OpenMP parallelization for station displacement
	interpolation and per-station file writes. Solver default is `0`
	(disabled); `run_phaseA_phaseB.sh` defaults this to `1`.
- `HERCULES_STATION_LAZY_OPEN=1|0`
	When `1`, station files are opened and headers written on first output
	write (lazy open) instead of during startup setup. Solver default is `1`.
	Set to `0` to restore eager station file creation during pre-setup.
- `HERCULES_STATION_SETUP_OMP=1|0`
	Enables or disables OpenMP parallelization in station startup setup
	(ownership scan + interpolation metadata prep). Solver default is `0`
	(disabled); `run_phaseA_phaseB.sh` defaults this to `1`.
- `HERCULES_PLANE_IO_OMP=1|0`
	Enables or disables OpenMP parallelization for per-strip plane
	interpolation before MPI send/write. Solver default is `0`
	(disabled); `run_phaseA_phaseB.sh` defaults this to `1`.
- `HERCULES_PLANE_SETUP_OMP=1|0`
	Enables or disables OpenMP parallelization in plane setup strip
	construction (`io_pool=0` old path). Solver default is `0`
	(disabled); `run_phaseA_phaseB.sh` defaults this to `1`.
- `HERCULES_OCTOR_CHUNK_DIAG=1|0`
	Enables one-line partition diagnostics for Octor chunked transfers
	(max chunk counts and bytes). Solver default is `0` (disabled);
	`run_phaseA_phaseB.sh` defaults this to `1`.

Build note:

- `USE_OPENMP` now defaults to `1` in `quake/forward/Makefile`.
	Use `USE_OPENMP=0` only if you want to disable OpenMP explicitly.
	Examples: `make` (OpenMP on by default) or `make USE_OPENMP=0`.

Two-phase script defaults:

- `run_phaseA_phaseB.sh` now auto-sets the following when they are unset:
	`OMP_PROC_BIND=close`, `OMP_PLACES=cores`,
	`HERCULES_ADJUST_ASSIGN_OMP=1`, `HERCULES_ADJUST_DIST_OMP=1`,
	`HERCULES_SPARSE_CLEAR_MODE=auto`, `HERCULES_STATION_IO_OMP=1`,
	`HERCULES_STATION_LAZY_OPEN=1`, `HERCULES_STATION_SETUP_OMP=1`,
	`HERCULES_PLANE_IO_OMP=1`, `HERCULES_PLANE_SETUP_OMP=1`,
	`HERCULES_OCTOR_CHUNK_DIAG=1`.
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
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export HERCULES_ADJUST_ASSIGN_OMP=1
export HERCULES_ADJUST_DIST_OMP=1
export HERCULES_SPARSE_CLEAR_MODE=auto
export HERCULES_STATION_IO_OMP=1
export HERCULES_STATION_LAZY_OPEN=1
export HERCULES_STATION_SETUP_OMP=1
export HERCULES_PLANE_IO_OMP=1
export HERCULES_PLANE_SETUP_OMP=1
export HERCULES_OCTOR_CHUNK_DIAG=1

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
