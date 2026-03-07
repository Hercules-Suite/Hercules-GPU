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
- `HERCULES_REUSE_SOURCE_FILES=1|0`
	When `1`, reuse existing `force_process.<rank>` files from
	`source_directory_output` instead of regenerating source forces.
	Use this for repeated runs with unchanged mesh/source inputs.
- `HERCULES_MESH_ONLY_PREP_SOURCE=1|0`
	When `1` and `mesh_only_mode=1`, run source generation during mesh-only
	phase (Phase A) before exiting.

Build note:

- Compile with `USE_OPENMP=1` to enable OpenMP code paths.

Two-phase job run (Phase A then Phase B):

- Script: `run_phaseA_phaseB.sh`
- Positional arguments: `PHASE_A_NP PHASE_B_NP [BASE_IN]`

Example (cluster case directory + external GPU binary):

```bash
cd /scratch/08458/wz_tacc/share/hercules/Hercules_Example/Test_case

export OMP_NUM_THREADS=64
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export HERCULES_ADJUST_ASSIGN_OMP=1
export HERCULES_ADJUST_DIST_OMP=1
export HERCULES_SPARSE_CLEAR_MODE=auto

PSOLVE_BIN=/work/08458/wz_tacc/vista/hercules_GPU_new/quake/forward/psolve \
/work/08458/wz_tacc/vista/hercules_GPU_new/run_phaseA_phaseB.sh \
288 4 \
/scratch/08458/wz_tacc/share/hercules/Hercules_Example/Test_case/inputfiles/Test_merge.in
```

Notes:

- With `BASE_IN`, the script auto-generates temporary Phase A/Phase B inputs
	by overriding `output_mesh`, `mesh_only_mode`, and `mesh_input_mode`.
- Phase B enforces one-rank-per-node via OpenMPI mapping inside the script.
- Script defaults for source files:
	Phase A sets `HERCULES_MESH_ONLY_PREP_SOURCE=0` by default.
	Script runs a dedicated Phase B prep source substep with Phase B rank layout
	(`mesh_input_mode=1`, `mesh_only_mode=1`) and
	`HERCULES_MESH_ONLY_PREP_SOURCE=1`, `HERCULES_REUSE_SOURCE_FILES=0`.
	Phase B sets `HERCULES_REUSE_SOURCE_FILES=1` (reuse files from Phase A).
	Disable this Phase B prep source substep with `PHASE_SOURCE_PREP=0` if needed.
- Host pinning for stable source reuse:
	By default (`PHASE_PIN_HOSTS=1`), the script pins Source Prep and
	Phase B to a deterministic host list (from `PBS_NODEFILE` or Slurm
	`SLURM_NODELIST`) so rank-to-host mapping stays stable across both
	`mpirun` launches.
	Set `PHASE_PIN_HOSTS=0` to disable auto pinning.
	Set `PHASE_HOSTFILE=/path/to/hosts` to force an explicit hostfile.
- Source file completeness check:
	By default (`PHASE_VERIFY_SOURCE_FILES=1`), the script verifies after
	the Phase B prep source substep that
	`force_process.0..(PHASE_B_NP-1)` all exist under
	`source_directory_output`. This protects mixed-rank runs (for example,
	`PHASE_A_NP=288`, `PHASE_B_NP=4`) from partial source file sets.
	Set `PHASE_VERIFY_SOURCE_FILES=0` to disable this check.

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
