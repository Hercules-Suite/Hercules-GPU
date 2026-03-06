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

PSOLVE_BIN=/work/08458/wz_tacc/vista/hercules_GPU_new/quake/forward/psolve \
/work/08458/wz_tacc/vista/hercules_GPU_new/run_phaseA_phaseB.sh \
288 4 \
/scratch/08458/wz_tacc/share/hercules/Hercules_Example/Test_case/inputfiles/Test_merge.in
```

Notes:

- With `BASE_IN`, the script auto-generates temporary Phase A/Phase B inputs
	by overriding `output_mesh`, `mesh_only_mode`, and `mesh_input_mode`.
- Phase B enforces one-rank-per-node via OpenMPI mapping inside the script.
