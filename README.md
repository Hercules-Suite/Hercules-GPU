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

- `HERCULES_SAFE_MODE=1`
	Disables aggressive `compute_adjust` optimizations and OpenMP paths for
	robustness checks.
- `HERCULES_ADJUST_ASSIGN_OMP=1|0`
	Enables or disables OpenMP parallelization in assignment adjust core.
- `HERCULES_ADJUST_DIST_OMP=1|0`
	Enables or disables OpenMP parallelization in parent-centric distribution
	adjust core.
- `HERCULES_ADJUST_DIST_PARENTCENTRIC=1|0`
	Enables parent-centric distribution algorithm (`1`, default) or falls back to
	legacy dangling-centric distribution (`0`).

Build note:

- Compile with `USE_OPENMP=1` to enable OpenMP code paths.
