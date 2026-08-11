For instructions on how to build, install, configure, and run Hercules, see the
[Hercules GPU Documentation](https://hercules-gpu-docs.liling-dpim.chatgpt.site).
The documentation source is maintained with the solver in [`docs-site/`](docs-site/README.md).

Some general description of the core components in this directory are given as follows.

| File           | Description |
| -------------- | ----------- |
| `octor`        | Hercules mesh generator/partitioner/load balancer module
| `quake`        | Hercules solver module
| `user.mk`      | User-specified compiler and linker settings. This should be the only file we need to change.
| `systemdef.mk` | generic system-specific compile and linker settings. Tailored by entries in user.mk
| `common.mk`    | system-independent settings

The forward solver separates host orchestration from CUDA launch code:
`quake/forward/psolve.cc` is compiled with the configured MPI C++ compiler,
while `quake/forward/psolve_gpu_launch.cu` is compiled with NVCC. Use a parallel
build such as `make -j 8` (or the scheduler-provided CPU count) so these objects
and the remaining solver modules compile concurrently. This also keeps an edit
limited to host orchestration from triggering the previously expensive NVCC
compilation of the entire solver source.

The nonlinear constitutive kernels make PTX assembly the most expensive
forward-solver build phase. Sources under `quake/forward/nonlinear/` therefore
instantiate independent production objects for the classic constitutive
families, Manzari--Dafalias, and PressureDependMultiYield03, plus small runtime
and material-point objects. The runtime selects one family kernel for the
globally configured `material_model`. This keeps the equations normally
inlined, avoids relocatable device code, and lets `make -j 8` compile the three
independent production kernels concurrently. The split also creates stable
boundaries for moving model-specific implementation details without returning
to a single production PTX-assembly job.
Prefer an explicit outer job count such as `-j 8` over bare `-j`, which gives
GNU Make no concurrency limit.

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
	Enables or disables OpenMP parallelization for host-side station setup,
	displacement/kinematic processing, acoustic kinematic reconstruction,
	formatting, and per-station file writes. Default is `1` (enabled). CUDA
	pressure, pore-pressure, strain, and stress interpolation/gather kernels are
	not controlled by this setting.
- `HERCULES_STATION_ROUTE_MODE=interval|all`
	Controls station ownership lookup during setup. Default is `interval`, which uses the octree partition interval to skip stations that cannot belong to the current rank before calling the local octree search. `all` restores the older behavior where every rank searches every station.
- `HERCULES_PLANE_IO_OMP=1|0`
	Enables or disables OpenMP parallelization for plane-strip construction and
	metadata setup in all output backends, reconstructed acoustic-plane
	kinematics on each timestep, and legacy host-side plane interpolation.
	Default is `1` (enabled). This setting does not control CUDA plane
	interpolation, packed-result device-to-host transfer, or file writing.
- `HERCULES_PLANE_IO_MODE=auto|hybrid|mpiio|legacy`
	Selects the plane-output backend. Default is `auto`, which selects legacy CPU
	interpolation followed by the rank-0 strip gather and buffered `fwrite`
	writer. `hybrid` explicitly selects fused GPU interpolation and one packed
	device-to-host transfer followed by that same rank-0 writer. `mpiio`
	selects the experimental distributed MPI-IO writer, in which each rank writes
	its owned strips through a precomputed file view. `legacy` explicitly selects
	the default CPU interpolation and rank-0 path. A configured dedicated
	I/O pool continues to select the existing I/O-pool backend. All backends
	preserve the `planedisplacements.N` binary layout. Plane nodes
	are excluded from the every-step displacement-output patch and are evaluated
	only at `output_planes_print_rate` for the GPU backends. Timers report
	`Plane GPU Interpolation`, `Plane Result D2H`, `Plane Host Patch` (only for
	reconstructed acoustic points), and either `Plane Root Gather` plus
	`Plane Root fwrite` or `Plane MPI-IO Write` inside `Print Planes`/`Solver I/O`.
- `HERCULES_PLANE_MPIIO_BATCH=<snapshots>`
	Sets the number of packed plane snapshots accumulated per rank before each
	collective MPI-IO write. Default is `64`; use `1` to restore one collective
	write per snapshot. A partial final batch is flushed before solver timing is
	reduced and again defensively when the files are closed. Buffer memory scales
	with the sum of the local per-plane payloads times this value; the run log
	reports the selected batch size and maximum buffer allocation per rank.
	Snapshots become visible in the output files when their batch is flushed;
	use `1` when immediate per-snapshot visibility is required. This setting has
	no effect unless `HERCULES_PLANE_IO_MODE=mpiio` is selected explicitly.
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
	With `HERCULES_GPU_MPI_BACKEND=cuda_aware` and a CUDA-aware MPI build,
	`boundary` sends those displacement boundary values directly from device
	memory and performs the hanging-node adjustment on the GPU, eliminating the
	host staging for that exchange.
	Station and plane output no longer force a full displacement host path.
	With the default scalable plane backend, plane interpolation reads
	device-resident displacements directly and does not enter the every-step
	host output patch. The legacy plane backend and station fallback paths patch
	only their interpolation nodes, including the separate tetrahedral station
	stencil used by VT topography. For
	correctness, `boundary` still falls back to `full` when checkpoint/DRM/
	nonlinear/fixed-base/geostatic/incident-plane-wave paths require full host
	displacement arrays. Completed station and plane output files are flushed
	and closed immediately after the solver time loop, before post-solver
	finalization. The run log prints `Displacement sync path:`, the boundary
	node count, and the sparse output patch node count when relevant.
- `HERCULES_PRESSURE_FORCE_SYNC_MODE=sparse|full` (or `1|0`)
	Selects acoustic pressure-force synchronization before the scalar MPI/adjust phase. Default is `sparse`, which gathers only the pressure nodes touched by communication schedules and hanging-node adjustment, exchanges on the existing host schedule path, and scatters only those force entries back to the GPU. `full` restores the legacy full `nharbored` pressure-force D2H/H2D copies.
- `HERCULES_PRESSURE_FIELD_SYNC_MODE=resident|full` (or `1|0`)
	Selects acoustic pressure-field host/device synchronization. Default is `resident`, which keeps pressure history arrays resident on the GPU, swaps their device pointers between steps, and patches only schedule/adjusted pressure nodes back to the current device field after host-side pressure sharing. Station output uses fused GPU interpolation, reconstructed planes use a sparse pressure-node patch, and checkpoint/restart forces the legacy `full` host-visible path.
	When both acoustic pressure synchronization modes are sparse/resident,
	`HERCULES_GPU_MPI_BACKEND=cuda_aware` sends both acoustic pressure arrays
	directly from GPU memory and applies the scalar hanging-node adjustment on
	the GPU. Acoustic station output is interpolated directly on the GPU: one
	fused kernel reads each local eight-node stencil from the two device
	histories and returns only pressure plus its three-component gradient. This
	supports every-step acoustic motion reconstruction without patching station
	nodal pressure arrays on the host. Reconstructed acoustic plane output gathers only
	the local plane pressure nodes needed for host-side pressure gradients, so
	it no longer forces a full acoustic pressure-field host path.
- `HERCULES_PORE_PRESSURE_FORCE_SYNC_MODE=sparse|full` (or `1|0`)
	Selects poroelastic pore-pressure force synchronization before the scalar
	MPI/adjust phase. Default is `sparse`, which gathers only nodes touched by
	communication schedules and hanging-node adjustment, exchanges through the
	existing host schedules, and scatters only those force entries back to the
	GPU. `full` restores the full `nharbored` D2H/H2D copies.
- `HERCULES_PORE_PRESSURE_FIELD_SYNC_MODE=resident|full` (or `1|0`)
	Selects poroelastic pore-pressure field synchronization. Default is
	`resident`, which keeps pore-pressure histories resident on the GPU and
	patches only MPI/adjust nodes after host-side pressure sharing. When station
	pressure output is enabled, a fused GPU interpolator returns one scalar per
	local poroelastic station directly from the resident field and nodal mass
	mask; it does not patch station stencils into a host field.
	`full` restores full `nharbored` field D2H/H2D copies.
- `HERCULES_GPU_MPI_BACKEND=host|cuda_aware`
	Selects the MPI exchange backend for GPU-resident scalar and vector fields.
	The default is `cuda_aware` when the linked Open MPI reports CUDA-buffer
	support through `MPIX_Query_cuda_support()`. This enables direct device
	packing, MPI messages, and GPU hanging-node adjustment for eligible
	poroelastic pressure, acoustic pressure, elastic displacement-boundary, and
	elastic force-contribution exchanges. The solver automatically uses the
	host schedules whenever checkpointing or a host-side correction/output path
	requires them. Set `host` explicitly to use the sparse/resident GPU
	gather-scatter path with the existing host MPI schedules.
	`HERCULES_PORE_PRESSURE_SYNC_BACKEND=host|cuda_aware` remains a legacy alias
	when `HERCULES_GPU_MPI_BACKEND` is unset. `run_phaseA_phaseB.sh` applies a
	Phase A-only `UCX_TLS=sm,self,tcp` default when the active Open MPI reports
	CUDA support and `UCX_TLS` is not already set, because Phase A is CPU meshing
	and large rank-count mesh builds should avoid mlx5 reliable-connection
	endpoint pressure. For Phase B, the wrapper defaults CUDA-aware Open MPI runs
	to `OMPI_MCA_pml=ucx`, unsets `OMPI_MCA_btl`, and sets
	`UCX_TLS=dc,sm,self,cuda_copy,cuda_ipc`. The scalable DC transport avoids the
	per-peer reliable-connection queue-pair pressure observed in large jobs,
	while `cuda_ipc` enables direct intra-node transfers between GPU ranks and
	`cuda_copy` provides the general CUDA-memory fallback. On systems whose UCX
	stack cannot open CUDA IPC handles reliably, set
	`PHASE_B_UCX_TLS=dc,sm,self,cuda_copy` to disable only the IPC transport.
	Unified mode applies this setting to its single MPI job from startup, so it
	also covers the larger Phase-A rank set. Set `PHASE_B_UCX_TLS` or
	`PHASE_B_OMPI_MCA_PML` only when intentionally overriding the Phase B
	default. A user-defined global `UCX_TLS` or `OMPI_MCA_pml` still takes
	precedence for both phases.
	Set `PHASE_B_KEEP_OMPI_MCA_BTL=1` only if a site MPI stack requires an
	explicit inherited BTL setting.

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
- `HERCULES_MESH_IO_MODE=mpiio|sharded|single`
	Selects Phase-A mesh storage. `sharded` lets every MPI rank write an etree
	shard concurrently and publishes a manifest for automatic Phase-B import;
	`mpiio` writes one flat, indexed archive collectively and Phase B reads its
	overlapping record ranges through MPI-IO; `single` retains the serialized
	rank-0 etree writer. Direct `psolve` runs default to `single` for
	compatibility. `run_phaseA_phaseB.sh` defaults to `sharded` in two-phase mode
	and `mpiio` in unified mode. The MPI-IO archive records the native Hercules
	mesh-record ABI and is therefore intended for compatible builds on the same
	platform family; its header rejects incompatible record sizes.
- `HERCULES_PHASE_MODE=two-phase|unified`
	Controls `run_phaseA_phaseB.sh` and defaults to `unified`. `two-phase`
	launches separate Phase-A and Phase-B MPI jobs. `unified` keeps the Phase-A
	MPI world alive, optionally writes the MPI-IO archive according to the input
	file's `output_mesh` value, and redistributes mesh leaves directly from the
	CPU ranks to the requested Phase-B rank count in memory. The Phase-B ranks per
	node are derived from `PHASE_B_NP / allocated_nodes`; no GPU count is
	hard-coded. Unified mode requires `BASE_IN`, an explicit Phase-A rank count,
	and `PHASE_B_OMP_NUM_THREADS` when it cannot be inferred from the allocation.
	On Linux, the selected ranks divide the node's allowed CPU affinity mask into
	disjoint blocks before enabling their Phase-B OpenMP threads.
- `HERCULES_UNIFIED_IDLE_POLL_USEC=100..1000000`
	Sets how often inactive Phase-A ranks enter MPI to test for their Phase-B
	completion message. The default is `10000` microseconds. Between tests they
	sleep instead of blocking in a world barrier, preventing inactive MPI ranks
	from aggressively polling on CPUs assigned to the active GPU ranks.
- `HERCULES_UNIFIED_MPI_YIELD_WHEN_IDLE=0|1`
	Controls the unified launch value of Open MPI's `mpi_yield_when_idle` and
	defaults to `0`. The active Phase-B ranks therefore retain aggressive MPI
	progress; low CPU use applies only to the inactive ranks through Hercules'
	targeted sleep/poll loop. Set `1` only for an intentional diagnostic run.
- `HERCULES_FORCE_GENERIC_DAMPING=0|1`
	Diagnostic/validation override for the viscoelastic damping kernels.
	Default unset/`0`: BKT/BKT2/BKT3/BKT3F use the fast
	hoisted-coefficient convolution kernels (bit-identical to leaving it
	unset). Set to `1` to force the generic per-element
	`kernelDampingCalcConv` for all BKT-family damping, bypassing the
	fast kernels. Use it to produce a generic-kernel reference from an
	otherwise-identical run and verify the fast kernels — they must match
	the generic kernel to single-precision tolerance for a
	spatially-uniform relaxation spectrum (the auto-detected common
	case). A one-time rank-0 line is printed to the monitor/stderr when
	the override is active so a forced run cannot be mistaken for a
	normal one. Read once at startup; production runs leave it unset.
- `HERCULES_SRC_DEBUG=0|1`
	Diagnostic logging for the source-application path (POINT, FORCE,
	SRFH, PLANE, PLANEWITHKINKS). Default unset/`0`: no output. Set to
	`1` to have the owning rank print, per source-application call,
	one `[SRC DEBUG <type>]` line with the source amplitude inputs
	(`mu`, `area`, `maxSlip`, `muArea`, `M0`, source local
	coordinates, `strike/dip/rake`, `forceVec`) plus eight further
	lines giving the 24 entries of the `nodalForce[8][3]` weighting
	array. Useful for validating reciprocity-based Green's-tensor
	workflows: lets you read off the **effective** strike Hercules
	applied after the UTM grid-convergence correction (which differs
	from the user-input `strike.in` value for SRFH). See
	`doc/reciprocity_workflow.md` for the full kernel-construction
	workflow.
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
- `OCTOR_SETCOM_ALL2ALL=0|1`
	Controls construction of Octor's communication manager during mesh
	extraction. The two-phase launcher defaults Phase A to `0`, using Octor's
	neighbor-discovery path; the all-to-all shortcut produced invalid hanging-node
	topology during extraction of a 2,304-rank production mesh. Mesh-input Phase B
	and unified in-memory reconstruction still default to `1` because their
	distributed ownership is already known. Set `1` explicitly for two-phase Phase
	A only after validating the generated mesh at the intended rank count. Unified
	Phase-A generation is unchanged unless this variable is explicitly exported
	for the whole unified launch.
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
# export HERCULES_PLANE_IO_MODE=hybrid
# export HERCULES_PLANE_IO_MODE=mpiio
# export HERCULES_PLANE_IO_MODE=legacy
# export HERCULES_PLANE_MPIIO_BATCH=64
# export HERCULES_SPARSE_UPLOAD_MODE=dense
# export HERCULES_SPARSE_UPLOAD_MODE=range
# export HERCULES_TOPO_FORCE_MODE=cpu
# export HERCULES_PRESSURE_FORCE_SYNC_MODE=full
# export HERCULES_PRESSURE_FIELD_SYNC_MODE=full
# Advanced/debug overrides:
# export HERCULES_CONTRIB_PACK_OMP=0
# export HERCULES_CONTRIB_UNPACK_OMP=0
# export HERCULES_CONTRIB_RECV_MODE=waitall
# export HERCULES_REPRO_MODE=1
# export HERCULES_MESH_IMPORT_VERIFY=1
# export HERCULES_FORCE_GENERIC_DAMPING=1
# export HERCULES_SRC_DEBUG=1
# export PHASE_B_MAP_MODE=numa

PSOLVE_BIN=<path to psolve> \
<path to run_phaseA_phaseB.sh> \
288 4 \
inputfiles/parameter.in
```

Notes:

- With `BASE_IN`, the script auto-generates temporary Phase A, Phase B, and
	unified inputs. The two-phase inputs override `output_mesh`, `mesh_only_mode`,
	and `mesh_input_mode`; the unified input preserves the base input's
	`output_mesh` behavior and sets `mesh_only_mode=0` and `mesh_input_mode=0`.
- `BASE_IN` must define `mesh_etree_output_file` for two-phase runs and for
	unified runs with `output_mesh=1`; it is optional for unified
	`output_mesh=0`. Generated inputs use the same resolved mesh path. Relative
	paths are resolved from the run directory, except when only the
	input-file-relative path already exists.
- When mesh output is enabled, the launcher creates the archive/manifest parent
	directory and prints the resolved `Mesh path`.
- The two-phase script defaults `HERCULES_MESH_IO_MODE=sharded`. In Phase A,
	each MPI rank writes `<mesh_etree_output_file>.partNNNNNN.etree` concurrently,
	then rank 0 publishes a small text manifest at `mesh_etree_output_file`.
	Phase B detects that manifest automatically and each rank reads only the shards
	that overlap its space-filling-curve interval. This removes the serialized
	rank-0 gather/writer from mesh output. Set `HERCULES_MESH_IO_MODE=single` to
	retain the legacy single-etree writer. Keep the manifest and all `.part*.etree`
	files together when moving or archiving a generated mesh.
- Unified mode is the default and avoids the Phase-A-to-Phase-B filesystem
	round trip. The script launches `PHASE_A_NP` ranks once, then the solver
	selects `PHASE_B_NP / allocated_nodes` local ranks per node, redistributes the
	mesh by SFC ownership, and runs the GPU phase on the resulting communicator.
	For example:

	`PHASE_B_OMP_NUM_THREADS=35 PSOLVE_BIN=/path/to/psolve ./run_phaseA_phaseB.sh 576 16 inputfiles/parameter.in`

	In unified mode, `output_mesh=1` writes one MPI-IO archive at
	`mesh_etree_output_file`; its directory/index is embedded in the header, so
	there are no `.part*.etree` files or external manifest. The archive write is
	synchronous and completes before the handoff. Set `output_mesh=0` in the input
	when only the immediate in-memory handoff is needed and no reusable archive
	should be written. Set `HERCULES_PHASE_MODE=two-phase` explicitly to use the
	separate-launch workflow.
	The timing report separates unified handoff setup, record packing, transfer,
	tree import, and reconstruction. Tree import is a subset of transfer and is
	displayed diagnostically without being double-counted in the accounted
	wall-clock total.
- Phase B mapping is controlled by `PHASE_B_MAP_MODE`.
	`auto` is the default and is recommended unless you are intentionally testing an alternative placement.
	For example, on NVL72 with `4` GPU ranks per node and `PHASE_B_OMP_NUM_THREADS=35`,
	`auto` typically uses a node-level mapping like `ppr:4:node:PE=35`,
	while `PHASE_B_MAP_MODE=numa` tries `ppr:2:numa:PE=35` when the local node exposes two NUMA domains.
- For multi-node Phase B runs, `PHASE_B_NP` must divide evenly by the allocated
	node count. This prevents accidental uneven GPU layouts such as 48 Phase B
	ranks on 10 four-GPU nodes; use 40 ranks on 10 GB200 nodes or request 12
	nodes for 48 ranks.

Nonlinear GPU execution (runtime environment):

- `HERCULES_NONLINEAR_BACKEND`
	Selects the nonlinear constitutive-update and internal-force backend.
	`gpu` is the default and uses the CUDA implementation for all current
	nonlinear material models, with double-precision constitutive state and
	device-resident force assembly. `cpu` explicitly selects the established
	CPU/GSL reference implementation. Setting `gpu` explicitly remains useful
	for reproducible launch environments. This setting changes execution, not
	the calibrated values in `parameters.in`. The selected nonlinear backend is
	reported to standard error during solver startup.
- `material_model=ManzariDafalias` is a device-first, pressure-dependent
	effective-stress sand model. Production execution currently requires the GPU
	backend. It uses adaptive explicit Modified Euler and reuses the established
	controls: `backbone_errTol` is the dimensionless local stress-increment error
	target, `no_substeps` caps accepted adaptive substeps, and `error_tolerance` bounds
	the absolute yield-surface residual after a plastic stress correction. The
	elastic predictor uses the one-sided admissibility condition
	`f_trial <= error_tolerance`. A trial increment that crosses from strictly
	inside to outside the yield surface is split at its elastic intersection;
	only the remainder enters Modified Euler. The same admissibility dispatch is
	applied to every adaptive substep: when its first-stage plastic multiplier is
	nonpositive and its elastic predictor remains admissible, Hercules accepts
	the elastic substep without evaluating a discontinuous plastic second stage.
	This avoids a false embedded-error singularity at an unloading/loading switch.
	Because alpha evolves within Modified Euler, every adaptive substep also
	reapplies the reversal-direction test and re-anchors `alpha_in` to the current
	alpha before evaluating the hardening modulus when a new reversal is found.
	This preserves the required positive loading-side hardening gap for a
	multiaxial direction switch inside one global increment.
	A consistency denominator with magnitude below `1e-10` follows the
	UW/OpenSees neutral-loading branch: stress and directional history remain
	unchanged while the local strain is assigned to plastic strain. Any other
	invalid first or second Modified-Euler stage is rejected and retried with a
	smaller fraction before Hercules reports failure.
	Final validation retains that one-sided elastic condition, so a negative
	yield function represents a valid state inside the surface; the absolute
	criterion remains confined to plastic stress correction.
	Elastic stages leave alpha, fabric, and plastic strain unchanged; these history
	variables evolve only under plastic flow.
	Plastic correction is limited to `min(no_substeps, 50)` iterations. Each
	corrected state is transactional: Hercules first tries a genuinely bracketed
	segment projection, then rejects the candidate and halves the constitutive
	substep if the absolute criterion remains unsatisfied. It aborts with GPU
	constitutive failure code 3 only after the minimum local fraction or another
	integration limit is reached, never by continuing from an off-surface or
	partially updated state. The failure report includes the limiting reason,
	retry counts, local error, effective mean stress, element coordinates, and
	the stress and strain-increment tensors. A minimum-local-error failure also
	reports both Modified-Euler raw/accepted plastic multipliers, consistency
	denominators, hardening moduli/gaps, and integration progress so a
	loading-switch discontinuity can be distinguished from an ordinary
	insufficient substep. `manzari_dafalias_recovery` selects the exhaustion
	policy. `strict` is the default and retains finite accepted/rejected-work
	limits without accepting a state above `backbone_errTol`. `force_accept`
	retains those limits but, at the minimum local fraction, commits a finite,
	stress-ratio-capped admissible state in the style of the OpenSees fallback.
	The stress cap is mapped through the elastic compliance so the additive
	strain split is preserved; alpha and alpha-in are recentered on the accepted
	stress while the finite candidate fabric is retained.
	`unlimited_substeps` disables the accepted-substep and rejected-attempt caps
	while retaining the strict error test. `opensees` combines force acceptance,
	uncapped work, and the fixed OpenSees normalized minimum fraction of `1e-6`.
	The other modes use a normalized local-fraction floor that adapts between
	`1e-9` and `1e-6` from the strain roundoff scale. Non-strict modes are
	opt-in robustness fallbacks: they can increase runtime or accept a state
	whose embedded error exceeds `backbone_errTol`, so the startup report names
	the active policy and production results require independent checking. This
	error-control floor
	is distinct from the physical effective-pressure cutoff
	`p_min=1e-4 P_atm`. Predictor and averaged-candidate crossings are first
	resolved by adaptive fraction reduction. At the minimum representable
	fraction, the remaining total strain is completed on an isotropic `p_min`
	state; plastic strain is updated consistently with the elastic stress change,
	and alpha/alpha-in are reset. This prevents tensile skeleton stress without
	returning a partially integrated global increment.
- `material_model=PressureDependMultiYield03` (alias `PDMY03`) is a
	device-first implementation of the standard PDMY03 automatically generated
	nested-surface formulation. It includes pressure-dependent elastic moduli,
	translated yield surfaces, phase-transformation contraction/dilation, and
	post-liquefaction-zone history. The update is explicit and deterministic; it
	does not use the Manzari--Dafalias Modified-Euler error/recovery loop.
	Hercules estimates the requested subincrements from the number of yield
	surfaces crossed and the shear/volumetric strain increment. `no_substeps`
	caps that request. If the request is larger, Hercules integrates exactly
	`no_substeps`, records the clamping diagnostic, and continues rather than
	treating the cap as constitutive nonconvergence. `backbone_errTol`,
	`error_tolerance`, and `manzari_dafalias_recovery` do not control this PDMY03
	integrator. The current interface supports the standard automatic curve with
	2--40 surfaces; custom user-supplied G/Gmax surface coordinates are not yet
	accepted.
	PDMY03 production execution currently requires the GPU backend. Its history
	storage is allocated only when PDMY03 is selected and uses the maximum surface
	count actually present on the local rank. Startup reports the global and
	maximum per-rank nonlinear counts, estimated host/device constitutive storage,
	and maximum-rank timings for classification, host allocation, constant and
	geostatic preparation, and GPU upload. Large cases classify nonlinear
	eligibility once in parallel and compact the cached flags in deterministic
	element order; large FEM cases also prepare independent PDMY03 element
	constants in parallel with the configured OpenMP threads. Startup statistics
	are reduced over the active solver communicator, so unified Phase-A-to-Phase-B
	handoffs do not wait on inactive mesh ranks. FEM
	setup does not perform the VT topography-element lookup; that lookup is needed
	only for the tetrahedral VT formulation. At the elastic-to-plastic transition, Hercules
	matches OpenSees by clamping tensile mean stress to zero and accepting an
	exactly zero state. This prevents a numerical initialization abort but is not
	a substitute for the physical in-situ effective stress required by a
	production soil model. Use `approximate_geostatic_state=yes`, a nonzero ramped
	geostatic loading stage, or a valid restart; otherwise startup warns that the
	physical initialization is absent.
- Manzari--Dafalias rows use 17 values in `material_properties_list`:
	`Vs_limit G0 e_init Mc c lambda_c e0 ksi P_atm m h0 ch nb A0 nd z_max cz`.
	Poisson's ratio and density come from the element's `Vp`, `Vs`, and `rho`.
	All other material families retain the established 16-column layout.
	There are no implicit calibrated Manzari--Dafalias constants: every row must
	provide all 17 values. The constants in the material-point and OpenSees
	fixtures are characterization inputs, not production defaults.
- PressureDependMultiYield03 rows use 22 values in
	`material_properties_list`:

	`Vs_limit Gref Kref friction_angle peak_shear_strain ref_pressure pressure_dependence phase_transform_angle loading_type ca cb cc cd ce da db dc surface_count liquefaction1 liquefaction2 atmospheric_pressure cohesion`

	Stress-like quantities use Pa and both angles use degrees. `liquefaction1`
	and `liquefaction2` are dimensionless OpenSees coefficients and must not be
	rescaled when changing stress units. `loading_type` is 0 for triaxial
	compression, 1 for triaxial extension, or 2 for direct simple shear.
	`surface_count` is an integer from 2 through 40. Continuous values are
	linearly interpolated by element `Vs`; the two categorical fields are selected
	from the nearest `Vs_limit` row. These inputs are case-specific calibration
	parameters and have no production defaults.
	`benchmark_inputs/nonlinear_material_point/opensees_pdmy03_cyclic_shear.tcl`
	provides an independent OpenSees 3.8 affine-brick characterization of the
	synthetic fixture. The companion comparator checks the tensor-shear path and
	all six stresses; the reference development run produced stress NRMSE about
	`1.9e-10` after 100 kPa elastic consolidation and the cyclic plastic path.
	`plot_opensees_pdmy03.py` writes PNG and PDF comparisons of shear response,
	effective mean stress, prescribed strain history, and component error.
- With `approximate_geostatic_state=yes`, Hercules uses the integrated 1-D
	density column stored with each mesh element, rather than a homogeneous
	local-density half-space, and applies `K0 = 1 - sin(phi)`. For an active
	BrickUP element, `alpha*rho_f*g*z_w` is removed from total overburden so the
	constitutive skeleton starts from effective stress. The BrickUP pressure
	unknown remains zero initially because it represents **excess** pore pressure.
	The small kinematic yield surface is centered on that pre-existing stress
	ratio so self weight is not interpreted as a new plastic increment.
	The static stress is subtracted again before dynamic force assembly, so this
	initialization does not inject an unbalanced gravity force. This remains an
	approximation and must be checked against the intended effective overburden
	profile before a production liquefaction analysis.
- Whenever nonlinear analysis overlaps BrickUP poroelastic elements, the Biot
	pressure traction remains explicit and the constitutive model owns the
	skeleton volumetric response. The linear undrained fast path
	(`lambda += alpha^2*M`) is automatically disabled even if requested. This is
	particularly important for Manzari--Dafalias and PDMY03, whose skeleton stress is
	effective stress, but it also avoids bypassing every other nonlinear
	constitutive volumetric update.
- The synthetic fixture includes a standalone OpenSees Modified-Euler driver
	and a characterization comparator under
	`benchmark_inputs/nonlinear_material_point/`. OpenSees runs that reference in
	kPa to preserve its low-pressure regularization scale; the comparator converts
	stress to Pa and does not treat implementation differences as an automatic
	pass/fail validation.
- Two canonical OpenSees characterizations are provided in the same directory.
	The drained monotonic triaxial driver maintains 100 kPa lateral confinement;
	Hercules replays the solved strain path, and radial-stress drift measures
	disagreement in strength and contraction/dilation tendency. Both drivers
	first consolidate with OpenSees' forced-elastic stage. The conversion tool
	seeds the recorded stress, strain, alpha, alpha-in, and fabric at the stage
	boundary and replays only the common elastoplastic history; it does not
	incorrectly pass the elastic consolidation strain through Hercules'
	pressure-dependent plastic stage. The drained driver uses one retained corner
	for all affine face constraints and records SSPbrick's actual center-point
	strain, stress, and material state. Obsolete version-5 output used an invalid
	multiple-retained-node constraint layout, while version 6 used a response
	prefix unsupported by SSPbrick. Version 8 also uses OpenSees' default elastic
	outer tangent because the optional accumulated tangent used by version 7
	became ill-conditioned during early drained plastic loading. Version 9 uses
	a force-residual convergence test for the stress-controlled lateral response
	and independently rejects any monotonic state whose radial material stress
	misses the confinement target by more than 1 Pa; version 8 could accept a
	small displacement correction with about 22 kPa of residual confinement
	error. Older versions are rejected by the current conversion and comparison
	tools. The
	constant-volume cyclic-simple-shear driver reports effective-stress reduction,
	the inferred excess pore pressure `delta_u=p'_0-p'`, `r_u`, and cyclic stress
	reversals. The latter validates constitutive undrained tendency only: it does
	not exercise Hercules U-P pressure degrees of freedom, storage, Darcy flow,
	drainage boundaries, or solid-fluid coupling.
- A coupled one-element benchmark now compares OpenSees
	`SSPbrickUP -lumped` with Hercules BrickUP under sealed, zero-permeability
	undrained triaxial compression. The conversion tool emits
	`material_point_brickup=yes`; the standalone GPU replay then executes the
	production BrickUP pressure-force and pressure-update kernels for the exact
	OpenSees affine strain history. Its independent CPU reference checks
	`p-p0=-M*(ev-ev0)`, with `M=Kf/n`, while the comparison script reports
	effective stress, total stress, pore pressure, `r_u`, confinement, storage,
	and material-state differences. OpenSees's row-summed pressure capacity and
	Hercules's lumped nodal storage therefore use the same representation. This
	isolates the common uniform sealed mode
	without mixing the two programs' global time integrators. The setup keeps
	SSPbrickUP's small `alpha` input (the pressure-field stabilization
	coefficient) distinct from the unit Biot coupling assumed by that element.
	Driver version 4 uses conventional symmetry supports, distributes each face
	traction over four nodes, prescribes the same axial displacement on all four
	top nodes, and rebuilds the transient analysis after opening the pressure
	DOFs. The coupled stage follows the official OpenSees Penalty/Newmark setup
	so the time-varying displacement constraint supplies the solid velocity that
	drives the U-P pressure equation. It also avoids the unsupported combination
	of a prescribed DOF and an `equalDOF` retained node under the Transformation
	handler. The validated version-4 acceptance bounds are effective-stress
	NRMSE at most `5.0e-5`, pore-pressure difference at most `1.0e-6 kPa`, and
	Hercules total radial-confinement error at most `2.0e-1 kPa`; a thresholded
	comparison prints an explicit `PASS`. These values characterize this exact
	100-kPa fixture rather than defining production constitutive tolerances. A
	companion plotting script compares deviatoric response, pore-pressure
	buildup, effective stress path, and confinement drift in PNG and PDF form.
	It does not validate nonuniform pressure, Darcy diffusion, drainage/open-pore
	boundaries, multi-element exchange, or the complete production time-stepping loop. The
		full commands and input-only controls are documented in
		`benchmark_inputs/nonlinear_material_point/README.md`.
- Two additional coupled one-element workflows use OpenSees
		`PressureDependMultiYield03` with `SSPbrickUP -lumped`. Both consolidate to
		100 kPa, close the pressure drainage constraints, enable plastic stage 1,
		and replay the same affine strain path through Hercules PDMY03 and the
		production BrickUP CUDA pressure kernels. One follows a reversing undrained
		triaxial path at constant total radial confinement. The other preserves the
		consolidated lateral strains, maintains total vertical confinement, and
		applies four cycles of x-z simple shear; its OpenSees reference reaches
		`r_u` approximately one. The shared comparator reports effective- and
		total-stress NRMSE, pore-pressure and `r_u` differences, and the relevant
		confinement drift. Acceptance thresholds are intentionally not hard-coded
		until both generated replays are reviewed on the target CUDA/MPI system.
		Driver version 2 explicitly supplies OpenSees 3.8.0's `mType` argument
		between the phase-transformation angle and `ca`. Version-1 reference CSVs
		are rejected because that omitted argument shifted every subsequent PDMY03
		constant and did not represent the intended material.
		Material-point CSV headers are generated from the same result-field table as
		the data rows, including all six PDMY03 surface/substep diagnostics. The
		PDMY03 BrickUP comparator also repairs the short-lived legacy header that
		omitted those names, so an existing successful replay remains usable.
- The short PDMY03 suite fixtures also exercise the BrickUP storage update.
		The cyclic fixture uses the `alpha=1` limit `M=Kf/n`; the zero/tension-onset
		fixture uses `alpha=0.8` and an explicit material-point grain modulus to
		check `1/M=(alpha-n)/Ks+n/Kf`. These are standalone numerical regressions,
		not production calibrations.
- For a new calibration of a model supported by both production backends,
	compare CPU and GPU nonlinear station histories, final yield statistics, and
	force/displacement norms before a production run. Manzari--Dafalias and
	PDMY03 are production-GPU-only; validate each with its shared host/CUDA
	material-point regression and an appropriate independent characterization,
	then run a GPU production case. Keep geostatic, tolerance, substep, stiffness-damping, and
	`material_properties_list` values case-specific.
- With `abc_type=cpml`, solid elements in the absorbing band are always
	linear elastic and are excluded from the CPU and GPU nonlinear constitutive
	maps, even when their `Vs` is covered by `material_properties_list`.
- With `include_acoustic=yes`, acoustic elements are always excluded from the
	CPU and GPU solid constitutive maps, including when the nonlinear property
	table starts at `Vs=0` or an Istanbul payload carries a nonzero `nl_id`.
- Before an element or production run, the standalone material-point mode can
	drive a prescribed tensor-strain history through both constitutive backends
	and fail on any state mismatch:
	`mpirun -np 1 psolve --nonlinear-material-point test.in`. See
	`benchmark_inputs/nonlinear_material_point/README.md` for the test-only
	constant layout and synthetic smoke fixture. A one-point result validates
	correctness, not GPU performance or physical calibration.
- CUDA-aware force exchange and boundary-only displacement synchronization can
	be combined with this path by setting `HERCULES_GPU_MPI_BACKEND=cuda_aware`
	and `HERCULES_DISP_SYNC_MODE=boundary`, subject to the usual run-time safety
	checks.
- The `nl_id` nonlinear-element selector is a legacy Istanbul-model
	classification and is active only with `Istanbul_velocity_model=yes`.
	Non-Istanbul runs ignore `nl_id` and select eligible elements from the `Vs`
	range covered by `material_properties_list`.

Mesh mode input flags (input file):

- `output_mesh`
	`1` writes the generated mesh to disk and is the default when the setting is
	omitted. In unified mode this write completes before the in-memory handoff. With
	`HERCULES_MESH_IO_MODE=sharded`, `mesh_etree_output_file` is the manifest
	path; with `mpiio`, it is the single archive path; with `single`, it is the
	legacy etree database path.
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
- `acoustic_interface_projection`
	Optional. `point` (default) uses the legacy one-point quadrature on each
	acoustic-elastic interface fragment. `face_mass` uses a 2x2 face quadrature
	to weakly project the normal velocity/pressure traction over each fragment
	before assembling nodal interface forces. This is an experimental diagnostic
	option intended to reduce grid-scale interface forcing in
	`acoustic_formulation=velocity_potential` runs. PML interface fragments
	currently retain the `point` path so their one-memory-variable CPML
	correction remains unchanged.
- `acoustic_interface_face_filter`
	Optional. `all` (default) couples every acoustic-elastic interface face.
	`bottom_only` couples only acoustic `Z+` faces, which correspond to the
	water-bottom interface in the depth-positive mesh. This is a diagnostic
	option for bathymetry/shoreline cases: if `bottom_only` is stable while
	`all` is not, the instability is likely coming from lateral wet/dry
	stair-step shoreline coupling rather than from the water-bottom coupling.
	`bathymetry` also couples only acoustic `Z+` faces, but reconstructs a
	local bathymetry normal from nearby acoustic-to-elastic transition depths
	and scales each projected bottom patch by the reconstructed surface area.
	This restores inclined-bottom normal coupling without reintroducing raw
	stair-step side-wall faces.
	The monitor reports candidate and actually coupled face/fragment counts by
	face direction, so this option also shows how much lateral interface area
	was filtered. In `bathymetry` mode it also reports how many bottom faces
	were reconstructed, how many fell back to flat, and the average/maximum
	projected area scale and bathymetry slope.
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
- `acoustic_formulation`
	Optional. Single controlling variable for scalar acoustic physics:
	`velocity_potential` (default), `pressure`, or `none` (`potential` is
	accepted as an alias for `velocity_potential`). `velocity_potential`
	enables the current validated potential formulation (`p = -dphi/dt`), which
	is the recommended acoustic-elastic coupling path. `pressure` enables the
	legacy pressure-variable acoustic solver. `none` is an explicit opt-out:
	it keeps acoustic material handling and interface classification available,
	but does not allocate scalar acoustic DOFs. `pressure` and
	`velocity_potential` require `include_acoustic=yes`; when
	`include_acoustic=no`, the acoustic module is inactive.
	Acoustic elements assemble scalar acoustic mass/stiffness instead of
	displacement mass/stiffness, while elastic elements continue to use
	displacement DOFs. On a flat model surface, acoustic scalar nodes on the top
	external domain face are constrained to enforce a pressure-release water
	surface. When `acoustic_interface_coupling=yes`, scalar acoustic
	formulations build an explicit acoustic-solid interface table. In octree
	meshes, nonconforming interface faces are split into area-weighted fragments
	with fragment centers and normals, and CUDA kernels apply the normal
	traction/acceleration coupling only on those fragments. Each fragment
	interpolates the scalar acoustic field on the acoustic element and
	displacement acceleration on the adjacent elastic element, so refined elastic
	neighbors can couple to coarser acoustic faces. Fragment subdivision follows
	`acoustic_interface_max_level` and finer local neighbors. When the
	neighboring element is local, the table uses that element's mesh-level
	acoustic/elastic classification, so mixed elements kept elastic still
	couple to adjacent acoustic elements. Subdivision has a per-face cap of
	`32 x 32` fragments to avoid accidental table blow-up.
	The legacy `acoustic_pressure_formulation` key is still accepted as a
	backward-compatibility alias only when `acoustic_formulation` is absent:
	`yes` maps to `velocity_potential`, and explicit `no` maps to `none`. New
	input files should use `acoustic_formulation`.
- `acoustic_interface_diagnostics_file`
	Optional acoustic-interface diagnostic output path. When this path is set
	and `acoustic_formulation` is `pressure` or `velocity_potential`, the solver
	writes one row per local acoustic-solid interface fragment after mesh setup.
	The file records the acoustic element id, face id, fragment area, fragment
	center in local and model coordinates, outward normal, element level/size,
	and adjacent elastic element local index/id and local fragment coordinates
	when the adjacent element is local. This file can be very large, so it is
	disabled when unset.
- `print_station_pressures`
	Optional. `yes` writes scalar pressure station histories when
	`acoustic_formulation` is `pressure` or `velocity_potential`, and writes
	BrickUP pore pressure when `include_poroelastic=yes`. Pressure is appended as a
	`Pressure(Pa)` column in the normal `station.N` files after any requested
	displacement, velocity, and acceleration columns. For stations whose
	containing element is acoustic, the motion columns are reconstructed from
	the pressure gradient; for elastic stations they continue to use the
	elastic displacement field. The pressure column is `nan` for stations whose
	containing element has neither acoustic pressure nor poroelastic pore
	pressure. CUDA-aware resident runs report the rank-reduced
	`Acoustic Station Gather` and `Pore-pressure Station Gather`
	kernel-plus-copy timings when those paths are active.
- `reconstruct_acoustic_plane_displacements`
	Optional. `yes` reconstructs displacement at acoustic output-plane points
	from the pressure gradient and writes those reconstructed values into the
	standard `planedisplacements.N` files. Elastic plane points continue to use
	the elastic displacement field, so mixed water/solid planes remain a single
	three-component displacement file. This requires `include_acoustic=yes` and
	`acoustic_formulation=pressure` or `acoustic_formulation=velocity_potential`.
	Default is `no`, which preserves the legacy plane-output behavior.
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
- `acoustic_poroelastic_interface`
	Optional. Default `no_flow`, which keeps an acoustic-poroelastic contact
	sealed for the pore fluid while the existing acoustic pressure traction
	still acts on the solid skeleton. Set to `open_pore` to add weak
	compliance-projected pressure continuity on detected acoustic-poroelastic
	interface fragments, driving pore pressure toward the adjacent acoustic
	pressure. The interface weight uses the adjacent poroelastic element's
	`1/M`, matching the heterogeneous interior pressure capacity. This requires
	`include_acoustic=yes`, the acoustic scalar formulation, and
	`include_poroelastic=yes`.
- `acoustic_poroelastic_pressure_penalty_factor`
	Optional. Default `1.0`. Multiplies the open-pore weak interface compliance
	used by `acoustic_poroelastic_interface=open_pore`. The target exchange
	uses the existing device-resident scalar schedule when CUDA-aware MPI is
	active, with sparse host staging as the fallback.

3D velocity-model quality factors:

- With `3D_velocity_model=yes`, the material directory may contain optional
	`qs.in` and `qp.in` files alongside `vs.in`, `vp.in`, and `rho.in`.
	Each file uses the same count-prefixed layout, contains one value for every
	entry in `z.in`, and follows the same `index.in` column/layer ordering.
	Hercules selects the same vertical layer at each of the four surrounding
	horizontal columns and interpolates `Qs` and `Qp` with the same shape
	functions used for `Vs`, `Vp`, and density. With parametric Q disabled,
	the BKT-family damping setup uses these interpolated values independently
	of the fallback `cvmdb_input_file` mode.
- `qs.in` and `qp.in` must either both be present or both be absent, and their
	value counts must match `z.in`. If both are absent, Hercules preserves the
	legacy behavior `Qs = 0.1*Vs` and `Qp = 2*Qs`. If only one is present or a
	count differs, initialization stops with an explicit error. Because Q values
	are stored in element material data, regenerate or re-correct the mesh after
	adding or changing these files.

Poroelastic U-P input flags (experimental V1):

- `include_poroelastic`
	Optional. `no` (default) keeps the existing elastic/acoustic solver
	behavior. `yes` enables the BrickUP-style U-P poroelastic path for saturated
	shallow solid soil. The skeleton may remain linear or use an eligible
	nonlinear constitutive model; nonlinear use always selects the explicit
	pressure-traction path described below.
- `poroelastic_formulation`
	Optional. Currently only `brick_up` is supported. The implementation
	uses the input `Vp`, `Vs`, and density as the **dry/drained skeleton**
	elastic stiffness and evolves a scalar pore pressure from volumetric strain
	increment plus optional Darcy diffusion. Do not provide saturated/undrained
	`Vp` unless a future velocity-convention conversion is added; otherwise the
	fluid stiffness can be counted twice.
	When the undrained fast path below is inactive, the
		displacement-pressure coupling uses one shared 2x2x2-Gauss BrickUP
		coupling operator, applied as `Q p` in the skeleton force path and as
		`alpha Q^T du` in the pore-pressure source path. The pressure equation is
		assembled in compliance form,
		`C dp = -alpha Q^T du - dt H p`, where each nodal capacity is
		`C_i = sum_e integral(N_i/M_e dV)` and `H` contains Darcy mobility
		`k/mu`. Keeping `1/M_e` inside the element capacity supports neighboring
		elements with different dry/drained moduli and therefore different Biot
		storage moduli; it also avoids moving one element's `M` onto a shared
		node. The GPU path is matrix-free:
	setup precomputes the
	poroelastic element mapper and topography-aware pressure-test weights,
	then CUDA kernels apply the coupling operators and pressure update. With a
	CUDA-aware MPI build, the default `HERCULES_GPU_MPI_BACKEND=cuda_aware`
	keeps eligible pore-pressure schedule exchange and hanging-node adjustment
	on the GPU; `host` restores the sparse host staging path. Pore-pressure
	station values are copied back only at station-output steps.
- `poroelastic_undrained_fast_path`
	Optional. `yes` (default) activates only when
	`poroelastic_surface_pressure_bc = undrained` and all three permeability
	values are exactly `0.0`, and nonlinear analysis is inactive. In that common
	linear validation/undrained limit, the
	solver folds `alpha^2 M` directly into the native elastic volumetric
	stiffness (`lambda += alpha^2 M`) and skips the mechanical `Q p`
	pressure-traction feedback. Pore pressure is still updated from `Q^T du`
	and remains available in station output as a diagnostic. Set this to `no`
	to use the older split `Q p` / `Q^T du` path for comparisons, or use
	nonzero permeability / `drained` pressure BC for the full split path.
	`acoustic_poroelastic_interface=open_pore` also disables this fast path so
	the pressure-continuity pore pressure feeds back mechanically through the
	split `Q p` term. Nonlinear analysis likewise disables it automatically: the
	configured value remains `yes`, but the effective runtime path is the
	explicit `Q p` coupling needed by the constitutive skeleton.
- `poroelastic_depth_max_m`
		Optional. Default `30.0`. Only solid, non-acoustic, non-PML elements
		within this depth below the local ground surface are eligible for
		poroelastic classification. When nonlinear analysis is enabled, startup
		prints a `Nonlinear/BrickUP map` summary showing how many nonlinear
		elements are active BrickUP elements and how many are above groundwater,
		below this depth limit, in CPML, or otherwise unmapped. Nonlinear elements
		below the limit remain single-phase nonlinear solids: their constitutive
		stress is assembled directly as total stress, with no pore-pressure degree
		of freedom or Biot feedback, and approximate geostatic initialization uses
		total overburden. This is an intentional dry/drained or calibrated
		total-stress fallback, not a coupled liquefaction calculation. Increase the
		depth limit whenever deeper saturated nonlinear soil is intended to build
		excess pore pressure.
- `surface_groundwatertable`
	Optional. `no` (default) means the groundwater level is the local
	ground surface, so eligible shallow soil is saturated from the
	surface down. `yes` reads a groundwater-level grid and only classifies
	elements whose centroids are below that level.
- `groundwaterlevel_file`
	Optional. Used only when `surface_groundwatertable=yes`. If unset,
	the solver reads `groundwaterlevel.in` from the 3D material model
	directory (`3D_velocity_model_directory`, normally
	`inputfiles/materialfiles`). If set to a bare filename, it is resolved
	relative to that materialfiles directory; if it contains `/`, it is
	treated as an explicit path. The file format mirrors `topography.in`:
	first value is grid spacing, followed by regular-grid values over the
	full model footprint in the same ordering as `topography.in`. Values
	are always positive depth below the local ground surface, so `100.0`
	means the water table is 100 m below ground everywhere.
- `poroelastic_surface_pressure_bc`
		Optional. `undrained` (default) leaves pore pressure unconstrained at
		the saturated surface. `drained` applies a weak, compliance-projected
		Dirichlet condition on the groundwater/free-drainage plane: the solver
		assembles a nodal drained-plane capacity from a 2x2 quadrature rule on the
		water-table cut plane inside the first saturated element, then blends
		the free pressure update toward zero using that capacity. The plane term
		uses the adjacent element's `1/M`, so it has the same units as the
		interior pressure capacity even when `M` varies spatially. This avoids the
		old nearest-node hard clamp when the water table cuts through an
		element.
- `poroelastic_surface_pressure_penalty_factor`
		Optional. Default `1.0`. Scales the weak drained-plane capacity used by
	`poroelastic_surface_pressure_bc = drained`; larger values enforce
	closer-to-zero pore pressure on the plane, while `0.0` disables the
	weak projection without changing the classification.
- `poroelastic_elastic_interface`
	Optional. Currently only `no_flow` is supported. Poroelastic and
	elastic regions share displacement nodes mechanically, while pore
	pressure exists only on poroelastic elements/nodes; Darcy diffusion is
	assembled only through poroelastic elements, giving a no-flow
	poroelastic-elastic boundary in V1.
- Poroelastic material constants:
	`poroelastic_porosity` (default `0.4`),
	`poroelastic_biot_alpha` (default `1.0`; values below `1.0` derive
	a local grain bulk modulus from the dry/drained skeleton relation
	`K_s = K_d / (1 - alpha)` and use
	`1/M = (alpha - porosity) / K_s + porosity / K_f`; `alpha = 1.0`
	uses the incompressible-grain limit `M = K_f / porosity`),
	`poroelastic_fluid_bulk_modulus_pa` (default `2.2e9`),
	`poroelastic_fluid_density_kgm3` (default `1000.0`),
	`poroelastic_fluid_viscosity_pa_s` (default `1.0e-3`), and
	`poroelastic_permeability_x_m2`,
	`poroelastic_permeability_y_m2`,
		`poroelastic_permeability_z_m2` (all default `0.0`, i.e. undrained
		no-diffusion response).

For nonlinear BrickUP runs, a positive `geostatic_loading_time_sec` with all
three permeability values exactly zero is rejected during startup. Such a
gravity ramp is undrained: it can place the applied overburden into the pore
pressure unknown instead of establishing the intended effective skeleton
stress. Use `approximate_geostatic_state = yes` with zero loading and cushion
times, or use a physically calibrated permeability and consolidation interval.
The approximate-state path initializes the constitutive effective stress and
keeps the pressure unknown as excess pressure.

CPML policy for poroelasticity: PML elements are kept acoustic or
elastic only. Any saturated poroelastic candidate inside the CPML zone is
left as elastic CPML and reported in the startup line
`Poroelastic BrickUP: ... saturated CPML candidates left as elastic CPML`,
matching the practical SPECFEM-style policy that avoids poroelastic PML in V1.

Scalar acoustic formulation — correctness fix & known limitation:

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
	at the absorbing boundary and logs the realized thickness to standard error
	at startup. The PML element counts and optional outer-clamp count use the
	same startup stream.
- `cpml_layers`
	Optional integer. Default `10`. Used only when `cpml_thickness_m` is
	`0` — "this many element layers deep" at the boundary, without needing
	to know the local mesh resolution.
- `cpml_partition_weight`
	Optional. Default `1.0`, which preserves the standard equal-leaf-count
	octree partition. Values `> 1.0` make the partitioner count future
	CPML elements as heavier work, so boundary-heavy ranks may receive
	fewer total elements. This is a performance tuning knob only: it does
	not refine/remesh the model, does not change `mesh.e`, and can be
	changed for Phase B runs that import an existing mesh. Use it only when
	the `CPML rank/load diagnostics` show that the same ranks have high
	`pml_total` and high `core_wait`; otherwise leave it at `1.0`.
	Moderate values (`1.25`-`1.5`) are the useful range seen so far.
	Larger values can reduce CPML compute skew but also increase MPI /
	patch imbalance, so they may slow the full run even when the PML
	kernel load looks better.
- `cpml_3d_material_extension`
	Optional. `yes` (default) or `no`. Used only with `abc_type=cpml`
	and `3D_velocity_model=yes`. When enabled, material queries made
	from inside the CPML are clamped to the inner PML box before reading
	the 3D velocity model, so Vs/Vp/rho and the acoustic/elastic water
	classification are extended normally into the absorber instead of
	letting a coarse `materialfiles` grid create a new coastline or
	bathymetry inside the PML. Corners clamp along every active PML
	normal. This changes sampled element material properties, so rebuild
	or re-correct an imported etree for the PML element values to change;
	Phase B imports that reuse baked properties keep those existing
	values. Set `no` to restore the legacy raw 3D material queries.
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
	Optional. `none` (default) | `rayleigh`. Whether material damping is
	applied inside elastic PML elements. `none` means material damping is
	skipped in the PML zone; non-PML elements still use the requested
	`type_of_damping`. For BKT-family runs (`bkt`, `bkt2`, `bkt3`,
	`bkt3f`), this is the supported policy: the interior uses BKT damping,
	while PML elements use native elastic stiffness plus the elastic CPML
	correction. The previous experimental BKT-in-PML path was removed
	because stretching the BKT memory variables inside CPML overdamped
	large heterogeneous runs to numerical zero. Therefore `bkt`, `bkt2`,
	`bkt3`, and `bkt3f` are not valid values for
	`cpml_pml_material_damping`; use `none` with BKT-family
	`type_of_damping`.

	`rayleigh` enables the Rayleigh material-damping CPML path and is valid
	only with `type_of_damping=rayleigh`. If `type_of_damping=rayleigh`,
	`abc_type=cpml`, and `cpml_pml_material_damping` is omitted, the solver
	auto-promotes it to `rayleigh`; explicit `none` is honored for
	debugging/regression with a warning.
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
	For `3D_velocity_model=yes`, `cpml_3d_material_extension=yes`
	(default) makes that classification come from the corresponding
	inner-PML material column/row, matching the usual flat-PML material
	extension used by spectral-element PML setups.
- Material damping on PML elements follows `cpml_pml_material_damping`
	(see that flag). With the default `none`, material damping is skipped
	on PML elements, but non-PML elements always use material damping as
	usual. This is the supported policy for BKT-family damping: the
	interior is viscoelastic, while the PML zone is elastic CPML.
	`rayleigh` is the only supported material-damping mode inside PML
	elements.
- **Viscoelastic damping runs on the GPU for all BKT-family models**
	(`bkt`, `bkt2`, `bkt3`, `bkt3f`). When the relaxation spectrum is
	spatially uniform — the usual case: fixed relaxation frequencies,
	per-element weights vary with Q — the solver uses fast convolution
	kernels with the spectral coefficients hoisted to scalars (computed
	once on the host, no per-element/-node `exp()`, no `g`-array
	traffic): `bkt`/`bkt2` share the 2-mechanism kernel (they differ
	only in the time-scale factor), `bkt3`/`bkt3f` use a 3-mechanism
	kernel. If the relaxation spectrum is spatially **heterogeneous**,
	the solver auto-detects this at setup (it scans every damping
	element) and falls back to a generic per-element kernel that is
	numerically identical (bit-for-bit for the uniform case). Both
	paths feed the interior BKT force correctly; BKT material damping is
	not applied inside PML elements. **Behavior change:**
	a heterogeneous-spectrum `bkt2` model now uses the generic kernel —
	previously it silently applied the first element's coefficients to
	all elements (a latent bug). Spatially-uniform `bkt2` (every
	validated run) is bit-identical to before and unaffected.
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
	This surface taper is independent of `cpml_3d_material_extension`:
	tapering makes the PML top surface flat, while
	`cpml_3d_material_extension` extends the 3D material/water mask
	through the PML volume.

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
# cpml_partition_weight:    leave at 1.0 unless CPML rank/load diagnostics justify tuning
```

Notes:

- The validated configuration is `cpml_kappa_max = 1`; the deep-corner
	vertical residual is reduced most reliably by a thicker PML
	(`cpml_layers` 30–50), not by `kappa`.
- `cpml_clamp_outer = yes` is recommended for long runs / thinner PMLs.
- `cpml_layers = 10` already beats the dashpot 2–5× on the horizontal
	components; 30–50 helps the deep-corner vertical and low-frequency
	absorption.
- `cpml_partition_weight` affects only how elements are assigned to MPI
	ranks. It is useful for diagnostics/tuning, but the default `1.0` is
	the safest production setting unless a specific case shows CPML-heavy
	ranks pacing the solve.
- For double-precision correctness checks build with `solver_float =
	double`; production single-precision runs are verified back-to-back.


Single-point body-force source (`type_of_source = force`):

In addition to the moment-tensor `point`, `plane`, `planewithkinks`, and
`srfh` source types, Hercules supports a single-point body-force source
analogous to SPECFEM3D's `FORCESOLUTION`. It applies a Cartesian force
`F(t) = (Fx, Fy, Fz) * stf(t)` at a single point inside one element,
distributed to the 8 nodes via the trilinear hex8 shape functions
`f_a^i = F_i * N_a(x,y,z)`. No moment tensor, strike, dip, or rake is
needed.

- `print_station_strains` / `print_station_stresses`
	Optional, independent (`yes`/`no`, default `no` → byte-identical to
	before). For **linear** analysis, append the 6 strain and/or 6 stress
	tensor components per station per output step — the linear-analysis
	counterpart of the strain/stress columns that `include_nonlinear_
	analysis=yes` already produces. Intended for reciprocity-based fault
	inversion (strain Green's tensors). Each is a **contiguous 6-column
	block** in the order `XX, YY, ZZ, XY, YZ, XZ`, appended to the normal
	`station.N` files after any displacement/velocity/acceleration/
	pressure columns (header labels `Eps_*` / `Sig_*`); if both flags are
	set, the strain block precedes the stress block. **Tensor strain
	convention:** `eps_ij = ½(∂u_i/∂x_j + ∂u_j/∂x_i)` for the
	off-diagonals (the true strain tensor, matching Hercules' existing
	nonlinear output — *not* engineering shear). **Strain (`Eps_*`)** is
	always the kinematic `ε = ½(∇u + ∇uᵀ)` of the current displacement
	`u(t)` — correct for both elastic and viscoelastic runs. In BKT-family
	runs the station gather kernel reads only the owning element's eight
	nodes from the device-resident displacement and element-node arrays,
	then applies the regular-hex or VT-tet strain interpolation directly
	on the GPU. The host `point_strain` / `compute_tetra_strain` path remains
	the fallback for non-BKT damping.
	**Stress (`Sig_*`) is *physically correct* under viscoelastic
	damping (BKT-family).** When `type_of_damping` is BKT/BKT2/BKT3/
	BKT3F, a small GPU gather kernel (`kernelStationVEStressGather`)
	computes the viscoelastic stress from Hercules' internal
	modified-displacement vectors using the same deviatoric/volumetric
	split the damping force operator uses internally:
	`σ_visco = 2μ·dev(ε_shearVector) + K·tr(ε_kappaVector)·I` with
	`K = λ + 2μ/3` (so off-diagonals are pure deviatoric `2μ·ε^s_ij`
	from `shearVector`, and the diagonal carries both the shear-channel
	and the bulk-channel SLS memory). The GPU kernel and its host transfer
	are request-aware: strain-only and stress-only output compute and copy
	only one packed 6-component tensor per station, while requesting both
	uses 12 packed values. One small `cudaMemcpyAsync` per station-output
	step ships those tensors to host; no bulk transfer of the full per-node
	arrays and no redundant host strain reconstruction. The cumulative
	kernel-plus-copy wait is reported separately as `Station Tensor Gather`
	inside the `Stations` and `Solver I/O` timers. Under non-BKT damping
	(`none`/`rayleigh`/`mass`) the GPU gather is skipped and `Sig_*`
	falls back to the elastic constitutive `σ = 2μ ε + λ tr(ε) I` on
	the kinematic strain (which IS the physical stress in that regime).
	**Caveat — one-step phase lag in `Sig_*`:** the per-step kick reads
	the previous step's `shearVector`/`kappaVector` (station output runs
	at the top of each timestep, before that step's damping conv
	kernel), so `Sig_*` is offset by one simulation `Δt` relative to
	`Eps_*`. Negligible for the typical reciprocity / kinematic-
	inversion sample rate and easily corrected post-hoc if needed.
	Acoustic/fluid stations (which have `print_station_pressures`
	instead) get no strain/stress columns. The GPU plumbing is fully
	gated: with both flags off and/or non-BKT damping, no extra device
	allocation, no kernel launch, no memcpy — byte-identical to before.

Use cases:

- **Reciprocity-based Green's tensor computation** for finite-fault slip
	inversion. Three runs with unit forces along x, y, z at the receiver
	location give the three strain Green's tensors needed to map an
	arbitrary moment-tensor source on the fault to the receiver.
- **Single-impulse studies** (e.g., hammer source, acoustic transducer
	calibration) where a body force is the physical excitation.

`source.in` keys (see `examples/test1/sourceforce/source.in` for a full
example):

```
type_of_source       = force
source_function_type = quadratic     # ramp | sine | quadratic | ricker | exponential | discrete
average_risetime_sec = 0.5           # (or ricker_Ts/ricker_Tp for ricker)
lonlat_or_cartesian  = 0             # 0 = lon/lat, 1 = cartesian
hypocenter_lat_deg   = 34.0
hypocenter_long_deg  = -118.5
hypocenter_depth_m   = 1000
force_x_n            = 0.0           # F_x in Newtons (north)
force_y_n            = 0.0           # F_y in Newtons (east)
force_z_n            = 1.0e15        # F_z in Newtons (downwards)
```

The source-time-function (`source_function_type`, `average_risetime_sec`,
`ricker_Ts`, `ricker_Tp`) keys are shared with the moment-tensor point
source. Force components are in Newtons in the simulation Cartesian
frame (`x = north`, `y = east`, `z = downwards`). Internally the FORCE
path reuses the same file-based per-step force pipeline as POINT — no
new MPI traffic, no new GPU plumbing.

`source_function_type = discrete` (tabulated waveform). For non-analytical
STFs — e.g., a recorded velocity pulse used as the excitation for a
reciprocity Green's function — drop a `forcefunction.in` file alongside
`source.in` in the source directory with this format:

```
<nt1>             # number of samples (int)
<tinit>           # initial delay before the pulse begins (s)
<dt>              # sampling interval (s)
<s[0]>            # nt1 unit-amplitude STF samples (double, one per line)
...
<s[nt1-1]>
```

The instantaneous force is `F_vec * interp(t - tinit, s, dt)`. Normalize
the samples to peak 1 and carry the physical amplitude in
`force_{x,y,z}_n`. Internally this reuses the existing `interpolate_linear`
path that SRFH already uses for `slipfunction.in`, so it is bit-identical
to the SRFH discrete-STF math.

Distributed compliant-base force source (`type_of_source = base_force`):

For OpenSees-style 1D site-response validation, Hercules also supports a
distributed base force source. It reads a base velocity history and applies
nodal forces on the deepest horizontal mesh-node plane:

```
F_i(t) = scale * rho_base * Vs_base * A_i * v(t) * direction_i
```

where `A_i` is the bottom-face tributary area assembled from local hex
faces. This is the force-source analogue of the Lysmer-Kuhlemeyer
compliant-base input often used in OpenSees. It precomputes the standard
`force_process.*` files during source initialization, so the solver hot path
is unchanged.

`source.in` keys:

```
type_of_source       = base_force
source_function_type = discrete

base_force_rho_kgm3      = 2500.0   # base/rock density
base_force_vs_mps        = 700.0    # base/rock shear velocity
base_force_direction_x   = 1.0      # x = north
base_force_direction_y   = 0.0      # y = east
base_force_direction_z   = 0.0      # z = down
base_force_scale         = 2.0      # use 2.0 for outcrop velocity input

# Optional. Omit to use the deepest mesh-node plane.
# base_force_depth_m     = 30000.0
# base_force_tolerance_m = 0.01
```

Place `baseforcefunction.in` beside `source.in`:

```
<nt1>             # number of velocity samples
<tinit>           # initial delay (s)
<dt>              # sampling interval (s)
<v[0]>            # base velocity sample in m/s
...
<v[nt1-1]>
```

Use `base_force_scale = 1.0` for an incident/within velocity history and
`2.0` for an outcrop velocity history, matching the usual compliant-base
conversion convention. The first implementation is intended for flat-bottom
validation columns; if `base_force_depth_m` is provided and no bottom faces
are found within `base_force_tolerance_m`, source generation aborts clearly.
