# Hercules GPU — Extended README

This is the comprehensive engineering reference: runtime controls, numerical
formulations, implementation details, input/output examples, and development
and validation history. For a concise operational guide and parameter index,
see [Technical README](Technical%20README.md). The
[documentation site](https://hercules-gpu-docs.liling-dpim.chatgpt.site) provides
the full input reference; its [source](../docs-site/README.md) is maintained
with the solver. Current parser and solver code under `quake/forward/` take
precedence over historical notes.

Use one `parameters.in` for simulation controls and `source.in` for
source-specific controls. Units are SI; solver X is north, Y is east, and
Z/depth is positive downward. Nonlinear examples require case-specific
calibration of geostatic state, tolerances, substeps, damping, and material
properties. Reported historical benchmarks apply to their stated cases and
versions; consolidating these documents does not validate new combinations.

## Contents

- [Domain sizing](#octree-compatible-domain-sizing) and [compilation](#host-and-cuda-compilation)
- [Runtime controls](#runtime-performance-controls) and [diagnostics](#advanced-and-diagnostic-controls)
- [Launcher defaults](#phase-launcher-defaults), [launching](#launch-a-simulation), and [mesh import](#mesh-generation-and-import)
- [Nonlinear constitutive models](#nonlinear-gpu-execution-and-constitutive-models)
- [Acoustic coupling](#acoustic-and-elastic-coupling), [basin grids](#basin-bottom-grid), and [3D material files](#3d-velocity-model-files)
- [Poroelasticity](#poroelastic-u-p-formulation) and [acoustic limitations](#scalar-acoustic-corrections-and-limitations)
- [Absorbing boundaries](#absorbing-boundaries-and-cpml)
- [Point forces](#single-point-body-force-source) and [compliant-base forces](#distributed-compliant-base-force-source)
- [Additional derivations and historical results](#additional-derivations-and-historical-results)
- [Station binary output](#single-file-station-binary-output)
- [Structured DRM box output](#structured-drm-box-output)

## Repository components

| File           | Description |
| -------------- | ----------- |
| `octor`        | Hercules mesh generator/partitioner/load balancer module
| `quake`        | Hercules solver module
| `user.mk`      | User-specified compiler and linker settings. This should be the only file we need to change.
| `systemdef.mk` | generic system-specific compile and linker settings. Tailored by entries in user.mk
| `common.mk`    | system-independent settings

## Octree-compatible domain sizing

Hercules domains must be compatible with the dyadic octree grid. Define

```math
\begin{aligned}
L_{\mathrm{north}} &= \mathtt{region\_length\_north\_m}, \\
L_{\mathrm{east}}  &= \mathtt{region\_length\_east\_m}, \\
L_{\mathrm{depth}} &= \mathtt{region\_depth\_deep\_m}
                     -\mathtt{region\_depth\_shallow\_m}, \\
L_{\max} &= \max(L_{\mathrm{north}},L_{\mathrm{east}},L_{\mathrm{depth}}).
\end{aligned}
```

At a common octree level $`N`$, every physical extent $`D`$ must be an integer
multiple of $`\Delta_N=L_{\max}/2^N`$. Equivalently,

```math
D=m\frac{L_{\max}}{2^N},\qquad m\in\mathbb{Z},\quad 0<D\le L_{\max},
```

or

```math
D=\sum_{i=0}^{N}a_i\frac{L_{\max}}{2^i},\qquad a_i\in\{0,1\}.
```

Use exact whole-metre extents with the current Octor path, which reduces the
physical dimensions using integer metre values.

For $`L_{\max}=\mathtt{region\_length\_east\_m}=1{,}000{,}000\,\mathrm{m}`$,
the dyadic sequence begins $`1{,}000{,}000`$, $`500{,}000`$, $`250{,}000`$,
$`125{,}000`$, $`62{,}500`$, and $`31{,}250\,\mathrm{m}`$. Thus
$`\mathtt{region\_length\_north\_m}=750{,}000\,\mathrm{m}`$ is compatible
because $`750{,}000=500{,}000+250{,}000`$, and a depth thickness of
$`187{,}500\,\mathrm{m}`$ is compatible because
$`187{,}500=125{,}000+62{,}500`$. If
$`\mathtt{region\_depth\_shallow\_m}=10{,}000\,\mathrm{m}`$, the corresponding
bottom coordinate is
$`\mathtt{region\_depth\_deep\_m}=197{,}500\,\mathrm{m}`$.

For topographic meshes, choose `maximum_octant_level` from the wavelength
target. Let $`V_{S,\min}^{\mathrm{model}}`$ be the smallest shear velocity in
the active solid material model. Choose the smallest integer $`N`$ satisfying

```math
\frac{
\max\!\left(
\mathtt{region\_length\_north\_m},
\mathtt{region\_length\_east\_m},
\mathtt{region\_depth\_deep\_m}-\mathtt{region\_depth\_shallow\_m}
\right)
}{2^N}
\le
\frac{
\max\!\left(V_{S,\min}^{\mathrm{model}},
\mathtt{simulation\_shear\_velocity\_min}\right)
}{
\mathtt{simulation\_wave\_max\_freq\_hz}\,
\mathtt{simulation\_node\_per\_wavelength}
},
```

then set $`\mathtt{maximum\_octant\_level}=N`$. The accepted input range is
0 through 30. For example, with a maximum extent of
$`1{,}000{,}000\,\mathrm{m}`$, $`V_{S,\min}^{\mathrm{model}}=500\,\mathrm{m/s}`$,
$`\mathtt{simulation\_shear\_velocity\_min}=400\,\mathrm{m/s}`$,
$`\mathtt{simulation\_wave\_max\_freq\_hz}=2\,\mathrm{Hz}`$, and
$`\mathtt{simulation\_node\_per\_wavelength}=8`$, the wavelength target is
$`31.25\,\mathrm{m}`$. Level 14 gives $`61.035\,\mathrm{m}`$, while level 15
gives $`30.518\,\mathrm{m}`$, so the smallest acceptable value is
$`\mathtt{maximum\_octant\_level}=15`$.
This explicit parameter is required only with `include_topography=yes`.
Acoustic or other feature-specific refinement rules can be more restrictive;
the chosen topography cap must not prevent the finest active rule.

## Host and CUDA compilation

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
families, Manzari--Dafalias, PressureDependMultiYield03, production RIVA-Sand,
and `RIVASand02`, plus small runtime and material-point
objects. The research successor also has its own material-point CUDA object
because its frozen restart state contains 138 values. The runtime selects one
family kernel for the globally configured `material_model`. This keeps the
equations normally inlined, avoids relocatable device code, and lets
`make -j 8` compile the independent kernels concurrently. The split also creates stable
boundaries for moving model-specific implementation details without returning
to a single production PTX-assembly job.
Prefer an explicit outer job count such as `-j 8` over bare `-j`, which gives
GNU Make no concurrency limit.

## Runtime performance controls

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
- `HERCULES_PLANE_IO_MODE=mpiio|hybrid|legacy`
	Selects the plane-output backend. Default is `mpiio`, including when the
	variable is unset or empty: fused GPU interpolation and a packed device-to-host
	transfer feed the distributed MPI-IO writer, in which each rank writes its
	owned strips through a precomputed file view. `hybrid` uses the same GPU
	interpolation followed by the rank-0 strip gather and buffered `fwrite`
	writer. `legacy` selects CPU interpolation and the rank-0 writer. A configured
	dedicated I/O pool (`IO_PES>0`) selects the existing legacy I/O-pool backend
	and logs the fallback if a GPU backend was requested or selected by default.
	All backends preserve the `planedisplacements.N` binary layout. Plane nodes
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
	no effect unless the selected backend is MPI-IO (the default with `IO_PES=0`).
- `HERCULES_PLANE_METADATA=0|1`
	Optional setup-only point audit; default is `0` (disabled, also when unset
	or empty). With `1`, each solver rank owning plane points writes
	`plane_point_metadata.RANK.csv` in the plane output directory after
	`solver_init()` has established CPML membership. This works with all plane
	backends and adds no per-step transfers or work. See
	[plane amplitude diagnostics](#diagnose-large-plane-displacements) for the
	columns and how to join them to existing displacement files.
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
	With a GPU plane backend (`mpiio`, the default, or `hybrid`),
	plane interpolation reads device-resident displacements directly and does
	not enter the every-step host output patch. The legacy plane backend and station fallback paths patch
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

## Advanced and diagnostic controls

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

### OpenMP build option

- `USE_OPENMP` now defaults to `1` in `quake/forward/Makefile`.
	Use `USE_OPENMP=0` only if you want to disable OpenMP explicitly.
	Examples: `make` (OpenMP on by default) or `make USE_OPENMP=0`.

## Phase launcher defaults

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

## Launch a simulation

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
# Plane output defaults to mpiio; choose an alternative only if needed:
# export HERCULES_PLANE_IO_MODE=hybrid
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

## Nonlinear GPU execution and constitutive models

- `HERCULES_NONLINEAR_BACKEND`
	Selects the nonlinear constitutive-update and internal-force backend.
	`gpu` is the default and uses the CUDA implementation for all current
	nonlinear material models, with double-precision constitutive state and
	device-resident force assembly. `cpu` explicitly selects the established
	CPU/GSL reference implementation. Setting `gpu` explicitly remains useful
	for reproducible launch environments. This setting changes execution, not
	the calibrated values in `parameters.in`. The selected nonlinear backend is
	reported to standard error during solver startup.
- `HERCULES_NONLINEAR_STATS_INTERVAL`
	Enables sampled Manzari-Dafalias GPU workload diagnostics when set to a
	positive integer `N`; unset or `0` is the default and adds no diagnostic work.
	Every `N` global time steps, Hercules reports globally reduced integration-
	point-slot counts classified as `zero_or_inactive`, `one_step`, and
	`multi_step`, together with total, mean, and maximum accepted constitutive
	substeps. For Forward Euler, `one_step` identifies an elastic-predictor return
	or a one-part plastic update. With the default fixed policy, a plastic point is
	`multi_step` whenever `no_substeps > 1`; with the opt-in demand policy, the
	classification reflects the requested count. For Modified Euler, a
	one-step point can be either
	elastic or a plastic update accepted in one substep. Zero/inactive slots also
	include elements skipped because their displacement increment is zero and the
	three unused slots of a five-point nonlinear-topography element. The sampled
	path clears only the existing diagnostic substep-count array, runs an extra
	device scan, and performs global MPI reductions; use it for diagnosis rather
	than production timing. All solver ranks must use the same interval.
- The specialized Manzari--Dafalias production kernel keeps only the committed
	constitutive state it consumes. It does not allocate or write the four legacy
	previous-state shadow arrays required by the classic yield-surface/backbone
	kernel. This removes 1,216 bytes of unused device state per nonlinear element
	and the same amount of device-memory write traffic per updated element and
	time step. Startup reports the omitted allocation. Its failure-only
	constitutive diagnostics use a per-thread shared-memory record so detailed
	failure reports do not enlarge the normal production thread's local stack.
	PDMY03 instead retains the legacy shadow allocation/writes and its compact
	thread-local diagnostic record. Those shadows do not enter the PDMY03
	constitutive equations, but controlled target-system production measurements
	favored the established generated kernel over their omission despite better
	static PTXAS spill totals. This model-specific execution choice does not alter
	the committed stress/history state.
- `material_model=ManzariDafalias` is a device-first, pressure-dependent
	effective-stress sand model. Production execution currently requires the GPU
	backend. `manzari_dafalias_integration_scheme=modified_euler` is the default
	and uses the established adaptive explicit Modified-Euler implementation.
	`manzari_dafalias_integration_scheme=forward_euler` selects the OpenSees
	integration-scheme-5 constitutive stage. After any elastic/yield intersection,
	Hercules applies that first-order stage to each fraction. The optional
	`manzari_dafalias_forward_euler_substeps` control defaults to `fixed`, for
	which `no_substeps` is the exact count over every plastic remainder. This
	preserves established inputs and integration-sensitivity studies. The opt-in
	`demand` policy treats `no_substeps` as a cap and sizes work with an MD-specific
	stress-ratio/deviatoric-strain heuristic. It is not the PDMY03 subdivision
	criterion, which also depends on nested-surface spacing, confinement,
	elastoplastic modulus, and volumetric strain. On one dense-sand cyclic path at
	30 m depth, changing from a fixed 322 substeps per step to a demanded mean of
	1.0 moved shear-stress NRMS against the Modified-Euler reference from 0.20% to
	0.44%; this is a path-specific characterization, not general validation.
	Setting `no_substeps=1` recovers OpenSees's original single-update behavior
	under either policy. Forward Euler and its subdivision policies are compiled
	as separate device entry functions, so ptxas allocates registers for each
	entry. While the scheme was a runtime branch, a Forward-Euler run inherited
	Modified Euler's allocation: 255 registers with 6,564 bytes of spill stores
	and 21,532 bytes of spill loads per thread, against zero spills for PDMY03 on
	the same build.
	Specializing the integrator lets each entry carry only its own working set. On
	the measured demand-policy build, Forward Euler's spill loads fell from 21,532
	to 2,108 bytes per thread and its spill stores from 6,564 to 508, while
	Modified Euler was unchanged.
	In that measured demand-policy build, both the demand Forward-Euler and
	Modified-Euler entries remained at the 255-register ceiling because ptxas
	allocates up to the cap and spills only the excess. Occupancy therefore
	remained 2 blocks of 128 threads (8 of 64 warps) per SM.
	`NL_MD_FE_MIN_BLOCKS_PER_SM` (default 2) and `NL_MD_ME_MIN_BLOCKS_PER_SM`
	(default 2) set that bound per scheme and are compile-time overridable, e.g.
	`-DNL_MD_FE_MIN_BLOCKS_PER_SM=3`. The useful value is empirical: raising it
	lowers the per-thread register budget to `65536/(block*N)` and trades more
	spilling for more latency hiding. The neutral value 2 preserves the measured
	allocation until target-system timing justifies a higher residency.
	Forward Euler has no embedded local-error estimate, adaptive retry, or stress
	correction. It therefore ignores `backbone_errTol` and
	`manzari_dafalias_recovery`, reports the resulting absolute yield drift, and
	is more sensitive to substep size. Hercules retains its finite-state and
	minimum-effective-pressure safeguards. Modified Euler remains the production
	default.
	For Modified Euler, `backbone_errTol` is the dimensionless local stress-increment error
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
- RIVA-Sand means **R**eversible--**I**rreversible **V**olumetric
  **A**ccumulation model for Sand. `material_model=RIVASand` selects the native binary64 port of the frozen
	Ottawa F65 RIVA-Sand pressure-dependent bounding-surface model. Production
	execution is CUDA-only. The compiled kernel retains 112 calibrated constants,
	four host-increment reversal controls, and all 39 logical history fields; the
	startup report prints parameter SHA-256
	`9585f0c155c9444885c5115e7753a1a5d97c783e2c685938a49a65435c9e8f83`.
	`no_substeps` is the exact number of equal Forward-Euler subdivisions per
	nonzero Hercules increment; an exact zero increment preserves history
	without subdivision. Reversal detection executes once before subdivision, so
	changing this number cannot create extra physical reversal events. Start
	with one and increase it only after a timestep/substep convergence study.
	`backbone_errTol` and `error_tolerance` remain required by the common
	nonlinear input block but do not control this fixed explicit update.
- RIVA-Sand uses one compressive geostatic-admission path for both analytical
	`approximate_geostatic_state=yes` stresses and the committed stress from a
	ramped `approximate_geostatic_state=no` solve. An admissible stress is
	unchanged. An over-bound stress is preserved, its mapping origins are
	recentered, and a closed-form radial rule prevents its normalized deviatoric
	stress from moving farther from the bounding cone. This adds no iterations
	or constitutive substeps. A physical mean effective stress at or below
	`p_min` is rejected instead of being silently changed at the start of
	shaking. At the first post-geostatic step Hercules initializes every active
	RIVA-Sand integration point and reports the global initialized count, admitted
	count, and maximum admitted cone ratio. The frozen 93-value numerical oracle
	remains unchanged; the committed admission selector is compared separately.
	The exact extension is specified in
	`benchmark_inputs/nonlinear_material_point/RIVA_SAND_GEOSTATIC_ADMISSION.md`.
	The internal residual cone-pressure research control remains zero and is not
	exposed in `parameters.in`.
- RIVASand and RIVASand02 use compact `Vs_limit Dr` rows by default.
  Optional positive G0 after Dr selects an explicit stiffness override:
  `Vs_limit Dr G0`. Their remaining material inputs use reference values.
  RIVASandCustom and RIVASandCustom02 require the full row:

  `Vs_limit Dr G0 M kd h m zeta e_max e_min Q R n_G`.

  Use G0=0 in every row for automatic stiffness, or positive G0 in every
  row for explicit stiffness; mixed modes and inconsistent row widths are
  rejected. G0 is dimensionless. The reference value 483.48301127222084
  preserves the former 48.976829 MPa stiffness at 101.3 kPa. Every numeric
  input must be finite; Dr and n_G lie in [0,1], G0>=0, M,kd,h,Q>0,
  m,zeta,R>=0, and e_max>e_min>=0 containing the anchor e=0.601.

  Automatic stiffness targets the actual mesh bulk-density times Vs squared,
  not Vs_limit and not buoyant density. After stress-preserving geostatic
  admission, each integration point is calibrated using its full initialized
  shear-modulus evaluator, including pressure guards and bias multipliers.
  Approximate=yes uses its admitted K0 field; ramped no uses the actually
  solved effective stress. The cached elastic scale is fixed throughout
  dynamics, so pressure-dependent degradation is retained. Bulk stiffness
  follows the existing nu=0.30 and does not also match mesh Vp. Near-zero
  effective stress must still pass admission; it is not silently changed.

  This adds one target-modulus double per element and one elastic-scale double
  per quadrature point; unrelated families allocate neither. Per-step updates
  load the cached scale but do not repeat the calibration. Checkpoints now
  store the per-point scale and require marker 4 (original), 8 (02 type 1),
  or 7 (02 type 2). Both `RIVASand02` and `RIVASandCustom02` accept the
  optional scalar `rivasand02_reversaltype = 1` in `parameters.in`; omission also means
  type 1. Type 1 compares adjacent valid increments. Set `rivasand02_reversaltype = 2`
  to retain the previous accumulated loading-branch rule. Detection occurs
  once per host increment before substeps. Only the integers 1 and 2 are
  accepted, and only for the two 02 names; it is not a material-row column.
  The material-point test input accepts the same key. Existing marker-7
  checkpoints require explicit `rivasand02_reversaltype = 2`; cross-type restarts and
  marker-6 or earlier 02 files are rejected. Retain the same model, mesh,
  physical inputs and reversal type on restart. Original RIVA-Sand inputs
  and checkpoints are unchanged.

  The material-point driver takes Dr [G0] or the twelve Custom values,
  without Vs_limit. Automatic mode additionally requires
  material_point_riva_density (kg/m3) and material_point_riva_vs (m/s).
  Full reference histories are retained using explicit reference G0.

  The 02 kernel remains the distinct revision-3 plastic-activity successor,
  with unchanged 138-value history and material-specific density correction.
  OpenSees BiasVolume and reversalLatch are not included. Both generations
  return effective skeleton stress; BrickUP owns pressure and its total-stress
  contribution. Automatic initial Vs matching is not a new cyclic calibration
  and does not establish full-mesh CUDA/MPI qualification.
  See `benchmark_inputs/nonlinear_material_point/RIVA_VELOCITY_STIFFNESS.md`.
- RIVA-Sand is compatible with BrickUP through explicit field ownership.
	RIVA-Sand is the effective-skeleton law: its compatibility state separates
	total volumetric strain into elastic confining, irreversible plastic, and
	reversible plastic parts, and reconstructs effective mean stress from the
	elastic confining part. BrickUP independently solves fluid mass balance,
	`d(zeta) = alpha*d(eps_v) + d(p_f)/M`, together with Darcy flow and the
	configured drainage boundaries. Equilibrium then uses
	`stress_total = stress_effective - alpha*p_f*I`. Hercules never copies the
	standalone diagnostic `1-p'/p'_anchor` or RIVA-Sand effective-pressure loss into
	`p_f`, so the same pore-pressure contribution is not generated twice. The
	nonlinear undrained fast path remains disabled, as for MD and PDMY03.
- Checkpoints preserve the historical displacement prefix and append a
	versioned coupled-state section containing both BrickUP pressure time levels
	and the complete RIVA-Sand integration-point state, including prior strain/stress.
	A RIVA-Sand/BrickUP restart rejects a legacy displacement-only checkpoint instead
	of silently resetting pore pressure or constitutive history. It also requires
	`RIVA_KERNEL_REVISION=2`, so a checkpoint created before geostatic admission
	and non-expansive regularization cannot silently resume under new equations.
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
	For the coupled SSPbrickUP/BrickUP triaxial and cyclic-shear workflows,
	`plot_opensees_pdmy03_brickup.py` writes stress-response, pore-pressure,
	effective-stress-path, and total-confinement-drift panels.
- The public Physics models guide now contains a complete `parameters.in`
	nonlinear block for every accepted `material_model`: `linear`,
	`vonMises_ep`, `vonMises_fa`, `vonMises_faM`, `vonMises_baE`,
	`vonMises_baH`, `vonMises_GQH`, `vonMises_MKZ`, `vonMises_RO`,
	`MohrCoulomb`, `DruckerPrager`, `ManzariDafalias`,
	`PressureDependMultiYield03`, `RIVASand`, `RIVASandCustom`,
	`RIVASand02`, and `RIVASandCustom02`. The displayed rows are syntax or regression
	fixtures, never implicit production calibrations. They use two identical
	property rows to make the nonlinear Vs interval explicit without adding an
	interpolation gradient.
- The documented pressure-dependent recipes for Mohr--Coulomb,
	Drucker--Prager, Manzari--Dafalias, PDMY03, and RIVA-Sand set
	`approximate_geostatic_state=yes`,
	`geostatic_loading_time_sec=0`, and `geostatic_cushion_time_sec=0`. This is
	the internally consistent approximate-K0 choice when a case has neither a
	calibrated ramp nor a valid restart. The parser rejects a positive ramp or
	cushion together with the approximate state. If a calibrated ramp is used,
	set the approximate flag to `no`; for a sealed, zero-permeability nonlinear
	BrickUP case, use the approximate effective-stress state, a permeable
	consolidation stage, or a restart rather than a positive undrained ramp.
- The same guide includes reviewable configuration blocks for acoustic--elastic,
	poroelastic--elastic, nonlinear BrickUP, and acoustic--poroelastic open-pore
	interfaces; lossless, Rayleigh, mass, material-Q BKT, and parametric-Q BKT
	damping; and baseline, attenuating, heterogeneous, acoustic, and topographic
	CPML operation. Keep these snippets synchronized with
	`docs-site/app/docs/content.ts`, the parameter reference, and the parser.
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
	particularly important for Manzari--Dafalias, PDMY03, and RIVA-Sand, whose skeleton stress is
	effective stress, but it also avoids bypassing every other nonlinear
	constitutive volumetric update.
- The synthetic fixture includes a standalone OpenSees driver for both
	Modified Euler and Forward Euler, plus a characterization comparator under
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
	force/displacement norms before a production run. Manzari--Dafalias,
	PDMY03, RIVA-Sand, and the intermediate-bias research successor are
	production-GPU-only; validate each with its shared host/CUDA
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

## Mesh generation and import

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

### Separate mesh and solver settings

- Phase A (generate mesh database):
	`output_mesh=1`, `mesh_only_mode=1`, `mesh_input_mode=0`
- Phase B (import mesh and solve):
	`output_mesh=0`, `mesh_only_mode=0`, `mesh_input_mode=1`

## Acoustic and elastic coupling

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
	Optional. `bathymetry` (default), `all`, or `bottom_only`. This applies when
	acoustic interface coupling is enabled. Explicit `all` couples every
	acoustic-elastic interface face using mesh-face normals.
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
	Omitting the parameter now selects `bathymetry`; an existing explicit `all`
	line still selects all faces. The `default` alias also now selects
	`bathymetry`. Other accepted aliases are `reconstructed_bathymetry` and
	`bathymetry_normal` for `bathymetry`, `none` for `all`, and `bottom` or
	`z_plus` for `bottom_only`. Values are case-insensitive.
	`bathymetry` omits lateral and upward-facing acoustic-solid interfaces, so
	check that this representation is appropriate for the model geometry.
	Stability with this filter does not by itself validate waveform accuracy.
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
	external domain face use the condition selected by
	`acoustic_free_surface_gravity`. When `acoustic_interface_coupling=yes`,
	scalar acoustic formulations build an explicit acoustic-solid interface table. In octree
	meshes, nonconforming interface faces are split into area-weighted fragments
	with fragment centers and normals, and CUDA kernels apply the normal
	traction/acceleration coupling only on those fragments. Each fragment
	interpolates the scalar acoustic field on the acoustic element and
	displacement acceleration on the adjacent elastic element, so refined elastic
	neighbors can couple to coarser acoustic faces. Fragment subdivision follows
	`acoustic_interface_max_level` and finer neighboring leaves, including
	remote leaves resolved on their MPI owner. The table uses the owner's
	corrected mesh-level acoustic/elastic classification, so mixed elements kept elastic still
	couple to adjacent acoustic elements. Subdivision has a per-face cap of
	`32 x 32` fragments to avoid accidental table blow-up.
	The legacy `acoustic_pressure_formulation` key is still accepted as a
	backward-compatibility alias only when `acoustic_formulation` is absent:
	`yes` maps to `velocity_potential`, and explicit `no` maps to `none`. New
	input files should use `acoustic_formulation`.
- `acoustic_free_surface_gravity`
	Optional yes/no control for the flat acoustic free surface. Default is `no`,
	which preserves the established pressure-release condition (`p=0`, or
	`phi[n+1]=phi[n]` for velocity potential). `yes` enables the linearized
	gravity condition `q_tt = g q_z`, with `g=9.81 m/s^2`, for either scalar
	formulation. Hercules applies it weakly by adding
	`integral_Gamma N_i N_j/(rho g)` to the acoustic surface mass; the Hex8 top
	face is lumped as `h^2/(4 rho g)` per node. This supports linear surface-
	gravity waves without adding a bulk `rho g` load, a moving mesh, or a
	hydrostatic prestress solve. The flat acoustic surface is at solver
	`Z=topographybase_zcoord` with topography, `Z=surface_shift_m` with buildings,
	or the octree top otherwise. A Hex8 top face must coincide with that datum;
	with topography its center must also be exposed at the topographic surface,
	so buried water beneath land is excluded. This supports flat water surfaces,
	not a general curved or cut-cell fluid free-surface operator. The same
	detection selects the pressure-release nodes when this option is `no`.
	The startup `Acoustic free-surface geometry` line reports the solver datum,
	exact exposed face count, and summed harbored surface-node copies (the last
	count includes MPI replicas). Where this face intersects a
	lateral acoustic CPML, the surface mass uses the tangential Jacobian
	`kappa_x*kappa_y`; the solver also adds the paired `x/y` mass-proportional
	damping and CFS restoring/memory terms on the four top-face nodes. The
	vertical CPML stretch is intentionally excluded because `z` is normal to the
	gravity surface. Corner PML faces include both tangential directions. This
	adds no per-step work to non-surface PML elements. The current face-average
	memory has a discretization inconsistency described under
	[scalar acoustic limitations](#scalar-acoustic-corrections-and-limitations). Changing this
	option affects phase-B solver setup but not mesh generation, so an existing
	compatible mesh etree can be reused.
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
- `print_station_pore_pressure_ratios`
	Optional; default `no`. With `include_poroelastic=yes` and
	`include_nonlinear_analysis=yes`, `yes` appends `Ru(-)` to every normal
	`station.N` row using the conventional experimental definition
	`r_u(t) = Delta_u(t) / sigma'_v0`. `Delta_u` is the BrickUP excess-pressure
	unknown (the same value written by the optional `Pressure(Pa)` column), and
	`sigma'_v0 > 0` is the initial vertical effective-stress magnitude at that
	station. With `approximate_geostatic_state=yes`, the denominator comes from
	Hercules' integrated 1-D total-overburden column, adjusted from the owning
	element center to the exact station depth, minus
	`poroelastic_biot_alpha * rho_f * g * z_w`. With a ramped geostatic phase,
	`Ru(-)` remains `nan` while that numerical reference does not yet exist.
	The last geostatic force evaluation fixes the denominator from the committed
	vertical skeleton stress at the owning nonlinear element's station-output
	quadrature point, before the first dynamic output row. This is the same
	constitutive state used by the nonlinear stress columns and to initialize
	RIVA-Sand before dynamics. The final value is
	written as `sigma_v0_eff_pa` in `station_mediums.txt` for auditability.
	Hercules rejects `print_station_pore_pressure_ratios=yes` when
	`approximate_geostatic_state=no` and both geostatic times are zero, because
	that combination supplies no physical initial effective-stress denominator.
	The valid alternatives are approximate initialization with zero loading and
	cushion times, or numerical geostatic initialization with a positive loading
	or cushion interval.
	The ratio switch also activates the sparse/fused pore-pressure station
	gather when `print_station_pressures=no`, so enabling only `Ru(-)` does not
	download the bulk pressure field. Hercules does not clip the ratio: dilation
	can give negative values and strong buildup can exceed one. Stations outside
	an active saturated BrickUP element, CPML stations, and zero/nonpositive
	reference-stress locations write `nan` rather than a misleading ratio.
- `reconstruct_acoustic_plane_displacements`
	Optional. `yes` reconstructs displacement at acoustic output-plane points
	from the pressure gradient and writes those reconstructed values into the
	standard `planedisplacements.N` files. Elastic plane points continue to use
	the elastic displacement field, so mixed water/solid planes remain a single
	three-component displacement file. This requires `include_acoustic=yes` and
	`acoustic_formulation=pressure` or `acoustic_formulation=velocity_potential`.
	Default is `no`, which preserves the legacy plane-output behavior.
	With `velocity_potential`, the stored scalar is first converted to physical
	pressure, `p[n] = -(phi[n] - phi[n-1])/dt`, before evaluating `-grad(p)/rho`.
	This matches station reconstruction in both formulations. The existing
	per-step host buffer swap retains the previous plane-node scalar values,
	so the resident GPU path needs no additional gather or device-to-host copy.
	Inside CPML, both station and plane reconstruction apply the inverse
	coordinate stretch to the pressure gradient before time integration. The
	output point keeps its own three gradient memories, using coefficients from
	its containing element's solver CPML profile; this adds no field exchange.
	The startup log reports the selected reconstruction formula. See the
	[acoustic correction note](#scalar-acoustic-corrections-and-limitations)
	when interpreting plane files from older executables.
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
	Optional ASCII station diagnostic output path; ignored with
	`station_binary=yes`. The solver writes this metadata
	file during station setup and rewrites it after a ramped geostatic reference
	refresh. It records each station's owning medium
	(`acoustic` or `elastic`), whether pressure is defined there, owner rank,
	element id/level, local coordinates, and material properties. If unset, it
	defaults to `output_stations_directory/station_mediums.txt`. The metadata
	includes `sigma_v0_eff_pa`, the positive station-specific denominator used
	by optional `Ru(-)` output; it is `nan` where the initial effective-stress
	reference itself is unavailable or nonpositive.
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
	attenuating Earth do not need it. The older flat-water benchmark reported
	`acoustic_quality_factor = 100` as a useful tradeoff; the
	[historical record](#historical-acoustic-experiments) gives its scope.
	That case-specific observation is not a general production calibration.
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

## Basin-bottom grid

Enable the basin in `parameters.in` and locate its grid origin geographically.
The coordinates below are illustrative; replace them with your grid's origin.
Keep a valid background `cvmdb_input_file` or profile for points outside the
basin. These settings explicitly select the basin replacement:

```ini
Istanbul_velocity_model    = no
3D_velocity_model          = no
basin_velocity_model       = yes
basin_input_file           = input/basin.in
basin_origin_latitude_deg  = 34.0000
basin_origin_longitude_deg = -118.5000
basin_material_type        = elastic
```

Relative file paths are resolved from the solver's working directory. The file
starts with `Nx Ny 1`, followed by exactly `Nx * Ny` rows of
`x y bottom_depth`. Use whitespace-separated numbers without column labels or
comments. The current lookup uses one bottom surface, so set the third header
value to `1` and use at least two points on each horizontal axis.

- All three columns are in metres. Solver X is north and Y is east, measured
  from the basin origin. Bottom depth is positive downward and uses the same
  reference as solver Z, rather than thickness below local topography.
- Start at `(0, 0)` with a positive uniform spacing and `dx = dy`. `readBasin`
  infers the common spacing from the first two X values.
- Write X fastest: loop over Y outside and X inside. For zero-based indices
  `i` and `j`, data row `j * Nx + i` contains
  `i * spacing`, `j * spacing`, and the bottom depth at that location.

For example, `input/basin.in` can contain this illustrative 3-by-3 grid with
1,000 m spacing and a 500 m bottom depth at its centre. This is a file-format
example, not a validated simulation case:

```text
3 3 1
0    0    0
1000 0    0
2000 0    0
0    1000 0
1000 1000 500
2000 1000 0
0    2000 0
1000 2000 0
2000 2000 0
```

**Current maximum-edge limitation:** `interpolateZ` accesses the next X/Y
points without clamping at the last row or column. Extend the grid beyond
the largest solver query in each direction so queried basin-relative
`X < (Nx - 1) * spacing` and `Y < (Ny - 1) * spacing`. Queries exactly on
either maximum edge can access invalid grid data.

`isInBasin` bilinearly interpolates the bottom depth and classifies a point
inside the horizontal footprint as basin material only when
`z < bottom_depth(x, y)`. Points on or below the bottom, or outside the
horizontal grid, use the background `cvmdb_input_file` or profile. A zero
bottom depth excludes points at nonnegative solver depths.

The grid supplies geometry. For `basin_material_type=elastic`,
`quake/forward/basin.cc` supplies the hard-coded depth-dependent material
formulas. For `basin_material_type=acoustic_water`, `include_acoustic=yes` is
required; `basin_water_vp` and `basin_water_rho` set the water properties and
`Vs`, `Qs`, and `Qp` are zero.

The existing [bowl basin file](../benchmark_inputs/bowl/bowl_basin_250m.in)
has header `162 162 1`, 250 m spacing, and the same X-fastest ordering.
The current reader and lookup are `readBasin`, `interpolateZ`, and `isInBasin`
in `quake/forward/psolve.cc`.

## 3D velocity-model files

- With `3D_velocity_model=yes`, `3D_velocity_model_directory` must contain
	`x.in`, `y.in`, `z.in`, `vs.in`, `vp.in`, `rho.in`, and `index.in`.
	Every file begins with an integer count followed by that many values.
	`x.in` and `y.in` contain at least two strictly increasing, uniformly
	spaced coordinates in metres. Solver X is north and solver Y is east; the
	axis spacings may differ from one another.
- Let `nx` and `ny` be the counts in `x.in` and `y.in`. Horizontal profiles
	are stored with Y as the outer loop and X as the inner loop. For zero-based
	axis indices `i` and `j`, the profile number is
	`p = j * nx + i`, giving `(x0,y0), (x1,y0), ...` before advancing to
	`y1`. The count in `index.in` must be `nx * ny`.
- `index.in` stores the zero-based start of each vertical profile in the
	concatenated arrays; it does not contain a final end offset. Profile `p`
	uses `[index[p], index[p+1])`, while the final profile uses
	`[index[last], number_of_z_samples)`. The first start is zero and all
	subsequent starts must increase.
- `z.in`, `vs.in`, `vp.in`, and `rho.in` have identical counts and aligned
	entries. Within each profile, Z is depth in metres, positive downward,
	strictly increasing, and unique. Profiles may have different depth values
	and layer counts. `Vs` and `Vp` are m/s; density is tonne/m3 and is
	multiplied by 1000 when Hercules constructs kg/m3 material density.
- The lookup selects the deepest listed Z not exceeding the query depth at
	each of the four surrounding horizontal columns, so vertical properties are
	piecewise constant and each Z is a layer-top depth. It then bilinearly
	interpolates the selected material values horizontally. The first sample is
	used above the shallowest listed depth and the last sample below the deepest
	listed depth. The public input-files guide contains a complete 2-by-2
	row-table-to-file conversion with unequal profile lengths.
- `script/convert_3d_velocity_model.py` converts a row-wise CSV with
	`x,y,z,vs,vp,rho` columns into the required count-prefixed files. Optional
	`qs,qp` columns must be provided together. The converter sorts profiles into
	the required Y-outer/X-inner order and validates the rectangular grid,
	horizontal spacing, vertical depth uniqueness, and material ranges:

	```bash
	python3 script/convert_3d_velocity_model.py model.csv \
	    --output-dir input/materialfiles
	```

### Quality-factor files

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

## Poroelastic U-P formulation

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

## Scalar acoustic corrections and limitations

- **Gravity-surface CPML consistency limitation.** Surface mass and damping
	are lumped at the four face nodes, while the current restoring/memory
	kernel uses a face average. A zero-mean face mode therefore receives nodal
	damping but no matching restoring/memory force. A nodal-memory candidate
	passed analytical CFS step-response checks and an 8 s shifted-surface
	GPU/MPI smoke test, but its reduced Izmir run produced nonfinite pressure
	at 4.833 s. It was rejected and is not included in the production solver.
	The failed experiment does not establish the cause of the original slow
	motion. See the [Izmir investigation](izmir_gravity_diagnostic_20261001.md)
	for scope, raw-output checks, and the separate CPML-thickness failure.
- **Coupled displacement synchronization.** Acoustic interface forcing reads
	new solid displacements only after MPI sharing and hanging-node assignment.
	The existing displacement exchange runs before the scalar update, using
	either CUDA-aware MPI or the host fallback; it is not repeated after the
	buffer swap. Earlier builds sampled unsynchronized displacements and could
	become unstable on adaptive meshes split across ranks even when a one-rank
	run was stable. The startup line `Acoustic interface displacement` reports
	this ordering. Compare one-rank and multi-rank wavefields when validating it;
	finite output and agreement between I/O backends alone do not establish
	stability.
- **Remote interface geometry.** Setup requests final neighbor leaf bounds
	and corrected material properties from their owning MPI ranks. Remote
	octree aggregates can be coarser than those leaves; treating an aggregate
	as the neighbor changed interface fragment quadrature with the partition.
	The startup `Acoustic interface geometry` line gives exact face/fragment
	counts. Batched discovery also supplies the remote elastic CPML geometry;
	its cache is released after solver initialization and adds no per-step
	communication. A recorded local shifted-surface test checked all 1024
	fragment centers and areas and compared one-rank and four-rank wavefields.
- **Shifted free surface and CPML output correction.** Earlier surface
	detection only recognized the octree's `Z=0` face. Embedded water at a
	nonzero topography datum therefore missed both pressure-release flags and
	gravity surface inertia, including the lateral CPML surface terms. Detection
	now follows the physical flat surface described above. This changes the
	acoustic wavefield in affected models; it requires a new simulation, not
	postprocessing of old output. In the Izmir diagnostic with datum 2812.5 m,
	the old log reported zero gravity/CPML surface faces despite water at the
	shifted surface. Check the new exact face count before interpreting results.
	Acoustic station/plane reconstruction also now uses the same gradient-side
	CPML stretch as the solver instead of the raw pressure gradient. For each
	axis, `psi[n]=b*psi[n-1]+a*grad(p[n])` and
	`grad_CPML(p[n])=grad(p[n])/kappa+psi[n]`; outside CPML this is an exact
	identity. The existing host reconstruction applies it at each solver step
	using setup-time coefficients, with no additional GPU transfer. The startup
	`Acoustic output CPML gradient correction` line counts affected observations.
	Recorded local host checks covered a compact potential pulse's integrated
	stretch gain `alpha/(kappa*alpha+d)` and zero residual velocity, including
	zero-alpha and unstretched axes. Separate checks covered translated surface
	selection, mass assembly, pressure-release flags and buried-water exclusion.
	A local CUDA/MPI smoke case checked surface enrollment, finite excited
	fields and station/plane agreement inside CPML. Those investigation tools
	and inputs are not distributed with the repository; their results do not
	validate the full Izmir wavefield.
- **CPML output interpretation.** Values inside the absorbing layer are
	computed in its stretched coordinates and attenuated by its artificial
	operator. Use interior observations for physical motion comparisons;
	correcting the reconstruction does not turn CPML samples into physical
	water or ground-motion observations.
- **Velocity-potential plane reconstruction correction.** Earlier plane
	output treated the stored velocity potential as pressure when
	`reconstruct_acoustic_plane_displacements=yes`. It therefore integrated
	`-grad(phi)/rho` as acceleration instead of differentiating the potential
	in time first. This introduces an incorrect sign/time integration and can
	produce artificial drift in acoustic plane output with any I/O backend.
	The corrected plane path uses the same physical-pressure conversion as
	stations. Pressure-formulation plane output and solid displacement
	interpolation retain their existing behavior; the wavefield evolution is
	not changed by this output correction.
	A recorded local host regression used a two-second potential pulse whose X velocity is
	`0.2*sin(pi*t/2)^2 m/s` during the pulse and zero afterward. Its final
	displacement is 0.2 m. The old plane code instead gives -0.2 m at 2 s and
	-0.6 m at 4 s, with a spurious -0.2 m/s residual velocity. The corrected
	code gives 0.2 m with zero residual velocity, agrees with the production
	station routine in all three components, and preserves the analytic
	constant-acceleration pressure case.
	These checks exercise production host reconstruction with synthetic fields;
	they do not validate an entire CUDA/MPI acoustic-gravity/CPML simulation.
	Regenerate affected acoustic plane output with a rebuilt solver. Agreement
	between two old I/O backends cannot detect this shared reconstruction error.
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

## Absorbing boundaries and CPML

Hercules supports two absorbing boundary conditions, selected at run
time. The default is the legacy first-order Lysmer dashpot. The newer
option is an unsplit Convolutional Perfectly Matched Layer (C-PML,
Komatitsch & Martin 2007 — the same second-order displacement form
SPECFEM3D's FE-PML uses). The C-PML is per-element, wraps the five
absorbing faces (`x_min`, `x_max`, `y_min`, `y_max`, `z_max`) and leaves
`z_min` as the free surface, and works for both elastic and
acoustic-pressure elements. Earlier CPS comparisons are summarized below and in the
[historical benchmark record](#historical-cpml-benchmark-record); those
results are specific to the recorded configuration and implementation.

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
	earlier reference setting. Values `> 1` improve grazing/evanescent absorption
	and are now consistently implemented, but still need re-validation;
	treat `2–3` as promising-but-unconfirmed and prefer a thicker PML
	(`cpml_layers` ~30–50) for the deep-corner vertical residual. Values
	above `3` are untested. `kappa = 1` keeps the PML↔interior interface
	exact.
- `cpml_alpha_max_hz`
	Optional. When `> 0`, sets the CFS frequency directly (Hz). When `< 0`
	or unset (default), it is derived from `simulation_wave_max_freq_hz`.
	The explicit value becomes `alpha_max=2*pi*cpml_alpha_max_hz`,
	whereas the fallback is `pi*simulation_wave_max_freq_hz`; see the
	[coefficient formulas](#cpml-coefficient-and-operator-details) before
	matching a source reference frequency.
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
	correction.

	`rayleigh` enables the Rayleigh material-damping CPML path and is valid
	only with `type_of_damping=rayleigh`. If `type_of_damping=rayleigh`,
	`abc_type=cpml`, and `cpml_pml_material_damping` is omitted, the solver
	auto-promotes it to `rayleigh`; explicit `none` is honored for
	debugging/regression with a warning.
- `cpml_clamp_outer`
	Optional. `yes` or `no` (default `yes`, effective with `abc_type=cpml`).
	Omitting this key enables the clamp; explicit `no` disables it. `yes` pins
	displacement to zero on the four outer lateral faces and bottom. It also clamps the
	acoustic scalar: pressure in pressure mode, or density-scaled velocity
	potential in velocity-potential mode. A zero potential does not require
	zero reconstructed water velocity, which depends on its spatial gradient.
	The top face is not selected by its Z coordinate, but its lateral rim
	can belong to a clamped side face. Pressure-release acoustic top nodes
	use their dedicated update before the clamp; gravity-enabled top nodes
	use the ordinary update and are clamped at that rim. `no` removes this
	mask for both media while retaining CPML absorption. Neither setting
	guarantees zero reflection from the outer termination. Compare interior
	waveforms with one setting changed at a time; smaller amplitudes alone
	do not establish which setting is more accurate.
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

### CPML execution and topography

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
	helper `script/taper_topography.py` does exactly this for
	any `topography.in`; the two example generators also take
	`--edge-taper`. (FEM is a strict no-op for flatten-PML; non-CPML /
	non-VT / flat runs are byte-identical. The
	`benchmark_inputs/topography_tilt_boundary` example is a deliberate
	worst-case test, not a model-setup template.)
	This surface taper is independent of `cpml_3d_material_extension`:
	tapering makes the PML top surface flat, while
	`cpml_3d_material_extension` extends the 3D material/water mask
	through the PML volume.

The reusable utilities and extracted Python examples are indexed in
[`script/README.md`](../script/README.md). The relocated taper follows
`point_elevation`: `flat_index = x_index * n_east + y_index`. Earlier
versions interpreted rectangular topography grids with transposed indices.

Using `script/taper_topography.py` (the recommended
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
	cd /path/to/Hercules_GPU
	python3 script/taper_topography.py path/to/topography.in --out topography_tapered.in
	```

	Then point Hercules at the tapered file — Hercules expects the
	filename `topography.in` in `topography_directory`, e.g.:

	```
	mkdir -p model/topo_tapered
	python3 script/taper_topography.py raw/topography.in \
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
	`--taper-zref`. For these square benchmark domains, output is identical
	to running the standalone utility on the un-tapered surface; the flag
	is off by default so the committed example files stay reproducible.

### CPML stability mechanisms

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

Earlier elastic half-space comparisons against CPS reported improved
late-coda absorption for the tested CPML configurations. The
[historical benchmark record](#historical-cpml-benchmark-record) retains the
receiver-by-receiver table and its implementation caveat. Keep receivers
outside the PML when judging physical seismograms, and distinguish the direct
arrival from the coda dominated by boundary reflections.

### CPML example from earlier validation

```
abc_type = cpml
cpml_clamp_outer = yes      # default: clamp the outer lateral faces and bottom
cpml_layers = 30            # raise to 30–50 to close the deep-corner Z residual
# cpml_kappa_max:           leave at 1.0 (validated); 2–3 promising but unconfirmed
# cpml_stiffness_mode:      leave at native (recommended)
# cpml_alpha_max_hz:        leave default unless the source peak freq differs
# cpml_mpml_cross_fraction / cpml_alpha_floor_fraction: leave at 0.1
# cpml_partition_weight:    leave at 1.0 unless CPML rank/load diagnostics justify tuning
```

Notes:

- The earlier reference configuration uses `cpml_kappa_max = 1`; the deep-corner
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


### Incident plane-wave controls

The separate legacy incident-wave forcing path is enabled by
`include_incident_planewaves=yes`. Its analytic Ricker displacement uses `fc`
as the central frequency in Hz and `Ts` as the time shift in seconds. These
are parsed by `quake/forward/drm_planewaves.cc`; the full conditional geometry
and wave inputs are listed in the
[DRM and wave parameter reference](../docs-site/app/docs/parameter-reference.ts).
They are independent of the structured DRM output recorder.

## Single-point body-force source

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

## Distributed compliant-base force source

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

## Additional derivations and historical results

The following sections explain the current operators and preserve
case-specific benchmark records. Current input descriptions above govern new
runs; historical measurements are labeled with their implementation scope.

### CPML coefficient and operator details

The implementation uses a gradient-side ADE memory to modify strain (or the
scalar acoustic gradient), plus divergence-side inertia, mass damping, and
complex-frequency-shift restoring/memory terms. The absorbing faces are
`x_min`, `x_max`, `y_min`, `y_max`, and `z_max`; `z_min` remains the free
surface. Elastic and acoustic PML elements have separate execution maps.

For layer thickness L, reference compressional speed Vp, target reflection R,
polynomial order n, and centroid distance xi into the layer, the coefficient
builder uses:

```text
r_j     = (xi_j / L)^n
D       = -(n + 1) * Vp * ln(R) / (2 * L)
d_j     = D * r_j
alpha_j = max(alpha_max * (1 - xi_j/L), alpha_floor_fraction * alpha_max)
```

Before cross-axis mixing, an inactive axis has zero profile and alpha.
`cpml_mpml_cross_fraction` supplies p for the M-PML coupling:

```text
d_eff_j = d_j + p * (d_k1 + d_k2)
r_mix_j = min(1, r_j + p * (r_k1 + r_k2))
kappa_j = 1 + (cpml_kappa_max - 1) * r_mix_j
b_j     = exp(-(d_eff_j/kappa_j + alpha_j) * dt)
a_j     = d_eff_j * (b_j - 1) / (kappa_j * (d_eff_j + kappa_j*alpha_j))
```

The code handles zero denominators explicitly. Axes receiving damping only
through cross-coupling also receive the alpha floor. These formulas are in
`compute_pml_profile_coeffs` in `quake/forward/psolve.cc`; n and p above denote
different controls. `cpml_alpha_max_hz>0` sets
`alpha_max=2*pi*cpml_alpha_max_hz`; otherwise the implemented fallback is
`alpha_max=pi*simulation_wave_max_freq_hz`. To match `pi*f_dominant` using the
explicit Hz input, use `f_dominant/2`, not `f_dominant`.

The gradient stretch `1/kappa_j` is paired with element inertia scaled by
`kappa_x*kappa_y*kappa_z`. Earlier strain-only kappa experiments lacked that
matched inertia; their results cannot establish the behavior of the current
matched formulation. The restoring term `+rho*d*alpha*u` and memory term
`-rho*d*alpha^2*zeta` use the same element-centroid displacement so they cancel
at DC. The earlier rejected stress-convolution experiment changed even the
direct arrival under reduced integration; it is not the current operator.

With `cpml_stiffness_mode=native`, the full native element stiffness is
retained and a center-evaluated correction based on `sigma_pml-sigma_raw` is
added. The native normalization uses `mu=9*c1` and `lambda=9*c2`. At an
unstretched, undamped point the correction vanishes, retaining the native
higher modes. `reduced` uses the standalone reduced operator. The independent
`cpml_stiffness_mode_acoustic` inherits the resolved elastic setting unless
overridden. These choices share the coefficient and divergence-side pipeline.

Membership and profiles use centroids: the first layer samples roughly
`xi=h/2`, which can leave a bounded interface residual.

### VT and BKT operator history

The older hill test exposed a mismatch between a full-cube GPU internal
force and the fill-fraction-scaled mass of a VT cut cell. CPU `TetraForcesBKT`
used a fill-scaled five-tetrahedron viscoelastic operator, while the old GPU
paths either applied the full-cube BKT force or skipped the cut-cell effective
force. Small cut-cell masses then made the explicit update unstable.

The GPU tetrahedral force path uses the existing BKT shear/bulk convolution
state, both `cube_part` layouts, and the five tetrahedra with their solid fill
fractions. Its material moduli follow `mu=rho*Vs^2` and
`kappa=rho*(Vp^2-4*Vs^2/3)`. Eligible cut cells are omitted from the full-cube
internal-force path to avoid counting both operators. Flat PML cells remain
on the cube-stiffness/CPML path. See `kernelTopoBktForceEffective` in
`quake/forward/psolve_gpu_launch.cu` and the current
`topography_bkt_gpu_stabilize` description above for activation rules.

The historical 25-station hill report recorded growth near 2.7 s and NaNs by
8.8 s before the fix. It reported convergence for VT+BKT2 under both dashpot
and CPML boundaries and reproduction of the old failure with the explicit
stabilization opt-out. Those are records of the earlier Vista runs, not tests
rerun during this documentation consolidation.

### Historical acoustic experiments

The scalar formulations exchange different time derivatives at the interface:
pressure coupling uses solid normal acceleration; velocity-potential coupling
uses normal velocity and recovers pressure as `p=-dphi/dt`. The earlier
flat-water tests reported that the latter removed a secular instability and
reduced excitation of a lossless-cavity pressure ring. The current
`acoustic_formulation` and optional gravity free-surface descriptions above
include subsequent interface changes.

For positive `acoustic_quality_factor=Q`, the implemented acoustic damping
coefficient is `beta=1/(Q*2*pi*f_ref*dt)`, from
`acoustic_stiffness_damping_beta`. In the modal interpretation of the discrete
force term, `zeta(omega)=beta*omega*dt/2` and
`Q(omega)=1/(beta*omega*dt)`. This is stiffness-proportional damping: it is
stronger at higher frequencies. It is independent of the elastic `theBBase`
coefficient, which an earlier implementation incorrectly used, leaving the
option ineffective for nondamped/BKT cases. Acoustic CPML elements skip this
interior regularizer.

The older flat-water, non-attenuating benchmark reported a roughly 4.6 Hz
pressure ring, estimated by `c/(4*h)`. Its Q sweep (500, 100, 50) reported
that Q=100 reduced pressure energy above 2 Hz by about 50–82%, improved
pressure correlation with CPS from about 0.68 to 0.86, and retained vertical
displacement correlation near 0.9965 with mean perturbation at most 0.5%.
Q=50 increased the displacement perturbation with little additional ring
reduction. These are case-specific historical measurements of a numerical
regularizer, not physical ocean attenuation values or a general calibration.

### Historical CPML benchmark record

The earlier Extended README recorded the following elastic half-space,
dip-slip source comparison against CPS. Entries are late-coda RMS misfit over
25–40 s; smaller is better. The reported run used `abc_type=cpml`,
`cpml_clamp_outer=yes`, approximately 50 layers, and `cpml_kappa_max=3` in
an **older strain-only implementation without matched kappa inertia**.
It is not evidence validating kappa=3, the native correction, or an arbitrary
current configuration. The separate example above uses different settings.

| Receiver | Lysmer | Earlier CFS-CPML (kappa=3) | Other-branch CPML |
| --- | ---: | ---: | ---: |
| R east, 12 km | 0.222 | 0.107 | 0.121 |
| R east, 16 km | 0.649 | 0.153 | 0.400 |
| Z up, 12 km | 0.561 | 0.171 | 0.886 |
| Z up, 16 km | 2.195 | 1.110 | 2.735 |

For R east at 16 km the old report also gave direct-arrival (10–22 s)
misfits 0.151 versus 0.154 for Lysmer, and whole-record (0–40 s) misfits
0.170 versus 0.366. It attributed much of the deep-corner vertical improvement
to layer thickness. Receiver labels above retain the original plotting
convention (including “Z up”); Hercules still stores solver Z positive down.

## Single-file station binary output

`station_binary` defaults to `no`, preserving the existing ASCII `station.N`
histories and station metadata. Set it to `yes` to write displacement only for
the existing `output_stations` list:

```ini
station_binary             = yes
station_output_binary_file = output/station_displacements.0
number_output_stations     = 2
output_stations_print_rate = 10

# latitude_deg longitude_deg depth_m
output_stations =
  34.0000 -118.5000   0.0
  34.0100 -118.4900 100.0
```

Coordinates remain latitude, longitude (degrees), and surface-relative depth
(metres), with the existing geographic conversion and building/topography
shifts. Every station must have exactly one owner. Input row order is preserved,
including duplicate rows. Binary mode uses the existing station displacement
interpolation, including cut-cell topography and scalar-acoustic reconstructed
motion; poroelastic samples contain skeleton displacement.

The file contains native-endian float64 values in metres, shaped
`[snapshot, station, XYZ]` in C order. XYZ means solver north, east, down.
There is no header, timestamp, coordinate table, or per-station file. ASCII
`print_station_*` flags and nonlinear diagnostic columns are ignored; binary
mode does not create `station_mediums.txt` or use `output_stations_directory`.
`station_output_binary_file` is required only in binary mode and is otherwise
ignored. Its path must contain 1 to 240 bytes, relative to the run directory
or absolute (`~` is not expanded). Parent directories are created automatically.

Sampling uses positive `output_stations_print_rate`: steps divisible by this
rate are saved before advancing, including step zero in a fresh run. Each
snapshot is `24 * number_output_stations` bytes. Every run truncates the file,
including checkpoint restarts; preserve earlier segments in a different path.
For a restart the first saved time is the first divisible step at or after the
restart step, multiplied by the solver timestep.

For the example above with solver dt=0.001 s:

```sh
python3 script/read_drm.py output/station_displacements.0 --points 2 --dt 0.01
```

The reader accepts either station counts (`--points`) or DRM geometry (`--box`).
To load station displacement with NumPy:

```python
import numpy as np

nstations = 2  # number_output_stations
raw = np.fromfile("output/station_displacements.0", dtype=np.float64)
if not raw.size or raw.size % (nstations * 3):
    raise ValueError("incomplete snapshot or wrong station count")
u = raw.reshape(-1, nstations, 3)
time = np.arange(u.shape[0]) * 10 * 0.001  # print rate * solver dt
# u[:, 0, :] is displacement for the first output_stations row.
```

The implementation reuses station interpolation and the existing sparse
station GPU downloads, followed by one collective MPI-IO write per saved
snapshot. Ranks without stations participate with zero data. Only local
sampled XYZ values are buffered; the station mode does not use
`drm_output_batch_size`. It can run alongside any number of structured DRM
boxes, with distinct output paths and independent sampling rates. The legacy
three-stage DRM replay workflow remains separate.

## Diagnose large plane displacements

Use `script/read_planes.py` to check a completed displacement file's layout,
finite values and vertical peak. Locate and compare motion using the actual
plane geometry, saved sample times and containing-element classification.
The largest displacements need not occur where two output backends differ
most. Keep raw histories when investigating slow motion.

The storage order is `[snapshot, along-strike point, down-dip point, XYZ]` in
C order, with native-endian float64 displacements in metres. Indices are
zero-based. For zero strike and dip, the grid axes are solver X/north and
Y/east; vector components always remain global X/north, Y/east, Z/down.
The saved interval is `output_planes_print_rate * simulation_delta_time_sec`,
not the solver time step. Geometry, timestamps and byte order are not recorded
in the headerless binary.

For actual containing-element classification, enable:

```bash
export HERCULES_PLANE_METADATA=1
```

With `output_planes_directory=outputfiles/planes`, each solver MPI rank owning
plane points writes `outputfiles/planes/plane_point_metadata.RANK.csv`. The
filename suffix is the rank, not the plane number. Each file can contain
points from multiple planes; ranks without points produce no file. The export
runs once after solver setup establishes CPML membership, before time
stepping, with every plane backend. It adds no per-step work and does not
write pressure, velocity or displacement histories.

Each CSV row has these columns:

```text
plane,i,j,rank,element,x_m,y_m,z_m,medium,cpml,vp_m_s,vs_m_s,rho_kg_m3,h_m,valid,topostrip,acoustic_reconstruction,pressure_defined
```

Element IDs are local to the recorded rank. Coordinates are the actual solver
query positions after plane rotation and topography/building shifts. `medium`
is `acoustic`, `solid`, `inactive`, or `unknown`; `solid` alone does not
identify its constitutive law. `cpml` is the actual element membership flag
(0 outside, 1 inside, or -1 unresolved), independent of medium. Material
speeds are m/s, density is kg/m³, and `h_m` is element width in metres.
`valid` records point validity, `topostrip` identifies tetrahedral topography
interpolation, and the last two flags identify configured acoustic
reconstruction and availability of its full hex pressure stencil. Acoustic
topography tetrahedra lack that stencil.

Join metadata to plane values using `(plane,i,j)`. Reject duplicate ownership
and keep missing ownership unknown. Do not mix rank files from different
runs. A geometric lateral boundary band is only a CPML candidate region;
actual membership can differ near the inner edge or bottom. A short run with
the same saved mesh, materials, topography and plane configuration can supply
metadata for existing output. Use fresh output directories because plane
setup truncates displacement files. A full-duration rerun is unnecessary for
this setup-only classification.

Acoustic-plane output integrates `-grad(p)/rho` every solver step. In
velocity-potential mode it first forms
`p[n] = -(phi[n] - phi[n-1])/dt`; older plane code omitted this conversion
(see the [correction note](#scalar-acoustic-corrections-and-limitations)).
For acoustic CPML points, reconstruction applies the containing element's
inverse gradient stretch and output-point memory before integration. Older
builds omitted that correction. Absorbing-layer samples remain artificial
CPML fields, not physical interior water or ground-motion observations.

## Acoustic pressure and motion diagnostics

For selected plane locations, add stations at the surface and nearby submerged
points. `output_stations` takes latitude, longitude and depth below the local
topographic surface. Convert solver X/north and Y/east using the model's
projection and origin, then verify the actual coordinates and acoustic
ownership in `station_mediums.txt`. Retain existing station IDs and record
any inward offset used at an exact domain boundary.

```ini
station_binary = no
print_station_velocities = yes
print_station_accelerations = yes
print_station_pressures = yes
output_stations_print_rate = 1
```

Read ASCII column headers rather than assuming optional fields have fixed
positions. Require complete finite samples, matching time axes and known
material density before comparing pressure and motion. Recheck ownership
after coarsening the mesh: a formerly acoustic observation can become solid.
After changing CPML thickness, recheck membership separately for each case.

Use separate case copies with the same executable, source, mesh controls,
time step, duration, damping, gravity and output settings. Compare the default
`cpml_clamp_outer=yes` with explicit `no`, changing only that control. CPML
remains active in both runs. A separate `acoustic_free_surface_gravity=no`
case changes the physical free-surface condition; it is not a CPML-only test.
Changing CPML thickness can also change material extension in the absorber.

For a flat gravity surface outside CPML, positive-down coordinates give:

```math
p=-\rho g u_z,\qquad g=9.81\;\mathrm{m/s^2}.
```

This is the linearized dynamic free-surface condition with the vertical sign
converted to Hercules coordinates; see the [MIT wave notes, section 1.2](https://web.mit.edu/fluids-modules/waves/www/material/chap-6.pdf).

The residual $`u_z+p/(\rho g)`$ is a diagnostic only for gravity-enabled,
non-CPML surface points. Lumped surface mass, element-gradient recovery and
the half-step pressure difference limit pointwise agreement. Integrated
motion checks also depend on saved sampling and ASCII precision. Station and
plane agreement cannot exclude a shared reconstruction error. Keep raw
signals; smaller amplitudes alone do not validate a boundary setting.

## Check gravity-wave dispersion and reconstruction

A uniform, lossless water column with a rigid bottom separates surface-wave
dispersion from motion recovered using the first element's gradient. This
reference excludes CPML, adaptive meshes, solid coupling and GPU/MPI effects.
For depth $`H`$, horizontal wavenumber $`k=2\pi/\lambda`$, water sound speed
$`c`$ and depth-positive coordinate $`z`$, the compressible reference has
$`\phi\propto\cosh[\mu(H-z)]`$ and:

```math
\omega^2=g\mu\tanh(\mu H),\qquad \mu^2=k^2-\omega^2/c^2.
```

This follows from the acoustic wave equation, the rigid-bottom condition,
and $`\phi_{tt}=g\phi_z`$ at the surface. Its incompressible limit is the
standard $`\omega^2=gk\tanh(kH)`$ relation in the
[MIT surface-gravity-wave notes](https://ocw.mit.edu/courses/12-802-wave-motion-in-the-ocean-and-the-atmosphere-spring-2008/fca533a5edb5f5d3e975ec4c2e63019f_MIT12_802S08_lec03.pdf).

The October 2026 local operator investigation extracted production Hex8
stiffness and lumped volume/surface masses and checked exact tensor
quadrature. With c=1500 m/s, depth 175.78125 m and wavelength 250 m,
h=21.97265625 m gave about −1.253% frequency error but surface-gradient motion
19.65% below the pressure-based value. At h=10.986328125 m these became
−0.316% and 11.58%. These are historical uniform-grid operator results,
not validation of an Izmir waveform or the complete CUDA/MPI solver.

A static reference integrates the potential equation after a unit sinusoidal
permanent bottom displacement. The exact equilibrium surface vertical
displacement is zero, while horizontal particle displacement amplitude is
$`1/\sinh(kH)`$ metres. The first-element gradient can retain a spurious
vertical offset that converges with refinement. Permanent horizontal particle
displacement alone is not proof of instability.

With zero initial fields, let $`w_n=\nabla\phi_n/\rho`$ (or its causal
CPML-filtered value with zero initial memory). The pressure difference and
motion integrator give the discrete identities:

```math
a_n=\frac{w_n-w_{n-1}}{\Delta t},\qquad
v_n=\frac{w_n+w_{n-1}}{2},\qquad
u_n=\Delta t\sum_{j=0}^{n-1}w_j.
```

Local checks of a compact pulse, a 100 s slow potential ramp and a CPML pulse
confirmed these identities for the recorded implementation. They isolate
output time integration from spatial recovery and field evolution; they do
not validate the potential field in a coupled simulation. Investigation
scripts and smoke-test inputs are retained locally, not distributed in this
repository. See the [recorded investigation](izmir_gravity_diagnostic_20261001.md).

Changing `simulation_wave_max_freq_hz` changes meshing and frequency-dependent
damping; it does not filter the source. Compare CPML variants within the same
configuration and keep the time step fixed unless testing time convergence.
Acoustic mesh resolution alone does not ensure gravity-wave resolution:
gravity wavelengths at a given frequency can be much shorter. Recheck actual
material and CPML membership after changing mesh or absorber settings.

## Compare plane output backends

Compare completed legacy, hybrid or MPI-IO output only when plane geometry,
point order, sampling times and restart segments match. Wait for both runs
to finish and flush output. Check file sizes and expected snapshot counts
independently; equally truncated files do not establish complete output.
A snapshot occupies `nx * ny * 3 * 8` bytes. Thus a 361-by-361 plane with
1,000 snapshots occupies 3,127,704,000 bytes; a 71-by-71 plane with the same
count occupies 120,984,000 bytes.

For finite values, an application-specific comparison may use:

```math
|u_{candidate}-u_{reference}| \leq \mathrm{atol}
  + \mathrm{rtol}\,|u_{reference}|.
```

The absolute tolerance is in metres; the relative tolerance is dimensionless
and scales with each reference value. Reject nonfinite values even when the
bytes match. Finite equality and byte identity differ for signed zeros.
Compare complete histories and record maximum errors and their locations,
not only global peaks. Agreement between backends cannot establish physical
accuracy, mesh/time convergence or correctness of a shared reconstruction.

## Structured DRM box output

`quake/forward/io_drm.cc` records any number of DRM boxes, with one flat binary
displacement file and one ASCII coordinate file per box, using a plane-style
input table:

```ini
number_output_drm_boxes       = 2
output_drm_boxes_print_rate   = 10
output_drm_boxes_directory    = outputfiles/drm_boxes_fullDomain

# top_center_x top_center_y exterior_lx exterior_ly exterior_depth h strike_deg
output_drm_boxes =
  1000.0 2000.0 100.0 100.0 50.0 1.0  0.0
  2000.0 3000.0 200.0 100.0 80.0 2.0 30.0
```

Each row has exactly seven finite values. Blank lines and `#` comments are
accepted. `number_output_drm_boxes` defaults to zero; a positive count requires
the table, a positive print rate, and an output directory. Files are
`drmdisplacements.0` with `drmcoordinates.0`, `drmdisplacements.1` with
`drmcoordinates.1`, etc., in row order. Coordinate output is automatic and
requires no additional parameter. Directories are
created automatically. The complete output path is limited to 240 bytes and
resolves from the run working directory (`~` is not expanded).

The first two fields (`drm_offset_x/y`) locate the **top center** in solver
X/north and Y/east metres from the domain origin. Lengths and depth describe
the **exterior**. For spacing h, the interior dimensions are
`(length_x - 2*h, length_y - 2*h, depth - h)`, sharing the same horizontal
center and top. Lengths must exceed 2h, depth must exceed h, and all three
exterior dimensions must be positive integer multiples of positive h.
Strike is degrees clockwise from north: for box-local (x,y),
`X = center_x + x*cos(strike) - y*sin(strike)` and
`Y = center_y + x*sin(strike) + y*cos(strike)`. Displacements remain in solver
XYZ components; they are not rotated into box-local components.

Each box has **ten complete structured face grids**: exterior front, back,
left, right, bottom, then interior front, back, left, right, bottom. There is
no top lid. Front/back mean negative/positive local X; left/right mean
negative/positive local Y. Define exterior interval counts nx=Lx/h, ny=Ly/h,
nz=depth/h. Every face includes the endpoints of each grid axis:

| Face | Grid rows, columns | Exterior shape | Interior shape |
| --- | --- | --- | --- |
| front / back | local Y, Z | (ny+1, nz+1) | (ny-1, nz) |
| left / right | local X, Z | (nx+1, nz+1) | (nx-1, nz) |
| bottom | local X, local Y | (nx+1, ny+1) | (nx-1, ny-1) |

Rows and columns increase in their positive local directions. Columns vary
fastest. **Shared edges/corners are repeated in adjoining faces**, so each face
can be reshaped to a rectangle. The setup log prints face names, shapes, and
zero-based point offsets. All ten faces share the same binary/ASCII file pair.
For the 100 × 100 × 50 m example at h=1 m, the interior is 98 × 98 × 49 m and
there are **60,406 samples / 1,449,744 bytes per snapshot**.

The first example spans exterior X=950..1050, Y=1950..2050, depth=0..50, and
interior X=951..1049, Y=1951..2049, depth=0..49 at zero strike on flat terrain.
A center of (0,0) would extend into negative X/Y; choose centers that keep the
entire rotated footprint inside the domain. Building/topography shifts follow
station/plane surface-relative depth conventions at each sample, so topography
can warp otherwise planar grids. Every point must have exactly one owning
solid element. Scalar-acoustic points are rejected; poroelastic samples record
skeleton displacement. Sampling does not refine the solver mesh.

Each snapshot concatenates the ten faces as native-endian float64 XYZ triples
in metres. The full file is `[time, concatenated_face_sample, XYZ]` with no
header or embedded timestamps. Keep the input box row to identify face shapes.

The matching `drmcoordinates.N` contains three whitespace-separated columns:
`X_m Y_m Z_m`. These are the actual solver sampling coordinates in metres:
X north, Y east, Z positive downward, after strike rotation and the same
building/topography depth shifts used for interpolation. Z is in the solver
coordinate frame, including those shifts. Ignoring lines beginning with `#`,
data row i (zero-based) gives the position of binary sample i in every
snapshot. The exterior/interior face ordering and repeated edge/corner rows
are identical to the binary layout; no points are sorted or deduplicated.
Values use 17 significant digits, sufficient to round-trip float64 coordinates,
with a decimal point independent of the host locale. The file begins with
comment lines describing the axes, ordering, point count, and column names.
For the first example on flat terrain with zero surface shift, its data begin:

```text
# X_m Y_m Z_m
950 1950 0
950 1950 1
950 1950 2
```

Solver rank zero writes and closes the coordinate file once during setup,
after validating point ownership. Each run overwrites it, including restarts,
alongside the binary file. It adds no per-step I/O or GPU transfers.

Inspect the first example (assuming solver dt=0.001 s) with:

```sh
python3 script/read_drm.py outputfiles/drm_boxes_fullDomain/drmdisplacements.0 --box 100 100 50 1 --dt 0.01
```

The reader reports each face's shape and offset and scans all values for finite
motion. For NumPy processing from the repository root:

```python
from pathlib import Path
import numpy as np
from script.read_drm import box_layout

layout = box_layout(100, 100, 50, 1)  # exterior dimensions and spacing
npoints = sum(rows * cols for _, rows, cols, _ in layout)
path = Path("outputfiles/drm_boxes_fullDomain/drmdisplacements.0")
xyz = np.loadtxt(path.with_name("drmcoordinates.0"), comments="#", ndmin=2)
if xyz.shape != (npoints, 3):
    raise ValueError("coordinate count does not match box dimensions")
nt, remainder = divmod(path.stat().st_size, npoints * 3 * 8)
if remainder or nt == 0:
    raise ValueError("incomplete snapshot or wrong dimensions")
u = np.memmap(path, dtype=np.float64, mode="r", shape=(nt, npoints, 3))
faces = {name: u[:, start:start + rows*cols, :].reshape(nt, rows, cols, 3)
         for name, rows, cols, start in layout}
face_xyz = {name: xyz[start:start + rows*cols].reshape(rows, cols, 3)
            for name, rows, cols, start in layout}
# faces["exterior_front"] -> (nt, 101, 51, 3)
# faces["interior_bottom"] -> (nt, 99, 99, 3)
# face_xyz["exterior_front"] -> (101, 51, 3), matching every time slice
```

Snapshots are taken before advancing steps divisible by
`output_drm_boxes_print_rate`, including step zero in fresh runs. The saved
interval is this rate times `simulation_delta_time_sec`. Every run truncates
its output, including restarts: use a new directory for a restart segment.
The first saved step is the first multiple of the print rate at or after the
restart step; multiply it by dt for the reader's `--first-time`.

Generated boxes use the plane GPU interpolation kernel, compact sampled-result
downloads, and precomputed collective MPI-IO views on solver ranks, including
empty ranks. There are no new field exchanges or full-field host copies. Each
file has independent buffers and sampling state. `drm_output_batch_size`
defaults to 16 snapshots per file; set 1 for writes at every saved step. A
partial final batch flushes before timing reductions. Buffer memory is
`24 * local_samples * batch_size` bytes per file, and
`3 * local_samples * batch_size` must fit an MPI int count.

Generated boxes work with `implement_drm=no` and are independent of stations,
planes, and the legacy `quake/forward/drm.cc` three-part replay workflow. The
standalone legacy `drm_boundary`, `drm_edgesize`, and `drm_offset_x/y` keys do
not define the new generated boxes. New displacement files cannot be passed
directly to the existing part2 reader. CPU geometry/reader and real MPI-IO
regressions are in `benchmark_inputs/drm_flat_output/README.md`; target-system
CUDA builds and station-history comparisons are still required before claiming
a physics/backend combination validated.
