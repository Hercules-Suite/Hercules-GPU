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
source-specific controls. Units are SI; at zero PROJ azimuth solver X is UTM
grid north and Y is grid east. Rotated runs use the local axes described below;
Z/depth is positive downward. Nonlinear examples require case-specific
calibration of geostatic state, tolerances, substeps, damping, and material
properties. Reported historical benchmarks apply to their stated cases and
versions; consolidating these documents does not validate new combinations.

## Contents

- [Domain rotation](#rotated-domains-and-material-grids)
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

## Rotated domains and material grids

In a PROJ build, required `region_azimuth_leftface_deg` is the finite angle
of solver +X clockwise from UTM grid north in the region origin's UTM zone.
Angles are normalized modulo 360 degrees. Zero retains X north/Y east;
+Y is always 90 degrees clockwise from +X and Z stays positive down.
`region_origin_latitude_deg` and `region_origin_longitude_deg` fix the
rotation pivot at solver `(0,0)`, not the centre of the rectangle.
`region_length_north_m` and `region_length_east_m` remain the local X and Y
side lengths; shallow/deep depths and octree sizing rules are unchanged.
UTM grid north is distinct from true geographic north.

Optional `3D_velocity_model_azimuth_deg` (finite degrees) defaults to
`region_azimuth_leftface_deg` when omitted, aligning the material axes with
the solver axes. An explicit value sets the material grid's +X orientation
independently in the same UTM zone; explicit `0` keeps UTM north/east axes.
The required `3D_velocity_model_origin_latitude_deg` and
`3D_velocity_model_origin_longitude_deg` locate its `(x=0,y=0)` point when
`3D_velocity_model=yes`. `x.in` and `y.in` contain metre offsets along its
own axes. Existing north/east grids require explicit `0`; a grid prepared
along the rotated solver axes can omit the material angle. The two origins may
differ. For domain angle $\theta$ and projected offsets $\Delta N,\Delta E$:

```math
\begin{aligned}
X &= \cos\theta\,\Delta N + \sin\theta\,\Delta E,\\
Y &= -\sin\theta\,\Delta N + \cos\theta\,\Delta E,\\
\Delta N &= \cos\theta\,X - \sin\theta\,Y,\\
\Delta E &= \sin\theta\,X + \cos\theta\,Y.
\end{aligned}
```

After subtracting the material origin expressed in solver coordinates,
rotate by the relative angle $\theta_{material}-\theta_{domain}$ into
material-grid axes. CPML material extension clamps the query in solver
coordinates first. Cached elevations at material columns use the inverse
rotation back into solver topography coordinates. Topography beyond its
horizontal coverage retains the existing edge-clamping behavior. Extend
material coverage past every intended query and its interpolation stencil;
the existing upper grid boundary remains exclusive.

Example additions to a complete case (illustrative, not scientifically
validated input):

```ini
region_origin_latitude_deg  = 37.791109
region_origin_longitude_deg = 26.383378
region_length_north_m       = 90000
region_length_east_m        = 90000
region_depth_shallow_m      = 0
region_depth_deep_m         = 39375
region_azimuth_leftface_deg = 30
cvmdb_input_file = profile
Istanbul_velocity_model = no
basin_velocity_model = no
3D_velocity_model = yes
3D_velocity_model_directory = input/materialfiles
3D_velocity_model_origin_latitude_deg  = 37.791109
3D_velocity_model_origin_longitude_deg = 26.383378
# Omit the next line for a material grid aligned with the domain.
3D_velocity_model_azimuth_deg = 0  # explicit override: UTM north/east grid
```

For these matching origins, solver `(1000,0)` samples material
`(866.0254,500)` metres with the explicit zero override. Omit the material
azimuth or set it to `30` for a grid prepared along the domain axes; the same
solver point then samples `(1000,0)`.
A valid background profile and the other required simulation/source inputs
are still necessary. Nonzero domain rotation currently requires the profile
background and `Istanbul_velocity_model=no`: legacy etree CVM and Istanbul
coordinate readers have not been migrated and are rejected. The general
3-D and basin models are supported; a PROJ basin grid stays aligned with
UTM north/east while queries rotate back into those axes.

Geographic stations, output-plane origins, point/SRFH source locations,
and planar-fault origins use the same mapping. Geographic fault strikes
are corrected from true north into solver axes, including UTM convergence.
Point/SRFH and kinked-fault strikes use quadrant-safe `atan2` directions.
The PROJ planar-fault path now uses projected origin coordinates instead of
the old constant metres-per-degree approximation. Cartesian source
coordinates and force components remain local. Output-plane strikes,
DRM-box strikes, and incident-wave directions are clockwise from solver +X.
Topography, groundwater, buildings, local DRM inputs, and directional
permeability use solver axes; prepare/resample these local grids for the
rotated domain because their files carry no independent geographic frame.

Output XYZ positions, displacement/velocity/acceleration vectors, and
strain/stress components stay in solver axes. PROJ writes
`coordinate_frame.txt` in each enabled station, plane, and structured DRM
output directory. It is overwritten at setup, including restarts, and
records `utm_zone`, `hemisphere`, `N0`, `E0`, `c`, and `s`, with
$c=\cos\theta$ and $s=\sin\theta$. Binary layouts are unchanged. For
`station_binary=yes`, the sidecar is in `output_stations_directory` even
when `station_output_binary_file` is elsewhere. A Python reader is:

```python
from pathlib import Path
frame = dict(line.split(" = ", 1) for line in
             Path("outputfiles/stations/coordinate_frame.txt").read_text().splitlines()
             if line and not line.startswith("#"))
c, s = float(frame["c"]), float(frame["s"])
def vector_to_utm(ux, uy, uz):
    return c*ux - s*uy, s*ux + c*uy, uz

def position_to_utm(x, y):
    return float(frame["N0"]) + c*x - s*y, float(frame["E0"]) + s*x + c*y
```

Convert tensor components with $T_{UTM}=A T_{solver} A^T$, using
$A=[[c,-s,0],[s,c,0],[0,0,1]]$. These are UTM grid components; conversion to
true north/east additionally requires meridian convergence. Inverse-project
positions using the recorded UTM zone and hemisphere for latitude/longitude.

Rebuild persisted meshes and regenerate source forces whenever either
origin or azimuth changes; previously stored material coefficients cannot
be corrected by relabeling output. Non-PROJ builds retain the existing
`domain_surface_corners` mapping and solver-aligned material axes when the
material angle is omitted or explicitly `0`; explicit nonzero overrides of
`3D_velocity_model_azimuth_deg` are rejected there. Local analytical geometry
and interpolation checks do not establish CUDA/MPI simulation validation.

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

For topographic meshes using `topography_surface_refinement=maximum_level`
(default), choose the forced surface target `maximum_octant_level` from the wavelength
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
This is not a global cap: the subsequent material wavelength rule can refine
further. Acoustic or other feature-specific rules can also require finer cells.

## Topography surface refinement

`topography_surface_refinement` accepts `maximum_level` (default) or
`wavelength`, case-insensitively. It is read only with `include_topography=yes`.
The default preserves forced refinement near the surface to
`maximum_octant_level` and the existing surface-size consistency check.
`wavelength` skips both and delegates sizing to the ordinary material rule,
as in the topography-off refinement path. Add this to an otherwise complete
topography case:

```ini
include_topography             = yes
topography_surface_refinement = wavelength
mesh_input_mode               = 0
# Keep maximum_octant_level and all other required topography/material inputs.
```

`maximum_octant_level` remains required (integer 0–30, no default) for geometry
classification/proximity bookkeeping even in wavelength mode. The startup log
reports the mode. Regenerate the mesh to change its sizing; importing an
existing mesh retains its stored cells.

For element edge length $`h`$, maximum frequency $`f_{\max}`$ and requested
`simulation_node_per_wavelength` $`n`$, the ordinary sizing targets are

```math
h \le \frac{V_S}{f_{\max} n}\quad\text{(solid)},\qquad
h \le \frac{V_P}{f_{\max} n\,r_a}\quad\text{(acoustic)},
```

where solid properties already reflect the existing shear-velocity floor and
$`r_a`$ is `acoustic_mesh_refinement_factor`. This option does not change
air/material classification, material interpolation, enabled building/DRM or
interface refinement, or octree 2:1 balancing. Those rules can produce finer
cells than the simple wavelength estimate.

Coarser surface cells reduce terrain resolution. FEM retains whole-cube
staircase geometry; VT computes cut volumes using the resulting element size.
Wavelength resolution alone does not establish surface-geometry accuracy.
Inspect surface representation, source/receiver ownership, actual timestep
limits and waveform convergence after rebuilding.

The decision is shared by elastic, viscoelastic, nonlinear elastoplastic and
poroelastic solids and supported acoustic–elastic, porous–water and
porous–elastic paths. Existing restrictions remain: acoustic field-specific
constraints require FEM with topography; local-matching rejects topography;
VT retains its existing pore-pressure constraint limitations. Acoustic free
surfaces still require supported flat, complete-face geometry; this option
does not add acoustic cut cells. Hydrostatic/geostatic initialization,
constitutive models, time integration, GPU kernels and field exchanges are
unchanged, but use the resulting mesh.

Local regression checks cover the parser, rank-broadcast serialization, and
production refinement/classification functions. They do not establish CUDA/MPI
propagation accuracy, coupled-physics accuracy or hydrostatic/geostatic
equilibrium for a coarser mesh. Those require case-specific comparisons.

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

### Scalar MPI schedules and profiling

`HERCULES_SCALAR_MPI_SCHEDULE=compact|full` selects the scalar communication
graph. `compact` is the default on both host and CUDA-aware MPI backends.
Setup retains nodes used by acoustic pressure/velocity potential, final active
BrickUP elements, open-pore target assembly, station/plane interpolation and
acoustic receiver recovery. It closes dependencies through the selected
field-specific hanging-node map (or the existing geometric fallback) and
propagates activity between owners and sharers. Neighbor global IDs and counts
are checked before filtering. A rank without local active cells still
participates when it owns or harbors a required degree of freedom.

The resulting separate acoustic, pore-pressure and open-pore-target schedules
retain original message/node order and duplicate entries. Empty peer messages
are removed; an entirely empty device schedule needs no stream synchronization.
Unbatched host scalar records contain one scalar plus the existing debug ID when enabled. Batched packets carry contiguous scalar segments using the setup-verified mappings.
Initialization exchanges and solid-vector schedules retain their full maps;
time levels, quadrature, constraint operators and physical boundary conditions
are unchanged. Host sparse GPU-transfer lists remain conservative supersets;
this optimization compacts MPI traffic, not all host/device transfers.
`full` restores the original general node mappings. Set both batching and
overlap to `0` as well to restore the sequential exchange path and legacy host
record format. All communication settings must be identical on all solver ranks.
Startup reports global original/retained node-entry and peer-schedule counts,
counting each contributor-to-owner edge once over dangling and anchored lists.

`HERCULES_SCALAR_MPI_BATCH=0|1` defaults to `1`. For sealed acoustic–poroelastic
coupling with the same transport backend for both fields, it assembles both
RHS vectors before exchanging them together, then updates and shares both
fields together. Including the unchanged four solid phases, the dependency
sequence goes from 12 to 8 exchange phases per step. Open-pore coupling must
finish the acoustic field first because the hydraulic target reads that new
field. It combines only the pore RHS and target contributions, reducing 14 to
12 phases. A packet concatenates each field's existing compact mapping; it
does not pad one field to the other field's nodes. Each shared neighbor receives
one message per phase. Neighbors used by only one field still receive that
field. Peer unpack order, duplicate-node handling and field-specific hanging
constraints are preserved. Batching works with host and CUDA-aware MPI and
with `compact` or `full` schedules. Different acoustic/pore backends retain
separate exchanges. Phase counts describe dependencies, not a runtime speedup.

`HERCULES_GPU_MPI_OVERLAP=0|1` defaults to `1`. On CUDA-aware scalar paths, the
anchored RHS exchange starts before an independent interior pressure update is
queued. Setup partitions all harbored nodes into two disjoint lists: nodes
that incoming anchored contributions can modify, and all remaining nodes.
The pore partition also includes incoming hydraulic-target destinations when
that target is batched. Dangling contributions and constraint distribution
finish before this split, so the interior update cannot read an unfinished
RHS. Receive completion, ordered unpack and boundary updates follow. These
launches use the same node-local formulas, including free-surface, CPML,
drainage and open-pore terms. No interior update overlap is attempted on the
host transport; its existing staging and fallback behavior remain available.
Consecutive host pore-RHS/target uploads protect their shared pinned staging
buffer with a completion event before the CPU rewrites it.

Both scalar batches and solid-vector device exchanges omit the immediate
host wait after GPU unpack when overlap is enabled. A CUDA event protects
receive-buffer reuse on the next call; MPI send completion separately protects
send-buffer reuse. Packing still waits for GPU readiness before MPI sends,
and subsequent field consumers remain ordered on the main CUDA stream.
The solver drains the final GPU tail before collecting solver timings.
There are no per-step buffer allocations or additional MPI collectives.
This implements overlap for scalar node updates, not a split of all solid
stiffness or nonlinear element kernels. Disable overlap with `0` to use
whole-field pressure updates and synchronous unpack; disable batching as well
for the previous sequential path. Neither switch alone skips diagnostic pore pressure;
how often it is updated is controlled separately below.

`HERCULES_PORE_PRESSURE_UPDATE=output|step` defaults to `output`.
It applies only to fresh linear runs where the undrained fast path is active:
zero permeability, undrained boundaries, no open-pore coupling and no nonlinear
analysis. Since fluid stiffening is already in the solid operator, pore pressure
is diagnostic and its linear increment equation telescopes from the initial
zero displacement/excess pressure. At a required state the solver reconstructs
`p = -M_lumped^-1 Q^T u` using the existing element quadrature, assembled mass,
MPI contribution maps and hanging-node constraints. It does not accumulate
unexchanged increments or allocate another displacement history. Floating-point
rounding can differ from the incremental recurrence.

All ranks reconstruct the state preceding an ASCII pressure station output,
even ranks without local stations. Only stations in poroelastic
elements trigger this work. Displacement-only binary stations and plane/DRM motion
do not consume pore pressure. The final two pressure states are also reconstructed.
Intervening steps skip pore RHS assembly, update and all four standalone pore
exchange phases. With sealed scalar batching, acoustic-only packets are used on
skipped steps and combined acoustic/pore packets on reconstruction steps; the
four combined phases themselves still occur for the acoustic field. Both host
and CUDA-aware MPI, full/compact maps, and batching/overlap controls remain usable.

`step` retains the original incremental update on every step. Checkpoint reads
or writes automatically select that path regardless of the request, preserving
both saved pressure time levels and the existing checkpoint format; no checkpoint
schema or restart-state migration is introduced. Nonlinear, diffusion, drained
and open-pore paths always retain every-step pressure feedback. Pore-pressure-ratio
(`Ru`) output currently requires nonlinear analysis and an initial effective-stress
reference, so it also retains every-step updates. The selected
update mode and output consumers are reported at startup, and reconstructed/skipped
step counts after the loop. This control must agree on all ranks. It is separate
from `poroelastic_undrained_fast_path=no`, which changes the mechanical formulation.

`HERCULES_GPU_MPI_PROFILE=0|1` defaults to `0`. Set `1` consistently on all
solver ranks to record device-helper and scalar-batch calls, send/receive bytes
and messages, and CPU elapsed times for posting receives, packing/waiting for
GPU readiness, posting sends, waiting for receives, unpacking, and waiting for
sends. End-of-loop output reports cumulative rank average/min/max for solid
force, displacement, acoustic RHS/field, pore RHS/field, open-pore target, and
combined acoustic+pore RHS/field or pore RHS+target packets when active. The
combined counters count actual packets and their concatenated bytes once.
It adds no per-step collective, allocation or GPU synchronization. With overlap,
unpack time measures enqueue work; receive-wait time can overlap an interior
GPU update. These CPU timings do not isolate kernel duration or network transfer.
Host schedule timers and dangling communication counters cover only the host helper;
batched host traffic is in the scalar-batch profile instead. Scalar exchanges
remain nested under `Compute Physics`, so top-level `Communication` is not total
MPI cost. The pipeline, RHS/field exchange, and interior/boundary enqueue timers
are reported separately; nested timings must not be added together.

For an identical-mesh comparison, preserve partition, timestep, source and
output cadence and run each combination of `compact|full` and
`host|cuda_aware`, for example:

```bash
export HERCULES_SCALAR_MPI_SCHEDULE=compact
export HERCULES_SCALAR_MPI_BATCH=1
export HERCULES_GPU_MPI_OVERLAP=1
export HERCULES_GPU_MPI_PROFILE=1
export HERCULES_GPU_MPI_BACKEND=cuda_aware
BASE_IN=/absolute/path/parameters.in PSOLVE_BIN=/absolute/path/quake/forward/psolve ./run_phaseA_phaseB.sh
```

This change requires target-system CUDA/MPI propagation and performance
comparison before claiming a speedup. Check standalone and coupled elastic,
acoustic and poroelastic cases, sealed/open-pore interfaces, refinement across
ranks, ranks without local active cells, both acoustic scalar formulations,
viscoelastic/nonlinear solid paths, station/plane pressure and motion outputs,
and restart. The linear undrained diagnostic-pressure update still runs each
step. Existing limits remain: poroelastic–elastic supports `no_flow`, saturated
CPML candidates remain elastic CPML, and VT pore constraints use the geometric
fallback. This change does not add poroelastic CPML or validate those limits.

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
	Default unset/`0`: allow the optimized BKT-family kernels, including
	the fused BKT2 path below. Set to `1` to force the generic per-element
	`kernelDampingCalcConv` for all BKT-family damping, bypassing the
	fast kernels. Use it to produce a generic-kernel reference from an
	otherwise-identical run and verify the fast kernels — they must match
	the generic kernel to single-precision tolerance for a
	spatially-uniform relaxation spectrum (the auto-detected common
	case). A one-time rank-0 line is printed to the monitor/stderr when
	the override is active so a forced run cannot be mistaken for a
	normal one. Read once at startup; production runs leave it unset.
- `HERCULES_BKT2_FUSED=0|1`
	Default unset/`1`: use shared histories and fused convolution/force
	for eligible BKT2 cube elements with device-resident force and displacement
	synchronization. Host-staged paths retain the reference kernels. `0` retains
	the previous separate
	hoisted-coefficient convolution and force kernels. The generic override
	above takes precedence. Both switches are read once; never switch history
	layouts during a run. See [BKT2 GPU execution](#bkt2-gpu-execution).
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
  store the per-point scale. Original models retain marker 4. Both
  `RIVASand02` and `RIVASandCustom02` now default to:

  ```ini
  rivasand02_reversaltype = 3
  rivasand02_biasvolume = 1
  ```

  These optional integer scalars belong in `parameters.in` or a material-point
  input, not in a material row. Invalid values, duplicates, and explicit use
  with other models are rejected. Reversal type 1 compares adjacent valid
  increments; type 2 uses accumulated loading-branch strain. Type 3 uses the
  unguarded full-host elastic stress-ratio predictor from OpenSees_RIVA
  `355d8f042`, evaluating
  $(\alpha_{trial}-\alpha_0):(\alpha_{trial}-\alpha)<0$ from the committed
  mechanical state. It applies at most one event, in the first substep.
  Before cyclic activation it retains the strain-based initialization rule.
  `reversalGuard` and `reversalLatch` are not exposed.

  BiasVolume modes follow that same OpenSees revision:

  | Mode | Behavior |
  | --- | --- |
  | 0 | Preserve the original bias-driven reversible-volume contribution. |
  | 1 (default) | Disable that contribution; retain separate phase-transformation volume. |
  | 2 | Preserve it with the upstream correction to its contractive mean part. |

  Mode 2 gates the correction by fixed bias, reference-pressure ratio,
  effective-pressure ratio, and accumulated plastic activity. It evolves
  `field_bias_mean_activity` once after the host increment, keeping that
  activity fixed during local substeps. This additional scalar is zero in
  modes 0/1 and is the 139th flattened state value. Mode 2 remains an optional
  research choice; selecting it does not establish field stability or convergence.
  None of these modes disables BrickUP fluid-pressure coupling or the separate
  phase-transformation mechanism.

  RS02 checkpoint marker = `10 + 3*BiasVolume + (reversaltype-1)`; the defaults
  use marker 15. Each of the nine combinations has its own marker. The new
  139-value layout rejects all older RS02 markers through 9 and changes to
  either selector on restart. Rerun initialization, keeping model, mesh, and
  material inputs unchanged. To reproduce an earlier response rule, explicitly
  select mode 0 and the earlier reversal type in a fresh run. The 257-value
  upstream reference calibration has SHA-256
  `16c8567fd2370c562899b9ec00cc2fc313be137f239cde5d0e22cc1cd3eee800`;
  this identifies the reference record, with runtime selector overrides
  recorded separately. User material-row values are unchanged.

  On a rejected RS02 GPU update, the solver writes
  `rivasand02_failure.rank<R>.step<S>.txt` in the run directory before aborting.
  This is a diagnostic snapshot, not a restart checkpoint. A repeat failure
  at the same rank and step overwrites that filename. The download and file
  write occur only after failure; no per-point history transfer is added to
  successful timesteps. The log also reports the full incoming tensor strain
  increment and cached `E_ref` in Pa. Non-finite strain increments are rejected,
  including NaNs that must not be mistaken for a zero increment.

  The text starts with `RIVASAND02_FAILURE_V2`, followed by `key=value` records:
  model, parameter hash, kernel/checkpoint revisions, rank, element, quadrature
  point, step, reversal type, BiasVolume mode, substeps, whether an update was attempted,
  initialization/admission flags, compatibility residual, increment, material,
  and state. For example, `increment=0 0 0 0 0 1e-6` uses physical tensor order
  `[xx yy zz xy yz xz]`. The eleven `material` entries are
  `h m M kd zeta e_max e_min Q R n_G E_ref`; `E_ref` uses Pa.
  `state_values=139` identifies the subsequent `state` array in the native
  `riva_ib_state_values` order. A rejected constitutive call preserves its
  input history; a failure of the subsequent compatibility check instead
  records the resulting state. Thus the printed pressure and residual alone
  do not identify the rejected trial. Investigation replay helpers remain
  local-only under `.local-diagnostics/`.

  The material-point driver takes Dr [G0] or the twelve Custom values,
  without Vs_limit. Automatic mode additionally requires
  material_point_riva_density (kg/m3) and material_point_riva_vs (m/s).
  Full reference histories are retained using explicit reference G0.

  The 02 kernel retains the plastic-activity and material-specific density
  corrections, with kernel revision 5 and the 139-value BiasVolume history.
  OpenSees reversalLatch is not included. Both generations
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
	Regenerate meshes made before the DRM Qs/Qp material-payload extension;
	current sharded manifests and MPI-IO archives use format version 2.

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
- `acoustic_interface_method`
	Optional, case-insensitive: `field-specific-constraints` (default for coupled
	scalar acoustics), `shared-geometric-constraints`, or `local-matching`. Select the geometric constraints,
	experimental field-specific constraints, or experimental automatic face
	matching during meshing. This is independent of quadrature, face filters,
	gravity, and CPML. See [details and limits](#interface-constraints-and-local-matching).
- `acoustic_interface_projection`
	Optional. `point` (default) uses the legacy one-point quadrature on each
	acoustic-elastic interface fragment. `face_mass` uses a 2x2 face quadrature
	to weakly project the normal velocity/pressure traction over each fragment
	before assembling nodal interface forces. This is an experimental diagnostic
	option intended to reduce grid-scale interface forcing in
	`acoustic_formulation=velocity_potential` runs. PML fragments retain the
	selected quadrature, with independent cofactor histories at each point.
- `acoustic_interface_face_filter`
	Optional. `all` (default), `bathymetry`, or `bottom_only`. This applies when
	acoustic interface coupling is enabled. `all` couples every
	acoustic-elastic interface face using mesh-face normals.
	`bottom_only` couples only acoustic `Z+` faces, which correspond to the
	water-bottom interface in the depth-positive mesh. This is a diagnostic
	option for bathymetry/shoreline cases. Differences between filters can help
	isolate the effect of lateral wet/dry stair-step shoreline coupling, but
	first verify that compared runs use identical, finite source forces.
	`bathymetry` also couples only acoustic `Z+` faces, but reconstructs a
	local bathymetry normal from nearby acoustic-to-elastic transition depths
	and scales each projected bottom patch by the reconstructed surface area.
	This restores inclined-bottom normal coupling without reintroducing raw
	stair-step side-wall faces.
	Omitting the parameter or using the `default` alias selects `all`.
	Existing explicit `bathymetry` and `bottom_only` settings retain their
	behavior. Other accepted aliases are `reconstructed_bathymetry` and
	`bathymetry_normal` for `bathymetry`, `none` for `all`, and `bottom` or
	`z_plus` for `bottom_only`. Values are case-insensitive.
	`bathymetry` omits lateral and upward-facing acoustic-solid interfaces, so
	check that this representation is appropriate for the model geometry.
	Stability with this filter does not by itself validate waveform accuracy.
	The source-corrected Izmir comparison completed 100 s with finite output
	for all three filters using identical generated forces; slow acoustic
	motion and reconstruction accuracy remain under investigation.
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
	traction and normal-motion coupling only on those fragments. Each fragment
	interpolates the scalar acoustic field on the acoustic element and solid
	normal acceleration (`pressure`) or velocity (`velocity_potential`) on the
	adjacent solid element, so refined solid
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
	hydrostatic prestress solve. With `seabed_geometry_file`, complete horizontal
	acoustic top faces must coincide with local topography over water, or solver
	Z=0 with `include_topography=no`. Separate
	flat water bodies may have different elevations. Without explicit seabed
	geometry, the flat acoustic surface is at solver
	`Z=topographybase_zcoord` with topography, `Z=surface_shift_m` with buildings,
	or the octree top otherwise. A Hex8 top face must coincide with that datum;
	with topography its center must also be exposed at the topographic surface,
	so buried water beneath land is excluded. This supports flat water surfaces,
	not a general curved or cut-cell fluid free-surface operator. The same
	detection selects the pressure-release nodes when this option is `no`.
	The startup `Acoustic free-surface geometry` line identifies local topography
	or the legacy solver datum, then reports the exact exposed face count and
	summed harbored surface-node copies (the last
	count includes MPI replicas). Where this face intersects a
	lateral acoustic CPML, the surface mass uses the tangential Jacobian
	`kappa_x*kappa_y`; the solver also adds the paired `x/y` mass-proportional
	damping and CFS restoring/memory terms on the four top-face nodes. The
	vertical CPML stretch is intentionally excluded because `z` is normal to the
	gravity surface. Corner PML faces include both tangential directions. This
	adds no per-step work to non-surface PML elements. Each top node has its
	own memory, consistent with the lumped surface mass. Changing this
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
- `acoustic_station_gradient`
	Optional, case-insensitive `linear` (default) or `element`. Applies to
	acoustic station motion in both `pressure` and `velocity_potential` modes.
	`linear` combines adjacent same-material water-cell derivatives to recover
	the pressure gradient at the actual station position before the existing
	acceleration and time integration. `element` keeps the historical single-cell
	Hex8 derivative. Pressure interpolation and solver propagation are unchanged.
	ASCII/binary station displacement and requested velocity/acceleration share
	this setting. Plane reconstruction does not use it; see
	[Acoustic station gradient recovery](#acoustic-station-gradient-recovery).
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
	file after solver setup finalizes BrickUP and CPML membership, and rewrites
	it after a ramped geostatic reference refresh. It records each station's
	owning medium (`acoustic`, `poroelastic` or `elastic`), pressure eligibility, owner rank,
	element id/level, local coordinates, and material properties. If unset, it
	defaults to `output_stations_directory/station_mediums.txt`. The metadata
	includes `sigma_v0_eff_pa`, the positive station-specific denominator used
	by optional `Ru(-)` output; it is `nan` where the initial effective-stress
	reference itself is unavailable or nonpositive.
	Poroelastic pressure is defined only when the owning element belongs to the
	final active BrickUP map and the receiver is saturated. Elastic-only CPML
	cells are labeled `elastic` with pressure undefined, even if shared nodes
	carry pressure from adjacent active soil. Their `Pressure(Pa)` and `Ru(-)`
	columns are `nan`; skeleton motion still uses the normal solid output path.
- `acoustic_quality_factor`
	Optional, dimensionless. Default `0`; omitting the parameter or setting
	it to `<= 0` disables interior acoustic damping. This default was verified
	in the documented flat-water CPS comparisons and MPI/backend consistency
	tests; see [validation scope](#interface-constraints-and-local-matching).
	Those results apply to the recorded configurations and durations.
	Positive values enable stiffness-proportional numerical damping referenced
	to `simulation_wave_max_freq_hz`, with stronger damping at higher frequencies.
	It works for `pressure` and `velocity_potential`, skips CPML elements, and is
	independent of material-profile `Qp`/`Qs`. CPML absorption and energy transfer
	into coupled solids remain separate from this control.
- `acoustic_poroelastic_interface`
	Optional. Default `no_flow`, which keeps an acoustic-poroelastic contact
	sealed for the pore fluid while the existing acoustic pressure traction
	still acts on the solid skeleton. Set to `open_pore` to add weak
	compliance-projected pressure continuity on detected acoustic-poroelastic
	interface fragments, driving pore pressure toward the adjacent acoustic
	pressure using four-point (2x2 Gauss) hydraulic quadrature, independently of
	`acoustic_interface_projection`. The interface weight uses the adjacent poroelastic element's
	`1/M`, matching the heterogeneous interior pressure capacity. This requires
	`include_acoustic=yes`, the acoustic scalar formulation, and
	`include_poroelastic=yes`.
- `acoustic_poroelastic_pressure_penalty_factor`
	Optional. Default `1.0`. Multiplies the open-pore weak interface compliance
	used by `acoustic_poroelastic_interface=open_pore`. The target exchange
	uses the existing device-resident scalar schedule when CUDA-aware MPI is
	active, with sparse host staging as the fallback.

### Open-pore hydraulic quadrature

`acoustic_poroelastic_interface=open_pore` always uses four-point (2x2 Gauss)
integration on each detected rectangular interface fragment. No new input is
required. `acoustic_interface_projection=point|face_mass` still selects only the
mechanical traction/motion quadrature. At each hydraulic point, the acoustic
Hex8 trace supplies pressure and the poroelastic Hex8 trace supplies the nodal
load weights. Each side uses its own fragment center and half widths, so
nonmatching elements sample the same physical points. Remote solid vertices
are resolved at setup for every open-pore path; missing trace geometry aborts
instead of falling back to acoustic nodes or centroid integration.

For fragment $f$, let $A_f$ be its existing effective area, $h_s$ the adjacent
solid edge length, $M_s$ its storage modulus, and $\gamma$ the existing positive
`acoustic_poroelastic_pressure_penalty_factor`. With four equal normalized
weights, the boundary capacity $b_j$ and pressure target $r_j$ receive

```math
c_f = \frac{\gamma A_f h_s}{2M_s}, \qquad
b_j \mathrel{+}= \frac{c_f}{4}\sum_{q=1}^{4}N_j^s(x_q), \qquad
r_j \mathrel{+}= \frac{c_f}{4}\sum_{q=1}^{4}N_j^s(x_q)
                 \sum_i N_i^a(x_q)p_i^a.
```

The setup capacity and runtime target share the same quadrature rule. For a
constant pressure, $r_j=p^a b_j$; on a rectangular fragment the capacity agrees
with the former midpoint capacity within roundoff. Varying pressure now changes
the nodal targets: on a matching face with $p^a=200+100\xi$ Pa, the two nodes
at $\xi=-1$ have $r_j/b_j=166.666\ldots$ Pa and those at $\xi=+1$ have
$233.333\ldots$ Pa, instead of assigning every node the 200 Pa centroid value.
For bilinear face traces and constant fragment compliance, 2x2 Gauss integrates
the shape-pressure product exactly on the represented rectangle.

The pressure formulation supplies its existing acoustic pressure time level;
velocity potential still uses $p_i^a=-(\phi_i^{\mathrm{Tm2}}-
\phi_i^{\mathrm{Tm1}})/\Delta t$. The pressure update, penalty scaling,
material parameters and time integration are unchanged. The kernel reads each
acoustic nodal value once and accumulates the four points before one atomic
addition per loaded solid node. It uses the existing CUDA-aware or host-staged
scalar exchange without extra per-step MPI rounds. The startup monitor names
the 2x2 Gauss hydraulic rule separately from mechanical projection.

This changes the spatial integration of the existing weak pressure-continuity
operator; it does not introduce a new fluid-inertia or Darcy interface law.
It uses the existing mesh-face geometry and area convention, without adding a
cut-face or curved-surface reconstruction. VT pressure constraints still retain
their documented geometric fallback. Earlier Horizon open-pore comparisons
below used centroid hydraulic integration and do not validate this upgrade.
Local analytical host checks cover the production capacity/target algebra,
constant and varying pressure, both acoustic formulations, six face orientations,
nonmatching fragments and reordered remote trace indices. CUDA/MPI propagation
and performance validation on the target system are still required.

### Interface constraints and local matching

`acoustic_interface_method` accepts three values. Omitting it selects
`field-specific-constraints` for coupled pressure or velocity-potential acoustics.
If acoustics, scalar acoustics or interface coupling is disabled, an omitted method
uses `shared-geometric-constraints`; explicit unsupported choices still abort.
Coupled VT topography requires explicit `shared-geometric-constraints`.
Selecting a constraint method does not change `acoustic_interface_projection`
(`point` or `face_mass`), `acoustic_interface_face_filter`,
`acoustic_free_surface_gravity`, or `abc_type`.

| Value | Behavior |
| --- | --- |
| `shared-geometric-constraints` | Acoustic and solid fields share the geometric hanging-node constraints. |
| `field-specific-constraints` | Default for coupled scalar acoustics. Separate maps retain only constraints required by an incident coarser element of the same field. Nodes constrained solely by the other material become independent for their own field. |
| `local-matching` | Experimental. During fresh generation, refine the coarser leaf at every detected acoustic–solid interface, rebalance, re-extract and reclassify until opposing faces match in size. Use geometric acoustic/solid constraints on the resulting mesh; pore pressure retains its independent map. |

Both alternatives accept either scalar acoustic formulation. In `pressure`
mode the interface applies pressure traction to the solid and solid normal
acceleration to the acoustic equation. In `velocity_potential` mode it obtains
traction from the potential time difference and drives the acoustic equation
with solid normal velocity. The field constraints and matching rules are the
same for both. Selecting pressure does not remove the need to check time-step,
mesh and interface accuracy; it does not change what `acoustic_station_gradient`
affects (station motion, not the pressure column). Velocity potential remains
the default.

Pressure-formulation validation (2026-10-08 UTC): a local 58-run CUDA/MPI
comparison exercised all three methods. All 58 runs used
`acoustic_quality_factor = 0`. Fourteen 110 s flat-water runs covered
coarse/fine meshes, both formulations and two half-time-step pressure controls.
Forty 4 s source-excited checks covered `point`/`face_mass`, gravity, finer water,
buildings combined with synthetic `vonMises_ep` nonlinear elastoplastic and
BrickUP open-pore coupling, one/eight ranks, both MPI backends and existing
pressure controls. Four 6 s real-terrain Izmir crop runs checked field-specific
FEM topography with gravity, eight/sixteen ranks and both backends. Matched
MPI/backend/formulation comparisons agreed below $10^{-6}$ aggregate relative
L2 in motion and acoustic pressure; available pore-pressure, strain, stress and
nonlinear station-state columns were also checked. Undefined elastic pressure
NaNs were excluded, not other output components.

Against the independent source-matched CPS reference, median submerged-receiver
pressure errors over $10 \le t < 100$ s were about 18.38% on the coarse mesh and
8.34% on the fine mesh for both newly enabled methods and both formulations.
Halving the coarse time step from 0.005 to 0.0025 s changed pressure error to
about 18.40%; changing the evolved variable did not improve pressure accuracy
in these zero-initialized tests. No simulated response was filtered, fitted or
shifted; the benchmark's prescribed source filter was retained. These local-only
fixtures are unavailable in a clean checkout. Gravity/module/topography checks
establish implementation consistency within the stated durations, not
independent gravity-wave accuracy, calibrated constitutive response, or
full-duration Izmir stability. Existing geometric restrictions below still apply.

The two alternatives require `include_acoustic=yes`,
`acoustic_interface_coupling=yes`, `acoustic_formulation=pressure` or
`acoustic_formulation=velocity_potential`,
and full Hex8 solids. Both alternatives accept `include_buildings=yes`,
`include_nonlinear_analysis=yes` and `include_poroelastic=yes`, separately or
together, subject to the existing restrictions of those modules. Nonlinear
geostatic initialization, tolerances, substeps, damping and
`material_properties_list` remain case-specific calibrated inputs. This change
does not select a constitutive model or supply production material values.
Field-specific constraints also
support `include_topography=yes` with `computation_method=fem`: retained
material cells still use regular Hex8 traces. Air cells (`Vp=-1`, zero density)
and carved voids carry no field, including when their MPI owner is remote.
Both `type_of_etree=sqd` and `full` retain their existing material-query meaning;
the constraint map uses the actual mesh cells, not a new terrain interpolation.
VT cut-tetrahedron topography remains unsupported by field-specific constraints.
Local matching still requires `include_topography=no`. Stepped or sloping
**material-defined** seabeds can be used within the full-Hex8 scope.
Unknown values and unsupported combinations abort with an explanatory message.

For an existing FEM-topography case such as Izmir, keep its topography datum,
directory, etree type and other physical inputs, and select:

```ini
include_topography = yes
computation_method = fem
acoustic_interface_method = field-specific-constraints
```

This choice does not remesh an imported mesh. Acoustic volume refinement remains
optional and also works with shared geometric constraints. Changing the volume
refinement factor requires fresh meshing to change the existing resolution.

Field maps resolve the actual incident leaves through batched MPI owner
queries. Every retained stencil must equal the corresponding same-field Hex8
trace; unsupported mixed-corner stencils abort rather than approximate it.
The same maps handle mass, damping, free-surface flags, force distribution and
field assignment. Device arrays are uploaded once at setup. Existing CUDA
adjustment kernels, CUDA-aware exchanges, host-staged fallback and owner
schedules remain in use. Independently of `acoustic_interface_method`, a separate
pore-pressure map is built for regular Hex8 and FEM-topography BrickUP models,
including standalone poroelastic–elastic runs with `include_acoustic=no`.
It uses final active BrickUP cells, excluding elastic-only CPML cells, and handles
pore-pressure capacity, drained and open-pore boundary weights, pressure RHS
distribution and pressure assignment on both MPI backends. Saturated solids
carry both displacement and pore pressure; dry solids and buildings still carry
solid displacement. The monitor reports acoustic, solid and active pore-pressure
keep/release/absent counts and identifies FEM topography when selected. Missing
remote replies still abort; only an owner-confirmed carved void is omitted.
Coupling uses the actual solid face vertices on both sides of MPI partitions.
Missing remote vertices are added as non-owned harbored nodes during setup;
their global IDs order the existing force-contribution and displacement-sharing
schedules. Both alternative methods use the solid trace for `point` and
`face_mass`. All three methods, including `shared-geometric-constraints`,
resolve actual remote solid vertices when `face_mass` is selected and use the
same four-point face quadrature on local and remote neighbors. They no longer
fall back to point integration on remote faces. Open-pore boundary mass and
water-pressure target loads also use actual solid vertices and their shape functions.
The [hydraulic projection](#open-pore-hydraulic-quadrature) uses four-point
quadrature independently of the mechanical `point`/`face_mass` selection. In particular,
independent fine-solid nodes must not be replaced by a coarse acoustic trace
when the neighboring solid element is remote. This also applies without
topography. The monitor reports remote fragments and additional harbored nodes.
The shared-geometric method retains its geometric hanging-node
constraints for acoustic and solid fields. Its `face_mass` correction
uses the same trace and global-ID schedule ordering as the alternatives,
including on ranks that only own nodes requested by another rank. Coupling-off
and non-acoustic cases do not build these traces.
`open_pore` also builds actual solid traces with
`shared-geometric-constraints + point`: hydraulic quadrature requires its
own solid-side load weights on remote faces. This uses the existing global-ID
MPI schedules. Sealed shared-geometric `point` coupling keeps its existing path.
This is a constrained experimental implementation,
not a general mortar method.

Horizon checks on 2026-10-04 covered an 8 s FEM case with the water surface
shifted below an air region, CPML, gravity on/off, both refinement directions,
`point`/`face_mass`, and mesh import. One versus eight GPU ranks and host versus
CUDA-aware MPI agreed within $2\times10^{-11}$ relative L2 across the tested
plane and station histories. The previous remote-solid fallback produced a
1.43% plane-history difference between one and eight ranks in the solid-finer
case. The corrected implementation also reproduced the original one-rank
topography-free operator and the shared-geometric baseline. These comparisons
check implementation consistency; they do not establish continuum accuracy or
long-time stability for arbitrary coastlines.
The same build imported the existing 1,179,872,281-element Izmir mesh on 512
Horizon GPUs and completed 20 s with FEM topography, gravity and the original
physical/source controls. All 18 stations recorded nonzero, finite motion
(largest displacement component 0.1384 m); station histories and generated
source forces passed output-integrity checks. Setup resolved 34,534 remote
solid fragments using 12,310 additional harbored nodes. This is a short-run
integration check, not validation of the full 100 s production simulation.

The subsequent module extension was checked on Horizon with 66 successful
candidate runs and six historical solver controls, yielding 49 accepted
physical-output comparisons. Four-second synthetic cases covered both methods,
buildings, yielding `vonMises_ep` solids, explicit BrickUP no-flow/open-pore
coupling, combined modules, both refinement directions, `point`/`face_mass`,
one/eight GPU ranks, host/CUDA-aware MPI, nonlinear CPU/GPU execution, drained
pore pressure with finite Darcy permeability, and matching-mesh import.
The final build also replayed the previous 8 s FEM/gravity baselines. Raw
responses were not filtered. Across the accepted comparisons, the largest
relative L2 differences were $2.02\times10^{-9}$ for an individual plane-motion
component and $2.50\times10^{-8}$ for a physical station column. The existing
CPU pre-arrival skip leaves the diagnostic `Fs` at its initialized zero where
GPU execution evaluates its negative elastic value; that diagnostic was
excluded from CPU/GPU comparisons, with no motion, pressure, strain or stress
column excluded. These synthetic materials are not production calibrations,
and these checks do not validate every constitutive model or long-time stability.

In the partition-crossing `local-matching + face_mass` control, the worst
plane-component one/eight-rank discrepancy dropped from 3.33% before the
remote-quadrature fix to $1.85\times10^{-13}$ relative L2 afterward. The new
one-rank result also reproduced the original one-rank operator. The final
follow-up change only guards the absent carving timer when printing timing
for imported building meshes; six import checks and the two FEM/gravity
regressions passed with that final build. Earlier failed or incorrectly
specified test attempts were excluded and replaced by corrected runs.

The later default-method correction explicitly compared
`shared-geometric-constraints + face_mass` on one/eight Horizon GPUs. In 4 s
flat-water cases, the worst plane-component relative L2 discrepancy changed
from 3.33% to $1.93\times10^{-13}$ on a matching interface, from 2.78% to
$1.81\times10^{-13}$ with the solid twice as fine, and from 1.54% to
$7.67\times10^{-14}$ with water twice as fine. The corrected one-rank results
reproduced the previous one-rank operator. These are raw, unfiltered outputs.
Checks also covered host/CUDA-aware exchanges, the unchanged default `point`
path, matching and nonmatching building/nonlinear/open-pore combinations,
8 s FEM/gravity cases, the legacy pressure formulation, coupling disabled,
and regressions of both alternative methods. The default-method correction
was not part of the earlier production-sized Izmir check.

One combined yielding/open-pore case with finer solids exceeded the strict
host/CUDA comparison tolerance in host performance mode: plane displacement
differed by $4.67\times10^{-8}$ relative L2, and the largest station-column
difference was $3.84\times10^{-5}$ in acceleration. Host repeats also varied;
the previous solver showed a host/CUDA difference on the same one-rank model.
With the existing `HERCULES_REPRO_MODE=strict`, host/CUDA and one/eight-rank
plane comparisons agreed within $1.5\times10^{-12}$, and station comparisons
within $4.6\times10^{-9}$. The original performance-mode failures are retained
as a repeatability limitation, not counted as strict passes. No host execution
defaults were changed by the face-integration correction.

Local matching inspects all detected water–solid faces, including those omitted
by a force face filter, and resolves refinement requests on their MPI owners.
There is no hardcoded seabed, band width, or extra target-level parameter.
`acoustic_interface_max_level` still refines mixed material cells and
`acoustic_mesh_refinement_factor` still controls the acoustic volume target.
Matching only adds refinement, and matching passes can propagate along an
interface and through octree balancing. Memory, setup and solver costs can
increase substantially. A 128-pass guard rejects failure to reach a matching
mesh. Reclassification of refined mixed cells can change the represented
interface geometry; compare material volumes and resolution as well as traces.
When buildings are enabled, each matching pass reapplies building carving after
refinement and balancing before mesh extraction and material reclassification.
Building meshing now also handles layered profiles: building/foundation material
and air classification precede the soil query, avoiding invalid negative-depth
profile queries above `surface_shift_m`. Air uses a positive meshing speed
(`Vp` for a water surface with `Vs=0`) before carving. Building and foundation
properties are selected using the actual cell depth. BrickUP ground and
groundwater depths use the building surface shift, so above-ground building
cells do not become saturated merely because their solver Z coordinate is
positive. Existing groundwater and depth-limit settings determine whether
below-ground foundation cells carry pore pressure.

Set `local-matching` before **Phase A**, regenerate the mesh and source caches,
and use the same setting for Phase B. An imported mesh is checked and rejected
if detected interface neighbors differ in size; import does not remesh.
Matching-pass counts and completion are printed in the monitor. This automatic
face-matching method differs from the earlier experimental fixed-width
refinement band around the benchmark seabed.

```ini
include_acoustic = yes
acoustic_formulation = velocity_potential
acoustic_interface_coupling = yes
acoustic_interface_method = field-specific-constraints
# Default: field-specific-constraints; alternatives: shared-geometric-constraints, local-matching
acoustic_interface_projection = face_mass
acoustic_interface_face_filter = all
acoustic_free_surface_gravity = yes
abc_type = cpml
cpml_clamp_outer = yes
include_topography = no
include_buildings = no
include_nonlinear_analysis = no
include_poroelastic = no
```

**Validation scope:** the preceding isolated, fixed-geometry sloping-seabed
prototype with separate field maps, `face_mass`, gravity and CPML remained
bounded for 80 s at two resolutions. That establishes observed boundedness,
not waveform accuracy. Separate dashpot comparisons did not establish improved
accuracy for those maps. Earlier band-matched sloping runs developed exponential
growth with CPML, including `face_mass`; stable flat-water results do not
establish general CPML stability. Neither experimental mode is a validated
replacement for every geometry. Keep source, timestep, absorber thickness in
metres, gravity and receivers fixed in comparisons, and recheck the timestep
after remeshing. Integration checks on Horizon (2026-10-03) exercised both
water-finer and solid-finer meshes with CPML, gravity and `face_mass`: all three
methods completed 8 s, matching removed the detected size mismatches, a matching
mesh imported successfully, and an incompatible mesh was rejected. Host and
CUDA-aware field-map histories agreed within $5.8\times10^{-9}$ relative L2
across displacement, velocity, acceleration and pressure. An 80 s replay of the
coarse sloping prototype through the new parameter reproduced its constraint
counts and pressure to $1.4\times10^{-9}$ relative L2. These are implementation
comparisons, not independent accuracy references. They predate the remote-solid
trace correction and did not establish MPI partition invariance. Investigation helpers and
their raw outputs remain local-only.

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
	spaced coordinates in metres along the material axes defined by
	`3D_velocity_model_azimuth_deg` (inherits `region_azimuth_leftface_deg`
	when omitted); the axis spacings may differ from one another.
- Let `nx` and `ny` be the counts in `x.in` and `y.in`. Horizontal profiles
	are stored with Y as the outer loop and X as the inner loop. For zero-based
	axis indices `i` and `j`, the profile number is
	`p = j * nx + i`, giving `(x0,y0), (x1,y0), ...` before advancing to
	`y1`. The count in `index.in` must be `nx * ny`.
- Material sample counts and profile start offsets use signed 64-bit integers.
	X/Y axis counts and the number of horizontal profiles remain limited to
	2,147,483,647. MPI initialization transfers large material arrays in bounded
	chunks; the count-prefixed text format is unchanged and old files remain valid.
	Build the updated reader before using models above 2,147,483,647 samples.
	Every solver rank holds the full model in host memory: approximately
	`32 * samples + 8 * profiles` bytes for Z, Vs, Vp, density and profile offsets,
	plus `8 * profiles` bytes for the SQD surface cache when needed, and another
	`16 * samples` bytes if Qs/Qp are supplied.
	For example, 3,003 by 6,003 columns with 121 levels require about 65.3 GiB
	per solver rank with the SQD cache and without Q files, before mesh/solver allocations. Plan the
	rank count per node accordingly; large-file support does not distribute the
	material arrays among ranks. Text parsing remains serial on solver rank zero.
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

### 3D material depths with topography

Without explicit seabed geometry, for `3D_velocity_model=yes` the selected vertical reference applies to all
aligned `z.in`, Vs, Vp, density and optional Q entries. With topography enabled,
`type_of_etree` is required (`sqd` or `full`, no default); it is independent of
the `computation_method=vt|fem` boundary treatment. Heights in `topography.in`
are positive upward; all material depths and model Z coordinates are metres,
positive downward. Write the files in the reference selected below:

| Mode | Meaning of `z.in` | Depth used at column $i$ |
| --- | --- | --- |
| `include_topography=yes`, `type_of_etree=sqd` | Layer-top depth below that column's own surface | $d_i=d+Z_s(q)-Z_s(i)=d+h_i-h_q$ |
| `include_topography=yes`, `type_of_etree=full` | Layer-top coordinate in fixed model Z | $Z_i=Z$; no terrain shift |
| `include_topography=no` | Existing fixed model Z reference | $Z_i=Z$ |

Here $q$ is the query location, $Z_s=B-h$ with
$B=\mathtt{topographybase\_zcoord}$, and the SQD caller supplies
$d=Z-Z_s(q)$. Existing building offsets are removed before the material query;
fixed model Z includes the existing mesh/database origin convention. It is
not automatically depth below sea level. For each column the SQD reader uses
the query's physical elevation:

```math
e_q=h_q-d,\qquad d_i=h_i-e_q=d+h_i-h_q.
```

Both $Z_s(q)$ and the cached column surfaces come from `point_elevation` on
the actual topography grid, in solver coordinates. The query surface is not
reconstructed by interpolating the four velocity-grid heights. This matters
for fine terrain on a coarse material grid and for translated/rotated grids.
The general 3-D convention samples at equal physical elevation; it does not
interpolate corresponding layers at equal burial depth. Terrain-following
interfaces can be encoded by column depths, but horizontal property blending
can smooth their discontinuities between columns.

For a concrete example, take $B=1000$ m and two columns with heights 0 and
100 m. A horizontal layer boundary at elevation $-40$ m has fixed
$Z=1040$ m. Encode that boundary as 40 and 140 m in SQD, or 1040 m in both
FULL columns. At the midpoint, a point at $Z=1030$ m has local burial depth
80 m. SQD samples the two profiles at 30 and 130 m; FULL samples both at
1030 m. Both select the upper layer. Switching the parameter alone without
converting the file depths changes the physical model.

Queries outside `xmin <= x < xmax` or `ymin <= y < ymax` retain the background
profile/CVM fallback, including the maximum grid edges. FULL and no-topography
queries below the smallest numeric `z.in` value also retain that fallback.
For SQD the raw query depth has a different origin from each profile, so no
global minimum-depth rejection is applied; each converted column query uses
the first/last material sample beyond that column's listed depth range.
Supply shallow coverage consistent with the intended model: this endpoint
extension is not an air-mask or a new material interface algorithm.

**Rebuild the solver and regenerate meshes with `mesh_input_mode=0` after
this correction.** The earlier general 3-D reader reversed the terrain-height
difference, used a surface interpolated on the velocity grid, and applied
column terrain shifts even in FULL mode. Imported element properties retain
those earlier values; new interface queries cannot repair an old mesh.
The setup monitor reports the selected 3-D depth convention.

The shared lookup serves refinement, element material assignment, nonlinear
geostatic density queries and acoustic interface reconstruction. Its chosen
layer supplies Vs, Vp, density and Q together. This affects elastic stiffness,
viscoelastic damping, Vs-indexed nonlinear elastoplastic/plastic and hysteretic
initialization, poroelastic material initialization, and acoustic/solid labels
used by acoustic–elastic, poroelastic–acoustic and poroelastic–elastic paths.
Water retains Vs=0 and solid Q weighting remains unchanged. Saturation mapping,
groundwater coordinates, yielding laws, interface constraints, time integration,
station interpolation, GPU kernels and both MPI exchange backends retain their
existing algorithms; their inputs can change after remeshing. Source and
receiver coordinate conventions are unchanged, but medium labels and output
waveforms can consequently change. CPML still clamps queries in solver
coordinates before the reader converts SQD depths at the clamped location.

Local analytical reader checks cover flat/sloping terrain, FULL/SQD equivalent
models, fine topography on coarse material grids, translated/rotated frames,
acoustic/solid transitions, Q selection, endpoint/fallback behavior and CPML
sampling. These ignored investigation fixtures are not shipped. They establish
material-coordinate correctness, not new wave-propagation validation: standalone
and coupled elastic/acoustic/poroelastic runs, applicable viscoelastic/nonlinear
cases and both GPU/MPI backends still require supercomputer comparisons.
Existing restrictions remain: field-specific coupled constraints require FEM
topography; VT uses the supported shared-geometric path. No new coupled
formulation or unsupported nonlinear/acoustic combination is enabled.

### Explicit seabed geometry

`seabed_geometry_file` enables geometry-defined water with the general 3-D
velocity model. Omit it to retain the legacy Vs-based convention. Negative
height alone never means water. Water depth and coverage come from the difference
between the outer surface in `topography.in` and the ground/bed surface in the
new file. Both contain elevations in the same datum, in metres positive upward.
A positive gap greater than `1e-8 m` is water; coincident surfaces are dry,
including dry land below sea level. There is no separate mask or water-level
input.

| Control | Requirement, default and meaning |
| --- | --- |
| `seabed_geometry_file` | Optional; omitted by default. Path relative to the solver working directory; grid format below. |
| `3D_velocity_model_depth_reference` | Default `legacy` without geometry. With geometry, explicitly choose `seabed`, `surface`, or `model_z`. |
| `seabed_water_vp_m_s` | Required with geometry, finite positive water Vp in m/s; no default. |
| `seabed_water_density_kg_m3` | Required with geometry, finite positive density in kg/m3; no default. Unlike `rho.in`, this is SI rather than tonne/m3. |

Explicit mode requires `3D_velocity_model=yes` and `include_acoustic=yes`.
`include_topography=no` is supported for a flat outer surface at solver Z=0. It rejects buildings and Istanbul, basin and
homogeneous-halfspace overrides, which otherwise bypass the shared reader.
Non-legacy depth references and water controls require the geometry file.
Use solid-only 3-D profiles: every Vs must exceed `1e-6 m/s`. Remove existing
water samples consistently from all aligned material arrays and rebuild
`index.in`. The solid dataset must contain its shallow solid endpoint at every
column; water is inserted from geometry before interpolation. Water uses
`Vs=0`, the supplied Vp/density and existing `acoustic_quality_factor` controls.
Solid Q interpolation, viscoelastic damping and material calibration retain
their existing roles.

With topography enabled, the geometry shares the origin, extent and spacing
of `topography.in`, in
solver X north/Y east axes, regardless of material-grid rotation. Let
`Nnorth=region_length_north_m/spacing+1` and similarly for `Neast`. Its text
format is:

1. One positive finite spacing dividing both horizontal domain lengths,
   matching `topography.in` when topography is enabled.
2. `Nnorth*Neast` finite ground/seabed elevations in metres, positive upward,
   ordered north outer/east inner (east fastest).

There are no count prefixes, coordinates, masks, comments or extra tokens.
One number per line is supported, exactly as in `topography.in`; spaces and
newlines are equivalent. The two elevation fields use bilinear interpolation
and nearest-edge extension outside the domain. At each query, water occupies
the gap between the surfaces when that gap exceeds the fixed `1e-8 m`
tolerance. Bed elevations above topography by more than that tolerance abort
setup; nodal differences within tolerance are snapped to the outer surface.
Dry queries also snap their bed reference to the outer surface.

With `include_topography=yes`, `topography.in` remains the **outer air/material surface**. It supplies the
local water surface wherever the bed is lower, and dry ground where both
surfaces coincide. A wet grid cell (any corner has positive water depth) must
have a flat outer surface below `topographybase_zcoord`: existing acoustic
free-surface operators require horizontal mesh faces. Separate flat water
bodies may have different elevations, with dry cells between them. Sloping
water surfaces are rejected at setup. Resolve shorelines on this grid; do not
put submerged bed elevations into the outer topographic surface.

For a 1000 m north by 1000 m east domain with 500 m spacing, this complete
`input/seabed.in` defines a flat bed 100 m below sea level:

```text
500
-100
-100
-100
-100
-100
-100
-100
-100
-100
```

The corresponding `topography.in` has the same spacing followed by nine zero
heights. Add this excerpt to the normal case configuration (illustrative water
properties, not a site calibration):

```ini
include_topography = yes
3D_velocity_model = yes
include_acoustic = yes
seabed_geometry_file = input/seabed.in
3D_velocity_model_depth_reference = seabed
seabed_water_vp_m_s = 1500
seabed_water_density_kg_m3 = 1025
```

For a flat ocean, use the same seabed file and water/solid controls but set:

```ini
include_topography = no
```

No `topography.in`, `topography_directory`, `topographybase_zcoord`,
`maximum_octant_level`, `type_of_etree` or `computation_method` is needed for
this mode; the topography-only controls are unused. The outer surface and
height datum are both zero: `bed_z = -bed_elevation`. Negative bed elevations
are underwater, zero is dry ground, and positive elevations beyond the fixed
tolerance are rejected. The seabed itself may slope. The uniform spacing must
divide both horizontal domain lengths into positive integer cell counts;
node counts and the exact number of file values are checked before use.
CVM mesh offsets must be zero, as normally obtained with
`cvmdb_input_file=profile`; translating a CVM subvolume is not supported in
this flat mode. The general 3-D dataset still supplies all solid properties.

The topography surface-refinement band and air/cut-cell operators are skipped.
Wavelength refinement, optional acoustic-interface refinement, local matching
when selected, and 2:1 balancing remain. In particular, `maximum_octant_level`
is not a global refinement cap with topography off. Adequate water-column,
seabed and gravity-wave resolution must still be established by convergence
checks. This mode provides no moving or sloping water surface.

With topography off, `surface` and `model_z` both mean depth below Z=0;
`seabed` means depth below the local bed. In the 100 m water example, a solid
layer 20 m below bed uses `z.in=120`, `120`, or `20`, respectively. Sources,
stations and planes retain the flat solver depth convention; a seabed
receiver is at depth 100 m, not zero.

The explicit reference is independent of `type_of_etree`. Both `full` and
`sqd` work: the caller convention is converted back to physical solver Z
before evaluating the selected column reference. For column $i$:

```math
z_i =
\begin{cases}
Z, & \texttt{model\_z},\\
Z-Z_s(i), & \texttt{surface},\\
Z-Z_b(i), & \texttt{seabed}.
\end{cases}
```

$Z_s=B-h_s$ is the outer surface, $Z_b=B-h_b$ the ground/bed and
$B=\mathtt{topographybase\_zcoord}$ with topography enabled (zero in flat mode). All `z.in` coordinates are metres,
positive downward. For $B=1000$, bed elevation $-100$, water level 0 and a
layer 20 m below bed, write 1120 (`model_z`), 120 (`surface`), or 20 (`seabed`).
On dry cells seabed-relative means below ground. At the query point, water
occupies $Z_s\le Z<Z_b$ where $Z_b-Z_s>10^{-8}\,\mathrm{m}$; the bed belongs to solid.
The actual fine geometry is used even with coarse/rotated velocity columns.

Solid columns retain piecewise-constant vertical lookup and bilinear horizontal
interpolation at equal physical elevation. First/last solid endpoints extend
beyond listed depths. A solid query outside the half-open velocity footprint
aborts rather than using an ambiguously referenced background model. Extend
velocity axes beyond maximum query coordinates. CPML clamps queries in solver
coordinates before geometry and material lookup; the bathymetry reconstruction
uses the same lateral clamp and reads bed depths directly.

Physics and runtime integration:

- Elastic, viscoelastic and nonlinear elastoplastic/plastic and hysteretic
  solids receive the selected solid properties; water is excluded from solid
  interpolation and retains the acoustic marker. Material laws are unchanged.
  Calibrated geostatic, tolerance, substep, damping and
  `material_properties_list` values remain case-specific.
- Porous depth limits are measured below ground/bed for both FULL and SQD.
  Beneath explicit water, the saturation/hydrostatic reference is the water
  surface from local topography, or Z=0 with topography off;
  `groundwaterlevel_file` is ignored there.
  Dry-cell groundwater depths remain relative to the ground. Keep the pore-fluid and water densities
  physically consistent. Nonlinear approximate geostatic initialization adds
  the exact water-column load and uses 100 midpoint samples for the solid
  overburden before the existing hydrostatic Biot-stress correction.
- Acoustic free-surface selection uses Z=0 with topography off, or local
  topography at each wet element with it enabled,
  including separate flat lakes/reservoirs at different elevations. Only
  complete horizontal mesh faces coincident with that surface are enrolled;
  crossing or submerged faces are excluded. Sloping/changing water surfaces
  and cut acoustic free surfaces remain unsupported. Existing gravity and
  pressure-release formulations retain their mass/force algorithms.
- The optional `acoustic_interface_face_filter=bathymetry` reads bed depths
  from geometry while retaining its existing Z+ bottom-face-only policy.
  Default all-face coupling and the acoustic–elastic, poroelastic–acoustic and
  poroelastic–elastic formulations retain existing constraints and algorithms.
  FEM/VT restrictions remain. Existing 27-point element classification gives
  mixed water/solid cells solid priority; resolve thin water layers and check
  mesh convergence. This is not a conforming seabed mesh.
- Rank zero reads the geometry and broadcasts identical arrays to solver
  ranks. All material/geometry work stays in initialization; time integration,
  device-resident GPU fields and both MPI exchange backends receive the existing
  element metadata without new per-step messages or transfers. Station/source
  depth conventions continue to reference the outer surface (Z=0 in flat mode): a station on the
  example seabed has surface-relative depth 100 m. Medium labels, interface maps
  and output waveforms can change after remeshing.

**Rebuild and regenerate with `mesh_input_mode=0` after changing these inputs.**
Importing a mesh does not update saved element material assignments. Geometry
is still loaded for interface/poroelastic setup on imported runs; it must match
the mesh. The setup monitor reports the explicit reference, water properties
and geometry tolerance.

A 2026-10-18 Horizon benchmark at solver commit `8c9a514` compared this
flat-ocean mode against independent CPS layered-medium Green functions. Two
nodes/eight GPUs used CUDA-aware MPI, velocity-potential acoustics and
field-specific constraints for 110 s, with a 937.5 m water layer over a uniform
elastic solid. The same physical 1 Hz source, 0.0025 s step and CPML controls
were used on 22,282,240- and 178,257,920-element meshes. Refined off-axis median raw
waveform L2 errors over 10–100 s were 2.89% displacement,
6.64% velocity and 9.80% pressure. No response filtering,
fitted gain or timing shift was applied. Some refined receiver peaks exceed the 5% target; the two resolutions do not establish mesh independence.
All 30 receivers and the generated source history were checked. This benchmark
covers a horizontal seabed, pressure-release surface and lossless materials;
it does not cover gravity waves, sloping-bed propagation or the host-MPI backend.

Local extracted C++ checks cover flat-ocean parser/geometry/material behavior,
all three depth references, sloping-bed normals, refinement dispatch and surface
mass assembly. Additional host tests exercise production porous pressure and
traction algebra, sealed-interface content conservation, open-pore quadrature,
local/remote constraint selection, station eligibility, solid Q conversion,
water exclusion from solid damping/nonlinearity, and prescribed hydrostatic and
approximate geostatic effective stress. CUDA scheduling and MPI transport are
mocked in these tests. Nonlinear skeletons retain explicit pore-pressure coupling.
A zero-permeability gravity-loading ramp is rejected by the existing solver.
Submerged ramped equilibrium and viscoelastic, poroelastic and nonlinear
elastoplastic/plastic or hysteretic standalone/coupled propagation still require
separate validation; an analytical initial-stress check is not a gravity-ramp
equilibration test. The ignored investigation harness and benchmark results are
not distributed in a clean checkout.

### Quality-factor files

- With `3D_velocity_model=yes`, the material directory may contain optional
	`qs.in` and `qp.in` files alongside `vs.in`, `vp.in`, and `rho.in`.
	Each file uses the same count-prefixed layout, contains one value for every
	entry in `z.in`, and follows the same `index.in` column/layer ordering.
	Hercules selects the same vertical layer at each of the four surrounding
	horizontal columns. For `Qs` and `Qp`, it excludes acoustic samples
	(`Vs <= 1e-6 m/s`) and renormalizes the remaining bilinear weights. An
	entirely solid stencil retains the previous interpolation arithmetic;
	an all-water query returns zero solid Q. Velocity, density, topography,
	and acoustic/solid classification are unchanged. With parametric Q disabled,
	the BKT-family damping setup uses these interpolated values independently
	of the fallback `cvmdb_input_file` mode.
- `qs.in` and `qp.in` must either both be present or both be absent, and their
	value counts must match `z.in`. Each solid sample must have finite,
	positive Qs and Qp. Acoustic placeholders may be zero but must be finite
	and nonnegative. These checks run as the files are read, even if solid
	attenuation is disabled. If both files are absent, Hercules preserves the
	legacy behavior `Qs = 0.1*Vs` and `Qp = 2*Qs`. If only one is present or a
	count differs, initialization stops with an explicit error.
- Raising a queried solid Vs to `simulation_shear_velocity_min` does **not**
	rescale explicitly supplied Qs or Qp. For example, `Vs=80 m/s`, `Qs=16`,
	and `Qp=32` retain Qs=16 and Qp=32 if Vs is raised to 120 m/s. Update the
	Q model explicitly if a velocity-dependent Q law is intended for the
	modified velocities. The optional parametric-Q path remains separate.

For either supplied Q, the mixed-stencil rule is

```math
Q_{\mathrm{solid}} = \frac{\sum_{i\in S} N_i Q_i}{\sum_{i\in S} N_i},
\qquad S=\{i:V_{s,i}>10^{-6}\ \mathrm{m/s},\ N_i>0\}.
```

For example, 10% weight on solid Qs=40 and 90% on water Qs=0 gives
solid Qs=40, not 4. This prevents water placeholders from introducing
artificially strong solid attenuation; it does not change the velocity mixture.

#### Setup-time attenuation validation

During the existing material-assignment loop, BKT-family setup requires finite,
positive solid Qs/Qp, velocities, density, and reference frequency, with a
positive bulk modulus. Unless `use_infinite_qk=yes`, it derives

```math
L=\frac{4V_s^2}{3V_p^2},\qquad
Q_\kappa=\frac{1-L}{1/Q_p-L/Q_s}.
```

With `use_infinite_qk=no`, a negative denominator beyond floating-point
roundoff now selects **Qk=1000 for that element only**. This is Hercules's
existing high-Q, lossless-bulk limit: all bulk relaxation weights, rates and
viscous coefficients are zero. It applies to `bkt`, `bkt2`, `bkt3`, and `bkt3f`.
Valid positive-bulk calculations and their existing high-Q cutoffs are
unchanged, as are Qs and the shear coefficients. Zero inverse bulk Q, including
cancellation within roundoff, continues to select lossless bulk. In contrast,
`use_infinite_qk=yes` suppresses bulk attenuation for every solid element;
both settings still require finite, positive solid Q inputs.

The setup diagnostic reports the global `negative-Qk bulk fallbacks` count.
That count excludes the roundoff/lossless limit and the explicit
`use_infinite_qk=yes` override. A nonzero count also prints the fallback value
and states that the original Qp is no longer reproduced in those elements.
No new input is needed: keep `use_infinite_qk=no` to retain valid bulk
attenuation elsewhere. This follows the high-Q fallback approach used in
[SPECFEM3D's tomography Q conversion](https://github.com/SPECFEM/specfem3d/blob/master/src/generate_databases/model_tomography.f90),
using Hercules's own lossless cutoff rather than SPECFEM's numerical constant.

For example, illustrative values `Vp=1500 m/s`, `Vs=1000 m/s`, `Qs=100`,
`Qp=200` imply a negative Qk and therefore select the fallback. With
`Vp=2000 m/s` and the same Vs/Qs/Qp, the valid derived `Qk=400` is retained.
These are examples, not calibrated material values. Setting lossless bulk
changes the effective P-wave attenuation, approximately

```math
\frac{1}{Q_{p,\mathrm{effective}}}=\frac{L}{Q_s},
\qquad Q_{p,\mathrm{effective}}=\frac{Q_s}{L}.
```

P waves still attenuate through the shear contribution. The usual
`simulation_velocity_profile_freq_hz` correction still runs with the selected
coefficients: shear dispersion is unchanged, and the fallback removes bulk
dispersion only in affected elements. It does not introduce a Vp/Vs clamp or
make an inconsistent supplied Qp/Qs pair physically equivalent to the original.

BKT2 additionally requires finite, nonnegative stored coefficients and positive
relaxed stiffness, `1-a0-a1 > 0`, in both shear and bulk branches. Active
relaxation mechanisms must have positive rates. Invalid values stop setup with
a global invalid-element count and one example giving rank, element center,
material values and Q values. The negative-Qk fallback does not bypass these
checks or repair nonfinite/nonpositive Q, invalid material properties, or
inadmissible BKT2 shear coefficients. The timestep is not changed by this
policy. These are admissibility checks, not a guarantee of stability for an
arbitrary timestep or model.

The checks reuse values already read or computed: no additional material
queries, mesh refinement, full mesh traversal, or time-step checks. Successful
material validation adds one small MPI reduction per material-assignment call.
Imported BKT2 meshes have their stored coefficients checked during the existing
critical-timestep traversal. They do not retain the original Q inputs, so that
diagnostic reports unavailable Q values as NaN. Regenerate a persisted mesh to
apply this fallback or corrected interpolation, or after changing Q files;
a coefficient check cannot recover the original Q or repair an older,
admissible but incorrectly interpolated Q model.

The fallback is assigned in the common solid material setup, including solids
in acoustic–elastic, poroelastic–elastic and poroelastic–acoustic cases. Pure
acoustic elements bypass BKT setup and remain lossless; none, Rayleigh and
mass damping are unaffected. Existing selection of BKT-active cubes and VT
cut cells, elastic CPML treatment, and exclusion of nonlinear elements from
the BKT force path remain unchanged. The policy does not alter poroelastic
pressure evolution or nonlinear elastoplastic/plastic and hysteretic laws.
CPU/GPU force and station-stress paths consume the same stored coefficients;
there are no new boundary/interface terms, GPU kernels or MPI field exchanges.
Host coefficient and phase-velocity comparisons cover all four BKT variants;
they do not validate full GPU/MPI wave propagation or any of the coupled or
nonlinear cases. Those still require platform benchmarks with the rebuilt
solver and regenerated mesh.

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
	pressure-traction feedback. Pore pressure remains available in station output as a diagnostic.
	`HERCULES_PORE_PRESSURE_UPDATE=output` (default) reconstructs it at
	required output states; `step` retains incremental `Q^T du` updates.
	Checkpoint/restart runs always retain every-step updates (see
	[pressure update modes](#scalar-mpi-schedules-and-profiling)). Set this to `no`
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
	On regular Hex8 meshes (including FEM topography), pore-pressure hanging-node
	constraints are built automatically from final active BrickUP elements,
	independent of acoustic enablement and `acoustic_interface_method`. A node
	constrained only by a coarser dry or elastic-only CPML cell remains independent
	for pore pressure. A coarser active BrickUP cell retains its exact Hex8 trace
	constraint; unsupported stencils abort setup. This applies to capacity,
	drained/open-pore boundary weights, pressure RHS and field assignment on both
	host and CUDA-aware MPI paths. Solid displacement keeps its selected map.
	The startup monitor reports pore-pressure keep/release/absent counts.
	VT cut-cell topography retains geometric pressure constraints and reports
	that limitation; the Hex8 change does not extend to cut-cell pressure traces.
	No new input switch is required. For example, `include_acoustic=no` with
	`include_poroelastic=yes` and `poroelastic_elastic_interface=no_flow` enables
	the independent pore map on a regular Hex8 mesh.
	The standalone extension requires new CUDA/MPI propagation and refinement
	comparisons on the target system; earlier acoustic-module benchmarks do not
	validate this newly enabled path.
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

- **Gravity-surface CPML nodal consistency.** The full second-order operator
	now evaluates the tangential product `S_x*S_y` separately at each of the
	four lumped top nodes, including its corner term. Earlier face-average
	memories omitted zero-mean face modes. Historical source-corrupted tests
	do not validate this formulation. Gravity-wave absorption and motion
	reconstruction require separate checks, since surface-gradient recovery
	error remains even with a consistent absorber.
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
	axis, the trapezoidal inverse-stretch section is applied to `grad(p[n])`;
	outside CPML this is an exact
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
second-order weak form with complex frequency shift and full Jacobian products). The C-PML is per-element, wraps the five
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
	or unset (default), the solver uses `cpml_layers ×` the largest non-air
	element edge touching any of the four lateral faces or the bottom face,
	taking an MPI-wide maximum. Air elements (`Vp < 0`) are excluded, including
	from the global non-air fallback used if no material element touches an
	absorbing face. Material extension and optional CPML partition weighting
	use the same exclusion rule. The realized thickness is logged to standard
	error at startup with `air excluded`; PML element counts and the optional
	outer-clamp count use the same startup stream.
- `cpml_layers`
	Optional integer. Default `10`. Used only when `cpml_thickness_m` is
	nonpositive. This multiplies one reference edge length to give a uniform
	physical thickness; finer parts of an adaptive boundary contain more than
	`cpml_layers` cells. For example, if a 90 km model has level-8 boundary air
	cells (351.5625 m) but its largest non-air boundary edge is level 9
	(175.78125 m), `cpml_layers=10` gives 1757.8125 m, excluding the air cells.
	This example does not establish the active-material maximum of a particular
	mesh; check its startup log. Automatic thickness can decrease after this
	fix, changing absorption and material-extension bounds. A positive
	`cpml_thickness_m` retains a specified physical width.
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
	normal in absolute model coordinates. For SQD topography, restore
	model Z from depth below the original local surface before clamping
	at the bottom, then convert back using the surface at the clamped
	X/Y location. Building surface offsets are also restored and removed.
	For example, a surface at model Z = 2000 m and a bottom-PML entrance
	at Z = 9000 m map a point at Z = 9500 m (relative depth 7500 m)
	to Z = 9000 m (relative depth 7000 m). At bottom/lateral corners,
	the final relative depth uses the clamped column's surface.
	This changes sampled element material properties, so rebuild
	or re-correct an imported etree after changing this setting or
	adopting the corrected topographic bottom/corner conversion;
	Phase B imports that reuse baked properties keep those existing
	values. Set `no` to restore the legacy raw 3D material queries.
- `cpml_reflection_target`
	Optional. Target theoretical reflection coefficient. Default `0.1`; valid `(0,1)`.
	Smaller values increase damping strength but do not guarantee smaller
	numerical reflection or better stability.
- `cpml_poly_order`
	Optional integer in `[1, 4]`. Default `2` (the standard
	K&M/SPECFEM3D choice). Controls how steeply the damping ramps up
	across the PML.
- `cpml_kappa_max`
	Optional. Default `10.0`; must be at least one. Higher values change the
	real stretch in all stiffness, inertia and surface products. Validate
	thickness, resolution and incidence angle together; increasing kappa is
	not a substitute for a resolved layer.
- `cpml_alpha_max_hz`
	Optional. Default: `simulation_wave_max_freq_hz/2` in Hz. Omitted or
	nonpositive inputs retain automatic selection (internal sentinel `-1`).
	A positive numeric value overrides the automatic frequency; arithmetic
	expressions are not evaluated in the input file. For `simulation_wave_max_freq_hz=4.36`,
	the automatic value is `2.18` Hz. The angular coefficient is
	`alpha_max=2*pi*cpml_alpha_max_hz`, giving `pi*simulation_wave_max_freq_hz`
	for automatic selection. This automatic behavior is unchanged; see the
	[coefficient formulas](#cpml-coefficient-and-operator-details) before
	matching a source reference frequency.
- `cpml_cfl_safety`
	Optional. Default `0.9`. The critical time step is multiplied by this
	factor when `abc_type=cpml`, for margin against the corner-element
	instability. `1.0` disables the reduction; `< 0.5` is rarely needed.
- `cpml_mpml_cross_fraction`
	Optional. Default **`0.0`**, changed from `0.1` with the complete operator.
	Valid range `[0,0.5]`. Nonzero values mix each damping and kappa profile
	into the other axes. This experimental stabilization breaks the separable
	coordinate mapping and can increase reflections substantially. Existing
	inputs explicitly setting `0.1` retain that behavior; remove the setting
	or write `0` for plain CPML.
- `cpml_acoustic_sponge_fraction`
	Optional. Default `0.3`; finite, dimensionless and nonnegative.
	Set `0` to disable this experimental acoustic bulk damping, applied only in CPML:

```math
a_{\mathrm{ac}}=\gamma_{\mathrm{ac}}(d_x+d_y+d_z),\qquad
\gamma_{\mathrm{ac}}=\texttt{cpml\_acoustic\_sponge\_fraction}.
```

	The acoustic volume inertia becomes `$J(s)(s^2+a_{\mathrm{ac}}s)$`;
	`$a_{\mathrm{ac}}$` has units of inverse seconds. This applies to pressure
	and velocity-potential acoustics. It adds no static force and does not
	change `acoustic_quality_factor`, which remains a separate attenuation
	control. Elastic and interior coefficients are unaffected by this option.
	The gravity-surface inertia retains its ordinary tangential stretching,
	with no additional sponge term; interface cofactors are also unchanged.
	Coupled waveforms can nevertheless change.

	This is independent of `cpml_mpml_cross_fraction`. Setup reuses existing inertia coefficients
	and GPU histories, without additional kernels, exchanges or per-step
	branches. Keep the default `0` for baseline validation. Nonzero values
	alter matching and require raw-waveform accuracy and long-time stability
	checks; they are not a universal remedy for heterogeneous CPML growth.
- `cpml_alpha_floor_fraction`
	Optional. Default `0.1`; valid `[0,1]`. Lower bound on alpha relative to
	its maximum. Zero and repeated poles are handled without divided-pole
	singularities; their acceptance does not establish long-time stability
	for every model.
- `cpml_pml_material_damping`
	Optional. `none` (default) | `rayleigh`. Whether material damping is
	applied inside elastic PML elements. `none` means material damping is
	skipped in the PML zone; non-PML elements still use the requested
	`type_of_damping`. For BKT-family runs (`bkt`, `bkt2`, `bkt3`,
	`bkt3f`), this is the supported policy: the interior uses BKT damping,
	while PML elements use the complete elastic CPML operator.

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
	* `native` — fully integrated trilinear Hex8 stiffness, including all
		tangential and mixed terms. Its unstretched limit equals the existing
		interior operator. Recommended.
	* `reduced` — one-point integration, retaining only linear modes.
		Rank deficient; kept as a diagnostic option.
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

### Topographic surface geometry and sampling

The topography grid defines a bilinear height field shared by mesh air/solid
classification, material-depth queries, VT volume sampling, sources, stations,
and planes. For file height $h(X,Y)$ and station depth $d$:

```math
Z_{\mathrm{surface}}(X,Y) = \mathtt{topographybase\_zcoord} - h(X,Y),
\qquad Z_{\mathrm{station}} = Z_{\mathrm{surface}} + d.
```

Solver Z is positive downward. `region_depth_deep_m` sets the domain extent;
it does not replace the surface datum. FEM retains whole Hex8 staircase
elements, so this consistency correction does not remove terrain-discretization
error or provide a fitted surface mesh.

Earlier builds used two triangles per grid cell in the mesh air/solid test,
but bilinear heights for material and receiver queries. On nonplanar terrain
those surfaces differ: a surface receiver could lie several fine cells above
the mesh's solid boundary. Its trace could be zero, or nonzero through a few
shared active nodes without being a valid material interpolation. Refinement
can expose this mismatch. The discrepancy predates the shifted acoustic
free-surface correction; changing the datum does not fix it.

**Regenerate affected topographic meshes with `mesh_input_mode=0` and rerun.**
Importing an older mesh preserves its stored air/solid labels, even with a
rebuilt executable. Existing output cannot recover samples missing for this
reason. Verify positive `Vp_mps` and `rho_kgm3` in the new `station_mediums.txt`
for land stations. Use `HERCULES_PLANE_METADATA=1` to generate fresh plane
ownership metadata; an older run's metadata is not evidence of current owners.

The final grid row/column uses its actual edge heights. Out-of-grid material
queries clamp to the nearest edge. For meshing proximity, `point_PlaneDist`
returns $(Z_{\mathrm{surface}}-Z)/\sqrt{1+h_X^2+h_Y^2}$: positive in air,
negative in material. This is a local tangent-normal offset, not an exact
shortest distance to a curved bilinear patch. It preserves the planar normal
distance and the existing refinement thresholds while aligning their zero
surface with all other topography queries. The startup message identifies this
shared bilinear geometry. No time-stepping exchange or GPU kernel is added.

Local regression on 2026-10-06 replayed the production classifier at four
air-owned Izmir stations and 4,806 zero-valued interior plane locations; all
containing fine cells became material with the shared surface. A separate
1 km terrain crop with homogeneous elastic material and a synthetic source
reproduced and removed the zero station trace and 431 missing histories in a
1,521-point plane. One/eight GPU ranks and host/CUDA-aware MPI agreed within
`2.3e-13` relative L2 in raw plane histories. An 8 s flat-water control with
gravity, CPML and field-specific constraints agreed with the previous build
within `3.7e-12` relative L2 in plane output. The local-only investigation
fixtures are not shipped. These checks establish geometry/output consistency,
not continuum accuracy or full-Izmir long-time stability.

### CPML execution and topography

- Each PML element is automatically classified as elastic or acoustic
	(the same test the acoustic-pressure path uses), so heterogeneous
	models and acoustic/elastic coupling need no extra configuration.
	For `3D_velocity_model=yes`, `cpml_3d_material_extension=yes`
	(default) makes that classification come from the corresponding
	inner-PML material column/row, matching the usual flat-PML material
	extension used by spectral-element PML setups. For SQD topography,
	query depth is measured below the local surface, while CPML bounds
	use mesh coordinates. The extension clamps the physical query point
	in mesh coordinates, including the bottom limit, then converts back
	to depth below the surface at the extended horizontal location.
	This handles a nonzero topography datum consistently; it does not
	by itself establish stability for a heterogeneous CPML model.
- Material damping on PML elements follows `cpml_pml_material_damping`
	(see that flag). With the default `none`, material damping is skipped
	on PML elements, but non-PML elements always use material damping as
	usual. This is the supported policy for BKT-family damping: the
	interior is viscoelastic, while the PML zone is elastic CPML.
	`rayleigh` is the only supported material-damping mode inside PML
	elements.
- **Viscoelastic damping runs on the GPU for all BKT-family models**
	(`bkt`, `bkt2`, `bkt3`, `bkt3f`). BKT2 cube elements with matching
	uniform shear/bulk relaxation rates use the shared-history fused path
	described in [BKT2 GPU execution](#bkt2-gpu-execution). For the separate
	convolution/force paths, when the relaxation spectrum is
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
	all elements (a latent bug). `HERCULES_BKT2_FUSED=0` retains the
	previous hoisted BKT2 path for comparisons with the fused implementation.
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

### CPML setup and validation requirements

Use `stiffness_calculation_method=effective` for the GPU elastic path.
Checkpoint reads and writes are rejected with CPML: the legacy checkpoint
format does not preserve absorber histories. Set `use_checkpoint=0` and
`checkpointing_rate=0`. Non-PML checkpoint behavior is unchanged.

The defaults select `cpml_stiffness_mode=native`, `cpml_kappa_max=10`,
`cpml_reflection_target=0.1`, `cpml_acoustic_sponge_fraction=0.3`,
automatic alpha in Hz equal to `simulation_wave_max_freq_hz/2`,
`cpml_mpml_cross_fraction=0`, and alpha floor `0.1`. These candidate settings
are not a universal stability or accuracy guarantee; slow gravity packets can
require a lower alpha. Explicit input values always take precedence. To reproduce
the previous defaults, explicitly set `acoustic_interface_method=shared-geometric-constraints`,
`cpml_kappa_max=1`, `cpml_reflection_target=0.001`, and
`cpml_acoustic_sponge_fraction=0`; automatic alpha is unchanged. Keep source and
receivers outside the absorbing layer, resolve the source bandwidth, and
compare with a larger domain at identical mesh spacing and timestep. Measure
P/S arrivals, surface waves, late coda, and the permanent displacement
separately; neither finite output nor a quieter trace proves matching.
The profile is sampled at element centers, so spatial refinement and a thicker
layer can reduce the entrance discretization error.

Historical CPS results below describe earlier incomplete operators and are
not acceptance evidence for the full formulation. Current implementation and
validation scope are described in [operator details](#cpml-coefficient-and-operator-details).

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
frame (local solver X/Y, Z down; UTM north/east at zero azimuth). Internally the FORCE
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

The interpolator holds the final sample at and after
`tinit + (nt1 - 1) * dt`. Its interval index is also bounded when
floating-point division rounds a time just below that endpoint to the last
sample index. For example, with 141 samples at 0.02 s, `2.8 / 0.02` can
round to 140 while `2.8 < 140 * 0.02`; reading sample 141 would be outside
the supplied history. Older builds could therefore generate an enormous
force or NaN from otherwise valid discrete input. Rebuild and regenerate
source-force files before repeating an affected simulation.

SRFH force histories are generated at the common `slipfunction.in` interval
and then interpolated at the solver timestep. Reducing
`simulation_delta_time_sec` does not change that source-generation interval.

Point and single-point body-force source generation writes all loaded nodes,
including when a long time history spans multiple memory batches. The binary
force-file header contains the complete node list, and each time sample uses
that same ordering. Increasing the simulation duration must preserve the
previously generated early-time forces. This batching is part of source setup;
it does not change the GPU time-stepping path.

### Source-force filtering

In `source.in`, `source_is_filtered` is required (`0` or `1`). When it is `1`,
`threshold_frequency` (cutoff in Hz) and `number_of_poles` (integer Butterworth
order) are required; both are ignored when filtering is disabled. For example:

```ini
source_is_filtered = 1
threshold_frequency = 1.0
number_of_poles = 14
```

Filtering is applied by the `plane`, `planewithkinks`, and `srfh` source-generation
paths. It uses the interval of the generated force samples: the common
`slipfunction.in` interval for SRFH, and the solver timestep for `plane` and
`planewithkinks`. An SRFH source sampled every 0.01 s uses
a 100 Hz sampling frequency even if `simulation_delta_time_sec = 0.005`.
Changing only the solver timestep must not change its physical filter cutoff.
Older builds instead used the solver interval for SRFH filtering, making the
effective cutoff `threshold_frequency * solver_dt / source_dt`. Rebuild and
regenerate source-force files after the correction; remove any compensation
for that old ratio and specify the desired physical cutoff directly.

Filtering differentiates the force, applies a zero-phase FFT Butterworth
magnitude, then integrates from time zero. The filtered response extends
before the nominal source start. If the source starts too early, truncating
this prehistory can leave a startup transient and change the final applied
force. Leave sufficient pre-source time and check generated force histories
and final force while increasing the buffer. The required duration depends
on cutoff, order, source history and the tolerated residual; it is not a
universal number of seconds. For SRFH, set the second value of each
`slipfunction.in` source block to its start time in seconds. A delay does not
make this zero-phase operation strictly causal relative to that start time.
Filtering changes the physical source spectrum, so choose a cutoff consistent
with both the desired bandwidth and mesh resolution.

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
base_force_direction_x   = 1.0      # solver +X (UTM north at zero domain azimuth)
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

The production helpers are in `quake/forward/cpml.h`. For Laplace variable
$s$, each axis has stretch $S_j$, and $J$ is their product:

```math
S_j(s)=\kappa_j+\frac{d_j}{s+\alpha_j},\qquad J(s)=S_xS_yS_z.
```

The elastic weak stiffness multiplies $C_{ijkl}$ by $J/(S_jS_l)$; the
scalar acoustic gradient term uses $J/S_j^2$. Thus an x-only PML includes
both inverse x stiffness and forward tangential stiffness. The full inertia
is $\rho J(s)s^2$. With Rayleigh material damping it becomes
$\rho J(s)(s^2+a s)$, while stiffness acts on the usual effective displacement
$u+b\dot u$. Edge/corner products are retained in every term.

For normalized depth $r_j=\xi_j/L$, reference speed $V_p$, target $R$ and
polynomial order $n$, the base profiles are:

```math
d_j=-\frac{(n+1)V_p\ln R}{2L}r_j^n,\qquad
\kappa_j=1+(\kappa_{\max}-1)r_j^n,\qquad
\alpha_j=\max\{\alpha_{\max}(1-r_j),f_\alpha\alpha_{\max}\}.
```

Inactive axes have $d=0$, $\kappa=1$, $\alpha=0$. Positive
`cpml_alpha_max_hz` gives $\alpha_{max}=2\pi f$; otherwise the fallback is
$\pi$ times `simulation_wave_max_freq_hz`. Optional cross mixing adds a
fraction of the other two damping profiles and normalized kappa profiles;
it defaults to zero because this is no longer a separable coordinate mapping.

A section $c+g/(s+p)$ uses a trapezoidal recurrence with one persistent value:

```math
\beta=\frac{\Delta t}{2+p\Delta t},\qquad
b=\frac{2-p\Delta t}{2+p\Delta t},\qquad
y_n=(c+g\beta)x_n+m_n,\qquad
m_{n+1}=b m_n+g\beta(1+b)x_n.
```

Products are cascades of these sections. Polynomial division splits
$J(s)(s^2+a s)$ into leading mass, velocity and restoring terms plus a
proper-rational remainder. The remainder uses nested denominators, avoiding
division by differences of poles. Inactive factors cancel before division;
equal alpha values and zero poles are well defined. The leading velocity
term retains Hercules' backward-difference damping convention; trapezoidal
ADE sections alone do not make the entire time integrator second order.

Native integration uses orthogonal Hex8 sign modes, algebraically equivalent
to eight-point Gauss integration for element-constant profiles. All modes
are stretched, and nodal inertia histories use the same lumping as the leading
mass. Interior kernels skip PML elements to avoid double counting. The
unmodified native normalization is $\mu\Delta t^2h=9c_1$ and
$\lambda\Delta t^2h=9c_2$. Reduced integration retains only linear modes.

An interface flux component normal to axis $j$ uses the cofactor $J/S_j$.
Both coupling directions use the same coefficients evaluated at the physical
interface-fragment center, independently of the adjacent solid element size
or MPI owner. Sampling each containing element separately can give different
cofactors for the same nonmatching fragment. Histories remain independent for
each direction, and all normal components are retained, including bathymetry
normals. `face_mass` retains four quadrature points in the PML, each with its
own histories. This corrects coefficient consistency across the interface;
it does not establish long-time stability for arbitrary heterogeneous CPML. A flat gravity
surface uses $S_xS_y s^2/(\rho g)$ at each of its four lumped nodes, including
the x/y corner product. Output reconstruction uses the same inverse stretch,
with a separate history per observation point.

Coefficients are prepared once on the host. Histories remain on the GPU and
feed the existing force assembly and CUDA-aware or fallback MPI exchange;
there is no new per-step exchange or bulk host transfer. Native bulk histories
reserve 204 `solver_float` values per elastic PML element and 60 per acoustic
PML element; gravity faces reserve eight, and each coupled face reserves 24
per side. This increases absorber memory relative to the earlier incomplete
operator. Kernel timing and total GPU memory must be checked for a target mesh.

### Complete CPML validation on Horizon (2026-10-02)

The development runs used double precision, CUDA 13.2/NVIDIA 26.5,
Open MPI, and native Hex8 integration. The prepared homogeneous case has
$V_p=2000$ m/s, $V_s=1000$ m/s, $\rho=2000$ kg/m³, $h=78.125$ m,
$\Delta t=0.005$ s, and a 781.25 m PML. The 20 km reference preserves mesh
spacing, source, and physical receiver locations from the 10 km case.
Results use the unchanged tabulated dislocation source, no amplitude fitting,
no time shifting, and an eighth-order zero-phase 1 Hz filter with DC retained.
Stations 20–24 inside the small-domain PML are excluded.

| 100 s displacement comparison, common 25 interior receivers | Median relative L2 | Worst relative L2 |
| --- | ---: | ---: |
| Earlier gradient-correction CPML, cross fraction 0.1 | 11.47% | 32.99% |
| Complete CPML, cross fraction 0 | 1.10% | 2.69% |
| Coupled flat water/elastic, gravity off, cross fraction 0 | 0.65% | 2.06% |
| Coupled flat water/elastic, gravity on, cross fraction 0 | 1.68% | 3.33% |

Eight additional elastic receivers near the sides, bottom and corners give
1.88% median and 4.98% worst error across all 33 interior receivers. Keeping
the old 0.1 cross fraction with the complete operator gave 185% median error;
this is why the default changed. Raw elastic displacement gives 1.11% median
and 3.15% worst error on the common set. Raw velocity has larger differences
(8.57% median, 20.04% worst); the input source is not band-limited to the mesh
frequency. These filtered results do not validate unresolved high frequencies.

Separate CUDA packet tests call the production operator helpers on structured
Hex8 meshes and compare against larger domains. P, SV and SH incidence at
0°, 30°, 60° and 80° gives 0.13–3.47% relative L2 error with a 16-cell layer.
A free-top Rayleigh packet gives 0.102% (0.071% with half the timestep),
and a slow gravity-surface packet gives 0.228% (0.375% with half the timestep).
The latter change is not monotone convergence; it bounds these two tested
steps only. Kappa=3 remained finite in the tested grazing and surface cases,
but is not uniformly more accurate. These harnesses test the shared operator,
not the entire MPI solver.

A 1,000 s, 200,000-step elastic run remained finite at all 38 stations,
including absorbing-layer observations. Every station's velocity RMS over
900–1,000 s was at most 0.513 times its RMS over 100–200 s. For a 30 s
coupled gravity run, one versus four CUDA-aware ranks differed by 1.1e-9 in
aggregate displacement L2; the host MPI fallback differed by 5.8e-5. Both
backends remained finite. The experimental face_mass projection remained
finite and differed from point projection by 0.15%; this is a compatibility
check, not an independent accuracy result for face_mass.

Host checks cover repeated/zero poles, full bulk/surface inertia transfer
functions, independent Gauss-integrated elastic/acoustic stiffness, and the
zero-damping limit against the existing native elastic operator (roundoff).
The CUDA compiler reports no register spills in the new PML kernels.
A matched 10 km timing run with 38 stations and planes disabled reported
14.13 s total solver wall time for the earlier executable and 17.16 s for
complete plain CPML (about 21% more). Maximum-rank solver compute rose from
4.73 s to 11.10 s. These single runs on different nodes are indicative, not a
performance distribution. The added memory and computation are the cost of
retaining the previously omitted terms; production-size profiling remains
necessary.
Checkpoint/restart, single precision, arbitrary heterogeneous materials,
curved surfaces and attenuation combinations are not validated by this suite.

Investigation harnesses and results are local-only under
`.local-diagnostics/cpml-completion-20261002/`; they are not shipped in a clean
checkout. Horizon inputs, executable and results are retained under
`/scratch/08458/wz_tacc/hercules/cpml-completion-20261002`. Principal jobs are
42115–42123 (domain comparisons), 42117 (P/S packets), 42129 (gravity packet),
42130 (long run), 42136–42139 (backend/rank/projection and Rayleigh checks),
42145 (matched old-binary timing), and 42152 (final-default executable check).

### BKT2 GPU execution

The default BKT2 cube path shares the two displacement-history recurrences
between shear and bulk, then evaluates their force in the same kernel. This
preserves the existing relaxation rates, Q-dependent weights, reference-frequency
correction, timestep, arithmetic precision, underflow threshold and force
operator. The optimized dispatch requires device-resident force and displacement
synchronization. Host-staged execution (including checkpoint-driven backend
fallback) retains the previous kernels: an initial host-backend fused comparison
exceeded the waveform tolerance, so that path is not enabled. No extra per-step
synchronization is added to make it eligible. An initialization scan requires
spatially uniform relaxation rates
and exact agreement between active shear and bulk spectra. In BKT/BKT2, a sector
needs relaxation history only when `a0` or `a1` is nonzero; the unused `a2/g2`
fields cannot make it active. Generated element records initialize all damping
fields, and two-mechanism material assignment explicitly zeros `a2/g2`. The
eligibility scan also ignores these unused fields in imported records. BKT3/BKT3F
still include their third mechanism in the scan.

Inactive sectors need no history and are excluded from rate comparisons, so
active, shear-only, bulk-only and inactive elements can coexist on a fused rank.
Different Q-dependent weights remain per element and are compatible with fusion
when the active rates match. The existing BKT2 `Qs >= 1000` and `Qk >= 1000`
cutoffs, and `use_infinite_qk`, retain their previous physical behavior; this
eligibility fix does not change those attenuation approximations. Other BKT
variants and nonmatching spectra retain the previous kernels.

The shared histories occupy the shear slots for eligible cubes. This halves
history update traffic when both sectors attenuate. Transformed shear/bulk
vectors stay local to the fused kernel unless station stresses are requested;
then they are saved every step so the existing previous-step station gather
remains correct, including decimated output. VT cut cells retain the original
convolution and tetrahedral force. CPML elements retain their elastic CPML
operator. This change adds no field exchange, host staging or per-step allocation.
The legacy GPU buffer allocations remain in place for compatibility; reduced
traffic must not be reported as reduced allocated GPU memory.

The fused kernel launches 128-thread blocks with a six-block occupancy target.
This deliberately permits compiler spills to increase the number of resident
warps. The setting was selected by alternating CUDA-event measurements on
Horizon GB200 GPUs; the unconstrained register allocation was slower than the
separate kernels. Performance depends on the GPU and compiler, so retain the
reference switch for comparisons on other systems.

Set `HERCULES_BKT2_FUSED=0` in the launch environment to compare with the old
two-kernel BKT2 path. `HERCULES_FORCE_GENERIC_DAMPING=1` overrides both optimized
paths. The setup log reports actual fused ranks and fused/fallback element slots,
then rank counts for each fallback reason: fusion disabled, generic override,
force or displacement synchronization not device-resident, nonuniform active
shear or bulk rates, and a shear/bulk rate mismatch. Reason counts may overlap.
Ranks with no damping slots are reported separately and excluded from fused and
fallback-reason rank counts. Slots describe mapper entries, which can include
inactive, acoustic or CPML elements; they are not a count of attenuating work.
For the fused path, `Force damping local wait only` includes convolution and
cube force plus any VT force; `Force damping conv wait only` covers only the
remaining legacy convolution subset. Compare total solver time or a GPU
profiler's combined force region, not individual timer names between paths.
Host wait timers are not isolated kernel timings, and per-rank maxima are not
additive. Atomic nodal accumulation can change floating-point summation order,
so full waveform comparisons require a precision-appropriate tolerance rather
than a universal bitwise-identity claim.

Horizon GB200 validation (CUDA 13.2, NVIDIA HPC SDK 26.5): an isolated CUDA
regression compared the unchanged two kernels with the fused kernel over 128
steps and 263 independently assembled elements. Forces and saved transformed
vectors matched bitwise in both double- and single-precision builds, including
different shear/bulk weights, inactive sectors, bulk-only history, interface
projection flags and PML skips. A 2,097,152-element shared-node assembly check
had a peak-scaled force difference of `3.12e-16` in double precision.

For that larger synthetic assembly with both sectors active, no CPML/VT and
station-vector stores disabled, four alternating 1,000-step CUDA-event timing
pairs after warm-up gave median combined damping times of `1.110192 ms`
(separate kernels) and `0.854177 ms` (fused): **1.30x faster**, or 23.1% less
damping-kernel time. This is a kernel benchmark, not a full-solver or Izmir
speedup. It does not establish production-scale waveform accuracy or stability.

The small-solver regression comprised 31 ASCII-output runs and 12 float64
binary-displacement runs against the unmodified `2e1492a` executable, including
the updated executable's reference mode. Cases used one or four GPU/MPI ranks,
heterogeneous Q, shear-only/bulk-only/infinite-Q sectors, acoustic--elastic
coupling, CPML, VT cut cells and decimated strain/stress output. The final host
backend check exercises the retained reference kernels. Maximum errors scaled
by each station/field group's peak were `7.51e-7` for ASCII outputs (acceptance
limit `2e-6`, accounting for the text precision) and `1.72e-13` for binary
displacements (limit `1e-10`). Undefined solid-station pressure retained its
documented `nan` marker; finite physical fields were checked separately.

Single-pass small-solver reference/fused time ratios ranged from `0.917` to
`1.039`, including slower VT timing. These runs show no consistent end-to-end
speedup; use a matched production benchmark before assigning an Izmir speedup.
The local-only harness and retained results are under
`.local-diagnostics/bkt2-shared-fused-20261007/` (unavailable in a clean checkout).
The corresponding Horizon directory is
`/scratch/08458/wz_tacc/hercules/bkt2-shared-fused-20261007/`; jobs 50973 and 51046
contain the final CUDA-kernel checks and final guarded solver comparisons,
respectively. The latter solver source includes the independent acoustic
receiver-halo ordering fix from `5fd3d7e`.

The inactive-sector eligibility fix was checked against unmodified `42b55a7`
on the same Horizon toolchain with 24 ASCII and eight float64-displacement runs.
One- and four-rank cases included mixed active/shear-only/bulk-only/inactive
layers, different active Q weights, CPML, VT, acoustic--elastic coupling and the
host backend, with fused, separate and generic dispatch comparisons. Eight
additional runs deliberately seeded unused `a2` weights in isolated test builds:
the old four-rank solver fell back on all ranks, while the fixed solver fused
on all four. Across all 40 runs, maximum peak-scaled errors were `2.33e-9` for
ASCII fields and `7.94e-14` for binary displacements, below the same acceptance
limits above. Setup reason counts matched the requested overrides/backends.
Extracted production-code checks also covered nonzero/NaN unused fields,
genuine active-rate mismatches and BKT3/BKT3F third-mechanism mismatches.
These checks establish small-case equivalence, not full-size Izmir timing or
512-rank eligibility. Local-only evidence is under
`.local-diagnostics/bkt2-inactive-20261007/` (unavailable in a clean checkout),
with Horizon build/test jobs 51728, 51739 and 51744 under
`/scratch/08458/wz_tacc/hercules/bkt2-inactive-20261007/`.

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
`acoustic_stiffness_damping_beta`, with `f_ref` set by
`simulation_wave_max_freq_hz` and `dt` the simulation timestep. In the modal
interpretation of the discrete force term, `zeta(omega)=beta*omega*dt/2` and
`Q(omega)=1/(beta*omega*dt)`. This is stiffness-proportional damping: it is
stronger at higher frequencies. It is independent of the elastic `theBBase`
coefficient, which an earlier implementation incorrectly used, leaving the
option ineffective for nondamped/BKT cases. Acoustic CPML elements skip this
interior regularizer.

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

## ASCII station buffering

With the default `station_binary=no`, `station_ascii_buffer_kib` sets the stdio
buffer size for each station owned by an MPI rank. This optional `parameters.in`
integer defaults to **64 KiB** (1 KiB = 1024 bytes). Accepted values are
1..2097151 (`INT_MAX / 1024`); zero, negative, fractional and overflowing values
are rejected during input parsing. The value is validated but unused with
`station_binary=yes`. Existing inputs retain the 64 KiB default.

```ini
# CPU buffer per ASCII station, shared by all its output columns.
station_ascii_buffer_kib = 64
```

Buffer memory per rank is `station_ascii_buffer_kib * 1024 * local_station_count`
bytes, plus small stream/container overhead. For example, 1,000 stations owned
by one rank use 62.5 MiB at the default or 15.625 MiB with a value of 16. A larger
buffer allows less frequent writes but uses more CPU RAM and delays live-file
visibility. A smaller buffer reduces memory use and writes more frequently.
Ranks without stations and binary mode allocate no ASCII station buffers.
Buffers are installed before writing file headers and remain valid until the
files close; allocation failure stops setup with an error.

The writer retains the existing ASCII formatting, timestamps, sample cadence
(`output_stations_print_rate * simulation_delta_time_sec`), and all requested
motion, pressure, strain/stress and nonlinear columns. There is no filtering,
averaging or downsampling. Writes occur as the byte buffer fills rather than
forcing a flush every ten samples. The number of rows per write depends on the
enabled columns and their formatted values; the buffer size is a byte limit, not a fixed
sample count. The policy covers acoustic, topography and ordinary solid station
writers. Binary station output uses its separate existing MPI-IO path.

At normal completion, the solver flushes each remaining partial buffer, checks
write errors, closes the station files and releases the buffer memory. Final
station output is included in the `Print Stations` and `Solver I/O` timers.
The startup monitor reports the selected ASCII buffering policy. No per-step
allocation, field exchange or GPU transfer is added.

Live station files may lag the current simulation time and may end in a partial
text row until more data is written. An interrupted or aborted run can lose the
unflushed tail; inspect completed files for the full requested history. Larger
writes reduce small-write frequency but cannot guarantee fast output on a
congested filesystem. Plane output and the ETA estimator are independent of
this station-buffering change.

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
`[snapshot, station, XYZ]` in C order. XYZ means solver-local X/Y, Z down; see `coordinate_frame.txt` for the PROJ mapping.
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
zero-based. For zero strike and dip, the grid axes are solver X/Y; vector components always remain solver-local X/Y, Z/down.
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
repository.

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

## Acoustic station gradient recovery

With `acoustic_station_gradient=linear`, each component uses the owning
Hex8 derivative and one adjacent water-cell derivative evaluated at the same
transverse station coordinates. For direction $a$, let $x_a$ be the station
coordinate and $c_0$, $c_1$ the owning and adjacent cell centers along that axis:

```math
t = \frac{x_a-c_0}{c_1-c_0}, \qquad
(\partial_a p)_{\mathrm{rec}} = (1-t)(\partial_a p)_0 + t(\partial_a p)_1.
```

The neighbor toward the receiver from its cell center is preferred. If it is
unavailable, the opposite neighbor gives a one-sided extrapolation. The formula
reproduces the gradient of any quadratic scalar polynomial sampled at the mesh
nodes, including mixed terms and unequal Cartesian cell sizes. It does not
change the approximation space or the physical time-step operators.

Setup resolves actual local or remote leaves and registers missing vertices as
non-owned harbored nodes before solver allocation. Existing anchored/dangling
field exchanges, after field-specific assignment, supply those samples on both
host and CUDA-aware paths. Sender and receiver MPI mappings are sorted by
global node ID for every interface method, including
`shared-geometric-constraints` with point projection. This ordering is required
when receiver-only vertices are appended after mesh extraction. Stencils are
sorted by global node ID and uploaded
once; the fused GPU station kernel returns only the existing four values per
acoustic station. There is no new per-step MPI exchange or full-field download.

Stencils only use acoustic cells with matching Vp and density (relative
comparison tolerance 1e-6). Solid, air, void and different water materials are
excluded. An unavailable neighbor on both sides retains the element derivative
for that component. CPML stations retain the original derivative and output
stretch memories, and an interior stencil excludes CPML neighbors. Startup
reports added halo nodes, recovered components, fallbacks and CPML stations.

This setting changes acoustic station motion, including binary displacement.
It leaves pressure interpolation, plane reconstruction, elastic/poroelastic
station motion, stresses/strains, sources and propagation operators unchanged.
Use `element` for comparisons to historical station output or a colocated plane
point using the original gradient. Recovery cannot repair a distorted nodal
field, suppress physical high-frequency content, or establish CPML accuracy.
A local 2026-10-07 gravity-packet comparison used the same 10 km model,
600 s duration, independent acoustic constraints and unfiltered analytical
reference for both settings. At the center station 20 m above the bed, the
coarse-water (39.0625 m) vertical-displacement L2 error over 0.1--150 s decreased
from 2.247% to 0.423%; the fine-water (19.53125 m) error decreased from 48.106%
to 0.239%. Surface-displacement median errors decreased from 5.856% to 0.987%
and 3.024% to 0.278%, respectively. Pressure differences between paired runs
were at most about 1e-7 Pa (below 8e-8% relative L2), and the late surface-pressure errors
remained 16.761% and 9.271% for those CPML settings. These are output-recovery
results for a specific smooth, lossless packet, not general propagation or
absorber tolerances. A separate 24-station check at heights 1--100 m above
the bed still found median vertical-displacement errors of 14.91% (coarse)
and 3.85% (fine) at 1 m, where the reference vertical motion is small; recovery
does not guarantee uniform relative accuracy arbitrarily close to the
interface. The local fixture and results are not in a clean checkout.

For the same 110-station, 600 s runs, rank-maximum acoustic station gathering
(kernel plus result copy) increased from 1.994 to 2.852 s on the coarse case
and 1.929 to 2.787 s on the fine case: about 0.86 s extra per run. Total solver
time changed from 105.16 to 111.64 s and 195.78 to 195.05 s, respectively.
These are single runs on different allocations, not repeated measurements
on fixed nodes, so their overall runtime differences cannot be attributed
solely to reconstruction. The cases added 192 halo nodes; communication uses
existing exchanges with a larger payload. Cost depends on station count,
placement and the frequency at which station gradients are evaluated.

Analytically initialized benchmarks must seed their initial receiver velocity
with the selected observation formula; ordinary runs initialized at rest need
no additional initialization.

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
local X/Y metres from the domain origin. Lengths and depth describe
the **exterior**. For spacing h, the interior dimensions are
`(length_x - 2*h, length_y - 2*h, depth - h)`, sharing the same horizontal
center and top. Lengths must exceed 2h, depth must exceed h, and all three
exterior dimensions must be positive integer multiples of positive h.
Strike is degrees clockwise from solver +X: for box-local (x,y),
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

The matching `drmcoordinates.N` contains eight whitespace-separated columns:
`X_m Y_m Z_m Vs_m_s Vp_m_s density_kg_m3 Qs Qp`. The first three are the actual
solver sampling coordinates in metres: local X/Y, Z positive downward, after
strike rotation and the same building/topography depth shifts used for
interpolation. Ignoring lines beginning with `#`, data row i (zero-based)
matches binary sample i in every snapshot. Face ordering and repeated
edge/corner rows are identical; no points are sorted or deduplicated.

The five appended columns describe the **initial containing mesh element**:
Vs and Vp in m/s, density in kg/m³, and dimensionless Qs and Qp. They are not
point queries of the original velocity model or interpolation across a material
interface. Vs/Vp/density include the solver's mesh-material adjustments,
including applicable velocity cutoffs, softening, and phase-velocity correction.
At an interface, the same unique owning element used for displacement sampling
supplies the properties. All points within that element share its values.

The Qp column reports the **selected nominal P-wave Q**, after the BKT
shear/bulk branch choices, rather than the input Qp. The Qs column reports the
selected shear Q; legacy `bkt` uses the selected table-entry Q values. Disabled
branches contribute zero inverse Q after negative-Qk fallback, explicit
`use_infinite_qk=yes`, and the existing high-Q cutoffs. The reported Qp is
computed with the final Vs/Vp printed in the same row, including any
phase-velocity correction:

```math
L=\frac{4V_s^2}{3V_p^2},\qquad
\frac{1}{Q_{p,\mathrm{selected}}}=
\frac{L}{Q_{s,\mathrm{selected}}}+\frac{1-L}{Q_{\kappa,\mathrm{selected}}}.
```

For a disabled branch, its inverse Q is exactly zero; the internal Qk=1000
sentinel is not treated as finite bulk loss in this calculation. For example,
with final Vs=1000 m/s, Vp=1500 m/s, selected Qs=100 and lossless bulk, the
exported Qp is **168.75**, even if the input Qp was 200. A valid case with final
Vs=1000 m/s, Vp=2000 m/s, Qs=100 and selected Qk=400 reports Qp=200. These are
analytical illustrations, not production calibration. No input-Qp, Qk or
fallback-flag columns are added.

`inf` denotes zero selected material loss: Qs is `inf` when the shear branch
is disabled, and Qp is `inf` when both branches are disabled. With
`type_of_damping=none`, both are `inf` for ordinary solid materials.
Rayleigh/mass damping and building/foundation overrides report `nan` because
there is no selected frequency-independent BKT Q to report. BKT inputs still
must pass the existing admissibility checks; this reporting calculation does
not repair invalid inputs or alter damping coefficients.

These are nominal material-attenuation values, not a frequency-resolved
measurement of the fitted BKT response (including frequency-dependent `bkt3f`)
or a total effective Q for evolving nonlinear elastoplastic/plastic and
hysteretic loss, CPML absorption, or poroelastic flow. In nonlinear or other
regions where the time-stepping path excludes BKT, they describe the selected
linear reference material only. Poroelastic points report the mesh
skeleton/drained material inputs, not fast/slow Biot wave speeds or attenuation;
density is the mesh `rho`, not the separately configured fluid density.

Elastic, viscoelastic, nonlinear-solid and poroelastic solid samples use the same
metadata path. Acoustic–elastic, poroelastic–acoustic and poroelastic–elastic runs
can sample supported solid elements; scalar-acoustic samples remain unsupported
and stop DRM setup. This change does not modify initialization physics,
boundary/interface conditions, time integration, station quantities, or GPU
field exchanges. Target-system standalone and coupled CUDA/MPI comparisons are
still required; CPU formatting and synthetic MPI checks do not validate physics.

Values use 17 significant digits and a locale-independent decimal point.
The header describes units, ordering, point count, column names and Q semantics.
For illustration only, a homogeneous BKT material with final Vs=300 m/s,
Vp=700 m/s, density=1800 kg/m³ and selected Qs=25, Qp=50 would start as follows for the first box on
flat terrain with zero surface shift (these are synthetic values, not defaults
or calibrated production properties):

```text
# X_m Y_m Z_m Vs_m_s Vp_m_s density_kg_m3 Qs Qp
950 1950 0 300 700 1800 25 50
950 1950 1 300 700 1800 25 50
950 1950 2 300 700 1800 25 50
```

Solver rank zero writes the table once during setup after point ownership
validation and a bounded-chunk MPI material gather. Each run overwrites it,
including restarts. There is no additional per-step I/O or GPU transfer.
Position-only readers should use `usecols=(0, 1, 2)`; readers expecting exactly
three columns must be updated. The displacement binary format is unchanged.

**Rebuild the solver and regenerate older meshes with `mesh_input_mode=0`.**
Selected Qs/Qp persist as two double-precision fields in the element and saved mesh
payload (16 additional host bytes per element). Single-etree and sharded payload
sizes are checked; sharded manifests and MPI-IO archives use format version 2.
Older meshes are rejected before records are imported. The fields travel with
the existing mesh payload through single/sharded/MPI-IO output and import,
local interface refinement, and unified rank redistribution; no separate
quality-factor cache or new configuration parameter is required. Regenerate
meshes from earlier development builds that saved input Qp, and regenerate
after changing damping settings; import preserves the saved material selection.

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
metadata = np.loadtxt(path.with_name("drmcoordinates.0"), comments="#", ndmin=2)
if metadata.shape != (npoints, 8):
    raise ValueError("coordinate/material table does not match box dimensions")
xyz = metadata[:, :3]
vs, vp, density, qs, qp = metadata[:, 3:].T
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

The DRM/station MPI regression executables require `OUTPUT_PATH EXPECTED_RANKS`
(`1` or `4`) and reject a mismatched `MPI_COMM_WORLD` before opening files.
CI uses matching Open MPI compiler/launcher commands on Ubuntu 24.04, verifies
that rejection leaves an existing output file intact, and compares all output
bytes from the one- and four-rank runs. The former Ubuntu MPICH package/Hydra
combination could start independent singleton worlds; green jobs before the
rank-count guard did not establish multi-rank coverage. These are host output
regressions, not CUDA propagation or physics validation.
