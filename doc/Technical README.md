# Hercules GPU — Technical README

This is the concise operational reference and index. Detailed explanations,
formulas, configuration examples, implementation notes, and historical results
live in [Extended README](Extended%20README.md). The
[parameter reference](https://hercules-gpu-docs.liling-dpim.chatgpt.site/docs/parameter-file)
contains the complete current input catalog; the parser under
`quake/forward/` is authoritative.

## Build and launch

Use C/C++, CUDA, a supported NVIDIA GPU, MPI, GSL, PROJ, zlib, bzip2, and GNU
Make. Put machine-specific compiler/library/CUDA settings in `user.mk`, starting
from the closest [example](../examples/). Build with an explicit parallel job
count, such as `make -j 8`. `USE_OPENMP=1` is the default; use
`make USE_OPENMP=0` only to disable it. See
[compilation details](Extended%20README.md#host-and-cuda-compilation).

Run the repository launcher from the case directory:

```sh
export HERCULES_ROOT=/path/to/Hercules_GPU
cd /path/to/case
export PSOLVE_BIN="$HERCULES_ROOT/quake/forward/psolve"
export BASE_IN="$PWD/parameters.in"
export PHASE_A_NP=64
export PHASE_B_NP=8
export PHASE_B_OMP_NUM_THREADS=8
"$HERCULES_ROOT/run_phaseA_phaseB.sh"
```

`BASE_IN` is the complete input filename. Adjust ranks and threads for the
allocation and GPU count. The default `HERCULES_PHASE_MODE=unified` meshes and
hands the mesh to solver ranks in one launch; select `two-phase` for separate
launches. See [launcher details](Extended%20README.md#launch-a-simulation).

| Explicit two-phase input | `output_mesh` | `mesh_only_mode` | `mesh_input_mode` |
| --- | ---: | ---: | ---: |
| Generate the reusable mesh | 1 | 1 | 0 |
| Import and solve | 0 | 0 | 1 |

## Configuration essentials

Use one `parameters.in`; source-specific controls stay in `source.in` and its
auxiliary files. Start from a current [benchmark case](../benchmark_inputs/).
Use SI units, solver X north, Y east, and Z/depth positive downward.

- Domain extents must be compatible with a common dyadic octree level. The
  depth extent is deep minus shallow depth. Set `maximum_octant_level` from
  the active refinement requirements when using topography; see
  [domain sizing](Extended%20README.md#octree-compatible-domain-sizing).
- Elastic/viscoelastic, acoustic–elastic, poroelastic, and nonlinear solid
  models have different requirements. Nonlinear capability includes
  elastoplastic/plastic and hysteretic response. Calibrate geostatic state,
  tolerances, substeps, damping, and `material_properties_list` for the case.
- For acoustic runs use `include_acoustic=yes` and `acoustic_formulation`;
  `velocity_potential` is the default scalar formulation. The older
  `acoustic_pressure_formulation` key is a compatibility alias only when the
  newer key is absent.
- CPML defaults to `cpml_pml_material_damping=none` for BKT-family interiors.
  `rayleigh` is the only supported material-damping option inside the PML and
  requires matching interior damping.
  Check [current CPML behavior](Extended%20README.md#absorbing-boundaries-and-cpml).

## Input-variable index

This index covers the current `parameters.in` controls discussed in the
Extended README. Required/conditional status, defaults, units, ranges, and
interactions are in the
[full parameter catalog](../docs-site/app/docs/parameter-reference.ts) and the
[Extended README](Extended%20README.md).

- **Domain geometry:** `region_depth_shallow_m`, `region_length_east_m`, `region_length_north_m`,
  `region_depth_deep_m`.
- **Time integration and spatial resolution:** `simulation_delta_time_sec`,
  `simulation_wave_max_freq_hz`, `simulation_node_per_wavelength`, `simulation_shear_velocity_min`.
- **Base material model:** `cvmdb_input_file`.
- **Optional regional material models:** `Istanbul_velocity_model`, `basin_velocity_model`,
  `basin_origin_latitude_deg`, `basin_origin_longitude_deg`, `basin_input_file`,
  `basin_material_type`, `basin_water_vp`, `basin_water_rho`, `3D_velocity_model`,
  `3D_velocity_model_directory`.
- **Mesh generation, persistence, and import:** `output_mesh`, `mesh_only_mode`, `mesh_input_mode`,
  `mesh_etree_output_file`.
- **Material damping:** `type_of_damping`.
- **Absorbing boundary and CPML:** `abc_type`, `cpml_thickness_m`, `cpml_layers`,
  `cpml_reflection_target`, `cpml_poly_order`, `cpml_kappa_max`, `cpml_alpha_max_hz`,
  `cpml_cfl_safety`, `cpml_mpml_cross_fraction`, `cpml_alpha_floor_fraction`,
  `cpml_partition_weight`, `cpml_clamp_outer`, `cpml_3d_material_extension`, `cpml_stiffness_mode`,
  `cpml_stiffness_mode_acoustic`, `cpml_pml_material_damping`.
  With `abc_type=cpml`, `cpml_clamp_outer` defaults to `yes`; set it explicitly
  to `no` to disable the outer displacement/acoustic-scalar clamp. The complete
  operator defaults to `cpml_mpml_cross_fraction=0`; explicit older `0.1`
  settings remain active. Use `native` Hex8 integration and the `effective`
  GPU stiffness path. CPML requires `use_checkpoint=0` and `checkpointing_rate=0`
  because the checkpoint format omits absorber histories. See the
  [formulation and memory cost](Extended%20README.md#cpml-coefficient-and-operator-details).
- **Acoustic material and acoustic–elastic coupling:** `include_acoustic`, `acoustic_formulation`,
  `acoustic_free_surface_gravity`, `acoustic_pressure_formulation`, `acoustic_interface_coupling`,
  `acoustic_interface_method`, `acoustic_interface_projection`, `acoustic_interface_face_filter`, `acoustic_interface_max_level`,
  `acoustic_mesh_refinement_factor`, `acoustic_interface_diagnostics_file`,
  `acoustic_quality_factor`.
  `acoustic_quality_factor` defaults to `0`, disabling the interior acoustic
  damping term; use this lossless setting for baseline and validation comparisons.
  Positive values are numerical damping trials: `500` gives weaker damping than
  `100`, which was useful in an older flat-water benchmark. These are not physical
  water-attenuation calibrations. See [damping controls](Extended%20README.md#acoustic-and-elastic-coupling)
  and the [historical scope](Extended%20README.md#historical-acoustic-experiments).
  `acoustic_interface_method` defaults to `shared-geometric-constraints`; experimental
  `field-specific-constraints` separates acoustic, solid and pore-pressure hanging-node maps, while
  `local-matching` refines the coarser interface neighbors during Phase A.
  Imported mismatched meshes are rejected. Both alternatives require coupled
  velocity-potential acoustics and full Hex8 solids. Buildings, eligible nonlinear
  constitutive models and BrickUP poroelasticity can be enabled with either method;
  existing module-specific restrictions and material calibration requirements apply.
  Both methods use actual remote solid traces, including `face_mass` quadrature
  and open-pore pressure loads. The default shared-geometric method also resolves
  actual remote solid traces when `face_mass` is selected, so local and remote
  faces use the same four-point quadrature; its `point` path is unchanged.
  Field-specific constraints also
  support FEM topography (`include_topography=yes`, `computation_method=fem`),
  excluding air and carved voids from the field maps. VT cut-tetrahedron topography
  remains unsupported; local matching still requires topography off. CPML, gravity and
  `face_mass` remain independent; local matching does not guarantee CPML
  stability. See [constraints and matching](Extended%20README.md#interface-constraints-and-local-matching).
  `acoustic_interface_face_filter` defaults to `all`: couple every detected
  interface face using mesh-face normals. The `default` alias also selects
  `all`. Explicit `bathymetry` uses reconstructed normals and area weights on
  Z+ bottom faces; `bottom_only` uses vertical normals on Z+ faces.
  See [interface geometry and aliases](Extended%20README.md#acoustic-and-elastic-coupling).
  The flat water surface follows `topographybase_zcoord` (or the building
  `surface_shift_m`); check the exact exposed-face count in the startup log.
  Coupled scalar updates wait for MPI displacement sharing and hanging-node
  assignment. See [acoustic corrections](Extended%20README.md#scalar-acoustic-corrections-and-limitations).
- **Poroelastic u-p formulation:** `include_poroelastic`, `poroelastic_formulation`,
  `poroelastic_depth_max_m`, `surface_groundwatertable`, `groundwaterlevel_file`,
  `poroelastic_surface_pressure_bc`, `poroelastic_surface_pressure_penalty_factor`,
  `poroelastic_elastic_interface`, `poroelastic_undrained_fast_path`, `poroelastic_porosity`,
  `poroelastic_biot_alpha`, `poroelastic_fluid_bulk_modulus_pa`, `poroelastic_fluid_density_kgm3`,
  `poroelastic_fluid_viscosity_pa_s`, `poroelastic_permeability_x_m2`,
  `poroelastic_permeability_y_m2`, `poroelastic_permeability_z_m2`,
  `acoustic_poroelastic_interface`, `acoustic_poroelastic_pressure_penalty_factor`.
  Building cases measure ground and groundwater depths from `surface_shift_m`;
  above-ground building cells remain outside the saturated region.
- **Station output:** `number_output_stations`, `station_binary`, `station_output_binary_file`,
  `output_stations_print_rate`, `output_stations_directory`, `output_stations`,
  `station_mediums_output_file`, `print_station_pressures`, `print_station_pore_pressure_ratios`,
  `print_station_strains`, `print_station_stresses`.
- **Plane output:** `output_planes_print_rate`, `reconstruct_acoustic_plane_displacements`.
  Acoustic station/plane output inside CPML includes the containing element's
  gradient stretch and point memory; use interior samples for physical motion.
  Acoustic reconstruction uses physical pressure in both formulations; for
  velocity potential it first computes `p = -(phi[n] - phi[n-1])/dt`.
  See the [correction and validation note](Extended%20README.md#scalar-acoustic-corrections-and-limitations)
  for plane output from older executables.
- **Topography:** `include_topography`, `topography_directory`, `maximum_octant_level`,
  `computation_method`, `topography_bkt_gpu_stabilize`.
- **Nonlinear and plastic controls:** `include_nonlinear_analysis`, `geostatic_loading_time_sec`,
  `geostatic_cushion_time_sec`, `approximate_geostatic_state`, `error_tolerance`, `backbone_errTol`,
  `material_model`, `material_properties_list`, `no_substeps`, `rivasand02_reversaltype`,
  `manzari_dafalias_integration_scheme`, `manzari_dafalias_forward_euler_substeps`,
  `manzari_dafalias_recovery`.
- **Standalone material-point test controls:** `material_point_riva_density`,
  `material_point_riva_vs`, `material_point_brickup`.
- **Domain Reduction Method:** `implement_drm`, `number_output_drm_boxes`,
  `output_drm_boxes_print_rate`, `output_drm_boxes_directory`, `output_drm_boxes`,
  `drm_output_batch_size`, `drm_edgesize`, `drm_offset_x`, `drm_offset_y`, `drm_boundary`.
- **Incident plane waves:** `include_incident_planewaves`, `Ts`, `fc`.

### Source controls

The source-specific names covered by the detailed force-source examples are:

`average_risetime_sec`, `base_force_depth_m`, `base_force_direction_x`, `base_force_direction_y`,
`base_force_direction_z`, `base_force_rho_kgm3`, `base_force_scale`, `base_force_tolerance_m`,
`base_force_vs_mps`, `force_x_n`, `force_y_n`, `force_z_n`, `hypocenter_depth_m`,
`hypocenter_lat_deg`, `hypocenter_long_deg`, `lonlat_or_cartesian`, `ricker_Tp`, `ricker_Ts`,
`source_function_type`, `type_of_source`, `source_is_filtered`,
`threshold_frequency`, `number_of_poles`.

See [point-force configuration](Extended%20README.md#single-point-body-force-source)
and [compliant-base input](Extended%20README.md#distributed-compliant-base-force-source).

`source_is_filtered` is required (`0` or `1`); cutoff `threshold_frequency`
(Hz) and integer Butterworth `number_of_poles` are required when enabled and
ignored otherwise. SRFH filtering uses the `slipfunction.in` sample interval,
independently of the solver timestep. Allow pre-source time for the zero-phase
filter and check startup and final force; there is no universal source delay.
See [source-force filtering](Extended%20README.md#source-force-filtering) for
sampling, regeneration after the SRFH cutoff correction, and timing details.

## Runtime-variable index

These environment controls retain their existing meanings. Consult
[runtime controls](Extended%20README.md#runtime-performance-controls),
[diagnostics](Extended%20README.md#advanced-and-diagnostic-controls), and
[nonlinear execution](Extended%20README.md#nonlinear-gpu-execution-and-constitutive-models)
for defaults, supported combinations, and performance implications.

- **Launcher and CPU placement:** `HERCULES_PHASE_MODE`, `HERCULES_UNIFIED_IDLE_POLL_USEC`,
  `HERCULES_UNIFIED_MPI_YIELD_WHEN_IDLE`, `OMP_NUM_THREADS`, `OMP_PLACES`, `OMP_PROC_BIND`,
  `PHASE_A_NP`, `PHASE_A_OMP_NUM_THREADS`, `PHASE_B_KEEP_OMPI_MCA_BTL`, `PHASE_B_MAP_MODE`,
  `PHASE_B_NP`, `PHASE_B_OMPI_MCA_PML`, `PHASE_B_OMP_NUM_THREADS`, `PHASE_B_UCX_TLS`.
- **Mesh import and partition diagnostics:** `HERCULES_MESH_IMPORT_TREE_MODE`,
  `HERCULES_MESH_IMPORT_VERIFY`, `HERCULES_MESH_IO_MODE`, `OCTOR_EXTRACT_PROFILE`,
  `OCTOR_IMPORT_PATH_CACHE`, `OCTOR_INIT_TASKS_PER_RANK`, `OCTOR_NODE_SORT_MODE`,
  `OCTOR_SETCOM_ALL2ALL`.
- **Station and plane output:** `HERCULES_PLANE_IO_MODE`, `HERCULES_PLANE_IO_OMP`,
  `HERCULES_PLANE_MPIIO_BATCH`, `HERCULES_PLANE_METADATA`, `HERCULES_STATION_IO_OMP`,
  `HERCULES_STATION_ROUTE_MODE`.
  Plane output defaults to `mpiio`; `hybrid` and `legacy` remain explicit choices.
  Dedicated I/O ranks (`IO_PES>0`) select the legacy I/O-pool backend.
  `HERCULES_PLANE_METADATA=1` optionally exports containing materials and actual
  CPML membership once after solver setup (default: `0`).
- **Nonlinear execution:** `HERCULES_NONLINEAR_BACKEND`, `HERCULES_NONLINEAR_STATS_INTERVAL`.
- **Field synchronization and MPI:** `HERCULES_CONTRIB_PACK_OMP`, `HERCULES_CONTRIB_RECV_MODE`,
  `HERCULES_CONTRIB_UNPACK_OMP`, `HERCULES_DISP_SYNC_MODE`, `HERCULES_GPU_MPI_BACKEND`,
  `HERCULES_PORE_PRESSURE_FIELD_SYNC_MODE`, `HERCULES_PORE_PRESSURE_FORCE_SYNC_MODE`,
  `HERCULES_PORE_PRESSURE_SYNC_BACKEND`, `HERCULES_PRESSURE_FIELD_SYNC_MODE`,
  `HERCULES_PRESSURE_FORCE_SYNC_MODE`.
- **Other performance and diagnostic controls:** `HERCULES_ADJUST_ASSIGN_OMP`,
  `HERCULES_ADJUST_DIST_OMP`, `HERCULES_ALLOW_MULTI_RANK_PER_GPU`, `HERCULES_FORCE_GENERIC_DAMPING`,
  `HERCULES_FULL_CLEAR_OMP`, `HERCULES_REPRO_MODE`, `HERCULES_SPARSE_CLEAR_MODE`,
  `HERCULES_SPARSE_UPLOAD_MODE`, `HERCULES_SRC_DEBUG`, `HERCULES_TOPO_FORCE_MODE`.

## Output and validation

`station_binary=no` is the default and preserves ASCII station histories.
To record displacement from the existing `output_stations` list in one file:

```ini
station_binary = yes
station_output_binary_file = output/station_displacements.0
```

Use `number_output_stations` and `output_stations_print_rate` as before. The
binary is native float64 `[snapshot][station][XYZ]` in metres, in station-row
order. See [station binary details](Extended%20README.md#single-file-station-binary-output).
Plane files are native float64 displacement grids; formats and readers are in
the [output guide](https://hercules-gpu-docs.liling-dpim.chatgpt.site/docs/output-files)
and [script index](../script/README.md).
`HERCULES_PLANE_METADATA=1` writes `plane_point_metadata.RANK.csv` once during
setup, with actual point coordinates, containing-element material and CPML
membership. See [plane metadata and interpretation](Extended%20README.md#diagnose-large-plane-displacements).
For acoustic observations, enable ASCII station velocity, acceleration and
pressure output and verify ownership in `station_mediums.txt`. See
[pressure and motion interpretation](Extended%20README.md#acoustic-pressure-and-motion-diagnostics)
and [gravity-wave recovery limits](Extended%20README.md#check-gravity-wave-dispersion-and-reconstruction).
The [Izmir gravity investigation](izmir_gravity_diagnostic_20261001.md)
records surface-gradient recovery errors, the earlier face-average CPML inconsistency,
and source-corrupted experiments with a corrected-source follow-up. Validate
generated source forces before interpreting a failed boundary experiment.
Small operator tests or finite field output alone do not validate CPML
absorption or reconstructed gravity-surface motion.

Validate the chosen physics, mesh/time resolution, and GPU/MPI path with an
appropriate analytical or before/after comparison. Documentation checks do not
validate solver changes; historical benchmark numbers apply to their recorded
cases and versions. Keep `parameters.in`, source inputs, and monitor logs with
results. [Extended README](Extended%20README.md) holds the engineering detail;
[`docs-site/`](../docs-site/README.md) holds the public guides and parser audit.

## Structured DRM box output

Each `output_drm_boxes` row gives top-center X/Y, exterior lengths X/Y, depth,
spacing h, and clockwise strike from north. Hercules writes one
`drmdisplacements.N` containing five exterior and five interior face grids;
interior dimensions are `(Lx-2h, Ly-2h, depth-h)`. Configure
`number_output_drm_boxes`, `output_drm_boxes_print_rate`, and
`output_drm_boxes_directory`. This sampling is independent of legacy DRM
replay and does not refine the mesh.

Each box also writes `drmcoordinates.N` once during setup: an ASCII table of
solver X/north, Y/east, Z/down in metres, including terrain/building shifts.
Its rows match binary sample order exactly, including repeated edges/corners.

See [the complete geometry, layout, and reader examples](Extended%20README.md#structured-drm-box-output).
