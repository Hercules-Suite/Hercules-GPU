# Hercules GPU — Technical README

This is the concise operational reference and index. Detailed explanations,
formulas, configuration examples, implementation notes, and historical results
live in [Extended README](Extended%20README.md). The
[parameter reference](https://hercules-suite.github.io/Hercules-GPU/docs/parameter-file/)
contains the complete current input catalog; the parser under
`quake/forward/` is authoritative.

## Build and launch

Use C/C++, CUDA, a supported NVIDIA GPU, MPI, GSL, PROJ, zlib, bzip2, and GNU
Make. Put machine-specific compiler/library/CUDA settings in `user.mk`, starting
from the closest [example](../examples/). Build with an explicit parallel job
count, such as `make -j 8`. `USE_OPENMP`: `1` (default) enables OpenMP; `0`
disables it. See
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
allocation and GPU count. `HERCULES_PHASE_MODE=unified` (default) meshes and
hands the mesh to solver ranks in one launch; `two-phase` uses separate
launches. See [launcher details](Extended%20README.md#launch-a-simulation).

| Explicit two-phase input | `output_mesh` | `mesh_only_mode` | `mesh_input_mode` |
| --- | ---: | ---: | ---: |
| Generate the reusable mesh | 1 | 1 | 0 |
| Import and solve | 0 | 0 | 1 |

## Configuration essentials

The startup **Some Input Data** block records the build's Git commit and tracked
tree state, frequency, `simulation_node_per_wavelength`, timestep, simulation
start/end times and duration, domain, source paths, and active topography, 3-D
material, seabed and output-directory settings. The revision is embedded when
building; source copies without Git metadata report `unknown`. Rebuild to
update it. See [build provenance and startup summary](Extended%20README.md#build-provenance-and-startup-summary).

Routine configuration, progress and timing go to stdout; warnings and errors
go to stderr. Stdout is line-buffered for batch jobs. Preserve both scheduler
logs plus `monitor_file`, which retains its existing monitor messages. See
[console output streams](Extended%20README.md#console-output-streams).

Use one `parameters.in`; source-specific controls stay in `source.in` and its
auxiliary files. Start from a current [benchmark case](../benchmark_inputs/).
Use SI units and Z/depth positive downward. At zero PROJ azimuth, solver
X/Y follow UTM north/east; rotated runs use local X/Y axes.

- Domain extents must be compatible with a common dyadic octree level. The
  depth extent is deep minus shallow depth. Set `maximum_octant_level` from
  the active refinement requirements when using topography; see
  [domain sizing](Extended%20README.md#octree-compatible-domain-sizing).
- Elastic/viscoelastic, acoustic–elastic, poroelastic, and nonlinear solid
  models have different requirements. Nonlinear capability includes
  elastoplastic/plastic and hysteretic response. Calibrate geostatic state,
  tolerances, substeps, damping, and `material_properties_list` for the case.
- For acoustic runs use `include_acoustic=yes` and `acoustic_formulation`;
  `velocity_potential` (default) is the scalar formulation. The older
  `acoustic_pressure_formulation` key is a compatibility alias only when the
  newer key is absent.
- CPML uses `cpml_pml_material_damping=none` (default) for BKT-family interiors.
  `rayleigh` is the only supported material-damping option inside the PML and
  requires matching interior damping.
  Check [current CPML behavior](Extended%20README.md#absorbing-boundaries-and-cpml).

## Input-variable index

Every row lists the supported choices or value shape; `(default)` follows the
value selected when the input is omitted. **No default** means the value must
be supplied when its feature is active. Numeric values are in SI units unless
stated otherwise; a table or file path has a format rather than a finite option
list. Defaults that depend on another option are qualified in the same row.
Example/calibration values are not defaults. The tables cover all 249 controls
in the [current parser catalog](../docs-site/app/docs/parameter-reference.ts),
including the separately identified material-point test inputs. Use the
[Extended README](Extended%20README.md) for detailed constraints and examples.

### Domain geometry

| Variable | Options / accepted values | Requirement / units |
| --- | --- | --- |
| `region_origin_latitude_deg` | Any valid latitude; no default | Required; Float (degrees) |
| `region_origin_longitude_deg` | Any valid longitude; no default | Required; Float (degrees) |
| `region_depth_shallow_m` | Usually >= 0; no default | Required; Float (m) |
| `region_length_east_m` | > 0; octree-compatible; no default | Required; Float (m) |
| `region_length_north_m` | > 0; octree-compatible; no default | Required; Float (m) |
| `region_depth_deep_m` | > region_depth_shallow_m; no default | Required; Float (m) |
| `region_azimuth_leftface_deg` | Finite degrees; no default | Required; Float (degrees) |
| `domain_surface_corners` | 4 rows of `longitude latitude` (degrees); required without PROJ; no default | Conditional; Table: 4 rows × 2 floats |

PROJ azimuth is clockwise from UTM grid north;
rotate about the origin and retain local X/Y lengths. See
[domain rotation and output conversion](Extended%20README.md#rotated-domains-and-material-grids).

### Time integration and spatial resolution

| Variable | Options / accepted values | Requirement / units |
| --- | --- | --- |
| `simulation_start_time_sec` | `0` (default); >= 0 | Optional; Float (s) |
| `simulation_end_time_sec` | >= simulation_start_time_sec; no default | Required; Float (s) |
| `simulation_delta_time_sec` | > 0; no default | Required; Float (s) |
| `simulation_wave_max_freq_hz` | > 0; no default | Required; Float (Hz) |
| `simulation_node_per_wavelength` | > 0; no default | Required; Integer |
| `simulation_shear_velocity_min` | > 0; no default | Required; Float (m/s) |

See the [detailed parameter entries](../docs-site/app/docs/parameter-reference.ts).

### Base material model

| Variable | Options / accepted values | Requirement / units |
| --- | --- | --- |
| `cvmdb_input_file` | `profile` or material-database file path; no default | Required; Path or keyword |
| `number_profile_layers` | > 0 when cvmdb_input_file=profile; no default | Conditional; Integer |
| `profile_layers` | N rows of `z_top Vp Vs rho Qp Qs`; m, m/s, m/s, kg/m³, dimensionless Q; no default | Conditional; Table: N rows × 6 floats |
| `the_threshold_Vp_over_Vs` | >= 0; no default | Required; Float |
| `softening_factor` | `0` (default); either 0 or > 1 | Optional; Float |
| `stiffness_calculation_method` | `effective` (default), `conventional` | Optional; String |
| `use_progressive_meshing` | `1` (default); >= 0 | Optional; Integer |

See the [detailed parameter entries](../docs-site/app/docs/parameter-reference.ts).

### Optional regional material models

| Variable | Options / accepted values | Requirement / units |
| --- | --- | --- |
| `Istanbul_velocity_model` | `no` (default), `yes` | Optional; Yes/no |
| `basin_velocity_model` | `no` (default), `yes` | Optional; Yes/no |
| `basin_origin_latitude_deg` | Float (degrees); Required when basin_velocity_model=yes; no default | Conditional; Float (degrees) |
| `basin_origin_longitude_deg` | Float (degrees); Required when basin_velocity_model=yes; no default | Conditional; Float (degrees) |
| `basin_input_file` | Path; Required when basin_velocity_model=yes; no default | Conditional; Path |
| `basin_material_type` | `elastic` (default), `acoustic_water` | Optional; String |
| `basin_water_vp` | `1500` (default); > 0 | Optional; Float (m/s) |
| `basin_water_rho` | `1000` (default); > 0 | Optional; Float (kg/m3) |
| `3D_velocity_model` | `no` (default), `yes` | Optional; Yes/no |
| `3D_velocity_model_directory` | Directory path; `inputfiles/materialfiles` (default) | Optional; Directory path |
| `3D_velocity_model_azimuth_deg` | Finite degrees; inherit `region_azimuth_leftface_deg` (default); nonzero explicit overrides require PROJ | Optional; Float (degrees) |
| `3D_velocity_model_origin_latitude_deg` | Float (degrees); Required when 3D_velocity_model=yes; no default | Conditional; Float (degrees) |
| `3D_velocity_model_origin_longitude_deg` | Float (degrees); Required when 3D_velocity_model=yes; no default | Conditional; Float (degrees) |

Material axes follow solver axes by
default; explicit `0` keeps UTM north/east axes with PROJ. Without PROJ,
omission or explicit `0` retains solver-aligned axes; explicit nonzero
overrides require PROJ. Rotated PROJ domains require a profile
background and disable the legacy Istanbul reader. Rebuild meshes/sources
when either frame changes; PROJ output directories include `coordinate_frame.txt`.
The 3-D material reader supports 64-bit sample counts/profile offsets with
chunked MPI initialization; full material arrays remain replicated per rank.
See [file layout and memory sizing](Extended%20README.md#3d-velocity-model-files).
Without explicit seabed geometry, `type_of_etree=sqd` uses `z.in` depths below each column's
surface and samples at equal physical elevation using the actual query-point
terrain. `full` uses fixed model Z without a terrain shift. Rebuild the solver
and regenerate meshes with `mesh_input_mode=0` to apply the corrected lookup;
see [depth conventions and example](Extended%20README.md#3d-material-depths-with-topography).

### Explicit seabed geometry (opt-in)

| Variable | Options / accepted values | Requirement / units |
| --- | --- | --- |
| `seabed_geometry_file` | File path; omitted (default), retaining legacy Vs-based classification | Optional; File path |
| `seabed_height_tolerance_m` | `0.01` (default); finite >= 0; `0` means exact comparisons | Optional; metres; requires `seabed_geometry_file` when supplied |
| `seabed_above_surface_policy` | `error` (default), `clip`; case-insensitive | Optional; requires `seabed_geometry_file` when supplied |
| `3D_velocity_model_depth_reference` | `legacy` (default) without seabed geometry, `model_z`, `surface`, `seabed`; explicit non-legacy choice required with geometry | Conditional; String |
| `seabed_water_vp_m_s` | Float (m/s); Required with seabed_geometry_file; finite > 0; no default | Conditional; Float (m/s) |
| `seabed_water_density_kg_m3` | Float (kg/m3); Required with seabed_geometry_file; finite > 0; no default | Conditional; Float (kg/m3) |

The optional seabed file contains spacing followed by elevations only; one number per
line is valid. Water depth and coverage come from the outer surface minus bed
(`seabed_height_tolerance_m`, `0.01` m (default)). With `include_topography=no`, the outer surface is
flat at solver Z=0 and no topography file is loaded. Spacing must divide both
horizontal domain lengths; bed elevations are negative below this zero datum.
The domain must have zero CVM mesh offsets (normally `cvmdb_input_file=profile`).
Topography refinement is skipped; wavelength/interface refinement and 2:1
balancing remain. Omitting the seabed file retains Vs-based water classification.
Set `3D_velocity_model_depth_reference` explicitly to `seabed`, `surface`,
or `model_z`; `legacy` (default) applies without the file. With topography on,
each explicit reference works with either `full` or `sqd`. With it off,
`surface` and `model_z` coincide, `seabed` remains bed-relative, and
`type_of_etree` and other topography-only controls are unused. Supply `seabed_water_vp_m_s` and
`seabed_water_density_kg_m3` (required, finite positive SI values). Requires
acoustics and solid-only 3-D profiles; buildings and other
material-model overrides are unsupported. With topography on, dry negative
elevations remain dry and separate flat water bodies may have different elevations.
Porous depth limits start at the bed; wet hydrostatic pressure and acoustic free
surfaces use Z=0 with topography off, or local topography with it on.
`seabed_above_surface_policy=clip` lowers bed elevations above the outer surface
to that surface in memory, reports the count and maximum lowering, and preserves
the input file. Corrected nodes are dry; valid water and the flat-water-surface
requirement remain. `error` (default) rejects crossings beyond the configured tolerance.
The same tolerance defines dry gaps, nodal snapping and the maximum corner-height
range in wet cells; it does not flatten the surface or change mesh-face alignment
checks. Use `1e-8` to restore the former tolerance. Larger tolerances can remove
real shallow water and change geometry-relative depths and loads.
Regenerate meshes and sources after changes. See
[file example, physics scope and limits](Extended%20README.md#explicit-seabed-geometry).

[`script/taper_topography.py`](../script/README.md#taper_topographypy) automatically
checks a sibling input `seabed.in` and repairs taper-induced inconsistencies
using the original wet/dry nodes. `--seabed PATH` selects another bed;
`--no-seabed` disables handling. `--topography-base` optionally supplies the
absolute datum. It retains wet beds and writes the paired files plus a summary
beside `--out`, preserving inputs even with `--force`. Validation uses a fixed
`1e-8 m` tolerance. For already tapered surfaces, the separate
[`repair_tapered_geometry.py`](../script/README.md#repair_tapered_geometrypy)
command remains available. See the usage guide for arguments and limits.

### Source directories

| Variable | Options / accepted values | Requirement / units |
| --- | --- | --- |
| `source_directory` | Directory path; Required for source-driven runs; no default | Conditional; Directory path |
| `source_directory_output` | Directory path; `outputfiles/srctmp` (default) | Optional; Directory path |

### Mesh generation, persistence, and import

| Variable | Options / accepted values | Requirement / units |
| --- | --- | --- |
| `output_mesh` | `0`, `1` (default) | Optional; Boolean integer |
| `mesh_only_mode` | `0` (default), `1` | Optional; Boolean integer |
| `mesh_input_mode` | `0` (default), `1` | Optional; Boolean integer |
| `mesh_etree_output_file` | File path; Required for mesh output or mesh import; no default | Conditional; File path |

See the [detailed parameter entries](../docs-site/app/docs/parameter-reference.ts).

### Material damping

| Variable | Options / accepted values | Requirement / units |
| --- | --- | --- |
| `type_of_damping` | `rayleigh`, `mass`, `none`, `bkt`, `bkt2`, `bkt3`, `bkt3f`; no default | Required; String |
| `simulation_velocity_profile_freq_hz` | `0` (default); 0 to simulation_wave_max_freq_hz | Optional; Float (Hz) |
| `use_parametric_q_factor` | `no`, `yes`; no default; Required for type_of_damping=bkt* | Conditional; Yes/no |
| `use_infinite_qk` | `no`, `yes`; no default; Required for type_of_damping=bkt* | Conditional; Yes/no |
| `parametric_q_factor_constant` | Float; Required when use_parametric_q_factor=yes; no default | Conditional; Float |
| `parametric_q_factor_alpha` | Float; Required when use_parametric_q_factor=yes; no default | Conditional; Float |
| `parametric_q_factor_beta` | Float; Required when use_parametric_q_factor=yes; no default | Conditional; Float |
| `the_threshold_damping` | >= 0; no default | Required; Float |
| `do_damping_statistics` | `0` (default), `1` | Optional; Boolean integer |

Optional 3-D `qs.in`/`qp.in` use solid-only interpolation at water/solid
stencils; supplied Q is independent of the Vs minimum. Setup rejects invalid
solid Q and inadmissible BKT2 coefficients without another mesh traversal.
With `use_infinite_qk=no`, only elements whose Qp/Qs imply negative bulk
attenuation use the existing `Qk=1000` lossless-bulk fallback. Valid bulk
attenuation and shear attenuation remain unchanged; setup reports the global
fallback count. The supplied Qp is no longer reproduced in fallback elements.
`use_infinite_qk=yes` suppresses bulk attenuation for all solids. Zero inverse
bulk Q remains lossless with either setting. Regenerate persisted meshes to
apply the policy; see
[Setup-time attenuation validation](Extended%20README.md#setup-time-attenuation-validation)
and [Quality-factor files](Extended%20README.md#quality-factor-files).

### Absorbing boundary and CPML

| Variable | Options / accepted values | Requirement / units |
| --- | --- | --- |
| `abc_type` | `dashpot` (default): solid Lysmer and scalar acoustic radiation damping; `cpml` | Optional; String |
| `cpml_thickness_m` | `0` (default), automatic thickness; positive physical width overrides `cpml_layers`; nonpositive values select automatic width | Optional; Float (m) |
| `cpml_layers` | Integer; `10` (default); must be positive when thickness is automatic | Optional; Integer |
| `cpml_reflection_target` | `0.01` (default); 0 < R < 1 | Optional; Float |
| `cpml_reference_vp` | Omitted or nonpositive number → automatic maximum Vp over PML elements (default); positive number overrides with one uniform reference speed | Optional; Float (m/s) |
| `cpml_poly_order` | `2` (default); 1 to 4 | Optional; Integer |
| `cpml_kappa_max` | `5` (default); >= 1 | Optional; Float |
| `cpml_alpha_max_hz` | Omitted or nonpositive number → `0.5` Hz (default), independent of `simulation_wave_max_freq_hz`; positive number overrides; `auto` is not an input token | Optional; Float (Hz) |
| `cpml_cfl_safety` | `0.9` (default); 0 < value <= 1 | Optional; Float |
| `cpml_mpml_cross_fraction` | `0` (default); 0 to 0.5 | Optional; Float |
| `cpml_acoustic_sponge_fraction` | `0` (default); finite and >= 0; 0 disables | Optional; Float (dimensionless) |
| `cpml_alpha_floor_fraction` | `0.1` (default); 0 to 1 | Optional; Float |
| `cpml_profile_slab` | `none` (default), `auto` | Optional; String (case-insensitive) |
| `cpml_partition_weight` | `1` (default); >= 1 | Optional; Float. CPML weighting applies only to the completed solver mesh; construction and mesh-only runs use ordinary partitions. Check startup GPU memory reports; this is not a memory cap. |
| `cpml_clamp_outer` | `no`, `yes` (default) | Optional; Yes/no |
| `cpml_3d_material_extension` | `no`, `yes` (default) | Optional; Yes/no |
| `cpml_stiffness_mode` | `native` (default), `reduced` | Optional; String |
| `cpml_stiffness_mode_acoustic` | `native`, `reduced`; inherit `cpml_stiffness_mode` (default) | Optional; String |
| `cpml_pml_material_damping` | `none` (default) except for Rayleigh interiors; `rayleigh` (default) with Rayleigh interiors when omitted; explicit `none` is retained | Optional; String |
| `cpml_gravity_in_bottom_band` | `abort` (default), `allow` (case-insensitive) | Optional; setup policy for gravity-surface acoustic CPML elements. Separately reports geometric bottom-band overlap and lateral M-PML cross-z damping; `allow` warns and continues for diagnostics. Known unstable bottom-band harness cases motivate the guard; no universal stability claim. |

With `dashpot`, scalar acoustic water now absorbs on the four outer sides
and bottom, using local `rho*Vp` impedance. Pressure and velocity-potential
formulations use centered damping; the top free surface and internal coupling
conditions retain their existing rules. No new input is required. This also
changes water boundaries in cases that omit `abc_type`; old dashpot water
results used reflecting outer faces. The condition targets outgoing acoustic
waves, permits oblique reflection, and is not a specialized surface-gravity-wave
absorber. It requires a `BOUNDARY` build, like the solid dashpot. See
[Acoustic dashpot boundary](Extended%20README.md#acoustic-dashpot-boundary)
for the formulation, applicability and validation limits.

Automatic CPML thickness is `cpml_layers` times the largest non-air element
edge on a lateral or bottom boundary (MPI-wide maximum); air is also
excluded from the global fallback. Finer boundary regions contain more
cells; use positive `cpml_thickness_m` to specify the physical width.
With `abc_type=cpml`, `cpml_clamp_outer=yes` (default); set it explicitly
to `no` to disable the outer displacement/acoustic-scalar clamp. The complete
operator uses `cpml_mpml_cross_fraction=0` (default); explicit older `0.1`
settings remain active. Use `native` Hex8 integration and the `effective`
GPU stiffness path. CPML requires `use_checkpoint=0` and `checkpointing_rate=0`
because the checkpoint format omits absorber histories. With
`acoustic_free_surface_gravity=yes`, the guard counts bottom-band surface
elements separately from surface elements outside that band with lateral
M-PML cross-z damping. Remove the first trigger by deepening the domain or
thinning the band; remove the second with `cpml_mpml_cross_fraction=0`.
`cpml_gravity_in_bottom_band=allow` bypasses the setup stop for diagnostics.
Disabling gravity changes the boundary physics. Startup also prints
`CPML acoustic-elastic interface:` fragment and PML-fragment counts with the
coupling assumption. Zero PML fragments produce an informational stdout line,
not a missing-correction warning: an interface may lie entirely outside the
PML even when water reaches it. A missing host side mask for enrolled PML
fragments produces a `CPML WARNING` on stderr. Both streams retain their
`monitor_file` copies. These host setup counts do not verify GPU execution
or long-time stability. See the
[formulation and memory cost](Extended%20README.md#cpml-coefficient-and-operator-details)
and [guard evidence and validation limits](Extended%20README.md#absorbing-boundaries-and-cpml).

With `cpml_profile_slab=auto` (opt-in) the damping profile is sampled per
octree-aligned slab of the coarsest non-air PML element edge, so every element
size in a slab shares one `d`/`alpha`/`kappa`; `cpml_reference_vp` keeps one
uniform reference speed. Startup prints the maximum per-cell stretch
`d_max*h/V` and warns above 10. See the
[formulation, profile sampling and memory cost](Extended%20README.md#cpml-coefficient-and-operator-details).

### Acoustic material and acoustic–elastic coupling

| Variable | Options / accepted values | Requirement / units |
| --- | --- | --- |
| `include_acoustic` | `no` (default), `yes` | Optional; Yes/no |
| `acoustic_formulation` | `velocity_potential` (default), `pressure`, `none`; `potential` is an alias for `velocity_potential` | Optional; String |
| `acoustic_free_surface_gravity` | `no` (default), `yes` | Optional; Yes/no |
| `acoustic_pressure_formulation` | Omitted (default), using `acoustic_formulation`; `yes` → `velocity_potential`, `no` → `none`; ignored if the newer key is present | Legacy; Yes/no |
| `acoustic_interface_coupling` | `no`, `yes` (default) | Optional; Yes/no |
| `acoustic_interface_method` | `field-specific-constraints` (default) for coupled scalar acoustics, `shared-geometric-constraints` (default) when acoustics, scalar acoustics or coupling is disabled, `local-matching` | Optional; String (case-insensitive) |
| `acoustic_interface_projection` | `point` (default), `face_mass`; aliases: `legacy` / `one_point` → `point`, `mass` / `mass_projected` / `weak` → `face_mass` | Optional; String |
| `acoustic_interface_face_filter` | `all` (default), `bathymetry`, `bottom_only`; aliases: `default` / `none` → `all`, `reconstructed_bathymetry` / `bathymetry_normal` → `bathymetry`, `bottom` / `z_plus` → `bottom_only` | Optional; String |
| `acoustic_interface_max_level` | `-1` (default); -1 to PIXELLEVEL | Optional; Integer octree level |
| `acoustic_mesh_refinement_factor` | `1` (default); >= 1 | Optional; Float |
| `acoustic_interface_diagnostics_file` | File path; omitted (default), no diagnostic file | Optional; File path |
| `acoustic_station_gradient` | `linear` (default), `element` | Optional; String (case-insensitive) |
| `acoustic_quality_factor` | Number; `0` (default); nonpositive disables interior damping, positive enables it | Optional; Float (dimensionless) |

The lossless acoustic-quality-factor default was verified in the documented
flat-water CPS comparisons and MPI/backend consistency tests; see the
[validation scope](Extended%20README.md#interface-constraints-and-local-matching).
`acoustic_station_gradient=linear` (default) recovers station motion at the
actual point using same-material neighboring water cells; `element` retains
the historical gradient. CPML and unavailable-neighbor fallbacks are described
in [the detailed reference](Extended%20README.md#acoustic-station-gradient-recovery).
`acoustic_interface_method=field-specific-constraints` (default) for coupled
scalar acoustics separating acoustic, solid and pore-pressure hanging-node maps.
When the method is omitted and acoustics, scalar acoustics or coupling is disabled,
shared geometric constraints remain active for acoustic and solid fields.
Pore-pressure constraints are independent as described below.
Explicit unsupported choices still abort.
The alternative `local-matching` refines the coarser interface neighbors during Phase A.
Imported mismatched meshes are rejected. Both alternatives require coupled
pressure or velocity-potential acoustics and full Hex8 solids. Buildings, eligible nonlinear
constitutive models and BrickUP poroelasticity can be enabled with either method;
existing module-specific restrictions and material calibration requirements apply.
Both methods use actual remote solid traces, including `face_mass` quadrature
and open-pore pressure loads. The shared-geometric method also resolves
actual remote solid traces when `face_mass` is selected, so local and remote
faces use the same four-point quadrature. `open_pore` coupling
also uses actual solid traces with shared-geometric `point` projection so
remote pressure loads reach independent pore nodes.
Field-specific constraints also
support FEM topography (`include_topography=yes`, `computation_method=fem`),
excluding air and carved voids from the field maps. VT cut-tetrahedron topography
requires explicit `shared-geometric-constraints` for coupled acoustics; local matching
still requires topography off. CPML, gravity and
`face_mass` remain independent; local matching does not guarantee CPML
stability. See [constraints and matching](Extended%20README.md#interface-constraints-and-local-matching).
`acoustic_interface_face_filter=all` (default) couples every detected
interface face using mesh-face normals. The `default` alias also selects
`all`. Explicit `bathymetry` uses reconstructed normals and area weights on
Z+ bottom faces; `bottom_only` uses vertical normals on Z+ faces.
See [interface geometry and aliases](Extended%20README.md#acoustic-and-elastic-coupling).
With explicit seabed geometry, flat water surfaces follow local topography
when enabled, or solver Z=0 when `include_topography=no`.
Otherwise they follow `topographybase_zcoord` (or building `surface_shift_m`).
Check the exact exposed-face count in the startup log.
Coupled scalar updates wait for MPI displacement sharing and hanging-node
assignment. See [acoustic corrections](Extended%20README.md#scalar-acoustic-corrections-and-limitations).

### Poroelastic u-p formulation

| Variable | Options / accepted values | Requirement / units |
| --- | --- | --- |
| `include_poroelastic` | `no` (default), `yes` | Optional; Yes/no |
| `poroelastic_formulation` | `brick_up` (default), alias `brickup`; only supported formulation | Optional; String |
| `poroelastic_depth_max_m` | `30` (default); >= 0 | Optional; Float (m) |
| `surface_groundwatertable` | `no` (default), `yes` | Optional; Yes/no |
| `groundwaterlevel_file` | File path; `3D_velocity_model_directory/groundwaterlevel.in` (default) when the groundwater grid is enabled | Optional; File path |
| `poroelastic_surface_pressure_bc` | `undrained` (default), `drained` | Optional; String |
| `poroelastic_surface_pressure_penalty_factor` | `1` (default); >= 0 | Optional; Float |
| `poroelastic_elastic_interface` | `no_flow` (default), alias `noflow`; only supported interface | Optional; String |
| `poroelastic_undrained_fast_path` | `yes` (default), `no`; automatically inactive with nonlinear analysis | Optional; Yes/no |
| `poroelastic_porosity` | `0.4` (default); 0 < value <= 1 | Optional; Float |
| `poroelastic_biot_alpha` | `1` (default); 0 to 1 | Optional; Float |
| `poroelastic_fluid_bulk_modulus_pa` | `2.2e9` (default); > 0 | Optional; Float (Pa) |
| `poroelastic_fluid_density_kgm3` | `1000` (default); > 0 | Optional; Float (kg/m3) |
| `poroelastic_fluid_viscosity_pa_s` | `1e-3` (default); > 0 | Optional; Float (Pa·s) |
| `poroelastic_permeability_x_m2` | `0` (default); >= 0 | Optional; Float (m2) |
| `poroelastic_permeability_y_m2` | `0` (default); >= 0 | Optional; Float (m2) |
| `poroelastic_permeability_z_m2` | `0` (default); >= 0 | Optional; Float (m2) |
| `acoustic_poroelastic_interface` | `no_flow` (default), `open_pore`; no-flow aliases `default`, `noflow`, `sealed`; open-pore aliases `open`, `pressure_continuity` | Optional; String |
| `acoustic_poroelastic_pressure_penalty_factor` | `1` (default); > 0 for open_pore | Optional; Float |

`open_pore` always uses four-point (2x2 Gauss) hydraulic integration for both
the pressure target and boundary capacity, independently of mechanical
`point`/`face_mass`. No new parameter is required; the penalty factor and
material inputs are unchanged. Earlier centroid-based benchmark results do
not validate this upgrade; target-system CUDA/MPI checks remain required.
See [hydraulic quadrature](Extended%20README.md#open-pore-hydraulic-quadrature).
Building cases measure ground and groundwater depths from `surface_shift_m`;
above-ground building cells remain outside the saturated region.
Regular Hex8 and FEM-topography models automatically use a separate pressure
constraint map from final active BrickUP cells, including `include_acoustic=no`.
Solid displacement retains its selected constraints; VT cut-cell pressure
retains geometric constraints with a startup notice. See
[poroelastic controls](Extended%20README.md#poroelastic-u-p-formulation).

### Monitor and distributed field output

| Variable | Options / accepted values | Requirement / units |
| --- | --- | --- |
| `monitor_file` | File path; `monitor.txt` (default) | Optional; File path |
| `simulation_output_rate` | `1000000` (default); > 0 | Optional; Integer time steps |
| `output_parallel` | `0` (default), `1` | Optional; Boolean integer |
| `output_displacement` | `0` (default), `1`; effective only when output_parallel=1 | Optional; Boolean integer |
| `output_displacement_file` | File path; Required when output_parallel=1 and output_displacement=1; no default | Conditional; File path |
| `output_velocity` | `0` (default), `1`; effective only when output_parallel=1 | Optional; Boolean integer |
| `output_velocity_file` | File path; Required when output_parallel=1 and output_velocity=1; no default | Conditional; File path |
| `output_stats_file` | File path; `output-stats.txt` (default) | Optional; File path |
| `output_debug` | `0` (default), `1`; effective only when output_parallel=1 | Optional; Boolean integer |
| `output_debug_file` | File path; `output-debug.txt` (default) | Optional; File path |

See the [detailed parameter entries](../docs-site/app/docs/parameter-reference.ts).

### Station output

| Variable | Options / accepted values | Requirement / units |
| --- | --- | --- |
| `number_output_stations` | `0` (default); >= 0 | Optional; Integer |
| `station_binary` | `no` (default), `yes` | Optional; Yes/no |
| `station_ascii_buffer_kib` | `64` (default); 1..2097151 (INT_MAX/1024) | Optional; Integer KiB (1024 bytes) |
| `station_output_binary_file` | File path (1 to 240 bytes); Required when station_binary=yes; ignored otherwise; no default | Conditional; File path (1 to 240 bytes) |
| `output_stations_print_rate` | Integer time steps; Required when number_output_stations>0; > 0; no default | Conditional; Integer time steps |
| `output_stations_directory` | Directory path; `outputfiles/stations` (default) | Optional; Directory path |
| `output_stations` | N rows of `latitude longitude depth`; degrees, degrees, m; required when station count > 0; no default | Conditional; Table: N rows × 3 floats |
| `station_mediums_output_file` | File path; `output_stations_directory/station_mediums.txt` (default) | Optional; File path |
| `print_station_velocities` | `no` (default), `yes` | Optional; Yes/no |
| `print_station_accelerations` | `no` (default), `yes` | Optional; Yes/no |
| `print_station_pressures` | `no` (default), `yes` | Optional; Yes/no |
| `print_station_pore_pressure_ratios` | `no` (default), `yes`; requires include_poroelastic=yes and include_nonlinear_analysis=yes | Optional; Yes/no |
| `print_station_strains` | `no` (default), `yes` | Optional; Yes/no |
| `print_station_stresses` | `no` (default), `yes` | Optional; Yes/no |

Pressure eligibility and `station_mediums.txt` are finalized after solver
setup. A poroelastic receiver must lie in saturated, active BrickUP material;
elastic-only CPML receivers have undefined pressure/`Ru` (`nan`), even when
adjacent active soil shares their nodes. Skeleton motion is unchanged.

### Plane output

| Variable | Options / accepted values | Requirement / units |
| --- | --- | --- |
| `number_output_planes` | `0` (default); >= 0 | Optional; Integer |
| `planes_input_file` | File path; current `parameters.in` (default) | Optional; File path |
| `output_planes_print_rate` | Integer time steps; Required when number_output_planes>0; no default | Conditional; Integer time steps |
| `output_planes_directory` | Directory path; `outputfiles/planes` (default) | Optional; Directory path |
| `output_planes` | N rows of `lat lon depth d_strike n_strike d_dip n_dip strike dip`; degrees and m, integer counts; no default | Conditional; Table: N rows × 9 values |
| `reconstruct_acoustic_plane_displacements` | `yes` (default for scalar acoustics), `no` (default otherwise) | Optional; Yes/no; explicit `no` disables reconstruction |

Plane points exactly on mesh-domain faces, edges, and corners use the adjacent
interior element without moving their interpolation coordinates. Setup rejects
any plane point without exactly one owner, for every output backend, before
any displacement file is created; planes are never trimmed automatically. See
[plane boundary ownership](Extended%20README.md#plane-points-on-domain-faces).

Acoustic station/plane output inside CPML includes the containing element's
gradient stretch and point memory; use interior samples for physical motion.
Acoustic reconstruction uses physical pressure in both formulations; for
velocity potential it first computes `p = -(phi[n] - phi[n-1])/dt`.
See the [correction and validation note](Extended%20README.md#scalar-acoustic-corrections-and-limitations)
for plane output from older executables.

### Checkpoints and selected mesh coordinates

| Variable | Options / accepted values | Requirement / units |
| --- | --- | --- |
| `use_checkpoint` | `0` (default), `1` | Optional; Boolean integer |
| `checkpointing_rate` | `0` (default); >= 0 | Optional; Integer time steps |
| `checkpoint_path` | Directory path; `outputfiles/checkpoints` (default) | Optional; Directory path |
| `mesh_coordinates_for_matlab` | `no` (default), `yes` | Optional; Yes/no |
| `mesh_coordinates_directory_for_matlab` | Directory path; `outputfiles/For_Matlab` (default) | Optional; Directory path |
| `mesh_corners_matlab` | min_x min_y max_x max_y min_z max_z; no default | Conditional; List of 6 floats (m) |

See the [detailed parameter entries](../docs-site/app/docs/parameter-reference.ts).

### Schedule and initialization diagnostics

| Variable | Options / accepted values | Requirement / units |
| --- | --- | --- |
| `print_matrix_k` | `no` (default), `yes` | Optional; Yes/no |
| `enable_timing_barriers` | `0` (default), `1` | Optional; Boolean integer |
| `forces_buffer_size` | `104857600` (default) | Optional; Integer bytes |
| `schedule_print_error_check` | `0` (default), `1` | Optional; Boolean integer |
| `schedule_print_stdout` | `0` (default), `1` | Optional; Boolean integer |
| `schedule_print_file` | `0` (default), `1` | Optional; Boolean integer |
| `stat_schedule_filename` | File path; `stat-sched.txt` (default) | Optional; File path |
| `stat_mesh_filename` | File path; `stat-mesh.txt` (default) | Optional; File path |
| `stat_topo_filename` | File path; `stat-topo.txt` (default) | Optional; File path |

See the [detailed parameter entries](../docs-site/app/docs/parameter-reference.ts).

### Topography

| Variable | Options / accepted values | Requirement / units |
| --- | --- | --- |
| `include_topography` | `no` (default), `yes` | Optional; Yes/no |
| `topography_directory` | Directory path; `inputfiles/topography` (default) | Optional; Directory path |
| `maximum_octant_level` | Integer 0 to 30; no default | Required with topography; forced surface target in `maximum_level` mode, geometry bookkeeping in both modes; not a global cap |
| `topography_surface_refinement` | `maximum_level` (default), `wavelength`; case-insensitive | Optional; read only with topography enabled; skip forced surface sizing and its consistency check with `wavelength`; ordinary Vs/Vp wavelength rules and 2:1 balancing remain |
| `computation_method` | `Required with topography`; no default; vt \| fem | Conditional; String |
| `topographybase_zcoord` | Float (m); Required with topography; > 0; no default | Conditional; Float (m) |
| `type_of_etree` | `Required with topography`; no default; sqd \| full; no default | Conditional; String |
| `min_enclosed_volume_perc` | `0.1` (default) | Optional; Float |
| `consider_nonlinear_topography` | `no` (default), `yes` | Optional; Yes/no |
| `topography_bkt_gpu_stabilize` | `auto` (default), `yes`, `no` | Optional; String |

Mesh air/solid classification and receiver/material queries share a bilinear
surface, `Z = topographybase_zcoord - height`. Regenerate meshes made with
the earlier triangular classifier before rerunning affected surface outputs;
see [surface sampling](Extended%20README.md#topographic-surface-geometry-and-sampling).
Use `topography_surface_refinement=wavelength` with `mesh_input_mode=0` to
regenerate a mesh without the forced surface target. Retain the required
`maximum_octant_level`. Coarser surface cells may reduce terrain accuracy;
see [surface refinement](Extended%20README.md#topography-surface-refinement)
for water sizing, coupling restrictions and validation limits.

### Buildings and push-down regions

| Variable | Options / accepted values | Requirement / units |
| --- | --- | --- |
| `include_buildings` | `no` (default), `yes` | Optional; Yes/no |
| `number_of_buildings` | Integer; Required when include_buildings=yes; >= 0; no default | Conditional; Integer |
| `number_of_pushdowns` | Integer; Required when include_buildings=yes; >= 0; no default | Conditional; Integer |
| `min_octant_size_m` | Float (m); Required when include_buildings=yes; >= 0; no default | Conditional; Float (m) |
| `surface_shift_m` | Float (m); Required when include_buildings=yes; >= 0; no default | Conditional; Float (m) |
| `consider_fixed_base` | `no`, `yes`; no default; Required when include_buildings=yes | Conditional; Yes/no |
| `fixedbase_input_dt` | Float (s); Required when consider_fixed_base=yes; no default | Conditional; Float (s) |
| `fixedbase_input_dir` | Directory path; Required when consider_fixed_base=yes; no default | Conditional; Directory path |
| `fixedbase_input_startindex` | Integer; Required when consider_fixed_base=yes; no default | Conditional; Integer |
| `fixedbase_input_sufix` | Required when consider_fixed_base=yes; no default | Conditional; String |
| `building_properties` | xmin xmax ymin ymax depth height Vp_b Vs_b rho_b Vp_f Vs_f rho_f; no default | Conditional; Table: N rows × 12 floats |
| `pushdown_properties` | xmin xmax ymin ymax depth height; no default | Conditional; Table: N rows × 6 floats |

See the [detailed parameter entries](../docs-site/app/docs/parameter-reference.ts).

### Nonlinear and plastic controls

| Variable | Options / accepted values | Requirement / units |
| --- | --- | --- |
| `include_nonlinear_analysis` | `no` (default), `yes` | Optional; Yes/no |
| `geostatic_loading_time_sec` | >= 0; required for nonlinear analysis; no default | Conditional; Float (s) |
| `geostatic_cushion_time_sec` | >= 0; required for nonlinear analysis; no default | Conditional; Float (s) |
| `approximate_geostatic_state` | `no`, `yes`; no default; Required for nonlinear analysis | Conditional; Yes/no |
| `error_tolerance` | > 0; required for nonlinear analysis; no default | Conditional; Float |
| `backbone_errTol` | > 0; required for nonlinear analysis; no default | Conditional; Float |
| `material_model` | `linear`, `vonMises_ep`, `vonMises_fa`, `vonMises_faM`, `vonMises_baE`, `vonMises_baH`, `vonMises_GQH`, `vonMises_MKZ`, `vonMises_RO`, `MohrCoulomb`, `DruckerPrager`, `ManzariDafalias`, `PressureDependMultiYield03` (alias `PDMY03`), `RIVASand`, `RIVASandCustom`, `RIVASand02`, `RIVASandCustom02`; no default | Conditional; String |
| `material_plasticity_type` | `rate_dependent`, `rate_independent`; no default | Conditional; String |
| `material_properties_count` | >= 1; no default | Conditional; Integer |
| `material_properties_list` | Most models: N × 16; ManzariDafalias: N × 17; PDMY03: N × 22; RIVASand: N × 2 or 3; RIVASandCustom: N × 13; RIVASand02: N × 2 or 3; RIVASandCustom02: N × 13; no default | Conditional; Model-specific table |
| `no_substeps` | >= 1; required for nonlinear analysis; no default | Conditional; Integer |
| `rivasand02_reversaltype` | `1`, `2`, `3` (default); RIVASand02 and RIVASandCustom02 only | Optional; Integer (dimensionless) |
| `rivasand02_biasvolume` | `0`, `1` (default), `2`; RIVASand02 and RIVASandCustom02 only | Optional; Integer (dimensionless) |
| `manzari_dafalias_integration_scheme` | `modified_euler` (default), `forward_euler`; ManzariDafalias only | Optional; String |
| `manzari_dafalias_forward_euler_substeps` | `fixed` (default), `demand`; ManzariDafalias forward_euler only | Optional; String |
| `manzari_dafalias_recovery` | `strict` (default), `force_accept`, `unlimited_substeps`, `opensees`; ManzariDafalias modified_euler only | Optional; String |
| `stiff_damp` | Float; Required for nonlinear analysis; no default | Conditional; Float |
| `freq_stiff_damp` | Positive frequency; required for nonlinear analysis; no default | Conditional; Float (Hz) |
| `tension_cutoff` | `no`, `yes`; no default; Required for nonlinear analysis | Conditional; Yes/no |

RS02 reversal types are 1 (previous increment), 2 (accumulated branch), and
3 (default), the unguarded host stress-ratio predictor. BiasVolume modes are 0
(original bias volume), 1 (default), disabling bias volume but retaining phase
volume, and 2 (gated mean correction). Restart requires both selectors
unchanged; all older 138-value RS02 checkpoints must be rerun. RS02 GPU failures save a diagnostic snapshot
as `rivasand02_failure.rank<R>.step<S>.txt` in the run directory. See [RS02 details](Extended%20README.md#nonlinear-gpu-execution-and-constitutive-models).

### Standalone material-point test controls

These keys belong only in the input to `psolve --nonlinear-material-point`.
They are not production simulation controls. The driver also requires the
shared `material_model`, `no_substeps`, `error_tolerance`, `backbone_errTol`,
and `tension_cutoff` controls listed above.

| Variable | Options / accepted values | Requirement / units |
| --- | --- | --- |
| `material_point_delta_time_sec` | > 0; material-point subcommand only; no default | Conditional; Float (s) |
| `material_point_steps` | > 0; material-point subcommand only; no default | Conditional; Integer |
| `material_point_constants` | 24 floats normally; 17 for ManzariDafalias; 21 for PDMY03; 1 or 2 for standard RIVA-Sand; 12 for Custom names; no default | Conditional; Model-specific list |
| `material_point_initial_strain` | 6 values in `XX YY ZZ XY YZ XZ` order; omitted (default) for unseeded eligible models; all five arrays required if any is supplied; RIVA-Sand requires a seed | Conditional; 6 floats |
| `material_point_initial_stress` | 6 values in `XX YY ZZ XY YZ XZ` order; omitted (default) for unseeded eligible models; all five arrays required if any is supplied; RIVA-Sand requires a seed | Conditional; 6 floats (Pa) |
| `material_point_initial_alpha` | 6 values in `XX YY ZZ XY YZ XZ` order; omitted (default) for unseeded eligible models; all five arrays required if any is supplied; RIVA-Sand requires a seed | Conditional; 6 floats |
| `material_point_initial_alpha_in` | 6 values in `XX YY ZZ XY YZ XZ` order; omitted (default) for unseeded eligible models; all five arrays required if any is supplied; RIVA-Sand requires a seed | Conditional; 6 floats |
| `material_point_initial_fabric` | 6 values in `XX YY ZZ XY YZ XZ` order; omitted (default) for unseeded eligible models; all five arrays required if any is supplied; RIVA-Sand requires a seed | Conditional; 6 floats |
| `material_point_strain_history` | epsilon_xx epsilon_yy epsilon_zz epsilon_xy epsilon_yz epsilon_xz; no default | Conditional; Table: N rows × 6 floats |
| `material_point_riva_density` | > 0; required for automatic RIVA-Sand material-point stiffness; no default | Conditional; Float (kg/m³) |
| `material_point_riva_vs` | > 0; required for automatic RIVA-Sand material-point stiffness; no default | Conditional; Float (m/s) |
| `material_point_compare_rtol` | `1e-6` (default); >= 0; material-point subcommand only | Optional; Float |
| `material_point_compare_atol` | `1e-9` (default); >= 0; material-point subcommand only | Optional; Float |
| `material_point_output_prefix` | Path prefix; `material-point` (default); material-point subcommand only | Optional; Path prefix |
| `material_point_brickup` | `no` (default), `yes`; material-point subcommand only | Optional; Yes/no |
| `material_point_brickup_porosity` | 0 < n <= 1; required when material_point_brickup=yes; no default | Conditional; Float |
| `material_point_brickup_fluid_bulk_modulus_pa` | > 0; required when material_point_brickup=yes; no default | Conditional; Float (Pa) |
| `material_point_brickup_biot_alpha` | 0 < alpha <= 1; required when material_point_brickup=yes; no default | Conditional; Float |
| `material_point_brickup_solid_bulk_modulus_pa` | > 0 when material_point_brickup_biot_alpha<1; material-point subcommand only; no default | Conditional; Float (Pa) |
| `material_point_initial_pore_pressure_pa` | Finite number; required when `material_point_brickup=yes`; no default | Conditional; Float (Pa) |

See the [detailed parameter entries](../docs-site/app/docs/parameter-reference.ts).

### Domain Reduction Method

| Variable | Options / accepted values | Requirement / units |
| --- | --- | --- |
| `implement_drm` | `no` (default), `yes` | Optional; Yes/no |
| `number_output_drm_boxes` | `0` (default) | Optional; Nonnegative integer |
| `output_drm_boxes_print_rate` | Required when number_output_drm_boxes > 0; no default | Conditional; Positive integer time steps |
| `output_drm_boxes_directory` | Directory path; Required when number_output_drm_boxes > 0; no default | Conditional; Directory path |
| `output_drm_boxes` | N rows of `top_center_x top_center_y length_x length_y depth spacing strike`; m and degrees; no default | Conditional; N rows × 7 floats |
| `drm_output_batch_size` | `16` (default) | Optional; Positive integer snapshots |
| `drm_directory` | Directory path; Required by the DRM module; no default | Conditional; Directory path |
| `which_drm_part` | `part0`, `part1`, `part2`; no default | Conditional; String |
| `drm_edgesize` | Float (m); Required when implement_drm=yes; no default | Conditional; Float (m) |
| `drm_offset_x` | Float (m); Required when implement_drm=yes; no default | Conditional; Float (m) |
| `drm_offset_y` | Float (m); Required when implement_drm=yes; no default | Conditional; Float (m) |
| `drm_print_rate` | Integer time steps; Required when implement_drm=yes; no default | Conditional; Integer time steps |
| `part1_delta_t` | Float (s); Required when implement_drm=yes; no default | Conditional; Float (s) |
| `part1_simulation_time` | Float (s); Required when implement_drm=yes; no default | Conditional; Float (s) |
| `drm_boundary` | Required when implement_drm=yes; no default | Conditional; List of 5 floats (m) |

See the [detailed parameter entries](../docs-site/app/docs/parameter-reference.ts).

### Incident plane waves

| Variable | Options / accepted values | Requirement / units |
| --- | --- | --- |
| `include_incident_planewaves` | `no` (default), `yes` | Optional; Yes/no |
| `type_of_wave` | `SV`, `P`; no default | Conditional; String |
| `DRMBox_Noelem_Halfwidth_EW` | Float (elements); Required for incident waves; no default | Conditional; Float (elements) |
| `DRMBox_Noelem_Halfwidth_NS` | Float (elements); Required for incident waves; no default | Conditional; Float (elements) |
| `DRMBox_Noelem_depth` | Float (elements); Required for incident waves; no default | Conditional; Float (elements) |
| `DRMBox_element_size` | Float (m); Required for incident waves; no default | Conditional; Float (m) |
| `DRM_xo` | Float (m); Required for incident waves; no default | Conditional; Float (m) |
| `DRM_yo` | Float (m); Required for incident waves; no default | Conditional; Float (m) |
| `fnc_type` | `RICKER`, `time_hist`; no default | Conditional; String |
| `planewave_strike` | Float (degrees); Required for incident waves; no default | Conditional; Float (degrees) |
| `planewave_Z_angle` | Float (degrees); Required for incident waves; no default | Conditional; Float (degrees) |
| `Ts` | Float (s); Required for incident waves; no default | Conditional; Float (s) |
| `fc` | Float (Hz); Required for incident waves; no default | Conditional; Float (Hz) |
| `Uo` | Float; Required for incident waves; no default | Conditional; Float |
| `no_components` | Normally 1 or 2; no default | Conditional; Integer |
| `ug_timestep` | Float (s); Required for incident waves; no default | Conditional; Float (s) |
| `ug_alongstrike` | File path; Required when fnc_type=time_hist; no default | Conditional; File path |
| `ug_alongnormal` | File path; Required when fnc_type=time_hist and no_components=2; no default | Conditional; File path |

See the [detailed parameter entries](../docs-site/app/docs/parameter-reference.ts).

### Homogeneous half-space

| Variable | Options / accepted values | Requirement / units |
| --- | --- | --- |
| `include_hmg_halfspace` | `no` (default), `yes` | Optional; Yes/no |
| `etreeBox_xo` | Float (m); Required when include_hmg_halfspace=yes; no default | Conditional; Float (m) |
| `etreeBox_yo` | Float (m); Required when include_hmg_halfspace=yes; no default | Conditional; Float (m) |
| `etreeBox_lengthEW` | Float (m); Required when include_hmg_halfspace=yes; no default | Conditional; Float (m) |
| `etreeBox_lengthNS` | Float (m); Required when include_hmg_halfspace=yes; no default | Conditional; Float (m) |
| `etreeBox_depth` | Float (m); Required when include_hmg_halfspace=yes; no default | Conditional; Float (m) |
| `hmgHS_Vs` | Float (m/s); Required when include_hmg_halfspace=yes; no default | Conditional; Float (m/s) |
| `hmgHS_Vp` | Float (m/s); Required when include_hmg_halfspace=yes; no default | Conditional; Float (m/s) |
| `hmgHS_rho` | Float (kg/m3); Required when include_hmg_halfspace=yes; no default | Conditional; Float (kg/m3) |

See the [detailed parameter entries](../docs-site/app/docs/parameter-reference.ts).

### Source controls

These controls belong in `source.in`. The table covers the force-source and
filtering controls discussed here; other source-specific geometry remains in
the detailed source documentation.

| Variable | Options / accepted values | Requirement / units |
| --- | --- | --- |
| `type_of_source` | `point`, `plane`, `planewithkinks`, `srfh`, `force`, `base_force`; no default | Required source family |
| `source_function_type` | `ramp`, `sine`, `quadratic`, `ricker`, `exponential`, `discrete`; no default | Required; `base_force` requires `discrete` |
| `srfh_depth_reference` | `surface` (default), `sea_level`; case-insensitive | Optional in `source.in`, SRFH only; interprets `coords.in` depths in m positive down. With topography: local outer surface or zero-elevation/sea-level datum. Independent of the material depth reference; see [SRFH source depths](Extended%20README.md#srfh-source-depth-reference). |
| `average_risetime_sec` | Positive duration; no default | Required for ramp, sine, quadratic and exponential; s |
| `ricker_Ts` | Time shift; no default | Required for Ricker; s |
| `ricker_Tp` | Positive characteristic period; no default | Required for Ricker; s |
| `lonlat_or_cartesian` | `0` (geographic), `1` (Cartesian); no default | Required for point/force location |
| `hypocenter_lat_deg` | Valid latitude; no default | Required with geographic point/force location; degrees |
| `hypocenter_long_deg` | Valid longitude; no default | Required with geographic point/force location; degrees |
| `hypocenter_x` | Solver X coordinate; no default | Required with Cartesian point/force location; m |
| `hypocenter_y` | Solver Y coordinate; no default | Required with Cartesian point/force location; m |
| `hypocenter_depth_m` | Depth positive downward; no default | Required for point/force source; m |
| `force_x_n` | Signed force component; no default | Required with `type_of_source=force`; N, solver X |
| `force_y_n` | Signed force component; no default | Required with `type_of_source=force`; N, solver Y |
| `force_z_n` | Signed force component; no default | Required with `type_of_source=force`; N, positive down |
| `base_force_rho_kgm3` | Positive density; no default | Required with `base_force`; kg/m³ |
| `base_force_vs_mps` | Positive shear-wave speed; no default | Required with `base_force`; m/s |
| `base_force_direction_x` | Signed component; no default | Required with `base_force`; combined vector must be nonzero and is normalized |
| `base_force_direction_y` | Signed component; no default | Required with `base_force`; dimensionless |
| `base_force_direction_z` | Signed component; no default | Required with `base_force`; dimensionless, positive down |
| `base_force_scale` | Nonnegative number; `1.0` (default) for incident/within history; use `2.0` for outcrop history | Optional; dimensionless |
| `base_force_depth_m` | Deepest mesh-node plane (default) when omitted or negative; nonnegative model depth selects a plane | Optional; m |
| `base_force_tolerance_m` | `max(1e-6 m, 0.01 × minimum mesh edge)` (default) when omitted or `-1`; positive number overrides | Optional; m |
| `source_is_filtered` | `0` (off), `1` (on); no default | Required; filtering applies to plane, planewithkinks and srfh |
| `threshold_frequency` | Positive cutoff frequency; no default | Required when filtering is enabled; Hz |
| `number_of_poles` | Positive integer Butterworth order; no default | Required when filtering is enabled |

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

These are environment variables, except `USE_OPENMP` (a Make variable) and
`HERCULES_ROOT` (an example shell variable). `(default)` marks omission
behavior; launcher, direct-solver and conditional defaults are distinguished.
External OpenMP/MPI/UCX values must be supported by the installed runtime.
The tables give canonical choices and the documented numeric aliases; see
[runtime controls](Extended%20README.md#runtime-performance-controls),
[diagnostics](Extended%20README.md#advanced-and-diagnostic-controls),
[launcher details](Extended%20README.md#launch-a-simulation), and
[nonlinear execution](Extended%20README.md#nonlinear-gpu-execution-and-constitutive-models)
for interactions and performance implications.

### Launcher, build and CPU placement

| Variable | Options / accepted values | Use / conditions |
| --- | --- | --- |
| `USE_OPENMP` | `1` (default), `0` | Make variable: enable/disable OpenMP |
| `HERCULES_ROOT` | Repository directory; no default | Convenience variable used by the launch example, not read by the solver |
| `BASE_IN` | Path to the complete `parameters.in`; no default | Required for unified mode; optional in two-phase mode when separate inputs are supplied |
| `PSOLVE_BIN` | Executable path; auto-detect `quake/forward/psolve`, then `psolve` beside the launcher (default) | Override when the binary is built elsewhere |
| `MPI_RUN` | MPI launcher command/path; `mpirun` (default) | Launcher executable |
| `PHASE_A_IN` | Input-file path; case `inputfiles/Test_A.in`, otherwise `Input/Test_A.in` under the launcher parent (default) | Two-phase input when `BASE_IN` is omitted |
| `PHASE_B_IN` | Input-file path; case `inputfiles/Test_B.in`, otherwise `Input/Test_B.in` under the launcher parent (default) | Two-phase input when `BASE_IN` is omitted |
| `PHASE_A_MPI_ARGS` | Space-separated MPI arguments; empty (default) | Additional two-phase meshing arguments |
| `PHASE_B_MPI_ARGS` | Space-separated MPI arguments; empty (default) | Additional two-phase solver arguments |
| `UNIFIED_MPI_ARGS` | Space-separated MPI arguments; `--bind-to none` (default) | Unified launch manages CPU affinity during handoff |
| `HERCULES_PHASE_MODE` | `unified` (default), `two-phase` | One MPI job with mesh handoff, or separate jobs |
| `HERCULES_UNIFIED_IDLE_POLL_USEC` | Integer 100..1000000; `10000` (default) | Inactive-rank polling interval, microseconds |
| `HERCULES_UNIFIED_MPI_YIELD_WHEN_IDLE` | `0` (default), `1` | Unified Open MPI progress policy |
| `PHASE_A_NP` | Positive integer; no default in unified mode; MPI-launcher rank selection (default) when omitted in two-phase mode | Meshing ranks; positional argument 1 overrides |
| `PHASE_B_NP` | Positive integer; allocated node count (default) if detected, otherwise `1` (default) | Solver ranks; positional argument 2 overrides; multi-node oversubscription must divide evenly |
| `PHASE_A_OMP_NUM_THREADS` | Positive integer; `1` (default) | Meshing threads per rank |
| `PHASE_B_OMP_NUM_THREADS` | Positive integer; inherit `OMP_NUM_THREADS`, otherwise derive CPUs per solver rank from the allocation (default) | Unified mode requires an explicit value if inference fails; two-phase can leave selection to the OpenMP runtime |
| `OMP_NUM_THREADS` | Positive thread count (or OpenMP nesting list for direct runs); runtime-selected when unset (default) | Phase-B fallback; launcher thread controls require a single positive integer |
| `OMP_PLACES` | `cores` (default) in launcher, `threads`, `sockets`, or a runtime-supported explicit place list | Direct runs use their OpenMP runtime default |
| `OMP_PROC_BIND` | `close` (default) in two-phase launcher, `spread`, `primary`, `true`, `false`, or a supported nesting list | Unified launcher forces `false` at startup, then assigns Phase-B affinity internally |
| `PHASE_B_MAP_MODE` | `auto` (default), `per-node`, `single-node`, `numa` | Two-phase placement; NUMA mode requires compatible ranks/node and NUMA domains |
| `PHASE_A_UCX_TLS` | UCX transport list; `sm,self,tcp` (default) with CUDA-aware Open MPI, otherwise inherited MPI selection | CPU-meshing transport override for two-phase mode |
| `PHASE_B_UCX_TLS` | UCX transport list; `dc,sm,self,cuda_copy,cuda_ipc` (default) with CUDA-aware Open MPI, otherwise inherited MPI selection | Also applies from startup in unified mode; omit `cuda_ipc` only for a deliberate site override |
| `PHASE_B_OMPI_MCA_PML` | PML supported by the MPI installation; `ucx` (default) with CUDA-aware Open MPI, otherwise inherited MPI selection | Phase-B/unified PML override |
| `PHASE_B_KEEP_OMPI_MCA_BTL` | `0` (default), `1` | Keep an inherited BTL configuration only when explicitly needed |
| `UCX_TLS` | UCX-supported comma-separated transport list; omitted (default), using phase-specific selection above | A global value overrides both phase-specific transport settings |
| `OMPI_MCA_pml` | PML supported by the MPI installation; omitted (default), using phase-specific selection above | A global value overrides `PHASE_B_OMPI_MCA_PML` |
| `OMPI_MCA_btl` | BTLs supported by the MPI installation; omitted (default) | CUDA-aware Phase-B/unified launch unsets an inherited value unless the keep flag is `1` |
| `IO_PES` | Nonnegative integer; `0` (default) | Dedicated plane-I/O ranks; positive values select the legacy I/O-pool backend; unified mode requires `0` |

### Mesh import and partition diagnostics

| Variable | Options / accepted values | Use / conditions |
| --- | --- | --- |
| `HERCULES_MESH_IMPORT_TREE_MODE` | `cursor` (default), `refine` | Saved-leaf cursor import or recursive refinement |
| `HERCULES_MESH_IMPORT_VERIFY` | `0` (default), `1` | Extra import-refinement verification; effective only in `refine` mode |
| `HERCULES_MESH_IO_MODE` | `mpiio` (default) with unified launcher, `sharded` (default) with two-phase launcher, `single` (default) for direct `psolve` | Storage format; explicit setting overrides the launcher choice |
| `OCTOR_EXTRACT_PROFILE` | `0` (default) for ordinary meshing, `1` (default) for imported/unified-reconstructed meshes, `2` | Off, stage timings, or detailed timings |
| `OCTOR_IMPORT_PATH_CACHE` | `1` (default), `0` (aliases `off`, `false`, `no`) | Reuse shared paths during cursor import |
| `OCTOR_INIT_TASKS_PER_RANK` | Positive integer; `10` (default); values above 8192 are clamped | Coarse bootstrap tasks; larger values require mesh validation |
| `OCTOR_NODE_SORT_MODE` | `radix` (default), `qsort` (alias `legacy`) | Morton-order node sorting |
| `OCTOR_SETCOM_ALL2ALL` | `0` (default) for mesh generation, `1` (default) for imported/unified-reconstructed meshes | Two-phase Phase A explicitly uses `0`; an explicit override applies throughout a unified job |

### Station and plane output

| Variable | Options / accepted values | Use / conditions |
| --- | --- | --- |
| `HERCULES_PLANE_IO_MODE` | `mpiio` (default), `hybrid`, `legacy` | GPU/MPI-IO, GPU/root writer, or CPU/root writer; `IO_PES>0` forces legacy pool |
| `HERCULES_PLANE_IO_OMP` | `1` (default), `0` | Host plane setup/reconstruction; does not control GPU interpolation |
| `HERCULES_PLANE_MPIIO_BATCH` | Integer 1..2147483647; `64` (default) | Snapshots buffered per collective plane write |
| `HERCULES_PLANE_METADATA` | `0` (default), `1` | Write material/CPML point audit once at setup |
| `HERCULES_STATION_IO_OMP` | `1` (default), `0` | Host station work; does not select CUDA gathers |
| `HERCULES_STATION_ROUTE_MODE` | `interval` (default), `all` | Ownership lookup restricted by partition intervals or on every rank |

### Nonlinear execution

| Variable | Options / accepted values | Use / conditions |
| --- | --- | --- |
| `HERCULES_NONLINEAR_BACKEND` | `gpu` (default), `cpu` | GPU constitutive updates or CPU reference; some material families require GPU at element scale |
| `HERCULES_NONLINEAR_STATS_INTERVAL` | Nonnegative integer; `0` (default) | Positive N samples Manzari-Dafalias workload every N steps |

### Field synchronization and MPI

| Variable | Options / accepted values | Use / conditions |
| --- | --- | --- |
| `HERCULES_CONTRIB_PACK_OMP` | `1` (default), `0` | OpenMP contribution-message packing |
| `HERCULES_CONTRIB_UNPACK_OMP` | `1` (default), `0` | OpenMP unpacking when destination mappings have no duplicates |
| `HERCULES_CONTRIB_RECV_MODE` | `waitsome` (default), `waitall` | Reproducibility modes force deterministic `waitall` |
| `HERCULES_DISP_SYNC_MODE` | `auto` (default), `full` / `0`, `boundary` / `1` | Unsafe boundary-only configurations fall back to full host visibility |
| `HERCULES_GPU_MPI_BACKEND` | `cuda_aware` (default) when MPI reports CUDA support, `host` (default) otherwise | Checkpoint/host-dependent paths may force host schedules |
| `HERCULES_SCALAR_MPI_SCHEDULE` | `compact` (default), `full` | Setup-built acoustic, pore-pressure and open-pore-target schedules on both backends; retains ownership, field constraints and output ghosts. `full` restores the general schedules for comparison. Same value on all ranks; [details](Extended%20README.md#scalar-mpi-schedules-and-profiling) |
| `HERCULES_SCALAR_MPI_BATCH` | `1` (default), `0` | Combine sealed acoustic/pore RHS and field exchanges on matching backends (12 to 8 total phases); combine open-pore RHS/target contributions (14 to 12). Host and CUDA-aware; preserves field maps and time levels. `0` keeps fields separate. Same value on all ranks; [details](Extended%20README.md#scalar-mpi-schedules-and-profiling) |
| `HERCULES_GPU_MPI_OVERLAP` | `1` (default), `0` | CUDA-aware scalar interior updates overlap anchored RHS messages; scalar/vector unpack uses receive-reuse events. Host transport retains synchronous staging. `0` uses whole-field updates and synchronous unpack; set batching to `0` too for the previous sequential path. Same value on all ranks. |
| `HERCULES_GPU_MPI_PROFILE` | `0` (default), `1` | Device-helper and scalar-batch counters/CPU times, including combined packets, reported as cumulative rank average/min/max. Adds no synchronization; overlapping unpack measures enqueue time. Same value on all ranks; [details](Extended%20README.md#scalar-mpi-schedules-and-profiling) |
| `HERCULES_PORE_PRESSURE_UPDATE` | `output` (default), `step` | Eligible fresh linear undrained fast-path runs reconstruct pore pressure for pressure station output and the final two states; other steps skip pore kernels/exchanges. Checkpoint/restart, nonlinear, diffusion, drainage and open-pore retain every-step updates. Ru output requires nonlinear analysis and therefore retains every-step updates. Same value on all ranks. |
| `HERCULES_PORE_PRESSURE_SYNC_BACKEND` | Omitted (default); `host`, `cuda_aware` | Legacy alias, used only when `HERCULES_GPU_MPI_BACKEND` is unset |
| `HERCULES_PRESSURE_FORCE_SYNC_MODE` | `sparse` / `1` (default), `full` / `0` | Acoustic scalar force exchange |
| `HERCULES_PRESSURE_FIELD_SYNC_MODE` | `resident` / `1` (default), `full` / `0` | Checkpoint restart forces full host-visible acoustic state |
| `HERCULES_PORE_PRESSURE_FORCE_SYNC_MODE` | `sparse` / `1` (default), `full` / `0` | Poroelastic scalar force exchange |
| `HERCULES_PORE_PRESSURE_FIELD_SYNC_MODE` | `resident` / `1` (default), `full` / `0` | Poroelastic field exchange and station gathers |

### Other performance and diagnostic controls

| Variable | Options / accepted values | Use / conditions |
| --- | --- | --- |
| `HERCULES_ADJUST_ASSIGN_OMP` | `1` (default), `0` | OpenMP assignment adjustment |
| `HERCULES_ADJUST_DIST_OMP` | `1` (default), `0` | OpenMP distribution adjustment |
| `HERCULES_ALLOW_MULTI_RANK_PER_GPU` | `0` (default), `1` | Allow intentional GPU oversubscription |
| `HERCULES_FORCE_GENERIC_DAMPING` | `0` (default), `1` | Force generic BKT-family kernels for comparison |
| `HERCULES_BKT2_FUSED` | `1` (default), `0` | Use fused kernels only for eligible device-resident BKT2 runs; generic override takes precedence |
| `HERCULES_FULL_CLEAR_OMP` | `1` (default), `0` | OpenMP full-array host force clearing |
| `HERCULES_REPRO_MODE` | `0` (default), `1`, `2` | Normal, deterministic receive ordering, or stricter OpenMP restrictions |
| `HERCULES_SPARSE_CLEAR_MODE` | `full` / `3` (default), `auto` / `2`, `range` / `0`, `nodes` / `1` | Host topography-force clearing |
| `HERCULES_SPARSE_UPLOAD_MODE` | `dense` / `1` (default), `auto` / `0`, `sparse` / `2`, `range` / `3` | Topography-force upload path |
| `HERCULES_SRC_DEBUG` | `0` (default), `1` | Per-source application diagnostic logging |
| `HERCULES_TOPO_FORCE_MODE` | `gpu` / `1` (default), `cpu` / `0` | Eligible VT effective-force assembly; FEM has no such per-step work |

Host-staged execution retains the reference BKT2 path. Setup reports fused-path
fallback reasons and excludes inactive sectors from rate comparisons; see
[BKT2 GPU execution](Extended%20README.md#bkt2-gpu-execution).

## Output and validation

`station_binary`: `no` (default) for ASCII, `yes` (single binary file).
`station_ascii_buffer_kib` sets the CPU buffer size per locally owned ASCII
station: optional integer KiB (1 KiB = 1024 bytes), **64** (default), range
1..2097151. It is validated but unused in binary mode. All requested samples and
columns are retained. Buffers are written when full and flushed at normal
completion; live files can lag the simulation. The final flush is included in
the station and solver I/O timers. See
[ASCII station buffering](Extended%20README.md#ascii-station-buffering).
To record displacement from the existing `output_stations` list in one file:

```ini
station_binary = yes
station_output_binary_file = output/station_displacements.0
```

Use `number_output_stations` and `output_stations_print_rate` as before. The
binary is native float64 `[snapshot][station][XYZ]` in metres, in station-row
order. See [station binary details](Extended%20README.md#single-file-station-binary-output).
Plane files are native float64 displacement grids; formats and readers are in
the [output guide](https://hercules-suite.github.io/Hercules-GPU/docs/output-files/)
and [script index](../script/README.md).
`HERCULES_PLANE_METADATA=1` writes `plane_point_metadata.RANK.csv` once during
setup, with actual point coordinates, containing-element material and CPML
membership. See [plane metadata and interpretation](Extended%20README.md#diagnose-large-plane-displacements).
For acoustic observations, enable ASCII station velocity, acceleration and
pressure output and verify ownership in `station_mediums.txt`. See
[pressure and motion interpretation](Extended%20README.md#acoustic-pressure-and-motion-diagnostics)
and [gravity-wave recovery limits](Extended%20README.md#check-gravity-wave-dispersion-and-reconstruction).
Validate generated source forces before interpreting a failed boundary experiment.
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
spacing h, and clockwise strike from solver +X. Hercules writes one
`drmdisplacements.N` containing five exterior and five interior face grids;
interior dimensions are `(Lx-2h, Ly-2h, depth-h)`. Configure
`number_output_drm_boxes`, `output_drm_boxes_print_rate`, and
`output_drm_boxes_directory`. This sampling is independent of legacy DRM
replay and does not refine the mesh.

Each box also writes `drmcoordinates.N` once during setup: an ASCII table of
`X_m Y_m Z_m Vs_m_s Vp_m_s density_kg_m3 Qs Qp`. Coordinates are solver-local
X/Y, Z/down in metres, including terrain/building shifts. Appended properties
come from the initial owning mesh element: speeds in m/s, density in kg/m³,
and dimensionless selected nominal Q values. Qp includes lossless-bulk fallback,
infinite-Qk overrides, high-Q cutoffs and final Vs/Vp; it is not the input Qp.
`inf` denotes zero selected loss; `nan` denotes no applicable scalar Q.
They do not represent evolving nonlinear loss or effective poroelastic attenuation.
Rows match binary samples exactly, including repeated edges/corners. Use
`numpy.loadtxt(path, usecols=(0, 1, 2))` for positions only. Rebuild and regenerate
older meshes with `mesh_input_mode=0`; mesh payloads now retain Qs/Qp.

See [the complete geometry, layout, and reader examples](Extended%20README.md#structured-drm-box-output).
