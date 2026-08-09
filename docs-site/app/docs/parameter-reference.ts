import type { DocPageContent, DocSection } from "./content";

export type ParameterRequirement = "Required" | "Conditional" | "Optional" | "Legacy";

export type ParameterDefinition = {
  name: string;
  requirement: ParameterRequirement;
  type: string;
  values: string;
  description: string;
};

type ParameterSection = {
  id: string;
  title: string;
  introduction?: string;
  parameters: ParameterDefinition[];
};

type ParameterPage = {
  slug: string;
  title: string;
  summary: string;
  sections: ParameterSection[];
};

const p = (
  name: string,
  requirement: ParameterRequirement,
  type: string,
  values: string,
  description: string,
): ParameterDefinition => ({ name, requirement, type, values, description });

const parameterPages: ParameterPage[] = [
  {
    slug: "parameters-domain-time",
    title: "Domain, time, and resolution parameters",
    summary: "Geographic placement, physical dimensions, simulation duration, time step, and wavelength-based mesh resolution.",
    sections: [
      {
        id: "domain",
        title: "Domain geometry",
        parameters: [
          p("region_origin_latitude_deg", "Required", "Float (degrees)", "Any valid latitude", "Latitude of the lower-left surface corner of the model domain."),
          p("region_origin_longitude_deg", "Required", "Float (degrees)", "Any valid longitude", "Longitude of the lower-left surface corner of the model domain."),
          p("region_depth_shallow_m", "Required", "Float (m)", "Usually >= 0", "Depth of the shallow domain boundary, positive downward."),
          p("region_length_east_m", "Required", "Float (m)", "> 0", "Domain length in the local east direction; this becomes the solver Y dimension."),
          p("region_length_north_m", "Required", "Float (m)", "> 0", "Domain length in the local north direction; this becomes the solver X dimension."),
          p("region_depth_deep_m", "Required", "Float (m)", "> region_depth_shallow_m", "Depth of the bottom boundary. The solver depth is deep minus shallow depth."),
          p("region_azimuth_leftface_deg", "Required", "Float (degrees)", "Angular azimuth", "Azimuth of the domain left face; also establishes the local horizontal frame used by source calculations."),
          p("domain_surface_corners", "Conditional", "Table: 4 rows × 2 floats", "longitude latitude", "Required only when Hercules is built without PROJ. Give the four surface corners clockwise, starting at the origin."),
        ],
      },
      {
        id: "time-resolution",
        title: "Time integration and spatial resolution",
        introduction: "For elastic material, the nominal element-size target is h <= Vs / (nodes per wavelength × maximum frequency). Acoustic elements use Vp and may apply an additional refinement factor.",
        parameters: [
          p("simulation_start_time_sec", "Optional", "Float (s)", "Default: 0; >= 0", "Simulation start time. It must not exceed simulation_end_time_sec."),
          p("simulation_end_time_sec", "Required", "Float (s)", ">= simulation_start_time_sec", "Simulation end time."),
          p("simulation_delta_time_sec", "Required", "Float (s)", "> 0", "Explicit time-step interval. Select it from the fastest wave speed and smallest element in the final model."),
          p("simulation_wave_max_freq_hz", "Required", "Float (Hz)", "> 0", "Highest frequency the mesh is intended to resolve; used with nodes per wavelength to set octree refinement."),
          p("simulation_node_per_wavelength", "Required", "Integer", "> 0", "Target number of mesh nodes per shortest modeled wavelength."),
          p("simulation_shear_velocity_min", "Required", "Float (m/s)", "> 0", "Lower bound used for elastic shear velocity and wavelength-based refinement."),
        ],
      },
    ],
  },
  {
    slug: "parameters-material-mesh",
    title: "Material, source-path, and mesh parameters",
    summary: "Material databases and profiles, optional basin and 3-D models, source locations, and mesh generation or import modes.",
    sections: [
      {
        id: "base-material",
        title: "Base material model",
        parameters: [
          p("cvmdb_input_file", "Required", "Path or keyword", "File path or profile", "Material database used for the background crustal structure. Set to profile to read number_profile_layers and profile_layers from this file."),
          p("number_profile_layers", "Conditional", "Integer", "> 0 when cvmdb_input_file=profile", "Number of rows in profile_layers."),
          p("profile_layers", "Conditional", "Table: N rows × 6 floats", "z_top(m) Vp(m/s) Vs(m/s) rho(kg/m3) Qp Qs", "Layered half-space definition, required when cvmdb_input_file=profile. Rows begin at their top depth and should be ordered from shallow to deep."),
          p("the_threshold_Vp_over_Vs", "Required", "Float", ">= 0", "Upper limit applied to Vp/Vs while correcting element material properties and computing the first Lamé parameter."),
          p("softening_factor", "Optional", "Float", "Default: 0; either 0 or > 1", "Material-compliance adjustment tied to element size. Zero disables softening; nonzero values require case-specific validation."),
          p("stiffness_calculation_method", "Optional", "String", "Default: effective; effective | conventional", "Finite-element stiffness operator. effective is the current recommended and default method; conventional is retained for compatibility."),
          p("use_progressive_meshing", "Optional", "Integer", "Default: 1; >= 0", "Progressive-meshing factor/iteration control. Current production inputs normally retain the default."),
        ],
      },
      {
        id: "velocity-models",
        title: "Optional regional material models",
        parameters: [
          p("Istanbul_velocity_model", "Optional", "Yes/no", "Default: no", "Enable the built-in Istanbul material-model replacement where it covers the domain."),
          p("basin_velocity_model", "Optional", "Yes/no", "Default: no", "Enable a basin geometry that replaces the background material inside its boundary."),
          p("basin_origin_latitude_deg", "Conditional", "Float (degrees)", "Required when basin_velocity_model=yes", "Latitude of the basin model's lower-left origin."),
          p("basin_origin_longitude_deg", "Conditional", "Float (degrees)", "Required when basin_velocity_model=yes", "Longitude of the basin model's lower-left origin."),
          p("basin_input_file", "Conditional", "Path", "Required when basin_velocity_model=yes", "Text file describing the basin boundary grid."),
          p("basin_material_type", "Optional", "String", "Default: elastic; elastic | acoustic_water", "Material assigned inside basin geometry. acoustic_water requires include_acoustic=yes."),
          p("basin_water_vp", "Optional", "Float (m/s)", "Default: 1500; > 0", "Compressional-wave speed used when basin_material_type=acoustic_water."),
          p("basin_water_rho", "Optional", "Float (kg/m3)", "Default: 1000; > 0", "Density used when basin_material_type=acoustic_water."),
          p("3D_velocity_model", "Optional", "Yes/no", "Default: no", "Enable the structured general 3-D velocity model."),
          p("3D_velocity_model_directory", "Optional", "Directory path", "Default: inputfiles/materialfiles", "Directory containing x.in, y.in, z.in, vs.in, vp.in, rho.in, index.in, and optional matched qs.in/qp.in files."),
          p("3D_velocity_model_origin_latitude_deg", "Conditional", "Float (degrees)", "Required when 3D_velocity_model=yes", "Latitude of the structured model's lower-left origin."),
          p("3D_velocity_model_origin_longitude_deg", "Conditional", "Float (degrees)", "Required when 3D_velocity_model=yes", "Longitude of the structured model's lower-left origin."),
        ],
      },
      {
        id: "source-paths",
        title: "Source directories",
        parameters: [
          p("source_directory", "Conditional", "Directory path", "Required for source-driven runs", "Directory containing source.in and its associated histories or fault files. It is read for ordinary runs and DRM part1, but not for DRM part0/part2 source bypasses."),
          p("source_directory_output", "Optional", "Directory path", "Default: outputfiles/srctmp", "Directory for generated source-force files and source-related temporary output."),
        ],
      },
      {
        id: "mesh-modes",
        title: "Mesh generation, persistence, and import",
        parameters: [
          p("output_mesh", "Optional", "Boolean integer", "Default: 1; 0 | 1", "Write the generated mesh. The exact artifact is a legacy etree, sharded manifest, or MPI-IO archive according to the runtime mesh I/O mode."),
          p("mesh_only_mode", "Optional", "Boolean integer", "Default: 0; 0 | 1", "When 1, stop after mesh/setup work instead of entering time stepping."),
          p("mesh_input_mode", "Optional", "Boolean integer", "Default: 0; 0 | 1", "When 1, import an existing mesh rather than generating one from the material and geometry models."),
          p("mesh_etree_output_file", "Conditional", "File path", "Required for mesh output or mesh import", "Mesh artifact path. In sharded mode this is the manifest; in MPI-IO mode the archive; in single mode the etree database."),
        ],
      },
    ],
  },
  {
    slug: "parameters-damping-cpml",
    title: "Damping and absorbing-boundary parameters",
    summary: "Material attenuation, damping limits, the legacy dashpot boundary, and the current CPML controls.",
    sections: [
      {
        id: "material-damping",
        title: "Material damping",
        parameters: [
          p("type_of_damping", "Required", "String", "rayleigh | mass | none | bkt | bkt2 | bkt3 | bkt3f", "Select the material damping formulation."),
          p("simulation_velocity_profile_freq_hz", "Optional", "Float (Hz)", "Default: 0; 0 to simulation_wave_max_freq_hz", "Reference frequency used by BKT-family velocity-profile calculations."),
          p("use_parametric_q_factor", "Conditional", "Yes/no", "Required for type_of_damping=bkt*", "Choose whether BKT-family Q values are generated from the three parametric-Q coefficients."),
          p("use_infinite_qk", "Conditional", "Yes/no", "Required for type_of_damping=bkt*", "Choose the BKT-family infinite-Qk treatment."),
          p("parametric_q_factor_constant", "Conditional", "Float", "Required when use_parametric_q_factor=yes", "Constant term in the parametric Q relation."),
          p("parametric_q_factor_alpha", "Conditional", "Float", "Required when use_parametric_q_factor=yes", "Alpha coefficient in the parametric Q relation."),
          p("parametric_q_factor_beta", "Conditional", "Float", "Required when use_parametric_q_factor=yes", "Beta coefficient in the parametric Q relation."),
          p("the_threshold_damping", "Required", "Float", ">= 0", "Upper limit on the element damping ratio."),
          p("do_damping_statistics", "Optional", "Boolean integer", "Default: 0; 0 | 1", "Print damping statistics during initialization."),
        ],
      },
      {
        id: "absorbing-boundary",
        title: "Absorbing boundary and CPML",
        introduction: "All cpml_* settings are read only when abc_type=cpml. Solid PML elements always use linear elasticity and are excluded from nonlinear constitutive updates, even when their Vs lies within material_properties_list. The present code accepts none or rayleigh for CPML material damping; BKT-family memory is intentionally not applied inside PML elements.",
        parameters: [
          p("abc_type", "Optional", "String", "Default: dashpot; dashpot | cpml", "Select the legacy first-order absorbing dashpot or the unsplit convolutional PML."),
          p("cpml_thickness_m", "Optional", "Float (m)", "Default: 0", "Physical PML thickness. If nonpositive, Hercules derives thickness from cpml_layers and the boundary element size."),
          p("cpml_layers", "Optional", "Integer", "Default: 10", "Number of boundary element layers used when cpml_thickness_m is not positive. At least one thickness specification must be positive."),
          p("cpml_reflection_target", "Optional", "Float", "Default: 1e-3; 0 < R < 1", "Target theoretical reflection coefficient used to derive maximum attenuation."),
          p("cpml_poly_order", "Optional", "Integer", "Default: 2; 1 to 4", "Polynomial order of the attenuation profile."),
          p("cpml_kappa_max", "Optional", "Float", "Default: 1; >= 1", "Maximum real-coordinate stretch. The validated conservative setting is 1; larger values require case validation."),
          p("cpml_alpha_max_hz", "Optional", "Float (Hz)", "Default: -1 (automatic)", "CFS alpha frequency. A positive value sets alpha_max=2πf; a negative value uses π times the simulation maximum frequency."),
          p("cpml_cfl_safety", "Optional", "Float", "Default: 0.9; 0 < value <= 1", "Safety multiplier applied to the CPML critical time-step estimate."),
          p("cpml_mpml_cross_fraction", "Optional", "Float", "Default: 0.1; 0 to 0.5", "Cross-axis M-PML stabilization fraction. Zero disables cross-axis mixing."),
          p("cpml_alpha_floor_fraction", "Optional", "Float", "Default: 0.1; 0 to 1", "Lower bound on each alpha profile as a fraction of alpha_max."),
          p("cpml_partition_weight", "Optional", "Float", "Default: 1; >= 1", "Relative partitioning weight assigned to CPML work; it changes load assignment, not the physical operator."),
          p("cpml_clamp_outer", "Optional", "Yes/no", "Default: no", "Clamp displacement on PML outer faces except the free surface, converting round-trip attenuation to a single-pass outer condition."),
          p("cpml_3d_material_extension", "Optional", "Yes/no", "Default: yes", "Extend the nearest interior 3-D material/water classification through the CPML rather than resampling outside the model footprint."),
          p("cpml_stiffness_mode", "Optional", "String", "Default: native; native | reduced", "Elastic PML stiffness operator. native retains the full effective stiffness plus a CPML correction; reduced is a comparison path."),
          p("cpml_stiffness_mode_acoustic", "Optional", "String", "Default: inherit cpml_stiffness_mode; native | reduced", "Independent stiffness-mode override for acoustic PML elements."),
          p("cpml_pml_material_damping", "Optional", "String", "Default: none; none | rayleigh", "Material damping inside elastic PML. With Rayleigh damping, an omitted value is automatically promoted to rayleigh; explicit none is retained with an instability warning."),
        ],
      },
    ],
  },
  {
    slug: "parameters-coupled-physics",
    title: "Acoustic and poroelastic parameters",
    summary: "Fluid elements, acoustic–solid projection, saturated u-p media, groundwater, and acoustic–poroelastic interface behavior.",
    sections: [
      {
        id: "acoustic",
        title: "Acoustic material and acoustic–elastic coupling",
        parameters: [
          p("include_acoustic", "Optional", "Yes/no", "Default: no", "Enable acoustic classification. With it enabled, material with Vp>0, rho>0, and Vs<=1e-6 is treated as acoustic."),
          p("acoustic_formulation", "Optional", "String", "Default: velocity_potential; velocity_potential | pressure | none", "Select scalar acoustic physics. velocity_potential is the recommended current path; potential is accepted as an alias."),
          p("acoustic_pressure_formulation", "Legacy", "Yes/no", "yes → velocity_potential; no → none", "Backward-compatibility alias used only when acoustic_formulation is absent. New inputs should not use it."),
          p("acoustic_interface_coupling", "Optional", "Yes/no", "Default: yes", "Apply acoustic–elastic normal traction/acceleration coupling on detected interface fragments."),
          p("acoustic_interface_projection", "Optional", "String", "Default: point; point | face_mass", "Choose one-point interface coupling or experimental weak face-mass projection. legacy/one_point and mass/weak/mass_projected are accepted aliases."),
          p("acoustic_interface_face_filter", "Optional", "String", "Default: all; all | bottom_only | bathymetry", "Select all acoustic–elastic faces, only depth-positive bottom faces, or bottom faces with reconstructed bathymetry normals."),
          p("acoustic_interface_max_level", "Optional", "Integer octree level", "Default: -1; -1 to PIXELLEVEL", "Maximum level for extra refinement of mixed acoustic/elastic octants; -1 disables this refinement rule."),
          p("acoustic_mesh_refinement_factor", "Optional", "Float", "Default: 1; >= 1", "Additional refinement multiplier for fully acoustic elements during fresh mesh generation."),
          p("acoustic_interface_diagnostics_file", "Optional", "File path", "Default: unset", "Write one diagnostic row per local acoustic–solid interface fragment. The file may be very large."),
          p("acoustic_quality_factor", "Optional", "Float", "Default: 0 (lossless)", "Positive Q adds a numerical, high-frequency-selective regularizer to interior acoustic elements; it is not a physical ocean-attenuation model."),
        ],
      },
      {
        id: "poroelastic",
        title: "Poroelastic u-p formulation",
        parameters: [
          p("include_poroelastic", "Optional", "Yes/no", "Default: no", "Enable the experimental linear BrickUP-style displacement–pore-pressure formulation for saturated shallow solids."),
          p("poroelastic_formulation", "Optional", "String", "Default: brick_up; brick_up only", "Poroelastic element formulation. brickup is accepted as an alias."),
          p("poroelastic_depth_max_m", "Optional", "Float (m)", "Default: 30; >= 0", "Maximum depth below local ground at which eligible solid, non-PML elements are classified as poroelastic."),
          p("surface_groundwatertable", "Optional", "Yes/no", "Default: no", "When no, saturation starts at the local ground surface. When yes, read a spatial groundwater-depth grid."),
          p("groundwaterlevel_file", "Optional", "File path", "Default: groundwaterlevel.in in the 3-D material directory", "Groundwater depth grid used when surface_groundwatertable=yes. Bare filenames resolve relative to 3D_velocity_model_directory."),
          p("poroelastic_surface_pressure_bc", "Optional", "String", "Default: undrained; undrained | drained", "Pore-pressure condition at the saturation surface. drained applies a weak zero-pressure projection."),
          p("poroelastic_surface_pressure_penalty_factor", "Optional", "Float", "Default: 1; >= 0", "Scale for the weak drained-surface pressure projection; zero disables that projection."),
          p("poroelastic_elastic_interface", "Optional", "String", "Default: no_flow; no_flow only", "Poroelastic–elastic hydraulic condition. Mechanical displacement remains shared while Darcy flux stops at the interface."),
          p("poroelastic_undrained_fast_path", "Optional", "Yes/no", "Default: yes", "Fold alpha²M into volumetric stiffness for undrained, zero-permeability cases while still updating diagnostic pore pressure."),
          p("poroelastic_porosity", "Optional", "Float", "Default: 0.4; 0 < value <= 1", "Connected pore-volume fraction."),
          p("poroelastic_biot_alpha", "Optional", "Float", "Default: 1; 0 to 1", "Biot effective-stress coefficient."),
          p("poroelastic_fluid_bulk_modulus_pa", "Optional", "Float (Pa)", "Default: 2.2e9; > 0", "Pore-fluid bulk modulus."),
          p("poroelastic_fluid_density_kgm3", "Optional", "Float (kg/m3)", "Default: 1000; > 0", "Pore-fluid density."),
          p("poroelastic_fluid_viscosity_pa_s", "Optional", "Float (Pa·s)", "Default: 1e-3; > 0", "Dynamic viscosity used in Darcy diffusion."),
          p("poroelastic_permeability_x_m2", "Optional", "Float (m2)", "Default: 0; >= 0", "Intrinsic permeability along solver X (north)."),
          p("poroelastic_permeability_y_m2", "Optional", "Float (m2)", "Default: 0; >= 0", "Intrinsic permeability along solver Y (east)."),
          p("poroelastic_permeability_z_m2", "Optional", "Float (m2)", "Default: 0; >= 0", "Intrinsic permeability along solver Z (positive depth)."),
          p("acoustic_poroelastic_interface", "Optional", "String", "Default: no_flow; no_flow | open_pore", "Hydraulic condition at acoustic–poroelastic contact. open_pore requires both formulations and weakly enforces pressure continuity."),
          p("acoustic_poroelastic_pressure_penalty_factor", "Optional", "Float", "Default: 1; > 0 for open_pore", "Multiplier for open-pore acoustic–poroelastic pressure-continuity coupling."),
        ],
      },
    ],
  },
  {
    slug: "parameters-output",
    title: "Output, checkpoint, and diagnostic parameters",
    summary: "Monitor files, distributed field output, stations, planes, checkpoints, selected mesh coordinates, and schedule diagnostics.",
    sections: [
      {
        id: "monitor-field-output",
        title: "Monitor and distributed field output",
        parameters: [
          p("monitor_file", "Optional", "File path", "Default: monitor.txt", "Run-progress, configuration, and timing log opened by rank 0."),
          p("simulation_output_rate", "Optional", "Integer time steps", "Default: 1000000; > 0", "Sampling interval for distributed mesh/displacement/velocity output."),
          p("output_parallel", "Optional", "Boolean integer", "Default: 0; 0 | 1", "Enable the distributed field-output subsystem and its summary statistics."),
          p("output_displacement", "Optional", "Boolean integer", "Default: 0; effective only when output_parallel=1", "Write distributed displacement fields."),
          p("output_displacement_file", "Conditional", "File path", "Required when output_parallel=1 and output_displacement=1", "Destination for distributed displacement output."),
          p("output_velocity", "Optional", "Boolean integer", "Default: 0; effective only when output_parallel=1", "Write distributed velocity fields."),
          p("output_velocity_file", "Conditional", "File path", "Required when output_parallel=1 and output_velocity=1", "Destination for distributed velocity output."),
          p("output_stats_file", "Optional", "File path", "Default: output-stats.txt", "Summary statistics destination when output_parallel=1."),
          p("output_debug", "Optional", "Boolean integer", "Default: 0; effective only when output_parallel=1", "Enable the distributed debug-output stream."),
          p("output_debug_file", "Optional", "File path", "Default: output-debug.txt", "Debug-output destination when output_parallel and output_debug are enabled."),
        ],
      },
      {
        id: "stations",
        title: "Station output",
        parameters: [
          p("number_output_stations", "Optional", "Integer", "Default: 0; >= 0", "Number of station rows Hercules reads from output_stations."),
          p("output_stations_print_rate", "Conditional", "Integer time steps", "Required when number_output_stations>0", "Station sampling interval."),
          p("output_stations_directory", "Optional", "Directory path", "Default: outputfiles/stations", "Directory containing station.N histories."),
          p("output_stations", "Conditional", "Table: N rows × 3 floats", "latitude(deg) longitude(deg) depth(m)", "Station coordinates, required when number_output_stations>0. Note the latitude–longitude order."),
          p("station_mediums_output_file", "Optional", "File path", "Default: output_stations_directory/station_mediums.txt", "One-time station ownership, medium, element, local-coordinate, and material metadata."),
          p("print_station_velocities", "Optional", "Yes/no", "Default: no", "Append three-component station velocity."),
          p("print_station_accelerations", "Optional", "Yes/no", "Default: no", "Append three-component acceleration; the current output path also provides velocity when acceleration is requested."),
          p("print_station_pressures", "Optional", "Yes/no", "Default: no", "Append acoustic pressure or poroelastic pore pressure. Requires an active scalar acoustic formulation or poroelastic physics."),
          p("print_station_strains", "Optional", "Yes/no", "Default: no", "Append the six symmetric strain components in XX, YY, ZZ, XY, YZ, XZ order."),
          p("print_station_stresses", "Optional", "Yes/no", "Default: no", "Append the six symmetric stress components in XX, YY, ZZ, XY, YZ, XZ order."),
        ],
      },
      {
        id: "planes",
        title: "Plane output",
        parameters: [
          p("number_output_planes", "Optional", "Integer", "Default: 0; >= 0", "Number of plane definitions to read."),
          p("planes_input_file", "Optional", "File path", "Default: current parameters.in", "File containing output_planes and output_planes_print_rate; it may be separated for large plane lists."),
          p("output_planes_print_rate", "Conditional", "Integer time steps", "Required when number_output_planes>0", "Plane sampling interval."),
          p("output_planes_directory", "Optional", "Directory path", "Default: outputfiles/planes", "Destination for plane metadata and planedisplacements.N files."),
          p("output_planes", "Conditional", "Table: N rows × 9 values", "lat lon depth d_strike n_strike d_dip n_dip strike(deg) dip(deg)", "Plane origins and regular sampling geometry, required when number_output_planes>0."),
          p("reconstruct_acoustic_plane_displacements", "Optional", "Yes/no", "Default: no", "Reconstruct scalar-acoustic plane motion; fully described with the acoustic parameters."),
        ],
      },
      {
        id: "checkpoints-mesh-selection",
        title: "Checkpoints and selected mesh coordinates",
        parameters: [
          p("use_checkpoint", "Optional", "Boolean integer", "Default: 0; 0 | 1", "Restart from checkpoint.in in checkpoint_path."),
          p("checkpointing_rate", "Optional", "Integer time steps", "Default: 0; >= 0", "Checkpoint interval. Zero disables checkpoint generation."),
          p("checkpoint_path", "Optional", "Directory path", "Default: outputfiles/checkpoints", "Directory used for checkpoint.in and rotating checkpoint.out0/checkpoint.out1 files."),
          p("mesh_coordinates_for_matlab", "Optional", "Yes/no", "Default: no", "Write mesh-coordinate data for a selected Cartesian subdomain."),
          p("mesh_coordinates_directory_for_matlab", "Optional", "Directory path", "Default: outputfiles/For_Matlab", "Destination for selected mesh-coordinate output."),
          p("mesh_corners_matlab", "Conditional", "List of 6 floats (m)", "min_x min_y max_x max_y min_z max_z", "Cartesian bounds required when mesh_coordinates_for_matlab=yes."),
        ],
      },
      {
        id: "diagnostics",
        title: "Schedule and initialization diagnostics",
        parameters: [
          p("print_matrix_k", "Optional", "Yes/no", "Default: no", "Print the element stiffness matrix for diagnostic inspection."),
          p("enable_timing_barriers", "Optional", "Boolean integer", "Default: 0", "Insert selected MPI barriers around timers. Useful for diagnosis but changes measured execution behavior."),
          p("forces_buffer_size", "Optional", "Integer bytes", "Default: 104857600", "Maximum in-memory buffer used while generating and staging source-force files."),
          p("schedule_print_error_check", "Optional", "Boolean integer", "Default: 0", "Enable additional communication-schedule error checks."),
          p("schedule_print_stdout", "Optional", "Boolean integer", "Default: 0", "Print communication-schedule details to standard output."),
          p("schedule_print_file", "Optional", "Boolean integer", "Default: 0", "Write communication-schedule details to schedule_info.txt."),
          p("stat_schedule_filename", "Optional", "File path", "Default: stat-sched.txt", "Communication-schedule statistics file."),
          p("stat_mesh_filename", "Optional", "File path", "Default: stat-mesh.txt", "Mesh statistics file."),
          p("stat_topo_filename", "Optional", "File path", "Default: stat-topo.txt", "Topography statistics file."),
        ],
      },
    ],
  },
  {
    slug: "parameters-topography-buildings",
    title: "Topography and building parameters",
    summary: "Surface geometry, VT/FEM topographic treatment, explicit buildings, foundation records, and push-down regions.",
    sections: [
      {
        id: "topography",
        title: "Topography",
        parameters: [
          p("include_topography", "Optional", "Yes/no", "Default: no", "Enable surface topography."),
          p("topography_directory", "Optional", "Directory path", "Default: inputfiles/topography", "Directory containing topography.in."),
          p("maximum_octant_level", "Conditional", "Integer octree level", "Required with topography; 0 to 30", "Maximum refinement level allowed by the topography mesher."),
          p("computation_method", "Conditional", "String", "Required with topography; vt | fem", "Topographic element treatment: volume tetrahedra (vt) or FEM."),
          p("topographybase_zcoord", "Conditional", "Float (m)", "Required with topography; > 0", "Reference depth coordinate below which the topographic surface is embedded."),
          p("type_of_etree", "Conditional", "String", "Required with topography; sqd | full", "Etree representation used by the topographic mesh path."),
          p("min_enclosed_volume_perc", "Optional", "Float", "Default: 0.1", "Minimum enclosed-volume fraction used by VT topographic element decisions."),
          p("consider_nonlinear_topography", "Optional", "Yes/no", "Default: no", "Apply nonlinear treatment to topographic elements; requires include_nonlinear_analysis=yes."),
          p("topography_bkt_gpu_stabilize", "Optional", "String", "Default: auto; auto | yes | no", "VT+BKT GPU stabilization. auto enables it only for the affected VT plus BKT-family configuration."),
        ],
      },
      {
        id: "buildings",
        title: "Buildings and push-down regions",
        parameters: [
          p("include_buildings", "Optional", "Yes/no", "Default: no", "Enable explicit building and foundation regions."),
          p("number_of_buildings", "Conditional", "Integer", "Required when include_buildings=yes; >= 0", "Number of rows in building_properties."),
          p("number_of_pushdowns", "Conditional", "Integer", "Required when include_buildings=yes; >= 0", "Number of rows in pushdown_properties."),
          p("min_octant_size_m", "Conditional", "Float (m)", "Required when include_buildings=yes; >= 0", "Minimum octant size used while resolving buildings."),
          p("surface_shift_m", "Conditional", "Float (m)", "Required when include_buildings=yes; >= 0", "Vertical surface shift applied by the building model and to associated source/station depths."),
          p("consider_fixed_base", "Conditional", "Yes/no", "Required when include_buildings=yes", "Enable prescribed fixed-base input histories."),
          p("fixedbase_input_dt", "Conditional", "Float (s)", "Required when consider_fixed_base=yes", "Sampling interval of the fixed-base input records."),
          p("fixedbase_input_dir", "Conditional", "Directory path", "Required when consider_fixed_base=yes", "Directory containing fixed-base input files."),
          p("fixedbase_input_startindex", "Conditional", "Integer", "Required when consider_fixed_base=yes", "Index of the first fixed-base record; successive buildings increment the index."),
          p("fixedbase_input_sufix", "Conditional", "String", "Required when consider_fixed_base=yes", "File prefix/suffix token used with fixedbase_input_startindex. The parameter name intentionally retains the historical sufix spelling."),
          p("building_properties", "Conditional", "Table: N rows × 12 floats", "xmin xmax ymin ymax depth height Vp_b Vs_b rho_b Vp_f Vs_f rho_f", "Building bounds, dimensions, and building/foundation material properties; required when number_of_buildings>0."),
          p("pushdown_properties", "Conditional", "Table: N rows × 6 floats", "xmin xmax ymin ymax depth height", "Push-down bounds and dimensions; required when number_of_pushdowns>0. Current code forces the supplied depth field to zero."),
        ],
      },
    ],
  },
  {
    slug: "parameters-nonlinear",
    title: "Nonlinear and plastic analysis parameters",
    summary: "Geostatic initialization, elastoplastic and hysteretic constitutive models, plasticity integration controls, and the material-property table.",
    sections: [
      {
        id: "nonlinear-controls",
        title: "Nonlinear and plastic controls",
        introduction: "With include_nonlinear_analysis=yes, eligible solid elements perform history-dependent constitutive updates, including yield evaluation, plastic-strain evolution, hardening or nonlinear-backbone state, and nonlinear station diagnostics. Every parameter on this page except the feature switch is required when the feature is enabled.",
        parameters: [
          p("include_nonlinear_analysis", "Optional", "Yes/no", "Default: no", "Enable nonlinear solid analysis, constitutive updates, plastic-strain evolution, and nonlinear station diagnostics in eligible elements. Acoustic elements are always excluded. Istanbul solid elements use the legacy nl_id classification; other solid elements ignore nl_id and use the Vs range covered by material_properties_list."),
          p("geostatic_loading_time_sec", "Conditional", "Float (s)", ">= 0; required for nonlinear analysis", "Duration used to establish geostatic loading."),
          p("geostatic_cushion_time_sec", "Conditional", "Float (s)", ">= 0; required for nonlinear analysis", "Post-loading cushion duration. Loading plus cushion may not exceed the simulation end time."),
          p("approximate_geostatic_state", "Conditional", "Yes/no", "Required for nonlinear analysis", "Use an approximate at-rest stress state. It must be no when either geostatic time is positive."),
          p("error_tolerance", "Conditional", "Float", "> 0; required for nonlinear analysis", "Constitutive convergence/error tolerance. For ManzariDafalias this is the maximum absolute yield-surface residual in Pa after plastic stress correction; the elastic predictor separately uses the one-sided condition f_trial <= error_tolerance."),
          p("backbone_errTol", "Conditional", "Float", "> 0; required for nonlinear analysis", "Tolerance used while matching nonlinear backbone behavior. ManzariDafalias reuses it as the dimensionless local error target for adaptive Modified Euler."),
          p("material_model", "Conditional", "String", "linear | vonMises_ep | vonMises_fa | vonMises_faM | vonMises_baE | vonMises_baH | vonMises_GQH | vonMises_MKZ | vonMises_RO | MohrCoulomb | DruckerPrager | ManzariDafalias", "Constitutive family: elastic diagnostic, von Mises/cyclic plasticity, nonlinear backbone response, frictional plasticity, or the pressure-dependent critical-state Manzari–Dafalias sand model. ManzariDafalias production updates currently require the GPU backend."),
          p("material_plasticity_type", "Conditional", "String", "rate_dependent | rate_independent", "Select the integration family for plastic constitutive response."),
          p("material_properties_count", "Conditional", "Integer", ">= 1", "Number of rows in material_properties_list."),
          p("material_properties_list", "Conditional", "Model-specific table", "Most models: N × 16; ManzariDafalias: N × 17", "Rows are ordered by increasing, distinct Vs_limit; the first and last limits define the inclusive nonlinear-element range and active values are linearly interpolated by Vs. Existing models use: Vs_limit alpha_cohesion K_phi strain_rate sensitivity hardening beta_dilatancy gamma0 psi m theta1 theta2 theta3 theta4 theta5 phi_RO. With vonMises_baE, alpha_cohesion is Su in Pa. ManzariDafalias uses: Vs_limit G0 e_init Mc c lambda_c e0 ksi P_atm m h0 ch nb A0 nd z_max cz. Hercules supplies no implicit calibrated ManzariDafalias constants; every row must provide all 17 values. Its Poisson ratio and density come from each element's Vp, Vs, and rho. Acoustic elements remain excluded even when the range starts at zero; CPML elements also remain linear."),
          p("no_substeps", "Conditional", "Integer", ">= 1; required for nonlinear analysis", "Constitutive substep control. For ManzariDafalias it caps accepted adaptive Modified-Euler substeps per global increment; rejected error-control trials are also finitely bounded."),
          p("stiff_damp", "Conditional", "Float", "Required for nonlinear analysis", "Stiffness-proportional damping coefficient used by the nonlinear path."),
          p("freq_stiff_damp", "Conditional", "Float (Hz)", "Required for nonlinear analysis", "Reference frequency associated with nonlinear stiffness damping."),
          p("tension_cutoff", "Conditional", "Yes/no", "Required for nonlinear analysis", "Remove tensile stress. yes is accepted only for MohrCoulomb and DruckerPrager."),
        ],
      },
      {
        id: "nonlinear-material-point-test",
        title: "Standalone material-point test controls",
        introduction: "These keys belong only in an input passed to psolve --nonlinear-material-point. They are not part of a production parameters.in simulation. The test also reads material_model, no_substeps, error_tolerance, backbone_errTol, and tension_cutoff from the same file.",
        parameters: [
          p("material_point_delta_time_sec", "Conditional", "Float (s)", "> 0; material-point subcommand only", "Time increment associated with each prescribed material-point strain row."),
          p("material_point_steps", "Conditional", "Integer", "> 0; material-point subcommand only", "Number of rows in material_point_strain_history."),
          p("material_point_constants", "Conditional", "Model-specific list", "24 floats normally; 17 for ManzariDafalias", "Existing fixtures use the documented direct nlconstants_t order. ManzariDafalias uses G0 nu e_init Mc c lambda_c e0 ksi P_atm m h0 ch nb A0 nd z_max cz. Stress-like entries use Pa and angles use radians."),
          p("material_point_initial_strain", "Conditional", "6 floats", "epsilon_xx epsilon_yy epsilon_zz epsilon_xy epsilon_yz epsilon_xz; ManzariDafalias only", "Optional committed initial tensor strain for a seeded ManzariDafalias replay. If any material-point initial-state array is present, all five are required. Off-diagonal entries are tensorial shear strains."),
          p("material_point_initial_stress", "Conditional", "6 floats (Pa)", "sigma_xx sigma_yy sigma_zz sigma_xy sigma_yz sigma_xz; ManzariDafalias only", "Optional committed initial effective-stress tensor for a seeded ManzariDafalias replay. Hercules uses tension-positive stress."),
          p("material_point_initial_alpha", "Conditional", "6 floats", "alpha_xx alpha_yy alpha_zz alpha_xy alpha_yz alpha_xz; ManzariDafalias only", "Optional committed ManzariDafalias backstress-ratio tensor. It must be supplied with the other four initial-state arrays."),
          p("material_point_initial_alpha_in", "Conditional", "6 floats", "alpha_in_xx alpha_in_yy alpha_in_zz alpha_in_xy alpha_in_yz alpha_in_xz; ManzariDafalias only", "Optional committed ManzariDafalias reversal-memory tensor. It must be supplied with the other four initial-state arrays."),
          p("material_point_initial_fabric", "Conditional", "6 floats", "fabric_xx fabric_yy fabric_zz fabric_xy fabric_yz fabric_xz; ManzariDafalias only", "Optional committed ManzariDafalias fabric tensor. It must be supplied with the other four initial-state arrays."),
          p("material_point_strain_history", "Conditional", "Table: N rows × 6 floats", "epsilon_xx epsilon_yy epsilon_zz epsilon_xy epsilon_yz epsilon_xz", "Prescribed tensor-strain history for the material-point subcommand. Off-diagonal entries are tensorial, not engineering, shear strains."),
          p("material_point_compare_rtol", "Optional", "Float", "Default: 1e-6; >= 0; material-point subcommand only", "Relative CPU/GPU state-comparison tolerance."),
          p("material_point_compare_atol", "Optional", "Float", "Default: 1e-9; >= 0; material-point subcommand only", "Absolute CPU/GPU state-comparison tolerance."),
          p("material_point_output_prefix", "Optional", "Path prefix", "Default: material-point; material-point subcommand only", "Prefix for the generated .cpu.csv and .gpu.csv state histories."),
          p("material_point_brickup", "Optional", "Yes/no", "Default: no; material-point subcommand only", "Enable the sealed uniform one-element BrickUP pressure regression alongside a ManzariDafalias strain replay. The GPU side executes the production pressure-force and pressure-update kernels; this is not a production parameters.in control."),
          p("material_point_brickup_porosity", "Conditional", "Float", "0 < n <= 1; required when material_point_brickup=yes", "Fluid-filled porosity used to form the sealed storage modulus M=Kf/n for the standalone SSPbrickUP/BrickUP comparison."),
          p("material_point_brickup_fluid_bulk_modulus_pa", "Conditional", "Float (Pa)", "> 0; required when material_point_brickup=yes", "Pore-fluid bulk modulus Kf for the standalone sealed BrickUP pressure update."),
          p("material_point_brickup_biot_alpha", "Conditional", "Float", "Currently exactly 1; required when material_point_brickup=yes", "Biot coefficient for the standalone SSPbrickUP/BrickUP comparison. The current benchmark restricts alpha to one to match the OpenSees element formulation."),
          p("material_point_initial_pore_pressure_pa", "Conditional", "Float (Pa)", "Finite; required when material_point_brickup=yes", "Committed pore pressure at the first seeded replay state. Generated drained-to-undrained comparisons use zero."),
        ],
      },
    ],
  },
  {
    slug: "parameters-drm-waves",
    title: "DRM, incident-wave, and half-space parameters",
    summary: "Domain Reduction Method stages, analytic or recorded incident plane waves, and the homogeneous half-space override.",
    sections: [
      {
        id: "drm",
        title: "Domain Reduction Method",
        parameters: [
          p("implement_drm", "Optional", "Yes/no", "Default: no", "Enable the three-part Domain Reduction Method workflow."),
          p("drm_directory", "Conditional", "Directory path", "Required by the DRM module", "Root directory containing part0, part1, and part2 DRM data. Although the main parser initializes a conventional path, the DRM module currently requires the key to be present."),
          p("which_drm_part", "Conditional", "String", "part0 | part1 | part2", "Select coordinate generation, response recording, or effective-force replay."),
          p("drm_edgesize", "Conditional", "Float (m)", "Required when implement_drm=yes", "Element edge size prescribed for the DRM box."),
          p("drm_offset_x", "Conditional", "Float (m)", "Required when implement_drm=yes", "DRM-box offset from the domain origin along local X/north."),
          p("drm_offset_y", "Conditional", "Float (m)", "Required when implement_drm=yes", "DRM-box offset from the domain origin along local Y/east."),
          p("drm_print_rate", "Conditional", "Integer time steps", "Required when implement_drm=yes", "DRM output sampling interval."),
          p("part1_delta_t", "Conditional", "Float (s)", "Required when implement_drm=yes", "Time step used by the part1 data that part2 will read; normally matches the part1 simulation delta time."),
          p("part1_simulation_time", "Conditional", "Float (s)", "Required when implement_drm=yes", "Duration of the part1 data set that part2 will read."),
          p("drm_boundary", "Conditional", "List of 5 floats (m)", "min_x min_y max_x max_y depth", "DRM boundary dimensions relative to the DRM offset box."),
        ],
      },
      {
        id: "incident-plane-waves",
        title: "Incident plane waves",
        parameters: [
          p("include_incident_planewaves", "Optional", "Yes/no", "Default: no", "Enable the DRM-box incident plane-wave forcing path."),
          p("type_of_wave", "Conditional", "String", "SV | P", "Incident wave type."),
          p("DRMBox_Noelem_Halfwidth_EW", "Conditional", "Float (elements)", "Required for incident waves", "DRM-box half-width in the east–west direction, expressed as an element count."),
          p("DRMBox_Noelem_Halfwidth_NS", "Conditional", "Float (elements)", "Required for incident waves", "DRM-box half-width in the north–south direction, expressed as an element count."),
          p("DRMBox_Noelem_depth", "Conditional", "Float (elements)", "Required for incident waves", "DRM-box depth expressed as an element count."),
          p("DRMBox_element_size", "Conditional", "Float (m)", "Required for incident waves", "Uniform DRM-box element edge size."),
          p("DRM_xo", "Conditional", "Float (m)", "Required for incident waves", "DRM-box local X origin."),
          p("DRM_yo", "Conditional", "Float (m)", "Required for incident waves", "DRM-box local Y origin."),
          p("fnc_type", "Conditional", "String", "RICKER | time_hist", "Incident-wave time function."),
          p("planewave_strike", "Conditional", "Float (degrees)", "Required for incident waves", "Horizontal propagation strike, converted internally to radians."),
          p("planewave_Z_angle", "Conditional", "Float (degrees)", "Required for incident waves", "Incidence angle relative to the solver Z convention, converted internally to radians."),
          p("Ts", "Conditional", "Float (s)", "Required for incident waves", "Time shift used by the analytic Ricker displacement."),
          p("fc", "Conditional", "Float (Hz)", "Required for incident waves", "Ricker central frequency."),
          p("Uo", "Conditional", "Float", "Required for incident waves", "Incident displacement-amplitude scale."),
          p("no_components", "Conditional", "Integer", "Normally 1 or 2", "Number of recorded horizontal components used by time_hist; 2 also requires ug_alongnormal."),
          p("ug_timestep", "Conditional", "Float (s)", "Required for incident waves", "Sampling interval of time-history ground motion; the current parser requires it even for RICKER."),
          p("ug_alongstrike", "Conditional", "File path", "Required when fnc_type=time_hist", "Count-prefixed displacement history along strike."),
          p("ug_alongnormal", "Conditional", "File path", "Required when fnc_type=time_hist and no_components=2", "Count-prefixed displacement history normal to strike."),
        ],
      },
      {
        id: "homogeneous-halfspace",
        title: "Homogeneous half-space",
        parameters: [
          p("include_hmg_halfspace", "Optional", "Yes/no", "Default: no", "Enable the homogeneous half-space material override in its etree-aligned box."),
          p("etreeBox_xo", "Conditional", "Float (m)", "Required when include_hmg_halfspace=yes", "Local X origin of the homogeneous box."),
          p("etreeBox_yo", "Conditional", "Float (m)", "Required when include_hmg_halfspace=yes", "Local Y origin of the homogeneous box."),
          p("etreeBox_lengthEW", "Conditional", "Float (m)", "Required when include_hmg_halfspace=yes", "East–west length of the homogeneous box."),
          p("etreeBox_lengthNS", "Conditional", "Float (m)", "Required when include_hmg_halfspace=yes", "North–south length of the homogeneous box."),
          p("etreeBox_depth", "Conditional", "Float (m)", "Required when include_hmg_halfspace=yes", "Depth of the homogeneous box."),
          p("hmgHS_Vs", "Conditional", "Float (m/s)", "Required when include_hmg_halfspace=yes", "Homogeneous half-space shear-wave speed."),
          p("hmgHS_Vp", "Conditional", "Float (m/s)", "Required when include_hmg_halfspace=yes", "Homogeneous half-space compressional-wave speed."),
          p("hmgHS_rho", "Conditional", "Float (kg/m3)", "Required when include_hmg_halfspace=yes", "Homogeneous half-space density."),
        ],
      },
    ],
  },
];

function toDocSection(section: ParameterSection): DocSection {
  return {
    id: section.id,
    title: section.title,
    paragraphs: section.introduction ? [section.introduction] : undefined,
    table: {
      headers: ["Parameter", "Requirement", "Type / format", "Default / accepted values", "Meaning"],
      rows: section.parameters.map((parameter) => [
        parameter.name,
        parameter.requirement,
        parameter.type,
        parameter.values,
        parameter.description,
      ]),
    },
  };
}

export const parameterReferencePages: Record<string, DocPageContent> = Object.fromEntries(
  parameterPages.map((page) => [
    page.slug,
    {
      eyebrow: "Parameter reference",
      title: page.title,
      summary: page.summary,
      sections: page.sections.map(toDocSection),
    },
  ]),
);

export const documentedParameters = parameterPages.flatMap((page) =>
  page.sections.flatMap((section) => section.parameters),
);

export const documentedParameterNames = documentedParameters.map((parameter) => parameter.name);

export const parameterReferenceSummary = parameterPages.map((page) => ({
  slug: page.slug,
  title: page.title,
  count: page.sections.reduce((total, section) => total + section.parameters.length, 0),
}));
