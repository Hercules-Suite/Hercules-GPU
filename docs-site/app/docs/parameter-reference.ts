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
        introduction: String.raw`Define \(L_{\mathrm{north}}=\mathtt{region\_length\_north\_m}\), \(L_{\mathrm{east}}=\mathtt{region\_length\_east\_m}\), \(L_{\mathrm{depth}}=\mathtt{region\_depth\_deep\_m}-\mathtt{region\_depth\_shallow\_m}\), and \(L_{\max}=\max(L_{\mathrm{north}},L_{\mathrm{east}},L_{\mathrm{depth}})\). Every extent \(D\) must be an integer multiple of \(L_{\max}/2^N\) for a common level \(N\). Use exact whole-metre extents with the current Octor path.`,
        parameters: [
          p("region_origin_latitude_deg", "Required", "Float (degrees)", "Any valid latitude", "Latitude of the lower-left surface corner of the model domain."),
          p("region_origin_longitude_deg", "Required", "Float (degrees)", "Any valid longitude", "Longitude of the lower-left surface corner of the model domain."),
          p("region_depth_shallow_m", "Required", "Float (m)", "Usually >= 0", "Depth of the shallow domain boundary, positive downward."),
          p("region_length_east_m", "Required", "Float (m)", "> 0; octree-compatible", String.raw`Domain length in the local east direction; this becomes solver Y. It must be an integer multiple of the common dyadic increment \(L_{\max}/2^N\).`),
          p("region_length_north_m", "Required", "Float (m)", "> 0; octree-compatible", String.raw`Domain length in the local north direction; this becomes solver X. It must be an integer multiple of the common dyadic increment \(L_{\max}/2^N\).`),
          p("region_depth_deep_m", "Required", "Float (m)", "> region_depth_shallow_m", "Depth of the bottom boundary. The octree-compatible vertical extent is deep minus shallow depth, not the deep coordinate alone."),
          p("region_azimuth_leftface_deg", "Required", "Float (degrees)", "Finite degrees; no default", "With PROJ, rotate solver +X clockwise from UTM grid north about the region origin; +Y follows at 90 degrees and Z stays down. Any finite angle is normalized modulo 360. Length north/east means local X/Y length. Nonzero rotation requires cvmdb_input_file=profile and Istanbul_velocity_model=no; general 3-D and basin models are supported. Without PROJ use consistent domain_surface_corners. Rebuild meshes and sources after changing the frame."),
          p("domain_surface_corners", "Conditional", "Table: 4 rows × 2 floats", "longitude latitude", "Required only when Hercules is built without PROJ. Give the four surface corners clockwise, starting at the origin."),
        ],
      },
      {
        id: "time-resolution",
        title: "Time integration and spatial resolution",
        introduction: String.raw`Let \(V_{S,\min}^{\mathrm{model}}\) be the smallest shear-wave velocity in the active solid material model and define \(V_{S,\min}^{\mathrm{effective}}=\max(V_{S,\min}^{\mathrm{model}},\mathtt{simulation\_shear\_velocity\_min})\). For elastic material, the nominal target is \(h\le V_{S,\min}^{\mathrm{effective}}/(\mathtt{simulation\_node\_per\_wavelength}\,\mathtt{simulation\_wave\_max\_freq\_hz})\). Acoustic elements use Vp and may apply an additional refinement factor.`,
        parameters: [
          p("simulation_start_time_sec", "Optional", "Float (s)", "Default: 0; >= 0", "Simulation start time. It must not exceed simulation_end_time_sec."),
          p("simulation_end_time_sec", "Required", "Float (s)", ">= simulation_start_time_sec", "Simulation end time."),
          p("simulation_delta_time_sec", "Required", "Float (s)", "> 0", "Explicit time-step interval. Select it from the fastest wave speed and smallest element in the final model."),
          p("simulation_wave_max_freq_hz", "Required", "Float (Hz)", "> 0", "Highest frequency the mesh is intended to resolve; used with simulation_node_per_wavelength to set octree refinement."),
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
          p("stiffness_calculation_method", "Optional", "String", "Default: effective; effective | conventional", "Finite-element stiffness operator. effective is the current recommended and default method; conventional is retained for compatibility. CPML requires effective for non-BKT runs."),
          p("use_progressive_meshing", "Optional", "Integer", "Default: 1; >= 0", "Progressive-meshing factor/iteration control. Current production inputs normally retain the default."),
        ],
      },
      {
        id: "velocity-models",
        title: "Optional regional material models",
        parameters: [
          p("Istanbul_velocity_model", "Optional", "Yes/no", "Default: no", "Enable the built-in Istanbul material-model replacement where it covers the domain."),
          p("basin_velocity_model", "Optional", "Yes/no", "Default: no", "Enable a basin geometry that replaces the background material inside its boundary."),
          p("basin_origin_latitude_deg", "Conditional", "Float (degrees)", "Required when basin_velocity_model=yes", "Latitude of the basin grid's (0, 0) origin; grid X is distance north and Y is distance east."),
          p("basin_origin_longitude_deg", "Conditional", "Float (degrees)", "Required when basin_velocity_model=yes", "Longitude of the basin grid's (0, 0) origin; grid X is distance north and Y is distance east."),
          p("basin_input_file", "Conditional", "Path", "Required when basin_velocity_model=yes", "Basin-bottom text grid: Nx Ny 1, then Nx × Ny rows of x y bottom_depth in metres, X fastest. Start at (0, 0), with Nx/Ny >= 2 and equal positive X/Y spacing; depth uses solver Z, positive downward. Relative paths use the solver working directory. Extend the grid beyond maximum query coordinates. See Input files and directories: Basin-bottom grid format for the example and edge limitation."),
          p("basin_material_type", "Optional", "String", "Default: elastic; elastic | acoustic_water", "Material assigned inside basin geometry. acoustic_water requires include_acoustic=yes."),
          p("basin_water_vp", "Optional", "Float (m/s)", "Default: 1500; > 0", "Compressional-wave speed used when basin_material_type=acoustic_water."),
          p("basin_water_rho", "Optional", "Float (kg/m3)", "Default: 1000; > 0", "Density used when basin_material_type=acoustic_water."),
          p("3D_velocity_model", "Optional", "Yes/no", "Default: no", "Enable the structured general 3-D velocity model."),
          p("seabed_geometry_file", "Optional", "File path", "Default: omitted (legacy Vs-based water)", "Enable explicit ground/bed geometry. Water coverage and depth are inferred from topography minus bed with a fixed 1e-8 m tolerance. Requires topography, 3D velocity model and acoustics; disables buildings and Istanbul, basin and homogeneous-halfspace overrides. Path is relative to the solver working directory. Text format: topography spacing followed by nodal ground/seabed elevations only; east varies fastest, and one number per line is valid. Same grid and vertical datum as topography.in. Bed above the outer surface is rejected beyond tolerance; coincident surfaces are dry. Wet grid cells require flat, mesh-face-aligned outer surfaces below topographybase_zcoord; separated water bodies may have different elevations. Solid-only material profiles required. Rebuild the mesh; see Input files: Explicit seabed geometry."),
          p("3D_velocity_model_depth_reference", "Conditional", "String", "Default: legacy; legacy | model_z | surface | seabed", "With seabed_geometry_file, explicitly choose model_z (fixed solver Z), surface (depth below each column's outer air/material boundary), or seabed (below bed; ground on dry cells). These references work with both full and sqd; the reader restores physical Z before converting each column. Non-legacy values require the geometry file. Without it, legacy retains type_of_etree conventions."),
          p("seabed_water_vp_m_s", "Conditional", "Float (m/s)", "Required with seabed_geometry_file; finite > 0", "Uniform acoustic Vp in the geometry-defined water volume. Vs is zero. Choose case-appropriate water properties; no implicit value is supplied."),
          p("seabed_water_density_kg_m3", "Conditional", "Float (kg/m3)", "Required with seabed_geometry_file; finite > 0", "Uniform water density in SI units, also used for nonlinear geostatic water load. This differs from rho.in, which stores tonne/m3. Keep pore-fluid density and hydrostatic assumptions physically consistent for poroelastic sediment."),
          p("3D_velocity_model_directory", "Optional", "Directory path", "Default: inputfiles/materialfiles", "Directory containing x.in, y.in, z.in, vs.in, vp.in, rho.in, index.in, and optional matched qs.in/qp.in files. Solid Q must be finite and positive; acoustic placeholders may be zero. Q interpolation excludes water samples and renormalizes solid weights; supplied Q is not rescaled by the Vs minimum."),
          p("3D_velocity_model_azimuth_deg", "Optional", "Float (degrees)", "Default: region_azimuth_leftface_deg; finite", "PROJ: material-grid +X azimuth clockwise from the domain UTM zone grid north. If omitted, inherits region_azimuth_leftface_deg so material and solver axes align. Explicit values override this default: 0 keeps the material grid aligned with UTM north/east. x.in/y.in are offsets in the material axes about the material origin. Queries and cached topography columns transform between frames; CPML clamps in solver coordinates first. Without PROJ, omission or explicit 0 retains solver-aligned material axes; explicit nonzero overrides require PROJ."),
          p("3D_velocity_model_origin_latitude_deg", "Conditional", "Float (degrees)", "Required when 3D_velocity_model=yes", "Latitude of the structured model's (x=0,y=0) origin."),
          p("3D_velocity_model_origin_longitude_deg", "Conditional", "Float (degrees)", "Required when 3D_velocity_model=yes", "Longitude of the structured model's (x=0,y=0) origin."),
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
          p("use_parametric_q_factor", "Conditional", "Yes/no", "Required for type_of_damping=bkt*", "When yes, compute Qs=C+alpha×Vs_km_s^beta from the three parametric-Q coefficients. When no, use material Q values when the selected profile or 3-D model supplies them."),
          p("use_infinite_qk", "Conditional", "Yes/no", "Required for type_of_damping=bkt*", "When yes, suppress bulk attenuation for all solids. When no, derive Qk from Qp, Qs, Vp, and Vs; negative bulk attenuation beyond roundoff selects Qk=1000 (lossless bulk) only for affected elements, with a global fallback count. Valid bulk and shear attenuation remain unchanged. Zero inverse bulk Q remains lossless. The supplied Qp is not reproduced in fallback elements. Solid Q inputs must still be finite and positive. Regenerate persisted meshes to apply this policy."),
          p("parametric_q_factor_constant", "Conditional", "Float", "Required when use_parametric_q_factor=yes", "Constant C in Qs=C+alpha×Vs_km_s^beta."),
          p("parametric_q_factor_alpha", "Conditional", "Float", "Required when use_parametric_q_factor=yes", "Multiplier alpha in Qs=C+alpha×Vs_km_s^beta."),
          p("parametric_q_factor_beta", "Conditional", "Float", "Required when use_parametric_q_factor=yes", "Exponent beta in Qs=C+alpha×Vs_km_s^beta."),
          p("the_threshold_damping", "Required", "Float", ">= 0", "Upper limit on the element damping ratio."),
          p("do_damping_statistics", "Optional", "Boolean integer", "Default: 0; 0 | 1", "Print damping statistics during initialization."),
        ],
      },
      {
        id: "absorbing-boundary",
        title: "Absorbing boundary and CPML",
        introduction: "All cpml_* settings are read only when abc_type=cpml. Solid PML elements always use linear elasticity and are excluded from nonlinear constitutive updates, even when their Vs lies within material_properties_list. The present code accepts none or rayleigh for CPML material damping; BKT-family memory is intentionally not applied inside PML elements. CPML requires effective GPU stiffness for non-BKT runs and disables checkpoint reads/writes.",
        parameters: [
          p("abc_type", "Optional", "String", "Default: dashpot; dashpot | cpml", "Select the legacy first-order absorbing dashpot or the unsplit convolutional PML."),
          p("cpml_thickness_m", "Optional", "Float (m)", "Default: 0", "Physical PML thickness. If nonpositive, Hercules uses cpml_layers times the largest non-air element edge touching a lateral or bottom boundary (MPI-wide maximum). Air is also excluded from the global fallback when no non-air boundary element exists."),
          p("cpml_layers", "Optional", "Integer", "Default: 10", "Multiplier of the largest non-air boundary element edge when cpml_thickness_m is not positive. Sets one physical thickness; finer regions contain more than this many cells. At least one thickness specification must be positive."),
          p("cpml_reflection_target", "Optional", "Float", "Default: 0.1; 0 < R < 1", "Target theoretical reflection coefficient used to derive maximum attenuation."),
          p("cpml_poly_order", "Optional", "Integer", "Default: 2; 1 to 4", "Polynomial order of the attenuation profile."),
          p("cpml_kappa_max", "Optional", "Float", "Default: 10; >= 1", "Maximum real-coordinate stretch. All stiffness, inertia and surface products include this stretch; larger values require case validation."),
          p("cpml_alpha_max_hz", "Optional", "Float (Hz)", "Default: simulation_wave_max_freq_hz / 2 (automatic; stored as -1)", "CFS alpha frequency in Hz. Omitted or nonpositive values select half the simulation maximum frequency, giving angular alpha_max=π×simulation_wave_max_freq_hz. A positive numeric value overrides this with alpha_max=2πf; the input parser does not evaluate arithmetic expressions."),
          p("cpml_cfl_safety", "Optional", "Float", "Default: 0.9; 0 < value <= 1", "Safety multiplier applied to the CPML critical time-step estimate."),
          p("cpml_mpml_cross_fraction", "Optional", "Float", "Default: 0; 0 to 0.5", "Experimental cross-axis damping/profile mixing. The complete operator defaults to plain CPML; nonzero values break separable matching and can increase reflection. Explicit older values remain active."),
          p("cpml_acoustic_sponge_fraction", "Optional", "Float (dimensionless)", "Default: 0.3; finite and >= 0; 0 disables", "Experimental additional acoustic volume mass damping a=gamma×(d_x+d_y+d_z) inside CPML, for pressure or velocity-potential acoustics. This control leaves gravity-surface, interface, elastic and interior coefficients unchanged. Independent of acoustic_quality_factor. Nonzero values change matching; validate raw reflection errors and long-time stability."),
          p("cpml_alpha_floor_fraction", "Optional", "Float", "Default: 0.1; 0 to 1", "Lower bound on each alpha profile as a fraction of alpha_max."),
          p("cpml_partition_weight", "Optional", "Float", "Default: 1; >= 1", "Relative partitioning weight assigned to CPML work; it changes load assignment, not the physical operator."),
          p("cpml_clamp_outer", "Optional", "Yes/no", "Default: yes", "With abc_type=cpml, zero elastic displacement and the acoustic scalar (pressure or velocity potential) on the outer lateral faces and bottom. The lateral rim of a gravity-enabled top surface is included. Pressure-release acoustic top nodes use their dedicated update first. Explicit no removes the clamp while retaining CPML absorption; neither setting guarantees a reflection-free outer termination."),
          p("cpml_3d_material_extension", "Optional", "Yes/no", "Default: yes", "With abc_type=cpml and 3D_velocity_model=yes, extend inner-PML material properties and water classification along each absorbing normal. Clamp in absolute model coordinates, accounting for SQD topography and building depth offsets at the bottom and corners. Rebuild or recorrect imported mesh properties after changing this setting or adopting the corrected coordinate conversion."),
          p("cpml_stiffness_mode", "Optional", "String", "Default: native; native | reduced", "Elastic PML stiffness operator. native integrates all Hex8 modes with full Jacobian products; reduced is rank deficient and intended only for diagnostics."),
          p("cpml_stiffness_mode_acoustic", "Optional", "String", "Default: inherit cpml_stiffness_mode; native | reduced", "Independent stiffness-mode override for acoustic PML elements."),
          p("cpml_pml_material_damping", "Optional", "String", "Default: none; none | rayleigh", "Material damping inside elastic PML. With Rayleigh damping, an omitted value is automatically promoted to rayleigh; explicit none is retained with a material-damping interface warning."),
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
          p("acoustic_free_surface_gravity", "Optional", "Yes/no", "Default: no", "Enable the linearized gravity condition at a flat acoustic free surface. With seabed_geometry_file, complete horizontal top faces must coincide with local topography over water; separate flat water bodies may have different elevations. Without explicit seabed geometry, the surface follows solver Z=topographybase_zcoord with topography, surface_shift_m with buildings, or the octree top otherwise; topographic face centers must also be exposed, excluding buried water beneath land. no applies pressure-release p=0; yes adds lumped weak surface mass integral N_i N_j/(rho g), with g=9.81 m/s^2. Lateral acoustic CPML includes the matched tangential x/y CFS stretch. Check the logged exact exposed-face count. Curved/cut-cell water surfaces, bulk gravity, and moving meshes are not provided by this option."),
          p("acoustic_pressure_formulation", "Legacy", "Yes/no", "yes → velocity_potential; no → none", "Backward-compatibility alias used only when acoustic_formulation is absent. New inputs should not use it."),
          p("acoustic_interface_coupling", "Optional", "Yes/no", "Default: yes", "Apply acoustic–elastic normal traction/acceleration coupling on detected interface fragments."),
          p("acoustic_interface_method", "Optional", "String (case-insensitive)", "Default: field-specific-constraints for coupled scalar acoustics; shared-geometric-constraints | field-specific-constraints | local-matching", "When omitted, uses shared-geometric-constraints if acoustics, scalar acoustics or interface coupling is disabled. Explicit unsupported choices still abort. Coupled VT topography must explicitly select shared-geometric-constraints. shared-geometric-constraints preserves shared geometric hanging-node constraints for acoustic and solid fields. Experimental field-specific-constraints uses separate acoustic, solid and pore-pressure maps, retaining only exact same-field Hex8 trace constraints; unsupported stencils abort. On regular Hex8 and FEM-topography meshes the pore-pressure map uses final active BrickUP cells, excluding elastic-only CPML cells, independently of this acoustic method and even with include_acoustic=no. VT cut-cell pressure retains geometric constraints with a startup notice. Experimental local-matching refines coarser leaves until every detected water–solid face has equal-size neighbors, then uses geometric acoustic/solid constraints. Generate Phase A with this setting; imported mismatched meshes abort. Both alternatives require include_acoustic=yes, acoustic_interface_coupling=yes, pressure or velocity_potential, and full Hex8 solids. Both support buildings, eligible nonlinear constitutive models and BrickUP poroelasticity, subject to existing module restrictions and case-specific calibration. Both use actual remote solid traces for mechanical quadrature and open-pore pressure loads. Field-specific constraints support FEM topography with include_topography=yes and computation_method=fem, excluding air cells and owner-confirmed carved voids. VT cut-tetrahedron topography is unsupported. Local matching requires include_topography=no and reapplies building carving during matching. They leave face_mass/point quadrature, face filters, gravity and CPML independently selectable. Matching can spread refinement and raise memory/runtime cost; it does not guarantee CPML stability or convergence."),
          p("acoustic_interface_projection", "Optional", "String", "Default: point; point | face_mass", "Choose one-point mechanical interface coupling or experimental weak face-mass projection. Open-pore hydraulic capacity and pressure targets always use four-point (2x2 Gauss) quadrature independently of this setting. All three constraint methods resolve actual remote solid vertices for face_mass, keeping four-point integration identical on local and remote faces. The shared-geometric point path is unchanged for sealed coupling; open_pore uses actual solid traces even with point projection so hydraulic weights reach the correct solid nodes. legacy/one_point and mass/weak/mass_projected are accepted aliases."),
          p("acoustic_interface_face_filter", "Optional", "String", "Default: all; all | bathymetry | bottom_only", "With acoustic interface coupling enabled, all couples every detected face using mesh-face normals. Explicit bathymetry couples only depth-positive Z+ bottom faces using reconstructed normals and surface-area weights; it omits side and Z− faces. Explicit bottom_only keeps only Z+ faces with vertical normals. Aliases: default and none select all; reconstructed_bathymetry and bathymetry_normal select bathymetry; bottom and z_plus select bottom_only. The default alias follows all. This geometry choice does not establish stability or waveform accuracy for every model."),
          p("acoustic_interface_max_level", "Optional", "Integer octree level", "Default: -1; -1 to PIXELLEVEL", "Maximum level for extra refinement of mixed acoustic/elastic octants; -1 disables this refinement rule."),
          p("acoustic_mesh_refinement_factor", "Optional", "Float", "Default: 1; >= 1", "Additional refinement multiplier for fully acoustic elements during fresh mesh generation."),
          p("acoustic_interface_diagnostics_file", "Optional", "File path", "Default: unset", "Write one diagnostic row per local acoustic–solid interface fragment. The file may be very large."),
          p("acoustic_station_gradient", "Optional", "String (case-insensitive)", "Default: linear; linear | element", "Reconstruct acoustic station motion from a gradient recovered at the receiver position using adjacent Cartesian water cells of the same Vp and density. Uses one-sided recovery when only one side is available; retains the original element gradient inside CPML or when neither side is suitable. Supports pressure and velocity_potential, ASCII and binary station displacement, and host/CUDA-aware MPI. element selects the historical single-element derivative. Changes station displacement, velocity and acceleration only; pressure interpolation, plane output and physical time stepping are unchanged. Does not repair interface-constraint or wavefield errors. Compare against a physical reference at fixed receiver positions."),
          p("acoustic_quality_factor", "Optional", "Float (dimensionless)", "Default: 0 (lossless)", "Omitting the parameter or setting Q<=0 disables interior acoustic damping. The default 0 was verified in the documented flat-water CPS comparisons and MPI/backend consistency tests, within their recorded configurations and durations. Positive Q enables a numerical, high-frequency-selective regularizer in interior acoustic elements."),
        ],
      },
      {
        id: "poroelastic",
        title: "Poroelastic u-p formulation",
        parameters: [
          p("include_poroelastic", "Optional", "Yes/no", "Default: no", "Enable the BrickUP-style displacement–pore-pressure formulation for saturated shallow solids. The skeleton may be linear or use an eligible nonlinear constitutive model."),
          p("poroelastic_formulation", "Optional", "String", "Default: brick_up; brick_up only", "Poroelastic element formulation. brickup is accepted as an alias."),
          p("poroelastic_depth_max_m", "Optional", "Float (m)", "Default: 30; >= 0", "Maximum depth below local ground (or the explicit seabed in either full or sqd) at which eligible solid, non-PML elements are classified as poroelastic. The ground datum includes surface_shift_m when buildings are enabled, excluding above-ground building cells; foundation cells follow the groundwater and depth limits. Nonlinear elements below the limit remain single-phase nonlinear solids: constitutive stress is assembled directly as total stress, without a pore-pressure degree of freedom or Biot feedback. This supports an intentional dry/drained or calibrated total-stress region but cannot reproduce coupled liquefaction below the limit. Startup reports the nonlinear-to-BrickUP intersection."),
          p("surface_groundwatertable", "Optional", "Yes/no", "Default: no", "When no, saturation starts at the local ground surface. When yes, read a spatial groundwater-depth grid."),
          p("groundwaterlevel_file", "Optional", "File path", "Default: groundwaterlevel.in in the 3-D material directory", "Groundwater depth grid used when surface_groundwatertable=yes. In explicit seabed mode, dry cells use ground-relative depths and wet cells instead use the water surface. Bare filenames resolve relative to 3D_velocity_model_directory."),
          p("poroelastic_surface_pressure_bc", "Optional", "String", "Default: undrained; undrained | drained", "Pore-pressure condition at the saturation surface. drained applies a weak zero-pressure projection assembled in the same local compliance units as the element storage capacity."),
          p("poroelastic_surface_pressure_penalty_factor", "Optional", "Float", "Default: 1; >= 0", "Scale for the weak drained-surface pressure projection; zero disables that projection."),
          p("poroelastic_elastic_interface", "Optional", "String", "Default: no_flow; no_flow only", "Poroelastic–elastic hydraulic condition. Mechanical displacement remains shared while Darcy flux stops at the interface. Regular Hex8 and FEM-topography models automatically build an independent pore-pressure constraint map from final active BrickUP cells, even with include_acoustic=no. Coarse dry/CPML cells do not constrain pore pressure; active coarse BrickUP cells retain exact Hex8 traces. Unsupported stencils abort. VT cut-cell pressure retains geometric constraints with a startup notice."),
          p("poroelastic_undrained_fast_path", "Optional", "Yes/no", "Default: yes for linear analysis; automatically inactive with nonlinear analysis", "Fold alpha²M into linear volumetric stiffness for eligible undrained, zero-permeability cases while still updating diagnostic pore pressure. Nonlinear skeletons always retain explicit Qp coupling so their constitutive volumetric and effective-stress response is not bypassed."),
          p("poroelastic_porosity", "Optional", "Float", "Default: 0.4; 0 < value <= 1", "Connected pore-volume fraction."),
          p("poroelastic_biot_alpha", "Optional", "Float", "Default: 1; 0 to 1", "Biot effective-stress coefficient. For alpha<1, Hercules derives the element grain modulus from K_s=K_d/(1-alpha), computes 1/M=(alpha-n)/K_s+n/K_f, and assembles the heterogeneous nodal pressure capacity as the sum of integral N/M over adjacent elements."),
          p("poroelastic_fluid_bulk_modulus_pa", "Optional", "Float (Pa)", "Default: 2.2e9; > 0", "Pore-fluid bulk modulus."),
          p("poroelastic_fluid_density_kgm3", "Optional", "Float (kg/m3)", "Default: 1000; > 0", "Pore-fluid density."),
          p("poroelastic_fluid_viscosity_pa_s", "Optional", "Float (Pa·s)", "Default: 1e-3; > 0", "Dynamic viscosity used in Darcy diffusion."),
          p("poroelastic_permeability_x_m2", "Optional", "Float (m2)", "Default: 0; >= 0", "Intrinsic permeability along solver X (grid north at zero domain azimuth)."),
          p("poroelastic_permeability_y_m2", "Optional", "Float (m2)", "Default: 0; >= 0", "Intrinsic permeability along solver Y (grid east at zero domain azimuth)."),
          p("poroelastic_permeability_z_m2", "Optional", "Float (m2)", "Default: 0; >= 0", "Intrinsic permeability along solver Z (positive depth)."),
          p("acoustic_poroelastic_interface", "Optional", "String", "Default: no_flow; no_flow | open_pore", "Hydraulic condition at acoustic–poroelastic contact. open_pore requires both formulations and weakly enforces pressure continuity using four-point (2x2 Gauss) integration of acoustic pressure against the solid trace, with the same quadrature for boundary capacity. This is independent of acoustic_interface_projection; the existing material inputs and penalty factor are retained."),
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
          p("number_output_stations", "Optional", "Integer", "Default: 0; >= 0", "Number of station rows Hercules reads from output_stations. Binary mode requires 1..INT_MAX/3 rows and exactly one owner for every station."),
          p("station_binary", "Optional", "Yes/no", "Default: no", "Select one binary displacement file for output_stations. Requires number_output_stations>0 and station_output_binary_file when yes. Uses existing station coordinate conversion and displacement interpolation, including topography and acoustic reconstruction. Ignores ASCII print_station_* switches, nonlinear columns, output_stations_directory, and station_mediums_output_file. No per-station files or metadata sidecar are created. Independent of structured DRM boxes."),
          p("station_ascii_buffer_kib", "Optional", "Integer KiB (1024 bytes)", "Default: 64; 1..2097151 (INT_MAX/1024)", "CPU stdio buffer per locally owned ASCII station, shared by all its columns. Memory per rank is value × 1024 × local station count bytes plus stream/container overhead. Writes when full and flushes at normal close, preserving every sample. Larger buffers reduce write frequency but delay live-file visibility. Always validated; unused with station_binary=yes. No ASCII buffers are allocated for binary output or ranks without stations."),
          p("station_output_binary_file", "Conditional", "File path (1 to 240 bytes)", "Required when station_binary=yes; ignored otherwise", "Headerless native-endian float64 [snapshot][station][XYZ] displacement in metres, preserving output_stations row order including duplicates. Parent directories are created; each run truncates the file, including restarts. One collective MPI-IO write per saved station snapshot. Paths resolve from the run directory without tilde expansion; use a path distinct from parameters.in and all DRM box output files."),
          p("output_stations_print_rate", "Conditional", "Integer time steps", "Required when number_output_stations>0; > 0", "Station sampling interval in both ASCII and binary modes. Saves steps divisible by this rate, including step zero in a fresh run."),
          p("output_stations_directory", "Optional", "Directory path", "Default: outputfiles/stations", "Directory containing ASCII station.N histories; ignored when station_binary=yes."),
          p("output_stations", "Conditional", "Table: N rows × 3 floats", "latitude(deg) longitude(deg) depth(m)", "Station coordinates, required when number_output_stations>0. Note the latitude–longitude order. Binary mode requires finite latitude in [-90,90], longitude in [-180,180], nonnegative depth, and converted coordinates inside the domain."),
          p("station_mediums_output_file", "Optional", "File path", "Default: output_stations_directory/station_mediums.txt", "Station ownership, medium, element, local-coordinate, and material metadata in ASCII mode, written after solver setup finalizes active BrickUP and CPML membership and refreshed after ramped geostatic initialization. Elastic-only CPML receivers have pressure undefined even when adjacent active soil shares their nodes. Ignored when station_binary=yes."),
          p("print_station_velocities", "Optional", "Yes/no", "Default: no", "Append three-component station velocity."),
          p("print_station_accelerations", "Optional", "Yes/no", "Default: no", "Append three-component acceleration; the current output path also provides velocity when acceleration is requested."),
          p("print_station_pressures", "Optional", "Yes/no", "Default: no", "Append acoustic pressure or poroelastic excess pore pressure in Pa. Requires an active scalar acoustic formulation or poroelastic physics. Poroelastic receivers must be saturated and belong to the final active BrickUP map; elastic-only CPML receivers output nan even when adjacent active soil shares their nodes. Skeleton motion is unchanged."),
          p("print_station_pore_pressure_ratios", "Optional", "Yes/no", "Default: no; requires include_poroelastic=yes and include_nonlinear_analysis=yes", "Append conventional experimental pore-pressure ratio Ru(-) = excess BrickUP pressure / positive initial vertical effective stress. Requires either approximate_geostatic_state=yes with zero geostatic loading/cushion times, or approximate_geostatic_state=no with a positive geostatic loading or cushion interval. A ramped geostatic run refreshes the denominator from committed skeleton stress before dynamics. Enable print_station_pressures separately for Pressure(Pa); the final denominator is recorded as sigma_v0_eff_pa in station_mediums.txt."),
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
          p("output_planes", "Conditional", "Table: N rows × 9 values", "lat lon depth d_strike n_strike d_dip n_dip strike(deg) dip(deg)", "Plane origins and regular sampling geometry, required when number_output_planes>0. Origins use geographic latitude/longitude; strike is clockwise from solver +X. Output vector/tensor components use solver axes. PROJ writes coordinate_frame.txt in the output directory."),
          p("reconstruct_acoustic_plane_displacements", "Optional", "Yes/no", "Default: no", "Integrate -grad(p)/rho at acoustic plane points. Requires include_acoustic=yes and pressure or velocity_potential formulation; the latter first forms p = -(phi[n] - phi[n-1])/dt, matching stations. Solid points retain elastic displacement interpolation."),
        ],
      },
      {
        id: "checkpoints-mesh-selection",
        title: "Checkpoints and selected mesh coordinates",
        parameters: [
          p("use_checkpoint", "Optional", "Boolean integer", "Default: 0; 0 | 1", "Restart from checkpoint.in in checkpoint_path. Must be 0 with CPML because the checkpoint format omits absorber histories."),
          p("checkpointing_rate", "Optional", "Integer time steps", "Default: 0; >= 0", "Checkpoint interval. Zero disables checkpoint generation and is required with CPML."),
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
          p("maximum_octant_level", "Conditional", "Integer octree level", "Required with topography; 0 to 30", String.raw`Maximum refinement level allowed by the topography mesher. Choose the smallest integer \(N\) satisfying \(\displaystyle \frac{\max(\mathtt{region\_length\_north\_m},\mathtt{region\_length\_east\_m},\mathtt{region\_depth\_deep\_m}-\mathtt{region\_depth\_shallow\_m})}{2^N}\le\frac{\max(V_{S,\min}^{\mathrm{model}},\mathtt{simulation\_shear\_velocity\_min})}{\mathtt{simulation\_wave\_max\_freq\_hz}\,\mathtt{simulation\_node\_per\_wavelength}}\), then set \(\mathtt{maximum\_octant\_level}=N\). Here \(V_{S,\min}^{\mathrm{model}}\) is derived from the active solid material inputs. Other active physics may require a finer cap.`),
          p("computation_method", "Conditional", "String", "Required with topography; vt | fem", "Topographic element treatment: volume tetrahedra (vt) or FEM."),
          p("topographybase_zcoord", "Conditional", "Float (m)", "Required with topography; > 0", "Reference depth coordinate below which the topographic surface is embedded."),
          p("type_of_etree", "Conditional", "String", "Required with topography; sqd | full; no default", "Legacy material depth reference (unless explicit seabed geometry selects 3D_velocity_model_depth_reference): sqd subtracts the local surface; full retains fixed model Z. For 3D_velocity_model=yes, sqd z.in entries are depths below each column's own surface and columns are sampled at equal physical elevation using the actual query-point topography; full z.in entries use fixed model Z with no terrain shift. Both are in metres, positive downward. Changing modes requires matching material depths. computation_method separately selects VT/FEM."),
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
        introduction: "With include_nonlinear_analysis=yes, eligible solid elements perform history-dependent constitutive updates, including yield evaluation, plastic-strain evolution, hardening or nonlinear-backbone state, and nonlinear station diagnostics. Parameters marked Conditional are required when the feature is enabled; the Manzari-Dafalias integration scheme and recovery policy are optional model-specific overrides. PressureDependMultiYield03 (PDMY03) uses its standard automatic nested surfaces and a separate 22-column property layout. RIVASand uses the frozen Ottawa F65 RIVA-Sand calibration and a two-column selection/relative-density table; RIVASandCustom uses the same production kernel with a 12-column table exposing ten existing physical constitutive parameters. RIVASand02 is a distinct research successor with its own reference calibration and two-column selection/relative-density boundary; RIVASandCustom02 exposes the same ten physical inputs through that successor kernel.",
        parameters: [
          p("include_nonlinear_analysis", "Optional", "Yes/no", "Default: no", "Enable nonlinear solid analysis, constitutive updates, plastic-strain evolution, and nonlinear station diagnostics in eligible elements. Acoustic elements are always excluded. Istanbul solid elements use the legacy nl_id classification; other solid elements ignore nl_id and use the Vs range covered by material_properties_list."),
          p("geostatic_loading_time_sec", "Conditional", "Float (s)", ">= 0; required for nonlinear analysis", "Duration used to establish geostatic loading. Nonlinear BrickUP analysis rejects a positive ramp with exactly zero permeability because that ramp is undrained and can transfer overburden into pore pressure instead of establishing the intended effective stress."),
          p("geostatic_cushion_time_sec", "Conditional", "Float (s)", ">= 0; required for nonlinear analysis", "Post-loading cushion duration. Loading plus cushion may not exceed the simulation end time."),
          p("approximate_geostatic_state", "Conditional", "Yes/no", "Required for nonlinear analysis", "Initialize an approximate K0 at-rest state from the integrated 1-D density overburden. Use yes for the documented MohrCoulomb, DruckerPrager, ManzariDafalias, PDMY03, RIVASand, RIVASandCustom, RIVASand02, and RIVASandCustom02 pressure-dependent examples when no calibrated ramp or restart supplies the in-situ state. Active BrickUP elements use effective overburden after subtracting alpha*rho_f*g*z_w; the pressure DOF starts at zero and represents excess pore pressure. It must be no when either geostatic time is positive."),
          p("error_tolerance", "Conditional", "Float", "> 0; required for nonlinear analysis", "Constitutive convergence/error tolerance. For ManzariDafalias Modified Euler this is the maximum absolute yield-surface residual in Pa after plastic stress correction; the elastic predictor separately uses the one-sided condition f_trial <= error_tolerance. Forward Euler reports its final yield drift but does not correct or reject the first-order state against this bound, matching OpenSees scheme 5."),
          p("backbone_errTol", "Conditional", "Float", "> 0; required for nonlinear analysis", "Tolerance used while matching nonlinear backbone behavior. ManzariDafalias Modified Euler reuses it as the dimensionless local error target. Forward Euler has no embedded error estimator and does not use it. PDMY03 and all four RIVA-Sand identifiers use deterministic fixed update rules and likewise do not use this tolerance."),
          p("material_model", "Conditional", "String", "linear | vonMises_ep | vonMises_fa | vonMises_faM | vonMises_baE | vonMises_baH | vonMises_GQH | vonMises_MKZ | vonMises_RO | MohrCoulomb | DruckerPrager | ManzariDafalias | PressureDependMultiYield03 (alias PDMY03) | RIVASand | RIVASandCustom | RIVASand02 | RIVASandCustom02", "Constitutive family: elastic diagnostic, von Mises/cyclic plasticity, nonlinear backbone response, frictional plasticity, the pressure-dependent critical-state Manzari–Dafalias sand model, the nested-surface PressureDependMultiYield03 model, the frozen-calibration production RIVA-Sand model, its custom-material interface, or the intermediate-bias successor and its custom interface. RIVASand02 and RIVASandCustom02 do not change either original production identifier. The former RIVASandIntermediateBiasResearch input name is retired. ManzariDafalias, PDMY03, and all four RIVA-Sand identifiers currently require the GPU backend for element-scale runs."),
          p("material_plasticity_type", "Conditional", "String", "rate_dependent | rate_independent", "Select the integration family for plastic constitutive response."),
          p("material_properties_count", "Conditional", "Integer", ">= 1", "Number of rows in material_properties_list."),
          p("material_properties_list", "Conditional", "Model-specific table", "Most models: N × 16; ManzariDafalias: N × 17; PDMY03: N × 22; RIVASand: N × 2 or 3; RIVASandCustom: N × 13; RIVASand02: N × 2 or 3; RIVASandCustom02: N × 13", "Rows are ordered by increasing, distinct Vs_limit; the first and last limits define the inclusive nonlinear-element range and active continuous values are linearly interpolated by Vs. Existing models use: Vs_limit alpha_cohesion K_phi strain_rate sensitivity hardening beta_dilatancy gamma0 psi m theta1 theta2 theta3 theta4 theta5 phi_RO. With vonMises_baE, alpha_cohesion is Su in Pa. The present vonMises_baH and vonMises_MKZ paths derive active backbone quantities from Vs and overburden, so their documented table rows primarily select the Vs interval. ManzariDafalias uses: Vs_limit G0 e_init Mc c lambda_c e0 ksi P_atm m h0 ch nb A0 nd z_max cz. PDMY03 uses: Vs_limit Gref Kref friction_angle peak_shear_strain ref_pressure pressure_dependence phase_transform_angle loading_type ca cb cc cd ce da db dc surface_count liquefaction1 liquefaction2 atmospheric_pressure cohesion. PDMY03 stress values use Pa and angles use degrees; liquefaction1 and liquefaction2 are dimensionless OpenSees coefficients and are not stress-unit conversions. loading_type is 0 (triaxial compression), 1 (triaxial extension), or 2 (direct simple shear), and surface_count is an integer from 2 through 40. Its categorical loading_type and surface_count come from the nearest Vs row rather than interpolation. RIVASand and RIVASand02 use Vs_limit Dr, optionally followed by G0; their other physical fields keep reference values. Custom names require Vs_limit Dr G0 M kd h m zeta e_max e_min Q R n_G. Each row occupies one line; all rows must have equal width. Missing G0 or G0=0 selects automatic stiffness from the actual mesh bulk density and Vs, not Vs_limit. Positive dimensionless G0 overrides it. Do not mix automatic and positive G0 rows. Continuous fields and explicit G0 are linearly interpolated by Vs. Automatic stiffness is calibrated per integration point from its admitted effective geostatic stress at dynamic activation, including initialized stiffness factors and pressure guards, then frozen and checkpointed. Require finite inputs, 0<=Dr<=1, G0>=0, M,kd,h,Q>0, m,zeta,R>=0, e_max>e_min>=0 containing e=0.601, and 0<=n_G<=1. Reference G0=483.48301127222084 preserves the previous 48.976829 MPa scale at 101.3 kPa. Bulk stiffness scales with shear stiffness at fixed nu=0.30, not mesh Vp. Reference hashes are 9585f0c155c9444885c5115e7753a1a5d97c783e2c685938a49a65435c9e8f83 for the original and 16c8567fd2370c562899b9ec00cc2fc313be137f239cde5d0e22cc1cd3eee800 for 02. Original models retain checkpoint marker 4. RS02 markers are 10 + 3*BiasVolume + (reversaltype-1), with default marker 15. All older RS02 markers through 9 must be rerun because the state now contains 139 values. Changing either selector on restart is rejected. Restart requires unchanged mesh and material inputs. Material-specific void-ratio bounds control both initial and reference density in 02. BiasVolume is selected separately with rivasand02_biasvolume; h_p0 and reversalLatch are not exposed. Acoustic elements remain excluded even when the range starts at zero; CPML elements also remain linear."),
          p("no_substeps", "Conditional", "Integer", ">= 1; required for nonlinear analysis", "Constitutive substep control. For ManzariDafalias Modified Euler it normally caps accepted adaptive substeps per global increment; rejected error-control trials are also finitely bounded. manzari_dafalias_recovery=unlimited_substeps or opensees disables those work caps, but no_substeps still bounds plastic stress-correction iterations. For Forward Euler with manzari_dafalias_forward_euler_substeps=fixed, it is the exact number of equal subdivisions over every plastic remainder. With demand, it is the upper bound on an MD-specific stress-ratio/deviatoric-strain request. no_substeps=1 reproduces OpenSees's original single-update scheme 5 under either policy. For PDMY03 it caps the separate deterministic surface/strain-based requested subincrements; a larger request is clamped and reported rather than treated as a constitutive failure. For RIVASand, RIVASandCustom, RIVASand02, and RIVASandCustom02 it is the exact number of equal fixed Forward-Euler substeps per nonzero Hercules increment; reversal detection still runs once per Hercules increment, and an exact zero increment preserves history without subdivision. Start with 1 and increase only after a timestep/substep convergence study."),
          p("rivasand02_reversaltype", "Optional", "Integer (dimensionless)", "Default: 3; 1 | 2 | 3; RIVASand02 and RIVASandCustom02 only", "Choose the RS02 reversal reference once for the simulation: 1 compares the current increment with the preceding valid increment; 2 compares it with accumulated loading-branch strain; 3 uses the full-increment elastic stress-ratio predictor from the committed mechanical state. All three check once per Hercules increment before material substeps. Omission selects type 3 for both standard and Custom02 inputs; invalid values, duplicate entries, and use with other models are rejected. The material-point runner accepts the same key. The retired reversalType and rivasand02_reversalType keys are rejected. Keep both reversal type and rivasand02_biasvolume unchanged on restart; all old 138-value RS02 checkpoints must be rerun. Type 3 ports the unguarded OpenSees_RIVA 355d8f042 rule: a negative (alpha_trial-alpha0):(alpha_trial-alpha) triggers at most one event in the first material substep. Before cyclic activation it retains the strain-based initialization rule. reversalGuard and reversalLatch are not exposed. This choice does not establish timestep independence."),
          p("rivasand02_biasvolume", "Optional", "Integer (dimensionless)", "Default: 1; 0 | 1 | 2; RIVASand02 and RIVASandCustom02 only", "Select the OpenSees_RIVA 355d8f042 BiasVolume behavior: 0 retains the original bias-driven reversible-volume contribution; 1 disables that contribution while retaining the separate phase-transformation volume; 2 retains it with the upstream gated correction to its contractive mean part. Mode 2 uses fixed-bias, reference-pressure, effective-pressure, and accumulated plastic-activity gates; it is an optional research choice, not a universally validated stability fix. Both field and material-point inputs accept this scalar; it is not a material-row column. Invalid values, duplicate entries, and use with other models are rejected. The 139th state value is field_bias_mean_activity; it stays zero in modes 0/1 and updates once after an accepted host increment in mode 2. Calibration fingerprint 16c8567fd2370c562899b9ec00cc2fc313be137f239cde5d0e22cc1cd3eee800 identifies the 257-value upstream reference record; runtime selector overrides are identified separately. Checkpoint marker = 10 + 3*mode + (reversaltype-1), so the defaults use marker 15. Changing either selector on restart, or loading an older 138-value RS02 checkpoint, is rejected. To reproduce earlier response rules from fresh initialization, explicitly set mode 0 and the earlier reversal type. BrickUP remains responsible for fluid pressure and total-stress assembly; this control does not disable its coupling."),
          p("manzari_dafalias_integration_scheme", "Optional", "String", "Default: modified_euler; modified_euler | forward_euler; ManzariDafalias only", "Select the Manzari-Dafalias explicit integration algorithm. modified_euler uses two stress-increment estimates, adaptive local-error control, stress correction, and the configured recovery policy. forward_euler applies the OpenSees scheme-5 first-order constitutive stage using the fixed or demand subdivision policy, reports the resulting absolute yield drift, and does not use backbone_errTol or manzari_dafalias_recovery. Modified Euler remains the production default because Forward Euler is conditionally stable and more substep-size sensitive."),
          p("manzari_dafalias_forward_euler_substeps", "Optional", "String", "Default: fixed; fixed | demand; ManzariDafalias forward_euler only", "Select how Forward Euler divides the plastic remainder after any elastic/yield intersection. fixed preserves the established behavior: no_substeps is the exact number of equal parts, which supports reproducible integration-sensitivity studies. demand treats no_substeps as a cap and requests work from an MD-specific trial stress-ratio excursion threshold of 0.05 and octahedral engineering deviatoric-strain threshold of 1.0e-5. The demand rule is not the PDMY03 criterion and is an opt-in performance heuristic; compare it with fixed and Modified Euler for the calibrated loading path before production use."),
          p("manzari_dafalias_recovery", "Optional", "String", "Default: strict; strict | force_accept | unlimited_substeps | opensees; ManzariDafalias modified_euler only", "Select the recovery policy after adaptive Modified-Euler trials are rejected. strict retains the finite no_substeps and attempt caps and never accepts a state above backbone_errTol. force_accept retains those caps but, at the minimum local fraction, commits a finite stress-ratio-capped admissible state in the style of OpenSees. unlimited_substeps disables the accepted-substep and rejected-attempt caps but retains the strict local-error test. opensees combines force acceptance, uncapped substeps/attempts, and OpenSees's fixed normalized minimum fraction of 1e-6. Non-strict modes are robustness fallbacks: they can increase cost or accept a state whose embedded error exceeds backbone_errTol, so production use must be reported and independently checked. Forward Euler requires the default strict value because it has no rejected-trial recovery loop."),
          p("stiff_damp", "Conditional", "Float", "Required for nonlinear analysis", "Stiffness-proportional damping coefficient used by the nonlinear path. Set 0 to disable this contribution; keep freq_stiff_damp positive."),
          p("freq_stiff_damp", "Conditional", "Float (Hz)", "Required for nonlinear analysis", "Reference frequency associated with nonlinear stiffness damping."),
          p("tension_cutoff", "Conditional", "Yes/no", "Required for nonlinear analysis", "Remove tensile stress. yes is accepted only for MohrCoulomb and DruckerPrager."),
        ],
      },
      {
        id: "nonlinear-material-point-test",
        title: "Standalone material-point test controls",
        introduction: "These keys belong only in an input passed to psolve --nonlinear-material-point. They are not part of a production parameters.in simulation. The test also reads material_model, no_substeps, error_tolerance, backbone_errTol, tension_cutoff, and the optional RS02 rivasand02_reversaltype and rivasand02_biasvolume from the same file.",
        parameters: [
          p("material_point_delta_time_sec", "Conditional", "Float (s)", "> 0; material-point subcommand only", "Time increment associated with each prescribed material-point strain row."),
          p("material_point_steps", "Conditional", "Integer", "> 0; material-point subcommand only", "Number of rows in material_point_strain_history."),
          p("material_point_constants", "Conditional", "Model-specific list", "24 floats normally; 17 for ManzariDafalias; 21 for PDMY03; 1 or 2 for standard RIVA-Sand; 12 for Custom names", "Existing fixtures use the documented direct nlconstants_t order. ManzariDafalias uses G0 nu e_init Mc c lambda_c e0 ksi P_atm m h0 ch nb A0 nd z_max cz. PDMY03 omits the production Vs_limit and uses Gref Kref friction_angle peak_shear_strain ref_pressure pressure_dependence phase_transform_angle loading_type ca cb cc cd ce da db dc surface_count liquefaction1 liquefaction2 atmospheric_pressure cohesion. RIVA-Sand standard names use Dr [G0], while Custom names require Dr G0 M kd h m zeta e_max e_min Q R n_G, omitting Vs_limit. Automatic stiffness (missing G0 or G0=0) requires material_point_riva_density and material_point_riva_vs because this driver has no mesh; it uses the supplied effective initial stress. The same bounds as material_properties_list apply. Stress-like entries use Pa. PDMY03 liquefaction1 and liquefaction2 are dimensionless. ManzariDafalias angles use radians; PDMY03 angles use degrees."),
          p("material_point_initial_strain", "Conditional", "6 floats", "epsilon_xx epsilon_yy epsilon_zz epsilon_xy epsilon_yz epsilon_xz; ManzariDafalias, PDMY03, RIVASand, or RIVASandCustom", "Optional committed initial tensor strain for a seeded pressure-dependent replay. At its elastic-to-plastic transition PDMY03 clamps tensile mean stress to zero and accepts an exactly zero state, matching OpenSees; production soil models should nevertheless begin from a physical in-situ effective stress. If any material-point initial-state array is present, all five are required; alpha, alpha_in, and fabric are zero placeholders for PDMY03 and both RIVA-Sand interfaces. Off-diagonal entries are tensorial shear strains."),
          p("material_point_initial_stress", "Conditional", "6 floats (Pa)", "sigma_xx sigma_yy sigma_zz sigma_xy sigma_yz sigma_xz; ManzariDafalias, PDMY03, RIVASand, or RIVASandCustom", "Optional committed initial effective-stress tensor for a seeded pressure-dependent replay. Hercules uses tension-positive stress. Both RIVA-Sand interfaces require this equilibrated seed and initialize their full native history from it."),
          p("material_point_initial_alpha", "Conditional", "6 floats", "alpha_xx alpha_yy alpha_zz alpha_xy alpha_yz alpha_xz; ManzariDafalias, PDMY03, RIVASand, or RIVASandCustom", "Optional committed ManzariDafalias backstress-ratio tensor. It must be supplied with the other four initial-state arrays. PDMY03 and both RIVA-Sand interfaces accept and ignore a zero placeholder because they initialize native history from stress."),
          p("material_point_initial_alpha_in", "Conditional", "6 floats", "alpha_in_xx alpha_in_yy alpha_in_zz alpha_xy alpha_yz alpha_xz; ManzariDafalias, PDMY03, RIVASand, or RIVASandCustom", "Optional committed ManzariDafalias reversal-memory tensor. It must be supplied with the other four initial-state arrays. PDMY03 and both RIVA-Sand interfaces accept and ignore a zero placeholder."),
          p("material_point_initial_fabric", "Conditional", "6 floats", "fabric_xx fabric_yy fabric_zz fabric_xy fabric_yz fabric_xz; ManzariDafalias, PDMY03, RIVASand, or RIVASandCustom", "Optional committed ManzariDafalias fabric tensor. It must be supplied with the other four initial-state arrays. PDMY03 and both RIVA-Sand interfaces accept and ignore a zero placeholder."),
          p("material_point_strain_history", "Conditional", "Table: N rows × 6 floats", "epsilon_xx epsilon_yy epsilon_zz epsilon_xy epsilon_yz epsilon_xz", "Prescribed tensor-strain history for the material-point subcommand. Off-diagonal entries are tensorial, not engineering, shear strains."),
          p("material_point_riva_density", "Conditional", "Float (kg/m³)", "> 0; required for automatic RIVA-Sand material-point stiffness", "Bulk mass density, not buoyant density. Used with material_point_riva_vs and supplied effective initial stress; this is a material-point subcommand input, not a mesh-run parameter."),
          p("material_point_riva_vs", "Conditional", "Float (m/s)", "> 0; required for automatic RIVA-Sand material-point stiffness", "Initial small-strain shear-wave velocity. The driver calibrates the initialized elastic shear modulus to density times Vs squared once before CPU/GPU replay."),
          p("material_point_compare_rtol", "Optional", "Float", "Default: 1e-6; >= 0; material-point subcommand only", "Relative CPU/GPU state-comparison tolerance."),
          p("material_point_compare_atol", "Optional", "Float", "Default: 1e-9; >= 0; material-point subcommand only", "Absolute CPU/GPU state-comparison tolerance."),
          p("material_point_output_prefix", "Optional", "Path prefix", "Default: material-point; material-point subcommand only", "Prefix for the generated .cpu.csv and .gpu.csv state histories."),
          p("material_point_brickup", "Optional", "Yes/no", "Default: no; material-point subcommand only", "Enable the sealed uniform one-element BrickUP pressure regression alongside a ManzariDafalias, PressureDependMultiYield03, RIVASand, or RIVASandCustom strain replay. The constitutive result is effective skeleton stress; the independent pressure history is combined only once through sigma_total=sigma_effective-alpha*p*I. The GPU side executes the production pressure-force and pressure-update kernels; this is not a production parameters.in control."),
          p("material_point_brickup_porosity", "Conditional", "Float", "0 < n <= 1; required when material_point_brickup=yes", "Fluid-filled porosity used to form the standalone sealed BrickUP storage modulus."),
          p("material_point_brickup_fluid_bulk_modulus_pa", "Conditional", "Float (Pa)", "> 0; required when material_point_brickup=yes", "Pore-fluid bulk modulus Kf for the standalone sealed BrickUP pressure update."),
          p("material_point_brickup_biot_alpha", "Conditional", "Float", "0 < alpha <= 1; required when material_point_brickup=yes", "Biot coefficient for the standalone BrickUP comparison. Alpha=1 uses M=Kf/n; alpha<1 also requires material_point_brickup_solid_bulk_modulus_pa and uses 1/M=(alpha-n)/Ks+n/Kf."),
          p("material_point_brickup_solid_bulk_modulus_pa", "Conditional", "Float (Pa)", "> 0 when material_point_brickup_biot_alpha<1; material-point subcommand only", "Solid-grain bulk modulus Ks used by the standalone arbitrary-alpha BrickUP storage regression. It is not a production parameters.in control."),
          p("material_point_initial_pore_pressure_pa", "Conditional", "Float (Pa)", "Finite; required when material_point_brickup=yes", "Committed pore pressure at the first seeded replay state. Generated drained-to-undrained comparisons use zero."),
        ],
      },
    ],
  },
  {
    slug: "parameters-drm-waves",
    title: "DRM, incident-wave, and half-space parameters",
    summary: "Single-file DRM displacement recording, Domain Reduction Method stages, incident plane waves, and the homogeneous half-space override.",
    sections: [
      {
        id: "drm",
        title: "Domain Reduction Method",
        parameters: [
          p("implement_drm", "Optional", "Yes/no", "Default: no", "Enable the three-part Domain Reduction Method workflow."),
          p("number_output_drm_boxes", "Optional", "Nonnegative integer", "Default: 0", "Number of generated DRM box rows. Zero disables generated boxes. Each row writes one drmdisplacements.N binary containing ten structured face grids and one drmcoordinates.N ASCII table in matching sample order. Independent of implement_drm, station output and planes."),
          p("output_drm_boxes_print_rate", "Conditional", "Positive integer time steps", "Required when number_output_drm_boxes > 0", "Common sampling interval for generated boxes. Snapshots are saved before advancing steps divisible by this rate, including step zero; saved interval is this rate times simulation_delta_time_sec."),
          p("output_drm_boxes_directory", "Conditional", "Directory path", "Required when number_output_drm_boxes > 0", "Created automatically. Each box writes drmdisplacements.N and drmcoordinates.N, numbered from zero in row order. Both are overwritten on startup, including restarts. Coordinates are automatic: comment headers followed by X_m Y_m Z_m Vs_m_s Vp_m_s density_kg_m3 Qs Qp. Coordinates are solver-local X/Y/down metres after rotation and building/topography shifts, preserving repeated edge/corner rows. Appended initial owning-element speeds are m/s, density kg/m³, and Qs/Qp dimensionless material/BKT inputs (nan when unavailable), not effective nonlinear or poroelastic attenuation. Position-only readers should select the first three columns. Rebuild and regenerate older mesh payloads with mesh_input_mode=0 to retain Qs/Qp. Rank zero writes the ASCII once at setup with 17 significant digits. The complete file path must contain at most 240 bytes. Paths resolve from the working directory; no tilde expansion."),
          p("output_drm_boxes", "Conditional", "N rows × 7 floats", "Required when number_output_drm_boxes > 0", "Each row: drm_offset_x drm_offset_y length_x length_y depth drm_edgesize strk. Offsets locate the top center in solver-local X/Y metres. Lengths and depth describe the exterior; interior dimensions are length_x−2h, length_y−2h, depth−h, with h=drm_edgesize. Lengths must exceed 2h, depth must exceed h, and all three dimensions must be integer multiples of positive h. Strike is finite degrees clockwise from solver +X. PROJ writes coordinate_frame.txt with the geographic mapping; binary components stay in solver axes. Each snapshot stores exterior then interior, with front (−local X), back (+local X), left (−local Y), right (+local Y), bottom; all faces include endpoints and retain shared-edge/corner samples. Rows advance local Y on front/back, local X on left/right/bottom; columns advance Z on walls and local Y on bottom. No top lid. Input rows follow the key on separate lines; blank lines and # comments are allowed. All samples must lie inside the domain and have one solid owner; topographic depth shifts follow plane conventions. No mesh refinement."),
          p("drm_output_batch_size", "Optional", "Positive integer snapshots", "Default: 16", "Snapshots buffered per collective MPI-IO write, per file, for output_drm_boxes. Set 1 to write every saved snapshot. Partial final batches are flushed at completion. For each file, local sample count × 3 × batch size must fit an MPI int count. Ignored when number_output_drm_boxes=0; does not control station binary output."),
          p("drm_directory", "Conditional", "Directory path", "Required by the DRM module", "Root directory containing part0, part1, and part2 DRM data. Although the main parser initializes a conventional path, the DRM module currently requires the key to be present."),
          p("which_drm_part", "Conditional", "String", "part0 | part1 | part2", "Select coordinate generation, response recording, or effective-force replay."),
          p("drm_edgesize", "Conditional", "Float (m)", "Required when implement_drm=yes", "Element edge size prescribed for the legacy DRM box. Generated output boxes specify their own spacing in column 6 of each output_drm_boxes row."),
          p("drm_offset_x", "Conditional", "Float (m)", "Required when implement_drm=yes", "Legacy DRM-box X/north offset: part0 adds it when exporting coordinates and part1 removes it on export. For output_drm_boxes, the row field with this name instead specifies the top-center X coordinate in metres from the solver domain origin."),
          p("drm_offset_y", "Conditional", "Float (m)", "Required when implement_drm=yes", "Legacy DRM-box Y/east offset: part0 adds it when exporting coordinates and part1 removes it on export. For output_drm_boxes, the row field with this name instead specifies the top-center Y coordinate in metres from the solver domain origin."),
          p("drm_print_rate", "Conditional", "Integer time steps", "Required when implement_drm=yes", "Legacy three-stage DRM sampling interval. Generated boxes use output_drm_boxes_print_rate; station output uses output_stations_print_rate."),
          p("part1_delta_t", "Conditional", "Float (s)", "Required when implement_drm=yes", "Time step used by the part1 data that part2 will read; normally matches the part1 simulation delta time."),
          p("part1_simulation_time", "Conditional", "Float (s)", "Required when implement_drm=yes", "Duration of the part1 data set that part2 will read."),
          p("drm_boundary", "Conditional", "List of 5 floats (m)", "Required when implement_drm=yes", "Legacy DRM boundary: min_x min_y max_x max_y depth. The independent generated recorder uses output_drm_boxes rows with top-center offsets and exterior dimensions; this standalone parameter does not define its geometry."),
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
          p("planewave_strike", "Conditional", "Float (degrees)", "Required for incident waves", "Horizontal propagation strike clockwise from solver +X, converted internally to radians."),
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
