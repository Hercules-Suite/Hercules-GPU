export type DocTable = {
  headers: string[];
  rows: string[][];
};

export type DocSection = {
  id: string;
  title: string;
  paragraphs?: string[];
  bullets?: string[];
  code?: string;
  codeLabel?: string;
  note?: string;
  table?: DocTable;
};

export type DocPageContent = {
  eyebrow: string;
  title: string;
  summary: string;
  sections: DocSection[];
};

export const docPages: Record<string, DocPageContent> = {
  "getting-started": {
    eyebrow: "Start here",
    title: "Run your first Hercules simulation",
    summary:
      "A short path from a clean checkout to a completed adaptive finite-element simulation.",
    sections: [
      {
        id: "before-you-begin",
        title: "Before you begin",
        paragraphs: [
          "Hercules is an MPI-based finite-element solver designed for large, three-dimensional wave-propagation problems. The GPU build expects a Linux HPC environment with a supported C/C++ toolchain, MPI, CUDA, and PROJ.",
        ],
        bullets: [
          "A CUDA-capable GPU for every active Phase B rank",
          "An MPI implementation compatible with your site launcher",
          "A PROJ installation for geographic coordinate transformations",
          "A batch allocation large enough for both mesh generation and solving",
        ],
        note: "Hercules intentionally stops when local GPU ranks exceed visible GPUs. Use oversubscription only for diagnostics.",
      },
      {
        id: "configure",
        title: "1. Configure the build",
        paragraphs: [
          "Copy an appropriate user configuration and update compiler, MPI, CUDA, and PROJ paths. Site-specific build choices belong in user.mk; systemdef.mk and common.mk supply shared rules.",
        ],
        codeLabel: "shell",
        code: "cp examples/user.mk-linux-mpich user.mk\n$EDITOR user.mk\nmake",
      },
      {
        id: "prepare-inputs",
        title: "2. Prepare a case",
        paragraphs: [
          "A run is configured by one parameters.in file plus the source, material, and geometry inputs it references. Start from a current benchmark_inputs case so the directory layout, units, and optional physics settings are explicit.",
        ],
        codeLabel: "typical case layout",
        code: "my-case/\n├── parameters.in\n├── sourcefiles/\n├── models/\n└── output/",
      },
      {
        id: "launch",
        title: "3. Launch",
        paragraphs: [
          "The unified workflow keeps the Phase A MPI world alive, builds or imports the adaptive mesh, and transfers leaves directly to the GPU solver ranks. It is the default scalable workflow in the current codebase.",
        ],
        codeLabel: "shell",
        code: "export BASE_IN=/path/to/my-case\nexport PHASE_A_NP=64\nexport PHASE_B_NP=8\nexport PHASE_B_OMP_NUM_THREADS=8\n\n./run_phaseA_phaseB.sh",
      },
      {
        id: "verify",
        title: "4. Verify the run",
        bullets: [
          "Confirm the selected phase mode, GPU mapping, and MPI backend in the monitor output.",
          "Check that the reported time step and total step count match the input files.",
          "Inspect station or plane outputs before scaling the case to a production allocation.",
          "Save the monitor log with benchmark results so performance choices remain reproducible.",
        ],
      },
    ],
  },
  installation: {
    eyebrow: "Start here",
    title: "Installation and build configuration",
    summary:
      "Configure Hercules for an HPC system while keeping machine-specific paths isolated in user.mk.",
    sections: [
      {
        id: "requirements",
        title: "Requirements",
        table: {
          headers: ["Component", "Purpose", "Notes"],
          rows: [
            ["C/C++ and CUDA toolchains", "Compile Octor and the GPU solver", "Use versions supported by the target CUDA toolkit"],
            ["MPI", "Distributed mesh and solver communication", "CUDA-aware support enables direct device-buffer exchange"],
            ["PROJ", "Geographic coordinate conversion", "Headers and libraries must both be visible"],
            ["zlib / bzip2", "Compressed inputs and outputs", "Usually supplied by the system environment"],
            ["Make", "Project build orchestration", "GNU Make is recommended"],
          ],
        },
      },
      {
        id: "user-mk",
        title: "Keep site settings in user.mk",
        paragraphs: [
          "user.mk is the supported customization surface. Set compiler wrappers, include directories, library directories, optimization flags, and CUDA architecture there instead of editing shared makefiles.",
        ],
        codeLabel: "user.mk — illustrative only",
        code: "CC          = mpicc\nCXX         = mpicxx\nNVCC        = nvcc\nPROJ_DIR    = /path/to/proj\nCUDA_DIR    = /path/to/cuda\n\nCFLAGS     += -O3\nCXXFLAGS   += -O3",
        note: "Exact variables differ across existing site configurations. Begin with the closest checked-in example and preserve your center's compiler-wrapper conventions.",
      },
      {
        id: "build",
        title: "Build and inspect",
        codeLabel: "shell",
        code: "make clean\nmake -j 8\n\n# Expected executables\nls octor/octor quake/forward/psolve",
        bullets: [
          "Treat missing PROJ symbols as a link-path problem before changing source.",
          "Ensure the host compiler used by nvcc matches the supported CUDA combination.",
          "Build on the target platform when cluster login and compute nodes differ materially.",
        ],
      },
    ],
  },
  architecture: {
    eyebrow: "Concepts",
    title: "How Hercules works",
    summary:
      "Hercules separates adaptive mesh construction from time integration, then joins them through a distributed octree representation.",
    sections: [
      {
        id: "components",
        title: "Core components",
        table: {
          headers: ["Component", "Responsibility"],
          rows: [
            ["Octor", "Adaptive octree mesh generation, partitioning, balancing, and mesh extraction"],
            ["Quake", "Finite-element assembly, sources, physics kernels, time stepping, communication, and output"],
            ["etree", "Persistent spatial database used by legacy mesh-storage paths"],
            ["Launcher", "Coordinates Phase A mesh work and Phase B GPU solving"],
          ],
        },
      },
      {
        id: "lifecycle",
        title: "Simulation lifecycle",
        bullets: [
          "Read parameters.in and the source, model, and output paths it defines.",
          "Refine the octree to satisfy wavelength, material, topography, and feature constraints.",
          "Balance and partition leaves across Phase A ranks.",
          "Persist or redistribute the mesh, depending on the selected phase mode.",
          "Construct local finite-element data and communication schedules on Phase B ranks.",
          "Advance elastic, acoustic, or poroelastic fields and write requested observations.",
        ],
      },
      {
        id: "adaptivity",
        title: "Why the octree matters",
        paragraphs: [
          "The octree concentrates resolution where material properties and modeled wavelengths require it. Neighbor constraints introduce hanging nodes, so Hercules maintains adjustment schedules that preserve a conforming finite-element solution across refinement boundaries and MPI partitions.",
        ],
        note: "The global domain uses a power-of-two octree dimension. The physical region of interest may occupy only part of that cube.",
      },
    ],
  },
  running: {
    eyebrow: "Run simulations",
    title: "Running Hercules at scale",
    summary:
      "Choose a phase workflow, map ranks deliberately, and make launch assumptions visible in every production log.",
    sections: [
      {
        id: "phase-modes",
        title: "Phase modes",
        table: {
          headers: ["Mode", "Mesh handoff", "Best fit"],
          rows: [
            ["unified", "Direct in-memory redistribution", "Production runs where one allocation covers both phases"],
            ["two-phase", "Sharded, MPI-IO, or single-file storage", "Independent mesh reuse or separate resource shapes"],
          ],
        },
      },
      {
        id: "unified",
        title: "Unified launch",
        codeLabel: "shell",
        code: "export HERCULES_PHASE_MODE=unified\nexport BASE_IN=/scratch/project/case\nexport PHASE_A_NP=2304\nexport PHASE_B_NP=256\nexport PHASE_B_OMP_NUM_THREADS=8\n\nsrun ./run_phaseA_phaseB.sh",
        paragraphs: [
          "Active GPU ranks divide each node's allowed CPU affinity mask into disjoint blocks. Inactive Phase A ranks use a targeted sleep-and-poll loop while Phase B advances, avoiding unnecessary CPU contention.",
        ],
      },
      {
        id: "two-phase",
        title: "Two-phase launch",
        paragraphs: [
          "Two-phase mode is useful when a mesh will be reused. The sharded backend lets Phase A ranks write concurrently; the MPI-IO backend creates one indexed archive; single retains the serialized compatibility path.",
        ],
        codeLabel: "shell",
        code: "export HERCULES_PHASE_MODE=two-phase\nexport HERCULES_MESH_IO_MODE=sharded\n./run_phaseA_phaseB.sh",
      },
    ],
  },
  "gpu-execution": {
    eyebrow: "Run simulations",
    title: "GPU and MPI execution",
    summary:
      "Keep fields resident, exchange device buffers directly when available, and move only the output data you need.",
    sections: [
      {
        id: "backend",
        title: "CUDA-aware communication",
        paragraphs: [
          "HERCULES_GPU_MPI_BACKEND selects host-staged or CUDA-aware exchanges. When Open MPI reports CUDA-buffer support, the current default enables direct device packing, MPI messaging, and GPU hanging-node adjustment for eligible elastic, acoustic, and poroelastic fields.",
        ],
        codeLabel: "shell",
        code: "export HERCULES_GPU_MPI_BACKEND=cuda_aware\nexport HERCULES_DISP_SYNC_MODE=boundary\nexport HERCULES_PRESSURE_FIELD_SYNC_MODE=resident\nexport HERCULES_PORE_PRESSURE_FIELD_SYNC_MODE=resident",
      },
      {
        id: "resident-fields",
        title: "Resident-field strategy",
        bullets: [
          "Elastic displacement history remains on the GPU; only boundary and adjustment nodes move when required.",
          "Acoustic pressure histories swap device pointers and patch only schedule or adjustment nodes.",
          "Poroelastic pore pressure uses the same sparse/resident pattern.",
          "Station interpolation runs on the GPU for supported fields and returns compact observations.",
        ],
      },
      {
        id: "transport",
        title: "Transport defaults",
        paragraphs: [
          "For large CUDA-aware Open MPI runs, the launcher prefers UCX and scalable DC transport for Phase B. cuda_copy supplies the general device-memory fallback, while cuda_ipc can accelerate intra-node transfers when the site stack handles CUDA IPC reliably.",
        ],
        note: "If CUDA IPC handles fail on a specific cluster, override only Phase B UCX transports to disable cuda_ipc while retaining cuda_copy.",
      },
      {
        id: "io",
        title: "Scalable plane output",
        paragraphs: [
          "Plane output supports legacy, hybrid GPU interpolation, and experimental distributed MPI-IO paths. MPI-IO batches packed snapshots so collective calls are amortized; a partial batch is flushed before timing reduction and close.",
        ],
        codeLabel: "shell",
        code: "export HERCULES_PLANE_IO_MODE=mpiio\nexport HERCULES_PLANE_MPIIO_BATCH=64",
      },
    ],
  },
  physics: {
    eyebrow: "Run simulations",
    title: "Physics models",
    summary:
      "Model elastic solids, acoustic fluids, saturated porous media, attenuation, topography, and absorbing boundaries in one adaptive framework.",
    sections: [
      {
        id: "models",
        title: "Supported formulations",
        table: {
          headers: ["Model", "Primary fields", "Typical use"],
          rows: [
            ["Elastic", "Three-component displacement", "Regional and site-scale seismic wave propagation"],
            ["Acoustic", "Scalar pressure", "Water layers and fluid domains"],
            ["Poroelastic u-p", "Solid displacement and pore pressure", "Saturated soils and fluid–solid interaction"],
          ],
        },
      },
      {
        id: "coupling",
        title: "Coupled interfaces",
        paragraphs: [
          "Material classification determines which element formulation advances each region. Interface terms transfer the physically appropriate traction, motion, and pressure information between neighboring acoustic, elastic, and poroelastic domains.",
        ],
        note: "Always record the physical assumptions used at coupled interfaces; a numerically stable run can still encode the wrong boundary condition.",
      },
      {
        id: "attenuation",
        title: "Attenuation and damping",
        bullets: [
          "BKT-family viscoelastic convolution kernels provide frequency-dependent attenuation.",
          "Rayleigh and mass-proportional damping remain available for compatible analyses.",
          "A diagnostic override can force the generic damping kernel to build a correctness reference.",
        ],
      },
      {
        id: "pml",
        title: "Perfectly matched layers",
        paragraphs: [
          "PML regions attenuate outgoing waves before they reach the truncation boundary. Profile order, target reflection coefficient, layer width, reference frequency, and material wave speeds jointly control the damping profile.",
        ],
      },
    ],
  },
  configuration: {
    eyebrow: "Run simulations",
    title: "Runtime configuration",
    summary:
      "Tune execution paths without recompiling. Defaults favor scalable GPU execution while explicit overrides remain available for diagnosis and comparison.",
    sections: [
      {
        id: "core",
        title: "Core execution controls",
        table: {
          headers: ["Variable", "Values", "Purpose"],
          rows: [
            ["HERCULES_PHASE_MODE", "unified | two-phase", "Select mesh/solver orchestration"],
            ["HERCULES_GPU_MPI_BACKEND", "host | cuda_aware", "Select field-exchange backend"],
            ["HERCULES_DISP_SYNC_MODE", "auto | full | boundary", "Control elastic host/device synchronization"],
            ["HERCULES_MESH_IO_MODE", "mpiio | sharded | single", "Select Phase A mesh storage"],
            ["HERCULES_PLANE_IO_MODE", "auto | hybrid | mpiio | legacy", "Select plane interpolation and writing"],
          ],
        },
      },
      {
        id: "pressure",
        title: "Pressure-field controls",
        table: {
          headers: ["Variable", "Recommended scalable setting"],
          rows: [
            ["HERCULES_PRESSURE_FORCE_SYNC_MODE", "sparse"],
            ["HERCULES_PRESSURE_FIELD_SYNC_MODE", "resident"],
            ["HERCULES_PORE_PRESSURE_FORCE_SYNC_MODE", "sparse"],
            ["HERCULES_PORE_PRESSURE_FIELD_SYNC_MODE", "resident"],
          ],
        },
      },
      {
        id: "reproducibility",
        title: "Reproducibility controls",
        paragraphs: [
          "HERCULES_REPRO_MODE selects increasingly strict deterministic paths. Reproducibility modes can disable OpenMP-sensitive optimizations and force deterministic receive handling, so use them for controlled comparisons rather than assuming production performance.",
        ],
        codeLabel: "shell",
        code: "# normal production path\nexport HERCULES_REPRO_MODE=0\n\n# deterministic diagnostic path\nexport HERCULES_REPRO_MODE=2",
      },
    ],
  },
  "parameter-file": {
    eyebrow: "Setting up simulations",
    title: "The parameters.in file",
    summary:
      "Current Hercules runs use one parameter file. It combines domain, numerical, physical, source, and output settings in a single readable configuration.",
    sections: [
      {
        id: "single-file",
        title: "One file per simulation",
        paragraphs: [
          "The solver accepts one positional input: the path to parameters.in. The file may have another name, but parameters.in is the established convention and is used throughout this documentation.",
          "Some legacy examples, scripts, and source comments still refer to separate numerical.in and physics.in files. That split is not the current psolve interface and should not be used when creating a new case.",
        ],
        codeLabel: "shell",
        code: "mpirun -np <ranks> /path/to/psolve parameters.in",
        note: "Keep parameters.in under version control with the exact source files, model data, executable revision, and monitor log used for a scientific result.",
      },
      {
        id: "syntax",
        title: "Basic syntax",
        paragraphs: [
          "A scalar parameter is written as a name, an equals sign, and one value. A parameter may also introduce a list or table whose values continue on following lines.",
        ],
        codeLabel: "parameters.in",
        code: "# Scalar values\noutput_mesh = 0\ntype_of_damping = none\nregion_origin_latitude_deg = 40.69146389\nsimulation_delta_time_sec = 1.0e-05\n\n# A multi-row value\nnumber_profile_layers = 2\nprofile_layers =\n     0  2880  1800  2200  100  100\n  1000  3600  2250  2300  100  100",
        bullets: [
          "Values may be integers, floating-point numbers, strings, paths, 0/1 booleans, yes/no switches, or multi-value records.",
          "Scientific notation is accepted for floating-point values.",
          "A line beginning with # is a comment and is ignored by Hercules.",
          "Whitespace around the name, equals sign, and values is flexible.",
          "Parameters may appear in any order, provided every required or conditionally required setting is present.",
        ],
      },
      {
        id: "required-conditional",
        title: "Required and conditional parameters",
        paragraphs: [
          "Many core settings are required. Others become required only when a related feature is enabled. For example, enabling an output or physical model also requires the paths and settings needed by that feature.",
          "For reproducible production cases, leave important choices explicit instead of relying on defaults that may change between revisions.",
        ],
        table: {
          headers: ["Status", "Meaning"],
          rows: [
            ["Required", "Every simulation must define the parameter."],
            ["Conditional", "Required when another feature or output is enabled."],
            ["Optional", "May be omitted; Hercules uses the documented default."],
          ],
        },
      },
      {
        id: "groups",
        title: "Parameter groups",
        table: {
          headers: ["Group", "What it controls"],
          rows: [
            ["Domain", "Origin, dimensions, depth range, azimuth, and coordinate conventions"],
            ["Time and resolution", "Start/end time, time step, maximum frequency, and nodes per wavelength"],
            ["Material", "Profiles, velocity models, density, damping, and constitutive settings"],
            ["Physics", "Elastic, acoustic, poroelastic, nonlinear, PML, and coupling choices"],
            ["Source", "Source directories, source representation, and incident waves"],
            ["Mesh", "Mesh generation/import mode, refinement inputs, and mesh output"],
            ["Output", "Monitor, station, plane, displacement, velocity, and checkpoint settings"],
            ["Optional features", "Topography, buildings, DRM, and site-specific model options"],
          ],
        },
      },
      {
        id: "minimal-structure",
        title: "Representative structure",
        paragraphs: [
          "A practical file is usually organized into commented blocks even though the parser does not require a particular order. The example below shows the shape of a file, not a complete simulation recipe.",
        ],
        codeLabel: "parameters.in",
        code: "# Domain\nregion_origin_latitude_deg  = 0.0\nregion_origin_longitude_deg = 0.0\nregion_length_east_m        = 1000.0\nregion_length_north_m       = 1000.0\nregion_depth_deep_m         = 1000.0\n\n# Time and mesh resolution\nsimulation_start_time_sec      = 0.0\nsimulation_end_time_sec        = 10.0\nsimulation_delta_time_sec      = 0.001\nsimulation_wave_max_freq_hz    = 2.0\nsimulation_node_per_wavelength = 8\n\n# Referenced inputs\nsource_directory = sourcefiles\ncvmdb_input_file = profile\n\n# Outputs\nmonitor_file = output/monitor.txt\nnumber_output_stations = 0\nnumber_output_planes   = 0",
      },
    ],
  },
  "input-files": {
    eyebrow: "Setting up simulations",
    title: "Input files and directories",
    summary:
      "parameters.in is the single configuration entry point. It references the source, material, topography, and other data required by the selected simulation features.",
    sections: [
      {
        id: "overview",
        title: "Required input structure",
        table: {
          headers: ["Input", "Contains"],
          rows: [
            ["parameters.in", "All numerical, physical, domain, source-path, and output configuration"],
            ["source files", "Point, force, rupture-plane, or incident-wave definitions and histories"],
            ["material models", "Layered, basin, etree, or other spatial velocity/density descriptions"],
            ["topography", "Surface samples or topographic configuration"],
          ],
        },
      },
      {
        id: "rules",
        title: "Path and case rules",
        bullets: [
          "Keep paths and units documented in parameters.in with nearby comments.",
          "Treat paths as relative to the run directory unless a site workflow requires absolute paths.",
          "Create output directories before launch unless the selected output path is documented to create them automatically.",
          "Archive the exact inputs with the monitor log and executable revision.",
        ],
        note: "Do not create separate numerical.in and physics.in files for a current psolve case. Their settings belong together in parameters.in.",
      },
      {
        id: "layout",
        title: "Typical case layout",
        paragraphs: [
          "Only the inputs enabled by parameters.in are required. Keeping each simulation self-contained makes it easier to move between systems and archive with a publication.",
        ],
        codeLabel: "directory tree",
        code: "case/\n├── parameters.in\n├── sourcefiles/\n│   └── source.in\n├── materialfiles/\n│   ├── x.in\n│   ├── y.in\n│   ├── z.in\n│   ├── vs.in\n│   ├── vp.in\n│   ├── rho.in\n│   └── index.in\n├── topography/\n│   └── topography.in\n└── output/",
      },
      {
        id: "sources",
        title: "Source input",
        paragraphs: [
          "The source file uses the same assignment and comment syntax as parameters.in. For a point source, it defines the source function, amplitude, position, and orientation. Hercules also supports force sources, rupture-plane descriptions, kinked planes, and incident plane waves.",
        ],
        codeLabel: "source.in — abbreviated point source",
        code: "type_of_source       = point\nsource_function_type = quadratic\naverage_risetime_sec = 0.5\nmoment_amplitude     = 1.0e15\n\nlonlat_or_cartesian  = 0\nhypocenter_lat_deg   = 38.419\nhypocenter_long_deg  = 27.129\nhypocenter_depth_m   = 500\nsource_strike_deg    = 0\nsource_dip_deg       = 90\nsource_rake_deg      = 0",
        note: "Validate source coordinates, units, orientation conventions, amplitude, and time history independently before scaling a simulation.",
      },
      {
        id: "velocity-model",
        title: "Three-dimensional velocity model",
        paragraphs: [
          "The structured 3D model is stored as seven files: x.in, y.in, z.in, vs.in, vp.in, rho.in, and index.in. x.in and y.in contain the unique horizontal coordinates. The z, velocity, and density arrays contain the vertical profiles; index.in identifies where each horizontal profile begins.",
        ],
        bullets: [
          "Write x and y coordinates in ascending order and without duplicates.",
          "For each horizontal point, write z coordinates in ascending order and without duplicates.",
          "Use metres for coordinates, metres per second for Vs and Vp, and the density units expected by the selected model workflow.",
          "The first value in each file records the number of entries that follow.",
        ],
      },
      {
        id: "topography",
        title: "Topography",
        paragraphs: [
          "When topography is enabled, topography_directory contains topography.in. Its first value is the regular spacing in the horizontal directions; subsequent values are surface elevations in metres.",
          "Samples begin at the lower-left domain corner. Complete the south-to-north sequence for one west-to-east position before moving to the next position.",
        ],
      },
    ],
  },
  "output-files": {
    eyebrow: "Reference",
    title: "Output files",
    summary:
      "Choose the smallest observation set that answers the scientific question, then select an output backend that scales with it.",
    sections: [
      {
        id: "outputs",
        title: "Available outputs",
        table: {
          headers: ["Output", "Use"],
          rows: [
            ["Stations", "Time histories at named or indexed observation points"],
            ["Planes", "Regularly sampled slices for wavefield visualization and analysis"],
            ["Mesh coordinates", "Spatial inspection of a full mesh or selected subdomain"],
            ["Checkpoints", "Restart from a saved solver state"],
            ["Monitor log", "Configuration, execution-path, timing, and diagnostic record"],
          ],
        },
      },
      {
        id: "stations",
        title: "Stations",
        paragraphs: [
          "GPU station paths interpolate supported elastic, acoustic, and poroelastic fields without making full nodal histories host-visible. Print rate controls temporal sampling; output quantities depend on the active physics and station configuration.",
        ],
      },
      {
        id: "planes",
        title: "Planes",
        paragraphs: [
          "Plane backends preserve the planedisplacements.N binary layout. Hybrid mode uses fused GPU interpolation and one packed result transfer; MPI-IO uses distributed views and configurable snapshot batches.",
        ],
      },
      {
        id: "checkpoints",
        title: "Checkpoints",
        paragraphs: [
          "Checkpointing intentionally activates host-visible compatibility paths for state that must be serialized. Measure the performance effect at the intended checkpoint rate instead of extrapolating from checkpoint-free runs.",
        ],
      },
    ],
  },
  validation: {
    eyebrow: "Reference",
    title: "Validation and benchmarks",
    summary:
      "Validate physical behavior, numerical convergence, and optimized execution paths independently before combining them in a production case.",
    sections: [
      {
        id: "strategy",
        title: "Validation strategy",
        bullets: [
          "Begin with homogeneous or layered analytical and semi-analytical references.",
          "Compare time series after documenting filtering, sampling, alignment, and component conventions.",
          "Run mesh and time-step convergence separately.",
          "Compare optimized GPU paths against a deterministic or generic reference when available.",
          "Record both waveform metrics and conservation or interface diagnostics.",
        ],
      },
      {
        id: "benchmarks",
        title: "Included benchmark families",
        table: {
          headers: ["Family", "Primary purpose"],
          rows: [
            ["CPS flat-water cases", "Acoustic–elastic coupling and attenuation comparisons"],
            ["Poroelastic u-p cases", "Wave speed, transfer functions, and saturated-layer response"],
            ["Topography tilt boundary", "Surface geometry and boundary handling"],
            ["Bowl basin", "Three-dimensional material structure and focusing"],
          ],
        },
      },
      {
        id: "performance",
        title: "Performance evidence",
        paragraphs: [
          "Report allocation shape, mesh size, time-step count, output rate, MPI and CUDA versions, active environment overrides, and per-stage timers. A speedup without this context is not reproducible evidence.",
        ],
      },
    ],
  },
  troubleshooting: {
    eyebrow: "Reference",
    title: "Troubleshooting",
    summary:
      "Start with the first explicit Hercules diagnostic, then separate launch, transport, physics, and I/O symptoms.",
    sections: [
      {
        id: "gpu-mapping",
        title: "More ranks than visible GPUs",
        paragraphs: [
          "The guard is intentional: each active local Phase B rank needs a GPU unless diagnostic oversubscription is explicitly enabled. Fix scheduler resource requests, rank placement, or CUDA_VISIBLE_DEVICES before weakening the guard.",
        ],
      },
      {
        id: "cuda-aware",
        title: "CUDA-aware MPI fails",
        bullets: [
          "Confirm the MPI library was built with CUDA support and that the runtime loads the same installation.",
          "Check UCX transport availability on compute nodes, not only login nodes.",
          "Disable cuda_ipc only when handle failures specifically implicate IPC.",
          "Use HERCULES_GPU_MPI_BACKEND=host as a diagnostic comparison, not an unexplained permanent fix.",
        ],
      },
      {
        id: "unstable",
        title: "Unstable or nonphysical wavefields",
        bullets: [
          "Verify material units, density, and wave-speed bounds.",
          "Check that the time step satisfies the discretization and fastest material speed.",
          "Inspect source orientation, scaling, and time function.",
          "Validate coupling and PML parameters in a smaller isolated model.",
        ],
      },
      {
        id: "output",
        title: "Missing or delayed plane output",
        paragraphs: [
          "MPI-IO plane snapshots become visible when their batch is flushed. Set HERCULES_PLANE_MPIIO_BATCH=1 when immediate per-snapshot visibility is required, accepting more collective I/O calls.",
        ],
      },
    ],
  },
};
