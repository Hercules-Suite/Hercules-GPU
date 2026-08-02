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
          "A run is defined by numerical and physics parameter files plus the model and source inputs they reference. Start from examples/simple or one of the benchmark_inputs cases so directory layout and units are explicit.",
        ],
        codeLabel: "typical case layout",
        code: "my-case/\n├── numerical.in\n├── physics.in\n├── sourcefiles/\n├── models/\n└── output/",
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
          "Read numerical, physics, source, model, and output configuration.",
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
  "input-files": {
    eyebrow: "Reference",
    title: "Input files",
    summary:
      "Keep numerical controls, physical assumptions, sources, and spatial models explicit and versioned with every case.",
    sections: [
      {
        id: "overview",
        title: "Input families",
        table: {
          headers: ["Input", "Contains"],
          rows: [
            ["numerical.in", "Domain dimensions, resolution, time stepping, partitioning, and output controls"],
            ["physics.in", "Material model, damping, analysis type, topography, DRM, and related physics"],
            ["source files", "Point, force, rupture-plane, or incident-wave definitions and histories"],
            ["material models", "Layered, basin, etree, or other spatial velocity/density descriptions"],
            ["topography", "Surface samples or topographic configuration"],
          ],
        },
      },
      {
        id: "rules",
        title: "Parameter-file rules",
        bullets: [
          "Use one assignment per line and preserve units in nearby comments.",
          "Treat paths as relative to the run directory unless a site workflow requires absolute paths.",
          "Keep optional settings explicit in production cases so future default changes are visible.",
          "Archive the exact inputs with the monitor log and executable revision.",
        ],
      },
      {
        id: "sources",
        title: "Sources",
        paragraphs: [
          "Hercules supports concentrated point and force sources, rupture-plane descriptions, kinked planes, and incident plane waves. Source coordinates, orientation, moment or force scaling, and time histories should be validated independently before a large simulation.",
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
