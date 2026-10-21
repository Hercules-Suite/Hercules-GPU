import { parameterReferencePages, parameterReferenceSummary } from "./parameter-reference";

export type DocTable = {
  headers: string[];
  rows: string[][];
};

export type DocSection = {
  id: string;
  title: string;
  paragraphs?: string[];
  equations?: string[];
  afterEquationsParagraphs?: string[];
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

type NonlinearExampleOptions = {
  model: string;
  columns: string;
  row: string;
  approximateGeostatic?: boolean;
  tensionCutoff?: boolean;
  errorTolerance?: string;
  backboneTolerance?: string;
  noSubsteps?: number;
  extraControls?: string[];
};

const classicMaterialColumns =
  "Vs_limit c_or_Su(Pa) phi(deg) strain_rate sensitivity hardening(Pa) dilatancy(deg) gamma0 psi m theta1 theta2 theta3 theta4 theta5 phi_RO";

function nonlinearModelExample({
  model,
  columns,
  row,
  approximateGeostatic = false,
  tensionCutoff = false,
  errorTolerance = "1.0e-8",
  backboneTolerance = "1.0e-8",
  noSubsteps = 100,
  extraControls = [],
}: NonlinearExampleOptions): string {
  return [
    "include_nonlinear_analysis = yes",
    "Istanbul_velocity_model = no",
    "geostatic_loading_time_sec = 0.0",
    "geostatic_cushion_time_sec = 0.0",
    `approximate_geostatic_state = ${approximateGeostatic ? "yes" : "no"}`,
    `material_model = ${model}`,
    "material_plasticity_type = rate_independent",
    `error_tolerance = ${errorTolerance}`,
    `backbone_errTol = ${backboneTolerance}`,
    `no_substeps = ${noSubsteps}`,
    "stiff_damp = 0.0",
    "freq_stiff_damp = 1.0",
    `tension_cutoff = ${tensionCutoff ? "yes" : "no"}`,
    ...extraControls,
    "material_properties_count = 2",
    "material_properties_list =",
    `# ${columns}`,
    ` 150 ${row}`,
    ` 450 ${row}`,
  ].join("\n");
}

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
          "The launcher is run_phaseA_phaseB.sh in the top level of the Hercules_GPU repository. Run it from the case directory so relative paths inside parameters.in resolve against that directory. BASE_IN is the path to the parameters.in file, not the case directory.",
        ],
        codeLabel: "shell",
        code: "cd /path/to/my-case\n\nexport PSOLVE_BIN=/path/to/Hercules_GPU/quake/forward/psolve\nexport BASE_IN=$PWD/parameters.in\nexport PHASE_A_NP=64\nexport PHASE_B_NP=8\nexport PHASE_B_OMP_NUM_THREADS=8\n\n/path/to/Hercules_GPU/run_phaseA_phaseB.sh",
        note: "PSOLVE_BIN may be omitted only when the executable is already at Hercules_GPU/quake/forward/psolve or Hercules_GPU/psolve. Setting the absolute path explicitly is recommended in job scripts.",
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
            ["GNU Scientific Library (GSL)", "Polynomial, root-finding, vector, and nonlinear solver routines", "Both GSL headers and the gsl/gslcblas libraries are required"],
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
        code: "CC       = mpicc\nCXX      = mpicxx\nLD       = mpicxx\nNVCC     = nvcc\n\nGSL_DIR  = /path/to/gsl\nPROJ_DIR = /path/to/proj\nCUDA_DIR = /path/to/cuda\n\nCFLAGS   += -O3 -I$(GSL_DIR)/include\nCXXFLAGS += -O3 -I$(GSL_DIR)/include\nCPPFLAGS += -I$(GSL_DIR)/include -I$(PROJ_DIR)/include -DPROJ\nNVCC     += -I$(GSL_DIR)/include\nLDFLAGS  += -L$(GSL_DIR)/lib -lgsl -lgslcblas\nLDFLAGS  += -L$(PROJ_DIR)/lib -lproj",
        note: "Exact variables differ across existing site configurations. Begin with the closest checked-in example and preserve your center's compiler-wrapper conventions.",
      },
      {
        id: "gsl",
        title: "Make GSL visible to every compiler stage",
        paragraphs: [
          "Hercules compiles source that includes gsl_poly.h, gsl_vector.h, gsl_multiroots.h, gsl_roots.h, and related headers. The final solver is linked with libgsl and libgslcblas.",
          "On a module-based HPC system, load the GSL module before compiling and running. Otherwise, add its include and library directories to user.mk. CUDA compilation also needs the GSL include directory because solver headers are visible to nvcc.",
        ],
        codeLabel: "module-based environment",
        code: "module load gsl\n\n# Verify that the active environment exposes GSL\nprintenv TACC_GSL_INC\nprintenv TACC_GSL_LIB",
        note: "TACC_GSL_INC and TACC_GSL_LIB are site-specific variables used by existing TACC configurations. Other systems may expose different module variables; use the actual include and library paths supplied by that system.",
      },
      {
        id: "build",
        title: "Build and inspect",
        paragraphs: [
          "Make embeds the full Git commit ID and tracked-tree state (clean, modified, or unavailable) in a small host-only object. The startup Some Input Data block prints them even when the executable is copied elsewhere; it never reads Git at runtime. A source copy without .git or a build without working Git reports unknown. Untracked/ignored files are excluded. Incremental builds refresh changed metadata without recompiling CUDA kernels; preserve local patches, compiler options and libraries separately because a commit ID alone does not identify them.",
          "Use an explicit parallel make job count appropriate for the login or build node. The large host-side solver orchestration in psolve.cc is compiled by the MPI C++ compiler while CUDA kernels and launch wrappers in psolve_gpu_launch.cu are compiled independently by nvcc; make -j 8 can run those jobs concurrently. Incremental builds also avoid rerunning nvcc when only host orchestration changes.",
          "PTX assembly is the expensive phase for nonlinear constitutive CUDA. Sources under quake/forward/nonlinear instantiate separate production kernels in nonlinear_gpu_classic.o, nonlinear_gpu_manzari_dafalias.o, nonlinear_gpu_pdmy03.o, nonlinear_gpu_riva_sand.o, and nonlinear_gpu_riva_sand_intermediate_bias_research.o, plus small runtime and material-point objects. The research successor has a dedicated material-point CUDA object for its 139-value state. The runtime launches only the family selected by material_model. Parallel make can therefore run independent NVCC/PTX-assembly jobs concurrently without relocatable device code, noinline annotations, or equation changes.",
        ],
        codeLabel: "shell",
        code: "make clean\nmake -j 8\n\n# Expected executables\nls octor/octor quake/forward/psolve",
        bullets: [
          "Use a numeric -j limit. Bare make -j permits unbounded outer build concurrency and can oversubscribe a shared build node.",
          "Treat missing gsl headers as an include-path problem and missing gsl_* symbols as a GSL link-path or link-order problem.",
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
          "Advance elastic, acoustic, or poroelastic fields; eligible nonlinear solids additionally update stress, plastic strain, and constitutive internal variables.",
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
      {
        id: "octree-domain-dimensions",
        title: "Octree-compatible domain dimensions",
        paragraphs: [
          String.raw`Hercules cannot use three arbitrary domain lengths. The physical north, east, and depth extents—and their maximum, \(L_{\max}\)—are defined from the exact parameters below. At octree level \(N\), the common dyadic increment is \(\Delta_N=L_{\max}/2^N\); every domain extent must be an integer multiple of that increment.`,
          String.raw`Equivalently, an extent \(D\) uses an integer multiplier \(m\), or a binary sum with coefficients \(a_i\in\{0,1\}\). This lets a shorter direction occupy a rectangular subset of the power-of-two root cube without breaking the octree grid.`,
        ],
        equations: [
          String.raw`\begin{aligned}
L_{\mathrm{north}} &= \mathtt{region\_length\_north\_m}, \\
L_{\mathrm{east}}  &= \mathtt{region\_length\_east\_m}, \\
L_{\mathrm{depth}} &= \mathtt{region\_depth\_deep\_m}-\mathtt{region\_depth\_shallow\_m}, \\
L_{\max} &= \max\!\left(L_{\mathrm{north}},L_{\mathrm{east}},L_{\mathrm{depth}}\right),
\end{aligned}`,
          String.raw`\begin{aligned}
D &= m\,\frac{L_{\max}}{2^N}, &&m\in\mathbb{Z},\quad 0<D\le L_{\max}, \\
D &= \sum_{i=0}^{N} a_i\,\frac{L_{\max}}{2^i}, &&a_i\in\{0,1\}.
\end{aligned}`,
        ],
        afterEquationsParagraphs: [
          String.raw`For \(L_{\max}=\mathtt{region\_length\_east\_m}=1{,}000{,}000\,\mathrm{m}\), the successive dyadic sizes are \(1{,}000{,}000\), \(500{,}000\), \(250{,}000\), \(125{,}000\), \(62{,}500\), and \(31{,}250\,\mathrm{m}\). The choice \(\mathtt{region\_length\_north\_m}=750{,}000\,\mathrm{m}\) is valid because \(750{,}000=500{,}000+250{,}000\). A depth thickness of \(187{,}500\,\mathrm{m}\) is valid because \(187{,}500=125{,}000+62{,}500\). If \(\mathtt{region\_depth\_shallow\_m}=10{,}000\,\mathrm{m}\), set \(\mathtt{region\_depth\_deep\_m}=197{,}500\,\mathrm{m}\) to obtain that thickness.`,
        ],
        codeLabel: "parameters.in — octree-compatible 1,000 km × 750 km × 187.5 km domain",
        code: "region_depth_shallow_m = 10000\nregion_depth_deep_m    = 197500\nregion_length_north_m  = 750000\nregion_length_east_m   = 1000000",
        note: "Use exact whole-metre extents with the current Octor path: it reduces the physical dimensions using integer metre values. The dyadic rule applies to the depth thickness, not to region_depth_deep_m by itself unless region_depth_shallow_m is zero.",
      },
    ],
  },
  running: {
    eyebrow: "Run simulations",
    title: "Running Hercules at scale",
    summary:
      "Use the repository launcher to connect one parameters.in file to mesh construction and GPU solving, with explicit executable and resource settings.",
    sections: [
      {
        id: "launcher",
        title: "Launcher location and working directory",
        paragraphs: [
          "run_phaseA_phaseB.sh is located in the top level of the Hercules_GPU repository, next to the quake and octor directories. It is executable and may be called by absolute path from any case directory.",
          "The launcher records the current directory as the run directory. Start it from the directory containing parameters.in when that file uses relative paths such as sourcefiles, models, or output.",
        ],
        codeLabel: "repository layout",
        code: "Hercules_GPU/\n├── run_phaseA_phaseB.sh\n├── quake/\n│   └── forward/\n│       └── psolve\n├── octor/\n└── ...",
      },
      {
        id: "launcher-inputs",
        title: "Launcher inputs",
        table: {
          headers: ["Setting", "Required", "Meaning"],
          rows: [
            ["PSOLVE_BIN", "When not auto-detected", "Absolute path to the executable psolve binary"],
            ["BASE_IN", "Unified mode", "Path to the complete parameters.in file; this is a file, not a directory"],
            ["PHASE_A_NP", "Unified mode", "Total MPI ranks launched for Phase A mesh work"],
            ["PHASE_B_NP", "Recommended", "Ranks that become active GPU solver ranks, normally nodes × GPUs per node"],
            ["PHASE_B_OMP_NUM_THREADS", "Unified mode when not inferred", "OpenMP threads assigned to each active Phase B rank"],
            ["MPI_RUN", "No", "MPI launcher command; defaults to mpirun"],
          ],
        },
        note: "PSOLVE_BIN is auto-detected only at Hercules_GPU/quake/forward/psolve or Hercules_GPU/psolve. A fresh source checkout has no psolve executable until Hercules is built.",
      },
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
        code: "cd /scratch/project/case\n\nexport HERCULES_PHASE_MODE=unified\nexport PSOLVE_BIN=/path/to/Hercules_GPU/quake/forward/psolve\nexport BASE_IN=$PWD/parameters.in\nexport PHASE_A_NP=2304\nexport PHASE_B_NP=256\nexport PHASE_B_OMP_NUM_THREADS=8\n\n/path/to/Hercules_GPU/run_phaseA_phaseB.sh",
        paragraphs: [
          "Unified is the default phase mode. The launcher reads BASE_IN, creates temporary Phase A, Phase B, and unified variants of the same parameter file, and adjusts the mesh-mode settings required by each stage.",
          "The MPI job is launched with PHASE_A_NP total ranks. PHASE_B_NP of those ranks become active solver ranks after mesh construction. Active GPU ranks divide each node's allowed CPU affinity mask into disjoint blocks while inactive ranks wait for Phase B to finish.",
        ],
      },
      {
        id: "positional",
        title: "Equivalent positional form",
        paragraphs: [
          "The first three launcher arguments override PHASE_A_NP, PHASE_B_NP, and BASE_IN respectively. Environment variables remain useful for PSOLVE_BIN and thread placement.",
        ],
        codeLabel: "shell",
        code: "cd /scratch/project/case\nexport PSOLVE_BIN=/path/to/Hercules_GPU/quake/forward/psolve\nexport PHASE_B_OMP_NUM_THREADS=8\n\n/path/to/Hercules_GPU/run_phaseA_phaseB.sh 2304 256 $PWD/parameters.in",
      },
      {
        id: "two-phase",
        title: "Two-phase launch",
        paragraphs: [
          "Two-phase mode is useful when a mesh will be reused. With BASE_IN set, the launcher derives separate mesh-only and solver-only parameter files from the same parameters.in. The sharded backend lets Phase A ranks write concurrently; the MPI-IO backend creates one indexed archive; single retains the serialized compatibility path.",
        ],
        codeLabel: "shell",
        code: "cd /scratch/project/case\n\nexport HERCULES_PHASE_MODE=two-phase\nexport HERCULES_MESH_IO_MODE=sharded\nexport PSOLVE_BIN=/path/to/Hercules_GPU/quake/forward/psolve\nexport BASE_IN=$PWD/parameters.in\nexport PHASE_A_NP=2304\nexport PHASE_B_NP=256\nexport PHASE_B_OMP_NUM_THREADS=8\n\n/path/to/Hercules_GPU/run_phaseA_phaseB.sh",
      },
      {
        id: "preflight",
        title: "Check the launcher summary",
        bullets: [
          "Confirm the printed psolve binary is the executable built for the target system.",
          "Confirm Base input points to the intended parameters.in file.",
          "Confirm the current run directory matches the base directory for relative input and output paths.",
          "Confirm Phase A ranks, Phase B ranks, OpenMP threads, phase mode, and mesh I/O mode match the allocation.",
          "On multi-GPU nodes, use a PHASE_B_NP value divisible by the allocated node count.",
        ],
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
        id: "scalar-schedules",
        title: "Scalar communication schedules and profiling",
        paragraphs: [
          "HERCULES_SCALAR_MPI_SCHEDULE=compact (default) builds separate acoustic, pore-pressure and open-pore-target schedules once at setup, on both host and CUDA-aware MPI backends. It retains active field nodes, ownership, hanging-node parents, interface traces and station/plane/receiver ghosts. Neighbor node order is checked before filtering. Ranks without local active elements still communicate required owned or shared nodes. Startup reports original and retained node entries and peer schedules; entirely empty device schedules skip synchronization. Initialization and solid-vector exchanges keep their original maps, and host sparse transfer lists remain conservative supersets.",
          "Set HERCULES_SCALAR_MPI_SCHEDULE=full to restore the original general scalar node mappings for comparison. Set batching and overlap to 0 as well to restore the previous sequential exchange path. The change preserves time levels, quadrature, field constraints and physical coupling; diagnostic pore pressure remains available at required output states. Target-system CUDA/MPI propagation and performance comparisons are required before claiming a speedup. Test standalone and coupled elastic/acoustic/poroelastic cases, sealed/open-pore and matching/nonmatching interfaces, both acoustic scalar formulations, applicable viscoelastic/nonlinear paths, outputs and restart. Existing poroelastic CPML exclusions and VT geometric constraint limitations remain.",
          "HERCULES_SCALAR_MPI_BATCH=1 (default) combines sealed acoustic/pore RHS and field exchanges when both fields use the same backend, reducing the total dependency sequence from 12 to 8 phases including solid exchanges. Open-pore coupling must finish the new acoustic field before computing its hydraulic target; only pore RHS and target contributions combine, reducing 14 to 12 phases. Each peer packet concatenates the existing field mappings without padding. Host and CUDA-aware MPI both support batching. These phase reductions are not measured speedups. Set 0 to keep fields separate.",
          "HERCULES_GPU_MPI_OVERLAP=1 (default) starts anchored scalar RHS messages, queues updates for nodes unaffected by incoming contributions, then completes communication and updates boundary nodes. The pore partition includes incoming hydraulic targets. Both partitions use the existing formulas. Scalar and solid-vector device unpack no longer wait immediately: CUDA events protect receive-buffer reuse and MPI completion protects send-buffer reuse. GPU packing still completes before MPI sends. There are no per-step allocations or extra collectives, and the final GPU tail is drained before solver timing ends. Host staging remains synchronous. This does not split all solid element kernels. Set overlap and batching to 0 for the previous sequential exchange path.",
          "HERCULES_PORE_PRESSURE_UPDATE=output (default) applies only to fresh linear undrained fast-path runs. Zero permeability and initial zero displacement/excess pressure let the diagnostic recurrence telescope to p = -M_lumped^-1 Q^T u, using the existing quadrature, mass and hanging constraints. All ranks reconstruct before requested ASCII pore-pressure station output and for the final two states; skipped steps perform no pore RHS/update/exchange. Only actual poroelastic station consumers trigger reconstruction. Displacement-only binary stations and plane/DRM motion do not need pore pressure. Sealed acoustic coupling uses acoustic-only packets between reconstructions and combined packets when pressure is required. Both MPI backends and full/compact maps are supported. Set step for the original incremental recurrence. Checkpoint reads/writes, nonlinear analysis, diffusion, drainage and open-pore coupling retain every-step pressure updates. Ru output requires nonlinear analysis and an initial effective-stress reference, so it always retains every-step updates. The checkpoint schema is unchanged. Startup and final counters report the effective path; rounding can differ from incremental accumulation.",
          "HERCULES_GPU_MPI_PROFILE=1 enables device-helper and scalar-batch calls, actual packet byte/message counts and CPU elapsed times; default 0 disables it. Combined packets have separate acoustic+pore RHS/field and pore RHS+target categories. End-of-loop output gives cumulative rank average/min/max. With overlap, unpack time measures enqueue work and receive-wait time can overlap GPU updates; these are not isolated kernel or network timings. Profiling adds no GPU synchronization or per-step collective. All communication settings must agree across solver ranks. Existing host schedule and dangling counters cover only the host helper; batched host traffic uses the scalar-batch profile. Scalar exchanges remain nested under Compute Physics rather than top-level Communication. Pipeline and interior/boundary enqueue timers are separate; do not sum nested timers.",
        ],
        codeLabel: "shell — repeat with full and with host on the same mesh and partition",
        code: "export HERCULES_SCALAR_MPI_SCHEDULE=compact\nexport HERCULES_GPU_MPI_PROFILE=1\nexport HERCULES_GPU_MPI_BACKEND=cuda_aware",
      },
      {
        id: "resident-fields",
        title: "Resident-field strategy",
        bullets: [
          "Elastic displacement history remains on the GPU; only boundary and adjustment nodes move when required.",
          "Nonlinear-solid stress, plastic strain, hardening/backbone state, and internal-force assembly remain on the GPU by default. Set HERCULES_NONLINEAR_BACKEND=cpu explicitly to select the CPU/GSL reference backend.",
          "Nonlinear station output gathers compact quadrature-point records only for local station elements; it does not download the full nonlinear state at each output step.",
          "Acoustic pressure histories swap device pointers and patch only schedule or adjustment nodes.",
          "Poroelastic pore pressure uses the same sparse/resident pattern.",
          "Station interpolation runs on the GPU for supported fields and returns compact observations.",
        ],
      },
      {
        id: "nonlinear-backend",
        title: "Nonlinear solid backend",
        paragraphs: [
          "The CUDA nonlinear constitutive and force path is the default when HERCULES_NONLINEAR_BACKEND is unset and keeps constitutive calculations in double precision. Set HERCULES_NONLINEAR_BACKEND=cpu explicitly for the established CPU/GSL reference models. ManzariDafalias, PressureDependMultiYield03, RIVASand, RIVASandCustom, RIVASand02, and RIVASandCustom02 are device-first and currently require the GPU backend for production-sized runs; their material-point regressions execute the same allocation-free constitutive cores on host and device. The selected backend is reported to standard output at startup.",
          "For a ManzariDafalias performance diagnosis, set HERCULES_NONLINEAR_STATS_INTERVAL=N to sample every N global time steps. The global report classifies integration-point slots as zero/inactive, one-step, or multi-step and gives total, mean, and maximum accepted constitutive substeps. With Forward Euler, one-step covers an elastic-predictor return or a one-part plastic update. Under the default fixed policy, a plastic point is multi-step whenever no_substeps is greater than one; under the opt-in demand policy, the classification reflects the requested count. The diagnostic adds a device scan and global MPI reductions on sampled steps, so leave it unset or zero for production timing and use the same value on every rank.",
          "The specialized ManzariDafalias production kernel omits the four legacy previous-state shadow arrays used by the classic model kernel. This removes 1,216 bytes of unused device allocation per nonlinear element and the corresponding writes on each updated time step. ManzariDafalias keeps its larger failure-only diagnostic record in per-thread shared memory. PDMY03 retains the legacy shadow allocation and writes together with its compact thread-local diagnostic record because controlled target-system production measurements favored that generated kernel, despite lower static PTXAS spill totals when the shadows were omitted. The PDMY03 shadows do not enter its constitutive equations, so this model-specific execution choice does not change committed stress or history.",
          "Fixed ManzariDafalias, PDMY03, and RIVA-Sand material inputs live in model-specific host and device arrays rather than the common constants of every nonlinear element. MD and PDMY03 allocate one interpolated material block per element only while that family is active. All four RIVA-Sand names keep one interpolated material block per element, including the setup-time stiffness cache derived from G0. Startup output reports active material storage separately from evolving history; unrelated nonlinear families pay none of this model-input memory.",
          "When abc_type=cpml, solid elements in the absorbing band remain linear elastic and are excluded from both CPU and GPU nonlinear constitutive maps regardless of their Vs. Nonlinearity is confined to the physical interior because the CPML correction assumes a linear time-invariant native operator.",
          "When include_acoustic=yes, acoustic elements are excluded from both CPU and GPU solid constitutive maps regardless of their Vs or legacy nl_id payload. They remain governed only by the selected acoustic formulation and never acquire plastic-strain or nonlinear solid history.",
          "Treat the GPU path as a separately validated execution backend. For models supported by both production backends, compare nonlinear station histories, final yield statistics, and force or displacement norms against a CPU run before using a new calibration. ManzariDafalias, PressureDependMultiYield03, RIVASand, RIVASandCustom, RIVASand02, and RIVASandCustom02 successors are GPU-only at element scale, so use their shared host/CUDA material-point regressions and an appropriate independent characterization before a GPU run.",
          "For production timing, verify both the total element count and the nonlinear count after regenerating any mesh whose velocity profile changed. Reusing an older mesh also reuses its stored Vs payload and can unintentionally select nearly the entire domain.",
        ],
        codeLabel: "shell",
        code: "export HERCULES_NONLINEAR_BACKEND=gpu\nexport HERCULES_GPU_MPI_BACKEND=cuda_aware\nexport HERCULES_DISP_SYNC_MODE=boundary",
        note: "The GPU backend changes execution only. Nonlinear geostatic, tolerance, damping, substep, and material_properties_list values remain case-specific physical/numerical inputs.",
      },
      {
        id: "nonlinear-material-point",
        title: "Material-point CPU/GPU regression",
        paragraphs: [
          "Before an element or production comparison, run a prescribed tensor-strain history through the CPU and CUDA constitutive routines with the standalone material-point mode. It writes complete CPU and GPU state histories and returns nonzero on a status, finite-value, or tolerance mismatch.",
          "The twenty-one-fixture suite covers EP, the FA/FAM quartic correctors, BAE and BAH, GQH, MKZ, Ramberg--Osgood, five ManzariDafalias paths including fixed and demand Forward-Euler policies, cyclic plus zero/tensile-onset PressureDependMultiYield03 paths, frozen, low-static-bias, ramped-geostatic-admission, plus custom-reference RIVA-Sand effective-stress/BrickUP replays, and the distinct intermediate-bias research successor. Existing families compare the CPU GSL-backed polynomial or scalar root solvers with CUDA; the device-first fixtures execute the same allocation-free implementation on host and device and record their internal-state and substep diagnostics.",
          "BAE, BAH, GQH, MKZ, and Ramberg--Osgood now enter their incremental bounding/backbone update directly. Those models decide loading and reversal from the strain increment and stored stress/history; the generic total-strain elastic predictor remains limited to conventional yield-surface models. A repeated, true zero strain increment preserves the converged stress and history.",
          "For these pressure-insensitive von Mises backbones, a purely volumetric increment updates only hydrostatic stress and preserves cyclic history instead of entering the scalar kappa solver. At untouched virgin state (kappa=100), deviatoric trial stress below 1.0e-3 Pa follows the incremental elastic asymptote until the backbone root is numerically resolvable. The onset rule is disabled after cyclic history is established, so it does not replace the model's unloading or reversal logic.",
          "This dispatch boundary follows the OpenSees MultiaxialCyclicPlasticity reference design: the cyclic integrator owns its zero-increment and loading/reversal decisions. The Hercules material-point suite checks that invariant independently on both CPU and CUDA backends.",
          "For a cross-code characterization, compare_opensees_bae.py maps the BAE constants into an OpenSees one-brick affine simple-shear test, converts tensorial epsilon_xy to engineering gamma_xy, and reports component stress errors plus integrated shear-work differences. The closest algorithmic mapping uses Ho=0, h=psi0 times G, and beta=1; beta=0.5 is available as an integration-sensitivity run.",
          "The standalone opensees_bae_synthetic.tcl fixture embeds the mapped constants and complete 20-step history. OpenSees can therefore run on a different system, with its CSV passed back through compare_opensees_bae.py using --opensees-csv. The 17-step legacy BAE suite fixture is not interchangeable with this dedicated cross-code history.",
          "The OpenSees driver uses a Path time series and prescribed top-face displacement, verifies the actual control-node displacement at every recorded step, and marks its CSV as driver version 2. This preserves displacement through the repeated zero point; using zero-increment load control there would instead release displacement control and re-equilibrate the brick. The comparison tool rejects older CSVs, accepts and reports penalty-enforcement shear-path error up to 1e-12, and retains the actual OpenSees strain in its detailed output.",
          "The OpenSees comparison is intentionally not described as an equivalence test. OpenSees uses coupled Newton iteration and strain splitting and resets trial mapping variables for a committed near-zero deviatoric increment, while Hercules uses a scalar bracketed root and preserves its complete state for a true zero increment. The comparison history includes a dwell and reload to expose that difference.",
          "Two canonical ManzariDafalias/OpenSees workflows extend the short regression: drained monotonic triaxial compression at constant lateral confinement, and constant-volume cyclic simple shear. The drained comparator reports radial-confinement drift and q/p so strength and contraction/dilation tendency are visible. The cyclic comparator reports effective-stress reduction, inferred excess pore pressure delta_u=p'_0-p', r_u, and stress reversals.",
          "Both canonical OpenSees drivers consolidate in a forced-elastic stage before enabling plasticity. The conversion tool seeds Hercules with the recorded stress, strain, alpha, alpha-in, and fabric at that boundary, then replays only the common elastoplastic history. This prevents the elastic consolidation strain from being reinterpreted through Hercules' pressure-dependent plastic stage and gives both codes the same initial effective stress.",
          "The drained driver uses one retained corner for all affine face constraints and records SSPbrick's actual center-point strain, stress, and material state. Obsolete version-5 output used a non-affine multiple-retained-node constraint layout, version 6 used a response prefix unsupported by SSPbrick, and version 7 selected an accumulated outer tangent that became ill-conditioned during early drained plastic loading. Version 8 restored OpenSees's default elastic outer tangent but could accept a small displacement correction while retaining about 22 kPa of lateral stress error. The current driver version 9 uses a force-residual convergence test and independently requires each recorded radial material stress to match the confinement target within 1 Pa.",
          "The cyclic workflow is a constitutive undrained characterization, not a coupled U-P validation. Its inferred pressure does not exercise Hercules pressure degrees of freedom, storage, permeability and Darcy flow, drainage boundaries, or solid-fluid coupling; those require a separate BrickUP benchmark.",
          "The PDMY03 fixture starts from 100 kPa isotropic effective confinement and reverses a tensor-shear path through its automatic nested surfaces. It compares host and CUDA stress, incremental plastic strain, active surface, translated surface-center norm, pressureD, surface crossings, requested subincrements, and clamping. It is a numerical regression, not an OpenSees equivalence test or a calibrated liquefaction validation.",
          "A companion OpenSees 3.8 affine-brick driver establishes the same 100 kPa confinement in its elastic stage, switches PDMY03 to plastic stage 1, and applies the same shear history. compare_opensees_pdmy03.py checks all six stress components and the tensor-shear path; the reference development run obtained stress NRMSE about 1.9e-10. plot_opensees_pdmy03.py writes PNG and PDF panels for shear response, effective mean stress, strain history, and cross-code differences. This characterizes the standard automatic-surface translation and PPZ arithmetic for the synthetic inputs.",
          "Two PDMY03 SSPbrickUP/BrickUP workflows consolidate a sealed unit element and enable plastic stage 1. One applies a reversing undrained triaxial path at constant total radial confinement; the other preserves the consolidated lateral strains, maintains total vertical confinement, and applies four cycles of x-z simple shear, driving the OpenSees reference r_u approximately to one. The generated Hercules replays invoke the production CUDA pressure-force and pressure-update kernels together with PDMY03. Their shared comparator reports effective- and total-stress NRMSE, pore-pressure and r_u differences, and confinement drift; the coupled plotter writes stress-response, pressure, effective-stress-path, and confinement-drift panels for either test. Driver version 2 explicitly supplies OpenSees 3.8.0's mType argument between the phase-transformation angle and ca; version-1 reference CSVs are rejected because the omitted argument shifted all later PDMY03 constants. The generator converts only stress-like values from the OpenSees kPa fixture to Hercules Pa; liquefaction1 and liquefaction2 remain dimensionless. Hercules also matches OpenSees's zero-initialized maximum-pressure history, which is first advanced by a committed plastic-stage update and controls PPZ damage after reversal. Current material-point CSV output labels all PDMY03 surface/substep diagnostics; the comparator repairs existing histories written with the earlier header that omitted those six labels. These are uniform zero-permeability coupled benchmarks; reviewed Horizon results are still required before fixing acceptance thresholds.",
          "The coupled SSPbrickUP/BrickUP benchmark fills the uniform sealed part of that gap. OpenSees runs a zero-permeability SSPbrickUP -lumped element under undrained triaxial compression; Hercules replays the exact affine strain history and invokes its production BrickUP pressure-force and pressure-update CUDA kernels on a unit element. Driver version 4 uses conventional symmetry supports, distributed face tractions, and identical prescribed axial displacement on all four top nodes. Its coupled stage follows the official OpenSees Penalty/Newmark pattern so constrained solid velocity drives the U-P pressure equation, while avoiding a prescribed equalDOF retained node. An independent CPU reference checks p-p0=-M(ev-ev0), with M=Kf/n. OpenSees's row-summed pressure capacity and Hercules's lumped nodal storage use the same representation. The validated fixture bounds are effective-stress NRMSE at most 5.0e-5, pore-pressure difference at most 1.0e-6 kPa, and Hercules total radial-confinement error at most 2.0e-1 kPa; a thresholded comparison prints PASS. The companion plot shows deviatoric stress, pressure buildup, effective stress path, and confinement drift. These bounds characterize the 100-kPa fixture rather than production tolerances. This deliberately isolates the common uniform mode and does not claim validation of Darcy diffusion, drainage/open-pore boundaries, nonuniform pressure, multi-element exchange, or the complete production time integrator.",
          "RIVA-Sand uses an explicit poromechanical split in its BrickUP fixture. RIVA-Sand evolves effective skeleton stress from total strain, including its elastic/plastic volumetric compatibility. BrickUP independently enforces fluid mass balance, d(zeta)=alpha*d(eps_v)+d(p_f)/M, and total stress is assembled once as sigma_total=sigma_effective-alpha*p_f*I. Hercules never injects RIVA-Sand's standalone effective-pressure-loss diagnostic into p_f. The fixture compares all 93 frozen flattened RIVA-Sand numerical history values and the separate committed geostatic-admission selector between host and CUDA, executes the production uniform BrickUP pressure kernels, reports the sealed fluid-balance residual, and writes effective stress, pressure, and combined total stress. Its zero-increment check also requires every RIVA-Sand history value and selector to remain unchanged. This verifies the implementation split and uniform sealed arithmetic, but a full stress-controlled coupled element and restart/reproduction run remain required before declaring the RIVA-Sand/BrickUP combination scientifically validated for production liquefaction.",
          "The intermediate-bias research successor has a separate material-point CUDA object so its larger kernel does not increase the established all-model regression translation unit. Its fixture starts from admitted K0 effective stress with nonzero static shear bias, compares all 138 flattened research state values and lifecycle selectors between the shared host/device core, checks exact zero-increment preservation, and independently invokes the uniform sealed BrickUP pressure replay. This is numerical implementation coverage for the research identifier; it is not production RIVA-Sand equivalence or physical validation.",
          "The fixture constants and strain path are test inputs, not a substitute for calibrated production properties. Off-diagonal inputs are tensorial shear strains, so epsilon_xy equals one half of engineering gamma_xy. A single point validates correctness only; use a many-element case for performance measurement.",
        ],
        codeLabel: "shell",
        code: "cd benchmark_inputs/nonlinear_material_point\n./run_suite.sh ../../quake/forward/psolve\n\nOpenSees opensees_manzari_dafalias_drained_triaxial.tcl opensees_md_drained.csv\npython3 generate_hercules_manzari_dafalias_path.py --opensees-csv opensees_md_drained.csv --output manzari_md_drained.in\nmpirun -np 1 ../../quake/forward/psolve --nonlinear-material-point manzari_md_drained.in\npython3 compare_opensees_manzari_dafalias_canonical.py --mode drained --hercules-csv hercules_md_drained_triaxial.cpu.csv --opensees-csv opensees_md_drained.csv\n\nOpenSees opensees_manzari_dafalias_sspbrickup_undrained_triaxial.tcl opensees_md_sspbrickup.csv\npython3 generate_hercules_manzari_dafalias_path.py --opensees-csv opensees_md_sspbrickup.csv --output manzari_md_brickup.in --output-prefix manzari_md_brickup\nmpirun -np 1 ../../quake/forward/psolve --nonlinear-material-point manzari_md_brickup.in\npython3 compare_opensees_manzari_dafalias_brickup.py --hercules-csv manzari_md_brickup.cpu.csv --opensees-csv opensees_md_sspbrickup.csv --output-dir opensees_md_sspbrickup_comparison --max-effective-stress-nrmse 5.0e-5 --max-pore-pressure-difference-kpa 1.0e-6 --max-total-radial-error-kpa 2.0e-1\npython3 plot_opensees_manzari_dafalias_brickup.py --comparison-csv opensees_md_sspbrickup_comparison/comparison.csv --output-prefix opensees_md_sspbrickup_comparison/hercules_vs_opensees\n\nOpenSees opensees_pdmy03_sspbrickup_undrained_triaxial.tcl opensees_pdmy03_sspbrickup.csv\npython3 generate_hercules_pdmy03_brickup_path.py --opensees-csv opensees_pdmy03_sspbrickup.csv --output pdmy03_brickup.in --output-prefix pdmy03_brickup\nmpirun -np 1 ../../quake/forward/psolve --nonlinear-material-point pdmy03_brickup.in\npython3 compare_opensees_pdmy03_brickup.py --hercules-csv pdmy03_brickup.cpu.csv --opensees-csv opensees_pdmy03_sspbrickup.csv --output-dir opensees_pdmy03_sspbrickup_comparison\n\nOpenSees opensees_pdmy03_sspbrickup_undrained_cyclic_shear.tcl opensees_pdmy03_sspbrickup_cyclic.csv\npython3 generate_hercules_pdmy03_brickup_path.py --opensees-csv opensees_pdmy03_sspbrickup_cyclic.csv --output pdmy03_brickup_cyclic.in --output-prefix pdmy03_brickup_cyclic\nmpirun -np 1 ../../quake/forward/psolve --nonlinear-material-point pdmy03_brickup_cyclic.in\npython3 compare_opensees_pdmy03_brickup.py --hercules-csv pdmy03_brickup_cyclic.cpu.csv --opensees-csv opensees_pdmy03_sspbrickup_cyclic.csv --output-dir opensees_pdmy03_sspbrickup_cyclic_comparison",
        note: "The supplied von Mises fixtures are synthetic numerical regression tests. Replace them with constants and histories derived from the calibrated case for scientific validation.",
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
          "Plane output defaults to mpiio when HERCULES_PLANE_IO_MODE is unset or empty. Fused GPU interpolation and a packed device-to-host transfer feed distributed MPI-IO writes. Set hybrid for GPU interpolation followed by rank-0 gather/fwrite, or legacy for CPU interpolation followed by rank-0 writing. All three preserve the planedisplacements.N binary layout.",
          "Dedicated plane I/O ranks (IO_PES>0) select the legacy I/O-pool backend and log the fallback when a GPU backend was requested or selected by default. Use IO_PES=0 for MPI-IO or hybrid output.",
          "MPI-IO batches 64 packed snapshots by default so collective calls are amortized; a partial batch is flushed before timing reduction and close. HERCULES_PLANE_MPIIO_BATCH sets a positive snapshot count. Files expose snapshots when the batch is flushed; use 1 for per-snapshot visibility. This setting has no effect on hybrid or legacy output.",
        ],
        codeLabel: "shell",
        code: "export IO_PES=0\n# Optional: make the solver defaults explicit\nexport HERCULES_PLANE_IO_MODE=mpiio\nexport HERCULES_PLANE_MPIIO_BATCH=64",
        note: "The former auto value is no longer accepted. Remove HERCULES_PLANE_IO_MODE=auto from existing job scripts to use the MPI-IO default, or set legacy to retain CPU interpolation.",
      },
    ],
  },
  physics: {
    eyebrow: "Run simulations",
    title: "Physics models",
    summary:
      "Configure linear and nonlinear/plastic solids, acoustic fluids, saturated porous media, material attenuation, coupled interfaces, and absorbing layers with complete parameters.in examples.",
    sections: [
      {
        id: "models",
        title: "Supported formulations",
        table: {
          headers: ["Model", "Primary fields", "Typical use"],
          rows: [
            ["Elastic", "Three-component displacement", "Regional and site-scale seismic wave propagation"],
            ["Nonlinear/plastic solid", "Displacement, stress, plastic strain, and constitutive state", "Yielding, permanent deformation, and nonlinear soil/site response"],
            ["Acoustic", "Scalar pressure", "Water layers and fluid domains"],
            ["Poroelastic u-p", "Solid displacement and pore pressure", "Saturated soils and fluid–solid interaction"],
          ],
        },
      },
      {
        id: "nonlinear-plasticity",
        title: "Nonlinear setup rules",
        paragraphs: [
          "Set include_nonlinear_analysis=yes and supply every nonlinear key shown in the examples below. With Istanbul_velocity_model=no, the first and last Vs_limit rows define the inclusive nonlinear-element range; elements outside that range remain linear. The two identical rows in these syntax examples therefore select solids with 150 <= Vs <= 450 m/s without introducing Vs-dependent interpolation.",
          "The legacy Istanbul nl_id classification controls element eligibility only when Istanbul_velocity_model=yes. Other runs ignore nl_id and use the Vs range covered by material_properties_list. For non-Istanbul vonMises_baE and vonMises_ep, c_or_Su supplies strength; the built-in nl_id constants remain limited to the Istanbul path.",
          "The displayed constants are synthetic regression or syntax values, not production soil calibrations. Replace strengths, angles, backbone coefficients, tolerances, substep limits, damping, and property rows with values justified for the case, then confirm the nonlinear-element count and selected backend in the startup report.",
          "Pressure-dependent examples—MohrCoulomb, DruckerPrager, ManzariDafalias, PDMY03, RIVASand, RIVASandCustom, RIVASand02, and RIVASandCustom02—set approximate_geostatic_state=yes. Hercules integrates the 1-D density column to initialize a K0 effective-stress field. In active BrickUP elements it subtracts hydrostatic Biot stress, while the pressure DOF starts at zero as excess pore pressure. Because the parser forbids combining that approximation with a positive geostatic loading or cushion time, those examples set both times to zero. A calibrated ramped geostatic stage or valid restart is an alternative, but then approximate_geostatic_state must be no.",
          "material_plasticity_type=rate_independent is the established production choice illustrated here. stiff_damp=0 disables the nonlinear stiffness-damping contribution while freq_stiff_damp remains positive. tension_cutoff=yes is accepted only for MohrCoulomb and DruckerPrager.",
        ],
        table: {
          headers: ["Family", "Available material_model values", "Initialization shown below"],
          rows: [
            ["Elastic diagnostic", "linear", "No initial stress"],
            ["von Mises / cyclic plasticity", "vonMises_ep, vonMises_fa, vonMises_faM", "No initial stress in the syntax fixture"],
            ["Nonlinear backbone", "vonMises_baE, vonMises_baH, vonMises_GQH, vonMises_MKZ, vonMises_RO", "Model-specific; BAH and MKZ examples use approximate overburden"],
            ["Frictional plasticity", "MohrCoulomb, DruckerPrager", "approximate_geostatic_state=yes"],
            ["Critical-state sand", "ManzariDafalias", "approximate_geostatic_state=yes"],
            ["Nested-surface liquefaction", "PressureDependMultiYield03 (PDMY03)", "approximate_geostatic_state=yes"],
            ["Bounding-surface liquefaction", "RIVASand, RIVASandCustom", "approximate_geostatic_state=yes"],
            ["Research successor", "RIVASand02, RIVASandCustom02", "approximate_geostatic_state=yes"],
          ],
        },
        note: "Each block is a nonlinear portion of parameters.in. Keep the domain, mesh, source, material-profile, output, and other required controls elsewhere in the same file.",
      },
      {
        id: "model-linear",
        title: "linear diagnostic model",
        paragraphs: [
          "material_model=linear keeps selected elements elastic while producing nonlinear-path diagnostics such as stress, strain, and invariants. Use it as a control run before enabling plastic evolution.",
        ],
        codeLabel: "parameters.in — synthetic diagnostic example",
        code: nonlinearModelExample({
          model: "linear",
          columns: classicMaterialColumns,
          row: "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
        }),
      },
      {
        id: "model-vonmises-ep",
        title: "vonMises_ep",
        paragraphs: [
          "The elastic-perfectly-plastic von Mises model reads its non-Istanbul shear strength from c_or_Su. It uses no isotropic hardening; the example makes the strength constant across the selected Vs range.",
        ],
        codeLabel: "parameters.in — synthetic syntax, not a calibration",
        code: nonlinearModelExample({
          model: "vonMises_ep",
          columns: classicMaterialColumns,
          row: "7.5e4 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
        }),
      },
      {
        id: "model-vonmises-fa",
        title: "vonMises_fa",
        paragraphs: [
          "The Frederick–Armstrong model uses c_or_Su for strength, gamma0 for the initial elastic strain scale, and psi for its kinematic hardening control. The remaining classic-table fields are inactive in this fixture.",
        ],
        codeLabel: "parameters.in — synthetic syntax, not a calibration",
        code: nonlinearModelExample({
          model: "vonMises_fa",
          columns: classicMaterialColumns,
          row: "1.0e5 0 0 0 0 0 2.0e-4 0.25 0 0 0 0 0 0 0",
        }),
      },
      {
        id: "model-vonmises-fam",
        title: "vonMises_faM",
        paragraphs: [
          "The modified Frederick–Armstrong model uses c_or_Su and psi and includes its load-reversal recalibration. gamma0 is not used by this formulation.",
        ],
        codeLabel: "parameters.in — synthetic syntax, not a calibration",
        code: nonlinearModelExample({
          model: "vonMises_faM",
          columns: classicMaterialColumns,
          row: "1.0e5 0 0 0 0 0 0 0.25 0 0 0 0 0 0 0",
        }),
      },
      {
        id: "model-vonmises-bae",
        title: "vonMises_baE",
        paragraphs: [
          "The Borja–Amies exponential model reads non-Istanbul strength from c_or_Su and its hardening law from psi and m. The sample psi and m come from the repository's synthetic regression fixture and are not general soil defaults.",
        ],
        codeLabel: "parameters.in — synthetic syntax, not a calibration",
        code: nonlinearModelExample({
          model: "vonMises_baE",
          columns: classicMaterialColumns,
          row: "7.5e4 0 0 0 0 0 0 0.9331 1.245 0 0 0 0 0 0",
        }),
      },
      {
        id: "model-vonmises-bah",
        title: "vonMises_baH",
        paragraphs: [
          "The Borja–Amies hyperbolic path derives its active strength and hardening quantities from Vs and overburden correlations. The classic-table values below are zeros because only the Vs limits select elements. The example initializes approximate overburden so the pressure-dependent correlation is not evaluated from an unstressed state.",
        ],
        codeLabel: "parameters.in — eligibility and initialization example",
        code: nonlinearModelExample({
          model: "vonMises_baH",
          columns: classicMaterialColumns,
          row: "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
          approximateGeostatic: true,
        }),
        note: "The embedded correlation is part of this model implementation; independently verify that it is appropriate for the site before using it.",
      },
      {
        id: "model-vonmises-gqh",
        title: "vonMises_GQH",
        paragraphs: [
          "The GQH backbone model uses c_or_Su together with theta1 through theta5. The illustrated coefficients exercise the implemented inversion but do not represent a calibrated modulus-reduction curve.",
        ],
        codeLabel: "parameters.in — synthetic syntax, not a calibration",
        code: nonlinearModelExample({
          model: "vonMises_GQH",
          columns: classicMaterialColumns,
          row: "5.0e4 0 0 0 0 0 0 0 0 0 1 1 1 1 0",
        }),
      },
      {
        id: "model-vonmises-mkz",
        title: "vonMises_MKZ",
        paragraphs: [
          "The current MKZ production path derives its reference strain and limiting strength from Vs and overburden and uses the implementation's beta and exponent. The property rows below therefore select the Vs interval; approximate geostatic initialization supplies the in-situ stress used by the correlation.",
        ],
        codeLabel: "parameters.in — eligibility and initialization example",
        code: nonlinearModelExample({
          model: "vonMises_MKZ",
          columns: classicMaterialColumns,
          row: "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
          approximateGeostatic: true,
        }),
      },
      {
        id: "model-vonmises-ro",
        title: "vonMises_RO",
        paragraphs: [
          "The Ramberg–Osgood model reads c_or_Su, alpha_RO from psi, eta_RO from m, and phi_RO from the final column. The values below are the synthetic material-point fixture values.",
        ],
        codeLabel: "parameters.in — synthetic syntax, not a calibration",
        code: nonlinearModelExample({
          model: "vonMises_RO",
          columns: classicMaterialColumns,
          row: "5.0e4 0 0 0 0 0 0 1.0 2.0 0 0 0 0 0 0.5",
        }),
      },
      {
        id: "model-mohr-coulomb",
        title: "MohrCoulomb",
        paragraphs: [
          "Mohr–Coulomb reads cohesion in Pa, friction and dilatancy angles in degrees, and the hardening modulus from the classic table. A physically meaningful in-situ stress is essential for frictional yielding, so the example enables approximate geostatic initialization and zeroes the mutually exclusive ramp times.",
        ],
        codeLabel: "parameters.in — synthetic frictional-plasticity syntax",
        code: nonlinearModelExample({
          model: "MohrCoulomb",
          columns: classicMaterialColumns,
          row: "5.0e4 32 0 0 0 5 0 0 0 0 0 0 0 0 0",
          approximateGeostatic: true,
          tensionCutoff: true,
        }),
      },
      {
        id: "model-drucker-prager",
        title: "DruckerPrager",
        paragraphs: [
          "Drucker–Prager uses the same cohesion, friction-angle, hardening, and dilatancy columns as Mohr–Coulomb and converts the angles to its smooth pressure-dependent cone. The example uses approximate geostatic initialization and enables the model's optional tension cutoff.",
        ],
        codeLabel: "parameters.in — synthetic frictional-plasticity syntax",
        code: nonlinearModelExample({
          model: "DruckerPrager",
          columns: classicMaterialColumns,
          row: "5.0e4 32 0 0 0 5 0 0 0 0 0 0 0 0 0",
          approximateGeostatic: true,
          tensionCutoff: true,
        }),
      },
      {
        id: "model-manzari-dafalias",
        title: "ManzariDafalias",
        paragraphs: [
          "Manzari–Dafalias is a pressure-dependent effective-stress critical-state sand model. Production runs currently require the GPU nonlinear backend. The default adaptive explicit Modified Euler scheme provides embedded local-error control; a fixed-substep Forward Euler scheme using the OpenSees scheme-5 constitutive stage is also available for controlled comparisons and increment-sensitivity studies. The 17-column table contains Vs_limit plus all sixteen calibrated model constants; Poisson's ratio and density come from each element.",
          "The example uses strict adaptive Modified Euler, approximate K0 initialization, zero ramp/cushion times, error_tolerance as a yield-residual bound in Pa, and backbone_errTol as the dimensionless local error target. With manzari_dafalias_integration_scheme=forward_euler, manzari_dafalias_forward_euler_substeps defaults to fixed and no_substeps is the exact number of equal parts over each plastic remainder. The opt-in demand policy treats no_substeps as a cap and uses an MD-specific stress-ratio and deviatoric-strain heuristic; it is not the PDMY03 nested-surface criterion and must be checked for the calibrated path. no_substeps=1 reproduces the original OpenSees single-update path under either policy. Forward Euler reports yield drift and does not use backbone_errTol or the recovery policy.",
        ],
        codeLabel: "parameters.in — synthetic Manzari–Dafalias fixture values",
        code: nonlinearModelExample({
          model: "ManzariDafalias",
          columns: "Vs_limit G0 e_init Mc c lambda_c e0 ksi P_atm(Pa) m h0 ch nb A0 nd z_max cz",
          row: "125 0.80 1.25 0.712 0.019 0.934 0.70 101300 0.01 7.05 0.968 1.10 0.704 3.50 4.0 600",
          approximateGeostatic: true,
          errorTolerance: "1.0e-1",
          backboneTolerance: "1.0e-4",
          extraControls: [
            "manzari_dafalias_integration_scheme = modified_euler",
            "manzari_dafalias_recovery = strict",
          ],
        }),
        note: "Hercules supplies no implicit Manzari–Dafalias calibration. Replace every displayed constant and validate drained and cyclic response independently.",
      },
      {
        id: "model-pdmy03",
        title: "PressureDependMultiYield03 (PDMY03)",
        paragraphs: [
          "PDMY03 is the device-first liquefaction model with standard automatically generated nested yield surfaces. Production runs currently require the GPU nonlinear backend. Its 22-column table contains Vs_limit plus 21 model values; stress-like quantities are Pa, angles are degrees, loading_type is 0, 1, or 2, and surface_count is an integer from 2 through 40.",
          "The example initializes the physical effective-stress field with approximate_geostatic_state=yes. no_substeps caps the deterministic surface- and strain-based subincrement request; backbone_errTol and error_tolerance are still parser-required but do not control the PDMY03 integrator.",
        ],
        codeLabel: "parameters.in — synthetic PDMY03 fixture values",
        code: nonlinearModelExample({
          model: "PressureDependMultiYield03",
          columns: "Vs_limit Gref(Pa) Kref(Pa) phi(deg) peak_gamma pref(Pa) d PT(deg) loading_type ca cb cc cd ce da db dc surfaces liq1 liq2 Patm(Pa) cohesion(Pa)",
          row: "9.0e7 2.2e8 32 0.10 1.01e5 0.50 26 0 0.067 0.50 0.23 1.0 0.0 0.06 3.0 0.27 20 1.0e4 0.01 1.01e5 1.73e3",
          approximateGeostatic: true,
          errorTolerance: "1.0e-1",
          backboneTolerance: "1.0e-4",
        }),
        note: "PDMY03 accepts the alias material_model=PDMY03. The shown constants characterize a regression fixture and are not production liquefaction parameters.",
      },
      {
        id: "model-riva-sand",
        title: "RIVA-Sand",
        paragraphs: [
          "RIVA-Sand means Reversible–Irreversible Volumetric Accumulation model for Sand; RIVASand is the input identifier. It is the native binary64 port of the frozen Ottawa F65 pressure-dependent bounding-surface liquefaction calibration. Production updates are CUDA-only. The kernel retains all 112 calibrated parameters, four host-increment reversal controls, and 39 logical history fields; the startup report prints calibration SHA-256 9585f0c155c9444885c5115e7753a1a5d97c783e2c685938a49a65435c9e8f83.",
          "RIVASand and RIVASand02 use compact Vs_limit Dr rows by default. The actual element density and mesh Vs set the target initial shear modulus rho*Vs^2; Vs_limit only selects/interpolates material rows. An optional third column G0 overrides automatic stiffness. G0 is dimensionless, positive for an explicit override; omitted G0 or G0=0 selects automatic stiffness. Standard names keep the other physical inputs at their reference values. Use one uniform row width and stiffness mode throughout the table. no_substeps remains the fixed explicit subdivision count.",
          "With BrickUP, RIVASand owns effective skeleton stress and its constitutive histories. Its volumetric compatibility is the skeleton relation between total volumetric strain, elastic confining strain, and reversible/irreversible plastic strain; it is not a second fluid-pressure degree of freedom. BrickUP separately enforces d(zeta)=alpha*d(eps_v)+d(p_f)/M, owns storage, Darcy flow and drainage, and contributes -alpha*p_f*I once to total stress. Hercules never copies RIVA-Sand's standalone inferred ru or effective-pressure loss into p_f.",
          "The reversible static-bias volume law uses a bounded low-bias continuation below one quarter of its frozen reference bias. A cubic smoothstep gate is raised to the third power, which cancels the calibrated inverse-square bias factor and makes the contribution vanish continuously as projected bias tends to zero; the original equation is unchanged at and above the onset. This is a contained equation revision, not a parameter recalibration. It removes the zero-bias singularity but does not establish validity at arbitrarily small cyclic amplitude or very low effective confinement; those regimes remain research limitations requiring independent characterization before production use.",
        ],
        codeLabel: "parameters.in — compact RIVA-Sand with mesh stiffness",
        code: nonlinearModelExample({
          model: "RIVASand",
          columns: "Vs_limit Dr",
          row: "0.66296296296296309",
          approximateGeostatic: true,
          errorTolerance: "1.0e-1",
          backboneTolerance: "1.0e-5",
          noSubsteps: 1,
        }),
        note: "Dr is the only soil input in this example. Append G0=483.48301127222084 to each row for the former stiffness calibration instead of mesh-derived stiffness. Automatic stiffness requires valid small-strain mesh Vs, bulk density, and compressive geostatic effective stress.",
      },
      {
        id: "model-riva-sand-intermediate-bias-research",
        title: "RIVASand02 — selectable reversal reference",
        paragraphs: [
          "RIVASand02 selects the native intermediate-bias research successor. It is a distinct internal material sourced from OpenSees_RIVA commit 637d207e867b863dee7b2bf14d9aeb75fdc32e7b, with the plastic-activity gate from commit 5e4ebd30cb5f43b677c61306c7f1fbbdcb24905c; it is not an alias for, replacement for, or retuning of production RIVASand or RIVASandCustom.",
          "RIVASand02 uses Vs_limit Dr [G0]; RIVASandCustom02 exposes Vs_limit Dr G0 M kd h m zeta e_max e_min Q R n_G. The same effective parameters produce the same response in either name. All active RIVA-Sand names store one interpolated material block, initial void ratio and target shear modulus per element, plus one cached elastic scale per integration point. Custom density bounds are used for both initial and reference density.",
          "The material returns effective skeleton stress. In a coupled BrickUP element, BrickUP remains the sole owner of pore pressure, storage, flow, drainage, and the one-time pore-pressure contribution to total stress. This research identifier should remain confined to verification or explicitly reported research runs until it has completed the project validation process.",
          "Its standard material-point fixture uses a dedicated CUDA object and the same allocation-free constitutive core on host and device. The runner compares the complete 139-value state contract, geostatic-admission and reversal selectors, exact zero-increment preservation, and the independent sealed uniform BrickUP pressure history.",
          "Both RS02 names default to rivasand02_reversaltype = 3 and rivasand02_biasvolume = 1 in parameters.in and material-point inputs. Reversal type 1 uses the preceding valid increment, type 2 uses accumulated loading-branch strain, and type 3 uses the full-host elastic stress-ratio predictor from the committed mechanical state, ported from OpenSees_RIVA 355d8f042 with reversalGuard=0. Detection occurs once per host increment before material substeps; type 3 applies at most one event in the first substep and retains the strain rule before cyclic activation. Explicit types 1 and 2 remain available. BiasVolume mode 0 preserves the original bias-driven reversible volume, mode 1 disables it while preserving the separate phase-transformation volume, and mode 2 applies the upstream gated correction to the contractive mean part. Mode 2 adds a committed plastic-activity history, updated once per host increment. These options do not establish timestep independence or field validation; reversalGuard and reversalLatch remain absent.",
          "Automatic stiffness uses the full initialized shear-modulus evaluator after stress-preserving geostatic admission, including pressure guards and initial bias multipliers. It uses approximate K0 stress for approximate_geostatic_state=yes and solved effective stress at dynamic activation for the ramped no workflow. The elastic scale is frozen during shaking and saved in checkpoints. Original models retain marker 4. RS02 checkpoint marker = 10 + 3*BiasVolume + (reversaltype-1); the defaults use marker 15. The 139-value state appends field_bias_mean_activity, which stays zero in modes 0/1. All old RS02 markers through 9 and changes to either selector on restart are rejected. Keep model, mesh, and material inputs unchanged. The retired RIVASandIntermediateBiasResearch input is not restored. To reproduce an earlier response from fresh initialization, explicitly select BiasVolume mode 0 and the earlier reversal type.",
          "On a rejected RS02 GPU update, the solver writes rivasand02_failure.rank<R>.step<S>.txt in the run directory before aborting. Its RIVASAND02_FAILURE_V2 records include both selectors, cached stiffness, incoming tensor strain increment, material properties, flags, and 139 state values. The snapshot is downloaded only after failure and is not a restart checkpoint; a repeated failure at the same rank and step overwrites it. The printed pressure and compatibility residual can describe the last accepted state rather than the rejected trial.",
        ],
        codeLabel: "parameters.in — research verification syntax",
        code: nonlinearModelExample({
          model: "RIVASand02",
          columns: "Vs_limit Dr",
          row: "0.663",
          approximateGeostatic: true,
          errorTolerance: "1.0e-1",
          backboneTolerance: "1.0e-5",
          noSubsteps: 1,
        }),
        note: "Omitting the selectors selects reversal type 3 and BiasVolume mode 1. Set rivasand02_reversaltype = 2 and rivasand02_biasvolume = 0 on separate lines to reproduce the earlier loading-branch behavior from fresh initialization. Neither selector is a material-row column. Material, geostatic, tolerance, substep and damping values remain case-specific.",
      },
      {
        id: "model-riva-sand-custom02",
        title: "RIVASandCustom02 — configurable intermediate-bias successor",
        paragraphs: [
          "Custom02 uses its own e_max/e_min for both the initial and reference relative densities. In particular, the shakedown reference void ratio e=0.601 must remain the same anchor when those bounds change. This applies to shakedown hardening, shear compliance, dense-state weighting, and the phase-volume transition; it is not a new calibration parameter. The cached reference density is used only when the material has the reference bounds.",
          "The material-specific density correction after 9203e02 is retained. The same reversal and BiasVolume defaults apply to Custom02. Markers 10 through 18 identify the nine reversal/BiasVolume combinations; default marker 15 identifies type 3 with mode 1. Older 02 checkpoints must be rerun.",
          "The column order and parameter meanings match RIVASandCustom, but the numerical kernel is RIVASand02, not the original RIVA-Sand kernel. Positive G0 values and continuous material fields are interpolated linearly by Vs. Set every G0 entry to zero for automatic mesh stiffness; never mix zero and positive G0 rows. Repeat the physical parameters across rows when only Dr should vary. Both 02 interfaces use the same effective-skeleton/BrickUP coupling and fixed no_substeps behavior.",
          "Require finite inputs, 0<=Dr<=1, G0>=0, M,kd,h,Q>0, m,zeta,R>=0, e_max>e_min>=0 containing e=0.601, and 0<=n_G<=1. Vs_limit is in m/s; all twelve following columns are dimensionless. G0=0 requests automatic stiffness, not zero modulus. M is supplied directly, not phi_cs. The non-exposed internal controls remain those of the 257-parameter reference record, SHA-256 16c8567fd2370c562899b9ec00cc2fc313be137f239cde5d0e22cc1cd3eee800. The hash identifies that reference record, not an arbitrary custom row. BiasVolume is a separate scalar control; h_p0 and reversalLatch remain unexposed.",
        ],
        codeLabel: "parameters.in — Custom02 reference-row equivalence",
        code: nonlinearModelExample({
          model: "RIVASandCustom02",
          columns: "Vs_limit Dr G0 M kd h m zeta e_max e_min Q R n_G",
          row: "0.663 483.48301127222084 1.25 1.125 122.44207260468994 0.945 0.025 0.78 0.51 10.0 1.5 0.65",
          approximateGeostatic: true,
          errorTolerance: "1.0e-1",
          backboneTolerance: "1.0e-5",
          noSubsteps: 1,
        }),
        note: "This public interface documentation is not an equation-level formulation manual. Both 02 names retain the existing research qualification limits; a shorter name does not qualify full GPU/MPI field simulations.",
      },
      {
        id: "model-riva-sand-custom",
        title: "RIVA-Sand custom material interface",
        paragraphs: [
          "RIVASandCustom exposes the complete Vs_limit Dr G0 M kd h m zeta e_max e_min Q R n_G row for the original generation, unlike the compact standard name. Positive G0 sets reference stiffness independently of Dr; G0=0 requests automatic mesh stiffness. G0=483.48301127222084 preserves the previous reference stiffness. M, kd, h, m, zeta, e_max, e_min, Q, R and n_G retain their existing meanings; no cyclic equation or calibration is replaced.",
          "Hercules linearly interpolates continuous inputs after Vs_limit. Use the same automatic or explicit G0 mode in all rows. Repeat the eleven physical fields across rows for a uniform soil calibration. Inputs must satisfy 0<=Dr<=1, G0>=0, M,kd,h,Q>0, m,zeta,R>=0, e_max>e_min>=0, and 0<=n_G<=1. The frozen shakedown law still refers to e=0.601, so that value must lie within the custom void-ratio bounds.",
          "Only the listed existing RIVA-Sand fields are configurable. Reversible compatibility, amplitude scaling, static-bias and ratchet laws, reversal detection, and the remaining empirical controls retain the frozen reference values and SHA. In particular this interface does not include the experimental pre-name V8.1 h_p0 contraction-rate parameter. A custom row needs independent calibration and validation; the supplied reference row is a parser/kernel equivalence fixture, not a new sand calibration.",
          "RIVA-Sand material data remains outside the common constants of unrelated nonlinear families. The automatic calibration runs only at first dynamic initialization, uses total bulk mass density rather than buoyant density, and matches initialized elastic shear stiffness. Bulk stiffness follows the existing fixed nu=0.30; this does not also fit the mesh Vp. Near-zero effective pressures remain subject to geostatic admission checks. Very large derived G0 should prompt review of the supplied velocity/stress profile, not silent clipping. Initial Vs matching does not establish improved strong-motion damping or field accelerations.",
        ],
        codeLabel: "parameters.in — custom RIVA-Sand reference-equivalence syntax",
        code: nonlinearModelExample({
          model: "RIVASandCustom",
          columns: "Vs_limit Dr G0 M kd h m zeta e_max e_min Q R n_G",
          row: "0.66296296296296309 483.48301127222084 1.25 1.125 122.44207260468994 0.945 0.025 0.78 0.51 10.0 1.5 0.65",
          approximateGeostatic: true,
          errorTolerance: "1.0e-1",
          backboneTolerance: "1.0e-5",
          noSubsteps: 1,
        }),
        note: "RIVASand and RIVASandCustom now accept the same full row. Use RIVASand02 or its Custom02 alias only to request the distinct successor.",
      },
      {
        id: "coupling",
        title: "Coupled-interface rules",
        paragraphs: [
          "Material classification determines which formulation advances each element. Acoustic–elastic faces exchange normal traction and normal acceleration. Poroelastic–elastic contact shares solid motion but stops Darcy flux. Acoustic–poroelastic contact can be no-flow or open-pore, with the latter weakly enforcing pressure continuity.",
          "For nonlinear BrickUP regions, the constitutive model owns the skeleton effective stress and the u-p operator supplies Biot pore-pressure traction. Hercules disables the zero-permeability undrained fast path whenever nonlinear analysis is active so alpha²M is not folded into a linear bulk stiffness. Nonlinear elements below poroelastic_depth_max_m are single-phase nonlinear solids, not coupled liquefaction elements.",
        ],
        note: "Check the startup classification counts and station_mediums.txt. A syntactically valid case can still encode the wrong hydraulic or mechanical interface assumption.",
      },
      {
        id: "coupling-acoustic-elastic",
        title: "Acoustic–elastic example",
        paragraphs: [
          "With include_acoustic=yes, a material sample with Vp>0, rho>0, and Vs<=1e-6 is acoustic. The following layered profile places lossless water above an elastic solid and enables the recommended velocity-potential coupling path.",
          "acoustic_interface_face_filter defaults to all: every detected acoustic–elastic interface is coupled using mesh-face normals. The default alias also selects all. Explicit bathymetry couples only downward Z+ water-bottom faces with reconstructed slope normals and surface-area weights, omitting side and upward Z− faces; bottom_only retains Z+ faces with vertical normals. All three filters completed the source-corrected 100 s Izmir comparison with finite output and identical generated forces; slow acoustic motion and reconstruction accuracy remain under investigation. Check that the selected geometry represents the intended water–solid boundary; stability alone does not establish waveform accuracy.",
          "With seabed_geometry_file, acoustic free surfaces follow local topography over water (flat solver Z=0 with include_topography=no): complete horizontal top faces must coincide with that surface, and separate flat water bodies may have different elevations. Without explicit seabed geometry, the flat acoustic free surface follows solver Z=topographybase_zcoord when topography is enabled, surface_shift_m with buildings, or the octree top otherwise. Topography faces must also be exposed at their centers, excluding buried water beneath land. The Acoustic free-surface geometry startup line identifies local topography or the legacy datum and reports the exact exposed face count; the node-copy count includes MPI replicas. This selects both pressure-release nodes and optional gravity surface inertia, including its lateral CPML terms. Earlier builds checked only the octree top and missed shifted water surfaces; affected cases need a new simulation. Curved or cut-cell water surfaces require a different operator. Acoustic interface forcing now waits for MPI displacement sharing and hanging-node assignment before its scalar update, using the existing CUDA-aware or host exchange. Earlier builds could become unstable on adaptive meshes split across ranks while a one-rank run stayed stable. The Acoustic interface displacement startup line reports this ordering. Setup also fetches actual neighbor leaf geometry and corrected materials from their owning MPI ranks; using remote octree aggregates previously changed interface subdivision with the partition. The Acoustic interface geometry line reports exact face and fragment counts. This batched lookup adds no per-step exchange. A recorded local shifted-surface case checked surface enrollment, bounded excited fields, and station/plane agreement at a CPML face. Compare one and four ranks as well as I/O and communication backends when validating a setup. The investigation fixture is not distributed with the repository; its historical results do not establish full wavefield validation.",
        ],
        codeLabel: "parameters.in — water over elastic solid",
        code: "include_acoustic = yes\nacoustic_formulation = velocity_potential\nacoustic_free_surface_gravity = no\nacoustic_interface_coupling = yes\nacoustic_interface_projection = point\nacoustic_interface_face_filter = all\nacoustic_interface_max_level = -1\nacoustic_mesh_refinement_factor = 1.0\nacoustic_quality_factor = 0\n\ncvmdb_input_file = profile\nnumber_profile_layers = 2\nprofile_layers =\n# z_top(m) Vp(m/s) Vs(m/s) rho(kg/m3) Qp Qs\n 0.0   1500.0    0.0 1000.0 1.0e6 1.0e6\n 100.0 2600.0 1400.0 2300.0 1.0e6 1.0e6",
        note: "The default pressure-release surface is appropriate for the existing hydroacoustic path. Set acoustic_free_surface_gravity=yes when linear surface-gravity waves are part of the modeled physics; lateral acoustic CPML then includes the tangential gravity-surface correction. The complete operator now uses independent nodal surface memories and the full tangential x/y product, including corners. Earlier candidate and doubled-CPML runs contained invalid source forces from a tabulated-source endpoint overread, so their nonfinite output cannot establish an operator instability. Regenerate and verify source forces with the corrected interpolator before comparing boundary or interface choices. Finite output alone does not establish CPML absorption or gravity-surface reconstruction accuracy. Use face_mass projection only as an experimental comparison path. A positive acoustic_quality_factor is a numerical high-frequency regularizer, not a physical ocean attenuation law.",
      },
      {
        id: "acoustic-quality-factor",
        title: "Acoustic quality-factor default and validation",
        paragraphs: [
          "acoustic_quality_factor defaults to 0, which disables interior acoustic damping. This default was verified in the documented flat-water CPS comparisons and MPI/backend consistency tests described below. Those results apply to the recorded configurations and durations.",
        ],
        codeLabel: "parameters.in — default acoustic quality factor",
        code: "acoustic_quality_factor = 0",
      },
      {
        id: "acoustic-interface-method",
        title: "Choose interface constraints or local mesh matching",
        paragraphs: [
          "acoustic_interface_method=field-specific-constraints is the default for coupled pressure or velocity-potential acoustics. An omitted method uses shared-geometric-constraints when acoustics, scalar acoustics or interface coupling is disabled; explicit unsupported choices still abort. Coupled VT topography must explicitly select shared-geometric-constraints. That option preserves the original geometric hanging-node constraints for acoustic and solid fields; regular Hex8 pore-pressure constraints are independent. Changing this method does not select a quadrature: acoustic_interface_projection=point or face_mass remains independent, as do the face filter, CPML and free-surface gravity.",
          "Experimental field-specific-constraints retains a hanging-node constraint only where a coarser element of that same field requires it. A water node constrained only by a coarse solid becomes an independent acoustic degree of freedom, and conversely for solid nodes next to coarser water. Setup resolves actual neighboring leaves across MPI partitions and checks retained weights against the same-field Hex8 trace. An unsupported mixed-corner stencil stops setup; it does not silently revert to the original map. Coupling uses actual solid face vertices for both point and face_mass, including remote solid neighbors. Missing remote vertices become non-owned harbored nodes, with global-ID ordering in the existing force-contribution and displacement-sharing schedules. This avoids replacing independent fine-solid nodes with the coarse water trace at MPI boundaries, with or without topography. The monitor reports keep/release/absent counts for each field, remote fragments and additional harbored nodes. The existing CUDA-aware and host-staged MPI paths use these maps without an additional per-step exchange round.",
          "Experimental local-matching refines the coarser side of every detected water–solid interface until opposing faces have the same element size. Each pass rebalances the mesh and rechecks the actual material classification. Same-material regions remain adaptive, although balancing and connected interface faces can spread refinement. This is automatic face matching, not a prescribed uniform-width refinement band. It uses geometric acoustic/solid constraints on the resulting mesh; regular Hex8 pore pressure retains its independent map. Set it before Phase A and regenerate the mesh and source caches. Phase B imports do not refine: an incompatible imported mesh is rejected. No extra matching level or band width is required. The monitor lists refinement passes and completion.",
          "A local 58-run pressure validation (2026-10-08 UTC), with acoustic_quality_factor=0 in every run, covered point/face_mass, finer water, gravity, combined buildings/nonlinear elastoplastic/BrickUP open-pore fixtures, one/eight MPI ranks, host/CUDA-aware exchanges, and a six-second FEM Izmir crop at eight/sixteen ranks. Fourteen 110-second flat-water comparisons used an independent source-matched CPS reference and raw simulated responses: median submerged pressure error was about 18.38% coarse and 8.34% fine with either formulation and either newly enabled method. Halving the coarse time step changed it to about 18.40%. Matched formulation and MPI/backend comparisons agreed below 1e-6 aggregate relative L2. The local-only fixtures are not shipped; the short gravity/module/crop comparisons establish implementation consistency, not independent gravity-wave accuracy, calibrated constitutive response, or full-duration Izmir stability.",
          "Both scalar formulations use the same constraint maps and mesh matching. Pressure mode couples solid normal acceleration into the acoustic equation; velocity-potential mode couples normal velocity and obtains pressure traction from the potential time difference. Enabling pressure does not establish an accuracy advantage or remove time-step and mesh checks. Velocity potential remains the default. acoustic_station_gradient still affects station motion only, not the pressure column.",
          "Both alternatives require include_acoustic=yes, acoustic_interface_coupling=yes, and acoustic_formulation=pressure or velocity_potential with full Hex8 solids. Buildings, eligible nonlinear constitutive models and BrickUP poroelasticity may be enabled separately or together, subject to their existing module restrictions. Nonlinear geostatic initialization, tolerances, substeps, damping and material_properties_list require case-specific calibration. Field-specific constraints also support FEM topography: use include_topography=yes and computation_method=fem while retaining the case's topography datum, directory and etree type. FEM retains regular Hex8 material cells; air cells and owner-confirmed carved voids contribute no field. Unresolved remote neighbors still abort. VT cut-tetrahedron topography remains unsupported, and local matching requires include_topography=no. Material-defined stepped or sloping seabeds remain allowed within the full-Hex8 scope. Local matching inspects all detected water–solid faces even when a force face filter omits some of them and reapplies building carving after each refinement and balance pass. Acoustic interface max level and acoustic volume refinement retain their separate roles.",
          "Regular Hex8 and FEM-topography BrickUP models automatically build an independent pore-pressure map from final active BrickUP cells, excluding elastic-only CPML cells. This also applies with include_acoustic=no and is independent of acoustic_interface_method. VT cut-cell topography retains geometric pressure constraints and reports that limitation at startup. Saturated solids carry displacement and pore pressure; dry solids and buildings carry displacement. The pore-pressure map handles capacity, drained and open-pore boundary weights, RHS distribution and pressure assignment on host and CUDA-aware MPI paths. Both alternative methods use actual remote solid vertices. All three methods, including shared-geometric-constraints, resolve actual remote solid vertices for face_mass and use the same four-point quadrature on local and remote faces. The shared-geometric method retains its geometric hanging-node constraints for acoustic and solid fields. open_pore also builds actual solid traces with shared-geometric point projection so hydraulic quadrature uses the solid-side shape weights on remote faces, including independent pore nodes. Sealed shared-geometric point coupling keeps its existing path. Trace-enabled ranks use the same global-ID schedule ordering, including ranks that only own requested nodes; coupling-off and non-acoustic cases skip trace setup. Open-pore boundary mass and water-pressure targets also use actual solid traces. Open-pore hydraulic projection uses four-point (2x2 Gauss) quadrature for both boundary capacity and the water-pressure target, independently of mechanical point/face_mass quadrature. These setup changes do not add per-step MPI exchange rounds.",
          "Building meshing handles layered profiles by classifying building/foundation cells and air before querying the soil profile. This avoids invalid negative-depth queries above surface_shift_m. Air uses a positive meshing speed before carving, including Vp when the surface is water with Vs=0. Building and foundation properties follow actual cell depth. BrickUP ground and groundwater depths honor surface_shift_m: above-ground building cells remain outside the saturated region, while foundation cells follow the existing groundwater and depth-limit settings. One-rank ownership queries also work without an MPI interval table.",
          "The preceding fixed-geometry sloping-seabed prototype with field-specific-constraints, face_mass, gravity and CPML remained bounded for 80 s at two resolutions. Its separate dashpot comparisons did not establish improved accuracy. Earlier band-matched sloping tests showed exponential growth with CPML, including face_mass; a stable flat-water comparison does not establish stability for arbitrary geometry. Neither experimental mode is a generally validated replacement, and bounded output alone is not a waveform-accuracy test. The production local-matching algorithm differs from that earlier prescribed band.",
          "Horizon integration checks on 2026-10-03 completed 8 s with all three methods, CPML, gravity and face_mass for both water-finer and solid-finer meshes. Matching meshes imported successfully; incompatible meshes and invalid settings were rejected. Host and CUDA-aware field-map histories agreed within 5.8e-9 relative L2. The new parameter reproduced the coarse sloping prototype over 80 s, with identical constraint counts and pressure differing by 1.4e-9 relative L2. These are implementation comparisons, not independent accuracy references.",
          "Those earlier checks predate the remote-solid trace correction and did not establish MPI partition invariance. On 2026-10-04, 8 s FEM-topography checks covered a shifted water surface with air above, CPML, gravity on/off, both refinement directions, point and face_mass, and mesh import. One versus eight GPU ranks and host versus CUDA-aware MPI agreed within 2e-11 relative L2 across tested plane and station histories. The previous remote-solid fallback gave a 1.43% plane-history discrepancy in the solid-finer case. The corrected code reproduced the original one-rank topography-free operator and shared-geometric baseline. These checks establish implementation consistency within that scope, not continuum accuracy or arbitrary-coastline long-time stability.",
          "The same build imported the existing 1,179,872,281-element Izmir mesh on 512 Horizon GPUs and completed 20 s with FEM topography, gravity and the original physical/source controls. All 18 stations recorded nonzero, finite motion (largest displacement component 0.1384 m); station histories and generated source forces passed integrity checks. Setup resolved 34,534 remote solid fragments using 12,310 additional harbored nodes. This is a short-run integration check, not validation of the full 100 s production simulation.",
          "The subsequent module extension completed 66 candidate runs and six historical solver controls on Horizon, with 49 accepted physical-output comparisons. Four-second synthetic cases covered both methods, buildings, yielding vonMises_ep solids, explicit BrickUP no-flow/open-pore and combined modules, both refinement directions, point/face_mass, one/eight GPU ranks, host/CUDA-aware MPI, nonlinear CPU/GPU execution, drained pressure with finite Darcy permeability, and matching-mesh import. The final build also replayed the earlier 8 s FEM/gravity baselines. Responses were unfiltered. Maximum relative L2 differences were 2.02e-9 per plane-motion component and 2.50e-8 per physical station column. The existing pre-arrival CPU skip causes a diagnostic Fs difference; that column was excluded only from CPU/GPU comparisons, with motion, pressure, strain and stress retained. The synthetic materials are not production calibrations, and this does not validate every constitutive law or long-time stability.",
          "The partition-crossing local-matching plus face_mass control improved from a 3.33% worst-component plane discrepancy between one and eight ranks to 1.85e-13 relative L2. The corrected one-rank result reproduced the original one-rank operator. A final follow-up guards an absent building-carving timer at shutdown for imported meshes; six import checks and the two saved FEM/gravity regressions passed on that final build. Failed or incorrectly specified attempts were excluded and replaced by corrected runs.",
          "The later shared-geometric-constraints plus face_mass correction explicitly compared one/eight Horizon GPUs on 4 s flat-water cases. Worst plane-component relative L2 discrepancies changed from 3.33% to 1.93e-13 for matching faces, from 2.78% to 1.81e-13 with solids twice as fine, and from 1.54% to 7.67e-14 with water twice as fine. Corrected one-rank results reproduced the previous one-rank operator. Raw, unfiltered checks also covered host/CUDA-aware MPI, the unchanged default point path, combined building/nonlinear/open-pore modules, FEM/gravity, legacy pressure acoustics, coupling disabled, and both alternative methods. This correction was not included in the earlier production-sized Izmir check.",
          "A combined yielding/open-pore case with finer solids exceeded the strict host/CUDA tolerance in host performance mode: 4.67e-8 plane-component relative L2 and 3.84e-5 in a station acceleration column. Host repeats varied, and the previous solver also showed a host/CUDA difference on the same one-rank model. With the existing HERCULES_REPRO_MODE=strict, host/CUDA and one/eight-rank comparisons agreed within 1.5e-12 for plane motion and 4.6e-9 for station columns. Performance-mode failures remain recorded as a repeatability limitation, not counted as strict passes; the correction does not change host execution defaults.",
        ],
        codeLabel: "parameters.in — experimental choice, other model inputs still required",
        code: "include_acoustic = yes\nacoustic_formulation = velocity_potential\nacoustic_interface_coupling = yes\nacoustic_interface_method = field-specific-constraints\n# Default: field-specific-constraints; alternatives: shared-geometric-constraints, local-matching\nacoustic_interface_projection = face_mass\nacoustic_interface_face_filter = all\nacoustic_free_surface_gravity = yes\nabc_type = cpml\ncpml_clamp_outer = yes\ninclude_topography = no\ninclude_buildings = no\ninclude_nonlinear_analysis = no\ninclude_poroelastic = no",
        note: "For a controlled comparison, keep source, timestep, absorber thickness in metres, gravity and output locations fixed. Inspect matching-pass counts and the resulting element count; local matching can increase memory and setup cost substantially. Recheck the stable timestep after remeshing.",
      },
      {
        id: "acoustic-station-gradient",
        title: "Recover acoustic motion at station coordinates",
        paragraphs: [
          "Acoustic station displacement, velocity and acceleration are reconstructed from pressure gradients. A single Hex8 derivative is constant along the differentiated direction within its element, so a near-bed station can show a large relative motion error even when the nodal potential is accurate. Refining the mesh changes the station's position within its containing cell and need not improve that pointwise derivative monotonically.",
          "acoustic_station_gradient=linear is the default. At setup it combines the containing cell derivative with an adjacent same-material water-cell derivative to estimate the gradient at the actual station coordinate. It prefers the neighbor toward the station from the containing cell center, and uses one-sided extrapolation when only the opposite neighbor is suitable. This recovery is exact for quadratic scalar fields sampled at Cartesian nodes, including unequal cell sizes; it is an output approximation, not a higher-order propagation scheme or a correction of constraint-induced field errors.",
          "No recovery crosses a water material jump (Vp or density), solid or void. If no suitable neighbor exists for a component, that component keeps the element derivative. CPML stations retain the existing element gradient and its inverse-stretch memory, and interior stencils exclude CPML neighbors. Remote stencil vertices use existing host or CUDA-aware field-sharing schedules. Sender and receiver mappings use global node ID order for every interface method, including shared-geometric-constraints with point projection. Weights are precomputed and GPU station gathers remain device-resident.",
          "Select element to reproduce the historical station gradient or to compare methods on an identical mesh and source. Pressure columns, plane reconstruction, elastic/poroelastic station processing, and the physical operators do not change. ASCII and binary acoustic station displacements use the same reconstructed motion. A station and a colocated acoustic plane point can therefore differ in linear mode. Check raw observables against a physical reference; the receiver-formula prediction is a separate diagnostic.",
        ],
        codeLabel: "parameters.in — acoustic station motion",
        code: "acoustic_station_gradient = linear\n# Historical single-element gradient: acoustic_station_gradient = element",
      },
      {
        id: "coupling-poroelastic",
        title: "Poroelastic–elastic example",
        paragraphs: [
          "This sealed, zero-permeability configuration treats saturated solids down to 50 m as linear BrickUP elements. The underlying solid remains elastic and the hydraulic interface is no-flow. For a drained surface, change poroelastic_surface_pressure_bc to drained and retain a positive penalty factor. On regular Hex8 and FEM-topography meshes, include_acoustic=no still builds the pore-pressure map automatically. Constraints required only by coarse dry or elastic-only CPML cells are released for pore pressure, while coarse active BrickUP cells retain their exact Hex8 traces. Solid displacement keeps its selected constraints. VT cut-cell models keep geometric pressure constraints with a startup notice. The standalone extension still requires target-system CUDA/MPI propagation and refinement comparisons; earlier acoustic-module checks do not validate this path.",
        ],
        codeLabel: "parameters.in — linear undrained BrickUP layer",
        code: "include_poroelastic = yes\nporoelastic_formulation = brick_up\nporoelastic_depth_max_m = 50.0\nsurface_groundwatertable = no\nporoelastic_surface_pressure_bc = undrained\nporoelastic_surface_pressure_penalty_factor = 1.0\nporoelastic_elastic_interface = no_flow\nporoelastic_undrained_fast_path = yes\nporoelastic_porosity = 0.40\nporoelastic_biot_alpha = 1.0\nporoelastic_fluid_bulk_modulus_pa = 2.2e9\nporoelastic_fluid_density_kgm3 = 1000.0\nporoelastic_fluid_viscosity_pa_s = 1.0e-3\nporoelastic_permeability_x_m2 = 0.0\nporoelastic_permeability_y_m2 = 0.0\nporoelastic_permeability_z_m2 = 0.0\nacoustic_poroelastic_interface = no_flow",
      },
      {
        id: "coupling-nonlinear-poroelastic",
        title: "Nonlinear PDMY03–BrickUP example",
        paragraphs: [
          "The following combines the complete PDMY03 nonlinear fixture with a sealed saturated layer. approximate_geostatic_state=yes supplies the initial skeleton effective stress, the pressure unknown begins at zero excess pressure, and poroelastic_undrained_fast_path=no makes the explicit nonlinear coupling policy visible even though Hercules disables it automatically.",
        ],
        codeLabel: "parameters.in — synthetic coupled syntax, not a calibration",
        code: `${nonlinearModelExample({
          model: "PressureDependMultiYield03",
          columns: "Vs_limit Gref Kref phi peak_gamma pref d PT loading_type ca cb cc cd ce da db dc surfaces liq1 liq2 Patm cohesion",
          row: "9.0e7 2.2e8 32 0.10 1.01e5 0.50 26 0 0.067 0.50 0.23 1.0 0.0 0.06 3.0 0.27 20 1.0e4 0.01 1.01e5 1.73e3",
          approximateGeostatic: true,
          errorTolerance: "1.0e-1",
          backboneTolerance: "1.0e-4",
        })}\n\ninclude_poroelastic = yes\nporoelastic_formulation = brick_up\nporoelastic_depth_max_m = 50.0\nsurface_groundwatertable = no\nporoelastic_surface_pressure_bc = undrained\nporoelastic_undrained_fast_path = no\nporoelastic_porosity = 0.40\nporoelastic_biot_alpha = 1.0\nporoelastic_fluid_bulk_modulus_pa = 2.2e9\nporoelastic_fluid_density_kgm3 = 1000.0\nporoelastic_fluid_viscosity_pa_s = 1.0e-3\nporoelastic_permeability_x_m2 = 0.0\nporoelastic_permeability_y_m2 = 0.0\nporoelastic_permeability_z_m2 = 0.0`,
        note: "A zero-permeability nonlinear BrickUP case must not use a positive geostatic loading ramp. Use the approximate state shown here, a calibrated permeable consolidation stage, or a valid restart.",
      },
      {
        id: "coupling-riva-sand-poroelastic",
        title: "RIVA-Sand–BrickUP example",
        paragraphs: [
          "RIVA-Sand is compatible with the BrickUP element through the standard effective-stress split. The RIVA-Sand kernel returns effective skeleton stress and evolves plastic volumetric compatibility. BrickUP independently advances excess fluid pressure from d(zeta)=alpha*d(eps_v)+d(p_f)/M plus Darcy flux, then contributes -alpha*p_f*I once during total-stress assembly. The constitutive model's standalone effective-pressure-loss diagnostic is not copied into the pressure degree of freedom, so it cannot be counted a second time.",
          "At the first nonlinear update, both approximate-K0 and ramped-gravity workflows pass their compressive effective stress through one admission path. An admissible state is unchanged. An over-bound state keeps its equilibrated stress, recenters its mapping origins, and may not move farther from the bounding cone; noncompressive physical stress is rejected. Hercules initializes every active RIVA-Sand integration point at this transition and reports the global initialized count, admitted count, and maximum admitted cone ratio. The frozen 93-value numerical oracle remains unchanged while the committed admission selector is compared separately. The internal residual cone pressure remains zero and is not a user input.",
          "The checkpoint extension preserves both BrickUP pressure time levels and every RIVA-Sand integration-point history field. RIVA-Sand restarts require both the frozen parameter SHA and equation-level RIVA_KERNEL_REVISION=2. Older equation revisions are rejected rather than being resumed under changed constitutive equations; BrickUP-only checkpoint compatibility is unchanged. A RIVA-Sand/BrickUP restart also refuses a legacy displacement-only checkpoint instead of silently resetting coupled state.",
        ],
        codeLabel: "parameters.in — frozen RIVA-Sand coupled syntax",
        code: `${nonlinearModelExample({
          model: "RIVASand",
          columns: "Vs_limit Dr",
          row: "0.66296296296296309",
          approximateGeostatic: true,
          errorTolerance: "1.0e-1",
          backboneTolerance: "1.0e-5",
          noSubsteps: 1,
        })}\n\ninclude_poroelastic = yes\nporoelastic_formulation = brick_up\nporoelastic_depth_max_m = 50.0\nsurface_groundwatertable = no\nporoelastic_surface_pressure_bc = undrained\nporoelastic_undrained_fast_path = no\nporoelastic_porosity = 0.40\nporoelastic_biot_alpha = 1.0\nporoelastic_fluid_bulk_modulus_pa = 2.2e9\nporoelastic_fluid_density_kgm3 = 1000.0\nporoelastic_fluid_viscosity_pa_s = 1.0e-3\nporoelastic_permeability_x_m2 = 0.0\nporoelastic_permeability_y_m2 = 0.0\nporoelastic_permeability_z_m2 = 0.0`,
        note: "The frozen material constants do not make a site automatically valid. Confirm the nonlinear/poroelastic map intersection, initial effective stress, timestep, pressure boundary, and coupled one-element response before production use.",
      },
      {
        id: "coupling-acoustic-poroelastic",
        title: "Acoustic–poroelastic open-pore example",
        paragraphs: [
          "Enable both formulations and select open_pore only when the acoustic and saturated porous regions physically share an open hydraulic contact. The positive penalty factor scales the existing weak pressure-continuity projection. No new parameter is required for four-point (2x2 Gauss) hydraulic integration; acoustic_interface_projection still controls mechanical traction/motion only.",
          "Each rectangular interface fragment samples four matching physical points using separate acoustic and solid natural-coordinate centers and half widths. Acoustic pressure is interpolated at each point and integrated against the solid shape functions. For fragment compliance c = penalty_factor * area * solid_edge_length / (2 * solid_storage_modulus), boundary capacity adds (c/4) * sum_q N_s(q) and the target adds (c/4) * sum_q N_s(q) * p_acoustic(q). Actual solid traces are required on local and remote neighbors; missing geometry aborts. Capacity and target share the same rule, preserving a constant-pressure target within roundoff. Bilinear trace products with constant fragment compliance are integrated exactly on the represented rectangle.",
          "For a matching face with p_acoustic = 200 + 100*xi Pa, the integrated target/capacity ratio is 166.666… Pa at the xi=-1 nodes and 233.333… Pa at xi=+1, instead of the former 200 Pa centroid target everywhere. Velocity potential retains the existing p = -(phi[Tm2] - phi[Tm1])/dt conversion. Material inputs, penalty scaling and the pressure time update are unchanged. The kernel reuses each acoustic nodal value across the four points and accumulates before nodal atomics. Existing host-staged and CUDA-aware scalar exchanges remain in use without extra per-step MPI rounds.",
          "This upgrades the spatial integration of the existing weak pressure-continuity operator, without adding fluid inertia, a new Darcy interface law, or cut-face/curved-surface reconstruction. Existing mesh-face geometry and effective areas still apply; VT retains geometric pressure constraints. Earlier Horizon open-pore benchmarks used centroid hydraulic integration and do not validate this upgrade. Local analytical host checks cover the production capacity/target algebra, six face orientations, constant and varying pressure, both acoustic formulations, nonmatching fragments and reordered remote trace indices. Target-system CUDA/MPI propagation and performance checks remain required. The startup monitor reports 2x2 Gauss hydraulic quadrature separately from mechanical projection.",
        ],
        codeLabel: "parameters.in — open-pore interface controls",
        code: "include_acoustic = yes\nacoustic_formulation = velocity_potential\nacoustic_interface_coupling = yes\n\ninclude_poroelastic = yes\nporoelastic_formulation = brick_up\nporoelastic_depth_max_m = 50.0\nporoelastic_surface_pressure_bc = undrained\nporoelastic_elastic_interface = no_flow\n\nacoustic_poroelastic_interface = open_pore\nacoustic_poroelastic_pressure_penalty_factor = 1.0",
      },
      {
        id: "attenuation",
        title: "Attenuation and damping rules",
        paragraphs: [
          "type_of_damping selects none, mass-proportional, Rayleigh, or a BKT-family viscoelastic convolution. Rayleigh and mass coefficients are derived from simulation_wave_max_freq_hz and scaled by each element's damping ratio. BKT-family models use material Qp/Qs or a parametric Qs relation and maintain convolution history.",
          "the_threshold_damping caps the element damping ratio, and do_damping_statistics=1 prints initialization statistics. Qp and Qs are the final two fields of profile_layers and are also available from the 3-D material-model Q files.",
        ],
      },
      {
        id: "damping-none-rayleigh-mass",
        title: "None, Rayleigh, and mass damping examples",
        paragraphs: [
          "Copy exactly one of the three type_of_damping selections into a case. The profile excerpt supplies Qp and Qs for the attenuating choices; none ignores them. Rayleigh includes both mass- and stiffness-proportional terms, while mass retains only the mass-proportional term.",
        ],
        codeLabel: "parameters.in — choose one damping selection",
        code: "# Lossless material response\ntype_of_damping = none\n\n# OR: Rayleigh material damping\n# type_of_damping = rayleigh\n\n# OR: mass-proportional material damping\n# type_of_damping = mass\n\nthe_threshold_damping = 0.05\ndo_damping_statistics = 1\n\n# z_top Vp Vs rho Qp Qs\ncvmdb_input_file = profile\nnumber_profile_layers = 1\nprofile_layers =\n 0.0 2000.0 1000.0 2000.0 200.0 100.0",
        note: "Uncomment only one type_of_damping line. Q values and the 0.05 cap are illustrative; use the attenuation model and quality factors supported by measurements for the site.",
      },
      {
        id: "damping-bkt-material-q",
        title: "BKT-family example with material Q",
        paragraphs: [
          "All BKT variants require use_parametric_q_factor and use_infinite_qk. With use_parametric_q_factor=no, profile or 3-D material Qp/Qs values are used. bkt and bkt2 use two relaxation mechanisms; bkt3 and bkt3f use three. Keep the selected variant consistent with the validation used for the project.",
          "BKT setup rejects nonfinite or nonpositive solid Q, invalid bulk modulus, and invalid velocities, density or reference frequency. With use_infinite_qk=no, negative derived Qk beyond roundoff instead selects Qk=1000 (lossless bulk) only for affected elements. Valid bulk and shear attenuation remain unchanged. All four BKT variants use this policy, and setup reports the global negative-Qk bulk fallbacks count. Zero inverse bulk Q remains lossless; use_infinite_qk=yes suppresses bulk attenuation for all solids.",
          "Fallback elements retain shear damping, so P waves still attenuate through their shear contribution; the original supplied Qp is no longer reproduced. For example, Vp=1500 m/s, Vs=1000 m/s, Qs=100 and Qp=200 select the fallback, while Vp=2000 m/s with the same Vs/Qs/Qp retains valid Qk=400. These are illustrative values. The existing velocity-profile frequency correction still runs: shear dispersion is unchanged, and affected elements have no bulk dispersion. No new Vp/Vs clamp is applied.",
          "BKT2 still checks finite, nonnegative coefficients and positive relaxed stiffness (1-a0-a1 > 0) in both branches. Remaining errors report a global count and one element location; the fallback does not repair invalid inputs or inadmissible shear coefficients and does not adjust the timestep. Validation reuses the existing input-reading, material-assignment and imported-mesh setup loops. It adds no material queries, mesh refinement, extra full mesh traversal, or time-step work. Imported BKT2 meshes have their stored coefficients checked, but must be regenerated to apply the fallback or changed Q interpolation. These admissibility checks do not establish stability for every timestep.",
          "The fallback uses common solid material setup, including solids in acoustic–elastic, poroelastic–elastic and poroelastic–acoustic cases. Pure acoustic elements remain lossless. Existing VT, CPML and nonlinear BKT eligibility is unchanged, as are poroelastic pressure evolution and nonlinear elastoplastic/plastic and hysteretic laws. Force and station-stress paths use the stored coefficients without new interface terms or GPU/MPI exchanges. Host coefficient comparisons cover all BKT variants; full GPU/MPI, coupled and nonlinear propagation still require platform benchmarks.",
        ],
        codeLabel: "parameters.in — BKT2 with explicit Qp and Qs",
        code: "type_of_damping = bkt2\nsimulation_velocity_profile_freq_hz = 1.0\nuse_parametric_q_factor = no\nuse_infinite_qk = no\nthe_threshold_damping = 0.05\ndo_damping_statistics = 1\n\ncvmdb_input_file = profile\nnumber_profile_layers = 1\nprofile_layers =\n# z_top Vp Vs rho Qp Qs\n 0.0 2000.0 1000.0 2000.0 200.0 100.0",
        note: "simulation_velocity_profile_freq_hz must not exceed simulation_wave_max_freq_hz. Change bkt2 to bkt, bkt3, or bkt3f only when that formulation has been validated for the case.",
      },
      {
        id: "damping-bkt2-gpu",
        title: "BKT2 GPU execution and comparison",
        paragraphs: [
          "BKT2 uses a fused convolution and force kernel for cube elements with device-resident force and displacement synchronization when an initialization scan confirms uniform, matching shear and bulk relaxation rates. Host-staged execution retains the reference kernels because its initial fused comparison exceeded the waveform tolerance. The two sectors share displacement histories while retaining their distinct Q-dependent weights. The equations, timestep and precision are unchanged. VT cut cells retain the existing tetrahedral path, CPML remains elastic, and other BKT variants or nonmatching spectra retain the previous kernels.",
          "Inactive sectors are excluded from rate comparisons: BKT2 needs relaxation history only where a0 or a1 is nonzero. Generated records initialize the unused third-mechanism a2/g2 fields. The BKT/BKT2 eligibility scan ignores these fields even in older mesh records. Active, shear-only, bulk-only and inactive elements can share a fused batch; different Q-dependent weights are compatible when the active rates match. BKT3/BKT3F still check their third mechanism. The existing BKT2 Qs >= 1000 and Qk >= 1000 attenuation cutoffs and use_infinite_qk behavior are unchanged.",
          "HERCULES_BKT2_FUSED defaults to 1. Set it to 0 in the launch environment to retain the previous separate BKT2 convolution and force kernels. HERCULES_FORCE_GENERIC_DAMPING=1 takes precedence and selects the generic reference. The setup log reports the selected paths and fallback reason rank counts: fusion disabled, generic override, force or displacement not device-resident, nonuniform active shear or bulk rates, and a shear/bulk rate mismatch. Reason counts may overlap. Ranks with no damping slots are counted separately. Element slots include inactive, acoustic or CPML entries and do not measure attenuating work. Station stress vectors are saved every step when requested, including with decimated station output. Buffer allocations are retained; the optimization reduces memory traffic, not allocated GPU memory.",
          "Compare total solver time with identical inputs and MPI/GPU settings. The fused work appears in Force damping local wait only; the convolution timer covers only the remaining legacy subset. These host wait timers are not isolated GPU kernel timings. Atomic force accumulation may change roundoff, so compare waveforms with a precision-appropriate tolerance. Small-case equivalence or timing does not establish a production Izmir speedup.",
        ],
        codeLabel: "Launch environment — use separate output directories for each run",
        code: "# Default optimized BKT2\nexport HERCULES_BKT2_FUSED=1\n\n# Reference run with the same parameters.in and executable\n# export HERCULES_BKT2_FUSED=0",
      },
      {
        id: "damping-bkt-parametric-q",
        title: "BKT-family example with parametric Q",
        paragraphs: [
          "With use_parametric_q_factor=yes, Hercules computes Qs=C+alpha×Vs_km_s^beta. Qp still comes from the material payload when a profile or 3-D Q model supplies it; use_infinite_qk=yes instead suppresses bulk attenuation by assigning the implementation's high-Qk limit.",
        ],
        codeLabel: "parameters.in — illustrative parametric-Q syntax",
        code: "type_of_damping = bkt3\nsimulation_velocity_profile_freq_hz = 1.0\nuse_parametric_q_factor = yes\nparametric_q_factor_constant = 20.0\nparametric_q_factor_alpha = 50.0\nparametric_q_factor_beta = 1.0\nuse_infinite_qk = no\nthe_threshold_damping = 0.05\ndo_damping_statistics = 1",
        note: "The coefficients merely demonstrate the parser and formula. Fit C, alpha, and beta to an appropriate regional attenuation model before production use.",
      },
      {
        id: "pml",
        title: "Perfectly matched layer rules",
        paragraphs: [
          "abc_type=dashpot selects the legacy Lysmer boundary. abc_type=cpml selects the unsplit convolutional PML on x_min, x_max, y_min, y_max, and z_max; z_min remains the free surface. CPML supports elastic and acoustic elements, while poroelastic and nonlinear candidates inside the PML are intentionally left as linear elastic CPML elements.",
          "Without a positive cpml_thickness_m, the thickness is cpml_layers times the largest non-air element edge touching any lateral or bottom boundary, using an MPI-wide maximum. Air elements (negative Vp) are excluded from both this scan and the global fallback used if no non-air element touches those faces. This is one physical thickness, not an exact cell count everywhere on an adaptive mesh. Material extension and optional CPML partition weighting use the same exclusion rule.",
          "Choose either a positive cpml_thickness_m or a positive cpml_layers count. The default example below uses kappa=1 and complete Hex8 CPML integration. Startup prints the realized thickness, reference wave speed, critical time-step estimate, and element counts; verify them before relying on absorption. With abc_type=cpml, cpml_clamp_outer defaults to yes: outer lateral and bottom nodes are clamped. Set cpml_clamp_outer=no explicitly to remove that clamp while keeping CPML active.",
          "The damping profile is a function of the normal coordinate only. With cpml_profile_slab=auto (opt-in) that coordinate is quantized into octree-aligned slabs as wide as the coarsest non-air PML element (MPI-wide maximum), so every element size inside a slab receives identical d, alpha and kappa and 2:1 size transitions inside the PML no longer sample d half a cell apart; none (default) retains per-element-centre sampling. A positive cpml_reference_vp replaces the automatic maximum PML Vp with one uniform reference speed. After enrollment the solver prints the sampling rule, the reference speed, and the maximum per-cell stretch d_max*h/V (V is Vs for elastic and Vp for acoustic elements) with the class, h and V of the responsible element, and warns when it exceeds 10 based on growth seen in some local operator tests. This is an empirical warning, not a stability bound. Check reflection and long-time stability when changing cpml_reflection_target, PML thickness (cpml_layers only for automatic thickness), or cpml_kappa_max. Increasing kappa does not change this diagnostic itself. These are setup-time changes; kernels and per-step work are unchanged.",
        ],
      },
      {
        id: "pml-partition-memory",
        title: "Completed-mesh partitioning and GPU memory",
        paragraphs: [
          "cpml_partition_weight defaults to 1 (equal leaf counts). With a weight above 1, ordinary partitions are still used throughout progressive construction, 2:1 balancing, building carving and local acoustic-interface matching. CPML weighting is deferred until the completed solver mesh: after unified handoff to the GPU ranks, during Phase B import, or by final re-extraction in a direct or equal-rank unified run. Mesh-only runs remain equally partitioned. The log identifies the stage, rank count, weight and completed-mesh thickness.",
          "Weighting changes ownership, not refinement or corrected material records. Direct re-extraction preserves Q and geostatic data and rebuilds local interface metadata before physics initialization. Existing porous/nonlinear CPML restrictions are unchanged. This avoids weighted intermediate partitions followed by more balancing, but does not establish a general repair of OCTOR anchor-search failures.",
          "GPU memory reports show maximum observed use and minimum free memory before CPML upload, after the main solver allocation, and before the first timestep after communication/source/physics initialization. Each includes the responsible solver/world rank, node/device and element/node/PML counts. The largest planned CPML maps, coefficients and histories (including gravity-surface state) are also reported before upload; this is only a subset of the required storage. Ordinary solver, interface, source, nonlinear/poroelastic and MPI buffers are additional.",
          "These are setup-only snapshots, not a GPU-memory cap or a guarantee of later peak use. Occupancy includes CUDA context and other processes, and later lazy allocations may add storage. Main solver setup allocation failures print requested bytes, CUDA error and free/total memory with rank/device identity; other modules retain their own failure handling. No diagnostics are added inside the timestep loop. Validate the final partition and target-device headroom before a large production run.",
        ],
      },
      {
        id: "pml-complete-operator",
        title: "Complete second-order operator",
        paragraphs: [
          "Each coordinate stretch is S_j(s)=kappa_j+d_j/(s+alpha_j). The full Jacobian J=S_x S_y S_z multiplies inertia; elastic weak stiffness uses J/(S_j S_k), and scalar acoustic stiffness uses J/S_j². This includes tangential stiffness and edge/corner products omitted by the earlier gradient-only correction. Native mode integrates every Hex8 mode and equals the interior stiffness when unstretched. Reduced mode remains rank deficient and is diagnostic only.",
          "Trapezoidal auxiliary equations implement the rational products without divided-pole singularities. Inertia memory is nodal, as is the leading lumped mass. Interface fluxes use J/S_normal with common coefficients sampled at the physical fragment center on both sides, independent of neighbor size or MPI owner. Side and quadrature histories remain independent. This prevents differing cofactors on a nonmatching face; it does not establish long-time stability for arbitrary heterogeneous CPML. The gravity surface uses the full tangential S_x S_y product at each top node. All histories stay on the GPU and use the existing force exchanges.",
          "The default cpml_mpml_cross_fraction is now 0. Nonzero mixing breaks the separable coordinate mapping and can increase reflections substantially; an explicit 0.1 in an older input remains active. Use stiffness_calculation_method=effective for non-BKT CPML runs. CPML rejects checkpoint reads and writes because the legacy format omits absorber histories: use_checkpoint and checkpointing_rate must both be zero.",
          "cpml_acoustic_sponge_fraction is an optional experimental acoustic bulk absorber and defaults to 0 (disabled); positive values enable it. It adds a=gamma*(d_x+d_y+d_z) to acoustic volume inertia for either scalar formulation, independently of acoustic_quality_factor. This option leaves gravity-surface stretching, interface cofactors and interior coefficients unchanged. It adds no static force and no GPU state or exchanges. Nonzero values change matching and require separate accuracy and long-time stability checks.",
          "Native histories require 204 solver scalars per elastic PML element and 60 per acoustic element, plus coefficients. Gravity faces add eight scalars and coupled faces add 24 per side. Verify device memory and timing on the target mesh. Earlier published CPS numbers concern incomplete operators and do not validate arbitrary settings of this formulation.",
        ],
      },
      {
        id: "pml-interface-diagnostics",
        title: "Interface coverage and gravity-surface diagnostics",
        paragraphs: [
          "CPML acoustic-elastic interface: reports global interface-fragment and PML-fragment counts, remote solid-side fragments and ranks with enrolled PML fragments but no host side mask. The assumption line names the active pressure or velocity_potential coupling, or says that coupling is disabled. The cofactor J/S_j is configured at the fragment centre on both sides, with six independent CFS histories per quadrature point per side; point uses one point and face_mass uses four. These are host setup counts, not verification of device histories or kernel execution.",
          "Zero PML fragments are informational, not a missing-correction warning. Water in a lateral PML band and an interface somewhere else do not establish an intersection: a solid inclusion may lie entirely inside the water domain, away from the PML. No correction is required by the centre-sampling rule when all fragment centres are outside. If an intersection was expected, inspect the fragment geometry and realized thickness. A nonzero count of ranks with enrolled PML fragments but no host side mask does produce a CPML WARNING. Routine diagnostics go to stdout and monitor_file; warnings go to stderr and monitor_file. The existing setup reduction grows from three counters to five, with no extra collective or per-step work.",
          "With acoustic_free_surface_gravity=yes, cpml_gravity_in_bottom_band separately counts surface elements inside the geometric bottom band and surface elements outside it with lateral M-PML cross-z damping. Deepen the domain or thin the band until the surface element centres lie outside it, using actual solver coordinates including the topography offset. Water depth alone does not determine clearance. Set cpml_mpml_cross_fraction=0 to remove lateral cross-z damping. The abort/allow policy is unchanged; allow is a diagnostic override. The gravity surface uses tangential S_x*S_y inertia, while the bulk uses S_x*S_y*S_z; the different factors alone do not prove an implementation error. Instability was observed in specific operator-harness configurations, not established for every flagged case. Switching to dashpots changes the absorber; disabling gravity changes the surface to pressure release.",
          "A local conforming-mesh harness stayed bounded with the cofactor filters at its tested profiles. Deliberately removing them produced about 0.45/s growth at 10.85 Hz for a strong profile; that differs from the reported 1.465 Hz Izmir mode and does not identify its cause. The harness omits hanging nodes, octree size transitions, fragment splitting, MPI and material damping. It does not validate arbitrary heterogeneous, poroelastic, viscoelastic or nonlinear coupled cases. The diagnostics change no operators, hydrostatic/geostatic initialization, GPU/MPI exchanges or output recovery; existing poroelastic and nonlinear exclusions inside CPML remain. Host diagnostic checks are not CUDA/MPI propagation validation.",
        ],
      },
      {
        id: "pml-horizon-validation",
        title: "Horizon validation, October 2026",
        paragraphs: [
          "The 100 s homogeneous 10 km versus 20 km comparison preserves mesh spacing, source and physical receivers. With native integration and zero cross-axis mixing, the median/worst 1 Hz displacement errors over the same 25 interior receivers fell from 11.47%/32.99% in the earlier implementation to 1.10%/2.69%. Eight added side/bottom/corner receivers give 1.88%/4.98% over 33 interior receivers. The comparison retains DC and uses no amplitude or time fit. Raw velocity errors remain larger because the source contains unresolved high frequencies.",
          "Coupled flat-water comparisons gave median/worst displacement errors of 0.65%/2.06% with gravity off and 1.68%/3.33% with gravity on. Separate CUDA tests of the production helpers cover P/S incidence through 80 degrees, Rayleigh packets and slow gravity waves. The full elastic solver remained finite through 1,000 s; late velocity RMS decreased at every recorded station. One/four-rank and host/CUDA-aware exchange comparisons remained finite. A matched single-run timing check increased total solver wall time from 14.13 to 17.16 seconds, about 21%; compute and absorber memory costs are higher. Detailed metrics, scope and retained local evidence are in doc/Extended README.md.",
          "These results cover double-precision homogeneous cases, not arbitrary heterogeneous meshes, attenuation combinations or curved free surfaces. These are historical comparisons with their recorded inputs, not validation of every current default combination. The investigation harnesses are local-only and unavailable in a clean checkout. CPML checkpoint/restart remains unsupported.",
        ],
      },
      {
        id: "pml-baseline",
        title: "Default CPML example",
        codeLabel: "parameters.in — current CPML defaults",
        code: "abc_type = cpml\ncpml_thickness_m = 0.0\ncpml_layers = 10\ncpml_reflection_target = 0.001\n# Automatic reference Vp = maximum Vp over PML elements\ncpml_reference_vp = -1.0\ncpml_poly_order = 2\ncpml_kappa_max = 1\n# Automatic alpha in Hz = simulation_wave_max_freq_hz / 2\ncpml_alpha_max_hz = -1.0\ncpml_cfl_safety = 0.9\ncpml_mpml_cross_fraction = 0\ncpml_acoustic_sponge_fraction = 0\ncpml_alpha_floor_fraction = 0.1\ncpml_profile_slab = none\ncpml_partition_weight = 1.0\ncpml_clamp_outer = yes\ncpml_3d_material_extension = yes\ncpml_stiffness_mode = native\ncpml_stiffness_mode_acoustic = native\ncpml_pml_material_damping = none",
        note: "Omitted cpml_alpha_max_hz or a nonpositive value uses simulation_wave_max_freq_hz / 2 in Hz; for 4.36 Hz this is 2.18 Hz. Enter a positive numeric value to override it, not an arithmetic expression. These defaults are not a universal stability or accuracy guarantee; slow gravity packets can require a lower alpha. Explicit settings in existing inputs still win. To reproduce the CPML defaults used before this change, set cpml_acoustic_sponge_fraction=0.3, cpml_reflection_target=0.1 and cpml_kappa_max=10. Automatic alpha and the acoustic interface method are unchanged by this defaults change. Rebuild the executable to adopt the new defaults; existing binaries retain their compiled defaults. A smaller reflection target increases attenuation strength but need not reduce numerical reflection.",
      },
      {
        id: "pml-material-damping",
        title: "CPML with Rayleigh or BKT attenuation",
        paragraphs: [
          "Rayleigh material damping may continue inside elastic PML elements. BKT memory variables are intentionally not applied inside the PML: use elastic CPML around a BKT-attenuating interior. These policies avoid a damping-operator mismatch while keeping the supported PML implementation explicit.",
        ],
        codeLabel: "parameters.in — choose one compatible pair",
        code: "# Rayleigh interior and Rayleigh-damped elastic CPML\ntype_of_damping = rayleigh\nabc_type = cpml\ncpml_pml_material_damping = rayleigh\n\n# OR: BKT-family interior and elastic CPML\n# type_of_damping = bkt2\n# use_parametric_q_factor = no\n# use_infinite_qk = no\n# abc_type = cpml\n# cpml_pml_material_damping = none",
        note: "Do not set cpml_pml_material_damping to a BKT value; the accepted values are none and rayleigh. With Rayleigh+CPML, omitting the setting auto-promotes it to rayleigh, but writing it explicitly makes the physical choice reviewable.",
      },
      {
        id: "pml-heterogeneous-topography",
        title: "CPML with 3-D materials, water, and topography",
        paragraphs: [
          "Keep cpml_3d_material_extension=yes when a finite 3-D material grid or water classification should be extended from the nearest interior column through the absorber. This avoids sampling outside the model footprint and creating a false coastline or basin boundary inside CPML corners.",
          "Material extension clamps in absolute model coordinates. With SQD topography, Hercules restores model Z from depth below the original surface before applying the bottom clamp, then converts back to depth below the surface at the clamped X/Y location. Building surface offsets are restored and removed in the same conversion. This also applies at bottom/lateral corners; a nonzero topographybase_zcoord must not delay the bottom extension.",
          "For example, with the surface at model Z = 2000 m and the bottom PML entrance at Z = 9000 m, a query at Z = 9500 m has an input depth of 7500 m. The extended material is sampled at Z = 9000 m, or depth 7000 m below that surface. At a corner, the final depth uses the surface at the clamped horizontal location.",
          "When topography reaches a lateral PML, taper it to the flat boundary datum before running. script/taper_topography.py provides the repository workflow. The material extension controls properties, while the taper controls geometry; both are needed for a clean heterogeneous PML interface.",
        ],
        codeLabel: "parameters.in — heterogeneous acoustic/elastic CPML additions",
        code: "3D_velocity_model = yes\ninclude_acoustic = yes\ninclude_topography = yes\n\nabc_type = cpml\ncpml_3d_material_extension = yes\ncpml_stiffness_mode = native\ncpml_stiffness_mode_acoustic = native\ncpml_kappa_max = 1.0\ncpml_clamp_outer = yes",
        note: "Changing cpml_3d_material_extension changes sampled PML material properties. Rebuild or recorrect an imported mesh when switching it or adopting the corrected topographic bottom/corner extension; a Phase B import otherwise retains the properties stored in that mesh.",
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
            ["HERCULES_SCALAR_MPI_SCHEDULE", "compact | full", "Scalar field schedules on both backends; default: compact; full restores general schedules"],
            ["HERCULES_SCALAR_MPI_BATCH", "0 | 1", "Combine compatible scalar packets on both backends; default: 1; preserves time levels"],
            ["HERCULES_GPU_MPI_OVERLAP", "0 | 1", "CUDA-aware scalar interior updates and receive-buffer reuse events; default: 1; host staging remains synchronous"],
            ["HERCULES_PORE_PRESSURE_UPDATE", "output | step", "Linear undrained diagnostic pressure; default: output; checkpoint/restart and ineligible physics retain step"],
            ["HERCULES_GPU_MPI_PROFILE", "0 | 1", "Device-helper and scalar-batch timing and packet counters; default: 0; cumulative rank average/min/max"],
            ["HERCULES_NONLINEAR_BACKEND", "cpu | gpu", "Select nonlinear constitutive and force backend; default: gpu"],
            ["HERCULES_NONLINEAR_STATS_INTERVAL", "0 | positive integer", "Sample global ManzariDafalias GPU substep workload; default: 0 (disabled)"],
            ["HERCULES_DISP_SYNC_MODE", "auto | full | boundary", "Control elastic host/device synchronization"],
            ["HERCULES_MESH_IO_MODE", "mpiio | sharded | single", "Select Phase A mesh storage"],
            ["HERCULES_PLANE_IO_MODE", "mpiio | hybrid | legacy", "Select plane interpolation and writing; default: mpiio (dedicated I/O ranks select legacy)"],
            ["HERCULES_PLANE_METADATA", "0 | 1", "Write plane point materials and actual CPML membership once after solver setup; default: 0"],
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
        title: "Complete parameter reference",
        paragraphs: [
          "The categorized pages below document all 249 unique keys currently consumed from parameters.in by psolve and its active feature modules. Each entry states whether it is required, conditional, optional, or a legacy compatibility alias, together with its type, accepted values or default, and scientific or operational meaning.",
        ],
        table: {
          headers: ["Reference page", "Documented entries"],
          rows: parameterReferenceSummary.map((group) => [group.title, String(group.count)]),
        },
        note: "Source-specific controls such as type_of_source, force_x_n, and source_function_type are intentionally not counted here: they are parsed from source.in, not parameters.in.",
      },
      {
        id: "minimal-structure",
        title: "Minimal elastic-profile skeleton",
        paragraphs: [
          "This skeleton contains the unconditional solver controls and a one-layer profile. It still needs a valid sourcefiles/source.in and appropriate scientific values; output and optional-physics blocks may then be added from the reference pages.",
        ],
        codeLabel: "parameters.in",
        code: "# Domain\nregion_origin_latitude_deg   = 0.0\nregion_origin_longitude_deg  = 0.0\nregion_depth_shallow_m       = 0.0\nregion_length_east_m         = 1000.0\nregion_length_north_m        = 1000.0\nregion_depth_deep_m          = 1000.0\nregion_azimuth_leftface_deg  = 0.0\n\n# Time and mesh resolution\nsimulation_end_time_sec        = 10.0\nsimulation_delta_time_sec      = 0.001\nsimulation_wave_max_freq_hz    = 2.0\nsimulation_node_per_wavelength = 8\nsimulation_shear_velocity_min  = 500.0\n\n# Material and damping\ntype_of_damping          = none\nthe_threshold_damping    = 0.05\nthe_threshold_Vp_over_Vs = 3.0\ncvmdb_input_file         = profile\nnumber_profile_layers    = 1\nprofile_layers =\n  0.0  3000.0  1700.0  2300.0  1000.0  1000.0\n\n# Source path; source physics stays in sourcefiles/source.in\nsource_directory = sourcefiles",
      },
    ],
  },
  ...parameterReferencePages,
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
          "Only the inputs enabled by parameters.in are required. Keeping each simulation self-contained makes it easier to move between systems and archive with a publication. The qs.in and qp.in quality-factor files are optional, but they must be supplied together when the three-dimensional material model provides Q values.",
        ],
        codeLabel: "directory tree",
        code: "case/\n├── parameters.in\n├── sourcefiles/\n│   └── source.in\n├── materialfiles/\n│   ├── x.in\n│   ├── y.in\n│   ├── z.in\n│   ├── vs.in\n│   ├── vp.in\n│   ├── rho.in\n│   ├── qs.in\n│   ├── qp.in\n│   └── index.in\n├── topography/\n│   └── topography.in\n└── output/",
        note: "For a 3-D material model, provide qs.in and qp.in as a matched pair with the same number of samples as z.in. Solid-only interpolation excludes water placeholders and renormalizes the remaining weights. Supplied Q is independent of the Vs minimum. If both files are absent, Hercules uses the legacy relations Qs=0.1 Vs and Qp=2 Qs.",
      },
      {
        id: "sources",
        title: "Source input",
        paragraphs: [
          "source_directory in parameters.in points to a directory containing source.in and any auxiliary histories or fault tables required by the selected source model. The source file uses the same assignment and comment syntax as parameters.in.",
          "Hercules currently implements six type_of_source values. The Source models page documents every source.in key, the auxiliary-file layouts, coordinate conventions, and complete examples.",
        ],
        table: {
          headers: ["type_of_source", "Representation", "Additional files"],
          rows: [
            ["point", "Single moment-tensor point source", "None"],
            ["force", "Single Cartesian body force", "forcefunction.in when the time function is discrete"],
            ["base_force", "Distributed compliant-base traction", "baseforcefunction.in"],
            ["plane", "Rectangular finite-fault plane", "slip.in and rake.in"],
            ["planewithkinks", "Vertical finite fault following a kinked surface trace", "slip.in and rake.in"],
            ["srfh", "Standard Rupture Format, Hercules layout", "coords, area, strike, dip, rake, slip, and slip-function files"],
          ],
        },
        note: "Incident DRM plane waves are configured in parameters.in and are not a type_of_source value.",
      },
      {
        id: "basin-configuration",
        title: "Configure a basin model",
        paragraphs: [
          "basin_input_file points to a plain-text grid of basin-bottom depths. Enable it in parameters.in with basin_velocity_model=yes and specify the geographic location of the grid's (0, 0) origin. Relative file paths are resolved from the solver's working directory.",
          "The latitude and longitude below are illustrative; replace them with your basin-grid origin. Keep a valid background cvmdb_input_file or profile for points outside the basin. The example disables the Istanbul and general 3-D replacements to select the basin path explicitly.",
        ],
        codeLabel: "parameters.in — basin settings",
        code: "Istanbul_velocity_model    = no\n3D_velocity_model          = no\nbasin_velocity_model       = yes\nbasin_input_file           = input/basin.in\nbasin_origin_latitude_deg  = 34.0000\nbasin_origin_longitude_deg = -118.5000\nbasin_material_type        = elastic",
      },
      {
        id: "basin-input-format",
        title: "Basin-bottom grid format",
        paragraphs: [
          "The first line is Nx Ny 1, followed by exactly Nx × Ny rows of x y bottom_depth. Use whitespace-separated numbers with no column labels or comments. Nx and Ny must each be at least 2. Set the third header value to 1 because the current lookup uses a single bottom surface.",
          "All three columns are in metres. X is north and Y is east, measured from the basin origin. Bottom depth is positive downward and uses the same depth reference as solver Z; it is not a thickness measured below local topography.",
        ],
        bullets: [
          "Start at (0, 0) and use a positive, uniform spacing with Δx = Δy. The reader infers the common spacing from the first two X values.",
          "Write X fastest, with Y as the outer loop and X as the inner loop. For zero-based indices i and j, data row = j * Nx + i, with x = i * spacing and y = j * spacing.",
          "The example below has 1,000 m spacing and a 500 m bottom depth at its centre. It illustrates the file format and is not a validated simulation case.",
        ],
        codeLabel: "input/basin.in — 3 × 3 grid",
        code: "3 3 1\n0    0    0\n1000 0    0\n2000 0    0\n0    1000 0\n1000 1000 500\n2000 1000 0\n0    2000 0\n1000 2000 0\n2000 2000 0",
        note: "Current maximum-edge limitation: the interpolator accesses the next X and Y grid points without clamping at the final row or column. Extend the grid beyond the largest solver query in each direction, keeping queried basin-relative X < (Nx - 1) * spacing and Y < (Ny - 1) * spacing. A query exactly on either maximum edge can access invalid grid data.",
      },
      {
        id: "basin-material-lookup",
        title: "Basin depth and material lookup",
        paragraphs: [
          "Hercules bilinearly interpolates bottom_depth from the four surrounding grid points. A point within the horizontal grid is inside the basin when z < bottom_depth(x, y). Points on or below the bottom, or outside the horizontal footprint, use the background cvmdb_input_file or profile. A zero bottom depth excludes points at nonnegative solver depths.",
          "The file supplies geometry. With basin_material_type=elastic, material properties come from the hard-coded depth-dependent formulas in quake/forward/basin.cc. With basin_material_type=acoustic_water and include_acoustic=yes, basin_water_vp and basin_water_rho supply water properties; Vs, Qs, and Qp are zero.",
          "The repository example benchmark_inputs/bowl/bowl_basin_250m.in has a 162 × 162 grid, 250 m spacing, and a 162 162 1 header. Its rows follow the same X-fastest ordering. The current reader and lookup are readBasin, interpolateZ, and isInBasin in quake/forward/psolve.cc.",
        ],
      },
      {
        id: "rotated-domain",
        title: "Rotated domains and material grids",
        paragraphs: [
          "With PROJ, region_azimuth_leftface_deg rotates solver +X clockwise from UTM grid north about the region origin. +Y is 90 degrees clockwise from +X; Z remains positive down. The north/east length parameters remain the local X/Y side lengths. This is grid north in the origin's UTM zone, not true north or a rotation about the rectangle's centre.",
          "If 3D_velocity_model_azimuth_deg is omitted, it inherits region_azimuth_leftface_deg, aligning the material axes with the solver axes. Explicit values override this default: set 0 to keep an existing north/east material grid geographically fixed while the solver domain rotates. Prepare x.in/y.in in the selected material frame. The material origin may differ from the domain origin.",
          "The example below is an input-format example, not a validated scientific case. With a 30-degree domain and zero-degree material grid, solver (X,Y)=(1000,0) lies 866.025 m north and 500 m east of the origin. When the material angle is omitted or explicitly set to 30 degrees, matching origins make the same point query material (x,y)=(1000,0).",
        ],
        codeLabel: "parameters.in — add to a complete case",
        code: "region_origin_latitude_deg  = 37.791109\nregion_origin_longitude_deg = 26.383378\nregion_length_north_m       = 90000\nregion_length_east_m        = 90000\nregion_depth_shallow_m      = 0\nregion_depth_deep_m         = 39375\nregion_azimuth_leftface_deg = 30\n\ncvmdb_input_file = profile\nIstanbul_velocity_model = no\nbasin_velocity_model = no\n3D_velocity_model = yes\n3D_velocity_model_directory = input/materialfiles\n3D_velocity_model_origin_latitude_deg  = 37.791109\n3D_velocity_model_origin_longitude_deg = 26.383378\n# Omit the next line for a material grid aligned with the domain.\n3D_velocity_model_azimuth_deg = 0  # explicit override: UTM north/east grid",
        bullets: [
          "Geographic stations, plane origins, point sources, SRFH sources, and planar-fault origins use the rotated mapping. Geographic source strikes are corrected from true north into solver axes. Cartesian source coordinates and body-force components remain solver-local. Plane and DRM-box strikes and incident-wave directions are measured from solver +X.",
          "Topography, groundwater, building coordinates, local DRM geometry, and directional permeability must be prepared in solver axes. These local-coordinate files have no geographic header and are not automatically resampled. Basin grids retain UTM north/east axes under PROJ. Legacy etree CVM and Istanbul readers reject nonzero domain rotation; use a valid profile background, optionally with the general 3-D or basin model.",
          "Material queries rotate into the material axes after CPML extension clamps in solver coordinates. Cached topographic elevations at material columns use the inverse transform. Topography outside its horizontal coverage retains the reader's existing edge clamping. Extend material coverage beyond all intended sample points, including interpolation stencils.",
          "All displacement, velocity, acceleration, stress, strain, and exported XYZ coordinates remain in solver axes. PROJ writes coordinate_frame.txt in station, plane, and DRM output directories, overwriting it at setup. It contains the UTM zone, hemisphere, origin N0/E0, and rotation c/s; binary layouts are unchanged. For station_binary=yes, this metadata is in output_stations_directory even if the binary file is elsewhere.",
          "Rebuild persisted meshes and regenerate source forces after changing either azimuth or origin. Non-PROJ builds retain their corner mapping and solver-aligned material grid when the material angle is omitted or explicitly 0; explicit nonzero overrides require PROJ. Local geometry/interpolation checks do not establish CUDA/MPI benchmark validation.",
        ],
      },
      {
        id: "rotated-output-reader",
        title: "Convert solver output to UTM north/east",
        paragraphs: [
          "Use the sidecar coefficients with any station, plane, or DRM vector. Add N0/E0 only for positions. UTM north/east components differ from true geographic north/east by meridian convergence; obtain latitude/longitude by inverse projection in the recorded UTM zone. For tensors use A T A-transpose, where A has rows (c,-s,0), (s,c,0), (0,0,1).",
        ],
        codeLabel: "Python — read coordinate_frame.txt and convert horizontal coordinates",
        code: "from pathlib import Path\n\nframe = dict(line.split(\" = \", 1) for line in\n             Path(\"outputfiles/stations/coordinate_frame.txt\").read_text().splitlines()\n             if line and not line.startswith(\"#\"))\nc, s = float(frame[\"c\"]), float(frame[\"s\"])\n\ndef vector_to_utm(ux, uy, uz):\n    return c*ux - s*uy, s*ux + c*uy, uz\n\ndef position_to_utm(x, y):\n    return float(frame[\"N0\"]) + c*x - s*y, float(frame[\"E0\"]) + s*x + c*y",
      },
      {
        id: "velocity-model-configuration",
        title: "Configure a three-dimensional velocity model",
        paragraphs: [
          "The general 3-D model replaces the background material inside its horizontal footprint. The origin parameters locate the model's (x=0,y=0) point geographically; x.in and y.in are metre offsets along the material axes set by 3D_velocity_model_azimuth_deg (inherits region_azimuth_leftface_deg when omitted).",
          "Keep the fallback cvmdb_input_file or profile valid because Hercules uses it outside the structured model. Do not enable the Istanbul or basin replacement for the same case unless their priority over the general 3-D model is intentional.",
        ],
        codeLabel: "parameters.in",
        code: "Istanbul_velocity_model                  = no\nbasin_velocity_model                     = no\n3D_velocity_model                         = yes\n3D_velocity_model_directory               = input/materialfiles\n3D_velocity_model_origin_latitude_deg     = 34.0000\n3D_velocity_model_origin_longitude_deg    = -118.5000",
        note: "The model lookup covers cells between horizontal grid points and requires x<x_max and y<y_max. Extend the material grid beyond every mesh location that should use the 3-D model instead of ending it exactly on the desired coverage boundary.",
      },
      {
        id: "velocity-model-files",
        title: "3-D model files and ordering",
        paragraphs: [
          "A convenient way to prepare a model is one row per material sample with columns x, y, z, Vs, Vp, and density. That row-wise CSV is only a preparation format: Hercules reads the separate count-prefixed files listed below. The column order in your preparation table is your choice, but the converted filenames and their meanings are fixed.",
          "Coordinates use the material-grid frame: x and y follow 3D_velocity_model_azimuth_deg (inherits region_azimuth_leftface_deg when omitted), and z is depth positive downward, all in metres. x.in and y.in must each contain at least two strictly increasing, uniformly spaced coordinates. The x and y spacings may differ, but the spacing within either axis must remain constant.",
          "Treat y as the outer loop and x as the inner loop. For zero-based indices i and j, profile_id = j * nx + i, where i is the x index and j is the y index. Therefore write every x profile at y0 before advancing to y1: (x0,y0), (x1,y0), ..., (x0,y1), (x1,y1), ....",
          "Within one x-y profile, write z from shallow to deep with no duplicate depths. Different profiles may use different z values and different numbers of layers. Append the matching Vs, Vp, density, and optional Q values in exactly the same order as z.in.",
        ],
        bullets: [
          "The first value in every file is a count, not model data. It must equal the number of values that follow.",
          "x.in has nx coordinates and y.in has ny coordinates; index.in has nx * ny zero-based profile starts.",
          "z.in, vs.in, vp.in, and rho.in must have the same total sample count. If used, qs.in and qp.in must both exist and have that same count.",
          "Each solid sample (Vs > 1e-6 m/s) needs finite, positive Qs and Qp. Acoustic placeholders may be zero, but must be finite and nonnegative. Q files are checked while reading, even when solid attenuation is disabled. At a mixed horizontal stencil, only solid samples contribute to Q, with renormalized weights: 10% solid Qs=40 plus 90% water Qs=0 gives solid Qs=40. Velocity and density interpolation are unchanged.",
          "Explicit Qs/Qp are not scaled when Vs is raised to simulation_shear_velocity_min. For example, Vs=80 m/s with Qs=16 and Qp=32 keeps those Q values if Vs is raised to 120 m/s. Regenerate persisted meshes after changing these files or to apply corrected Q interpolation; imported coefficient checks cannot recover the original Q inputs.",
          "index.in contains starts only—do not append a final end offset. The end of the last profile is the total sample count from z.in.",
          "Material sample counts and profile start offsets are signed 64-bit integers; axis counts and the number of horizontal profiles remain limited to 2,147,483,647. Rebuild Hercules with this reader for larger sample counts. Existing count-prefixed files remain compatible, and MPI initialization broadcasts large arrays in bounded chunks.",
          "The complete material model is replicated in host memory on every solver rank. Budget approximately 32 bytes per material sample plus 8 bytes per horizontal profile, with another 8 bytes per profile for the SQD surface cache when needed; optional Qs/Qp add 16 bytes per sample. A 3,003 by 6,003 grid with 121 depths needs about 65.3 GiB per solver rank with the SQD cache and without Q files, before mesh and solver allocations. Text parsing runs on solver rank zero; 64-bit indexing does not distribute the material storage.",
        ],
        table: {
          headers: ["File", "Units", "Contents"],
          rows: [
            ["x.in", "m", "Material-grid X coordinates, ascending and uniformly spaced"],
            ["y.in", "m", "Material-grid Y coordinates, ascending and uniformly spaced"],
            ["z.in", "m", "Depth levels, concatenated in index.in profile order"],
            ["vs.in", "m/s", "Shear-wave speed for every z.in entry; zero is permitted for acoustic material"],
            ["vp.in", "m/s", "Compressional-wave speed for every z.in entry"],
            ["rho.in", "tonne/m3", "Density for every z.in entry; Hercules multiplies by 1000 to obtain kg/m3"],
            ["index.in", "zero-based index", "One vertical-profile start for each x-y column; count must equal nx times ny"],
            ["qs.in and qp.in", "dimensionless", "Optional matched quality factors; provide both or neither, with counts equal to z.in"],
          ],
        },
        note: "For profile p, use samples [index[p], index[p+1]) in every concatenated array. For the last profile, use [index[last], total_sample_count). These are zero-based, half-open ranges.",
      },
      {
        id: "velocity-model-example",
        title: "Worked 2 × 2 velocity-model example",
        paragraphs: [
          "Start with the readable table below. It has two x coordinates, two y coordinates, and four vertical profiles. The profiles deliberately have unequal layer counts: 4, 4, 3, and 4 samples. Rows are already sorted with y outermost, x next, and z innermost.",
          "Hercules does not read this CSV directly. Convert its columns into the files in 3D_velocity_model_directory as shown in the next two steps.",
        ],
        codeLabel: "preparation table — hercules-velocity-model.csv",
        code: [
          "x,y,z,vs,vp,rho",
          "0.0,0.0,0,1709.769,4188.062,1.850",
          "0.0,0.0,1000,1894.495,4640.546,1.898",
          "0.0,0.0,3000,2430.156,5952.642,2.020",
          "0.0,0.0,17400,3400.000,8328.265,2.197",
          "500.0,0.0,0,1743.876,4271.606,1.859",
          "500.0,0.0,1000,1920.505,4704.257,1.905",
          "500.0,0.0,5000,2625.958,6432.257,2.060",
          "500.0,0.0,17400,3400.000,8328.265,2.197",
          "0.0,500.0,0,1577.978,3865.241,1.814",
          "0.0,500.0,1000,1840.914,4509.300,1.885",
          "0.0,500.0,6000,2714.210,6648.430,2.077",
          "500.0,500.0,0,1609.918,3943.478,1.823",
          "500.0,500.0,1000,1861.275,4559.174,1.890",
          "500.0,500.0,5000,2712.018,6643.060,2.076",
          "500.0,500.0,17400,3400.000,8328.265,2.197",
        ].join("\n"),
        note: "Density in this preparation table is in tonne/m3, matching rho.in. Hercules multiplies it by 1000 when constructing material density in kg/m3.",
      },
      {
        id: "velocity-model-example-axes",
        title: "Step 1: write x.in, y.in, and index.in",
        paragraphs: [
          "x.in and y.in contain the unique axis values, not one coordinate per CSV row. Their first lines are therefore both 2.",
          "The profile order is (0,0), (500,0), (0,500), (500,500). The first profile starts at sample 0. After its 4 samples the next starts at 4; after another 4 it starts at 8; after the 3-sample third profile the last starts at 11. That produces index.in = [0, 4, 8, 11].",
        ],
        codeLabel: "input/materialfiles — copy each labeled group to its own file",
        code: [
          "# x.in",
          "2",
          "0.0",
          "500.0",
          "",
          "# y.in",
          "2",
          "0.0",
          "500.0",
          "",
          "# index.in",
          "4",
          "0",
          "4",
          "8",
          "11",
        ].join("\n"),
        note: "The # filename labels above separate the examples; do not copy those comment lines into the numeric files.",
      },
      {
        id: "velocity-model-example-ranges",
        title: "How index.in maps profiles to samples",
        paragraphs: [
          "Apply each half-open range to z.in, vs.in, vp.in, rho.in, and to qs.in and qp.in when present. The next start closes the current range; the z.in count of 15 closes the final range.",
        ],
        table: {
          headers: ["profile_id", "Coordinates (x,y)", "Start", "Sample range", "Layer count"],
          rows: [
            ["0", "(0, 0)", "index[0] = 0", "[0, 4)", "4"],
            ["1", "(500, 0)", "index[1] = 4", "[4, 8)", "4"],
            ["2", "(0, 500)", "index[2] = 8", "[8, 11)", "3"],
            ["3", "(500, 500)", "index[3] = 11", "[11, 15)", "4"],
          ],
        },
      },
      {
        id: "velocity-model-example-arrays",
        title: "Step 2: write the concatenated material arrays",
        paragraphs: [
          "There are 15 CSV rows, so each file begins with 15. z.in is the z column in row order, including the restart at zero when a new x-y profile begins. The Vs, Vp, and density arrays must remain aligned sample-for-sample with those depths.",
        ],
        codeLabel: "input/materialfiles — copy each labeled group to its own file",
        code: [
          "# z.in",
          "15",
          "0", "1000", "3000", "17400",
          "0", "1000", "5000", "17400",
          "0", "1000", "6000",
          "0", "1000", "5000", "17400",
          "",
          "# vs.in (m/s)",
          "15",
          "1709.769", "1894.495", "2430.156", "3400.000",
          "1743.876", "1920.505", "2625.958", "3400.000",
          "1577.978", "1840.914", "2714.210",
          "1609.918", "1861.275", "2712.018", "3400.000",
          "",
          "# vp.in (m/s)",
          "15",
          "4188.062", "4640.546", "5952.642", "8328.265",
          "4271.606", "4704.257", "6432.257", "8328.265",
          "3865.241", "4509.300", "6648.430",
          "3943.478", "4559.174", "6643.060", "8328.265",
          "",
          "# rho.in (tonne/m3)",
          "15",
          "1.850", "1.898", "2.020", "2.197",
          "1.859", "1.905", "2.060", "2.197",
          "1.814", "1.885", "2.077",
          "1.823", "1.890", "2.076", "2.197",
        ].join("\n"),
        note: "The # filename labels are not file content. This example omits qs.in and qp.in; when both are absent, Hercules uses the legacy relations Qs=0.1 Vs and Qp=2 Qs. Regenerate or re-correct a persisted mesh after changing any material or quality-factor file.",
      },
      {
        id: "velocity-model-lookup",
        title: "How Hercules samples the model",
        paragraphs: [
          "At each of the four horizontal grid columns surrounding a query point, Hercules selects the deepest listed z that does not exceed the local depth. Values are therefore piecewise constant with depth within each profile; a listed z acts as the top depth of that material interval. Above the first listed depth, Hercules uses the first sample, and below the last listed depth it uses the last sample.",
          "Hercules then bilinearly interpolates Vs, Vp, density, and optional Qs/Qp across those four columns. The general 3-D reader samples the columns at the same physical elevation. Without explicit seabed geometry, type_of_etree selects the vertical reference described below; computation_method independently selects VT or FEM boundary treatment.",
          "In legacy mode, queries outside xmin <= x < xmax or ymin <= y < ymax use the background model. With full or topography disabled, a query Z below the smallest numeric z.in value also uses the background. In sqd, depths have different surface origins, so endpoint extension is applied separately after converting each column depth; there is no global minimum-depth rejection.",
        ],
      },
      {
        id: "velocity-model-topography-depths",
        title: "3-D material depths with full and sqd topography",
        paragraphs: [
          "In legacy mode, with include_topography=yes and type_of_etree=sqd, z.in stores positive-downward layer-top depths below each column's own terrain surface. The caller supplies d = Z - Z_surface(query). The reader uses d_i = d + h_i - h_query, where heights h are positive upward. h_query comes directly from the actual topography grid at the solver query point, including when that grid is finer than, offset from, or rotated relative to the velocity grid.",
          "With type_of_etree=full, z.in stores layer tops in fixed model Z coordinates and every column is queried at Z, with no column terrain shift. The fixed model origin includes the existing mesh-origin convention; existing building offsets are removed before material lookup. With topography disabled, the fixed-coordinate lookup is retained. Changing full to sqd requires converting z.in to the corresponding reference, not just changing the parameter.",
          "For example, topographybase_zcoord=1000 m and terrain heights 0 and 100 m give surface Z values 1000 and 900 m. A horizontal interface at elevation -40 m is at Z=1040 m: encode its layer tops as 40 and 140 m in sqd, or 1040 m in both columns in full. At the midpoint, Z=1030 m is 80 m below the surface; sqd queries depths 30 and 130 m, and full queries Z=1030 m. Both select the material above the interface.",
          "The chosen general 3-D convention is interpolation at equal physical elevation. It does not interpolate corresponding layers at equal burial depth. A terrain-following layer can be encoded through column depths, but a discontinuity between columns may be blended by the horizontal property interpolation.",
          "Rebuild the solver and regenerate material meshes with mesh_input_mode=0 after this correction. Older general 3-D lookups reversed the terrain offset, reconstructed the query surface on the velocity grid, and shifted full-mode columns. Imported element materials retain those old values. The setup monitor now reports the selected 3-D depth convention.",
          "Validation covers local analytical material queries and CPML sampling, not new GPU/MPI wave-propagation benchmarks. The correction selects Vs, Vp, density, and optional Qs/Qp together and therefore affects elastic, viscoelastic, nonlinear elastoplastic and poroelastic initialization, plus acoustic/solid classification and coupled interfaces. Existing VT/FEM and coupled-formulation restrictions still apply.",
        ],
        equations: [String.raw`d_i=d+Z_{\mathrm{surface}}(q)-Z_{\mathrm{surface}}(i)=d+h_i-h_q\quad(\mathrm{sqd}),\qquad Z_i=Z\quad(\mathrm{full}).`],
      },
      {
        id: "explicit-seabed-geometry",
        title: "Explicit seabed geometry",
        paragraphs: [
          "Set seabed_geometry_file to supply a ground/seabed elevation grid. Omitting it preserves the existing Vs-based acoustic convention. Water coverage and depth come from topography minus bed: a gap greater than seabed_height_tolerance_m (default 0.01 m) is water; coincident surfaces are dry. Elevation sign never determines water: dry land can be below zero, and a lake can be above zero. No separate water-level parameter or wet/dry mask is needed.",
          "With include_topography=no, the outer surface is flat at solver Z=0. No topography.in, topographybase_zcoord, maximum_octant_level, computation_method or type_of_etree is required for this mode; these topography controls do not affect it. The seabed may slope, with elevations negative below the flat zero datum. Positive bed elevations beyond tolerance are rejected by default, or lowered to zero with seabed_above_surface_policy=clip; zero is dry ground. The geometry spacing must divide both horizontal domain lengths. CVM mesh offsets must be zero; cvmdb_input_file=profile normally provides this. Wavelength refinement, optional acoustic interface refinement and 2:1 balancing remain active, but there is no forced topography surface-refinement band. This does not establish sufficient water-column or gravity-wave resolution: check mesh convergence.",
          "With include_topography=yes, topography.in remains the outer air/material boundary: water-surface elevation on wet cells, ground elevation on dry cells. The new file shares its domain, origin, spacing and node ordering: solver X north, Y east, east varying fastest. Both files use the same vertical datum and elevations positive upward. Bed elevations above topography by more than seabed_height_tolerance_m abort setup under seabed_above_surface_policy=error (default); nodal differences within tolerance are snapped to topography under either policy. Wet grid cells require flat outer surfaces below topographybase_zcoord. Distinct water bodies may have different elevations, separated by dry cells; sloping water surfaces are rejected because the acoustic free-surface operator requires horizontal faces.",
          "Set seabed_height_tolerance_m in parameters.in to adjust the shared geometry tolerance: default 0.01 m (1 cm), finite and >= 0; requires seabed_geometry_file when supplied. Gaps <= tolerance are dry, nodal differences within tolerance snap the bed to the outer surface, and wet-cell maximum-minus-minimum surface elevation must be <= tolerance. All ranks use it in interpolated queries as well as loading. Zero means exact comparisons; 1e-8 restores the former tolerance. Larger values can remove real shallow water and change material reference depths, porous eligibility and hydrostatic/geostatic loads. It does not flatten topography or change the separate mesh-face alignment check (one quarter of an octree tick); horizontal acoustic free-surface faces remain required. Rebuild meshes and sources after changing it.",
          "Opt in with seabed_above_surface_policy=clip in parameters.in to lower above-surface bed elevations to the outer surface inside Hercules. This requires seabed_geometry_file; error is the default and both values are case-insensitive. Rank zero reports the lowered-node count (including positive excesses within tolerance) and maximum lowering in metres, then broadcasts the corrected grid before material and interface setup. Original input files are unchanged. Coincident nodes are dry; valid water is retained except for the existing tolerance snapping. Clipping does not make sloping water surfaces valid. model_z keeps fixed layer coordinates, surface keeps the outer-surface datum, and seabed follows the corrected bed. Porous ground/groundwater depths and hydrostatic/geostatic loads use the corrected geometry. Rebuild meshes and sources; this does not repair existing mesh archives or checkpoints and adds no time-stepping work. Local host checks compare against manually corrected geometry and exercise material/Q, phase, depth-reference and initialization paths; full CUDA/MPI propagation validation of clipping remains outstanding.",
          "The geometry file contains one positive spacing value (matching topography.in when enabled), followed by Nnorth × Neast finite ground/bed elevations (metres, positive upward), ordered east-fast. Dimensions are inferred from the domain and spacing. There are no count prefixes, coordinates, masks, comments or extra tokens. One number per line is valid; spaces and newlines are equivalent. Both surfaces are bilinear, and water coverage is inferred at the query point. Outside queries extend the nearest grid edge.",
          "For a 1000 m × 1000 m domain and 500 m geometry spacing, the file below has nine bed elevations. The flat-ocean example needs no topography.in. To use include_topography=yes instead, supply its required controls and nine zero heights in topography.in. Water occupies elevations from −100 m up to 0 m; the bed itself belongs to solid material. Water Vp and density below are illustrative, not a site calibration.",
        ],
        codeLabel: "parameters.in excerpt and geometry file (save separately)",
        code: "# parameters.in excerpt\ninclude_topography = no\n3D_velocity_model = yes\ninclude_acoustic = yes\nseabed_geometry_file = input/seabed.in\n# Optional: clip bed/surface crossings in memory (default: error)\nseabed_above_surface_policy = clip\nseabed_height_tolerance_m = 0.01\n3D_velocity_model_depth_reference = seabed\nseabed_water_vp_m_s = 1500\nseabed_water_density_kg_m3 = 1025\n\n# input/seabed.in — omit this label from the file\n500\n-100\n-100\n-100\n-100\n-100\n-100\n-100\n-100\n-100",
        bullets: [
          "Provide solid-only 3D profiles: every vs.in sample must exceed 1e-6 m/s. Remove water entries consistently from z.in, vs.in, vp.in, rho.in, optional Q files and index.in. Water is supplied by geometry before any solid interpolation, with Vs=0 and the configured Vp and SI density. Acoustic damping still uses acoustic_quality_factor; solid Q never leaks into water.",
          "Explicitly select 3D_velocity_model_depth_reference=seabed, surface or model_z. seabed measures layer-top depth below each column's bed (ground on dry cells); surface measures below the outer boundary; model_z uses fixed solver Z. All are positive downward and work with either full or sqd when topography is enabled. With topography off, type_of_etree is unused: model_z and surface coincide at the zero datum, while seabed remains bed-relative. For bed elevation −100 m, water level 0 and topographybase_zcoord=1000, a solid layer 20 m below bed has z.in=20, 120 or 1120 respectively.",
          "Columns are still sampled at equal physical elevation. Vertical endpoint extension uses the first/last solid sample; horizontal interpolation uses only solid profiles. Solid queries outside the half-open material footprint abort instead of falling back to a background with an ambiguous reference. Extend the material grid beyond the largest query coordinate. CPML clamps in solver coordinates before water classification and solid sampling.",
          "Porous sediment depth limits are measured below the bed in both full and sqd; submerged hydrostatic pressure is referenced to the outer water surface (Z=0 with topography off). groundwaterlevel_file remains ground-relative on dry cells and is ignored beneath explicit water. Nonlinear approximate geostatic setup adds the exact water load and quadratures the solid density column before the existing Biot hydrostatic correction. Nonlinear elastoplastic/plastic and hysteretic calibrations, geostatic controls, tolerances, substeps, damping and material_properties_list remain case-specific.",
          "Acoustic free-surface selection uses the flat Z=0 surface with topography off, or local topography at each wet element with it on, including separate flat lakes at different elevations. A complete horizontal mesh face must coincide with that surface; crossing and submerged faces are excluded. This feature does not introduce cut acoustic free surfaces. The optional bathymetry face filter takes bed depths from the explicit grid, retaining its existing bottom-face-only coupling policy. Default all-face coupling and existing acoustic–elastic, poroelastic–acoustic and poroelastic–elastic formulations remain unchanged.",
          "Rebuild the solver and regenerate mesh_input_mode=0 after changing either elevation grid, depth reference or water properties; imported element properties are not repaired by new interface queries. Geometry is read once by rank zero and broadcast to solver ranks. No new per-step kernels, transfers or MPI messages are introduced. Station/source depth coordinates still use the outer surface (flat Z=0 with topography off), so a seabed station in the example has surface-relative depth 100 m.",
          "This mode requires acoustics and the general 3D model; topography is optional for a flat zero-datum outer surface; buildings, Istanbul, basin and homogeneous-halfspace overrides are rejected. Sloping or changing water surfaces are unsupported. Existing FEM/VT and coupled-formulation restrictions still apply. The 27-point element classifier retains solid priority for mixed water/solid cells: resolve thin water layers and compare mesh convergence; explicit input geometry is not a conforming seabed mesh.",
        ],
        note: "A 2026-10-18 Horizon benchmark at solver commit 8c9a514 compared explicit flat-seabed/topography-off acoustics with independent CPS Green functions: two nodes/eight GPUs, CUDA-aware MPI, velocity potential, field-specific constraints, 110 s and fixed 0.0025 s dt. Coarse/fine meshes had 22,282,240/178,257,920 elements, with the same 1 Hz source and CPML controls. Refined off-axis median raw waveform L2 errors over 10–100 s were 2.89% displacement, 6.64% velocity and 9.80% pressure, without output filtering or fitted gains/shifts. Some refined receiver peaks exceed the 5% target; the two resolutions do not establish mesh independence. All 30 receivers and source history were checked. This wave benchmark covers a horizontal seabed, gravity off and lossless materials; it does not cover gravity waves or host-MPI exchange. Local extracted checks also cover flat initialization, input rejection, all depth references, phase separation, sloping-bed normals, surface assembly, refinement, porous pressure/traction and no-flow conservation, open-pore quadrature, constraint/station maps, solid Q isolation, nonlinear dispatch and hydrostatic/approximate geostatic stress. These setup checks are not a CUDA/MPI propagation benchmark; GPU scheduling and MPI transport are mocked. Nonlinear porous skeletons retain explicit pressure coupling, and zero-permeability gravity-loading ramps are rejected. Submerged ramped equilibrium, sloping-bed propagation and porous, viscoelastic and nonlinear plastic/hysteretic standalone/coupled propagation still need separate supercomputer validation. The local investigation harness and results are ignored and unavailable in a clean checkout.",
      },
      {
        id: "velocity-model-converter",
        title: "Convert a CSV model automatically",
        paragraphs: [
          "The repository includes script/convert_3d_velocity_model.py, a complete converter based on the script pattern in the original Hercules documentation. It uses only the Python standard library, so Pandas and NumPy are not required.",
          "The CSV header must contain x,y,z,vs,vp,rho. It may also contain both qs and qp; providing only one is an error. The converter sorts axis coordinates and profile depths, writes profiles in the required Y-outer/X-inner order, and rejects incomplete rectangular grids, nonuniform horizontal spacing, duplicate depths, non-finite numbers, and invalid material values.",
        ],
        codeLabel: "shell — run from the Hercules_GPU repository",
        code: "python3 script/convert_3d_velocity_model.py \\\n  hercules-velocity-model.csv \\\n  --output-dir input/materialfiles\n\n# Expected summary for the worked example:\n# wrote 2 x 2 profiles and 15 samples to input/materialfiles",
        bullets: [
          "Always writes x.in, y.in, index.in, z.in, vs.in, vp.in, and rho.in.",
          "Writes qs.in and qp.in only when both optional CSV columns are present.",
          "Creates the output directory when necessary and replaces files of the same names.",
          "Run the structural validator in the next section on the generated directory before mesh generation.",
        ],
      },
      {
        id: "velocity-model-converter-script",
        title: "CSV converter script",
        paragraphs: [
          "The checked-in script is reproduced below so the format and conversion logic remain visible with the example. Save it as convert_3d_velocity_model.py if you are using the documentation independently of the repository.",
        ],
        codeLabel: "Python 3 — script/convert_3d_velocity_model.py",
        code: String.raw`#!/usr/bin/env python3
"""Convert a row-wise CSV model to Hercules count-prefixed 3-D model files."""

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

REQUIRED_COLUMNS = ("x", "y", "z", "vs", "vp", "rho")
QUALITY_COLUMNS = ("qs", "qp")

def require_uniform_axis(values, name):
    if len(values) < 2:
        raise ValueError(f"{name}: need at least two unique coordinates")
    spacing = values[1] - values[0]
    for left, right in zip(values, values[1:]):
        if not math.isclose(right-left, spacing, rel_tol=1e-9, abs_tol=1e-8):
            raise ValueError(f"{name}: coordinates are not uniformly spaced")

def write_counted(path, values):
    with path.open("w", encoding="utf-8") as output:
        output.write(f"{len(values)}\n")
        for value in values:
            output.write(f"{value if isinstance(value, int) else format(value, '.17g')}\n")

def read_csv_rows(path):
    with path.open(newline="", encoding="utf-8-sig") as source:
        reader = csv.DictReader(source)
        if reader.fieldnames is None:
            raise ValueError("CSV is missing a header row")
        fieldnames = [name.strip().lower() for name in reader.fieldnames]
        reader.fieldnames = fieldnames
        missing = [name for name in REQUIRED_COLUMNS if name not in fieldnames]
        if missing:
            raise ValueError(f"CSV is missing columns: {', '.join(missing)}")
        has_qs, has_qp = "qs" in fieldnames, "qp" in fieldnames
        if has_qs != has_qp:
            raise ValueError("CSV must contain both qs and qp, or neither")
        columns = REQUIRED_COLUMNS + (QUALITY_COLUMNS if has_qs else ())
        rows = []
        for line_number, raw in enumerate(reader, start=2):
            try:
                row = {name: float(raw[name]) for name in columns}
            except (KeyError, TypeError, ValueError) as error:
                raise ValueError(f"line {line_number}: fields must be numeric") from error
            if not all(math.isfinite(value) for value in row.values()):
                raise ValueError(f"line {line_number}: values must be finite")
            if row["z"] < 0 or row["vs"] < 0 or row["vp"] <= 0 or row["rho"] <= 0:
                raise ValueError(f"line {line_number}: invalid z, Vs, Vp, or density")
            if has_qs and (row["qs"] <= 0 or row["qp"] <= 0):
                raise ValueError(f"line {line_number}: Qs and Qp must be > 0")
            rows.append(row)
    if not rows:
        raise ValueError("CSV contains no model rows")
    return rows, columns

def convert(csv_path, output_dir):
    rows, columns = read_csv_rows(csv_path)
    x_values = sorted({row["x"] for row in rows})
    y_values = sorted({row["y"] for row in rows})
    require_uniform_axis(x_values, "x")
    require_uniform_axis(y_values, "y")

    profiles = defaultdict(list)
    for row in rows:
        profiles[(row["y"], row["x"])].append(row)
    expected = {(y, x) for y in y_values for x in x_values}
    if expected.difference(profiles):
        raise ValueError("CSV does not contain a complete rectangular x-y grid")

    starts, ordered_rows = [], []
    for y_value in y_values:          # Y is the outer loop.
        for x_value in x_values:      # X is the inner loop.
            profile = sorted(profiles[(y_value, x_value)], key=lambda row: row["z"])
            depths = [row["z"] for row in profile]
            if any(b <= a for a, b in zip(depths, depths[1:])):
                raise ValueError(f"duplicate depth at x={x_value:g}, y={y_value:g}")
            starts.append(len(ordered_rows))
            ordered_rows.extend(profile)

    output_dir.mkdir(parents=True, exist_ok=True)
    write_counted(output_dir / "x.in", x_values)
    write_counted(output_dir / "y.in", y_values)
    write_counted(output_dir / "index.in", starts)
    for name in columns:
        if name not in ("x", "y"):
            write_counted(output_dir / f"{name}.in",
                          [row[name] for row in ordered_rows])
    print(f"wrote {len(x_values)} x {len(y_values)} profiles and "
          f"{len(ordered_rows)} samples to {output_dir}")

def main():
    parser = argparse.ArgumentParser(
        description="Convert CSV to Hercules 3-D velocity-model .in files"
    )
    parser.add_argument("csv_file", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    convert(args.csv_file, args.output_dir)

if __name__ == "__main__":
    main()`,
        note: "The displayed version is intentionally compact; the checked-in script reports the coordinates of missing profiles and performs additional header validation. Both produce the same Hercules file ordering.",
      },
      {
        id: "velocity-model-validation",
        title: "Validate 3-D model files with Python",
        paragraphs: [
          "Run this inexpensive structural check before launching mesh generation. It verifies counts, horizontal regularity, profile starts, vertical ordering, material ranges, and the optional Q-file pair.",
          "The standalone script/validate_3d_velocity_model.py accepts command-line arguments; run the command shown below without its leading comment marker. script/README.md lists the utilities and their dependencies.",
        ],
        codeLabel: "Python — script/validate_3d_velocity_model.py",
        code: "# Runnable version from the repository root:\n# python3 script/validate_3d_velocity_model.py input/materialfiles\n\nfrom pathlib import Path\nimport numpy as np\n\nroot = Path(\"input/materialfiles\")\n\ndef counted(name, dtype=float):\n    raw = np.loadtxt(root / name, dtype=dtype, ndmin=1)\n    n, values = int(raw[0]), raw[1:]\n    if values.size != n:\n        raise ValueError(f\"{name}: header says {n}, found {values.size}\")\n    return values\n\nx = counted(\"x.in\")\ny = counted(\"y.in\")\nz = counted(\"z.in\")\nstarts = counted(\"index.in\", int).astype(int)\nvs, vp, rho = (counted(name) for name in (\"vs.in\", \"vp.in\", \"rho.in\"))\n\nfor name, axis in ((\"x.in\", x), (\"y.in\", y)):\n    if axis.size < 2 or np.any(np.diff(axis) <= 0):\n        raise ValueError(f\"{name}: need at least two increasing values\")\n    if not np.allclose(np.diff(axis), np.diff(axis)[0]):\n        raise ValueError(f\"{name}: spacing is not uniform\")\nif starts.size != x.size * y.size or starts[0] != 0 or np.any(np.diff(starts) <= 0):\n    raise ValueError(\"index.in does not define one ordered profile per x-y column\")\nif not (z.size == vs.size == vp.size == rho.size):\n    raise ValueError(\"z, Vs, Vp, and density counts differ\")\nfor column, (begin, end) in enumerate(zip(starts, np.r_[starts[1:], z.size])):\n    if begin >= end or np.any(np.diff(z[begin:end]) <= 0):\n        raise ValueError(f\"vertical profile {column} is empty or not increasing\")\nif np.any(vs < 0) or np.any(vp <= 0) or np.any(rho <= 0):\n    raise ValueError(\"invalid Vs, Vp, or density value\")\nqs_exists, qp_exists = (root / \"qs.in\").exists(), (root / \"qp.in\").exists()\nif qs_exists != qp_exists:\n    raise ValueError(\"provide both qs.in and qp.in, or neither\")\nif qs_exists and (counted(\"qs.in\").size != z.size or counted(\"qp.in\").size != z.size):\n    raise ValueError(\"Q-file counts must match z.in\")\nprint(f\"valid model: {x.size} x {y.size} columns, {z.size} material samples\")",
      },
      {
        id: "topography-configuration",
        title: "Configure topography",
        paragraphs: [
          "topography_directory contains one required file, topography.in. The topographic spacing must divide both region_length_north_m and region_length_east_m exactly. maximum_octant_level and topographybase_zcoord are scientific mesh choices, not universal defaults.",
          "The block below is illustrative for a 40 km square domain with relief below 800 m. Recompute the refinement level from the material-model minimum shear velocity, simulation_shear_velocity_min, simulation_wave_max_freq_hz, and simulation_node_per_wavelength; choose a positive base coordinate that safely contains the full elevation range.",
        ],
        codeLabel: "parameters.in — illustrative values",
        code: "include_topography             = yes\ntopography_directory           = input/topography\ncomputation_method             = vt\nmaximum_octant_level           = 9\ntopographybase_zcoord          = 1250.0\ntype_of_etree                  = full\nmin_enclosed_volume_perc       = 0.1\nconsider_nonlinear_topography  = no\ntopography_bkt_gpu_stabilize   = auto",
        note: "computation_method=vt uses the volume-tetrahedra topographic treatment; fem retains whole-cube staircase behavior. type_of_etree=full uses fixed model Z; sqd uses surface-relative material depths. For the general 3-D files, follow the depth-reference example above. Nonlinear topography additionally requires include_nonlinear_analysis=yes.",
      },
      {
        id: "topography-wavelength-refinement",
        title: "Opt in to wavelength-based surface sizing",
        paragraphs: [
          "topography_surface_refinement=maximum_level is the default and preserves forced refinement near the surface to maximum_octant_level. Set topography_surface_refinement=wavelength to skip that target and its legacy surface-size consistency check. The option is case-insensitive and is read only with include_topography=yes. Startup output reports the selected mode.",
          "The ordinary wavelength rule then sizes solid cells from Vs, simulation_wave_max_freq_hz and simulation_node_per_wavelength, with the existing shear-velocity floor. Water uses Vp and acoustic_mesh_refinement_factor. Air/material classification, other enabled building/DRM/interface refinement and octree 2:1 balancing remain in place; wavelength mode therefore does not promise a uniform mesh or an exact element count.",
          "Keep maximum_octant_level in parameters.in: it is still required for topographic geometry classification/proximity bookkeeping. It is not a global cap on wavelength refinement. Changing this option requires a fresh mesh; mesh_input_mode=0 regenerates it, whereas importing a stored mesh retains its sizing.",
          "This trades forced terrain resolution for fewer surface elements where wavelengths allow them. FEM retains its whole-cube staircase surface; VT samples cut volumes at the resulting element size. Neither guarantees accurate terrain representation from wavelength resolution alone. Check the represented surface, receiver ownership, timestep limits and waveform convergence before using a coarser production mesh.",
          "The shared meshing decision applies to elastic, viscoelastic, nonlinear elastoplastic and poroelastic solids, and to supported acoustic–elastic, porous–water and porous–elastic combinations. Existing formulation restrictions remain: acoustic field-specific constraints require FEM when topography is enabled; local-matching rejects topography; VT retains its existing pore-pressure constraint limitations. Acoustic free surfaces still require the supported flat, complete-face geometry. Hydrostatic/geostatic initialization and constitutive/time-integration paths are unchanged but operate on the new mesh.",
        ],
        codeLabel: "parameters.in — opt-in addition to a valid topography case",
        code: "include_topography             = yes\ntopography_surface_refinement = wavelength\nmesh_input_mode               = 0\n# Retain maximum_octant_level and the remaining topography/material inputs.",
        note: "Local regression checks cover parsing, rank-broadcast serialization and the production refinement/classification functions. This is not CUDA/MPI propagation, coupled-physics accuracy or hydrostatic/geostatic equilibrium validation of a coarser mesh.",
      },
      {
        id: "topography-surface-consistency",
        title: "Surface geometry and receiver sampling",
        paragraphs: [
          "Hercules uses one bilinear topography grid for mesh air/solid classification, material depths, VT volume sampling, sources, stations, and planes. If the interpolated height is h, solver Z = topographybase_zcoord - h. A station depth d is placed at Z + d; d=0 samples the local surface. region_depth_deep_m controls the domain extent, not this surface datum. FEM still uses whole Hex8 staircase elements; sharing the height field does not turn it into a fitted surface mesh.",
          "Earlier builds used a triangular surface for the mesh air/solid test while receivers and material queries used bilinear heights. Nonplanar grid cells could therefore put surface receivers inside air elements, yielding zero traces or interpolation from only a few active nodes. A finer mesh can expose this discrepancy even when a coarser run appeared usable. This predates the shifted acoustic free-surface correction.",
          "Regenerate topographic meshes made with the earlier classifier using mesh_input_mode=0, then rerun the simulation. Importing an old mesh retains its stored air/solid labels; rebuilding only the executable cannot repair those labels or recover missing output samples. Check the new station_mediums.txt for positive Vp and density at land stations. For plane ownership checks, generate fresh metadata with HERCULES_PLANE_METADATA=1; metadata from an earlier mesh cannot establish current ownership.",
          "At the last grid row or column, interpolation uses that edge's heights. Material queries outside the grid use the nearest edge. The mesh proximity test divides the vertical offset by the local bilinear normal length; its sign classifies air versus material, while its magnitude is a tangent-normal offset rather than an exact shortest distance to a curved patch.",
          "The 2026-10-06 local regression reproduced an all-zero station in a 1 km real-terrain crop with a synthetic elastic source. Regenerating the mesh with the shared surface restored that station and all 431 missing histories in a 1,521-point surface plane. Raw plane histories agreed within 2.3e-13 relative L2 across one/eight GPU ranks and host/CUDA-aware MPI. An 8 s flat-water control with gravity, CPML and field-specific constraints agreed with the previous build within 3.7e-12 relative L2 in plane output. These local-only fixtures are not distributed with the repository. They test geometry and output consistency, not continuum accuracy or full-model long-time stability.",
        ],
      },
      {
        id: "topography-maximum-level",
        title: "Choose maximum_octant_level",
        paragraphs: [
          String.raw`\(\mathtt{maximum\_octant\_level}\) is required only when \(\mathtt{include\_topography=yes}\). Let \(V_{S,\min}^{\mathrm{model}}\) denote the smallest shear-wave velocity in the active solid material model; this is derived from the selected material inputs, not read from a separate parameters.in key.`,
          String.raw`For the default topography_surface_refinement=maximum_level, choose the smallest integer \(N\) that satisfies the wavelength inequality below as a baseline for the forced surface target. The accepted range is 0 through 30. The subsequent material wavelength rule can refine further; this parameter is not a global refinement cap. In wavelength mode, the forced target and its legacy surface-size consistency check are skipped.`,
        ],
        equations: [
          String.raw`\frac{
\max\!\left(
\mathtt{region\_length\_north\_m},\,
\mathtt{region\_length\_east\_m},\,
\mathtt{region\_depth\_deep\_m}-\mathtt{region\_depth\_shallow\_m}
\right)
}{2^N}
\;\le\;
\frac{
\max\!\left(V_{S,\min}^{\mathrm{model}},\mathtt{simulation\_shear\_velocity\_min}\right)
}{
\mathtt{simulation\_wave\_max\_freq\_hz}\,
\mathtt{simulation\_node\_per\_wavelength}
}`,
          String.raw`\begin{aligned}
\frac{\max(500,400)}{2\times 8} &= 31.25\ \mathrm{m}, \\
\frac{1{,}000{,}000}{2^{14}} &= 61.035\ \mathrm{m} > 31.25\ \mathrm{m}, \\
\frac{1{,}000{,}000}{2^{15}} &= 30.518\ \mathrm{m} \le 31.25\ \mathrm{m}.
\end{aligned}`,
        ],
        afterEquationsParagraphs: [
          String.raw`For the worked example, the maximum domain extent is \(1{,}000{,}000\,\mathrm{m}\), \(V_{S,\min}^{\mathrm{model}}=500\,\mathrm{m/s}\), \(\mathtt{simulation\_shear\_velocity\_min}=400\,\mathrm{m/s}\), \(\mathtt{simulation\_wave\_max\_freq\_hz}=2\,\mathrm{Hz}\), and \(\mathtt{simulation\_node\_per\_wavelength}=8\). Level 14 gives \(61.035\,\mathrm{m}\) and is too coarse; level 15 gives \(30.518\,\mathrm{m}\), so \(\mathtt{maximum\_octant\_level}=15\) is the smallest acceptable value.`,
        ],
        note: "This is a shear-wave baseline for selecting the default forced surface target. Terrain geometry, acoustic refinement and other enabled feature constraints may require finer elements. With wavelength mode, assess geometry accuracy separately from wave resolution.",
      },
      {
        id: "topography-format",
        title: "topography.in format and ordering",
        paragraphs: [
          "The first value is the common x-y grid spacing in metres. It is followed by one elevation in metres for every inclusive grid point, so the required count is (region_length_north_m / spacing + 1) times (region_length_east_m / spacing + 1).",
          "The current solver lookup uses flat_index = x_index * n_east + y_index. Begin at the southwest corner; for each fixed x position moving south to north, write all y positions from west to east before advancing x. This explicit x-major order prevents an asymmetric surface from being transposed.",
          "For a 1000 m by 1000 m domain with 500 m spacing, the nine elevations below represent z(0,0), z(0,500), z(0,1000), z(500,0), and so on.",
        ],
        codeLabel: "topography.in — 3 × 3 grid",
        code: "500.0\n0.0\n10.0\n20.0\n50.0\n60.0\n70.0\n100.0\n110.0\n120.0",
        note: "topography.in is parsed as whitespace-separated numbers. Do not include explanatory labels in the actual file, and ensure the number of elevations is exact.",
      },
      {
        id: "topography-generation",
        title: "Generate a regular topography grid with Python",
        paragraphs: [
          "This example generates a Gaussian hill using the same x-major flattening implemented by point_elevation. Replace the analytical expression with resampled DEM elevations while preserving the grid, order, and units.",
          "The standalone script/generate_topography.py accepts command-line arguments; run the command shown below without its leading comment marker. script/README.md lists the utilities and their dependencies.",
        ],
        codeLabel: "Python — script/generate_topography.py",
        code: "# Runnable version from the repository root:\n# python3 script/generate_topography.py --lx 40000 --ly 40000 --spacing 100 --height 800 --sigma 3000 --out input/topography/topography.in\n\nimport numpy as np\n\nlength_north = 40_000.0\nlength_east = 40_000.0\nspacing = 100.0\nheight = 800.0\nsigma = 3_000.0\ncenter_x = center_y = 20_000.0\n\nnx = int(round(length_north / spacing)) + 1\nny = int(round(length_east / spacing)) + 1\nif not np.isclose((nx - 1) * spacing, length_north) or not np.isclose((ny - 1) * spacing, length_east):\n    raise ValueError(\"spacing must divide both domain lengths exactly\")\nx = np.arange(nx) * spacing\ny = np.arange(ny) * spacing\nX, Y = np.meshgrid(x, y, indexing=\"ij\")  # array shape: (north, east)\nelevation = height * np.exp(-((X-center_x)**2 + (Y-center_y)**2) / (2*sigma**2))\n\nwith open(\"topography.in\", \"w\") as output:\n    output.write(f\"{spacing:g}\\n\")\n    np.savetxt(output, elevation.reshape(-1, order=\"C\"), fmt=\"%.6f\")\nprint(f\"wrote {nx} x {ny} = {elevation.size} elevations\")",
      },
      {
        id: "topography-validation",
        title: "Validate and inspect topography.in",
        paragraphs: [
          "The check below rejects a spacing/count mismatch, non-finite elevations, and an invalid base coordinate. It also reports perimeter relief, which should be essentially flat when CPML is used for a production model.",
          "The standalone script/validate_topography.py accepts command-line arguments; run the command shown below without its leading comment marker. script/README.md lists the utilities and their dependencies.",
        ],
        codeLabel: "Python — script/validate_topography.py",
        code: "# Runnable version from the repository root:\n# python3 script/validate_topography.py input/topography/topography.in --lx 40000 --ly 40000 --base 1250\n\nimport numpy as np\n\nlength_north = 40_000.0\nlength_east = 40_000.0\ntopography_base = 1_250.0\nraw = np.loadtxt(\"input/topography/topography.in\")\nspacing, values = raw[0], raw[1:]\nnx = int(round(length_north / spacing)) + 1\nny = int(round(length_east / spacing)) + 1\nif not np.isclose((nx-1)*spacing, length_north) or not np.isclose((ny-1)*spacing, length_east):\n    raise ValueError(\"topographic spacing does not divide the domain\")\nif values.size != nx * ny:\n    raise ValueError(f\"expected {nx*ny} elevations, found {values.size}\")\nif not np.isfinite(values).all():\n    raise ValueError(\"topography contains NaN or infinite values\")\nelevation = values.reshape(nx, ny)  # [north index, east index]\nif topography_base <= elevation.max():\n    raise ValueError(\"topographybase_zcoord must exceed the highest elevation\")\nperimeter = np.r_[elevation[0,:], elevation[-1,:], elevation[:,0], elevation[:,-1]]\nprint(f\"grid={nx}x{ny}, range={elevation.min():.3f}..{elevation.max():.3f} m\")\nprint(f\"perimeter range={np.ptp(perimeter):.6g} m\")",
        note: "For CPML runs, taper real relief smoothly to one flat reference within an interior buffer that becomes fully flat before the PML inner edge. The solver's flatten-PML route is a stability backstop, not a substitute for preprocessing the surface.",
      },
      {
        id: "topography-cpml",
        title: "Why topography must be tapered before C-PML",
        paragraphs: [
          "C-PML is a regular, axis-aligned absorbing layer with a flat top. With computation_method=vt, Hercules can force topography-intersecting lateral PML elements back onto a flat cube-and-CPML path. This flatten-PML route prevents an unstable element treatment, but it does not remove the artificial scattering produced when substantial relief is cut off abruptly at the PML inner edge.",
          "For a production model, begin the transition inside the physical domain. Leave the deep-interior elevations unchanged, reduce their deviation from a chosen reference elevation through a smooth raised-cosine ramp, and provide an exactly flat outer collar that covers the complete PML thickness with extra margin.",
          "The repository helper applies z_out = z_ref + factor × (z_in - z_ref). The factor is one in the interior, follows a C1 raised-cosine through the taper ramp, and is exactly zero throughout the outer flat collar. The helper now uses the solver's X-outer/Y-inner topography order, including rectangular grids; earlier versions used transposed point indices.",
        ],
        note: "Surface tapering and cpml_3d_material_extension solve different problems. The taper flattens the PML's top surface; material extension controls the 3-D material or water mask inside the PML volume.",
      },
      {
        id: "topography-taper-command",
        title: "Use script/taper_topography.py",
        paragraphs: [
          "Run the dependency-free Python 3 helper from the Hercules_GPU repository. It reads the grid spacing from the first line, checks the supplied domain dimensions and point count, and stages another file in the same topography.in format. Create the output directory first because the script does not create parent directories. A seabed.in beside the input surface triggers automatic paired checking and repair; --seabed PATH selects another bed and --no-seabed explicitly disables this behavior.",
          "This example is configured for the bundled 40 km × 40 km benchmark. Replace every dimensional value with the actual case values. Set the PML thickness from the realized thickness printed by Hercules to standard output at startup rather than assuming the benchmark value applies to another mesh. The PML element counts and optional outer-clamp count are reported on that same startup stream.",
        ],
        codeLabel: "shell — run from the Hercules_GPU repository",
        code: "mkdir -p input/topography_pml\n\npython3 script/taper_topography.py \\\n  input/topography_raw/topography.in \\\n  --out input/topography_pml/topography.in \\\n  --lx 40000 \\\n  --ly 40000 \\\n  --buffer 6000 \\\n  --flat-margin 1500 \\\n  --pml-thickness 780 \\\n  --z-ref auto",
        note: "Existing outputs require --force. Paired mode never replaces either input, even with --force; choose a separate output directory. Surface-only mode permits in-place replacement with --force. Output symlinks are refused. Hercules requires the final surface filename to be topography.in.",
      },
      {
        id: "topography-taper-options",
        title: "Choose the taper dimensions",
        table: {
          headers: ["Option", "Meaning", "How to choose it"],
          rows: [
            ["--lx / --ly", "North-south and east-west domain lengths in metres", "Match region_length_north_m and region_length_east_m; the file spacing must divide both exactly"],
            ["--buffer", "Total inward band occupied by the flat collar plus the cosine ramp", "Make it wide enough for a gradual transition; it must exceed --flat-margin, and opposite tapers must not overlap"],
            ["--flat-margin", "Exactly flat collar measured inward from every lateral boundary", "Set it at least as wide as the realized PML, with additional margin"],
            ["--pml-thickness", "Thickness used by the helper's safety warning", "Match the startup-reported PML thickness; this option does not configure Hercules"],
            ["--z-ref", "Elevation of the flat collar", "Use a fixed elevation such as 0.0, or auto for the mean of the first interior ring when a DEM has a large elevation offset"],
            ["--seabed", "Original bed file", "Omit to detect seabed.in beside the input surface; supply a path for other names/locations"],
            ["--no-seabed", "Disable seabed detection/check/repair", "Off by default; explicit surface-only operation, mutually exclusive with --seabed"],
            ["--topography-base", "Optional topographybase_zcoord in metres", "Finite and nonnegative; positive values check wet surfaces against the solver top plane. Omission skips only this absolute-datum check"],
            ["--force", "Permit replacing existing outputs", "Paired mode never overwrites inputs; only surface-only mode permits in-place input replacement"],
          ],
        },
        note: "The script rejects 2 × buffer greater than or equal to the shorter domain extent because the tapers from opposite boundaries would overlap.",
      },
      {
        id: "topography-taper-enable",
        title: "Enable the tapered surface and verify it",
        paragraphs: [
          "Point topography_directory at the directory containing the generated topography.in. The helper's --pml-thickness value is only a consistency check, so configure C-PML independently in parameters.in.",
        ],
        codeLabel: "parameters.in — matching the example above",
        code: "abc_type                    = cpml\ncpml_thickness_m            = 780.0\ninclude_topography           = yes\ncomputation_method           = vt\ntopography_directory         = input/topography_pml",
        bullets: [
          "Require the helper's perimeter-flatness result to end in OK with max |z_boundary - z_ref| = 0.",
          "Confirm the helper reports the expected grid dimensions, flat-collar width, ramp width, reference elevation, and number of modified points.",
          "At solver startup, compare the realized PML thickness with --flat-margin and confirm the latter is larger with reasonable slack.",
          "Inspect the tapered surface or several boundary-to-interior profiles before a large production run; the outer collar should be constant, the ramp smooth, and the central model unchanged.",
        ],
        note: "Do not use the boundary-reaching topography_tilt_boundary benchmark as a production setup. It deliberately exercises the flatten-PML backstop; real surfaces should be tapered first.",
      },
      {
        id: "repair-tapered-geometry",
        title: "Pair tapered topography with explicit seabed geometry",
        paragraphs: [
          "script/taper_topography.py checks the tapered candidate automatically when seabed.in is beside the input. It copies an already valid bed byte-for-byte. If the taper creates water gaps, bed crossings or invalid wet cells, the shared repair implementation sets originally dry beds to the tapered surface and retains original wet bed tokens. Invalid original geometry, changed original water levels or invalid repaired cells stop the command before publishing outputs. Both scripts must remain together when copying these utilities to a cluster.",
          "Integrated output is the named --out surface plus sibling seabed.in and repair-summary.json. Status is valid or repaired, with candidate-problem counts, actual changes and hashes. No seabed found means surface-only output. --seabed PATH overrides discovery and --no-seabed disables it. Supply --topography-base for the absolute top-plane check; omission is reported and skips only that check. Outputs are staged before publication, but the multi-file update is not crash-atomic; keep the destination idle until success.",
          "For a surface already tapered, script/repair_tapered_geometry.py remains available as the standalone command below. Originally dry beds receive the exact tapered-surface token; originally wet beds retain their original token without re-rounding. The tapered surface is copied byte-for-byte. This prevents artificial water gaps caused by moving or rounding land elevations while leaving the bed unchanged.",
          "The Python 3 standard-library helper requires all seven command-line arguments below, with no defaults. Supply common grid spacing, origin, ordering and vertical datum; elevations are metres positive upward, X/north outer and Y/east inner. Domain lengths and topographybase_zcoord are in metres. Create the output parent first. Existing output paths are refused, with no --force option; original inputs are read-only.",
          "Both workflows use the helper's fixed 1e-8 m tolerance, independent of parameters.in and the solver's configurable geometry tolerance. Original bed crossings, changed wet-node water levels, sloping original or repaired wet cells, invalid values and incorrect counts cause failure before publishing a pair. This preserves nodal classification, not necessarily an identical interpolated coastline. It does not flatten water surfaces or check acoustic mesh-face alignment.",
        ],
        codeLabel: "shell — replace example dimensions and datum with your case",
        code: "python3 script/repair_tapered_geometry.py \\\n  --original-surface inputfiles/geometry/topography.in \\\n  --original-bed inputfiles/geometry/seabed.in \\\n  --tapered-surface inputfiles/geometry_tapered/topography.in \\\n  --out-dir inputfiles/geometry_tapered_consistent \\\n  --lx 40000 --ly 40000 --topography-base 1250",
        bullets: [
          "Outputs are topography.in, seabed.in and repair-summary.json. The summary records wet/dry counts, changed dry-bed nodes, maximum dry-bed adjustment, wet-cell counts and input/output SHA-256 hashes.",
          "Set topography_directory = inputfiles/geometry_tapered_consistent and seabed_geometry_file = inputfiles/geometry_tapered_consistent/seabed.in. Keep the other required water/material controls, regenerate meshes and mesh-dependent sources, and review source/receiver placement.",
          "The helper streams rows; memory scales with row width and the longest input line. Geometry-relative materials and hydrostatic/geostatic loads follow the resulting surface/bed pair. File validation does not establish simulation accuracy or stability. See script/README.md for formats, all arguments, overwrite behavior and limits.",
        ],
      },
    ],
  },
  "source-models": {
    eyebrow: "Setting up simulations",
    title: "Source models and source.in",
    summary:
      "Configure all source representations implemented by the current Hercules solver, from one double-couple or body force to distributed finite-fault and compliant-base inputs.",
    sections: [
      {
        id: "source-directory",
        title: "Source directory contract",
        paragraphs: [
          "Set source_directory in parameters.in to the directory that contains source.in. Relative paths are resolved from the run directory. Hercules reads the source description, generates per-rank force_process files in source_directory_output, and consumes those generated forces during time stepping.",
          "source.in selects exactly one type_of_source. Auxiliary filenames are fixed for each model, so keep the documented names and place the files beside source.in.",
        ],
        codeLabel: "parameters.in and directory layout",
        code: "# parameters.in\nsource_directory        = sourcefiles\nsource_directory_output = outputfiles/srctmp\n\n# case layout\nsourcefiles/\n├── source.in\n└── model-specific auxiliary files",
        note: "source_directory_output contains generated temporary force files, not user-authored source input. Use a separate writable directory for each concurrent run.",
      },
      {
        id: "supported-models",
        title: "Supported source models",
        table: {
          headers: ["Value", "Physical representation", "Primary use"],
          rows: [
            ["point", "One moment-tensor point source", "Compact earthquake source specified by moment, location, strike, dip, and rake"],
            ["force", "One Cartesian body force", "Reciprocity, Green's functions, impulse tests, or physical force excitation"],
            ["base_force", "Traction distributed over the bottom mesh plane", "Compliant-base and one-dimensional site-response input"],
            ["plane", "Rectangular finite fault divided into subfault cells", "Kinematic rupture with spatial slip and rake"],
            ["planewithkinks", "Vertical finite fault following a polyline trace", "Strike-slip rupture on a segmented surface trace"],
            ["srfh", "Collection of independently described rupture points", "Converted SRF-style finite-fault models with tabulated slip histories"],
          ],
        },
      },
      {
        id: "common-parameters",
        title: "Parameters common to every source",
        table: {
          headers: ["Key", "Requirement", "Meaning"],
          rows: [
            ["source_is_filtered", "Required; 0 or 1", "Enable the source-force Butterworth low-pass filter."],
            ["threshold_frequency", "Required when source_is_filtered=1", "Low-pass cutoff frequency in hertz."],
            ["number_of_poles", "Required when source_is_filtered=1", "Integer Butterworth filter order."],
            ["type_of_source", "Required", "One of point, force, base_force, plane, planewithkinks, or srfh."],
            ["source_function_type", "Required", "Time-function family used by the chosen source model."],
          ],
        },
        note: "Existing examples often retain threshold_frequency and number_of_poles when filtering is disabled. They are ignored when source_is_filtered=0.",
      },
      {
        id: "source-filter-sampling",
        title: "Source-filter sampling and start time",
        paragraphs: [
          "The source-force filter is applied by the plane, planewithkinks, and srfh source-generation paths. It uses the sampling interval of the generated force history: the common slipfunction.in interval, not simulation_delta_time_sec, for SRFH, and the solver timestep for plane and planewithkinks. Changing the solver timestep must not change an SRFH filter's physical cutoff.",
          "Older SRFH builds used the solver timestep to label source-frequency bins. If source samples were 0.01 s apart and the solver timestep was 0.005 s, a configured 1 Hz cutoff acted at 0.5 Hz. Rebuild and regenerate source-force files after updating; use the intended physical cutoff without compensating for the timestep ratio.",
          "Filtering differentiates the force history, applies a zero-phase FFT Butterworth magnitude, and integrates the retained result from time zero. Its response extends before the nominal source start. Starting a source at time zero can truncate that response, create a startup transient, and change the final applied force. Provide enough pre-source time and compare generated force histories as the buffer increases; there is no universal delay in seconds. In slipfunction.in, the second value of each source block is its start time in seconds.",
        ],
        codeLabel: "source.in: cutoff in physical hertz",
        code: "source_is_filtered = 1\nthreshold_frequency = 1.0\nnumber_of_poles = 14\n# Choose source start times with enough room for the filtered prehistory.",
        note: "Filtering changes the source spectrum. Choose a cutoff that the mesh resolves, retain the intended physical bandwidth, and check final force and startup behavior. Delaying a source does not make the zero-phase filter strictly causal relative to its nominal start.",
      },
      {
        id: "time-functions",
        title: "Source-time functions",
        paragraphs: [
          "Analytical functions scale from zero according to the selected rise-time or Ricker parameters. Tabulated discrete histories use linear interpolation and hold the final sample after the tabulated interval. Endpoint indexing is bounded even when floating-point rounding places a time just below the endpoint in the final sample index.",
        ],
        table: {
          headers: ["source_function_type", "Additional keys or file", "Behavior"],
          rows: [
            ["ramp", "average_risetime_sec", "Linear rise to the final amplitude."],
            ["sine", "average_risetime_sec", "Smooth sine-shaped rise to the final amplitude."],
            ["quadratic", "average_risetime_sec", "Piecewise-quadratic smooth rise."],
            ["exponential", "average_risetime_sec", "Asymptotic exponential rise."],
            ["ricker", "ricker_Ts and ricker_Tp", "Ricker pulse shifted by Ts; Tp sets its time scale and an approximate characteristic frequency of 1/Tp."],
            ["discrete", "Model-specific history file", "Tabulated samples with nt1, initial delay, and sample interval."],
          ],
        },
        note: "The discrete path is fully defined for force, base_force, and srfh. A point source does not load a discrete history file; use an analytical function for type_of_source=point.",
      },
      {
        id: "point-source",
        title: "Point moment-tensor source",
        paragraphs: [
          "The point model represents a double-couple source at one location. Supply either scalar moment moment_amplitude in N m or moment_magnitude Mw. When moment_magnitude is present Hercules computes M0 = 10^(1.5 Mw + 9.1) N m.",
          "Point and single-point body-force histories preserve every loaded node when source generation needs multiple memory batches. The binary force file stores one complete node list followed by all time samples in that same node order; a longer simulation must not change the earlier source samples.",
        ],
        table: {
          headers: ["Key", "Requirement", "Meaning"],
          rows: [
            ["moment_amplitude", "Required when moment_magnitude is absent", "Scalar seismic moment M0 in N m."],
            ["moment_magnitude", "Alternative to moment_amplitude", "Moment magnitude Mw; if both forms are present, the magnitude path is selected."],
            ["lonlat_or_cartesian", "Required; 0 or 1", "0 selects geographic coordinates; 1 selects local Cartesian coordinates."],
            ["hypocenter_lat_deg / hypocenter_long_deg", "Required when coordinate mode is 0", "Geographic hypocentre latitude and longitude in degrees."],
            ["hypocenter_x / hypocenter_y", "Required when coordinate mode is 1", "Local x and y coordinates in metres."],
            ["hypocenter_depth_m", "Required", "Depth in metres, positive downward."],
            ["source_strike_deg", "Required", "Fault strike in degrees clockwise from true north for geographic sources, or solver +X for Cartesian sources."],
            ["source_dip_deg", "Required", "Fault dip in degrees."],
            ["source_rake_deg", "Required", "Slip rake in degrees."],
          ],
        },
        codeLabel: "source.in — point source",
        code: "source_is_filtered   = 0\ntype_of_source       = point\nsource_function_type = ricker\nricker_Tp            = 0.40\nricker_Ts            = 0.60\n\nmoment_amplitude     = 1.0e15\nlonlat_or_cartesian  = 1\nhypocenter_x         = 1000.0\nhypocenter_y         = 1000.0\nhypocenter_depth_m   = 800.0\nsource_strike_deg    = 0.0\nsource_dip_deg       = 90.0\nsource_rake_deg      = 0.0",
      },
      {
        id: "force-source",
        title: "Single-point body-force source",
        paragraphs: [
          "The force model applies F(t) = (Fx, Fy, Fz) s(t) at one point and distributes it to the containing eight-node element with trilinear shape functions. The simulation frame uses local solver X/Y and Z positive downward; X/Y coincide with UTM north/east at zero domain azimuth.",
        ],
        table: {
          headers: ["Key", "Requirement", "Meaning"],
          rows: [
            ["lonlat_or_cartesian", "Required; 0 or 1", "Select geographic or local Cartesian position."],
            ["hypocenter_lat_deg / hypocenter_long_deg", "Required for mode 0", "Geographic source latitude and longitude."],
            ["hypocenter_x / hypocenter_y", "Required for mode 1", "Local source x and y in metres."],
            ["hypocenter_depth_m", "Required", "Source depth in metres, positive downward."],
            ["force_x_n", "Required", "Northward Cartesian force component in newtons."],
            ["force_y_n", "Required", "Eastward Cartesian force component in newtons."],
            ["force_z_n", "Required", "Downward Cartesian force component in newtons."],
          ],
        },
        codeLabel: "source.in — analytical force",
        code: "source_is_filtered   = 0\ntype_of_source       = force\nsource_function_type = quadratic\naverage_risetime_sec = 0.5\n\nlonlat_or_cartesian  = 0\nhypocenter_lat_deg   = 34.0\nhypocenter_long_deg  = -118.5\nhypocenter_depth_m   = 1000.0\nforce_x_n            = 0.0\nforce_y_n            = 0.0\nforce_z_n            = 1.0e15",
      },
      {
        id: "force-history",
        title: "Discrete forcefunction.in",
        paragraphs: [
          "For source_function_type=discrete, forcefunction.in contains a unit or otherwise normalized waveform. Hercules linearly interpolates the samples and multiplies them by force_x_n, force_y_n, and force_z_n.",
        ],
        codeLabel: "forcefunction.in — five samples",
        code: "5       # nt1; must be at least 2\n0.0     # tinit, seconds\n0.01    # dt, seconds; must be positive\n0.0\n0.5\n1.0\n0.5\n0.0",
        note: "A peak-normalized history keeps physical amplitude in the force components. Non-unit sample scaling is mathematically allowed but should be documented explicitly.",
      },
      {
        id: "base-force",
        title: "Distributed compliant-base force",
        paragraphs: [
          "The base_force model converts a base velocity history to nodal force using F_i(t) = scale rho Vs A_i v(t) direction_i, where A_i is the tributary bottom-face area. It is intended for flat-bottom validation columns and compliant-base site-response input.",
        ],
        table: {
          headers: ["Key", "Requirement", "Meaning"],
          rows: [
            ["source_function_type", "Required: discrete", "base_force accepts only a tabulated velocity history."],
            ["base_force_rho_kgm3", "Required; positive", "Base or bedrock density in kg/m3."],
            ["base_force_vs_mps", "Required; positive", "Base or bedrock shear-wave speed in m/s."],
            ["base_force_direction_x/y/z", "Required; nonzero vector", "Propagation/polarization direction in the solver-local X/Y, Z-down frame; Hercules normalizes the vector."],
            ["base_force_scale", "Optional; default 1", "Use 1 for an incident/within motion and commonly 2 for an outcrop motion."],
            ["base_force_depth_m", "Optional", "Target base depth in metres; omission selects the deepest mesh-node plane."],
            ["base_force_tolerance_m", "Optional", "Positive depth-matching tolerance; omission uses 1% of the minimum mesh edge."],
          ],
        },
        codeLabel: "source.in — compliant base",
        code: "source_is_filtered        = 0\ntype_of_source            = base_force\nsource_function_type      = discrete\nbase_force_rho_kgm3       = 2200.0\nbase_force_vs_mps         = 450.0\nbase_force_direction_x    = 1.0\nbase_force_direction_y    = 0.0\nbase_force_direction_z    = 0.0\nbase_force_scale          = 1.0\n# base_force_depth_m      = 30000.0\n# base_force_tolerance_m  = 0.01",
        note: "baseforcefunction.in uses the same nt1, tinit, dt, samples layout as forcefunction.in, but its samples are base velocity in m/s rather than a unit force-time function.",
      },
      {
        id: "finite-plane",
        title: "Rectangular finite-fault plane",
        paragraphs: [
          "The plane source divides a rectangular fault into cells_along_strike by cells_down_dip subfaults. Rupture time follows distance from the hypocentre divided by extended_average_rupture_velocity; each time window then adds its value from time_windows.",
        ],
        table: {
          headers: ["Key", "Requirement", "Meaning"],
          rows: [
            ["number_of_time_windows", "Required", "Number of slip/rake matrix blocks."],
            ["time_windows", "Required list", "One additional delay in seconds per time window."],
            ["extended_cell_size_along_strike_m", "Required", "Subfault-cell length along strike."],
            ["extended_cell_size_down_dip_m", "Required", "Subfault-cell length down dip."],
            ["extended_cells_along_strike / extended_cells_down_dip", "Required", "Number of fault cells in each direction."],
            ["extended_lat_deg / extended_long_deg / extended_depth_m", "Required", "Geographic location and depth of the plane origin."],
            ["extended_strike_deg / extended_dip_deg", "Required", "Fault-plane orientation in degrees. Under PROJ, strike is clockwise from true north at the fault origin and is transformed into solver axes."],
            ["extended_hypocenter_along_strike_m", "Required", "Hypocentre distance from the origin along strike."],
            ["extended_hypocenter_down_dip_m", "Required", "Hypocentre distance from the origin down dip."],
            ["extended_average_rupture_velocity", "Required", "Average rupture speed in m/s."],
            ["extended_isminimumedgeautomatic", "Required; 0 or 1", "1 derives source discretization from the mesh; 0 uses extended_minimum_edge_m."],
            ["extended_minimum_edge_m", "Required when automatic=0", "Manual minimum source-discretization edge in metres."],
          ],
        },
        codeLabel: "source.in — 2 × 2 plane",
        code: "source_is_filtered                     = 0\ntype_of_source                         = plane\nsource_function_type                   = quadratic\naverage_risetime_sec                   = 0.6\nnumber_of_time_windows                 = 1\ntime_windows                           = 0.0\nextended_cell_size_along_strike_m     = 1000.0\nextended_cell_size_down_dip_m         = 1000.0\nextended_cells_along_strike           = 2\nextended_cells_down_dip               = 2\nextended_lat_deg                      = 34.3870\nextended_long_deg                     = -118.5976\nextended_depth_m                      = 5000.0\nextended_strike_deg                   = 122.0\nextended_dip_deg                      = 40.0\nextended_hypocenter_along_strike_m    = 1000.0\nextended_hypocenter_down_dip_m        = 1000.0\nextended_average_rupture_velocity     = 3000.0\nextended_isminimumedgeautomatic       = 1",
      },
      {
        id: "plane-matrices",
        title: "Plane slip.in and rake.in matrices",
        paragraphs: [
          "slip.in and rake.in each contain number_of_time_windows consecutive matrices. Every matrix has extended_cells_down_dip rows and extended_cells_along_strike values per row. Slip values are metres; rake values are degrees.",
        ],
        codeLabel: "one 2 × 2 window",
        code: "# slip.in\n0.40 0.60\n0.25 0.10\n\n# rake.in\n90.0 90.0\n90.0 90.0",
        note: "For multiple windows, append the next complete matrix immediately after the previous one in both files, using the same window order as time_windows.",
      },
      {
        id: "kinked-plane",
        title: "Kinked vertical finite fault",
        paragraphs: [
          "planewithkinks follows a polyline surface trace. Hercules computes the along-strike cell size from total trace length divided by extended_cells_along_strike and fixes dip to 90 degrees; strike is derived locally from each trace segment.",
          "It uses the plane model's time windows, down-dip spacing, cell counts, hypocentre distances, rupture velocity, automatic/manual minimum-edge control, and slip.in/rake.in matrices. It replaces the rectangular origin and strike with the kink list.",
        ],
        table: {
          headers: ["Key", "Requirement", "Meaning"],
          rows: [
            ["extended_number_of_kinks", "Required; at least 2", "Number of surface-trace vertices."],
            ["extended_kinks", "Required pairs", "longitude latitude pairs in trace order."],
            ["extended_depth_m", "Required", "Top depth of the vertical fault."],
            ["extended_cell_size_down_dip_m", "Required", "Down-dip cell size in metres."],
            ["extended_cells_along_strike / extended_cells_down_dip", "Required", "Fault grid dimensions."],
          ],
        },
        codeLabel: "source.in — trace excerpt",
        code: "type_of_source                       = planewithkinks\nsource_function_type                 = exponential\naverage_risetime_sec                 = 0.5\nnumber_of_time_windows               = 1\ntime_windows                         = 0.0\nextended_cell_size_down_dip_m       = 200.0\nextended_cells_along_strike         = 3\nextended_cells_down_dip             = 20\nextended_depth_m                    = 0.0\nextended_hypocenter_along_strike_m  = 1000.0\nextended_hypocenter_down_dip_m      = 2000.0\nextended_average_rupture_velocity   = 3000.0\nextended_isminimumedgeautomatic     = 1\nextended_number_of_kinks            = 4\nextended_kinks =\n  -117.50 34.29\n  -117.24 34.15\n  -116.77 33.99\n  -116.48 33.92",
        note: "extended_dip_deg, extended_strike_deg, extended_lat_deg, extended_long_deg, and extended_cell_size_along_strike_m are not used by planewithkinks in the current solver path.",
      },
      {
        id: "srfh",
        title: "SRFH finite-fault source",
        paragraphs: [
          "srfh reads one record per rupture point from synchronized auxiliary files. Use source_function_type=discrete; each point carries its own onset time and tabulated slip history.",
          "Hercules generates SRFH force histories at the common slipfunction.in sample interval, then interpolates those forces at solver time steps. Reducing simulation_delta_time_sec does not change the source-history sample interval. A cumulative slip history normally ends at its final slip fraction, which is held after tinit + (nt1 - 1) * dt.",
        ],
        table: {
          headers: ["File or key", "Entries per point", "Meaning"],
          rows: [
            ["number_of_point_sources", "One value in source.in", "Number of point-source records read from every auxiliary file."],
            ["coords.in", "longitude latitude depth_m", "Geographic location; depth is positive downward."],
            ["area.in", "area_m2", "Subfault area."],
            ["strike.in", "strike_deg", "Strike in degrees."],
            ["dip.in", "dip_deg", "Dip in degrees."],
            ["rake.in", "rake_deg", "Rake in degrees."],
            ["slip.in", "slip_m", "Final slip amplitude."],
            ["slipfunction.in", "nt1, tinit, dt, then nt1 samples", "Tabulated normalized slip-time history for each point."],
          ],
        },
        codeLabel: "source.in — SRFH controller",
        code: "source_is_filtered      = 1\nthreshold_frequency     = 4.5\nnumber_of_poles         = 14\ntype_of_source          = srfh\nsource_function_type    = discrete\nnumber_of_point_sources = 2",
        note: "All SRFH points must use the same dt in slipfunction.in. Hercules aborts when adjacent point records have different sample intervals.",
      },
      {
        id: "srfh-files",
        title: "SRFH auxiliary-file example",
        paragraphs: [
          "Each scalar file contains exactly number_of_point_sources values in the same point order. slipfunction.in concatenates one variable-length block per point.",
        ],
        codeLabel: "two-point layout",
        code: "# coords.in\n-118.50 34.00 1000.0\n-118.49 34.00 1000.0\n\n# area.in       strike.in  dip.in  rake.in  slip.in\n# Each item below represents the corresponding two-line scalar file.\narea:   1.0e6, 1.0e6\nstrike: 0.0,   0.0\ndip:    90.0,  90.0\nrake:   0.0,   0.0\nslip:   0.50,  0.35\n\n# slipfunction.in: block for point 0, then point 1\n3\n0.00\n0.01\n0.0\n1.0\n0.0\n3\n0.05\n0.01\n0.0\n1.0\n0.0",
        note: "The labels in this compact illustration explain the separate scalar files; do not place the words area, strike, dip, rake, or slip inside the actual numeric files.",
      },
      {
        id: "validation",
        title: "Source validation checklist",
        bullets: [
          "Confirm source coordinates fall inside the domain and use the intended geographic or Cartesian convention.",
          "Check positive-down Z and the local solver X/Y frame before interpreting force or station components.",
          "Plot every analytical or tabulated source-time function at its native sample interval.",
          "Verify scalar moment, force, slip, area, density, velocity, and time units independently.",
          "For finite faults, check matrix dimensions and point ordering before starting a large run.",
          "Archive source.in and every auxiliary file with parameters.in and the monitor log.",
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
          headers: ["Output", "Primary files", "Format and use"],
          rows: [
            ["Stations", "station.N; station_mediums.txt, or station_output_binary_file", "ASCII histories by default; station_binary=yes writes one float64 displacement file"],
            ["Planes", "planedisplacements.N", "Raw float64 displacement snapshots on a regular plane grid"],
            ["DRM boxes", "drmdisplacements.N", "One raw float64 displacement file per ten-face structured box"],
            ["Selected mesh", "mesh_coordinates.RANK; mesh_data.RANK", "Binary coordinates and material properties for a Cartesian subdomain"],
            ["Checkpoints", "checkpoint.out0; checkpoint.out1; checkpoint.in", "Binary rotating solver states and restart input"],
            ["Monitor", "monitor_file, default monitor.txt", "Configuration, execution-path, timing, and diagnostic record"],
            ["Distributed fields", "Configured output paths", "Optional parallel mesh, displacement, velocity, statistics, and debug streams"],
          ],
        },
      },
      {
        id: "station-configuration",
        title: "Configure station output",
        paragraphs: [
          "output_stations contains latitude, longitude, and depth in that order. With station_binary=no (the default), Hercules writes one ASCII station.N file for each row. The sampling interval is output_stations_print_rate multiplied by simulation_delta_time_sec.",
          "ASCII station writes default to a 64 KiB buffer per locally owned station, allocated once during setup. Set station_ascii_buffer_kib to change this size: an optional integer in KiB (1024 bytes), default 64, accepted range 1..2097151 (INT_MAX/1024). The value is validated but unused with station_binary=yes. Buffers write when full and flush at normal completion, replacing the former flush every ten samples. Every requested sample and column is retained without filtering, averaging, or downsampling; timestamps and text formatting stay unchanged. This applies to acoustic, topography, and solid station histories. The number of rows per write depends on the enabled columns and their formatted values.",
          "CPU buffer memory per rank is station_ascii_buffer_kib × 1024 × locally owned stations bytes, plus stream/container overhead; all output columns share that station's buffer. For 1,000 local stations, the default uses 62.5 MiB; setting 16 uses 15.625 MiB. Smaller buffers save memory and write more frequently; larger buffers allow less frequent writes but delay live-file visibility. Binary output and ranks without stations allocate no ASCII station buffers.",
          "Live station files can lag the simulation and may temporarily end in a partial row. An interrupted or aborted run can lose the unflushed tail. Normal completion checks write errors and includes the final flush in the Print Stations and Solver I/O timers. The startup monitor reports the buffer size. Binary station output uses its separate MPI-IO path; plane output and ETA estimation are unchanged.",
          "Displacement is always present. Velocity is appended when velocity or acceleration output is requested; acceleration, pressure, pore-pressure ratio, strain, and stress columns are then appended according to their switches. Strain and stress use XX, YY, ZZ, XY, YZ, XZ order.",
          "Station pressure eligibility is finalized after solver setup from the active BrickUP map and receiver saturation. Elastic-only CPML receivers are labeled elastic with pressure undefined, even when they share nodes with active poroelastic soil. Their Pressure(Pa) and Ru(-) values are nan; skeleton motion retains the normal solid interpolation. station_mediums.txt is written after this final classification and refreshed again when a ramped geostatic reference becomes available.",
          "For nonlinear BrickUP runs, print_station_pore_pressure_ratios=yes appends the conventional experimental ratio Ru(-) = Δu/σ′v0. Δu is the same excess pore pressure used by the Pressure(Pa) column. With approximate geostatic initialization, σ′v0 is derived at the exact station depth from the integrated total-overburden column minus Biot-weighted hydrostatic pressure. With a ramped geostatic phase, Ru remains nan during the ramp because its numerical reference does not yet exist. On the last geostatic force evaluation, Hercules fixes σ′v0 from the committed vertical skeleton stress at the owning nonlinear element's station-output quadrature point—the same constitutive state used by its nonlinear stress columns and by RIVA-Sand initialization—before the first dynamic output row. The final denominator is recorded as sigma_v0_eff_pa in station_mediums.txt; other undefined ratios are also nan. Hercules rejects the ratio request when approximate_geostatic_state=no and both geostatic times are zero, because that configuration supplies no initial effective-stress reference. The ratio is not clipped, so dilation may produce negative values and values may exceed one.",
        ],
        codeLabel: "parameters.in — two stations",
        code: "number_output_stations                 = 2\noutput_stations_print_rate             = 10\noutput_stations_directory              = output/stations\nstation_ascii_buffer_kib              = 64\nstation_mediums_output_file            = output/stations/station_mediums.txt\nprint_station_velocities               = yes\nprint_station_accelerations            = yes\nprint_station_pressures                = no\n# Set yes with nonlinear BrickUP to append Ru(-).\nprint_station_pore_pressure_ratios     = no\nprint_station_strains                  = no\nprint_station_stresses                 = no\n\n# latitude_deg longitude_deg depth_m\noutput_stations =\n  34.0000 -118.5000 0.0\n  34.0100 -118.4900 100.0",
        note: "With topography enabled, station depth is shifted to the local topographic surface. Check station_mediums.txt to confirm the owning medium, MPI rank, element, local coordinates, Vp, Vs, density, and sigma_v0_eff_pa reference. Ru(-) is nan for acoustic, elastic, unsaturated, CPML, surface-zero-stress, or otherwise inadmissible station references.",
      },
      {
        "id": "station-binary",
        "title": "Write all station displacements to one binary file",
        "paragraphs": [
          "Set station_binary=yes and station_output_binary_file to record the existing output_stations list in one file. station_binary defaults to no, which preserves the existing station.N output functions. Binary mode requires number_output_stations>0 and a positive output_stations_print_rate. Coordinates still come from output_stations as latitude longitude depth, with the existing geographic conversion and surface-relative building/topography shifts.",
          "The file is headerless native-endian float64 [snapshot, station, XYZ] displacement in metres, where XYZ means solver-local X/Y, Z down. PROJ writes the geographic mapping in coordinate_frame.txt in output_stations_directory. Input row order is preserved, including duplicate rows. Every station must have exactly one owner. The existing station interpolation handles solid/skeleton displacement, topography cut cells, and reconstructed scalar-acoustic motion.",
          "Binary mode writes displacement only. It ignores print_station_* switches, nonlinear diagnostic columns, output_stations_directory, and station_mediums_output_file. No per-station text files or station metadata sidecar are created. The filename is required only with station_binary=yes and is ignored otherwise. Paths contain 1 to 240 bytes, resolve from the run directory without tilde expansion, and have parent directories created automatically."
        ],
        "codeLabel": "parameters.in — one displacement file, using the existing station table",
        "code": "station_binary             = yes\nstation_output_binary_file = output/station_displacements.0\nnumber_output_stations     = 2\noutput_stations_print_rate = 10\n\n# latitude_deg longitude_deg depth_m\noutput_stations =\n  34.0000 -118.5000   0.0\n  34.0100 -118.4900 100.0",
        "note": "One snapshot is 24 × number_output_stations bytes. Sampling uses output_stations_print_rate × simulation_delta_time_sec. Step zero is saved in a fresh run, followed by steps divisible by the print rate before they advance. Each run truncates station_output_binary_file, including checkpoint restarts; preserve earlier segments in a different path."
      },
      {
        "id": "read-station-binary",
        "title": "Read the station binary",
        "paragraphs": [
          "Use script/read_drm.py output/station_displacements.0 --points 2 --dt 0.01 for the configuration above with solver dt=0.001 s. This shared flat-displacement reader checks complete snapshots and finite values; --points is the station count. Retain parameters.in to identify coordinates and timing. Set --first-time for a restart segment to the first sampled step times solver dt.",
          "The station binary path reuses existing sparse station GPU downloads and writes one collective MPI-IO snapshot at each requested step, including participation by ranks without stations. It is independent of drm_output_batch_size. Stations can run alongside multiple DRM boxes with separate sampling rates; their output paths must be distinct."
        ],
        "codeLabel": "Python — station order matches output_stations",
        "code": "import numpy as np\n\nnstations = 2  # number_output_stations\nraw = np.fromfile(\"output/station_displacements.0\", dtype=np.float64)\nif not raw.size or raw.size % (nstations * 3):\n    raise ValueError(\"incomplete snapshot or wrong station count\")\nu = raw.reshape(-1, nstations, 3)\ntime = np.arange(u.shape[0]) * 10 * 0.001  # print rate * solver dt\n# u[:, 0, :] is displacement for the first output_stations row."
      },
      {
        id: "nonlinear-station-output",
        title: "Interpret nonlinear and plastic station output",
        paragraphs: [
          "When include_nonlinear_analysis=yes and station_binary=no, Hercules appends nonlinear solid-state columns to station files. These contain paired strain and stress tensor components plus volumetric strain/stress, yield-function state Fs, accumulated inelastic measure Kappa, one-dimensional strain/stress diagnostics Gamma1d and Tao1d, and GGmax.",
          "Use these histories to identify yield onset, hysteresis, stiffness reduction, and permanent/plastic response. Read the station header instead of assuming fixed column offsets because velocity, acceleration, pressure, pore-pressure ratio, linear strain/stress, and nonlinear columns are controlled independently.",
        ],
        codeLabel: "nonlinear station columns",
        code: "Epsilon_XX Sigma_XX  Epsilon_YY Sigma_YY  Epsilon_ZZ Sigma_ZZ\nEpsilon_KK Sigma_KK  Epsilon_XY Sigma_XY  Epsilon_YZ Sigma_YZ\nEpsilon_XZ Sigma_XZ  Fs Kappa Gamma1d Tao1d GGmax",
      },
      {
        id: "read-stations",
        title: "Read station.N with Python",
        paragraphs: [
          "The header records the actual columns. This example matches the configuration above: time, three displacement components, three velocity components, and three acceleration components.",
          "The standalone script/read_stations.py accepts command-line arguments; run the command shown below without its leading comment marker. script/README.md lists the utilities and their dependencies.",
        ],
        codeLabel: "Python — script/read_stations.py",
        code: "# Runnable version from the repository root:\n# python3 script/read_stations.py output/stations/station.0\n\nimport numpy as np\n\nstation = np.loadtxt(\"output/stations/station.0\", comments=\"#\")\ntime = station[:, 0]\ndisplacement = station[:, 1:4]\nvelocity = station[:, 4:7]\nacceleration = station[:, 7:10]\n\nprint(f\"{len(time)} samples, dt_out={np.diff(time).mean():.6g} s\")\nprint(\"peak |u| =\", np.linalg.norm(displacement, axis=1).max())",
        note: "Do not hard-code later column offsets when pressure, strain, stress, or nonlinear output is enabled. Read the header and document the exact column selection used in post-processing.",
      },
      {
        id: "plane-configuration",
        title: "Configure plane output",
        paragraphs: [
          "Each output_planes row defines a regular grid as latitude, longitude, depth, spacing along strike, number along strike, spacing down dip, number down dip, strike, and dip. The origin is the first grid point.",
          "Every backend preserves the planedisplacements.N raw binary layout: float64 x, y, z displacement triples for all grid points at one output time, followed by the next output time.",
        ],
        codeLabel: "parameters.in — one horizontal plane",
        code: "number_output_planes       = 1\nplanes_input_file          = parameters.in\noutput_planes_print_rate   = 20\noutput_planes_directory    = output/planes\n\n# lat lon depth  d_strike n_strike  d_dip n_dip  strike dip\noutput_planes =\n  34.0000 -118.5000 0.0  100.0 101  100.0 81  0.0 0.0",
      },
      {
        id: "read-planes",
        title: "Read planedisplacements.N with Python",
        paragraphs: [
          "The file has no header, so the grid dimensions must come from the exact parameters.in used for the run. In the reshaped array below, the indices are output time, along-strike point, down-dip point, and Cartesian component.",
          "The standalone script/read_planes.py accepts command-line arguments; run the command shown below without its leading comment marker. script/README.md lists the utilities and their dependencies.",
        ],
        codeLabel: "Python — script/read_planes.py",
        code: "# Runnable version from the repository root:\n# python3 script/read_planes.py output/planes/planedisplacements.0 --nx 101 --ny 81 --dt 0.02\n\nimport numpy as np\n\nn_strike = 101\nn_dip = 81\nprint_rate = 20\nsimulation_dt = 0.001\n\nraw = np.fromfile(\"output/planes/planedisplacements.0\", dtype=np.float64)\nvalues_per_snapshot = n_strike * n_dip * 3\nif raw.size % values_per_snapshot:\n    raise ValueError(\"file size does not match the configured plane grid\")\n\nu = raw.reshape(-1, n_strike, n_dip, 3)\ntime = np.arange(u.shape[0]) * print_rate * simulation_dt\nvertical = u[..., 2]\nprint(u.shape, time[-1], np.abs(vertical).max())",
        note: "reconstruct_acoustic_plane_displacements=yes reconstructs motion at scalar-acoustic plane points from -grad(p)/rho. With acoustic_formulation=velocity_potential, it first computes p = -(phi[n] - phi[n-1])/dt, matching station reconstruction. Elastic points continue to use the elastic displacement field in the same file.",
      },
      {
        "id": "drm-box-configuration",
        "title": "Record multiple DRM boxes",
        "paragraphs": [
          "Use number_output_drm_boxes, output_drm_boxes_print_rate, output_drm_boxes_directory, and the output_drm_boxes table, following the plane-output style. Each row writes one headerless native-float64 drmdisplacements.N file plus a matching ASCII drmcoordinates.N file, numbered from zero in row order. Coordinate output is automatic; no extra parameter is needed. Each file contains ten structured faces: front, back, left, right, and bottom on both exterior and interior sides. There is no top lid.",
          "Each row has seven values: drm_offset_x drm_offset_y length_x length_y depth drm_edgesize strk. The offsets locate the top center in solver-local X/Y metres from the domain origin. The three dimensions describe the exterior box. With h=drm_edgesize, the interior dimensions are (length_x − 2h, length_y − 2h, depth − h), centered on the same X/Y location with the same top elevation. Spacing must be positive; lengths must exceed 2h and depth must exceed h; all dimensions must be integer multiples of h.",
          "Strike rotates the box clockwise from solver +X in degrees: for box-local horizontal coordinates (x,y), global X = center_x + x*cos(strk) − y*sin(strk), and global Y = center_y + x*sin(strk) + y*cos(strk). Zero strike aligns local X/Y with solver X/Y. The output displacement components remain solver X/Y/Z; they are not rotated into box-local components.",
          "This independent recorder works with implement_drm=no. It requires no drm_boundary or standalone drm_offset_x/y/drm_edgesize keys. Sampling does not refine the mesh. All generated points must lie inside the simulation domain and have exactly one owning solid element. A box centered at (0,0) extends to negative X/Y at zero strike; choose a center that keeps the full rotated footprint inside the domain."
        ],
        "codeLabel": "parameters.in — two boxes; all distances in metres, strike in degrees",
        "code": "number_output_drm_boxes       = 2\noutput_drm_boxes_print_rate   = 10\noutput_drm_boxes_directory    = outputfiles/drm_boxes_fullDomain\n\n# center_x center_y length_x length_y depth spacing strike\noutput_drm_boxes =\n  1000.0 2000.0 100.0 100.0 50.0 1.0  0.0\n  2000.0 3000.0 200.0 100.0 80.0 2.0 30.0\n\n# Optional: buffered snapshots per file (default 16).\ndrm_output_batch_size = 16",
        "note": "Each box produces a file pair in outputfiles/drm_boxes_fullDomain: drmdisplacements.0 with drmcoordinates.0, then drmdisplacements.1 with drmcoordinates.1. On flat terrain, the first exterior box spans X=950..1050, Y=1950..2050, depth=0..50; its interior spans X=951..1049, Y=1951..2049, depth=0..49. Building/topography depth shifts follow plane conventions at each point; topography can therefore warp the otherwise planar grids. The geometry is case-specific and must be checked against the simulation domain."
      },
      {
        "id": "drm-box-face-layout",
        "title": "Ten structured face grids in each snapshot",
        "paragraphs": [
          "Each snapshot concatenates the exterior front, back, left, right, bottom, then the interior front, back, left, right, bottom. Front/back mean negative/positive box-local X; left/right mean negative/positive box-local Y. All coordinates increase along each grid axis. Each face includes both endpoints of each axis, and shared edges/corners are repeated in each adjoining face so every face remains rectangular.",
          "For exterior interval counts nx=length_x/h, ny=length_y/h, nz=depth/h, use the shapes below. Face data are C-order [row, column, XYZ], with the column index varying fastest. A snapshot contains N = 2*(ny+1)*(nz+1) + 2*(nx+1)*(nz+1) + (nx+1)*(ny+1) exterior samples plus the same formula with nx−2, ny−2, nz−1 for the interior. Faces are observation grids, not the solver octree mesh.",
          "For a 100 × 100 × 50 m exterior box with h=1 m, the interior is 98 × 98 × 49 m. The ten faces contain 60,406 samples, or 1,449,744 bytes per displacement snapshot. The setup log prints each face name, row/column counts, and zero-based point offset. All ten faces share one binary/ASCII file pair."
        ],
        "table": {
          "headers": [
            "Face (exterior, then interior)",
            "Grid axes",
            "Exterior rows × columns",
            "Interior rows × columns"
          ],
          "rows": [
            [
              "front (−local X)",
              "local Y, Z",
              "(ny+1) × (nz+1)",
              "(ny−1) × nz"
            ],
            [
              "back (+local X)",
              "local Y, Z",
              "(ny+1) × (nz+1)",
              "(ny−1) × nz"
            ],
            [
              "left (−local Y)",
              "local X, Z",
              "(nx+1) × (nz+1)",
              "(nx−1) × nz"
            ],
            [
              "right (+local Y)",
              "local X, Z",
              "(nx+1) × (nz+1)",
              "(nx−1) × nz"
            ],
            [
              "bottom",
              "local X, local Y",
              "(nx+1) × (ny+1)",
              "(nx−1) × (ny−1)"
            ]
          ]
        },
        "note": "The file stores solver-local X/Y, Z/down displacement in metres, in native-endian float64 triples. The companion ASCII table records the actual sample positions. Keep the input box row to identify face shapes and offsets. No timestamps, velocities, accelerations, pressure, or stress are embedded."
      },
      {
        "id": "drm-coordinate-file",
        "title": "ASCII coordinates and materials for every DRM sample",
        "paragraphs": [
          "drmcoordinates.N has comment headers beginning with #, followed by eight whitespace-separated columns X_m Y_m Z_m Vs_m_s Vp_m_s density_kg_m3 Qs Qp. Coordinates are solver-local X/Y, Z/positive-down in metres after strike rotation and the same building/topography depth shifts used for interpolation. Z includes those shifts in the solver coordinate frame. Data row i (zero-based, excluding comments) matches binary sample i in every snapshot, including repeated edge/corner rows; points are never sorted or deduplicated.",
          "The appended properties come from the initial containing mesh element, using the same unique owner as displacement sampling. Vs/Vp are in m/s and density is in kg/m³, including mesh-material adjustments and applicable phase-velocity correction. Properties are not interpolated across interfaces or queried anew from the original model. Poroelastic points report skeleton/drained mesh inputs and mesh density, not fast/slow Biot wave speeds or separately configured fluid density.",
          "Qp reports the selected nominal P-wave Q after negative-Qk fallback, use_infinite_qk and high-Q cutoffs; it is not the input Qp. Qs reports the selected shear Q, including table-entry selection for legacy bkt. Using the final exported Vs/Vp (after any phase-velocity correction), L = (4/3)(Vs/Vp)² and 1/Qp = L/Qs + (1−L)/Qk. Disabled branches contribute exactly zero inverse Q; the internal Qk=1000 sentinel does not contribute finite bulk loss. For final Vs=1000 m/s, Vp=1500 m/s, Qs=100 and lossless bulk, the file reports Qp=168.75 even if the input was 200. No input-Qp, Qk or fallback-flag columns are added.",
          "inf means zero selected material loss: Qs is inf for disabled shear, and Qp is inf when both branches are disabled. With type_of_damping=none both are inf for ordinary solids. Rayleigh/mass damping and building/foundation overrides produce nan because no frequency-independent BKT Q is selected. These nominal values do not measure the frequency-dependent fitted BKT response, evolving nonlinear elastoplastic/plastic and hysteretic loss, CPML absorption or poroelastic attenuation. Where BKT is excluded from time stepping, the values describe only the selected linear reference material. BKT input admissibility checks still apply.",
          "Elastic, viscoelastic, nonlinear and poroelastic solid samples share this metadata path, including supported solid samples in acoustic–elastic, poroelastic–acoustic and poroelastic–elastic runs. Scalar-acoustic samples remain unsupported. CPU and synthetic MPI checks do not validate standalone or coupled CUDA physics; target-system comparisons remain necessary.",
          "Rebuild the solver and regenerate older meshes with mesh_input_mode=0: Selected Qs/Qp persist in element and saved mesh payloads (16 additional host bytes per element), including unified rank redistribution. Sharded manifests and MPI-IO archives use format version 2; incompatible payloads are rejected before import. Regenerate earlier development meshes that stored input Qp, and regenerate after changing damping settings; import preserves the saved material selection. No new parameter is needed. The displacement binary is unchanged; use usecols=(0, 1, 2) for position-only readers.",
          "The header gives the point count, axes and row convention. Values use 17 significant digits with a locale-independent decimal point, allowing exact float64 round trips. Solver rank zero writes and closes the ASCII file once during setup after point ownership validation and a bounded-chunk MPI material gather. Every run overwrites it, including restarts; keep each coordinate/binary pair together. This adds no per-step I/O or GPU transfers."
        ],
        "codeLabel": "Synthetic material illustration for box 0 — values are not defaults or production calibration",
        "code": "# X_m Y_m Z_m Vs_m_s Vp_m_s density_kg_m3 Qs Qp\n950 1950 0 300 700 1800 25 50\n950 1950 1 300 700 1800 25 50\n950 1950 2 300 700 1800 25 50",
        "note": "Load all eight columns with numpy.loadtxt, which skips # comments. The first three columns are positions and the last five are materials. The row count must equal the number of binary samples per snapshot. Use the same face offsets to reshape positions, materials and displacements."
      },
      {
        "id": "read-drm",
        "title": "Read a DRM box as ten arrays",
        "paragraphs": [
          "script/read_drm.py accepts --box length_x length_y depth spacing and reports the shape and offset of every face. It scans completed files for complete snapshots and finite displacement values using the Python standard library. Its box_layout helper also provides the face slices for NumPy processing, as shown below.",
          "The saved interval is output_drm_boxes_print_rate multiplied by simulation_delta_time_sec. A fresh run saves step zero, then steps divisible by the print rate before advancing those steps. With batching enabled, an in-progress file may be missing its buffered tail. Each run truncates its output files, including on restart; choose a different directory for a restart segment and use --first-time for the time of its first saved step."
        ],
        "codeLabel": "From the repository root — reshape each face independently",
        "code": "# python3 script/read_drm.py outputfiles/drm_boxes_fullDomain/drmdisplacements.0 --box 100 100 50 1 --dt 0.01\nfrom pathlib import Path\nimport numpy as np\nfrom script.read_drm import box_layout\n\nlayout = box_layout(100, 100, 50, 1)\nnpoints = sum(rows * cols for _, rows, cols, _ in layout)\npath = Path(\"outputfiles/drm_boxes_fullDomain/drmdisplacements.0\")\nmetadata = np.loadtxt(path.with_name(\"drmcoordinates.0\"), comments=\"#\", ndmin=2)\nif metadata.shape != (npoints, 8):\n    raise ValueError(\"coordinate/material table does not match box dimensions\")\nxyz = metadata[:, :3]\nvs, vp, density, qs, qp = metadata[:, 3:].T\nnt, remainder = divmod(path.stat().st_size, npoints * 3 * 8)\nif remainder or nt == 0:\n    raise ValueError(\"incomplete snapshot or wrong box dimensions\")\nu = np.memmap(path, dtype=np.float64, mode=\"r\", shape=(nt, npoints, 3))\nfaces = {name: u[:, offset:offset + rows*cols, :].reshape(nt, rows, cols, 3)\n         for name, rows, cols, offset in layout}\nface_xyz = {name: xyz[offset:offset + rows*cols].reshape(rows, cols, 3)\n            for name, rows, cols, offset in layout}\nprint(face_xyz[\"exterior_front\"].shape) # (101, 51, 3)\nprint(faces[\"exterior_front\"].shape)  # (nt, 101, 51, 3)\nprint(faces[\"interior_bottom\"].shape) # (nt, 99, 99, 3)\ntime = np.arange(nt) * (10 * 0.001)  # use your print rate and solver dt",
        "note": "For a restart, the first saved step is the first multiple of the print rate at or after the restart step. Add first_saved_step * simulation_delta_time_sec to time. The reader dimensions are exterior dimensions; center offsets and strike affect positions but not face shapes or the global XYZ component convention."
      },
      {
        "id": "drm-flat-output",
        "title": "DRM box execution details",
        "paragraphs": [
          "The recorder in quake/forward/io_drm.cc uses the existing plane GPU interpolation kernel and collective MPI-IO on solver ranks. Only compact sampled XYZ results are downloaded at requested steps, with independent file views and buffers per box. It introduces no field exchanges and is independent of the plane output backend and plane I/O pool. drm_output_batch_size defaults to 16 snapshots per file; set 1 to write every snapshot. Partial final batches flush on normal completion.",
          "Missing, multiply owned, or scalar-acoustic points stop setup. Poroelastic samples represent skeleton displacement. Close sample spacing does not refine the mesh. The existing three-stage implementation remains in quake/forward/drm.cc, controlled by implement_drm, which_drm_part, drm_directory, drm_print_rate, and its legacy geometry parameters. The new observation files cannot be fed directly to the existing part2 reader. Full CUDA/MPI station-history comparisons remain necessary before calling a physics/backend combination validated.",
          "Host output regressions are documented in benchmark_inputs/drm_flat_output/README.md. Both MPI test executables require OUTPUT_PATH EXPECTED_RANKS (1 or 4), and check MPI_COMM_WORLD before opening output. CI uses matching Open MPI compiler/launcher commands on Ubuntu 24.04, checks that a rejected launch leaves existing output untouched, and compares every byte from one- and four-rank runs. Earlier green jobs without the rank-count guard did not establish multi-rank coverage because independent singleton processes could race on the same files. These tests do not validate CUDA propagation or physics."
        ],
        "note": "Paths resolve from the run working directory; use absolute paths when needed, without ~ expansion. Generated displacement files, coordinate files, and the station binary must use distinct paths. Station and plane outputs remain controlled by their own settings."
      },
      {
        id: "animate-planes",
        title: "Animate horizontal-plane vertical velocity",
        paragraphs: [
          "The script/animations_mpi.py renderer uses manual arguments only; it does not read parameters.in or any other configuration file. Supply --plane-file, --nx, --ny, and --dt for the chosen binary. --dt is the interval between saved snapshots in seconds, equal to simulation_delta_time_sec times output_planes_print_rate; supply that product directly. --dx and --dy set grid spacing in metres (defaults: 500 and 325). Install NumPy, xarray, joblib, and PyGMT with its GMT dependencies in the plotting environment.",
          "This renderer assumes strike=0 (modulo 360) and dip=0; the headerless binary cannot verify the angles. Its C-order array is [time, along-strike point, down-dip point, component]. Rows advance in solver X and columns in solver Y: a 201 by 126 plane with 500 m and 325 m spacing plots X upward over 100 km and Y rightward over 40.625 km. Do not transpose or mirror this matrix. Coordinates are offsets from the plane origin, not latitude/longitude or an inferred hypocenter location.",
          "Vertical velocity is positive downward in cm/s, calculated with central differences inside the record and one-sided differences at the endpoints. The color table saturates outside the selected limits; the velocity arrays and reported extrema remain unclipped. This derivative of sampled displacement need not exactly match velocity computed at the solver time step.",
          "Each worker closes its PyGMT session before deleting its temporary directory and removes the matching exit callback; reused workers start a new session for their next chunk. GMT and Python temporary files, including download locks, are isolated per worker. GMT_DATA_UPDATE_INTERVAL = off disables remote catalogue refresh before PyGMT is imported, since this renderer uses local displacement data and the installed vik colour table. Older script versions could report Not available in classic mode or Module 'end' failed with status code 30 during shutdown, or gmt_data_server.txt.download and Bad record counter errors during catalogue refresh. These diagnostics do not change already saved PNGs; the updated script fixes their causes and still propagates rendering failures.",
        ],
        codeLabel: "Shell — select a completed plane file; dt=0.02 assumes snapshots every 0.02 s",
        code: "python /path/to/Hercules_GPU/script/animations_mpi.py \\\n  --plane-file planes_fullDomain/planedisplacements.0 \\\n  --nx 201 --ny 126 --dx 500 --dy 325 --dt 0.02 --check-only\n\n# Open MPI: two allocated nodes, at least four cores per rank.\n# --workers is per MPI rank. Use a fresh output directory for each run.\nmpirun -np 2 --map-by ppr:1:node:PE=4 --bind-to core --report-bindings \\\n  python /path/to/Hercules_GPU/script/animations_mpi.py \\\n  --plane-file planes_fullDomain/planedisplacements.0 \\\n  --nx 201 --ny 126 --dx 500 --dy 325 --dt 0.02 \\\n  --workers 4 --output-dir Animations_velocity_run1 --velocity-limit 10",
        note: "Relative paths resolve from the current working directory. Each worker retains one memory mapping and has its own temporary GMT directory and color table. Frame names encode time order even when workers finish out of order. --check-only checks complete snapshot sizes and samples the first, middle, and last velocity frames; it does not validate the source location or all snapshots. For restarted output, set --first-time to the time of its first stored snapshot (default: 0 s); restarts cannot be detected from the binary. Input byte order must match the machine reading the native float64 file.",
      },
      {
        id: "compare-plane-outputs",
        title: "Compare plane output backends",
        paragraphs: [
          "Compare every saved XYZ displacement only when geometry, point order, sample times and restart segments match. Plane binaries contain native-endian float64 metres in [snapshot, along-strike point, down-dip point, XYZ] order; the headerless files cannot verify geometry or time metadata.",
          "Check expected file sizes independently: equally truncated outputs do not establish complete output. A snapshot occupies nx*ny*3*8 bytes. A 361-by-361 plane with 1,000 snapshots occupies 3,127,704,000 bytes; a 71-by-71 plane with that count occupies 120,984,000 bytes. Wait for both runs to finish and flush their output.",
          "For finite values, an application-specific tolerance can use abs(candidate-reference) <= atol + rtol*abs(reference), with atol in metres and dimensionless rtol. NaN and infinity always fail, even in byte-identical files. Record complete-history errors and their locations, not only differences between global peaks. Agreement between backends cannot establish physical accuracy or correctness of a shared reconstruction.",
        ],
      },
      {
        id: "plane-amplitude-diagnostics",
        title: "Plane-point metadata and motion interpretation",
        paragraphs: [
          "Use script/read_planes.py to inspect completed plane layout, finite values and vertical peak motion. The largest displacements need not coincide with the largest differences between output backends. Indices are zero-based and components are solver-local X/Y, Z/down in metres. The saved interval is output_planes_print_rate times simulation_delta_time_sec.",
          "HERCULES_PLANE_METADATA=1 enables a setup-only audit after solver_init establishes CPML membership. Each solver MPI rank owning plane points writes plane_point_metadata.RANK.csv in output_planes_directory. The suffix is the rank, not the plane number; each file can contain several planes. Ranks without points write no file. This works with every plane backend and adds no per-step work. The export contains no pressure, velocity or displacement histories.",
          "Rows contain plane,i,j,rank,element,x_m,y_m,z_m,medium,cpml,vp_m_s,vs_m_s,rho_kg_m3,h_m,valid,topostrip,acoustic_reconstruction,pressure_defined. Element IDs are local to the rank. Coordinates include plane rotation and topography/building shifts. Speeds are m/s, density is kg/m³ and h_m is element width in metres. Medium is acoustic, solid, inactive or unknown; solid alone does not identify its constitutive law. CPML membership is 0 outside, 1 inside, or -1 unresolved. The final flags record point validity, tetrahedral topography interpolation, configured acoustic reconstruction and availability of its full hex pressure stencil.",
          "A geometric boundary band is only a CPML candidate region. Join metadata using (plane,i,j), reject duplicate ownership and keep missing ownership unknown. Do not mix rank files from different runs. A short run with the same saved mesh, materials, topography and plane configuration can supply metadata for existing displacement output. Use fresh output directories because plane setup truncates displacement files.",
          "Acoustic plane displacement is reconstructed by integrating -grad(p)/rho every solver step. In velocity-potential mode, p = -(phi[n] - phi[n-1])/dt is formed before taking the pressure gradient. The existing host buffer swap retains the previous gathered field at plane nodes, so this conversion adds no GPU transfers. The startup log prints the reconstruction formula. Inside CPML, plane and station reconstruction now apply the containing element's inverse gradient stretch and a separate memory at each observation point before integrating acceleration. Outside CPML this is an exact identity. The Acoustic output CPML gradient correction log line counts affected stations and plane points. Absorbing-layer values remain artificial CPML fields; use interior observations for physical motion comparisons. Older builds omitted this output memory correction.",
          "Earlier plane output treated velocity potential as pressure, introducing an incorrect sign/time integration and artificial acoustic displacement drift shared by legacy, hybrid, and MPI-IO output. Rebuild and regenerate affected acoustic plane files. The correction matches the existing station formulation and changes output reconstruction rather than wavefield evolution. A host test of a two-second potential pulse with X velocity 0.2*sin(pi*t/2)^2 gives the expected final 0.2 m displacement and zero residual velocity; the old code gives -0.6 m at 4 s and continues drifting. The test also checks all three components against the production station routine and retains the pressure-mode constant-acceleration solution. This is host reconstruction coverage, not validation of a complete CUDA/MPI acoustic-gravity/CPML case.",
        ],
        codeLabel: "Enable the setup-only export before launching Hercules",
        code: "export HERCULES_PLANE_METADATA=1",
      },
      {
        id: "acoustic-station-diagnostics",
        title: "Compare acoustic pressure and motion",
        paragraphs: [
          "Add stations at selected plane locations and nearby submerged points. output_stations takes latitude, longitude and depth below the local topographic surface. Use the model projection and origin, then verify actual XYZ coordinates and acoustic ownership in station_mediums.txt. Retain station IDs and record any inward offset at an exact boundary.",
          "Use ASCII output with pressure, velocity and acceleration enabled. Read its column headers, require complete finite samples and matching time axes, and verify density before comparing pressure and motion. Recheck material ownership after coarsening and CPML membership after changing layer thickness; a formerly acoustic point can become solid or move into CPML.",
          "Save separate cases with the same executable, source, mesh controls, time step, duration, damping, gravity and output settings. Compare the default cpml_clamp_outer=yes with explicit no, changing only that control. CPML remains active in both runs. A separate acoustic_free_surface_gravity=no run changes the physical surface condition and is not a CPML-only comparison.",
          "For flat gravity-surface points outside CPML, the continuum positive-down relation is p = −ρg uZ. The residual uZ + p/(ρg) is a diagnostic, not an acceptance threshold. Lumped surface mass, element-gradient recovery and staggered pressure timing limit this check. Saved sampling and ASCII precision affect integrated-motion checks. Station/plane agreement cannot exclude a shared reconstruction error. Keep raw histories; smaller amplitudes alone do not validate a boundary setting.",
        ],
        codeLabel: "parameters.in — record acoustic station pressure and motion",
        code: "station_binary = no\nprint_station_velocities = yes\nprint_station_accelerations = yes\nprint_station_pressures = yes\noutput_stations_print_rate = 1",
      },
      {
        id: "gravity-mode-diagnostics",
        title: "Gravity-wave dispersion and reconstruction limits",
        paragraphs: [
          "A uniform, lossless water-column reference can check gravity-wave dispersion without CPML, adaptive meshes, solid coupling or GPU/MPI execution. For depth H, horizontal wavenumber k and water sound speed c, a rigid bottom gives omega² = g mu tanh(mu H), with mu² = k² − omega²/c². Frequency accuracy and surface motion recovered from the first element gradient converge differently.",
          "The October 2026 local operator investigation extracted the production Hex8 stiffness and lumped volume/surface masses and verified independent exact element quadrature. For a 250 m wavelength in 175.78125 m depth, h=21.97265625 m gave −1.253% frequency error and surface-gradient motion 19.65% below the pressure-based value. Halving h gave −0.316% and 11.58%. These historical uniform-grid results do not validate an Izmir waveform.",
          "After unit sinusoidal permanent bottom displacement, the equilibrium surface vertical displacement is zero but the first-element gradient can retain a spurious vertical offset. Horizontal particle displacement can remain physically nonzero. Separate local checks of a compact pulse and 100 s slow potential ramp confirmed the output time-integration identities for the recorded implementation. They do not validate the evolved potential field. Investigation scripts and smoke-test inputs are retained locally, not distributed with the repository.",
          "Lowering simulation_wave_max_freq_hz changes mesh and frequency-dependent damping settings; it does not filter the source. Retaining the original source also retains frequencies beyond the reduced mesh resolution. Gravity wavelengths can be much shorter than acoustic wavelengths at the same frequency. After changing resolution or CPML thickness, recheck materials and actual containing-element membership before comparing histories.",
        ],
      },
      {
        id: "plane-animation-colours",
        title: "Smooth the animation colours",
        paragraphs: [
          "The renderer uses a continuous vik colour palette and defaults to --plot-upsample 4 --plot-interpolation linear. It interpolates velocity onto a finer display grid before plotting, reducing visible blocks between native samples. A 201 by 126 plane becomes 801 by 501 display points. --dpi controls exported PNG resolution; increasing DPI alone does not add intermediate grid values.",
          "Set --plot-upsample 1 to display the native grid, or 8 for a finer display. The factor must be a positive integer; display-point count grows approximately with its square, increasing memory and rendering work per worker. --plot-interpolation cubic produces smoother gradients but can overshoot locally; GMT clipping limits values to the input's global range. The default linear interpolation avoids overshoot.",
          "Only the display grid is interpolated. Displacement files, original velocity samples, matrix orientation, and reported extrema and saturation percentages are unchanged. Interpolation does not increase the simulation's physical resolution. Broad solid red or blue areas can instead indicate saturation: increase --velocity-limit to distinguish larger amplitudes, keeping the same limit across frames. --velocity-limit 50 displays the range ±50 cm/s.",
        ],
        codeLabel: "Options to append to the animation command (also the defaults)",
        code: "--plot-upsample 4 --plot-interpolation linear --dpi 300",
      },
      {
        id: "plane-animation-affinity",
        title: "Set Open MPI affinity for animation workers",
        paragraphs: [
          "Check mpirun --version first: these flags are for Open MPI. Intel MPI and MPICH use different affinity controls. On Linux, child processes inherit CPU affinity, so binding a Python rank to one core can leave all its joblib workers competing on that core. Set both mapping and binding explicitly when launching the animation.",
          "For two exclusively allocated nodes with at least 100 available cores each, use one rank per node and --bind-to none as below. ppr:1:node controls rank placement; --bind-to none adds no Open MPI CPU binding. Scheduler-enforced CPU sets still apply. Each rank starts 100 rendering workers, for 200 workers across the two nodes. Workers are not individually pinned. Match -np to the node count and --workers to the cores available to each rank.",
          "For four workers per rank, the preceding example uses --map-by ppr:1:node:PE=4 --bind-to core. The rank and its workers share a four-core group under the default core-based PE interpretation. For two such ranks per node, use ppr:2:node:PE=4, allocate at least eight cores per node, and adjust the total -np accordingly. Mapping and binding do not request resources from the scheduler.",
        ],
        codeLabel: "Open MPI — two exclusively allocated nodes with at least 100 cores each",
        code: "mpirun -np 2 \\\n  --map-by ppr:1:node --bind-to none --report-bindings \\\n  python /path/to/Hercules_GPU/script/animations_mpi.py \\\n  --plane-file planes_fullDomain/planedisplacements.0 \\\n  --nx 201 --ny 126 --dx 500 --dy 325 --dt 0.02 \\\n  --workers 100 --output-dir Animations_velocity_run1",
        note: "Use --report-bindings to inspect each launched rank's binding. These examples describe placement and CPU availability; they do not establish an optimal worker count or validate rendering performance on a particular cluster.",
      },
      {
        id: "plane-animation-affinity-check",
        title: "Check the CPU set visible to Python",
        paragraphs: [
          "On Linux, run the same placement policy with os.sched_getaffinity(0) to report the parent rank's allowed logical CPUs on each host. Substitute the PE=4 policy to check that alternative. Logical CPU counts may exceed physical core counts when simultaneous multithreading is enabled; this does not show individual worker placement.",
        ],
        codeLabel: "Open MPI — Linux CPU-affinity diagnostic",
        code: "mpirun -np 2 --map-by ppr:1:node --bind-to none --report-bindings \\\n  python -c 'import os, socket; cpus = sorted(os.sched_getaffinity(0)); print(socket.gethostname(), \"allowed logical CPUs:\", len(cpus), cpus, flush=True)'",
      },
      {
        id: "checkpoints-write",
        title: "Write rotating checkpoints",
        paragraphs: [
          "A positive checkpointing_rate writes every selected time step and alternates between checkpoint.out0 and checkpoint.out1. Only the two most recent checkpoint states are retained.",
        ],
        codeLabel: "parameters.in",
        code: "use_checkpoint     = 0\ncheckpointing_rate = 5000\ncheckpoint_path    = output/checkpoints",
        note: "CPML rejects checkpoint reads and writes because its histories are absent from the legacy format. Other checkpoint output serializes solver state and activates host-visible compatibility paths. Include its cost when measuring production performance.",
      },
      {
        id: "checkpoints-restart",
        title: "Restart from a checkpoint",
        paragraphs: [
          "Choose one completed rotating file, copy it to checkpoint.in in checkpoint_path, and enable use_checkpoint. The checkpoint header records the MPI rank count; a restart with a different rank count aborts.",
        ],
        codeLabel: "shell and parameters.in",
        code: "# Select a completed checkpoint without modifying the original.\ncp output/checkpoints/checkpoint.out1 output/checkpoints/checkpoint.in\n\n# In parameters.in:\nuse_checkpoint     = 1\ncheckpointing_rate = 5000\ncheckpoint_path    = output/checkpoints",
        note: "Use the same mesh, MPI rank count, precision, and compatible Hercules executable. Keep checkpoint.out0 and checkpoint.out1 until the restarted run has been validated.",
      },
      {
        id: "mesh-output",
        title: "Export and read a selected mesh subdomain",
        paragraphs: [
          "mesh_coordinates_for_matlab writes one pair of binary files per rank that owns selected elements. mesh_coordinates.RANK contains eight records per hexahedral element; mesh_data.RANK contains one material record per element.",
          "The standalone script/read_mesh.py accepts command-line arguments; run the command shown below without its leading comment marker. script/README.md lists the utilities and their dependencies.",
        ],
        codeLabel: "Python — script/read_mesh.py",
        code: "# Runnable version from the repository root:\n# python3 script/read_mesh.py output/mesh/mesh_coordinates.0 output/mesh/mesh_data.0\n\nimport numpy as np\n\ncoord_dtype = np.dtype([\n    (\"geid\", np.int64), (\"x\", np.float64),\n    (\"y\", np.float64), (\"z\", np.float64),\n])\nmaterial_dtype = np.dtype([\n    (\"geid\", np.int64), (\"Vs\", np.float32),\n    (\"Vp\", np.float32), (\"rho\", np.float32),\n])\ncoords = np.fromfile(\"output/mesh/mesh_coordinates.0\", dtype=coord_dtype)\nmaterials = np.fromfile(\"output/mesh/mesh_data.0\", dtype=material_dtype)\nassert coords.size == 8 * materials.size",
        note: "Enable the export with mesh_coordinates_for_matlab=yes, set mesh_coordinates_directory_for_matlab, and define mesh_corners_matlab as min_x min_y max_x max_y min_z max_z in metres.",
      },
      {
        id: "monitor-output",
        title: "Monitor and distributed field output",
        paragraphs: [
          "monitor_file is the first record to preserve with a scientific result: it reports the parsed configuration, execution paths, progress, and timing. output_parallel separately enables distributed mesh-field output and its statistics; use it only when full fields are scientifically necessary because file volume grows rapidly.",
          "Routine Hercules configuration, progress and timing go to stdout; warnings and errors go to stderr. This includes CPML/geometry/plane-output reports, nonlinear and hydrostatic/geostatic initialization, and checkpoint success reports. Stdout is line-buffered for batch jobs. The existing monitor_file messages are retained: routine monitor entries are mirrored to stdout and monitor diagnostics to stderr; direct module messages are not automatically copied to monitor_file. Preserve both scheduler logs and monitor_file. Keep Slurm --output and --error destinations separate, and avoid 2>&1 or MPI stream-merging options when separate logs are needed. Rank-zero guards remain unchanged; ordering across MPI ranks is not guaranteed. Module commands, MPI launchers and external libraries may still print non-error information to stderr.",
          "Some Input Data includes the embedded Hercules Git commit and Build tracked tree state, maximum frequency (Hz), dt (s), the parsed integer simulation_node_per_wavelength, simulation start/end times and configured duration (end minus start) in seconds, region origin latitude/longitude, east/north lengths and deep boundary. It lists source_directory and source_directory_output; active topography adds topographybase_zcoord, type_of_etree, maximum_octant_level, computation_method and topography_directory. A 3-D model adds 3D_velocity_model_directory and 3D_velocity_model_depth_reference, including legacy when omitted; explicit seabed geometry adds its file, water Vp and density. Configured ASCII stations and planes show their output directories with effective defaults; binary stations do not use output_stations_directory. Topography controls and source/station paths are echoed before their modules initialize, retaining their later validation; missing required entries appear as <not set>. This is a configuration record, not proof that a run completed.",
        ],
        codeLabel: "parameters.in",
        code: "monitor_file           = output/monitor.txt\noutput_parallel        = 1\nsimulation_output_rate = 100\noutput_displacement    = 1\noutput_displacement_file = output/displacement\noutput_velocity        = 1\noutput_velocity_file   = output/velocity\noutput_stats_file      = output/output-stats.txt",
        note: "Station and plane output usually provide a much smaller observation set. Estimate expected samples and storage before enabling full distributed fields.",
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
          "For nonlinear/plastic cases, validate yield onset, hysteresis, accumulated plastic deformation, and sensitivity to constitutive substeps and tolerances.",
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
