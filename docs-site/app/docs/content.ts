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
          "Use an explicit parallel make job count appropriate for the login or build node. The large host-side solver orchestration in psolve.cc is compiled by the MPI C++ compiler while CUDA kernels and launch wrappers in psolve_gpu_launch.cu are compiled independently by nvcc; make -j 8 can run those jobs concurrently. Incremental builds also avoid rerunning nvcc when only host orchestration changes.",
          "PTX assembly is the expensive phase for nonlinear constitutive CUDA. Sources under quake/forward/nonlinear instantiate separate production kernels in nonlinear_gpu_classic.o, nonlinear_gpu_manzari_dafalias.o, nonlinear_gpu_pdmy03.o, nonlinear_gpu_riva_sand.o, and nonlinear_gpu_riva_sand_intermediate_bias_research.o, plus small runtime and material-point objects. The research successor has a dedicated material-point CUDA object for its 138-value frozen state. The runtime launches only the family selected by material_model. Parallel make can therefore run independent NVCC/PTX-assembly jobs concurrently without relocatable device code, noinline annotations, or equation changes.",
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
          "The CUDA nonlinear constitutive and force path is the default when HERCULES_NONLINEAR_BACKEND is unset and keeps constitutive calculations in double precision. Set HERCULES_NONLINEAR_BACKEND=cpu explicitly for the established CPU/GSL reference models. ManzariDafalias, PressureDependMultiYield03, RIVASand, RIVASandCustom, RIVASand02, and RIVASandCustom02 are device-first and currently require the GPU backend for production-sized runs; their material-point regressions execute the same allocation-free constitutive cores on host and device. The selected backend is reported to standard error at startup.",
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
        title: "RIVASand02 — branch-reversal research successor",
        paragraphs: [
          "RIVASand02 selects the native intermediate-bias research successor. It is a distinct internal material sourced from OpenSees_RIVA commit 637d207e867b863dee7b2bf14d9aeb75fdc32e7b, with the plastic-activity gate from commit 5e4ebd30cb5f43b677c61306c7f1fbbdcb24905c; it is not an alias for, replacement for, or retuning of production RIVASand or RIVASandCustom.",
          "RIVASand02 uses Vs_limit Dr [G0]; RIVASandCustom02 exposes Vs_limit Dr G0 M kd h m zeta e_max e_min Q R n_G. The same effective parameters produce the same response in either name. All active RIVA-Sand names store one interpolated material block, initial void ratio and target shear modulus per element, plus one cached elastic scale per integration point. Custom density bounds are used for both initial and reference density.",
          "The material returns effective skeleton stress. In a coupled BrickUP element, BrickUP remains the sole owner of pore pressure, storage, flow, drainage, and the one-time pore-pressure contribution to total stress. This research identifier should remain confined to verification or explicitly reported research runs until it has completed the project validation process.",
          "Its standard material-point fixture uses a dedicated CUDA object and the same allocation-free constitutive core on host and device. The runner compares the complete 138-value state contract, geostatic-admission and reversal selectors, exact zero-increment preservation, and the independent sealed uniform BrickUP pressure history.",
          "On the branch-reversal research branch, kernel revision 4 adds the loading-branch reversal reference ported from OpenSees_RIVA e97c12fbc. It retains completed-half-cycle plastic activity in both backbones. The new reference prevents a smooth combined strain path from losing reversals when the host step is refined. This changes event handling and remains a research successor, with CPU/GPU and field validation required before acceptance. Automatic stiffness uses the full initialized shear-modulus evaluator after stress-preserving geostatic admission, including pressure guards and initial bias multipliers. It uses approximate K0 stress for approximate_geostatic_state=yes and the solved effective stress at dynamic activation for the ramped no workflow. The elastic scale is then frozen, not recomputed during shaking, and saved in checkpoints. Both 02 names use native checkpoint marker 7 for the changed reversal-reference state, and original names marker 4; older checkpoints must be rerun. Keep the same model, mesh and material inputs on restart. The retired RIVASandIntermediateBiasResearch input is not restored. OpenSees BiasVolume and reversalLatch remain absent from this Hercules port.",
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
        note: "Do not use the production RIVASand identifier to request this successor, and do not treat these research results as production RIVA-Sand results.",
      },
      {
        id: "model-riva-sand-custom02",
        title: "RIVASandCustom02 — configurable intermediate-bias successor",
        paragraphs: [
          "Custom02 uses its own e_max/e_min for both the initial and reference relative densities. In particular, the shakedown reference void ratio e=0.601 must remain the same anchor when those bounds change. This applies to shakedown hardening, shear compliance, dense-state weighting, and the phase-volume transition; it is not a new calibration parameter. The cached reference density is used only when the material has the reference bounds.",
          "The material-specific density correction after 9203e02 is retained. Velocity-calibrated stiffness introduces marker 6 for both 02 names, replacing marker 5 and earlier formats. This marker does not identify the separate OpenSees BiasVolume options.",
          "The column order and parameter meanings match RIVASandCustom, but the numerical kernel is RIVASand02, not the original RIVA-Sand kernel. Positive G0 values and continuous material fields are interpolated linearly by Vs. Set every G0 entry to zero for automatic mesh stiffness; never mix zero and positive G0 rows. Repeat the physical parameters across rows when only Dr should vary. Both 02 interfaces use the same effective-skeleton/BrickUP coupling and fixed no_substeps behavior.",
          "Require finite inputs, 0<=Dr<=1, G0>=0, M,kd,h,Q>0, m,zeta,R>=0, e_max>e_min>=0 containing e=0.601, and 0<=n_G<=1. Vs_limit is in m/s; all twelve following columns are dimensionless. G0=0 requests automatic stiffness, not zero modulus. M is supplied directly, not phi_cs. The non-exposed internal controls remain those of the 249-parameter reference record, SHA-256 3c17e962e64d3a3b29d4797399380dfc85669a5b9fd20b187af359f57aa5b56e. The hash identifies that reference record, not an arbitrary custom row. Neither h_p0 nor OpenSees BiasVolume/reversalLatch flags are part of this Hercules interface.",
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
        ],
        codeLabel: "parameters.in — water over elastic solid",
        code: "include_acoustic = yes\nacoustic_formulation = velocity_potential\nacoustic_free_surface_gravity = no\nacoustic_interface_coupling = yes\nacoustic_interface_projection = point\nacoustic_interface_face_filter = all\nacoustic_interface_max_level = -1\nacoustic_mesh_refinement_factor = 1.0\nacoustic_quality_factor = 0\n\ncvmdb_input_file = profile\nnumber_profile_layers = 2\nprofile_layers =\n# z_top(m) Vp(m/s) Vs(m/s) rho(kg/m3) Qp Qs\n 0.0   1500.0    0.0 1000.0 1.0e6 1.0e6\n 100.0 2600.0 1400.0 2300.0 1.0e6 1.0e6",
        note: "The default pressure-release surface is appropriate for the existing hydroacoustic path. Set acoustic_free_surface_gravity=yes when linear surface-gravity waves are part of the modeled physics; lateral acoustic CPML then includes the matched tangential gravity-surface correction. Use face_mass projection only as an experimental comparison path. A positive acoustic_quality_factor is a numerical high-frequency regularizer, not a physical ocean attenuation law.",
      },
      {
        id: "coupling-poroelastic",
        title: "Poroelastic–elastic example",
        paragraphs: [
          "This sealed, zero-permeability configuration treats saturated solids down to 50 m as linear BrickUP elements. The underlying solid remains elastic and the hydraulic interface is no-flow. For a drained surface, change poroelastic_surface_pressure_bc to drained and retain a positive penalty factor.",
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
          "Enable both formulations and select open_pore only when the acoustic and saturated porous regions physically share an open hydraulic contact. The positive penalty factor scales the weak pressure-continuity projection.",
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
        ],
        codeLabel: "parameters.in — BKT2 with explicit Qp and Qs",
        code: "type_of_damping = bkt2\nsimulation_velocity_profile_freq_hz = 1.0\nuse_parametric_q_factor = no\nuse_infinite_qk = no\nthe_threshold_damping = 0.05\ndo_damping_statistics = 1\n\ncvmdb_input_file = profile\nnumber_profile_layers = 1\nprofile_layers =\n# z_top Vp Vs rho Qp Qs\n 0.0 2000.0 1000.0 2000.0 200.0 100.0",
        note: "simulation_velocity_profile_freq_hz must not exceed simulation_wave_max_freq_hz. Change bkt2 to bkt, bkt3, or bkt3f only when that formulation has been validated for the case.",
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
          "Choose either a positive cpml_thickness_m or a positive cpml_layers count. The baseline below uses the validated conservative kappa=1 native-stiffness path. Startup prints the realized thickness, reference wave speed, critical time-step estimate, and element counts; verify them before relying on absorption.",
        ],
      },
      {
        id: "pml-baseline",
        title: "Baseline CPML example",
        codeLabel: "parameters.in — conservative CPML configuration",
        code: "abc_type = cpml\ncpml_thickness_m = 0.0\ncpml_layers = 10\ncpml_reflection_target = 1.0e-3\ncpml_poly_order = 2\ncpml_kappa_max = 1.0\ncpml_alpha_max_hz = -1.0\ncpml_cfl_safety = 0.9\ncpml_mpml_cross_fraction = 0.1\ncpml_alpha_floor_fraction = 0.1\ncpml_partition_weight = 1.0\ncpml_clamp_outer = no\ncpml_3d_material_extension = yes\ncpml_stiffness_mode = native\ncpml_stiffness_mode_acoustic = native\ncpml_pml_material_damping = none",
        note: "A smaller reflection target is not automatically better: it increases attenuation strength and can tighten stability. Prefer additional physical thickness before changing kappa from the validated value of 1.",
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
          "When topography reaches a lateral PML, taper it to the flat boundary datum before running. benchmark_inputs/taper_topography.py provides the repository workflow. The material extension controls properties, while the taper controls geometry; both are needed for a clean heterogeneous PML interface.",
        ],
        codeLabel: "parameters.in — heterogeneous acoustic/elastic CPML additions",
        code: "3D_velocity_model = yes\ninclude_acoustic = yes\ninclude_topography = yes\n\nabc_type = cpml\ncpml_3d_material_extension = yes\ncpml_stiffness_mode = native\ncpml_stiffness_mode_acoustic = native\ncpml_kappa_max = 1.0\ncpml_clamp_outer = no",
        note: "Changing cpml_3d_material_extension changes sampled PML material properties. Rebuild or recorrect an imported mesh when switching it; a Phase B import otherwise retains the properties stored in that mesh.",
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
            ["HERCULES_NONLINEAR_BACKEND", "cpu | gpu", "Select nonlinear constitutive and force backend; default: gpu"],
            ["HERCULES_NONLINEAR_STATS_INTERVAL", "0 | positive integer", "Sample global ManzariDafalias GPU substep workload; default: 0 (disabled)"],
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
        title: "Complete parameter reference",
        paragraphs: [
          "The categorized pages below document all 208 unique keys currently consumed from parameters.in by psolve and its active feature modules. Each entry states whether it is required, conditional, optional, or a legacy compatibility alias, together with its type, accepted values or default, and scientific or operational meaning.",
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
        note: "For a 3-D material model, provide qs.in and qp.in as a matched pair with the same number of samples as z.in. If both files are absent, Hercules uses the legacy relations Qs=0.1 Vs and Qp=2 Qs.",
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
        id: "velocity-model-configuration",
        title: "Configure a three-dimensional velocity model",
        paragraphs: [
          "The general 3-D model replaces the background material inside its horizontal footprint. The origin parameters locate the model's lower-left corner geographically; x.in is distance north and y.in is distance east from that origin.",
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
          "Coordinates use the solver convention: x is distance north, y is distance east, and z is depth positive downward, all in metres. x.in and y.in must each contain at least two strictly increasing, uniformly spaced coordinates. The x and y spacings may differ, but the spacing within either axis must remain constant.",
          "Treat y as the outer loop and x as the inner loop. For zero-based indices i and j, profile_id = j * nx + i, where i is the x index and j is the y index. Therefore write every x profile at y0 before advancing to y1: (x0,y0), (x1,y0), ..., (x0,y1), (x1,y1), ....",
          "Within one x-y profile, write z from shallow to deep with no duplicate depths. Different profiles may use different z values and different numbers of layers. Append the matching Vs, Vp, density, and optional Q values in exactly the same order as z.in.",
        ],
        bullets: [
          "The first value in every file is a count, not model data. It must equal the number of values that follow.",
          "x.in has nx coordinates and y.in has ny coordinates; index.in has nx * ny zero-based profile starts.",
          "z.in, vs.in, vp.in, and rho.in must have the same total sample count. If used, qs.in and qp.in must both exist and have that same count.",
          "index.in contains starts only—do not append a final end offset. The end of the last profile is the total sample count from z.in.",
        ],
        table: {
          headers: ["File", "Units", "Contents"],
          rows: [
            ["x.in", "m", "Unique northward coordinates, ascending and uniformly spaced"],
            ["y.in", "m", "Unique eastward coordinates, ascending and uniformly spaced"],
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
          "Hercules then bilinearly interpolates Vs, Vp, density, and optional Qs/Qp across those four columns. With topography enabled, the lookup adjusts the column depths relative to the local interpolated ground surface before selecting layers.",
        ],
      },
      {
        id: "velocity-model-converter",
        title: "Convert a CSV model automatically",
        paragraphs: [
          "The repository includes benchmark_inputs/convert_3d_velocity_model.py, a complete converter based on the script pattern in the original Hercules documentation. It uses only the Python standard library, so Pandas and NumPy are not required.",
          "The CSV header must contain x,y,z,vs,vp,rho. It may also contain both qs and qp; providing only one is an error. The converter sorts axis coordinates and profile depths, writes profiles in the required Y-outer/X-inner order, and rejects incomplete rectangular grids, nonuniform horizontal spacing, duplicate depths, non-finite numbers, and invalid material values.",
        ],
        codeLabel: "shell — run from the Hercules_GPU repository",
        code: "python3 benchmark_inputs/convert_3d_velocity_model.py \\\n  hercules-velocity-model.csv \\\n  --output-dir input/materialfiles\n\n# Expected summary for the worked example:\n# wrote 2 x 2 profiles and 15 samples to input/materialfiles",
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
        codeLabel: "Python 3 — benchmark_inputs/convert_3d_velocity_model.py",
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
        ],
        codeLabel: "Python",
        code: "from pathlib import Path\nimport numpy as np\n\nroot = Path(\"input/materialfiles\")\n\ndef counted(name, dtype=float):\n    raw = np.loadtxt(root / name, dtype=dtype, ndmin=1)\n    n, values = int(raw[0]), raw[1:]\n    if values.size != n:\n        raise ValueError(f\"{name}: header says {n}, found {values.size}\")\n    return values\n\nx = counted(\"x.in\")\ny = counted(\"y.in\")\nz = counted(\"z.in\")\nstarts = counted(\"index.in\", int).astype(int)\nvs, vp, rho = (counted(name) for name in (\"vs.in\", \"vp.in\", \"rho.in\"))\n\nfor name, axis in ((\"x.in\", x), (\"y.in\", y)):\n    if axis.size < 2 or np.any(np.diff(axis) <= 0):\n        raise ValueError(f\"{name}: need at least two increasing values\")\n    if not np.allclose(np.diff(axis), np.diff(axis)[0]):\n        raise ValueError(f\"{name}: spacing is not uniform\")\nif starts.size != x.size * y.size or starts[0] != 0 or np.any(np.diff(starts) <= 0):\n    raise ValueError(\"index.in does not define one ordered profile per x-y column\")\nif not (z.size == vs.size == vp.size == rho.size):\n    raise ValueError(\"z, Vs, Vp, and density counts differ\")\nfor column, (begin, end) in enumerate(zip(starts, np.r_[starts[1:], z.size])):\n    if begin >= end or np.any(np.diff(z[begin:end]) <= 0):\n        raise ValueError(f\"vertical profile {column} is empty or not increasing\")\nif np.any(vs < 0) or np.any(vp <= 0) or np.any(rho <= 0):\n    raise ValueError(\"invalid Vs, Vp, or density value\")\nqs_exists, qp_exists = (root / \"qs.in\").exists(), (root / \"qp.in\").exists()\nif qs_exists != qp_exists:\n    raise ValueError(\"provide both qs.in and qp.in, or neither\")\nif qs_exists and (counted(\"qs.in\").size != z.size or counted(\"qp.in\").size != z.size):\n    raise ValueError(\"Q-file counts must match z.in\")\nprint(f\"valid model: {x.size} x {y.size} columns, {z.size} material samples\")",
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
        note: "computation_method=vt uses the volume-tetrahedra topographic treatment; fem retains whole-cube staircase behavior. Nonlinear topography additionally requires include_nonlinear_analysis=yes.",
      },
      {
        id: "topography-maximum-level",
        title: "Choose maximum_octant_level",
        paragraphs: [
          String.raw`\(\mathtt{maximum\_octant\_level}\) is required only when \(\mathtt{include\_topography=yes}\). Let \(V_{S,\min}^{\mathrm{model}}\) denote the smallest shear-wave velocity in the active solid material model; this is derived from the selected material inputs, not read from a separate parameters.in key.`,
          String.raw`Choose the smallest integer \(N\) that satisfies the exact wavelength inequality below, then set \(\mathtt{maximum\_octant\_level}=N\). The accepted range is 0 through 30. A smaller value prevents the topography mesher from reaching the resolution required by \(\mathtt{simulation\_wave\_max\_freq\_hz}\) and \(\mathtt{simulation\_node\_per\_wavelength}\).`,
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
        note: "This is the shear-wave-controlled baseline. Acoustic refinement and other enabled feature constraints can require a finer mesh; ensure the topography cap does not block the most restrictive active refinement rule.",
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
        ],
        codeLabel: "Python",
        code: "import numpy as np\n\nlength_north = 40_000.0\nlength_east = 40_000.0\nspacing = 100.0\nheight = 800.0\nsigma = 3_000.0\ncenter_x = center_y = 20_000.0\n\nnx = int(round(length_north / spacing)) + 1\nny = int(round(length_east / spacing)) + 1\nif not np.isclose((nx - 1) * spacing, length_north) or not np.isclose((ny - 1) * spacing, length_east):\n    raise ValueError(\"spacing must divide both domain lengths exactly\")\nx = np.arange(nx) * spacing\ny = np.arange(ny) * spacing\nX, Y = np.meshgrid(x, y, indexing=\"ij\")  # array shape: (north, east)\nelevation = height * np.exp(-((X-center_x)**2 + (Y-center_y)**2) / (2*sigma**2))\n\nwith open(\"topography.in\", \"w\") as output:\n    output.write(f\"{spacing:g}\\n\")\n    np.savetxt(output, elevation.reshape(-1, order=\"C\"), fmt=\"%.6f\")\nprint(f\"wrote {nx} x {ny} = {elevation.size} elevations\")",
      },
      {
        id: "topography-validation",
        title: "Validate and inspect topography.in",
        paragraphs: [
          "The check below rejects a spacing/count mismatch, non-finite elevations, and an invalid base coordinate. It also reports perimeter relief, which should be essentially flat when CPML is used for a production model.",
        ],
        codeLabel: "Python",
        code: "import numpy as np\n\nlength_north = 40_000.0\nlength_east = 40_000.0\ntopography_base = 1_250.0\nraw = np.loadtxt(\"input/topography/topography.in\")\nspacing, values = raw[0], raw[1:]\nnx = int(round(length_north / spacing)) + 1\nny = int(round(length_east / spacing)) + 1\nif not np.isclose((nx-1)*spacing, length_north) or not np.isclose((ny-1)*spacing, length_east):\n    raise ValueError(\"topographic spacing does not divide the domain\")\nif values.size != nx * ny:\n    raise ValueError(f\"expected {nx*ny} elevations, found {values.size}\")\nif not np.isfinite(values).all():\n    raise ValueError(\"topography contains NaN or infinite values\")\nelevation = values.reshape(nx, ny)  # [north index, east index]\nif topography_base <= elevation.max():\n    raise ValueError(\"topographybase_zcoord must exceed the highest elevation\")\nperimeter = np.r_[elevation[0,:], elevation[-1,:], elevation[:,0], elevation[:,-1]]\nprint(f\"grid={nx}x{ny}, range={elevation.min():.3f}..{elevation.max():.3f} m\")\nprint(f\"perimeter range={np.ptp(perimeter):.6g} m\")",
        note: "For CPML runs, taper real relief smoothly to one flat reference within an interior buffer that becomes fully flat before the PML inner edge. The solver's flatten-PML route is a stability backstop, not a substitute for preprocessing the surface.",
      },
      {
        id: "topography-cpml",
        title: "Why topography must be tapered before C-PML",
        paragraphs: [
          "C-PML is a regular, axis-aligned absorbing layer with a flat top. With computation_method=vt, Hercules can force topography-intersecting lateral PML elements back onto a flat cube-and-CPML path. This flatten-PML route prevents an unstable element treatment, but it does not remove the artificial scattering produced when substantial relief is cut off abruptly at the PML inner edge.",
          "For a production model, begin the transition inside the physical domain. Leave the deep-interior elevations unchanged, reduce their deviation from a chosen reference elevation through a smooth raised-cosine ramp, and provide an exactly flat outer collar that covers the complete PML thickness with extra margin.",
          "The repository helper applies z_out = z_ref + factor × (z_in - z_ref). The factor is one in the interior, follows a C1 raised-cosine through the taper ramp, and is exactly zero throughout the outer flat collar.",
        ],
        note: "Surface tapering and cpml_3d_material_extension solve different problems. The taper flattens the PML's top surface; material extension controls the 3-D material or water mask inside the PML volume.",
      },
      {
        id: "topography-taper-command",
        title: "Use benchmark_inputs/taper_topography.py",
        paragraphs: [
          "Run the dependency-free Python 3 helper from the Hercules_GPU repository. It reads the grid spacing from the first line, checks the supplied domain dimensions and point count, and writes another file in the same topography.in format. Create the output directory first because the script does not create parent directories.",
          "This example is configured for the bundled 40 km × 40 km benchmark. Replace every dimensional value with the actual case values. Set the PML thickness from the realized thickness printed by Hercules to standard error at startup rather than assuming the benchmark value applies to another mesh. The PML element counts and optional outer-clamp count are reported on that same startup stream.",
        ],
        codeLabel: "shell — run from the Hercules_GPU repository",
        code: "mkdir -p input/topography_pml\n\npython3 benchmark_inputs/taper_topography.py \\\n  input/topography_raw/topography.in \\\n  --out input/topography_pml/topography.in \\\n  --lx 40000 \\\n  --ly 40000 \\\n  --buffer 6000 \\\n  --flat-margin 1500 \\\n  --pml-thickness 780 \\\n  --z-ref auto",
        note: "The input is kept unchanged. If the output already exists, the script stops; add --force only when replacing that output intentionally. Hercules requires the final filename to be topography.in.",
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
            ["--force", "Permit overwriting the selected output", "Use only for an intentional rerun; it also permits in-place replacement when input and output are identical"],
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
        id: "time-functions",
        title: "Source-time functions",
        paragraphs: [
          "Analytical functions scale from zero according to the selected rise-time or Ricker parameters. Tabulated discrete histories are linearly interpolated at solver time steps.",
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
            ["source_strike_deg", "Required", "Fault strike in degrees clockwise from north."],
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
          "The force model applies F(t) = (Fx, Fy, Fz) s(t) at one point and distributes it to the containing eight-node element with trilinear shape functions. The simulation frame is x north, y east, and z positive downward.",
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
            ["base_force_direction_x/y/z", "Required; nonzero vector", "Propagation/polarization direction in the x-north, y-east, z-down frame; Hercules normalizes the vector."],
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
            ["extended_strike_deg / extended_dip_deg", "Required", "Fault-plane orientation in degrees."],
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
          "Check positive-down z, x-north, and y-east signs before interpreting force or station components.",
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
            ["Stations", "station.N; station_mediums.txt", "ASCII time histories and station ownership/material metadata"],
            ["Planes", "planedisplacements.N", "Raw float64 displacement snapshots on a regular plane grid"],
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
          "output_stations contains latitude, longitude, and depth in that order. Hercules writes one ASCII station.N file for each row. The sampling interval is output_stations_print_rate multiplied by simulation_delta_time_sec.",
          "Displacement is always present. Velocity is appended when velocity or acceleration output is requested; acceleration, pressure, pore-pressure ratio, strain, and stress columns are then appended according to their switches. Strain and stress use XX, YY, ZZ, XY, YZ, XZ order.",
          "For nonlinear BrickUP runs, print_station_pore_pressure_ratios=yes appends the conventional experimental ratio Ru(-) = Δu/σ′v0. Δu is the same excess pore pressure used by the Pressure(Pa) column. With approximate geostatic initialization, σ′v0 is derived at the exact station depth from the integrated total-overburden column minus Biot-weighted hydrostatic pressure. With a ramped geostatic phase, Ru remains nan during the ramp because its numerical reference does not yet exist. On the last geostatic force evaluation, Hercules fixes σ′v0 from the committed vertical skeleton stress at the owning nonlinear element's station-output quadrature point—the same constitutive state used by its nonlinear stress columns and by RIVA-Sand initialization—before the first dynamic output row. The final denominator is recorded as sigma_v0_eff_pa in station_mediums.txt; other undefined ratios are also nan. Hercules rejects the ratio request when approximate_geostatic_state=no and both geostatic times are zero, because that configuration supplies no initial effective-stress reference. The ratio is not clipped, so dilation may produce negative values and values may exceed one.",
        ],
        codeLabel: "parameters.in — two stations",
        code: "number_output_stations                 = 2\noutput_stations_print_rate             = 10\noutput_stations_directory              = output/stations\nstation_mediums_output_file            = output/stations/station_mediums.txt\nprint_station_velocities               = yes\nprint_station_accelerations            = yes\nprint_station_pressures                = no\n# Set yes with nonlinear BrickUP to append Ru(-).\nprint_station_pore_pressure_ratios     = no\nprint_station_strains                  = no\nprint_station_stresses                 = no\n\n# latitude_deg longitude_deg depth_m\noutput_stations =\n  34.0000 -118.5000 0.0\n  34.0100 -118.4900 100.0",
        note: "With topography enabled, station depth is shifted to the local topographic surface. Check station_mediums.txt to confirm the owning medium, MPI rank, element, local coordinates, Vp, Vs, density, and sigma_v0_eff_pa reference. Ru(-) is nan for acoustic, elastic, unsaturated, CPML, surface-zero-stress, or otherwise inadmissible station references.",
      },
      {
        id: "nonlinear-station-output",
        title: "Interpret nonlinear and plastic station output",
        paragraphs: [
          "When include_nonlinear_analysis=yes, Hercules appends nonlinear solid-state columns to station files. These contain paired strain and stress tensor components plus volumetric strain/stress, yield-function state Fs, accumulated inelastic measure Kappa, one-dimensional strain/stress diagnostics Gamma1d and Tao1d, and GGmax.",
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
        ],
        codeLabel: "Python",
        code: "import numpy as np\n\nstation = np.loadtxt(\"output/stations/station.0\", comments=\"#\")\ntime = station[:, 0]\ndisplacement = station[:, 1:4]\nvelocity = station[:, 4:7]\nacceleration = station[:, 7:10]\n\nprint(f\"{len(time)} samples, dt_out={np.diff(time).mean():.6g} s\")\nprint(\"peak |u| =\", np.linalg.norm(displacement, axis=1).max())",
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
        ],
        codeLabel: "Python",
        code: "import numpy as np\n\nn_strike = 101\nn_dip = 81\nprint_rate = 20\nsimulation_dt = 0.001\n\nraw = np.fromfile(\"output/planes/planedisplacements.0\", dtype=np.float64)\nvalues_per_snapshot = n_strike * n_dip * 3\nif raw.size % values_per_snapshot:\n    raise ValueError(\"file size does not match the configured plane grid\")\n\nu = raw.reshape(-1, n_strike, n_dip, 3)\ntime = np.arange(u.shape[0]) * print_rate * simulation_dt\nvertical = u[..., 2]\nprint(u.shape, time[-1], np.abs(vertical).max())",
        note: "reconstruct_acoustic_plane_displacements=yes reconstructs motion at scalar-acoustic plane points; elastic points continue to use the elastic displacement field in the same file.",
      },
      {
        id: "animate-planes",
        title: "Animate horizontal-plane vertical velocity",
        paragraphs: [
          "The repository-level animations_mpi.py reads the exact run parameters, including a separate planes_input_file when configured. It derives the snapshot interval from simulation_delta_time_sec times output_planes_print_rate and checks complete float64 snapshots before rendering. Install NumPy, xarray, joblib, and PyGMT with its GMT dependencies in the plotting environment.",
          "This renderer accepts strike=0 (modulo 360) and dip=0. Its C-order array is [time, along-strike point, down-dip point, component]. Rows advance in solver X and columns in solver Y: a 201 by 126 plane with 500 m and 325 m spacing plots X upward over 100 km and Y rightward over 40.625 km. Do not transpose or mirror this matrix. Coordinates are offsets from the plane origin, not latitude/longitude or an inferred hypocenter location.",
          "Vertical velocity is positive downward in cm/s, calculated with central differences inside the record and one-sided differences at the endpoints. The color table saturates outside the selected limits; the velocity arrays and reported extrema remain unclipped. This derivative of sampled displacement need not exactly match velocity computed at the solver time step.",
        ],
        codeLabel: "Shell — run from the simulation case directory after the solver finishes",
        code: "python /path/to/Hercules_GPU/animations_mpi.py --parameters parameters.in --check-only\n\n# --workers is per MPI rank; match the allocation's CPUs per task.\n# Use a fresh output directory to avoid mixing frames from different runs.\nsrun python /path/to/Hercules_GPU/animations_mpi.py \\\n  --parameters parameters.in --workers 4 \\\n  --output-dir Animations_velocity_run1 --velocity-limit 10\n\n# If the binary was moved, add: --plane-file planes_fullDomain/planedisplacements.0",
        note: "Relative paths resolve from the case working directory, as in Hercules. Each worker retains one memory mapping and has its own temporary GMT directory and color table. Frame names encode time order even when workers finish out of order. --check-only samples the first, middle, and last velocity frames; it does not validate the source location or all snapshots. Restarted output requires --first-time with the time of its first stored snapshot. Input byte order must match the machine reading the native float64 file.",
      },
      {
        id: "checkpoints-write",
        title: "Write rotating checkpoints",
        paragraphs: [
          "A positive checkpointing_rate writes every selected time step and alternates between checkpoint.out0 and checkpoint.out1. Only the two most recent checkpoint states are retained.",
        ],
        codeLabel: "parameters.in",
        code: "use_checkpoint     = 0\ncheckpointing_rate = 5000\ncheckpoint_path    = output/checkpoints",
        note: "Checkpoint output serializes solver state and activates host-visible compatibility paths. Include its cost when measuring production performance.",
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
        ],
        codeLabel: "Python",
        code: "import numpy as np\n\ncoord_dtype = np.dtype([\n    (\"geid\", np.int64), (\"x\", np.float64),\n    (\"y\", np.float64), (\"z\", np.float64),\n])\nmaterial_dtype = np.dtype([\n    (\"geid\", np.int64), (\"Vs\", np.float32),\n    (\"Vp\", np.float32), (\"rho\", np.float32),\n])\ncoords = np.fromfile(\"output/mesh/mesh_coordinates.0\", dtype=coord_dtype)\nmaterials = np.fromfile(\"output/mesh/mesh_data.0\", dtype=material_dtype)\nassert coords.size == 8 * materials.size",
        note: "Enable the export with mesh_coordinates_for_matlab=yes, set mesh_coordinates_directory_for_matlab, and define mesh_corners_matlab as min_x min_y max_x max_y min_z max_z in metres.",
      },
      {
        id: "monitor-output",
        title: "Monitor and distributed field output",
        paragraphs: [
          "monitor_file is the first record to preserve with a scientific result: it reports the parsed configuration, execution paths, progress, and timing. output_parallel separately enables distributed mesh-field output and its statistics; use it only when full fields are scientifically necessary because file volume grows rapidly.",
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
