import { parameterReferencePages, parameterReferenceSummary } from "./parameter-reference";

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
        codeLabel: "shell",
        code: "make clean\nmake -j 8\n\n# Expected executables\nls octor/octor quake/forward/psolve",
        bullets: [
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
          "The CUDA nonlinear constitutive and force path is the default when HERCULES_NONLINEAR_BACKEND is unset and keeps constitutive calculations in double precision. Set HERCULES_NONLINEAR_BACKEND=cpu explicitly for the established CPU/GSL reference models. ManzariDafalias is device-first and currently requires the GPU backend for production runs; its material-point regression executes the same allocation-free core on host and device. The selected backend is reported to standard error at startup.",
          "When abc_type=cpml, solid elements in the absorbing band remain linear elastic and are excluded from both CPU and GPU nonlinear constitutive maps regardless of their Vs. Nonlinearity is confined to the physical interior because the CPML correction assumes a linear time-invariant native operator.",
          "When include_acoustic=yes, acoustic elements are excluded from both CPU and GPU solid constitutive maps regardless of their Vs or legacy nl_id payload. They remain governed only by the selected acoustic formulation and never acquire plastic-strain or nonlinear solid history.",
          "Treat the GPU path as a separately validated execution backend: compare nonlinear station histories, final yield statistics, and force or displacement norms against a CPU run before using a new material calibration in production.",
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
          "The eleven-fixture suite covers EP, the FA/FAM quartic correctors, BAE and BAH, GQH, MKZ, Ramberg--Osgood, ManzariDafalias Modified Euler, and seeded ManzariDafalias first yield. Existing families compare the CPU GSL-backed polynomial or scalar root solvers with CUDA; the MD fixtures execute one allocation-free implementation on host and device and record alpha, alpha_in, fabric, void ratio, plastic multiplier, absolute yield residual, and adaptive-substep diagnostics.",
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
          "The coupled SSPbrickUP/BrickUP benchmark fills the uniform sealed part of that gap. OpenSees runs a zero-permeability SSPbrickUP -lumped element under undrained triaxial compression; Hercules replays the exact affine strain history and invokes its production BrickUP pressure-force and pressure-update CUDA kernels on a unit element. Driver version 2 uses conventional symmetry supports, distributed face tractions, and identical prescribed axial displacement on all four top nodes; it avoids combining a prescribed DOF with an equalDOF retained node under the Transformation handler. An independent CPU reference checks p-p0=-M(ev-ev0), with M=Kf/n. OpenSees's row-summed pressure capacity and Hercules's lumped nodal storage use the same representation. The comparator reports effective and total stress, pore pressure, r_u, radial confinement, storage residual, and material state. This deliberately isolates the common uniform mode and does not claim validation of Darcy diffusion, drainage/open-pore boundaries, nonuniform pressure, multi-element exchange, or the complete production time integrator.",
          "The fixture constants and strain path are test inputs, not a substitute for calibrated production properties. Off-diagonal inputs are tensorial shear strains, so epsilon_xy equals one half of engineering gamma_xy. A single point validates correctness only; use a many-element case for performance measurement.",
        ],
        codeLabel: "shell",
        code: "cd benchmark_inputs/nonlinear_material_point\n./run_suite.sh ../../quake/forward/psolve\n\nOpenSees opensees_manzari_dafalias_drained_triaxial.tcl opensees_md_drained.csv\npython3 generate_hercules_manzari_dafalias_path.py --opensees-csv opensees_md_drained.csv --output manzari_md_drained.in\nmpirun -np 1 ../../quake/forward/psolve --nonlinear-material-point manzari_md_drained.in\npython3 compare_opensees_manzari_dafalias_canonical.py --mode drained --hercules-csv hercules_md_drained_triaxial.cpu.csv --opensees-csv opensees_md_drained.csv\n\nOpenSees opensees_manzari_dafalias_sspbrickup_undrained_triaxial.tcl opensees_md_sspbrickup.csv\npython3 generate_hercules_manzari_dafalias_path.py --opensees-csv opensees_md_sspbrickup.csv --output manzari_md_brickup.in --output-prefix manzari_md_brickup\nmpirun -np 1 ../../quake/forward/psolve --nonlinear-material-point manzari_md_brickup.in\npython3 compare_opensees_manzari_dafalias_brickup.py --hercules-csv manzari_md_brickup.cpu.csv --opensees-csv opensees_md_sspbrickup.csv",
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
      "Model linear and nonlinear/plastic solids, acoustic fluids, saturated porous media, attenuation, topography, and absorbing boundaries in one adaptive framework.",
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
        title: "Nonlinear and plastic solid response",
        paragraphs: [
          "Set include_nonlinear_analysis=yes to activate history-dependent constitutive updates in eligible solid elements. Hercules evaluates stress against the selected yield or backbone relation, evolves plastic strain and internal variables, and can retain permanent deformation after unloading.",
          "The material_model setting selects the constitutive family. material_plasticity_type selects rate-dependent or rate-independent integration, while no_substeps and the nonlinear tolerances control the constitutive update within each global time step.",
          "The legacy Istanbul nl_id classification controls element eligibility only when Istanbul_velocity_model=yes. Other runs ignore nl_id and select eligible elements from the Vs range covered by material_properties_list.",
          "For non-Istanbul vonMises_baE, material_properties_list supplies Su through alpha_cohesion and the exponential hardening parameters through psi and m; non-Istanbul vonMises_ep takes Su from alpha_cohesion. The built-in nl_id constants for these models remain limited to Istanbul_velocity_model=yes.",
          "ManzariDafalias uses adaptive explicit Modified Euler. backbone_errTol is its local stress-increment error target, no_substeps caps accepted adaptive substeps, and error_tolerance bounds the absolute yield-surface residual after plastic stress correction. Its elastic predictor retains the one-sided admissibility condition f_trial <= error_tolerance. A trial increment crossing from strictly inside to outside is split at its elastic yield intersection, and only the remainder enters Modified Euler. Elastic stages preserve alpha, fabric, and plastic strain; these history variables evolve only under plastic flow. Plastic correction uses at most min(no_substeps, 50) iterations; if the absolute criterion remains unsatisfied, Hercules reports constitutive failure code 3 and aborts rather than continuing from the off-surface state. With approximate_geostatic_state=yes, the approximate K0 stress initializes its effective-stress state and centers the small kinematic yield surface on the pre-existing stress ratio, without adding the static offset to dynamic force assembly.",
        ],
        table: {
          headers: ["Family", "Available material_model values", "Role"],
          rows: [
            ["Elastic diagnostic", "linear", "Compute nonlinear diagnostic quantities while keeping elements elastic for comparison"],
            ["von Mises / cyclic plasticity", "vonMises_ep, vonMises_fa, vonMises_faM", "Elastoplastic yielding with the selected hardening formulation"],
            ["Nonlinear backbone models", "vonMises_baE, vonMises_baH, vonMises_GQH, vonMises_MKZ, vonMises_RO", "Nonlinear and hysteretic soil response using the selected backbone law"],
            ["Pressure-dependent plasticity", "MohrCoulomb, DruckerPrager", "Frictional/cohesive yielding with optional tension cutoff"],
            ["Critical-state sand plasticity", "ManzariDafalias", "Pressure-dependent effective-stress response with dilatancy, fabric evolution, and cyclic mobility state"],
          ],
        },
        codeLabel: "parameters.in — feature selection (additional nonlinear controls are required)",
        code: "include_nonlinear_analysis = yes\nmaterial_model = vonMises_ep\nmaterial_plasticity_type = rate_independent",
        note: "A plastic production case also requires calibrated geostatic, tolerance, substep, damping, and model-specific material_properties_list values. Most models use 16 columns; ManzariDafalias uses 17. Validate the selected constitutive model against an appropriate reference before a large run.",
      },
      {
        id: "nonlinear-material-properties",
        title: "Non-Istanbul BAE property table",
        paragraphs: [
          "Rows are ordered by increasing Vs_limit. For non-Istanbul models, the first and last limits define the inclusive nonlinear range and Hercules linearly interpolates active properties for intermediate Vs. BAE reads Su from alpha_cohesion and its exponential hardening parameters from psi and m.",
          "The example brackets a homogeneous Vs=1000 m/s profile. Repeated values make its active BAE properties constant across the selected range; replace them with calibrated strength and modulus-reduction parameters for scientific use.",
        ],
        codeLabel: "parameters.in — synthetic syntax example, not a calibration",
        code: "Istanbul_velocity_model = no\nmaterial_model = vonMises_baE\nmaterial_plasticity_type = rate_independent\nmaterial_properties_count = 2\nmaterial_properties_list =\n# Vs   Su(Pa) phi rate sensitivity hardening dilatancy gamma0 psi m   theta1 theta2 theta3 theta4 theta5 phi_RO\n  900  7.5e4  0   0    0           0         0          0      0.2 1.2 0      0      0      0      0      0\n 1100  7.5e4  0   0    0           0         0          0      0.2 1.2 0      0      0      0      0      0",
        note: "Su, psi, and m are case-specific constitutive inputs and cannot be inferred from Vs alone. Confirm the selected nonlinear-element count in the setup output before timing or production runs.",
      },
      {
        id: "coupling",
        title: "Coupled interfaces",
        paragraphs: [
          "Material classification determines which element formulation advances each region. Interface terms transfer the physically appropriate traction, motion, and pressure information between neighboring acoustic, solid, and poroelastic domains; eligible solid elements may use linear or nonlinear/plastic constitutive response.",
          "For a poroelastic ManzariDafalias region, the constitutive model returns effective skeleton stress and the u-p pressure operator supplies Biot pore-pressure traction separately. Hercules therefore disables the zero-permeability undrained fast path for this model instead of folding alpha²M into a linear bulk stiffness.",
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
            ["HERCULES_NONLINEAR_BACKEND", "cpu | gpu", "Select nonlinear constitutive and force backend; default: gpu"],
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
          "The categorized pages below document all 206 unique keys currently consumed from parameters.in by psolve and its active feature modules. Each entry states whether it is required, conditional, optional, or a legacy compatibility alias, together with its type, accepted values or default, and scientific or operational meaning.",
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
          "Only the inputs enabled by parameters.in are required. Keeping each simulation self-contained makes it easier to move between systems and archive with a publication.",
        ],
        codeLabel: "directory tree",
        code: "case/\n├── parameters.in\n├── sourcefiles/\n│   └── source.in\n├── materialfiles/\n│   ├── x.in\n│   ├── y.in\n│   ├── z.in\n│   ├── vs.in\n│   ├── vp.in\n│   ├── rho.in\n│   └── index.in\n├── topography/\n│   └── topography.in\n└── output/",
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
          "Every file begins with an integer count followed by exactly that many values. x.in and y.in contain at least two strictly increasing, uniformly spaced coordinates. z.in, vs.in, vp.in, and rho.in contain concatenated vertical profiles; index.in gives the zero-based start of each profile.",
          "Horizontal profile order is y outer and x inner: (x0,y0), (x1,y0), then the remaining x positions before advancing to y1. Within each profile, z is strictly increasing and positive downward. Different horizontal columns may contain different z levels and layer counts.",
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
      },
      {
        id: "velocity-model-example",
        title: "Minimal 2 × 2 velocity-model example",
        paragraphs: [
          "This model has four horizontal columns and two depth levels per column. The separators identify different files and are not part of the numeric file contents.",
        ],
        codeLabel: "input/materialfiles",
        code: "# x.in\n2\n0.0\n1000.0\n\n# y.in\n2\n0.0\n1000.0\n\n# index.in — profiles: (x0,y0), (x1,y0), (x0,y1), (x1,y1)\n4\n0\n2\n4\n6\n\n# z.in\n8\n0.0\n1000.0\n0.0\n1000.0\n0.0\n1000.0\n0.0\n1000.0\n\n# vs.in (m/s)\n8\n500.0\n800.0\n550.0\n850.0\n600.0\n900.0\n650.0\n950.0\n\n# vp.in (m/s)\n8\n1000.0\n1600.0\n1100.0\n1700.0\n1200.0\n1800.0\n1300.0\n1900.0\n\n# rho.in (tonne/m3)\n8\n2.00\n2.20\n2.02\n2.22\n2.04\n2.24\n2.06\n2.26",
        note: "If qs.in and qp.in are absent, Hercules uses the legacy relations Qs=0.1 Vs and Qp=2 Qs. Regenerate or re-correct a persisted mesh after changing any material or quality-factor file.",
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
          "The block below is illustrative for a 40 km square domain with relief below 800 m. Recompute the refinement level from the case's minimum element size, and choose a positive base coordinate that safely contains the full elevation range.",
        ],
        codeLabel: "parameters.in — illustrative values",
        code: "include_topography             = yes\ntopography_directory           = input/topography\ncomputation_method             = vt\nmaximum_octant_level           = 9\ntopographybase_zcoord          = 1250.0\ntype_of_etree                  = full\nmin_enclosed_volume_perc       = 0.1\nconsider_nonlinear_topography  = no\ntopography_bkt_gpu_stabilize   = auto",
        note: "computation_method=vt uses the volume-tetrahedra topographic treatment; fem retains whole-cube staircase behavior. Nonlinear topography additionally requires include_nonlinear_analysis=yes.",
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
          "This example is configured for the bundled 40 km × 40 km benchmark. Replace every dimensional value with the actual case values. Set the PML thickness from the realized thickness printed by Hercules at startup rather than assuming the benchmark value applies to another mesh.",
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
          "Displacement is always present. Velocity is appended when velocity or acceleration output is requested; acceleration, pressure, strain, and stress columns are then appended according to their switches. Strain and stress use XX, YY, ZZ, XY, YZ, XZ order.",
        ],
        codeLabel: "parameters.in — two stations",
        code: "number_output_stations       = 2\noutput_stations_print_rate   = 10\noutput_stations_directory    = output/stations\nstation_mediums_output_file  = output/stations/station_mediums.txt\nprint_station_velocities     = yes\nprint_station_accelerations  = yes\nprint_station_pressures      = no\nprint_station_strains        = no\nprint_station_stresses       = no\n\n# latitude_deg longitude_deg depth_m\noutput_stations =\n  34.0000 -118.5000 0.0\n  34.0100 -118.4900 100.0",
        note: "With topography enabled, station depth is shifted to the local topographic surface. Check station_mediums.txt to confirm the owning medium, MPI rank, element, local coordinates, Vp, Vs, and density.",
      },
      {
        id: "nonlinear-station-output",
        title: "Interpret nonlinear and plastic station output",
        paragraphs: [
          "When include_nonlinear_analysis=yes, Hercules appends nonlinear solid-state columns to station files. These contain paired strain and stress tensor components plus volumetric strain/stress, yield-function state Fs, accumulated inelastic measure Kappa, one-dimensional strain/stress diagnostics Gamma1d and Tao1d, and GGmax.",
          "Use these histories to identify yield onset, hysteresis, stiffness reduction, and permanent/plastic response. Read the station header instead of assuming fixed column offsets because velocity, acceleration, pressure, linear strain/stress, and nonlinear columns are controlled independently.",
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
