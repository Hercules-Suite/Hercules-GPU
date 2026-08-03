export type NavItem = {
  slug: string;
  title: string;
  description: string;
};

export type NavGroup = {
  label: string;
  items: NavItem[];
};

export const navGroups: NavGroup[] = [
  {
    label: "Introduction",
    items: [
      {
        slug: "architecture",
        title: "Overview",
        description: "Purpose, numerical method, components, and simulation lifecycle.",
      },
      {
        slug: "installation",
        title: "Installation",
        description: "Compilers, MPI, CUDA, PROJ, and build settings.",
      },
      {
        slug: "getting-started",
        title: "A first simulation",
        description: "Build Hercules and run a small, complete example.",
      },
    ],
  },
  {
    label: "Running Hercules",
    items: [
      {
        slug: "running",
        title: "Running Hercules",
        description: "Unified and two-phase launch workflows.",
      },
      {
        slug: "gpu-execution",
        title: "GPU & MPI execution",
        description: "CUDA-aware MPI, resident fields, and scalable I/O.",
      },
      {
        slug: "configuration",
        title: "Runtime configuration",
        description: "Runtime controls and performance-oriented defaults.",
      },
    ],
  },
  {
    label: "Setting up simulations",
    items: [
      {
        slug: "parameter-file",
        title: "Parameter file",
        description: "Syntax and organization of the single parameters.in file.",
      },
      {
        slug: "physics",
        title: "Physical models",
        description: "Elastic, nonlinear plasticity, acoustic, poroelastic, damping, and PML.",
      },
      {
        slug: "input-files",
        title: "Input files and directories",
        description: "Source, material, topography, and auxiliary inputs.",
      },
      {
        slug: "source-models",
        title: "Source models",
        description: "Point, force, finite-fault, kinked-fault, and SRFH sources.",
      },
    ],
  },
  {
    label: "Parameter reference",
    items: [
      {
        slug: "parameters-domain-time",
        title: "Domain, time & resolution",
        description: "Coordinates, dimensions, duration, time step, and wavelength resolution.",
      },
      {
        slug: "parameters-material-mesh",
        title: "Material, source & mesh",
        description: "Profiles, regional models, source paths, and mesh modes.",
      },
      {
        slug: "parameters-damping-cpml",
        title: "Damping & CPML",
        description: "Attenuation, damping limits, dashpots, and CPML controls.",
      },
      {
        slug: "parameters-coupled-physics",
        title: "Acoustic & poroelastic",
        description: "Scalar acoustics, saturated media, groundwater, and interfaces.",
      },
      {
        slug: "parameters-output",
        title: "Output & diagnostics",
        description: "Stations, planes, checkpoints, field output, and diagnostics.",
      },
      {
        slug: "parameters-topography-buildings",
        title: "Topography & buildings",
        description: "Surface geometry, structures, fixed-base input, and push-downs.",
      },
      {
        slug: "parameters-nonlinear",
        title: "Nonlinear & plastic analysis",
        description: "Geostatic state, elastoplastic models, plasticity controls, and material tables.",
      },
      {
        slug: "parameters-drm-waves",
        title: "DRM & incident waves",
        description: "DRM stages, plane waves, histories, and homogeneous half-space.",
      },
    ],
  },
  {
    label: "Output and validation",
    items: [
      {
        slug: "output-files",
        title: "Output files",
        description: "Stations, planes, checkpoints, meshes, and monitor logs.",
      },
      {
        slug: "validation",
        title: "Validation and benchmarks",
        description: "Benchmarks and correctness checks.",
      },
      {
        slug: "troubleshooting",
        title: "Troubleshooting",
        description: "Common launch, GPU, MPI, and I/O problems.",
      },
    ],
  },
];

export const allNavItems = navGroups.flatMap((group) => group.items);

export function navItemForSlug(slug: string) {
  return allNavItems.find((item) => item.slug === slug);
}
