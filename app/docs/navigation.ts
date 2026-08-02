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
        description: "Elastic, acoustic, poroelastic, damping, and PML.",
      },
      {
        slug: "input-files",
        title: "Input files and directories",
        description: "Source, material, topography, and auxiliary inputs.",
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
