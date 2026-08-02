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
    label: "Start here",
    items: [
      {
        slug: "getting-started",
        title: "Quick start",
        description: "Build Hercules and run a first simulation.",
      },
      {
        slug: "installation",
        title: "Installation",
        description: "Compilers, MPI, CUDA, PROJ, and build settings.",
      },
      {
        slug: "architecture",
        title: "How Hercules works",
        description: "Octor, Quake, adaptive meshes, and the solver lifecycle.",
      },
    ],
  },
  {
    label: "Run simulations",
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
        slug: "physics",
        title: "Physics models",
        description: "Elastic, acoustic, poroelastic, damping, and PML.",
      },
      {
        slug: "configuration",
        title: "Configuration",
        description: "Runtime controls and performance-oriented defaults.",
      },
    ],
  },
  {
    label: "Reference",
    items: [
      {
        slug: "input-files",
        title: "Input files",
        description: "Numerical, physics, source, and model inputs.",
      },
      {
        slug: "output-files",
        title: "Output files",
        description: "Stations, planes, checkpoints, and meshes.",
      },
      {
        slug: "validation",
        title: "Validation",
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
