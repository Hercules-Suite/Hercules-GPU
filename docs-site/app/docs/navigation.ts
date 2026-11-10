import { setupGuides } from "./setup-guides";

export type NavItem = {
  children?: NavItem[];
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
      { slug: "simulation-setup", title: "Setup contents", description: "A step-by-step directory of simulation setup guides." },
      { slug: "parameter-file", title: "Parameter file", description: "Syntax and organization of the single parameters.in file." },
      ...setupGuides.filter((guide) => !guide.parent).map((guide) => ({
        slug: guide.slug, title: guide.title, description: guide.description,
        children: setupGuides.filter((child) => child.parent === guide.slug).map((child) => ({ slug: child.slug, title: child.title, description: child.description })),
      })),
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
        description: "Structured DRM output boxes, DRM stages, incident waves, and half-space.",
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

export const allNavItems = navGroups.flatMap((group) => group.items.flatMap((item) => [item, ...(item.children ?? [])]));

export function navItemForSlug(slug: string) {
  return allNavItems.find((item) => item.slug === slug);
}
