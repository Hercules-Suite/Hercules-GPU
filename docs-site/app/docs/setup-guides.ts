import type { DocPageContent, DocSection } from "./content";

type SectionGroup = { title: string; ids: string[] };
export type SetupGuide = {
  slug: string;
  title: string;
  description: string;
  source: string;
  parent?: string;
  groups: SectionGroup[];
};
const group = (title: string, ids: string): SectionGroup => ({ title, ids: ids.split(" ") });

// Each existing section has one canonical home. The builder below checks that
// every section from the three original guides is assigned exactly once.
export const setupGuides: SetupGuide[] = [
  {
    slug: "input-files",
    title: "Input files & directories",
    description: "Case layout, file paths, and a directory of input guides.",
    source: "input-files",
    groups: [group("Case essentials", "overview rules layout sources")],
  },
  {
    slug: "coordinates",
    title: "Coordinates & rotations",
    description: "Domain and material axes, geographic coordinates, and output rotation.",
    source: "input-files",
    parent: "input-files",
    groups: [group("Coordinates", "rotated-domain rotated-output-reader")],
  },
  {
    slug: "velocity-models",
    title: "3-D velocity models",
    description: "File formats, profile ordering, depth references, conversion, and validation.",
    source: "input-files",
    parent: "input-files",
    groups: [group("Model inputs", "velocity-model-configuration velocity-model-files"), group("Depth references", "velocity-model-depth-reference velocity-model-topography-depths velocity-model-lookup"), group("Worked example", "velocity-model-example velocity-model-example-axes velocity-model-example-ranges velocity-model-example-arrays"), group("Conversion & validation", "velocity-model-converter velocity-model-converter-script velocity-model-validation")],
  },
  {
    slug: "basin-models",
    title: "Basin models",
    description: "Basin-bottom grids and the material assigned inside them.",
    source: "input-files",
    parent: "input-files",
    groups: [group("Basin inputs", "basin-configuration basin-input-format basin-material-lookup")],
  },
  {
    slug: "topography",
    title: "Topography & mesh",
    description: "Surface grids, base alignment, refinement, and receiver geometry.",
    source: "input-files",
    parent: "input-files",
    groups: [group("Configure the surface", "topography-configuration topography-format topography-base-alignment"), group("Choose the mesh", "topography-maximum-level topography-wavelength-refinement topography-surface-consistency"), group("Prepare & check", "topography-generation topography-validation")],
  },
  {
    slug: "seabed-geometry",
    title: "Seabed & water geometry",
    description: "Wet/dry classification, clipping, tolerances, and solid-profile depths.",
    source: "input-files",
    parent: "input-files",
    groups: [group("Geometry", "explicit-seabed-geometry")],
  },
  {
    slug: "topography-tapering",
    title: "Tapering & geometry repair",
    description: "Prepare topography and seabed together near absorbing boundaries.",
    source: "input-files",
    parent: "input-files",
    groups: [group("Tapering workflow", "topography-cpml topography-taper-command topography-taper-options topography-taper-enable repair-tapered-geometry")],
  },
  {
    slug: "physics",
    title: "Physical models",
    description: "Choose a formulation and find its setup and coupling guide.",
    source: "physics",
    groups: [group("Choose the physics", "models coupling")],
  },
  {
    slug: "nonlinear-solids",
    title: "Nonlinear & hysteretic solids",
    description: "Nonlinear setup, von Mises, hysteretic, and frictional models.",
    source: "physics",
    parent: "physics",
    groups: [group("Setup", "nonlinear-plasticity model-linear"), group("Von Mises models", "model-vonmises-ep model-vonmises-fa model-vonmises-fam model-vonmises-bae model-vonmises-bah"), group("Hysteretic models", "model-vonmises-gqh model-vonmises-mkz model-vonmises-ro"), group("Frictional models", "model-mohr-coulomb model-drucker-prager")],
  },
  {
    slug: "soil-models",
    title: "Pressure-dependent soils",
    description: "Manzari–Dafalias, PDMY03, and the RIVA-Sand model family.",
    source: "physics",
    parent: "physics",
    groups: [group("Constitutive models", "model-manzari-dafalias model-pdmy03 model-riva-sand"), group("RIVA-Sand variants", "model-riva-sand-intermediate-bias-research model-riva-sand-custom02 model-riva-sand-custom")],
  },
  {
    slug: "acoustic-elastic",
    title: "Acoustic–elastic coupling",
    description: "Water/solid interfaces, constraints, acoustic Q, and motion recovery.",
    source: "physics",
    parent: "physics",
    groups: [group("Configure coupling", "coupling-acoustic-elastic acoustic-quality-factor acoustic-interface-method"), group("Receiver output", "acoustic-station-gradient")],
  },
  {
    slug: "poroelastic-coupling",
    title: "Poroelastic coupling",
    description: "Elastic, nonlinear, and acoustic interfaces with saturated materials.",
    source: "physics",
    parent: "physics",
    groups: [group("Solid & pore-pressure coupling", "coupling-poroelastic coupling-nonlinear-poroelastic coupling-riva-sand-poroelastic"), group("Water & sediment", "coupling-acoustic-poroelastic")],
  },
  {
    slug: "attenuation",
    title: "Attenuation & damping",
    description: "Rayleigh and BKT damping, material Q, and GPU execution.",
    source: "physics",
    parent: "physics",
    groups: [group("Choose damping", "attenuation damping-none-rayleigh-mass"), group("BKT-family setup", "damping-bkt-material-q damping-bkt-parametric-q damping-bkt2-gpu")],
  },
  {
    slug: "absorbing-boundaries",
    title: "CPML & absorbing boundaries",
    description: "Boundary setup, partitioning, operator details, and validation limits.",
    source: "physics",
    parent: "physics",
    groups: [group("Boundary setup", "acoustic-dashpot pml pml-baseline pml-material-damping pml-heterogeneous-topography"), group("Implementation & checks", "pml-partition-memory pml-complete-operator pml-interface-diagnostics pml-horizon-validation")],
  },
  {
    slug: "source-models",
    title: "Source models",
    description: "Choose a source, configure its time function, and check its inputs.",
    source: "source-models",
    groups: [group("Source essentials", "source-directory supported-models common-parameters"), group("Time function & checks", "source-filter-sampling time-functions validation")],
  },
  {
    slug: "point-force-sources",
    title: "Point & force sources",
    description: "Moment tensors, point forces, force histories, and compliant-base forcing.",
    source: "source-models",
    parent: "source-models",
    groups: [group("Source types", "point-source force-source force-history base-force")],
  },
  {
    slug: "finite-fault-sources",
    title: "Plane & kinked faults",
    description: "Finite-fault geometry and slip/rake matrices.",
    source: "source-models",
    parent: "source-models",
    groups: [group("Finite-fault inputs", "finite-plane plane-matrices kinked-plane")],
  },
  {
    slug: "srfh-sources",
    title: "SRFH sources & depths",
    description: "Finite-fault files and surface versus sea-level source depths.",
    source: "source-models",
    parent: "source-models",
    groups: [group("SRFH inputs", "srfh srfh-depth-reference srfh-files")],
  },
];

const shortTitles: Record<string, string> = {
  "velocity-model-configuration": "Configuration", "velocity-model-files": "Files & ordering",
  "velocity-model-example": "Example overview", "velocity-model-example-axes": "Axes & index.in",
  "velocity-model-example-ranges": "Profile index ranges", "velocity-model-example-arrays": "Material arrays",
  "velocity-model-lookup": "Column sampling", "velocity-model-depth-reference": "Depth reference options", "velocity-model-topography-depths": "Legacy: full & sqd",
  "velocity-model-converter": "Convert CSV inputs", "velocity-model-converter-script": "Converter script", "velocity-model-validation": "Validate model files",
  "topography-configuration": "Configuration", "topography-format": "topography.in format", "topography-base-alignment": "Base & water-surface alignment",
  "topography-maximum-level": "Maximum octant level", "topography-wavelength-refinement": "Wavelength refinement", "topography-surface-consistency": "Surface & receivers",
  "topography-generation": "Generate a grid", "topography-validation": "Validate the surface",
  "topography-cpml": "Why taper the surface?", "topography-taper-command": "Run taper_topography.py", "topography-taper-options": "Taper dimensions", "topography-taper-enable": "Enable & verify", "repair-tapered-geometry": "Repair paired seabed geometry",
  "model-riva-sand-intermediate-bias-research": "RIVASand02", "model-riva-sand-custom02": "RIVASandCustom02", "model-riva-sand-custom": "Custom material interface",
  "coupling-acoustic-elastic": "Setup example", "acoustic-quality-factor": "Acoustic quality factor", "acoustic-interface-method": "Interface methods & constraints", "acoustic-station-gradient": "Recover station motion",
  "pml": "CPML setup rules", "pml-baseline": "Default example", "pml-material-damping": "Material damping in CPML", "pml-heterogeneous-topography": "3-D materials & topography", "pml-partition-memory": "Partitioning & memory", "pml-complete-operator": "Second-order operator", "pml-interface-diagnostics": "Interface diagnostics", "pml-horizon-validation": "Validation scope",
  "srfh-depth-reference": "Source depth reference", "srfh-files": "Auxiliary-file example",
};

function splitSeabed(section: DocSection): DocSection[] {
  const p = section.paragraphs ?? [], b = section.bullets ?? [];
  if (p.length !== 7 || b.length !== 7) throw new Error("Update the seabed guide partition after changing its source content");
  return [
    { id: section.id, title: "Water coverage and wet/dry classification", tocTitle: "Wet/dry classification", tocGroup: "Geometry", paragraphs: [p[0]] },
    { id: "seabed-topography-modes", title: "With and without topography", tocGroup: "Geometry", paragraphs: p.slice(1, 3) },
    { id: "seabed-height-tolerance", title: "Height tolerance", tocGroup: "Geometry", paragraphs: [p[3]] },
    { id: "seabed-clipping", title: "Above-surface policy: error or clip", tocTitle: "Error or clip", tocGroup: "Geometry", paragraphs: [p[4]] },
    { id: "seabed-file-format", title: "File format and flat-ocean example", tocTitle: "File format & example", tocGroup: "Geometry", paragraphs: p.slice(5), code: section.code, codeLabel: section.codeLabel },
    { id: "seabed-material-depths", title: "Solid profiles and depth references", tocTitle: "Profiles & depth references", tocGroup: "Materials & setup", bullets: b.slice(0, 3) },
    { id: "seabed-coupled-physics", title: "Initialization and acoustic surfaces", tocTitle: "Initialization & free surface", tocGroup: "Materials & setup", bullets: b.slice(3, 5) },
    { id: "seabed-rebuild", title: "Rebuilding and supported configurations", tocTitle: "Rebuild & supported modes", tocGroup: "Materials & setup", bullets: b.slice(5) },
    { id: "seabed-validation", title: "Validation scope and limits", tocTitle: "Validation scope", tocGroup: "Validation", note: section.note },
  ];
}

export function organizeSetupGuides(authored: Record<string, DocPageContent>): Record<string, DocPageContent> {
  const pages = { ...authored };
  const assigned = new Set<string>();
  for (const guide of setupGuides) {
    const sections = guide.groups.flatMap(({ title, ids }) => ids.map((id) => {
      const key = `${guide.source}#${id}`;
      if (assigned.has(key)) throw new Error(`Section assigned twice: ${key}`);
      const section = authored[guide.source]?.sections.find((item) => item.id === id);
      if (!section) throw new Error(`Missing section: ${key}`);
      assigned.add(key);
      return { ...section, tocGroup: title, tocTitle: shortTitles[id] };
    }));
    const children = setupGuides.filter((item) => item.parent === guide.slug);
    pages[guide.slug] = {
      eyebrow: "Setting up simulations", title: guide.title, summary: guide.description,
      parent: guide.parent ?? "simulation-setup",
      sections: [
        ...(children.length ? [{ id: "topic-guides", title: "Explore the guides", tocGroup: "Start here", links: children.map((item) => ({ href: `/docs/${item.slug}`, title: item.title, description: item.description })) }] : []),
        ...(guide.slug === "seabed-geometry" ? splitSeabed(sections[0]) : sections),
      ],
    };
  }
  for (const source of ["input-files", "physics", "source-models"]) {
    for (const section of authored[source].sections) {
      if (!assigned.has(`${source}#${section.id}`)) throw new Error(`Unassigned section: ${source}#${section.id}`);
    }
    pages[source].legacySections = setupGuides.filter((guide) => guide.source === source && guide.slug !== source).flatMap((guide) => guide.groups.flatMap(({ ids }) => ids.map((id) => ({
      id, title: authored[source].sections.find((section) => section.id === id)!.title,
      href: `/docs/${guide.slug}#${id}`,
    }))));
  }
  pages["simulation-setup"] = {
    eyebrow: "Setting up simulations", title: "Set up a simulation",
    summary: "Start with the case structure, then choose the material, geometry, physics, and source inputs your simulation needs.",
    sections: [
      { id: "case-files", title: "1. Organize your case", links: [
        { href: "/docs/parameter-file", title: "Parameter file", description: "Syntax, required inputs, and a minimal configuration." },
        { href: "/docs/input-files", title: "Input files & directories", description: "Case layout and the files referenced by parameters.in." },
      ] },
      { id: "material-geometry", title: "2. Define materials and geometry", links: setupGuides.filter((guide) => guide.parent === "input-files").map((guide) => ({ href: `/docs/${guide.slug}`, title: guide.title, description: guide.description })) },
      { id: "physics-sources", title: "3. Choose physics and sources", links: [
        { href: "/docs/physics", title: "Physical models", description: "Constitutive models, coupling, damping, and absorbing boundaries." },
        { href: "/docs/source-models", title: "Source models", description: "Source types, time functions, finite faults, and depth references." },
      ] },
      { id: "run-check", title: "4. Run and check", links: [
        { href: "/docs/running", title: "Launch the simulation", description: "Mesh generation, unified execution, and two-phase runs." },
        { href: "/docs/validation", title: "Validation and benchmarks", description: "Check the numerical and physical behavior of your configuration." },
      ] },
    ],
  };
  return pages;
}
