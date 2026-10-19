import assert from "node:assert/strict";
import { existsSync, mkdtempSync, readdirSync, readFileSync, rmSync, statSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { spawnSync } from "node:child_process";
import { fileURLToPath } from "node:url";
import test from "node:test";

const staticExport = process.env.DOCS_STATIC_EXPORT === "1";
const basePath = staticExport ? (process.env.NEXT_PUBLIC_DOCS_BASE_PATH ?? "/Hercules_GPU") : "";
const exportDirectory = fileURLToPath(new URL("../out/", import.meta.url));

function filesUnder(directory) {
  return readdirSync(directory, { withFileTypes: true }).flatMap((entry) => {
    const path = join(directory, entry.name);
    return entry.isDirectory() ? filesUnder(path) : [path];
  });
}

async function render(path = "/") {
  if (staticExport) {
    const file = join(exportDirectory, path, "index.html");
    return existsSync(file)
      ? new Response(readFileSync(file), { headers: { "content-type": "text/html" } })
      : new Response("Not found", { status: 404 });
  }
  const workerUrl = new URL("../dist/server/index.js", import.meta.url);
  workerUrl.searchParams.set("test", `${process.pid}-${Date.now()}`);
  const { default: worker } = await import(workerUrl.href);

  return worker.fetch(
    new Request(`http://localhost${path}`, {
      headers: { accept: "text/html" },
    }),
    {
      ASSETS: { fetch: async () => new Response("Not found", { status: 404 }) },
    },
    { waitUntil() {}, passThroughOnException() {} },
  );
}

test("documents compact scalar MPI schedules and device profiling limits", async () => {
  const html = await (await render("/docs/gpu-execution")).text();
  assert.match(html, /id="scalar-schedules"/);
  assert.match(html, /HERCULES_SCALAR_MPI_SCHEDULE=compact \(default\)/);
  assert.match(html, /HERCULES_SCALAR_MPI_SCHEDULE=full/);
  assert.match(html, /HERCULES_GPU_MPI_PROFILE=1/);
  assert.match(html, /HERCULES_SCALAR_MPI_BATCH=1 \(default\)/);
  assert.match(html, /HERCULES_GPU_MPI_OVERLAP=1 \(default\)/);
  assert.match(html, /HERCULES_PORE_PRESSURE_UPDATE=output \(default\)/);
  assert.match(html, /Checkpoint reads\/writes/);
  assert.match(html, /final two states/);
  assert.match(html, /12 to 8 phases/);
  assert.match(html, /14 to 12 phases/);
  assert.match(html, /CUDA events protect receive-buffer reuse/);
  assert.match(html, /Host staging remains synchronous/);
  assert.match(html, /Ranks without local active elements/);
  assert.match(html, /host and CUDA-aware MPI backends/);
  assert.match(html, /cumulative rank average\/min\/max/);
  assert.match(html, /Target-system CUDA\/MPI propagation and performance comparisons/);
  assert.match(html, /cover only the host helper/);
  const controls = await (await render("/docs/configuration")).text();
  assert.match(controls, /HERCULES_SCALAR_MPI_SCHEDULE/);
  assert.match(controls, /HERCULES_GPU_MPI_PROFILE/);
  assert.match(controls, /HERCULES_SCALAR_MPI_BATCH/);
  assert.match(controls, /HERCULES_GPU_MPI_OVERLAP/);
  assert.match(controls, /HERCULES_PORE_PRESSURE_UPDATE/);
});

test("renders the Hercules documentation home", async () => {
  const response = await render();
  assert.equal(response.status, 200);
  assert.match(response.headers.get("content-type") ?? "", /^text\/html\b/i);

  const html = await response.text();
  assert.match(html, /<title>Hercules Documentation<\/title>/i);
  assert.match(html, /parallel, octree-based finite-element code/);
  assert.match(html, /soil liquefaction/i);
  assert.match(html, /saturated poroelastic media/i);
  assert.match(html, /acoustic domains/i);
  assert.match(html, /acoustic–poroelastic coupled interfaces/i);
  assert.match(html, /Hercules emblem formed from adaptive octree cells and a seismic waveform/);
  assert.match(html, /one.*parameters\.in.*file/i);
  assert.match(html, /Search documentation/);
  assert.match(html, /hercules-octree-mesh-mark\.png\?v=graded-centered-256/);
  assert.ok(html.includes(`src="${basePath}/hercules-octree-mesh-mark.png?v=graded-centered-256"`));
  assert.doesNotMatch(html, /_vinext\/image\?url=.*hercules-octree-mesh-mark/);
  assert.ok(existsSync(join(fileURLToPath(new URL("../public/", import.meta.url)), "hercules-octree-mesh-mark.png")));
  assert.doesNotMatch(html, /codex-preview|react-loading-skeleton/);
});

test("exports complete pages with working repository-prefixed links and assets", { skip: !staticExport }, () => {
  const origin = new URL(process.env.DOCS_SITE_URL ?? "https://zwyll.github.io/Hercules_GPU/").origin;
  const navigation = readFileSync(new URL("../app/docs/navigation.ts", import.meta.url), "utf8");
  const slugs = [...navigation.matchAll(/slug: "([^"]+)"/g)].map((match) => match[1]);
  for (const slug of slugs) {
    assert.ok(existsSync(join(exportDirectory, "docs", slug, "index.html")), `Missing exported page: ${slug}`);
  }
  assert.ok(existsSync(join(exportDirectory, "404.html")), "Static hosting needs a 404 page");

  for (const file of filesUnder(exportDirectory).filter((path) => path.endsWith(".html"))) {
    const html = readFileSync(file, "utf8");
    const documentPath = `${basePath}/${file.slice(exportDirectory.length).replace(/index\.html$/, "")}`;
    for (const tag of html.match(/<(?:a|link|img|script)\b[^>]*>/g) ?? []) {
      const reference = tag.match(/\b(?:href|src)="([^"]+)"/)?.[1]?.replaceAll("&amp;", "&");
      if (!reference || reference.startsWith("#")) continue;
      const url = new URL(reference, `${origin}${documentPath}`);
      if (url.origin !== origin) continue;
      assert.ok(url.pathname.startsWith(`${basePath}/`), `URL escapes the repository path: ${reference} in ${file}`);
      const relativePath = decodeURIComponent(url.pathname.slice(basePath.length + 1));
      let target = join(exportDirectory, relativePath);
      assert.ok(existsSync(target), `Missing link/asset target: ${reference} in ${file}`);
      if (statSync(target).isDirectory()) target = join(target, "index.html");
      assert.ok(existsSync(target), `Missing directory index: ${reference} in ${file}`);
      if (url.hash && target.endsWith(".html")) {
        const fragment = decodeURIComponent(url.hash.slice(1));
        assert.ok(readFileSync(target, "utf8").includes(`id="${fragment}"`), `Missing section: ${reference}`);
      }
    }
  }
  const home = readFileSync(join(exportDirectory, "index.html"), "utf8");
  assert.ok(home.includes(`content="${origin}${basePath}/og.png"`), "Social image must use the published site URL");
  assert.doesNotMatch(home, /localhost:3000|_vinext\/image|_next\/image/);
});

test("builds the octree emblem from square mesh cells and square waveform terminals", () => {
  const logoSource = readFileSync(
    fileURLToPath(new URL("../assets/hercules-octree-mark.svg", import.meta.url)),
    "utf8",
  );
  const cells = [...logoSource.matchAll(
    /<rect data-octree-cell="true"[^>]*x="([0-9.]+)" y="([0-9.]+)" width="([0-9.]+)" height="([0-9.]+)"/g,
  )].map(([, x, y, width, height]) => ({
    x: Number(x),
    y: Number(y),
    width: Number(width),
    height: Number(height),
  }));

  assert.equal(cells.length, 54);
  for (const { width, height } of cells) {
    assert.equal(width, height, `octree cell ${width} × ${height} is not square`);
  }

  for (let index = 0; index < cells.length; index += 1) {
    for (let otherIndex = index + 1; otherIndex < cells.length; otherIndex += 1) {
      const first = cells[index];
      const second = cells[otherIndex];
      const horizontalOverlap = Math.min(first.x + first.width, second.x + second.width)
        - Math.max(first.x, second.x);
      const verticalGap = Math.max(
        first.y - (second.y + second.height),
        second.y - (first.y + first.height),
        0,
      );
      if (horizontalOverlap > 0 && verticalGap <= 8) {
        const gradingRatio = Math.max(first.width, second.width) / Math.min(first.width, second.width);
        assert.ok(gradingRatio <= 2.2, `adjacent octree cells have an abrupt ${gradingRatio}:1 transition`);
      }
    }
  }

  assert.match(logoSource, /x="364" y="287" width="23" height="23"/);
  assert.match(logoSource, /x="391" y="287" width="23" height="23"/);

  const terminals = [...logoSource.matchAll(
    /<rect (?=[^>]*data-wave-terminal="true")[^>]*y="([0-9.]+)"[^>]*width="([0-9.]+)" height="([0-9.]+)"/g,
  )];
  assert.equal(terminals.length, 2);
  for (const [, y, width, height] of terminals) {
    assert.equal(width, height, `waveform terminal ${width} × ${height} is not square`);
    assert.equal(Number(y) + Number(height) / 2, 256, "waveform terminal is not vertically centered");
  }
  assert.match(
    logoSource,
    /<rect data-octree-cell="true" data-wave-terminal="true"/,
    "waveform terminals must be octree mesh cells",
  );
  assert.match(
    logoSource,
    /d="M148 256 .* 364 256"/,
    "waveform baseline must meet both terminal centers at y = 256",
  );
});

test("server-renders a structured documentation page", async () => {
  const response = await render("/docs/getting-started");
  assert.equal(response.status, 200);

  const html = await response.text();
  assert.match(html, /Run your first Hercules simulation/);
  assert.match(html, /Before you begin/);
  assert.match(html, /Current Hercules GPU code path/);
});

test("documents octree-compatible domain sizing and refinement level selection", async () => {
  const architectureResponse = await render("/docs/architecture");
  assert.equal(architectureResponse.status, 200);
  const architecture = await architectureResponse.text();
  assert.match(architecture, /Octree-compatible domain dimensions/);
  assert.match(architecture, /class="mathExpression mathInline"/);
  assert.match(architecture, /class="mathExpression mathDisplay"/);
  assert.match(architecture, /class="katex-display"/);
  assert.match(architecture, /region_depth_deep_m/);
  assert.match(architecture, /region_depth_shallow_m/);
  assert.match(architecture, /region_length_north_m/);
  assert.match(architecture, /region_length_east_m/);

  const topographyResponse = await render("/docs/input-files");
  assert.equal(topographyResponse.status, 200);
  const topography = await topographyResponse.text();
  assert.match(topography, /Choose maximum_octant_level/);
  assert.match(topography, /simulation_shear_velocity_min/);
  assert.match(topography, /simulation_wave_max_freq_hz/);
  assert.match(topography, /simulation_node_per_wavelength/);
  assert.match(topography, /maximum_octant_level/);
  assert.match(topography, /61\.035.*30\.518/s);
  assert.match(topography, /Opt in to wavelength-based surface sizing/);
  assert.match(topography, /topography_surface_refinement=wavelength/);
  assert.match(topography, /not a global cap on wavelength refinement/);
  assert.match(topography, /Hydrostatic\/geostatic initialization/);

  const parametersResponse = await render("/docs/parameters-domain-time");
  assert.equal(parametersResponse.status, 200);
  const parameters = await parametersResponse.text();
  assert.match(parameters, /common dyadic increment/);
  assert.match(parameters, /deep minus shallow depth/);

  const topographyParametersResponse = await render("/docs/parameters-topography-buildings");
  assert.equal(topographyParametersResponse.status, 200);
  const topographyParameters = await topographyParametersResponse.text();
  assert.match(topographyParameters, /region_length_north_m/);
  assert.match(topographyParameters, /region_length_east_m/);
  assert.match(topographyParameters, /region_depth_deep_m/);
  assert.match(topographyParameters, /region_depth_shallow_m/);
  assert.match(topographyParameters, /simulation_shear_velocity_min/);
  assert.match(topographyParameters, /simulation_wave_max_freq_hz/);
  assert.match(topographyParameters, /simulation_node_per_wavelength/);
  assert.doesNotMatch(topographyParameters, /north length,east length,deep-shallow thickness/);
  assert.doesNotMatch(topographyParameters, /maximum frequency×nodes per wavelength/);
});

test("searches complete documentation text and parameter names", async () => {
  const { matchesSearchQuery, normalizeSearchText } = await import(
    "../app/docs/search-utils.mjs"
  );
  const searchable = normalizeSearchText(
    "Pressure-dependent models require approximate_geostatic_state = yes.",
  );
  assert.equal(matchesSearchQuery(searchable, "geostatic pressure"), true);
  assert.equal(matchesSearchQuery(searchable, "approximate_geostatic_state"), true);
  assert.equal(matchesSearchQuery(searchable, "unrelated keyword"), false);

  const assetsDirectory = fileURLToPath(new URL(staticExport ? "../out/_next/static/" : "../dist/client/assets/", import.meta.url));
  const clientJavaScript = filesUnder(assetsDirectory)
    .filter((path) => path.endsWith(".js"))
    .map((path) => readFileSync(path, "utf8"))
    .join("\n");
  assert.match(clientJavaScript, /approximate_geostatic_state/);
  assert.match(clientJavaScript, /taper_topography\.py/);
  assert.match(clientJavaScript, /PressureDependMultiYield03/);
});

test("documents the current single parameter-file interface", async () => {
  const response = await render("/docs/parameter-file");
  assert.equal(response.status, 200);

  const html = await response.text();
  assert.match(html, /The parameters\.in file/);
  assert.match(html, /psolve parameters\.in/);
  assert.match(html, /not the current psolve interface/);
  assert.match(html, /all 249 unique keys currently consumed from parameters\.in/);
  assert.match(html, /Source-specific controls.*source\.in/);
});

test("renders the categorized parameter reference", async () => {
  const response = await render("/docs/parameters-coupled-physics");
  assert.equal(response.status, 200);

  const html = await response.text();
  assert.match(html, /Acoustic and poroelastic parameters/);
  assert.match(html, /acoustic_free_surface_gravity/);
  assert.match(html, /Z=topographybase_zcoord/);
  assert.match(html, /excluding buried water beneath land/);
  assert.match(html, /complete horizontal top faces must coincide with local topography over water/);
  assert.match(html, /acoustic_interface_face_filter/);
  const faceFilterRow = html.split("</tr>").find((row) => row.includes("acoustic_interface_face_filter"));
  assert.ok(faceFilterRow, "the acoustic interface filter must appear in the parameter reference");
  assert.match(faceFilterRow, /Default: all/);
  assert.match(faceFilterRow, /all couples every detected face/);
  assert.match(faceFilterRow, /default alias follows all/);
  const methodRow = html.split("</tr>").find((row) => row.includes("acoustic_interface_method"));
  assert.ok(methodRow);
  assert.match(methodRow, /Default: field-specific-constraints for coupled scalar acoustics/);
  assert.match(methodRow, /When omitted, uses shared-geometric-constraints if acoustics, scalar acoustics or interface coupling is disabled/);
  assert.match(methodRow, /Explicit unsupported choices still abort/);
  assert.match(methodRow, /Coupled VT topography must explicitly select shared-geometric-constraints/);
  assert.match(methodRow, /local-matching/);
  assert.match(methodRow, /imported mismatched meshes abort/);
  assert.match(methodRow, /Field-specific constraints support FEM topography/);
  assert.match(methodRow, /computation_method=fem/);
  assert.match(methodRow, /owner-confirmed carved voids/);
  assert.match(methodRow, /separate acoustic, solid and pore-pressure maps/);
  assert.match(methodRow, /Both support buildings, eligible nonlinear constitutive models and BrickUP poroelasticity/);
  assert.match(methodRow, /Local matching requires include_topography=no/);
  assert.match(methodRow, /does not guarantee CPML stability/);
  const projectionRow = html.split("</tr>").find((row) => row.includes("acoustic_interface_projection"));
  assert.ok(projectionRow);
  assert.match(projectionRow, /All three constraint methods resolve actual remote solid vertices for face_mass/);
  assert.match(projectionRow, /shared-geometric point path is unchanged/);
  const qualityRow = html.split("</tr>").find((row) => row.split("</td>")[0].includes("acoustic_quality_factor"));
  assert.ok(qualityRow);
  assert.match(qualityRow, /Float \(dimensionless\)/);
  assert.match(qualityRow, /Default: 0 \(lossless\)/);
  assert.match(qualityRow, /disables interior acoustic damping/);
  assert.match(qualityRow, /verified in the documented flat-water CPS comparisons and MPI\/backend consistency tests/);
  assert.match(qualityRow, /within their recorded configurations and durations/);
  assert.doesNotMatch(qualityRow, /\b(?:100|500)\b|trade.?off/i);
  assert.match(html, /poroelastic_fluid_bulk_modulus_pa/);
  assert.match(html, /Optional|Legacy/);
  assert.match(html, /parameterReferenceMain/);
  assert.match(html, /data-label="Parameter"/);
});

test("parallelizes host, CUDA launch, and nonlinear PTX assembly", async () => {
  const sourceRoot = fileURLToPath(new URL("../../quake/forward/", import.meta.url));
  if (!existsSync(sourceRoot)) return;

  const makefile = readFileSync(`${sourceRoot}Makefile`, "utf8");
  const launchHeader = readFileSync(`${sourceRoot}psolve_gpu_launch.h`, "utf8");
  assert.match(makefile, /PSOLVE_OBJECTS\s*=\s*psolve\.o psolve_gpu_launch\.o/);
  assert.match(makefile, /psolve\.o:\s*psolve\.cc[\s\S]*?\$\(CXX\)[^\n]*-c \$< -o \$@/);
  assert.match(
    makefile,
    /psolve_gpu_launch\.o:\s*psolve_gpu_launch\.cu[\s\S]*?\$\(NVCC\)[^\n]*-c \$< -o \$@/,
  );
  assert.match(
    makefile,
    /NONLINEAR_GPU_OBJECTS\s*=\s*nonlinear\/nonlinear_gpu_runtime\.o[\s\S]*?nonlinear\/nonlinear_gpu_pdmy03\.o/,
  );
  assert.match(
    makefile,
    /nonlinear\/nonlinear_gpu_%\.o:\s*nonlinear\/nonlinear_gpu_%\.cu[\s\S]*?\$\(NVCC\)[^\n]*-c \$< -o \$@/,
  );
  const gpu = readFileSync(
    `${sourceRoot}nonlinear/nonlinear_gpu_implementation.cuh`,
    "utf8",
  );
  const runtime = readFileSync(`${sourceRoot}nonlinear/nonlinear_gpu_runtime.cu`, "utf8");
  const classic = readFileSync(`${sourceRoot}nonlinear/nonlinear_gpu_classic.cu`, "utf8");
  const md = readFileSync(
    `${sourceRoot}nonlinear/nonlinear_gpu_manzari_dafalias.cu`,
    "utf8",
  );
  const pdmy03 = readFileSync(`${sourceRoot}nonlinear/nonlinear_gpu_pdmy03.cu`, "utf8");
  const materialPoint = readFileSync(
    `${sourceRoot}nonlinear/nonlinear_gpu_material_point.cu`,
    "utf8",
  );
  assert.match(
    gpu,
    /#ifndef NL_GPU_MATERIAL_POINT_ONLY\s+__device__ static inline void nl_record_error[\s\S]*?__device__ static tensor_t nl_tetra_strain[\s\S]*?#endif\s+__device__ static inline double nl_gqh_value/,
  );
  assert.match(gpu, /#ifndef NL_GPU_PRODUCTION_ONLY\s+__global__ static void kernel_nonlinear_material_point_history/);
  assert.match(gpu, /#if defined\(NL_GPU_CLASSIC_ONLY\)[\s\S]*?__device__ static inline void nl_accumulate_stress_force/);
  assert.match(runtime, /#define NL_GPU_RUNTIME_ONLY/);
  assert.match(classic, /#define NL_GPU_CLASSIC_ONLY/);
  assert.match(md, /#define NL_GPU_MANZARI_DAFALIAS_ONLY/);
  assert.match(pdmy03, /#define NL_GPU_PDMY03_ONLY/);
  assert.match(materialPoint, /#define NL_GPU_MATERIAL_POINT_ONLY/);
  assert.match(
    gpu,
    /NL_GPU_PDMY03_ONLY[\s\S]*?ok = pdmy03_update[\s\S]*?NL_GPU_MANZARI_DAFALIAS_ONLY[\s\S]*?ok = nl_update_manzari_dafalias[\s\S]*?ok = nl_update_material/,
  );
  assert.match(
    gpu,
    /g_config\.material_model == MANZARI_DAFALIAS[\s\S]*?nonlinear_gpu_launch_manzari_dafalias[\s\S]*?g_config\.material_model == PRESSURE_DEPEND_MULTIYIELD_03[\s\S]*?nonlinear_gpu_launch_pdmy03[\s\S]*?nonlinear_gpu_launch_classic/,
  );
  assert.ok(
    launchHeader.indexOf('#include "quake_util.h"') <
      launchHeader.indexOf('#include "topography.h"'),
    "noyesflag_t must be defined before topography and launch declarations",
  );

  const html = await (await render("/docs/installation")).text();
  assert.match(html, /psolve\.cc/);
  assert.match(html, /psolve_gpu_launch\.cu/);
  assert.match(html, /make -j 8/);
  assert.match(html, /nonlinear_gpu_pdmy03\.o/);
  assert.match(html, /PTX assembly/);
});

test("documents acoustic station recovery and its historical fallback", async () => {
  const html = await (await render("/docs/parameters-coupled-physics")).text();
  const row = html.split("</tr>").find((part) => part.includes("acoustic_station_gradient"));
  assert.ok(row);
  assert.match(row, /Default: linear/);
  assert.match(row, /linear \| element/);
  assert.match(row, /inside CPML/);
  assert.match(row, /plane output/);
  const guide = await (await render("/docs/physics")).text();
  assert.match(guide, /Sender and receiver mappings use global node ID order for every interface method/);
  const source = new URL("../../quake/forward/psolve.cc", import.meta.url);
  if (existsSync(source)) {
    const solver = readFileSync(source, "utf8");
    assert.match(solver, /acoustic_station_gradient\[64\]\s*=\s*"linear"/);
    for (const mode of ["linear", "element"])
      assert.ok(solver.includes(`strcasecmp(acoustic_station_gradient, "${mode}")`));
    assert.match(solver, /int_message\[66\] = Param\.acousticStationGradientRecovered/);
    assert.match(solver, /Param\.acousticStationGradientRecovered = int_message\[66\]/);
  }
});

test("documents every parameters.in key parsed by the current source tree", () => {
  const referencePath = fileURLToPath(new URL("../app/docs/parameter-reference.ts", import.meta.url));
  const reference = readFileSync(referencePath, "utf8");
  const documented = [...reference.matchAll(/\bp\("([^"]+)"/g)].map((match) => match[1]);
  assert.equal(new Set(documented).size, 251);
  assert.equal(documented.length, 251, "parameter definitions must not be duplicated");
  assert.ok(documented.includes("topography_surface_refinement"));
  assert.match(reference, /p\("topography_surface_refinement", "Optional", "String \(case-insensitive\)", "Default: maximum_level; maximum_level \| wavelength"/);
  assert.ok(documented.includes("cpml_acoustic_sponge_fraction"));
  assert.ok(documented.includes("rivasand02_biasvolume"));
  assert.match(reference, /p\("rivasand02_biasvolume", "Optional", "Integer \(dimensionless\)", "Default: 1; 0 \| 1 \| 2/);
  assert.match(reference, /p\("rivasand02_reversaltype", "Optional", "Integer \(dimensionless\)", "Default: 3; 1 \| 2 \| 3/);
  assert.ok(documented.includes("station_binary"));
  assert.ok(documented.includes("station_ascii_buffer_kib"));
  assert.match(reference, /p\("station_ascii_buffer_kib", "Optional", "Integer KiB \(1024 bytes\)", "Default: 64; 1\.\.2097151/);
  assert.ok(documented.includes("station_output_binary_file"));
  assert.ok(!documented.includes("drm_points_file"));
  assert.ok(!documented.includes("drm_output_file"));
  assert.match(reference, /p\("cpml_clamp_outer", "Optional", "Yes\/no", "Default: yes"/);
  assert.match(reference, /p\("acoustic_quality_factor", "Optional", "Float \(dimensionless\)", "Default: 0 \(lossless\)"/);
  assert.match(reference, /p\("acoustic_interface_face_filter", "Optional", "String", "Default: all;/);

  const sourceRoot = fileURLToPath(new URL("../../quake/forward/", import.meta.url));
  if (!existsSync(sourceRoot)) return;

  const solver = readFileSync(`${sourceRoot}psolve.cc`, "utf8");
  const topographySolver = readFileSync(`${sourceRoot}topography.cc`, "utf8");
  assert.match(topographySolver, /char surface_refinement\[64\] = "maximum_level"/);
  for (const mode of ["maximum_level", "wavelength"]) {
    assert.ok(topographySolver.includes(`strcasecmp(surface_refinement, "${mode}")`));
  }
  assert.match(solver, /\.theAcousticQualityFactor\s*=\s*0\.0\s*,/);
  assert.ok(documented.includes("acoustic_interface_method"));
  assert.match(solver, /acoustic_interface_method\[64\]\s*=\s*"field-specific-constraints"/);
  for (const method of ["shared-geometric-constraints", "field-specific-constraints", "local-matching"]) {
    assert.ok(solver.includes(`strcasecmp(acoustic_interface_method, "${method}")`));
  }
  assert.match(solver, /get_topo_meth\(\) != FEM/);
  assert.match(solver, /VT cut-tetrahedron constraints/);
  assert.match(solver, /acoustic_interface_method=local-matching does not support/);
  const interfaceGuard = solver.slice(
    solver.indexOf("if (Param.theAcousticInterfaceMethod != ACOUSTIC_INTERFACE_ORIGINAL &&"),
    solver.indexOf('strcpy(Param.theCheckPointingDirOut'),
  );
  assert.doesNotMatch(interfaceGuard, /Param\.include(?:Buildings|NonlinearAnalysis|Poroelastic)/);
  assert.doesNotMatch(interfaceGuard, /Param\.theAcousticVelocityPotential/);
  assert.match(interfaceGuard, /Param\.acousticPressureFormulation != YES/);
  assert.match(reference, /Both alternatives require include_acoustic=yes, acoustic_interface_coupling=yes, pressure or velocity_potential, and full Hex8 solids/);
  assert.match(reference, /Both support buildings, eligible nonlinear constitutive models and BrickUP poroelasticity/);
  assert.match(reference, /separate acoustic, solid and pore-pressure maps/);
  assert.match(solver, /\.theCpmlClampOuter\s*=\s*YES\b/);
  assert.match(solver, /strcasecmp\(cpml_clamp_outer, "no"\) == 0\)\s*\{\s*Param\.theCpmlClampOuter = NO;/);
  assert.match(solver, /\.theAcousticInterfaceFaceFilter\s*=\s*ACOUSTIC_INTERFACE_FACE_FILTER_ALL/);
  assert.match(solver, /acoustic_interface_face_filter\[64\]\s*=\s*"all"/);
  assert.match(solver, /int acousticInterfaceFaceFilter\s*=\s*ACOUSTIC_INTERFACE_FACE_FILTER_ALL/);
  const filterParser = solver.slice(
    solver.indexOf('if ((strcasecmp(acoustic_interface_face_filter, "all")'),
    solver.indexOf('if ((strcasecmp(acoustic_poroelastic_interface, "no_flow")'),
  );
  const bathymetryBranch = filterParser.slice(filterParser.indexOf('} else if ((strcasecmp(acoustic_interface_face_filter, "bathymetry")'));
  assert.match(bathymetryBranch, /"reconstructed_bathymetry"[\s\S]*?"bathymetry_normal"[\s\S]*?ACOUSTIC_INTERFACE_FACE_FILTER_BATHYMETRY/);
  assert.doesNotMatch(bathymetryBranch, /"default"/);
  assert.match(filterParser.slice(0, filterParser.indexOf("} else if")), /"default"[\s\S]*?ACOUSTIC_INTERFACE_FACE_FILTER_ALL/);

  const parserFiles = [
    "psolve.cc",
    "topography.cc",
    "nonlinear.cc",
    "buildings.cc",
    "drm.cc",
    "io_drm.cc",
    "drm_planewaves.cc",
    "meshformatlab.cc",
    "io_planes.cc",
    "io_planes_back.cc",
    "qmesh.cc",
  ];
  const parsed = new Set([
    "source_directory",
    "output_displacement_file",
    "output_velocity_file",
    "output_stats_file",
    "output_debug_file",
    "output_planes",
  ]);
  const parserPattern = /(?:\b(?:hu_)?parsetext(?:_v)?|\bhu_config_get_[A-Za-z0-9_]+|\bparsedarray2?|\briva_read_reversal_type|\briva_read_bias_volume|\bread_boxes)\s*\(\s*[^,]+,\s*"([^"]+)"/gs;
  for (const filename of parserFiles) {
    const source = readFileSync(`${sourceRoot}${filename}`, "utf8");
    for (const match of source.matchAll(parserPattern)) parsed.add(match[1]);
  }
  parsed.delete("4D_output_file");
  parsed.delete("drm_numberofelements");

  const missing = [...parsed].filter((name) => !documented.includes(name)).sort();
  const extra = documented.filter((name) => !parsed.has(name)).sort();
  assert.deepEqual(missing, [], `undocumented parsed parameters: ${missing.join(", ")}`);
  assert.deepEqual(extra, [], `documented names not found in the parser audit: ${extra.join(", ")}`);
});

test("documents elevation-only seabed geometry and depth-reference contract", async () => {
  const html = await (await render("/docs/input-files")).text();
  const section = html.match(/<section\b[^>]*id="explicit-seabed-geometry"[\s\S]*?<\/section>/)?.[0];
  assert.ok(section);
  for (const name of ["seabed_geometry_file", "seabed_above_surface_policy", "3D_velocity_model_depth_reference", "seabed_water_vp_m_s", "seabed_water_density_kg_m3"]) {
    assert.ok(section.includes(name), `missing seabed input example: ${name}`);
  }
  assert.match(section, /solid-only/i);
  assert.match(section, /east-fast/);
  assert.match(section, /topography minus bed/);
  assert.match(section, /1e-8 m/);
  assert.match(section, /No separate water-level parameter or wet\/dry mask/);
  assert.match(section, /One number per line is valid/);
  assert.match(section, /Distinct water bodies may have different elevations/);
  assert.match(section, /sloping water surfaces are rejected/);
  assert.doesNotMatch(html, /seabed_water_level_m/);
  assert.match(section, /both full and sqd/);
  assert.match(section, /20, 120 or 1120/);
  assert.match(section, /dry land can be below zero/);
  assert.match(section, /1000 m × 1000 m/);
  assert.match(section, /Porous sediment depth limits are measured below the bed/);
  assert.match(section, /supercomputer validation/);
  assert.match(section, /include_topography=no/);
  assert.match(section, /flat at solver Z=0/);
  assert.match(section, /No topography.in/);
  assert.match(section, /no forced topography surface-refinement band/);
  assert.match(section, /spacing must divide both horizontal domain lengths/);
  assert.match(section, /CVM mesh offsets must be zero/);
  assert.match(section, /model_z and surface coincide/);
  assert.match(section, /These setup checks are not a CUDA\/MPI propagation benchmark/);
  assert.match(section, /independent CPS Green functions/);
  assert.match(section, /Submerged ramped equilibrium/);
  assert.match(section, /GPU scheduling and MPI transport are mocked/);
  assert.match(section, /seabed_above_surface_policy=error \(default\)/);
  assert.match(section, /seabed_above_surface_policy=clip/);
  assert.match(section, /lowered-node count/);
  assert.match(section, /maximum lowering in metres/);
  assert.match(section, /Original input files are unchanged/);
  assert.match(section, /full CUDA\/MPI propagation validation of clipping remains outstanding/);
  const sourceRoot = "../quake/forward/";
  const solver = readFileSync(`${sourceRoot}psolve.cc`, "utf8");
  const geometry = readFileSync(`${sourceRoot}topography.cc`, "utf8");
  assert.match(solver, /char depth_reference\[64\] = "legacy"/);
  assert.match(solver, /char seabed_policy\[64\] = "error"/);
  for (const policy of ["error", "clip"]) {
    assert.ok(solver.includes(`strcasecmp(seabed_policy, "${policy}")`));
  }
  assert.match(solver, /MPI_Bcast\(&Param\.theSeabedClipAboveSurface, 1, MPI_INT, 0, comm_solver\)/);
  assert.match(solver, /has_seabed_policy == 1/);
  assert.match(solver, /topo_seabed_init\(Global\.myID, Param\.theSeabedGeometryFile,\s*Param\.theSeabedClipAboveSurface\)/);
  assert.match(solver, /topo_seabed_init_flat\(Global\.myID, Param\.theSeabedGeometryFile,\s*Param\.theDomainX, Param\.theDomainY,\s*Param\.theSeabedClipAboveSurface\)/);
  assert.match(geometry, /water_depth < -SEABED_HEIGHT_TOLERANCE_M && !clip_above_surface/);
  assert.match(geometry, /topo_seabed_init\(myID, path, clip_above_surface\)/);
  assert.match(solver, /have_buildings == YES \|\| Param\.the3DDepthReference == 0/);
  assert.match(solver, /MPI_Bcast\(Param\.theSeabedGeometryFile/);
  assert.match(geometry, /MPI_Bcast\(theSeabedInfo, ntp, MPI_DOUBLE/);
  const seabedParser = solver.slice(solver.indexOf('  char depth_reference[64]'), solver.indexOf('  parsetext(fp, "source_directory_output"'));
  assert.doesNotMatch(seabedParser, /haveTopography != YES/);
  assert.match(solver, /topo_seabed_init_flat\(Global\.myID, Param\.theSeabedGeometryFile/);
  assert.match(geometry, /void topo_seabed_init_flat/);
  assert.match(geometry, /MPI_Bcast\(counts, 2, MPI_INT/);
  assert.match(geometry, /thebase_zcoord = 0\.0;\s+theEtreeType = FULL;/);
  assert.doesNotMatch(solver + geometry, /seabed_water_level_m|theSeabedWaterLevel|theSeabedWet|topo_seabed_water_surface_z/);
});

test("documents the default acoustic quality factor and its validation scope", async () => {
  const html = await (await render("/docs/physics")).text();
  const section = html.match(/<section\b[^>]*id="acoustic-quality-factor"[\s\S]*?<\/section>/)?.[0];
  assert.ok(section);
  assert.match(section, /acoustic_quality_factor defaults to 0/);
  assert.match(section, /disables interior acoustic damping/);
  assert.match(section, /verified in the documented flat-water CPS comparisons and MPI\/backend consistency tests/);
  assert.match(section, /recorded configurations and durations/);
  assert.doesNotMatch(section, /\b(?:100|500)\b|trade.?off/i);
  assert.match(html, /acoustic_quality_factor=0 in every run/);
});

test("explains experimental constraints, remeshing and CPML limitations", async () => {
  const html = await (await render("/docs/physics")).text();
  assert.match(html, /id="acoustic-interface-method"/);
  assert.match(html, /acoustic_interface_method = field-specific-constraints/);
  assert.match(html, /Phase B imports do not refine/);
  assert.match(html, /unsupported mixed-corner stencil stops setup/);
  assert.match(html, /exponential growth with CPML/);
  assert.match(html, /not a prescribed uniform-width refinement band/);
  assert.match(html, /Field-specific constraints also support FEM topography/);
  assert.match(html, /acoustic_formulation=pressure or velocity_potential/);
  assert.match(html, /Pressure mode couples solid normal acceleration/);
  assert.match(html, /Unresolved remote neighbors still abort/);
  assert.match(html, /VT cut-tetrahedron topography remains unsupported/);
});

test("wires the optional acoustic gravity free surface through setup and CUDA", () => {
  const sourceRoot = fileURLToPath(new URL("../../quake/forward/", import.meta.url));
  if (!existsSync(sourceRoot)) return;

  const psolve = readFileSync(`${sourceRoot}psolve.cc`, "utf8");
  const header = readFileSync(`${sourceRoot}psolve.h`, "utf8");
  const kernel = readFileSync(`${sourceRoot}kernel.cu`, "utf8");
  assert.match(psolve, /parsetext\(fp, "acoustic_free_surface_gravity"/);
  assert.match(psolve, /gravitySurfaceMass[\s\S]*?edgesize[\s\S]*?4\.0[\s\S]*?rho/);
  assert.match(psolve, /gravity_surface_kprod = surface->inertia.mass/);
  assert.match(psolve, /surface->inertia = cpml_inertia[\s\S]*?2, 0\.0, dt/);
  assert.match(psolve, /psolve_launch_pml_acoustic_gravity_surface/);
  assert.match(header, /ACOUSTIC_GRAVITY_ACCELERATION_M_S2 9\.81/);
  assert.match(kernel, /pressureFreeSurfaceNodeFlagsDevice[\s\S]*?acousticFreeSurfaceGravity == 0/);
  assert.match(kernel, /kernelPMLAcousticGravitySurfaceStep/);
});

test("documents the phase launcher contract", async () => {
  const response = await render("/docs/running");
  assert.equal(response.status, 200);

  const html = await response.text();
  assert.match(html, /run_phaseA_phaseB\.sh/);
  assert.match(html, /PSOLVE_BIN/);
  assert.match(html, /BASE_IN/);
  assert.match(html, /Path to the complete parameters\.in file/);
  assert.match(html, /file, not a directory/);
  assert.match(html, /Hercules_GPU/);
  assert.doesNotMatch(html, /GPU_new/);
});

test("lists GSL as a required build dependency", async () => {
  const response = await render("/docs/installation");
  assert.equal(response.status, 200);

  const html = await response.text();
  assert.match(html, /GNU Scientific Library \(GSL\)/);
  assert.match(html, /-lgsl -lgslcblas/);
  assert.match(html, /GSL include directory/);
});

test("documents nonlinear and plastic simulation capability", async () => {
  const physicsResponse = await render("/docs/physics");
  assert.equal(physicsResponse.status, 200);
  const physics = await physicsResponse.text();
  assert.match(physics, /Nonlinear\/plastic solid/);
  assert.match(physics, /plastic strain/i);
  assert.match(physics, /vonMises_ep/);
  assert.match(physics, /MohrCoulomb/);
  assert.match(physics, /DruckerPrager/);
  assert.match(physics, /ManzariDafalias/);
  assert.match(physics, /PressureDependMultiYield03/);
  assert.match(physics, /RIVASand/);
  assert.match(physics, /RIVASand02/);
  assert.match(physics, /RIVASandCustom02/);
  assert.match(physics, /Custom02 uses its own e_max\/e_min for both the initial and reference relative densities/);
  assert.match(physics, /Markers 10 through 18 identify the nine reversal\/BiasVolume combinations; default marker 15 identifies type 3 with mode 1/);
  assert.match(physics, /retired RIVASandIntermediateBiasResearch input is not restored/);
  assert.match(physics, /not an alias for, replacement for, or retuning of production RIVASand/);
  assert.match(physics, /BrickUP remains the sole owner of pore pressure/);
  assert.match(physics, /frozen Ottawa F65/);
  assert.match(physics, /standard automatically generated nested yield surfaces/);
  assert.match(physics, /adaptive explicit Modified Euler/);
  assert.match(physics, /integrates the 1-D density column/);
  assert.match(physics, /pressure DOF starts at zero as excess pore pressure/);
  assert.match(physics, /disables the zero-permeability undrained fast path whenever nonlinear analysis is active/);
  assert.match(physics, /legacy Istanbul nl_id classification/);
  assert.match(physics, /Vs range covered by material_properties_list/);
  assert.match(physics, /non-Istanbul vonMises_baE/);
  assert.match(physics, /built-in nl_id constants/);
  assert.match(physics, /Nonlinear setup rules/);
  for (const model of [
    "linear",
    "vonMises_ep",
    "vonMises_fa",
    "vonMises_faM",
    "vonMises_baE",
    "vonMises_baH",
    "vonMises_GQH",
    "vonMises_MKZ",
    "vonMises_RO",
    "MohrCoulomb",
    "DruckerPrager",
    "ManzariDafalias",
    "PressureDependMultiYield03",
    "RIVASand",
    "RIVASandCustom",
    "RIVASand02",
    "RIVASandCustom02",
  ]) {
    assert.match(physics, new RegExp(`material_model = ${model}`));
  }
  for (const pressureDependentModel of [
    "MohrCoulomb",
    "DruckerPrager",
    "ManzariDafalias",
    "PressureDependMultiYield03",
    "RIVASand",
    "RIVASandCustom",
    "RIVASand02",
    "RIVASandCustom02",
  ]) {
    assert.match(
      physics,
      new RegExp(
        `approximate_geostatic_state = yes[\\s\\S]{0,500}material_model = ${pressureDependentModel}`,
      ),
    );
  }
  assert.match(physics, /The displayed constants are synthetic regression or syntax values/);
  assert.match(physics, /Acoustic–elastic example/);
  assert.match(physics, /acoustic_interface_projection = point/);
  assert.match(physics, /Poroelastic–elastic example/);
  assert.match(physics, /Nonlinear PDMY03–BrickUP example/);
  assert.match(physics, /RIVA-Sand–BrickUP example/);
  assert.match(physics, /BrickUP independently advances excess fluid pressure/);
  assert.match(physics, /contributes -alpha\*p_f\*I once/);
  assert.match(physics, /acoustic_poroelastic_interface = open_pore/);
  assert.match(physics, /type_of_damping = rayleigh/);
  assert.match(physics, /type_of_damping = mass/);
  assert.match(physics, /type_of_damping = bkt2/);
  assert.match(physics, /BKT setup rejects nonfinite or nonpositive solid Q/);
  assert.match(physics, /Qk=1000 \(lossless bulk\) only for affected elements/);
  assert.match(physics, /Valid bulk and shear attenuation remain unchanged/);
  assert.match(physics, /global negative-Qk bulk fallbacks count/);
  assert.match(physics, /original supplied Qp is no longer reproduced/);
  assert.match(physics, /must be regenerated to apply the fallback/);
  assert.match(physics, /positive relaxed stiffness/);
  assert.match(physics, /no material queries, mesh refinement, extra full mesh traversal, or time-step work/);
  assert.match(physics, /Imported BKT2 meshes have their stored coefficients checked/);
  assert.match(physics, /BKT2 GPU execution and comparison/);
  assert.match(physics, /Host-staged execution retains the reference kernels/);
  assert.match(physics, /HERCULES_BKT2_FUSED=0/);
  assert.match(physics, /HERCULES_FORCE_GENERIC_DAMPING=1 takes precedence/);
  assert.match(physics, /Inactive sectors are excluded from rate comparisons/);
  assert.match(physics, /BKT3\/BKT3F still check their third mechanism/);
  assert.match(physics, /fallback reason rank counts/);
  assert.match(physics, /Reason counts may overlap/);
  assert.match(physics, /not allocated GPU memory/);
  assert.match(physics, /Atomic force accumulation may change roundoff/);
  assert.match(physics, /Qs=C\+alpha×Vs_km_s\^beta/);
  assert.match(physics, /Default CPML example/);
  assert.match(physics, /cpml_reflection_target = 0\.1/);
  assert.match(physics, /cpml_pml_material_damping = rayleigh/);
  assert.match(physics, /script\/taper_topography\.py/);

  const parameterResponse = await render("/docs/parameters-nonlinear");
  assert.equal(parameterResponse.status, 200);
  const parameters = await parameterResponse.text();
  assert.match(parameters, /Nonlinear and plastic analysis parameters/);
  assert.match(parameters, /include_nonlinear_analysis/);
  assert.match(parameters, /material_plasticity_type/);
  assert.match(parameters, /inclusive nonlinear-element range/);
  assert.match(parameters, /alpha_cohesion is Su in Pa/);
  assert.match(parameters, /plastic-strain evolution/);
  assert.match(parameters, /ManzariDafalias: N × 17/);
  assert.match(parameters, /PDMY03: N × 22/);
  assert.match(parameters, /RIVASand: N × 2 or 3/);
  assert.match(parameters, /RIVASandCustom: N × 13/);
  assert.match(parameters, /RIVASand02: N × 2 or 3/);
  assert.match(parameters, /RIVASandCustom02: N × 13/);
  assert.match(parameters, /distinct research successor/);
  assert.match(parameters, /Vs_limit Dr G0 M kd h m zeta e_max e_min Q R n_G/);
  assert.match(parameters, /BiasVolume is selected separately with rivasand02_biasvolume; h_p0 and reversalLatch are not exposed/);
  assert.match(parameters, /9585f0c155c9444885c5115e7753a1a5d97c783e2c685938a49a65435c9e8f83/);
  assert.match(parameters, /surface_count is an integer from 2 through 40/);
  assert.match(parameters, /Vs_limit G0 e_init Mc c lambda_c/);
  assert.match(parameters, /maximum absolute yield-surface residual/);
  assert.match(parameters, /documented MohrCoulomb, DruckerPrager, ManzariDafalias, PDMY03, RIVASand, RIVASandCustom, RIVASand02, and RIVASandCustom02/);

  const gpuResponse = await render("/docs/gpu-execution");
  assert.equal(gpuResponse.status, 200);
  const gpu = await gpuResponse.text();
  assert.match(gpu, /HERCULES_NONLINEAR_BACKEND/);
  assert.match(gpu, /HERCULES_NONLINEAR_STATS_INTERVAL/);
  assert.match(gpu, /one-step[\s\S]*multi-step/);
  assert.match(gpu, /ManzariDafalias production kernel omits the four legacy previous-state shadow arrays/);
  assert.match(gpu, /1,216 bytes of unused device allocation per nonlinear element/);
  assert.match(gpu, /ManzariDafalias keeps its larger failure-only diagnostic record in per-thread shared memory/);
  assert.match(gpu, /PDMY03 retains the legacy shadow allocation and writes/);
  assert.match(gpu, /CUDA nonlinear constitutive and force path is the default/);
  assert.match(gpu, /GPU by default/);
  assert.match(gpu, /constitutive calculations in double precision/);
  assert.match(gpu, /CPU\/GSL reference/);
  assert.match(gpu, /Material-point CPU\/GPU regression/);
  assert.match(gpu, /run_suite\.sh/);
  assert.match(gpu, /FA\/FAM quartic correctors/);
  assert.match(gpu, /twenty-one-fixture suite/);
  assert.match(gpu, /device-first and currently require the GPU backend/);
  assert.match(gpu, /incremental bounding\/backbone update directly/);
  assert.match(gpu, /true zero strain increment preserves/);
  assert.match(gpu, /purely volumetric increment updates only hydrostatic stress/);
  assert.match(gpu, /deviatoric trial stress below 1\.0e-3 Pa/);
  assert.match(gpu, /OpenSees MultiaxialCyclicPlasticity/);
  assert.match(gpu, /compare_opensees_bae\.py/);
  assert.match(gpu, /opensees_bae_synthetic\.tcl/);
  assert.match(gpu, /--opensees-csv/);
  assert.match(gpu, /17-step legacy BAE suite fixture is not interchangeable/);
  assert.match(gpu, /Path time series and prescribed top-face displacement/);
  assert.match(gpu, /marks its CSV as driver version 2/);
  assert.match(gpu, /rejects older CSVs/);
  assert.match(gpu, /penalty-enforcement shear-path error up to 1e-12/);
  assert.match(gpu, /retains the actual OpenSees strain/);
  assert.match(gpu, /one-brick affine simple-shear test/);
  assert.match(gpu, /integrated shear-work differences/);
  assert.match(gpu, /Ho=0, h=psi0 times G, and beta=1/);
  assert.match(gpu, /dwell and reload/);
  assert.match(gpu, /tensorial shear strains/);
  assert.match(gpu, /compact quadrature-point records/);
  assert.match(gpu, /does not download the full nonlinear state/);
});

test("defaults nonlinear execution to GPU with an explicit CPU override", async () => {
  const response = await render("/docs/configuration");
  assert.equal(response.status, 200);
  assert.match(await response.text(), /default: gpu/);

  const sourceRoot = fileURLToPath(new URL("../../quake/forward/", import.meta.url));
  if (!existsSync(sourceRoot)) return;

  const cpu = readFileSync(`${sourceRoot}nonlinear.cc`, "utf8");
  const gpu = readFileSync(`${sourceRoot}nonlinear/nonlinear_gpu_implementation.cuh`, "utf8");
  const selector = gpu.slice(
    gpu.indexOf("int nonlinear_gpu_requested"),
    gpu.indexOf("int nonlinear_gpu_active"),
  );

  assert.match(selector, /g_backend_request = 1/);
  assert.match(selector, /parse_boolean_cpu\(value\)[\s\S]*g_backend_request = 0/);
  assert.match(gpu, /fprintf\(stderr,\s*"Nonlinear backend: GPU \(default backend/);
  assert.match(cpu, /fprintf\(stderr,\s*"Nonlinear backend: CPU \(explicit/);
  assert.match(await (await render("/docs/gpu-execution")).text(),
    /reported to standard error at startup/);
});

test("documents MPI-IO as the default plane backend and the I/O-pool fallback", async () => {
  const response = await render("/docs/gpu-execution");
  assert.equal(response.status, 200);
  const html = await response.text();
  assert.match(html, /Plane output defaults to mpiio/);
  assert.match(html, /HERCULES_PLANE_IO_MODE is unset or empty/);
  assert.match(html, /select the legacy I\/O-pool backend/);
  assert.match(html, /MPI-IO batches 64 packed snapshots by default/);
  assert.match(html, /The former auto value is no longer accepted/);
  assert.match(html, /set legacy to retain CPU interpolation/);

  const configuration = await render("/docs/configuration");
  assert.equal(configuration.status, 200);
  const config = await configuration.text();
  assert.match(config, /mpiio \| hybrid \| legacy/);
  assert.match(config, /default: mpiio \(dedicated I\/O ranks select legacy\)/);
  assert.doesNotMatch(config, /auto \| hybrid \| mpiio \| legacy/);
});

test("keeps CPML elements out of nonlinear constitutive maps", async () => {
  const sourceRoot = fileURLToPath(new URL("../../quake/forward/", import.meta.url));
  if (!existsSync(sourceRoot)) return;

  const cpu = readFileSync(`${sourceRoot}nonlinear.cc`, "utf8");
  const psolve = readFileSync(`${sourceRoot}psolve.cc`, "utf8");
  const classifier = cpu.slice(
    cpu.indexOf("int isThisElementNonLinear"),
    cpu.indexOf("double interpolate_property_value"),
  );

  assert.match(classifier, /cpml_pml_element\(eindex\)[\s\S]*return NO/);
  assert.ok(
    classifier.indexOf("cpml_pml_element(eindex)") <
      classifier.indexOf("edata->Vs < theVsLimits[0]"),
    "CPML exclusion must precede nonlinear material-range selection",
  );
  assert.match(psolve, /int cpml_pml_element\(int32_t eindex\)/);
  assert.match(psolve, /g_cpmlBuild\.face_mask != NULL/);
  assert.match(psolve, /g_cpmlHostElementBits\[bit_index >> 3\]/);
  assert.match(psolve, /CPML membership is unavailable before nonlinear/);
  assert.match(psolve, /fprintf\(stderr,[\s\S]{0,80}"CPML: derived thickness/);
  assert.match(psolve, /fprintf\(stderr,[\s\S]{0,80}"CPML: %lld elastic/);
  assert.match(psolve, /fprintf\(stderr,[\s\S]{0,80}"CPML: cpml_clamp_outer=yes/);
  assert.doesNotMatch(psolve, /fprintf\(stdout,[\s\S]{0,80}"CPML: derived thickness/);
  assert.doesNotMatch(psolve, /fprintf\(stdout,[\s\S]{0,80}"CPML: %lld elastic/);
  assert.doesNotMatch(psolve, /fprintf\(stdout,[\s\S]{0,80}"CPML: cpml_clamp_outer=yes/);

  const gpu = await (await render("/docs/gpu-execution")).text();
  assert.match(gpu, /solid elements in the absorbing band remain linear elastic/);
  assert.match(gpu, /excluded from both CPU and GPU nonlinear constitutive maps/);

  const parameters = await (await render("/docs/parameters-damping-cpml")).text();
  assert.match(parameters, /Solid PML elements always use linear elasticity/);
  assert.match(parameters, /excluded from nonlinear constitutive updates/);
  const clampRow = parameters.split("</tr>").find((row) => row.includes("cpml_clamp_outer"));
  assert.ok(clampRow, "the CPML clamp must appear in the parameter reference");
  assert.match(clampRow, /Default: yes/);
  assert.match(clampRow, /Explicit no removes the clamp/);
});

test("keeps acoustic elements out of nonlinear constitutive maps", async () => {
  const sourceRoot = fileURLToPath(new URL("../../quake/forward/", import.meta.url));
  if (!existsSync(sourceRoot)) return;

  const cpu = readFileSync(`${sourceRoot}nonlinear.cc`, "utf8");
  const psolve = readFileSync(`${sourceRoot}psolve.cc`, "utf8");
  const classifier = cpu.slice(
    cpu.indexOf("int isThisElementNonLinear"),
    cpu.indexOf("double interpolate_property_value"),
  );

  assert.match(classifier, /element_is_acoustic\(edata\)[\s\S]*return NO/);
  assert.ok(
    classifier.indexOf("element_is_acoustic(edata)") <
      classifier.indexOf("theIstanbulVelocityModel == YES"),
    "acoustic exclusion must precede Istanbul nl_id selection",
  );
  assert.ok(
    classifier.indexOf("element_is_acoustic(edata)") <
      classifier.indexOf("edata->Vs < theVsLimits[0]"),
    "acoustic exclusion must precede nonlinear material-range selection",
  );
  assert.match(psolve, /int element_is_acoustic\(const edata_t\* edata\)/);
  assert.match(psolve, /Param\.includeAcoustic == YES/);
  assert.match(psolve, /Nonlinear\/acoustic policy: acoustic elements are excluded/);

  const gpu = await (await render("/docs/gpu-execution")).text();
  assert.match(gpu, /acoustic elements are excluded from both CPU and GPU solid constitutive maps/);
  assert.match(gpu, /never acquire plastic-strain or nonlinear solid history/);

  const parameters = await (await render("/docs/parameters-nonlinear")).text();
  assert.match(parameters, /Acoustic elements are always excluded/);
  assert.match(parameters, /Acoustic elements remain excluded even when the range starts at zero/);
});

test("dispatches cyclic backbone models before the conventional predictor", () => {
  const sourceRoot = fileURLToPath(new URL("../../quake/forward/", import.meta.url));
  if (!existsSync(sourceRoot)) return;

  const cpu = readFileSync(`${sourceRoot}nonlinear.cc`, "utf8");
  const gpu = readFileSync(`${sourceRoot}nonlinear/nonlinear_gpu_implementation.cuh`, "utf8");
  const cpuUpdate = cpu.slice(
    cpu.indexOf("void material_update"),
    cpu.indexOf("/* computes the derivatives of the flow potential"),
  );
  const gpuUpdate = gpu.slice(
    gpu.indexOf("__device__ static int nl_update_material"),
    gpu.indexOf("__global__ static void kernel_nonlinear_material_point_history"),
  );

  assert.ok(cpuUpdate.indexOf("MatUpd_BA_GSL_v2") < cpuUpdate.indexOf("Fs_pr < Tol_sigma"));
  assert.ok(gpuUpdate.indexOf("nl_update_backbone") < gpuUpdate.indexOf("predictor < 1.0e-3"));
  assert.match(cpu, /A bounding-surface update is incremental[\s\S]*De\.xx == 0\.0/);
  assert.match(gpu, /The model is incremental[\s\S]*de\.xx == 0\.0/);
  assert.match(cpu, /DSdev_measure <= 64\.0 \* DBL_EPSILON \* stress_scale/);
  assert.match(gpu, /dsdev_measure <= 64\.0 \* DBL_EPSILON \* stress_scale/);
  assert.match(cpu, /kappa_n == 100\.0[\s\S]*NONLINEAR_BACKBONE_STRESS_FLOOR_PA/);
  assert.match(gpu, /kappa_n == 100\.0[\s\S]*NONLINEAR_BACKBONE_STRESS_FLOOR_PA/);
});

test("wires Manzari-Dafalias integration schemes through parser, CUDA state, and u-p", () => {
  const sourceRoot = fileURLToPath(new URL("../../quake/forward/", import.meta.url));
  const benchmarkRoot = fileURLToPath(
    new URL("../../benchmark_inputs/nonlinear_material_point/", import.meta.url),
  );
  if (!existsSync(sourceRoot)) return;

  const cpu = readFileSync(`${sourceRoot}nonlinear.cc`, "utf8");
  const gpu = readFileSync(`${sourceRoot}nonlinear/nonlinear_gpu_implementation.cuh`, "utf8");
  const kernel = readFileSync(`${sourceRoot}kernel.cu`, "utf8");
  const header = readFileSync(`${sourceRoot}nonlinear.h`, "utf8");
  const gpuHeader = readFileSync(`${sourceRoot}nonlinear_gpu.h`, "utf8");
  const psolve = readFileSync(`${sourceRoot}psolve.cc`, "utf8");
  const generator = readFileSync(
    `${benchmarkRoot}generate_hercules_manzari_dafalias_path.py`,
    "utf8",
  );
  const fixedForwardEulerFixture = readFileSync(
    `${benchmarkRoot}manzari_dafalias_forward_euler_synthetic.in`,
    "utf8",
  );
  const demandForwardEulerFixture = readFileSync(
    `${benchmarkRoot}manzari_dafalias_forward_euler_demand_synthetic.in`,
    "utf8",
  );
  const materialPointSuite = readFileSync(
    `${benchmarkRoot}run_suite.sh`,
    "utf8",
  );
  const canonicalComparator = readFileSync(
    `${benchmarkRoot}compare_opensees_manzari_dafalias_canonical.py`,
    "utf8",
  );
  const drainedDriver = readFileSync(
    `${benchmarkRoot}opensees_manzari_dafalias_drained_triaxial.tcl`,
    "utf8",
  );
  const undrainedDriver = readFileSync(
    `${benchmarkRoot}opensees_manzari_dafalias_undrained_cyclic_shear.tcl`,
    "utf8",
  );
  const brickupDriver = readFileSync(
    `${benchmarkRoot}opensees_manzari_dafalias_sspbrickup_undrained_triaxial.tcl`,
    "utf8",
  );
  const brickupComparator = readFileSync(
    `${benchmarkRoot}compare_opensees_manzari_dafalias_brickup.py`,
    "utf8",
  );
  const brickupPlot = readFileSync(
    `${benchmarkRoot}plot_opensees_manzari_dafalias_brickup.py`,
    "utf8",
  );

  assert.match(header, /MANZARI_DAFALIAS/);
  assert.match(cpu, /MATERIAL_MODEL_IS_RIVA_SAND_FAMILY\(theMaterialModel\) \? 13 : 16/);
  assert.match(cpu, /theManzariDafaliasProperties\[15\]/);
  assert.match(cpu, /parsetext\(fp, "manzari_dafalias_integration_scheme"/);
  assert.match(cpu, /parsetext\(fp, "manzari_dafalias_forward_euler_substeps"/);
  assert.match(cpu, /parsetext\(fp, "manzari_dafalias_recovery"/);
  assert.match(cpu, /force_accept[\s\S]*unlimited_substeps[\s\S]*opensees/);
  assert.match(gpuHeader, /NL_MD_RECOVERY_STRICT[\s\S]*NL_MD_RECOVERY_OPENSEES/);
  assert.match(gpuHeader, /NL_MD_INTEGRATION_MODIFIED_EULER[\s\S]*NL_MD_INTEGRATION_FORWARD_EULER/);
  assert.match(gpuHeader, /NL_MD_FORWARD_EULER_SUBSTEPS_FIXED[\s\S]*NL_MD_FORWARD_EULER_SUBSTEPS_DEMAND/);
  assert.match(gpu, /nl_update_manzari_dafalias/);
  assert.match(header, /manzari_dafalias_constants_t \*md_materials/);
  assert.match(gpu, /state\.md_materials\+nl_eindex/);
  assert.match(gpu, /Manzari-Dafalias material properties/);
  assert.match(gpu, /Fixed subdivision remains the default/);
  assert.match(gpu, /not the[\s\S]*PDMY03 surface-spacing criterion/);
  // The integrator is specialized per scheme and Forward-Euler subdivision
  // policy so ptxas allocates one path per entry function.
  assert.match(gpu, /template <nl_md_integration_scheme_t MD_SCHEME,[\s\S]*nl_md_forward_euler_substeps_t MD_FE_SUBSTEPS>/);
  assert.match(gpu, /nl_update_manzari_dafalias_scheme<MD_SCHEME,[\s\S]*MD_FE_SUBSTEPS>\(/);
  assert.match(gpu, /if \(MD_SCHEME == NL_MD_INTEGRATION_FORWARD_EULER\)/);
  assert.match(gpu, /NL_MD_INTEGRATION_FORWARD_EULER,[\s\S]*NL_MD_FORWARD_EULER_SUBSTEPS_FIXED/);
  assert.match(gpu, /NL_MD_INTEGRATION_FORWARD_EULER,[\s\S]*NL_MD_FORWARD_EULER_SUBSTEPS_DEMAND/);
  assert.match(gpu, /NL_MD_INTEGRATION_MODIFIED_EULER,[\s\S]*NL_MD_FORWARD_EULER_SUBSTEPS_FIXED/);
  // Per-scheme occupancy target: ptxas fills to the 255-register ceiling
  // unless a residency bound asks for less.
  assert.match(gpu, /__launch_bounds__\(NL_GPU_UPDATE_BLOCK_SIZE,/);
  assert.match(gpu, /NL_MD_FE_MIN_BLOCKS_PER_SM/);
  assert.match(gpu, /#define NL_MD_FE_MIN_BLOCKS_PER_SM 2/);
  assert.match(gpu, /NL_MD_ME_MIN_BLOCKS_PER_SM/);
  assert.match(gpu, /substep < requested/);
  assert.match(gpu, /md_substep_request\([\s\S]{0,200}?cfg->no_substeps\)/);
  assert.match(gpu, /NL_MD_FE_RATIO_STEP[\s\S]*NL_MD_FE_STRAIN_STEP/);
  assert.match(gpu, /cfg->md_integration_scheme == NL_MD_INTEGRATION_FORWARD_EULER/);
  assert.match(gpu, /enforce_md_yield_tolerance[\s\S]*NL_MD_INTEGRATION_MODIFIED_EULER[\s\S]*model_error > cfg\.error_tolerance/);
  assert.match(gpu, /cfg->backbone_error_tolerance/);
  assert.match(gpu, /\*accepted_substeps >= cfg->no_substeps/);
  assert.match(gpu, /cfg->error_tolerance/);
  assert.match(gpu, /if \(fabs\(f\) < tolerance\)/);
  assert.match(gpu, /trial_f <= cfg->error_tolerance/);
  assert.match(gpu, /first\.dgamma == 0\.0 &&[\s\S]*predictor_f <= cfg->error_tolerance/);
  assert.match(gpu, /Do not evaluate a plastic second[\s\S]*loading switch/);
  assert.match(gpu, /alpha evolves during accepted Modified-Euler[\s\S]*alpha_in_work = alpha_work/);
  assert.match(gpu, /first\.hardening_gap < 0\.0 &&[\s\S]*first\.loading_numerator > 0\.0/);
  assert.match(gpu, /diagnostics->correction_residual = correction_residual/);
  assert.match(gpu, /fmax\(\*yield_error, fmax\(final_f, 0\.0\)\)/);
  assert.doesNotMatch(gpu, /fmax\(\*yield_error, fabs\(final_f\)\)/);
  assert.match(gpu, /out\.dalpha = md_zero\(\)/);
  assert.match(gpu, /md_elastic_yield_intersection/);
  assert.match(gpu, /md_stress_segment_projection/);
  assert.match(gpu, /A failed yield[\s\S]*correction must not modify committed history/);
  assert.match(gpu, /dT = fmax\(0\.5\*dT, dT_min\)/);
  assert.match(gpu, /fabs\(denominator\) < 1\.0e-10[\s\S]*neutral-loading branch/i);
  assert.match(gpu, /out\.neutral = 1/);
  assert.match(gpu, /first\.dgamma == 0\.0 && !first\.neutral/);
  assert.match(gpu, /md_force_accept_candidate/);
  assert.match(gpu, /!unlimited_substeps && attempts > 4LL\*cfg->no_substeps/);
  assert.match(gpu, /cfg->md_recovery_mode == NL_MD_RECOVERY_OPENSEES[\s\S]*1\.0e-6/);
  assert.match(gpu, /strain_roundoff =\s*64\.0\*DBL_EPSILON/);
  assert.match(gpu, /fmax\(1\.0e-9,\s*strain_roundoff\/fmax\(md_norm\(de\), DBL_MIN\)\)/);
  assert.match(gpu, /NL_MD_FAILURE_STRESS_CORRECTION_MIN_STEP/);
  assert.match(gpu, /md_minimum_mean_stress/);
  assert.match(gpu, /md_complete_at_pressure_floor/);
  assert.match(gpu, /candidate_mean_stress < minimum_mean_stress/);
  assert.match(gpu, /"Forward-Euler" : "Modified-Euler"/);
  assert.match(gpu, /hardening_gap=/);
  assert.match(gpu, /raw_dgamma=/);
  assert.match(gpu, /HERCULES_NONLINEAR_STATS_INTERVAL/);
  assert.match(gpu, /kernel_nonlinear_workload_stats/);
  assert.match(gpu, /zero_or_inactive_points/);
  assert.match(gpu, /multi_step_fraction/);
  assert.match(gpu, /md_diagnostics_shared\[NL_GPU_UPDATE_BLOCK_SIZE\]/);
  assert.match(gpu, /allocate_previous_state_shadows[\s\S]*NL_ALLOC_COPY\(strains1, strains1\)/);
  assert.match(gpu, /defined\(NL_GPU_CLASSIC_ONLY\) \|\| defined\(NL_GPU_PDMY03_ONLY\)[\s\S]*state\.pstrains1\[nl_eindex\]/);
  assert.match(psolve, /!nonlinear_requires_explicit_pore_pressure\(\)/);
  assert.match(psolve, /Qp traction remains explicit/);
  assert.match(cpu, /return theNonlinearFlag != 0/);
  assert.match(cpu, /effective_overburden =\s*\(double\)edata->sigma_0 - pore_stress/);
  assert.match(cpu, /ecp->sigmaZ_st = effective_overburden/);
  assert.match(cpu, /pore-pressure unknown stores excess pressure/);
  assert.match(gpu, /c->sigmaZ_st\+c->sigmaZ_gradient\*h\*lz\*0\.5/);
  assert.match(gpu, /const tensor_t total_predictor = nl_add\(stress_predictor, static_stress\)/);
  assert.match(psolve, /poroelastic_geostatic_stress_correction/);
  assert.match(psolve, /Param\.thePoroelasticBiotAlpha \* Param\.thePoroelasticFluidDensity \* 9\.81/);
  assert.match(cpu, /below ["\s]*poroelastic_depth_max_m[\s\S]*remain nonlinear ["\s]*single-phase solids/);
  assert.match(cpu, /constitutive stress is assembled ["\s]*directly as total stress/);
  assert.ok(existsSync(`${benchmarkRoot}manzari_dafalias_modified_euler_synthetic.in`));
  assert.ok(existsSync(`${benchmarkRoot}manzari_dafalias_forward_euler_synthetic.in`));
  assert.ok(existsSync(`${benchmarkRoot}manzari_dafalias_forward_euler_demand_synthetic.in`));
  assert.ok(existsSync(`${benchmarkRoot}manzari_dafalias_seeded_onset_synthetic.in`));
  assert.ok(existsSync(`${benchmarkRoot}manzari_dafalias_pressure_floor_synthetic.in`));
  assert.ok(existsSync(`${benchmarkRoot}opensees_manzari_dafalias_synthetic.tcl`));
  assert.ok(existsSync(`${benchmarkRoot}compare_opensees_manzari_dafalias.py`));
  assert.ok(existsSync(`${benchmarkRoot}opensees_manzari_dafalias_drained_triaxial.tcl`));
  assert.ok(existsSync(`${benchmarkRoot}opensees_manzari_dafalias_undrained_cyclic_shear.tcl`));
  assert.ok(existsSync(`${benchmarkRoot}opensees_manzari_dafalias_sspbrickup_undrained_triaxial.tcl`));
  assert.ok(existsSync(`${benchmarkRoot}generate_hercules_manzari_dafalias_path.py`));
  assert.ok(existsSync(`${benchmarkRoot}compare_opensees_manzari_dafalias_canonical.py`));
  assert.ok(existsSync(`${benchmarkRoot}compare_opensees_manzari_dafalias_brickup.py`));
  assert.ok(existsSync(`${benchmarkRoot}plot_opensees_manzari_dafalias_brickup.py`));
  assert.match(generator, /material_point_initial_stress/);
  assert.match(generator, /--recovery-mode/);
  assert.match(generator, /--integration-scheme/);
  assert.match(generator, /--forward-euler-substeps/);
  assert.match(generator, /--compare-rtol/);
  assert.match(generator, /forward_euler[\s\S]*1\.0e-5/);
  assert.match(generator, /f"\{args\.integration_scheme\}"/);
  assert.match(generator, /forward_euler[\s\S]*--recovery-mode=strict/);
  assert.match(generator, /manzari_dafalias_recovery = \{args\.recovery_mode\}/);
  assert.match(generator, /manzari_dafalias_forward_euler_substeps = /);
  assert.doesNotMatch(fixedForwardEulerFixture, /manzari_dafalias_forward_euler_substeps/);
  assert.match(demandForwardEulerFixture, /manzari_dafalias_forward_euler_substeps = demand/);
  assert.match(materialPointSuite, /validate_md_forward_euler_substeps fixed/);
  assert.match(materialPointSuite, /validate_md_forward_euler_substeps demand/);
  assert.match(materialPointSuite, /\$substeps_column != cap/);
  assert.match(materialPointSuite, /\$substeps_column < cap/);
  assert.match(generator, /"alpha_in": tuple/);
  assert.match(generator, /"md_drained_triaxial": \("9", "complete"\)/);
  assert.match(generator, /"md_sspbrickup_undrained_triaxial": \("4", "complete"\)/);
  assert.match(canonicalComparator, /"drained": \("md_drained_triaxial", "monotonic", "9"\)/);
  assert.match(cpu, /material_point_initial_alpha_in/);
  assert.match(canonicalComparator, /def replay_rows/);
  assert.match(drainedDriver, /equalDOF 7 3 1 2/);
  assert.match(drainedDriver, /set stress \[eleResponse 1 stress3D6\]/);
  assert.match(drainedDriver, /set strain \[eleResponse 1 strain3D6\]/);
  assert.match(drainedDriver, /set state \[eleResponse 1 state\]/);
  assert.doesNotMatch(drainedDriver, /eleResponse 1 material 1/);
  assert.match(drainedDriver, /1 0 1 1\.0e-4 1\.0e-7/);
  assert.match(drainedDriver, /0\.5\*\[lindex \$strain 3\]/);
  assert.match(drainedDriver, /OpenSees drained deformation is not affine/);
  assert.match(drainedDriver, /test NormUnbalance 1\.0e-3 100 0/);
  assert.match(drainedDriver, /OpenSees drained confinement residual exceeds/);
  assert.match(generator, /maximum_confinement_error_kpa > 1\.0e-3/);
  assert.match(undrainedDriver, /set stress \[eleResponse 1 material 1 stress\]/);
  assert.match(brickupDriver, /element SSPbrickUP[\s\S]*-lumped/);
  assert.doesNotMatch(brickupDriver, /^\s*equalDOF/m);
  assert.match(brickupDriver, /sp 5 3 1\.0[\s\S]*sp 8 3 1\.0/);
  assert.match(brickupDriver, /wipeAnalysis[\s\S]*constraints Penalty 1\.0e18 1\.0e18[\s\S]*integrator Newmark/);
  assert.match(brickupDriver, /test NormDispIncr 1\.0e-9 100 0/);
  assert.match(brickupDriver, /set pressureStabilization 1\.5e-9/);
  assert.match(brickupDriver, /\$voidRatio \$pressureStabilization -lumped/);
  assert.match(brickupDriver, /set pressureValues[\s\S]*nodeVel \$node 4/);
  assert.match(brickupDriver, /OpenSees SSPbrickUP sealed-storage residual/);
  assert.match(cpu, /material_point_brickup_fluid_bulk_modulus_pa/);
  assert.match(kernel, /poroelastic_gpu_uniform_element_history/);
  assert.match(kernel, /kernelPoroelasticPressureForce<<<1, 1>>>/);
  assert.match(kernel, /kernelPoroelasticPressureUpdate<<<1, 8>>>/);
  assert.match(brickupComparator, /maximum_pore_pressure_difference_kpa/);
  assert.match(brickupComparator, /PASS: all requested SSPbrickUP\/BrickUP acceptance/);
  assert.match(brickupPlot, /hercules_pore_pressure_pa/);
  assert.match(brickupPlot, /hercules_total_sig_xx_pa/);
  assert.match(brickupPlot, /with_suffix\("\.png"\)/);
  assert.match(brickupPlot, /with_suffix\("\.pdf"\)/);
});

test("wires PressureDependMultiYield03 through parser, shared core, CUDA state, and fixture", () => {
  const sourceRoot = fileURLToPath(new URL("../../quake/forward/", import.meta.url));
  const benchmarkRoot = fileURLToPath(
    new URL("../../benchmark_inputs/nonlinear_material_point/", import.meta.url),
  );
  if (!existsSync(sourceRoot)) return;

  const cpu = readFileSync(`${sourceRoot}nonlinear.cc`, "utf8");
  const gpu = readFileSync(`${sourceRoot}nonlinear/nonlinear_gpu_implementation.cuh`, "utf8");
  const header = readFileSync(`${sourceRoot}nonlinear.h`, "utf8");
  const core = readFileSync(`${sourceRoot}nonlinear_pdmy03.h`, "utf8");
  const psolve = readFileSync(`${sourceRoot}psolve.cc`, "utf8");
  const pressureKernel = readFileSync(`${sourceRoot}kernel.cu`, "utf8");
  const suite = readFileSync(`${benchmarkRoot}run_suite.sh`, "utf8");
  const pdmyBrickupComparator = readFileSync(
    `${benchmarkRoot}compare_opensees_pdmy03_brickup.py`,
    "utf8",
  );
  const pdmyBrickupPlot = readFileSync(
    `${benchmarkRoot}plot_opensees_pdmy03_brickup.py`,
    "utf8",
  );
  const pdmyBrickupGenerator = readFileSync(
    `${benchmarkRoot}generate_hercules_pdmy03_brickup_path.py`,
    "utf8",
  );
  const pdmyBrickupTriaxialDriver = readFileSync(
    `${benchmarkRoot}opensees_pdmy03_sspbrickup_undrained_triaxial.tcl`,
    "utf8",
  );
  const pdmyBrickupCyclicDriver = readFileSync(
    `${benchmarkRoot}opensees_pdmy03_sspbrickup_undrained_cyclic_shear.tcl`,
    "utf8",
  );
  const cyclicFixture = readFileSync(
    `${benchmarkRoot}pdmy03_cyclic_shear_synthetic.in`,
    "utf8",
  );
  const onsetFixture = readFileSync(
    `${benchmarkRoot}pdmy03_zero_tension_onset_synthetic.in`,
    "utf8",
  );

  assert.match(header, /PRESSURE_DEPEND_MULTIYIELD_03/);
  assert.match(header, /#define PDMY03_MAX_SURFACES 40/);
  assert.match(header, /pdmy03_constants_t \*pdmy03_materials/);
  assert.match(cpu, /PDMY03_PROPERTY_COUNT 21/);
  assert.match(cpu, /MATERIAL_MODEL_IS_RIVA_SAND_FAMILY\(theMaterialModel\) \? 13 : 16/);
  assert.match(cpu, /nearest_property_value/);
  assert.match(cpu, /pdmy03_prepare_constants/);
  assert.match(cpu, /PDMY03 initialization complete:/);
  assert.match(gpu, /state\.pdmy03_materials\+nl_eindex/);
  assert.match(gpu, /PDMY03 material properties/);
  assert.match(
    cpu,
    /MPI_Reduce\(&myNonlinElementsCount,[\s\S]{0,160}MPI_SUM, 0, comm_solver\)/,
  );
  assert.match(
    cpu,
    /MPI_Reduce\(local_phase_seconds,[\s\S]{0,160}MPI_MAX, 0, comm_solver\)/,
  );
  assert.match(cpu, /myNonlinElementClassification\[eindex\] = selected/);
  assert.match(
    cpu,
    /#pragma omp parallel for schedule\(static\) if\(myMesh->lenum >= 100000\)/,
  );
  assert.match(
    cpu,
    /#pragma omp parallel for if\(parallel_pdmy03_initialization\)[\s\S]{0,1000}const int32_t eindex = myNonlinElementsMapping\[nl_eindex\]/,
  );
  assert.match(
    cpu,
    /get_topo_nonlin_flag\(\)[\s\S]{0,100}get_topo_meth\(\) == VT[\s\S]{0,100}isTopoElement/,
  );
  assert.match(core, /pdmy03_requested_substeps/);
  assert.match(cpu, /material_point_result_names\[i\]/);
  assert.match(pdmyBrickupComparator, /PDMY03_DIAGNOSTIC_FIELDS/);
  assert.match(pdmyBrickupComparator, /legacy_header/);
  assert.match(pdmyBrickupPlot, /TRIAXIAL_TEST/);
  assert.match(pdmyBrickupPlot, /CYCLIC_TEST/);
  assert.match(pdmyBrickupPlot, /Effective stress path/);
  assert.match(pdmyBrickupComparator, /DRIVER_VERSION = "2"/);
  assert.match(pdmyBrickupGenerator, /DRIVER_VERSION = "2"/);
  assert.match(pdmyBrickupGenerator, /"20 10\.0 0\.01 1\.013e5 1\.73e3"/);
  assert.doesNotMatch(pdmyBrickupGenerator, /"20 1\.0e4 0\.01/);
  assert.match(pdmyBrickupTriaxialDriver, /phaseTransform mType ca cb/);
  assert.match(pdmyBrickupTriaxialDriver, /26\.0 0/);
  assert.match(pdmyBrickupTriaxialDriver, /undrained_triaxial 2/);
  assert.match(pdmyBrickupCyclicDriver, /phaseTransform mType ca cb/);
  assert.match(pdmyBrickupCyclicDriver, /26\.0 0/);
  assert.match(pdmyBrickupCyclicDriver, /undrained_cyclic_shear 2/);
  assert.match(core, /c->failure_ratio\*cone\/[\s\n ]*\(double\)c->surface_count/);
  assert.match(core, /requested > maximum_substeps \? maximum_substeps : requested/);
  assert.match(core, /pdmy03_update_ppz/);
  assert.match(core, /pdmy03_update_active_surface/);
  assert.match(core, /if \(c == 0\.0\) return 0\.0/);
  assert.match(core, /return low;/);
  assert.match(core, /pdmy03_mean\(\*stress\) > 0\.0/);
  assert.match(core, /state->maximum_pressure = 0\.0/);
  assert.match(core, /surface <= active && surface < c->surface_count/);
  assert.match(core, /trial, current, centers/);
  assert.match(gpu, /pdmy03_surface_stride/);
  assert.match(gpu, /PDMY03 history allocation/);
  assert.match(gpu, /nl_record_pdmy03_error/);
  assert.match(gpu, /nonlinear_pdmy03_host_history/);
  assert.ok(existsSync(`${benchmarkRoot}pdmy03_cyclic_shear_synthetic.in`));
  assert.ok(existsSync(`${benchmarkRoot}pdmy03_zero_tension_onset_synthetic.in`));
  assert.ok(existsSync(`${benchmarkRoot}opensees_pdmy03_cyclic_shear.tcl`));
  assert.ok(existsSync(`${benchmarkRoot}compare_opensees_pdmy03.py`));
  assert.ok(existsSync(`${benchmarkRoot}plot_opensees_pdmy03.py`));
  assert.ok(existsSync(`${benchmarkRoot}opensees_pdmy03_sspbrickup_undrained_triaxial.tcl`));
  assert.ok(existsSync(`${benchmarkRoot}opensees_pdmy03_sspbrickup_undrained_cyclic_shear.tcl`));
  assert.ok(existsSync(`${benchmarkRoot}generate_hercules_pdmy03_brickup_path.py`));
  assert.ok(existsSync(`${benchmarkRoot}compare_opensees_pdmy03_brickup.py`));
  assert.ok(existsSync(`${benchmarkRoot}plot_opensees_pdmy03_brickup.py`));
  assert.match(cyclicFixture, /material_point_brickup_biot_alpha = 1\.0/);
  assert.match(onsetFixture, /material_point_brickup_biot_alpha = 0\.8/);
  assert.match(onsetFixture, /material_point_brickup_solid_bulk_modulus_pa/);
  assert.match(psolve, /nodalPoreVolume \/ storageM/);
  assert.match(cpu, /Nonlinear\/BrickUP map:/);
  assert.match(psolve, /geostatic_loading_time_sec=%g with zero permeability/);
  assert.match(pressureKernel, /const double mobility\[3\] = \{mobilityX, mobilityY, mobilityZ\}/);
  assert.match(pressureKernel, /double rhs = -qtDeltaU\[j\] - deltaT \* diffusionWeak\[j\]/);
  assert.match(suite, /pdmy03_cyclic_shear_synthetic\.in/);
});

test("wires both RIVA-Sand 02 interfaces without changing the original family", () => {
  const sourceRoot = fileURLToPath(new URL("../../quake/forward/", import.meta.url));
  const benchmarkRoot = fileURLToPath(
    new URL("../../benchmark_inputs/nonlinear_material_point/", import.meta.url),
  );
  if (!existsSync(sourceRoot)) return;
  const header = readFileSync(`${sourceRoot}nonlinear.h`, "utf8");
  const cpu = readFileSync(`${sourceRoot}nonlinear.cc`, "utf8");
  const gpu = readFileSync(`${sourceRoot}nonlinear/nonlinear_gpu_implementation.cuh`, "utf8");
  const point = readFileSync(
    `${sourceRoot}nonlinear/nonlinear_gpu_riva_sand_intermediate_bias_material_point.cu`, "utf8",
  );
  const suite = readFileSync(`${benchmarkRoot}run_suite.sh`, "utf8");
  const originalFixture = readFileSync(
    `${benchmarkRoot}riva_sand_intermediate_bias_research_synthetic.in`, "utf8",
  );
  const customFixture = readFileSync(`${benchmarkRoot}riva_sand_custom02_equivalence.in`, "utf8");
  assert.match(header, /RIVA_SAND, RIVA_SAND_CUSTOM,\s*RIVA_SAND_02, RIVA_SAND_CUSTOM_02\s*\} materialmodel_t/);
  assert.match(header, /#define MATERIAL_MODEL_IS_RIVA_SAND\(model\) \\\n\s*\(\(model\) == RIVA_SAND \|\| \(model\) == RIVA_SAND_CUSTOM\)/);
  assert.match(header, /#define MATERIAL_MODEL_IS_RIVA_SAND_02\(model\) \\\n\s*\(\(model\) == RIVA_SAND_02 \|\| \(model\) == RIVA_SAND_CUSTOM_02\)/);
  assert.match(header, /#define MATERIAL_MODEL_IS_RIVA_SAND_CUSTOM\(model\) \\\n\s*\(\(model\) == RIVA_SAND_CUSTOM \|\| \(model\) == RIVA_SAND_CUSTOM_02\)/);
  for (const name of ["RIVASand02", "RIVASandCustom02"]) {
    assert.match(cpu, new RegExp(`strcasecmp\\(material_model, "${name}"\\)`));
    assert.match(cpu, new RegExp(`strcasecmp\\(name, "${name}"\\)`));
  }
  assert.doesNotMatch(cpu, /RIVASandIntermediateBiasResearch|RIVA_SAND_INTERMEDIATE_BIAS_RESEARCH/);
  const customGuard = "MATERIAL_MODEL_IS_RIVA_SAND_CUSTOM(theMaterialModel)";
  const rows = cpu.slice(cpu.indexOf("static riva_material_parameters_t riva_material_at_row"),
    cpu.indexOf("/* Categorical material-table fields"));
  assert.equal(rows.split(customGuard).length - 1, 0);
  assert.equal((rows.match(/riva_material_set_G0/g) ?? []).length, 2);
  for (const field of ["M", "kd", "h", "m", "zeta", "e_max", "e_min", "Q", "R", "n_G"])
    assert.equal((rows.match(new RegExp(`material\\.${field} =`, "g")) ?? []).length, 2);
  assert.match(cpu, /const size_t material_count =\s*\(size_t\)myNonlinElementsCount/);
  assert.doesNotMatch(cpu, /myNonlinSolver->riva_materials\[0\] =/);
  assert.match(cpu, /case RIVA_SAND_02:\s*case RIVA_SAND_CUSTOM_02:/);
  assert.match(cpu, /MATERIAL_MODEL_IS_RIVA_SAND_FAMILY\(material_model\)\s*\? MATERIAL_POINT_RIVA_CONSTANT_COUNT/);
  const load = cpu.slice(cpu.indexOf("static void material_point_load_constants"),
    cpu.indexOf("static void material_point_load_constants") + 2500);
  assert.match(load, /riva_material_set_G0\(&reference, riva_material, values\[1\]\)/);
  ["M", "kd", "h", "m", "zeta", "e_max", "e_min", "Q", "R", "n_G"].forEach((field, i) => {
    assert.match(load, new RegExp(`riva_material->${field} = values\\[${i + 2}\\]`));
  });
  assert.match(gpu, /point_material=state\.riva_materials\[nl_eindex\]/);
  assert.match(gpu, /MATERIAL_MODEL_IS_RIVA_SAND_02\(config->material_model\) &&\s*count > 0\) \{\s*const size_t material_count =\s*\(size_t\)count;/);
  assert.match(gpu, /MATERIAL_MODEL_IS_RIVA_SAND_FAMILY\(g_config\.material_model\) &&\s*!g_riva_geostatic\.reported/);
  assert.match(gpu, /MATERIAL_MODEL_IS_RIVA_SAND_02\(g_config\.material_model\)\) \{\s*nonlinear_gpu_launch_riva_sand_intermediate_bias_research/);
  assert.match(gpu, /g_active && MATERIAL_MODEL_IS_RIVA_SAND_02\(g_config\.material_model\)\)\s*flags \|= NL_CHECKPOINT_RIVA_INTERMEDIATE_BIAS_RESEARCH_STATE/);
  assert.match(point, /!MATERIAL_MODEL_IS_RIVA_SAND_02\(config->material_model\)/);
  assert.match(customFixture, /material_model = RIVASandCustom02/);
  assert.match(customFixture, /Dr G0 M kd h m zeta e_max e_min Q R n_G/);
  assert.equal(customFixture.slice(customFixture.indexOf("material_point_initial_strain")),
    originalFixture.slice(originalFixture.indexOf("material_point_initial_strain")));
  assert.match(suite, /riva_sand_custom02_equivalence\.in/);
  assert.match(suite, /validate_riva_custom02_equivalence/);
  assert.match(suite, /RIVASandCustom02 reference row differs from RIVASand02/);
});

test("wires frozen and custom RIVA-Sand through CUDA, BrickUP, and checkpoints", () => {
  const sourceRoot = fileURLToPath(new URL("../../quake/forward/", import.meta.url));
  const benchmarkRoot = fileURLToPath(
    new URL("../../benchmark_inputs/nonlinear_material_point/", import.meta.url),
  );
  if (!existsSync(sourceRoot)) return;

  const cpu = readFileSync(`${sourceRoot}nonlinear.cc`, "utf8");
  const header = readFileSync(`${sourceRoot}nonlinear.h`, "utf8");
  const gpuHeader = readFileSync(`${sourceRoot}nonlinear_gpu.h`, "utf8");
  const core = readFileSync(`${sourceRoot}nonlinear_riva_sand.h`, "utf8");
  const researchCore = readFileSync(
    `${sourceRoot}nonlinear_riva_sand_intermediate_bias_research.h`, "utf8",
  );
  const gpu = readFileSync(`${sourceRoot}nonlinear/nonlinear_gpu_implementation.cuh`, "utf8");
  const gpuUnit = readFileSync(`${sourceRoot}nonlinear/nonlinear_gpu_riva_sand.cu`, "utf8");
  const researchMaterialPointUnit = readFileSync(
    `${sourceRoot}nonlinear/nonlinear_gpu_riva_sand_intermediate_bias_material_point.cu`,
    "utf8",
  );
  const checkpoint = readFileSync(`${sourceRoot}io_checkpoint.cc`, "utf8");
  const makefile = readFileSync(`${sourceRoot}Makefile`, "utf8");
  const fixture = readFileSync(`${benchmarkRoot}riva_sand_brickup_synthetic.in`, "utf8");
  const customFixture = readFileSync(`${benchmarkRoot}riva_sand_custom_equivalence.in`, "utf8");
  const lowBiasFixture = readFileSync(`${benchmarkRoot}riva_sand_low_bias_synthetic.in`, "utf8");
  const geostaticFixture = readFileSync(
    `${benchmarkRoot}riva_sand_ramped_geostatic_admission_synthetic.in`,
    "utf8",
  );
  const researchFixture = readFileSync(
    `${benchmarkRoot}riva_sand_intermediate_bias_research_synthetic.in`,
    "utf8",
  );
  const suite = readFileSync(`${benchmarkRoot}run_suite.sh`, "utf8");
  const nativeDriver = readFileSync(`${benchmarkRoot}riva_sand_native_driver.cc`, "utf8");
  const nativeStateTest = readFileSync(`${benchmarkRoot}riva_sand_native_state_test.cc`, "utf8");
  const geostaticSpecification = readFileSync(
    `${benchmarkRoot}RIVA_SAND_GEOSTATIC_ADMISSION.md`,
    "utf8",
  );
  const frozenSha = "9585f0c155c9444885c5115e7753a1a5d97c783e2c685938a49a65435c9e8f83";
  const commonConstants = header.match(
    /typedef struct nlconstants_t \{([\s\S]*?)\} nlconstants_t;/,
  );

  assert.match(header, /RIVA_SAND_CUSTOM/);
  assert.match(header, /RIVA_SAND_02/);
  assert.match(
    header,
    /#define MATERIAL_MODEL_IS_RIVA_SAND\(model\)[\s\S]*RIVA_SAND_CUSTOM\)[\s\S]*#define MATERIAL_MODEL_IS_RIVA_SAND_02/,
  );
  assert.match(header, /MATERIAL_MODEL_IS_RIVA_SAND_FAMILY/);
  assert.match(header, /riva_material_parameters_t/);
  assert.ok(commonConstants);
  assert.doesNotMatch(
    commonConstants[1],
    /manzari_dafalias_constants_t|pdmy03_constants_t|riva_material|riva_void_ratio/,
  );
  assert.match(header, /typedef struct nlsolver_t \{[\s\S]*md_materials[\s\S]*pdmy03_materials[\s\S]*riva_materials[\s\S]*riva_void_ratios/);
  assert.match(cpu, /RIVASand[\s\S]*RIVA-Sand/);
  assert.match(cpu, /RIVASandCustom[\s\S]*RIVA-Sand-Custom/);
  assert.match(cpu, /RIVASand02/);
  assert.match(cpu, /RIVA_IB_PARAMETER_SHA256/);
  assert.doesNotMatch(cpu, /PJLiq|PJ[-_ ]?Liq/);
  assert.match(cpu, /MATERIAL_MODEL_IS_RIVA_SAND_FAMILY\(theMaterialModel\) \? 13 : 16/);
  assert.match(cpu, /relative density must be[\s\S]*dimensionless range \[0, 1\]/);
  assert.match(core, /riva_void_ratio_from_relative_density/);
  assert.match(core, /return material->e_max-relative_density\*\(material->e_max-material->e_min\)/);
  assert.match(core, /riva_reference_material_parameters/);
  assert.match(core, /riva_material_parameters_valid/);
  assert.match(core, /riva_update_material/);
  assert.match(cpu, /"%s\/BrickUP coupling/);
  assert.match(cpu, /RIVA-Sand evolves effective ["\s]*skeleton stress/);
  assert.match(cpu, /BrickUP independently solves fluid ["\s]*mass balance and applies sigma=sigma'-alpha\*p\*I/);
  assert.match(cpu, /RIVA-Sand effective-pressure ["\s]*loss is not copied into the pressure DOF/);
  assert.match(core, new RegExp(frozenSha));
  assert.match(core, /#define RIVA_STATE_VALUE_COUNT 93/);
  assert.match(core, /riva_geostatic_cone_ratio/);
  assert.match(core, /#define RIVA_KERNEL_REVISION 2u/);
  assert.match(core, /p\.p_residual=0\.0/);
  assert.match(core, /riva_initialize_geostatic_material/);
  assert.match(core, /riva_admit_geostatic_state/);
  assert.match(core, /riva_regularize_deviator/);
  assert.match(core, /bias_onset=0\.25\*p->bias_reversible_volume_reference_bias/);
  assert.match(core, /bias_factor \*= gate\*gate\*gate/);
  assert.match(core, /riva_host_reversal/);
  assert.match(core, /riva_forward_euler/);
  assert.match(core, /riva_state_values/);
  assert.match(gpuHeader, /NL_RIVA_STATE_VALUE_COUNT 93/);
  assert.match(gpuHeader, /NL_RIVA_INTERMEDIATE_BIAS_STATE_VALUE_COUNT 139/);
  assert.match(researchCore, /#define RIVA_IB_KERNEL_REVISION 5u/);
  assert.equal((cpu.match(/riva_read_reversal_type\(fp, "rivasand02_reversaltype"/g) ?? []).length, 2);
  assert.match(cpu, /theRivaReversalType=3/);
  assert.equal((cpu.match(/riva_read_bias_volume\(fp, "rivasand02_biasvolume"/g) ?? []).length, 2);
  assert.match(cpu, /theRivaBiasVolume=1/);
  assert.match(cpu, /int_message\[16\] = theRivaBiasVolume/);
  assert.match(cpu, /theRivaBiasVolume\s*= int_message\[16\]/);
  assert.match(cpu, /gpuConfig\.riva_bias_volume = theRivaBiasVolume/);
  assert.match(cpu, /config\.riva_bias_volume = riva_bias_volume/);
  assert.match(gpu, /riva_ib_reference_parameters\(1000\.0,config->riva_bias_volume\)/);
  assert.match(researchMaterialPointUnit, /riva_ib_reference_parameters\(1000\.0,config->riva_bias_volume\)/);
  assert.match(researchCore, /bias_volume!=1/);
  assert.match(researchCore, /field_bias_mean_correction_enabled=bias_volume==2/);
  assert.match(researchCore, /values\[i\+\+\]=s->field_bias_mean_activity/);
  assert.match(gpu, /RIVASAND02_FAILURE_V2/);
  assert.match(readFileSync(`${sourceRoot}nonlinear_riva_sand_input.h`, "utf8"), /riva_read_02_selector\(fp,key,is_02,out,1,0,2\)/);

  assert.match(cpu, /int_message\[15\] = theRivaReversalType/);
  assert.match(cpu, /MPI_Bcast\(int_message,\s*17,/);
  assert.match(cpu, /theRivaReversalType\s*= int_message\[15\]/);
  assert.match(cpu, /gpuConfig\.riva_reversal_type = theRivaReversalType/);
  assert.match(cpu, /config\.riva_reversal_type = riva_reversal_type/);
  assert.match(gpuHeader, /int32_t riva_reversal_type/);
  assert.match(gpu, /&riva_info,cfg\.riva_reversal_type/);
  assert.match(researchMaterialPointUnit, /&state, &stress, &info, config->riva_reversal_type/);
  assert.match(researchCore, /riva_update_info_t \*info,int32_t reversal_type=3/);

  assert.match(researchCore, /riva_ib_host_stress_ratio_predictor/);
  assert.match(researchCore, /riva_ib_host_stress_ratio_reversal/);
  assert.match(researchCore, /reversal_type>=1 && reversal_type<=3/);
  assert.match(researchCore, /\(objective \|\| host_ratio\) && reversal && i==0/);
  assert.match(readFileSync(`${sourceRoot}nonlinear_riva_sand_input.h`, "utf8"), /riva_read_02_selector\(fp,key,is_02,out,3,1,3\)/);
  assert.match(researchCore, /#define RIVA_IB_PARAMETER_COUNT 257/);
  assert.match(researchCore, /#define RIVA_IB_LOGICAL_STATE_COUNT 64/);
  assert.match(researchCore, /#define RIVA_IB_STATE_VALUE_COUNT 139/);
  assert.match(researchCore, /16c8567fd2370c562899b9ec00cc2fc313be137f239cde5d0e22cc1cd3eee800/);
  assert.equal(
    (researchCore.match(/state\.ep_half_last=old->base\.ep_eq_since_reversal;/g) ?? []).length,
    2, "both actual research backbones must capture completed activity",
  );
  assert.match(researchCore, /state->ep_half_last=0\.0/);
  assert.match(researchCore, /s->ep_half_last=0\.0/);
  assert.match(researchCore, /values\[i\+\+\]=s->ep_half_last/);
  assert.match(researchCore, /target=riva_ib_mul_rn\(target,plastic_gate\)/);
  assert.match(researchCore, /s->ep_half_last\/p->bias_reversible_volume_ep_ref/);
  assert.match(researchCore, /#define RIVA_IB_CUSTOM_KERNEL_REVISION 5u/);
  assert.match(researchCore, /10\+3\*bias_volume\+reversal_type-1/);
  assert.match(gpu, /header\.riva_kernel_revision = riva_ib_checkpoint_revision\(\s*g_config\.material_model == RIVA_SAND_CUSTOM_02,\s*g_config\.riva_reversal_type,g_config\.riva_bias_volume\)/);
  assert.match(gpu, /!riva_ib_checkpoint_compatible\(header\.riva_kernel_revision,[\s\S]*?g_config\.material_model == RIVA_SAND_CUSTOM_02,\s*g_config\.riva_reversal_type,g_config\.riva_bias_volume\)/);
  assert.match(researchCore, /riva_ib_reference_relative_density\(\s*const riva_ib_parameters_t \*p,const riva_material_parameters_t \*m\)/);
  assert.match(researchCore, /if \(m->e_max==p->base\.e_max && m->e_min==p->base\.e_min\)\s*return p->reference_relative_density_value/);
  assert.equal((researchCore.match(/p->reference_relative_density_value/g) ?? []).length,
    1, "only the material-bound-aware helper may use the shared cached density");
  assert.match(researchCore, /state_bytes!=sizeof\(riva_ib_state_t\)/);
  assert.match(researchCore, /parameter_hash\[i\]!=expected\[i\]/);
  assert.match(
    gpuHeader,
    /nonlinear_gpu_riva_intermediate_bias_research_material_point_history/,
  );
  assert.match(gpu, /MATERIAL_MODEL_IS_RIVA_SAND\(g_config\.material_model\)[\s\S]*nonlinear_gpu_launch_riva_sand/);
  assert.match(gpu, /riva_initialize_geostatic_material/);
  assert.match(gpu, /"%s geostatic admission: initialized/);
  assert.match(gpu, /RIVA-Sand intermediate-bias research/);
  assert.match(gpu, /maximum admitted cone ratio/);
  assert.match(gpuHeader, /riva_geostatic_admitted/);
  assert.match(cpu, /RIVASand geostatic-admission mismatch/);
  assert.match(cpu, /RIVASand02 reversal-selector/);
  assert.match(gpu, /riva_update_material\(riva_parameters, riva_material/);
  assert.match(gpu, /host_state->riva_materials[\s\S]*host_state->riva_void_ratios/);
  assert.match(gpu, /local_material=%.1f MiB/);
  assert.match(gpu, /full_effective_stress/);
  assert.match(gpuUnit, /NL_GPU_RIVA_SAND_ONLY/);
  assert.match(makefile, /nonlinear\/nonlinear_gpu_riva_sand\.o/);
  assert.match(
    makefile,
    /nonlinear\/nonlinear_gpu_riva_sand_intermediate_bias_research\.o/,
  );
  assert.match(
    makefile,
    /nonlinear\/nonlinear_gpu_riva_sand_intermediate_bias_material_point\.o/,
  );
  assert.match(
    researchMaterialPointUnit,
    /riva_ib_state_values\(state,[\s\n ]*result->riva_intermediate_bias_state_values\)/,
  );
  assert.match(
    researchMaterialPointUnit,
    /#include "psolve\.h"[\s\S]*#include "nonlinear_gpu\.h"/,
  );
  assert.match(researchMaterialPointUnit, /riva_ib_zero_increment/);
  assert.match(researchMaterialPointUnit, /riva_ib_update_material/);
  assert.match(
    researchMaterialPointUnit,
    /nonlinear_riva_intermediate_bias_research_host_history/,
  );
  assert.match(fixture, /material_model = RIVASand/);
  assert.match(fixture, /0\.66296296296296309/);
  assert.match(fixture, /material_point_brickup = yes/);
  assert.match(fixture, /material_point_brickup_biot_alpha = 1\.0/);
  assert.match(suite, /riva_sand_brickup_synthetic\.in/);
  assert.match(lowBiasFixture, /material_model = RIVASand/);
  assert.match(lowBiasFixture, /-1\.8e4 -1\.8e4 -4\.0e4 1\.0e3 0 0/);
  assert.match(suite, /riva_sand_low_bias_synthetic\.in/);
  assert.match(geostaticFixture, /-2\.0e3 -2\.0e3 -4\.0e4/);
  assert.match(suite, /riva_sand_ramped_geostatic_admission_synthetic\.in/);
  assert.match(customFixture, /material_model = RIVASandCustom/);
  assert.match(customFixture, /Dr G0 M kd h m zeta e_max e_min Q R n_G/);
  assert.match(suite, /riva_sand_custom_equivalence\.in/);
  assert.match(suite, /RIVASandCustom reference row differs from frozen RIVA-Sand/);
  assert.match(
    researchFixture,
    /material_model = RIVASand02/,
  );
  assert.match(researchFixture, /material_point_steps = 18/);
  assert.match(researchFixture, /material_point_brickup = yes/);
  assert.match(
    suite,
    /riva_sand_intermediate_bias_research_synthetic\.in/,
  );
  assert.match(nativeDriver, /RIVA_STANDALONE/);
  assert.doesNotMatch(nativeDriver, /expected.*stress/i);
  assert.match(nativeStateTest, /frozen\/custom equivalence, zeta exposure/);
  assert.match(nativeStateTest, /relative-density conversion failed/);
  assert.match(nativeStateTest, /riva_geostatic_cone_ratio/);
  assert.match(nativeStateTest, /fluid mass balance before pressure is combined with[\s\S]*the[\s\S]*effective skeleton stress exactly once/);
  assert.match(geostaticSpecification, /rho_geo = \|\|alpha\|\|\/b/);
  assert.match(geostaticSpecification, /b_limit = max\(b_new, \|\|alpha_old\|\|\)/);
  assert.match(geostaticSpecification, /separate comparison field rather than a 94th oracle value/);

  assert.match(checkpoint, /nonlinear_gpu_checkpoint_write/);
  assert.match(checkpoint, /nonlinear_gpu_checkpoint_read/);
  assert.match(gpu, /g_nonlinear_checkpoint_magic/);
  assert.match(gpu, /pore_pressure_tm1/);
  assert.match(gpu, /pore_pressure_tm2/);
  assert.match(gpu, /g_nl\.riva_states/);
  assert.match(gpu, /RIVA_PARAMETER_SHA256/);
  assert.match(gpu, /riva_kernel_revision = RIVA_G0_CHECKPOINT_REVISION/);
  assert.match(gpu, /header\.riva_kernel_revision != RIVA_G0_CHECKPOINT_REVISION/);
  assert.match(gpu, /Checkpoint lacks a compatible BrickUP\/%s ["\s]*extension/);
});

test("ships the BAE virgin-onset regression fixture", () => {
  const benchmarkRoot = fileURLToPath(
    new URL("../../benchmark_inputs/nonlinear_material_point/", import.meta.url),
  );
  const suite = readFileSync(`${benchmarkRoot}run_suite.sh`, "utf8");
  const fixture = readFileSync(
    `${benchmarkRoot}vonmises_bae_virgin_onset_synthetic.in`,
    "utf8",
  );

  assert.match(suite, /vonmises_bae_virgin_onset_synthetic\.in/);
  assert.match(fixture, /Step 0 is purely volumetric/);
  assert.match(fixture, /crosses that floor/);
  assert.match(fixture, /another purely volumetric update/);
});

test("ships the OpenSees BAE characterization fixture and driver", () => {
  const benchmarkRoot = fileURLToPath(
    new URL("../../benchmark_inputs/nonlinear_material_point/", import.meta.url),
  );
  const driver = readFileSync(`${benchmarkRoot}compare_opensees_bae.py`, "utf8");
  const fixture = readFileSync(
    `${benchmarkRoot}vonmises_bae_opensees_synthetic.in`,
    "utf8",
  );
  const standalone = readFileSync(
    `${benchmarkRoot}opensees_bae_synthetic.tcl`,
    "utf8",
  );
  const readme = readFileSync(`${benchmarkRoot}README.md`, "utf8");

  assert.match(driver, /gamma_xy = 2\.0 \* strain\[3\]/);
  assert.match(driver, /hardening_scale = fixture\.psi0 \* fixture\.shear/);
  assert.match(driver, /element stdBrick/);
  assert.match(driver, /OPEN_STRESS_FIELDS[\s\S]*sigma_xy[\s\S]*sigma_yz[\s\S]*sigma_xz/);
  assert.match(driver, /--opensees-csv/);
  assert.match(fixture, /material_point_steps = 20/);
  assert.match(fixture, /repeated zero-strain point/);
  assert.match(standalone, /nDMaterial MultiaxialCyclicPlasticity/);
  assert.match(standalone, /element stdBrick/);
  assert.match(standalone, /set gammaHistory/);
  assert.match(standalone, /timeSeries Path/);
  assert.match(standalone, /sp 3 1 1\.0/);
  assert.match(standalone, /set actualGamma \[nodeDisp 3 1\]/);
  assert.match(standalone, /driver_version,domain_time/);
  assert.doesNotMatch(standalone, /integrator LoadControl 0\.0/);
  assert.doesNotMatch(standalone, /integrator DisplacementControl/);
  assert.match(standalone, /set stress \[eleResponse 1 material 1 stress\]/);
  assert.match(readme, /characterization comparison, not an equivalence proof/i);
});

test("uses sparse CUDA state gathering for nonlinear station output", () => {
  const sourceRoot = fileURLToPath(new URL("../../quake/forward/", import.meta.url));
  if (!existsSync(sourceRoot)) return;

  const psolve = readFileSync(`${sourceRoot}psolve.cc`, "utf8");
  const gpu = readFileSync(`${sourceRoot}nonlinear/nonlinear_gpu_implementation.cuh`, "utf8");
  const finish = psolve.slice(
    psolve.indexOf("static void solver_finish_nonlinear_gpu"),
    psolve.indexOf("static void solver_compute_force_planewaves"),
  );

  assert.match(finish, /nonlinear_gpu_sync_station_state\(solver\)/);
  assert.doesNotMatch(finish, /nonlinear_gpu_sync_host_state\(solver\)/);
  assert.match(gpu, /kernel_nonlinear_gather_station_state/);
  assert.match(gpu, /g_station\.count\s*\*\s*sizeof\(\*g_station\.host\)/);
  assert.match(gpu, /record->nonlinear_index/);
});

test("scopes legacy von Mises constants to Istanbul inputs", () => {
  const sourcePath = fileURLToPath(new URL("../../quake/forward/nonlinear.cc", import.meta.url));
  if (!existsSync(sourcePath)) return;

  const source = readFileSync(sourcePath, "utf8");
  const bae = source.slice(
    source.indexOf("case VONMISES_BAE:"),
    source.indexOf("case VONMISES_BAH:"),
  );
  const ep = source.slice(
    source.indexOf("case VONMISES_EP:"),
    source.indexOf("case DRUCKERPRAGER:"),
  );

  assert.match(bae, /if \( theIstanbulVelocityModel == YES \)/);
  assert.match(bae, /else \{\s*ecp->c\s*= get_cohesion\(elementVs\);/);
  assert.match(bae, /ecp->psi0 = interpolate_property_value\(elementVs, thePsi\);/);
  assert.match(bae, /ecp->m\s*= interpolate_property_value\(elementVs, theM\);/);
  assert.match(ep, /if \( theIstanbulVelocityModel == YES \)/);
  assert.match(ep, /else \{\s*ecp->c = get_cohesion\(elementVs\);/);
});

test("documents every supported source model and its auxiliary inputs", async () => {
  const response = await render("/docs/source-models");
  assert.equal(response.status, 200);

  const html = await response.text();
  for (const sourceType of ["point", "force", "base_force", "plane", "planewithkinks", "srfh"]) {
    assert.match(html, new RegExp(`>${sourceType}<`));
  }
  const sourceKeys = [
    "source_is_filtered", "threshold_frequency", "number_of_poles",
    "type_of_source", "source_function_type", "average_risetime_sec",
    "ricker_Ts", "ricker_Tp", "moment_magnitude", "moment_amplitude",
    "lonlat_or_cartesian", "hypocenter_lat_deg", "hypocenter_long_deg",
    "hypocenter_x", "hypocenter_y", "hypocenter_depth_m",
    "source_strike_deg", "source_dip_deg", "source_rake_deg",
    "force_x_n", "force_y_n", "force_z_n", "base_force_rho_kgm3",
    "base_force_vs_mps", "base_force_direction_x", "base_force_direction_y",
    "base_force_direction_z", "base_force_scale", "base_force_depth_m",
    "base_force_tolerance_m", "number_of_time_windows", "time_windows",
    "extended_cell_size_down_dip_m", "extended_isminimumedgeautomatic",
    "extended_depth_m", "extended_cells_along_strike",
    "extended_cells_down_dip", "extended_hypocenter_along_strike_m",
    "extended_hypocenter_down_dip_m", "extended_average_rupture_velocity",
    "extended_cell_size_along_strike_m", "extended_lat_deg",
    "extended_long_deg", "extended_strike_deg", "extended_dip_deg",
    "extended_minimum_edge_m", "extended_number_of_kinks",
    "extended_kinks", "number_of_point_sources",
  ];
  for (const key of sourceKeys) assert.match(html, new RegExp(key));
  assert.match(html, /forcefunction\.in/);
  assert.match(html, /baseforcefunction\.in/);
  assert.match(html, /slipfunction\.in/);
  assert.match(html, /hold the final sample/);
  assert.match(html, /Reducing simulation_delta_time_sec does not change the source-history sample interval/);
  assert.match(html, /Source-filter sampling and start time/);
  assert.match(html, /common slipfunction\.in interval, not simulation_delta_time_sec/);
  assert.match(html, /Changing the solver timestep must not change an SRFH filter/);
  assert.match(html, /there is no universal delay in seconds/);
  assert.match(html, /change the final applied force/);
  assert.match(html, /extended_number_of_kinks/);
  assert.match(html, /number_of_point_sources/);
  assert.match(html, /local solver X\/Y and Z positive downward/);
});

test("provides configuration and reader examples for output files", async () => {
  const response = await render("/docs/output-files");
  assert.equal(response.status, 200);

  const html = await response.text();
  assert.match(html, /np\.loadtxt/);
  assert.match(html, /np\.fromfile/);
  assert.match(html, /raw\.reshape/);
  for (const name of ["read_stations", "read_planes", "read_drm", "read_mesh"]) {
    assert.ok(html.includes(`script/${name}.py`));
  }
  assert.match(html, /script\/animations_mpi\.py/);
  assert.match(html, /compare-plane-outputs/);
  assert.match(html, /3,127,704,000 bytes/);
  assert.match(html, /NaN and infinity always fail/);
  assert.match(html, /equally truncated outputs/);
  assert.doesNotMatch(html, /(?:compare_plane_outputs|diagnose_plane_output|diagnose_acoustic_stations|check_acoustic_gravity_mode)\.py/);
  assert.match(html, /does not read parameters\.in or any other configuration file/);
  assert.match(html, /--nx 201 --ny 126 --dx 500 --dy 325 --dt 0\.02/);
  assert.match(html, /mpirun -np 2/);
  assert.match(html, /--map-by ppr:1:node --bind-to none --report-bindings/);
  assert.match(html, /--map-by ppr:1:node:PE=4 --bind-to core/);
  assert.match(html, /--workers 100/);
  assert.match(html, /Scheduler-enforced CPU sets still apply/);
  assert.match(html, /Workers are not individually pinned/);
  assert.match(html, /Mapping and binding do not request resources from the scheduler/);
  assert.match(html, /mpirun --version/);
  assert.match(html, /os\.sched_getaffinity\(0\)/);
  assert.match(html, /Do not transpose or mirror this matrix/);
  assert.match(html, /40\.625 km/);
  assert.match(html, /--check-only/);
  assert.match(html, /--first-time/);
  assert.match(html, /reported extrema remain unclipped/);
  assert.match(html, /--plot-upsample 4 --plot-interpolation linear/);
  assert.match(html, /--plot-upsample 1/);
  assert.match(html, /801 by 501 display points/);
  assert.match(html, /--plot-interpolation cubic/);
  assert.match(html, /keeping the same limit across frames/);
  assert.match(html, /closes its PyGMT session before deleting its temporary directory/);
  assert.match(html, /GMT_DATA_UPDATE_INTERVAL = off/);
  assert.match(html, /checkpoint\.out1.*checkpoint\.in/s);
  assert.match(html, /mesh_coordinates\.0/);
  assert.match(html, /station_mediums\.txt/);
  assert.match(html, /print_station_pore_pressure_ratios/);
  assert.match(html, /Ru\(-\)/);
  assert.match(html, /Δu\/σ′v0/);
  assert.match(html, /sigma_v0_eff_pa/);
  assert.match(html, /ramped geostatic phase/);
  assert.match(html, /Ru remains nan during the ramp/);
  assert.match(html, /committed vertical skeleton stress/);
  assert.match(html, /output_displacement_file/);
});

test("documents plane amplitude diagnostics without equating geometry with CPML membership", async () => {
  const response = await render("/docs/output-files");
  assert.equal(response.status, 200);
  const html = await response.text();
  assert.match(html, /plane,i,j,rank,element,x_m,y_m,z_m,medium,cpml/);
  assert.match(html, /The suffix is the rank, not the plane number/);
  assert.match(html, /contains no pressure, velocity or displacement histories/);
  assert.match(html, /only a CPML candidate region/);
  assert.match(html, /HERCULES_PLANE_METADATA=1/);
  assert.match(html, /plane_point_metadata\.RANK\.csv/);
  assert.match(html, /fresh output directories/);
  assert.match(html, /keep missing ownership unknown/);
  assert.match(html, /no per-step work/);
  assert.match(await (await render("/docs/configuration")).text(), /HERCULES_PLANE_METADATA/);
});

test("converts velocity potential to pressure before reconstructing acoustic planes", async () => {
  const response = await render("/docs/output-files");
  assert.equal(response.status, 200);
  const html = await response.text();
  assert.match(html, /p = -\(phi\[n\] - phi\[n-1\]\)\/dt/);
  assert.match(html, /conversion adds no GPU transfers/);
  assert.match(html, /old code gives -0\.6 m at 4 s/);
  assert.match(html, /not validation of a complete CUDA\/MPI/);
  assert.match(html, /inverse gradient stretch/);
  assert.match(html, /separate memory at each observation point/);
  assert.doesNotMatch(html, /routine has no CPML memory correction/);
  const models = await (await render("/docs/physics")).text();
  assert.match(models, /Acoustic free-surface geometry/);
  assert.match(models, /acoustic free surfaces follow local topography over water/);
  assert.match(models, /Without explicit seabed geometry/);
  assert.match(models, /missed shifted water surfaces/);
  assert.match(models, /MPI displacement sharing and hanging-node assignment/);
  assert.match(models, /Compare one and four ranks/);
  assert.match(models, /Acoustic interface geometry/);
  assert.match(models, /independent nodal surface memories/);
  assert.match(models, /full tangential x\/y product, including corners/);
  assert.match(models, /Earlier candidate and doubled-CPML runs contained invalid source forces/);
  assert.match(models, /nonfinite output cannot establish an operator instability/);

  const sourceRoot = fileURLToPath(new URL("../../quake/forward/", import.meta.url));
  if (!existsSync(sourceRoot)) return;
  const solver = readFileSync(`${sourceRoot}psolve.cc`, "utf8");
  assert.match(solver, /planes_update_acoustic_kinematics\([\s\S]*?Param\.theDeltaTSquared,\s*Param\.theAcousticVelocityPotential\)/);
  assert.match(solver, /Acoustic plane reconstruction:/);
  // The sparse plane gather still patches the raw current scalar; pointer
  // swaps retain its previous value instead of transferring another field.
  const sync = solver.slice(solver.indexOf("static void solver_sync_acoustic_plane_pressure_values(mysolver_t* solver) {"),
                            solver.indexOf("static void solver_init_pore_pressure_station_sync(mysolver_t* solver) {"));
  assert.equal((sync.match(/cudaMemcpyAsync\(/g) ?? []).length, 1);
  assert.doesNotMatch(sync, /pressureTm2Device/);
});

test("documents controlled acoustic pressure and motion diagnostics", async () => {
  const response = await render("/docs/output-files");
  const html = await response.text();
  assert.match(html, /acoustic-station-diagnostics/);
  assert.match(html, /print_station_pressures = yes/);
  assert.match(html, /output_stations_print_rate = 1/);
  assert.match(html, /Recheck material ownership after coarsening/);
  assert.match(html, /shared reconstruction error/);
  assert.match(html, /CPML remains active in both runs/);
  assert.match(html, /gravity-mode-diagnostics/);
  assert.match(html, /Investigation scripts and smoke-test inputs are retained locally/);
  assert.match(html, /without CPML, adaptive meshes, solid coupling or GPU\/MPI execution/);
  assert.match(html, /does not filter the source/);
});

test("documents structured DRM boxes and generic station binary output", async () => {
  const response = await render("/docs/output-files");
  const html = await response.text();
  for (const key of ["station_binary", "station_ascii_buffer_kib", "station_output_binary_file", "drm_output_batch_size",
    "drm_print_rate", "drm_boundary", "drm_edgesize", "drm_offset_x", "drm_offset_y",
    "number_output_drm_boxes", "output_drm_boxes_print_rate", "output_drm_boxes_directory",
    "output_drm_boxes"]) {
    assert.ok(html.includes(key), `missing DRM configuration key ${key}`);
  }
  assert.match(html, /input row order is preserved/i);
  assert.match(html, /including duplicate rows/);
  assert.match(html, /exterior front, back, left, right, bottom, then the interior/);
  assert.match(html, /shared edges\/corners are repeated/);
  assert.match(html, /60,406/);
  assert.match(html, /1,449,744 bytes/);
  assert.match(html, /top center/);
  assert.match(html, /length_x − 2h, length_y − 2h, depth − h/);
  assert.match(html, /clockwise from solver \+X/);
  assert.match(html, /drmdisplacements\.0/);
  assert.match(html, /drmdisplacements\.1/);
  assert.match(html, /drmcoordinates\.0/);
  assert.match(html, /drmcoordinates\.1/);
  assert.match(html, /X_m Y_m Z_m Vs_m_s Vp_m_s density_kg_m3 Qs Qp/);
  assert.match(html, /initial containing mesh element/);
  assert.match(html, /produce nan/);
  assert.match(html, /selected nominal P-wave Q/);
  assert.match(html, /Qp=168\.75/);
  assert.match(html, /inf means zero selected material loss/);
  assert.match(html, /No input-Qp, Qk or fallback-flag columns/);
  assert.match(html, /format version 2/);
  assert.match(html, /metadata.shape != \(npoints, 8\)/);
  assert.match(html, /vs, vp, density, qs, qp/);
  assert.match(html, /Data row i .* matches binary sample i/);
  assert.match(html, /17 significant digits/);
  assert.match(html, /once during setup after point ownership validation/);
  assert.match(html, /np\.loadtxt/);
  assert.doesNotMatch(html, /No extra face files or coordinate sidecars/);
  assert.match(html, /box_layout\(100, 100, 50, 1\)/);
  assert.match(html, /a box centered at \(0,0\) extends to negative X\/Y/i);
  assert.match(html, /cannot be fed directly to the existing part2 reader/);
  assert.match(html, /Each run truncates station_output_binary_file/);
  assert.match(html, /scalar-acoustic points stop setup/);
  assert.match(html, /does not refine the mesh/);
  assert.match(html, /shape=\(nt, npoints, 3\)/);
  assert.match(html, /station_binary defaults to no/);
  assert.match(html, /64 KiB buffer per locally owned station/);
  assert.match(html, /default 64, accepted range 1\.\.2097151/);
  assert.match(html, /validated but unused with station_binary=yes/);
  assert.match(html, /62\.5 MiB; setting 16 uses 15\.625 MiB/);
  assert.match(html, /Every requested sample and column is retained without filtering/);
  assert.match(html, /includes the final flush in the Print Stations and Solver I\/O timers/);
  assert.match(html, /interrupted or aborted run can lose the unflushed tail/);
  assert.match(html, /raw.reshape\(-1, nstations, 3\)/);
  assert.match(html, /No per-station text files or station metadata sidecar/);
  assert.match(html, /reconstructed scalar-acoustic motion/);
  assert.match(html, /OUTPUT_PATH EXPECTED_RANKS/);
  assert.match(html, /check MPI_COMM_WORLD before opening output/);
  assert.match(html, /did not establish multi-rank coverage/);
  assert.doesNotMatch(html, /drm_points_file|drm_output_file|CSV recorder/);
});

test("reads structured DRM faces with exterior dimensions and rejects degenerate interiors", () => {
  const reader = fileURLToPath(new URL("../../script/read_drm.py", import.meta.url));
  if (!existsSync(reader)) return;
  const dir = mkdtempSync(join(tmpdir(), "hercules-drm-box-reader-"));
  try {
    const file = join(dir, "drmdisplacements.0");
    // Analytically: exterior 131 + interior 63 = 194 grid samples for 4×6×3, h=1.
    const payload = new Float64Array(2 * 194 * 3);
    writeFileSync(file, Buffer.from(payload.buffer));
    const run = (dimensions) => spawnSync("python3", [reader, file, "--box", ...dimensions,
      "--dt", "0.01"], { encoding: "utf8" });
    const result = run(["4", "6", "3", "1"]);
    assert.equal(result.status, 0, result.stderr);
    assert.match(result.stdout, /layout=\(2, 194, 3\)/);
    assert.match(result.stdout, /exterior_front: shape=\(2, 7, 4, 3\), point_offset=0/);
    assert.match(result.stdout, /interior_front: shape=\(2, 5, 3, 3\), point_offset=131/);
    assert.match(result.stdout, /interior_bottom: shape=\(2, 3, 5, 3\), point_offset=179/);
    const bad = run(["2", "6", "3", "1"]);
    assert.notEqual(bad.status, 0);
    assert.match(bad.stderr, /exterior lengths must exceed/);
    assert.notEqual(run(["4", "6", "3", "0.7"]).status, 0);
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
});

test("reads complete flat DRM snapshots and rejects incomplete or nonfinite data", () => {
  const reader = fileURLToPath(new URL("../../script/read_drm.py", import.meta.url));
  if (!existsSync(reader)) return;
  const dir = mkdtempSync(join(tmpdir(), "hercules-drm-reader-"));
  try {
    const file = join(dir, "drm.bin");
    const payload = new Float64Array([0, 0, 0, 1, 2, 3, 4, -5, 6, 7, 8, -9]);
    writeFileSync(file, Buffer.from(payload.buffer));
    const run = () => spawnSync("python3", [reader, file, "--points", "2", "--dt", "0.02",
      "--first-time", "0.1"], { encoding: "utf8" });
    const result = run();
    assert.equal(result.status, 0, result.stderr);
    assert.match(result.stdout, /layout=\(2, 2, 3\)/);
    assert.match(result.stdout, /time range=0.1\.\.0.12 s/);
    assert.match(result.stdout, /7.0 8.0 9.0/);
    writeFileSync(file, Buffer.from(payload.buffer).subarray(0, 95));
    assert.match(run().stderr, /complete 48-byte snapshots/);
    payload[11] = Number.NaN;
    writeFileSync(file, Buffer.from(payload.buffer));
    assert.match(run().stderr, /snapshot 1, point 1 contains non-finite/);
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
});

test("reports the conventional BrickUP station pore-pressure ratio", () => {
  const sourceRoot = fileURLToPath(new URL("../../quake/forward/", import.meta.url));
  if (!existsSync(sourceRoot)) return;

  const psolve = readFileSync(`${sourceRoot}psolve.cc`, "utf8");
  const nonlinear = readFileSync(`${sourceRoot}nonlinear.cc`, "utf8");
  const header = readFileSync(`${sourceRoot}psolve.h`, "utf8");
  assert.match(psolve, /parsetext\(fp, "print_station_pore_pressure_ratios"/);
  assert.match(psolve, /print_station_pore_pressure_ratios=yes requires[\s\S]{0,120}include_poroelastic=yes[\s\S]{0,120}include_nonlinear_analysis=yes/);
  assert.match(psolve, /print_station_pore_pressure_ratios=yes requires an initial[\s\S]{0,240}approximate_geostatic_state=no with a positive geostatic loading/);
  assert.match(psolve, /hydrostaticBiotStress = Param\.thePoroelasticBiotAlpha \*[\s\S]{0,120}Param\.thePoroelasticFluidDensity \* 9\.81 \* waterDepth/);
  assert.match(psolve, /effectiveVerticalStress =[\s\S]{0,100}totalVerticalStress - hydrostaticBiotStress/);
  assert.match(psolve, /get_geostatic_total_time\(\) > 0\.0[\s\S]{0,80}return NAN/);
  assert.match(psolve, /\? pressure \/ sigmaV0[\s\S]{0,80}: NAN/);
  assert.match(psolve, /sigma_v0_eff_pa/);
  assert.match(psolve, /step == get_geostatic_final_step\(\)/);
  assert.match(psolve, /station_rewrite_mediums_file_from_current_state\(1\)/);
  assert.match(psolve, /return \(Param\.printStationPressures == YES\) \|\|[\s\S]{0,100}\(Param\.printStationPorePressureRatios == YES\)/);
  assert.match(nonlinear, /-myNonlinSolver->stresses\[nl_eindex\]\.qp\[4\]\.zz/);
  assert.match(header, /initial_vertical_effective_stress/);
});

test("documents independent BrickUP constraints and final station eligibility", async () => {
  const coupling = await (await render("/docs/physics")).text();
  assert.match(coupling, /include_acoustic=no still builds the pore-pressure map automatically/);
  assert.match(coupling, /VT cut-cell models keep geometric pressure constraints/);
  assert.match(coupling, /standalone extension still requires target-system CUDA\/MPI/);
  assert.match(coupling, /open_pore also builds actual solid traces with shared-geometric point projection/);
  const output = await (await render("/docs/output-files")).text();
  assert.match(output, /Station pressure eligibility is finalized after solver setup/);
  assert.match(output, /Elastic-only CPML receivers are labeled elastic with pressure undefined/);
  assert.match(output, /Their Pressure\(Pa\) and Ru\(-\) values are nan/);
  const reference = await (await render("/docs/parameters-coupled-physics")).text();
  assert.match(reference, /Coarse dry\/CPML cells do not constrain pore pressure/);
});

test("documents four-point open-pore hydraulic integration and its validation scope", async () => {
  const physics = await (await render("/docs/physics")).text();
  assert.match(physics, /four-point \(2x2 Gauss\) hydraulic integration/);
  assert.match(physics, /acoustic_interface_projection still controls mechanical traction\/motion only/);
  assert.match(physics, /Capacity and target share the same rule/);
  assert.match(physics, /166\.666… Pa/);
  assert.match(physics, /233\.333… Pa/);
  assert.match(physics, /Earlier Horizon open-pore benchmarks used centroid hydraulic integration/);
  assert.match(physics, /Target-system CUDA\/MPI propagation and performance checks remain required/);
  const reference = await (await render("/docs/parameters-coupled-physics")).text();
  assert.match(reference, /Open-pore hydraulic capacity and pressure targets always use four-point/);
  assert.match(reference, /same quadrature for boundary capacity/);
});

test("documents the basin-bottom grid format and renders a usable numeric example", async () => {
  const response = await render("/docs/input-files");
  assert.equal(response.status, 200);
  const html = await response.text();
  assert.match(html, /id="basin-configuration"/);
  assert.match(html, /basin_velocity_model\s*= yes/);
  assert.match(html, /basin_input_file\s*= input\/basin\.in/);
  assert.match(html, /basin_origin_latitude_deg/);
  assert.match(html, /basin_origin_longitude_deg/);
  assert.match(html, /X is north and Y is east/);
  assert.match(html, /same depth reference as solver Z/);
  assert.match(html, /Δx = Δy/);
  assert.match(html, /Write X fastest/);
  assert.match(html, /no column labels or comments/);
  assert.match(html, /Current maximum-edge limitation/);
  assert.match(html, /without clamping/);
  assert.match(html, /z &lt; bottom_depth\(x, y\)/);
  assert.match(html, /basin_material_type=acoustic_water/);
  assert.match(html, /benchmark_inputs\/bowl\/bowl_basin_250m\.in/);

  const formatSection = html.match(/<section[^>]*id="basin-input-format"[^>]*>([\s\S]*?)<\/section>/)?.[1];
  assert.ok(formatSection, "basin format section is rendered");
  const example = formatSection.match(/<pre><code>([\s\S]*?)<\/code><\/pre>/)?.[1];
  assert.ok(example, "basin example is a copyable code block");
  const [header, ...rows] = example.trim().split("\n").map((line) => line.trim().split(/\s+/).map(Number));
  const [nx, ny, nz] = header;
  assert.deepEqual(header, [3, 3, 1]);
  assert.equal(rows.length, nx * ny * nz);
  for (const [index, row] of rows.entries()) {
    assert.equal(row.length, 3);
    assert.ok(row.every(Number.isFinite));
    assert.equal(row[0], (index % nx) * 1000, "X varies fastest from zero");
    assert.equal(row[1], Math.floor(index / nx) * 1000, "Y uses the same spacing");
    assert.equal(row[2], index === 4 ? 500 : 0, "only the centre has nonzero bottom depth");
  }

  const parameterResponse = await render("/docs/parameters-material-mesh");
  assert.equal(parameterResponse.status, 200);
  const parameters = await parameterResponse.text();
  assert.match(parameters, /Basin-bottom text grid: Nx Ny 1/);
  assert.match(parameters, /Basin-bottom grid format for the example and edge limitation/);
});

test("provides complete examples for 3-D material and topography inputs", async () => {
  const response = await render("/docs/input-files");
  assert.equal(response.status, 200);

  const html = await response.text();
  assert.match(html, /3D_velocity_model_directory/);
  assert.match(html, /rho\.in.*multiplies by 1000/s);
  assert.match(html, /qs\.in and qp\.in/);
  assert.match(html, /rho\.in.*qs\.in.*qp\.in.*index\.in/s);
  assert.match(html, /provide qs\.in and qp\.in as a matched pair/i);
  assert.match(html, /Solid-only interpolation excludes water placeholders/);
  assert.match(html, /Supplied Q is independent of the Vs minimum/);
  assert.match(html, /10% solid Qs=40 plus 90% water Qs=0 gives solid Qs=40/);
  assert.match(html, /Acoustic placeholders may be zero, but must be finite and nonnegative/);
  assert.match(html, /CSV is only a preparation format|CSV directly/);
  assert.match(html, /profile_id = j \* nx \+ i/);
  assert.match(html, /y as the outer loop and x as the inner loop/i);
  assert.match(html, /index\.in contains starts only/);
  assert.match(html, /index\.in = \[0, 4, 8, 11\]/);
  assert.match(html, /\[8, 11\)/);
  assert.match(html, /unequal layer counts: 4, 4, 3, and 4/);
  assert.match(html, /piecewise constant with depth/);
  assert.match(html, /bilinearly interpolates Vs, Vp, density/);
  assert.match(html, /3-D material depths with full and sqd topography/);
  assert.match(html, /d_i = d \+ h_i - h_query/);
  assert.match(html, /actual topography grid at the solver query point/);
  assert.match(html, /every column is queried at Z, with no column terrain shift/);
  assert.match(html, /sqd queries depths 30 and 130 m, and full queries Z=1030 m/);
  assert.match(html, /there is no global minimum-depth rejection/);
  assert.match(html, /does not interpolate corresponding layers at equal burial depth/);
  assert.match(html, /Imported element materials retain those old values/);
  assert.match(html, /not new GPU\/MPI wave-propagation benchmarks/);
  assert.match(html, /script\/convert_3d_velocity_model\.py/);
  for (const name of ["validate_3d_velocity_model", "generate_topography", "validate_topography"]) {
    assert.ok(html.includes(`script/${name}.py`));
  }
  assert.match(html, /csv\.DictReader/);
  assert.match(html, /--output-dir/);
  assert.match(html, /both qs and qp/);
  assert.match(html, /complete rectangular x-y grid/);
  assert.match(html, /index\.in does not define one ordered profile/);
  assert.match(html, /topography\.in — 3 × 3 grid/);
  assert.match(html, /flat_index = x_index \* n_east \+ y_index/);
  assert.match(html, /one bilinear topography grid for mesh air\/solid classification/);
  assert.match(html, /Regenerate topographic meshes.*mesh_input_mode=0/);
  assert.match(html, /region_depth_deep_m controls the domain extent, not this surface datum/);
  assert.match(html, /np\.meshgrid/);
  assert.match(html, /values\.reshape\(nx, ny\)/);
  assert.match(html, /perimeter range/);
  assert.match(html, /script\/taper_topography\.py/);
  assert.match(html, /--flat-margin 1500/);
  assert.match(html, /--pml-thickness.*does not configure Hercules/s);
  assert.match(html, /max \|z_boundary - z_ref\| = 0/);
  assert.match(html, /cpml_3d_material_extension/);
});

test("documents the required full/sqd material-coordinate convention", async () => {
  const html = await (await render("/docs/parameters-topography-buildings")).text();
  const row = html.split("</tr>").find((part) => part.includes(">type_of_etree</td>"));
  assert.ok(row);
  assert.match(row, /Required with topography; sqd \| full; no default/);
  assert.match(row, /sqd subtracts the local surface; full retains fixed model Z/);
  assert.match(row, /actual query-point topography/);
  assert.match(row, /Both are in metres, positive downward/);
  assert.match(row, /Changing modes requires matching material depths/);
});

test("converts a row-wise CSV into Hercules 3-D model files", () => {
  const repositoryRoot = fileURLToPath(new URL("../../", import.meta.url));
  const script = join(repositoryRoot, "script", "convert_3d_velocity_model.py");
  assert.ok(existsSync(script), "the documented converter exists under script/");

  const temporaryRoot = mkdtempSync(join(tmpdir(), "hercules-3d-model-"));
  try {
    const csvPath = join(temporaryRoot, "model.csv");
    const outputDir = join(temporaryRoot, "materialfiles");
    writeFileSync(
      csvPath,
      [
        "x,y,z,vs,vp,rho",
        "500,500,1000,900,1800,2.2",
        "0,0,1000,800,1600,2.1",
        "500,0,0,550,1100,2.0",
        "0,500,0,600,1200,2.0",
        "500,500,0,650,1300,2.0",
        "0,0,0,500,1000,2.0",
        "500,0,1000,850,1700,2.1",
        "0,500,1000,875,1750,2.1",
      ].join("\n"),
      "utf8",
    );

    const conversion = spawnSync(
      "python3",
      [script, csvPath, "--output-dir", outputDir],
      { encoding: "utf8" },
    );
    assert.equal(conversion.status, 0, conversion.stderr);
    assert.match(conversion.stdout, /wrote 2 x 2 profiles and 8 samples/);
    assert.equal(readFileSync(join(outputDir, "x.in"), "utf8"), "2\n0\n500\n");
    assert.equal(readFileSync(join(outputDir, "y.in"), "utf8"), "2\n0\n500\n");
    assert.equal(readFileSync(join(outputDir, "index.in"), "utf8"), "4\n0\n2\n4\n6\n");
    assert.equal(
      readFileSync(join(outputDir, "z.in"), "utf8"),
      "8\n0\n1000\n0\n1000\n0\n1000\n0\n1000\n",
    );
    assert.ok(existsSync(join(outputDir, "vs.in")));
    assert.ok(existsSync(join(outputDir, "vp.in")));
    assert.ok(existsSync(join(outputDir, "rho.in")));
    assert.equal(existsSync(join(outputDir, "qs.in")), false);
    assert.equal(existsSync(join(outputDir, "qp.in")), false);
  } finally {
    rmSync(temporaryRoot, { recursive: true, force: true });
  }
});

test("tapers rectangular topography in solver X-major order", () => {
  const repositoryRoot = fileURLToPath(new URL("../../", import.meta.url));
  const script = join(repositoryRoot, "script", "taper_topography.py");
  const temporaryRoot = mkdtempSync(join(tmpdir(), "hercules-topography-"));
  try {
    const source = join(temporaryRoot, "topography.in");
    const output = join(temporaryRoot, "tapered.in");
    const values = Array.from({ length: 63 }, (_, k) => 100 * Math.floor(k / 9) + k % 9);
    writeFileSync(source, [100, ...values].join("\n") + "\n");
    const result = spawnSync("python3", [
      script, source, "--out", output, "--lx", "600", "--ly", "800",
      "--buffer", "250", "--flat-margin", "50", "--pml-thickness", "25", "--z-ref", "0",
    ], { encoding: "utf8" });
    assert.equal(result.status, 0, result.stderr);
    const [spacing, ...tapered] = readFileSync(output, "utf8").trim().split(/\s+/).map(Number);
    assert.equal(spacing, 100);
    assert.equal(tapered.length, 63);
    assert.equal(tapered[3 * 9 + 4], 304, "deep interior retains its elevation");
    assert.equal(tapered[1 * 9 + 4], 15.2304, "north ramp has the analytical raised-cosine value");
    for (let i = 0; i < 7; i++) {
      assert.equal(tapered[i * 9], 0);
      assert.equal(tapered[i * 9 + 8], 0);
    }
    for (let j = 0; j < 9; j++) {
      assert.equal(tapered[j], 0);
      assert.equal(tapered[6 * 9 + j], 0);
    }
  } finally {
    rmSync(temporaryRoot, { recursive: true, force: true });
  }
});

test("preserves explicit G0 and caches velocity-derived RIVA-Sand stiffness", () => {
  const root = fileURLToPath(new URL("../../quake/forward/", import.meta.url));
  const cpu = readFileSync(root + "nonlinear.cc", "utf8");
  const core = readFileSync(root + "nonlinear_riva_sand.h", "utf8");
  const gpu = readFileSync(root + "nonlinear/nonlinear_gpu_implementation.cuh", "utf8");
  assert.match(cpu, /#define RIVA_PROPERTY_COUNT 12/);
  assert.match(cpu, /RIVA_PROPERTY_DR = 0,\s*RIVA_PROPERTY_G0,/);
  assert.match(cpu, /#define MATERIAL_POINT_RIVA_CONSTANT_COUNT 12/);
  assert.match(core, /RIVA_REFERENCE_G0 483\.48301127222084/);
  assert.match(core, /material->E_ref\*pow\(ratio,material->n_G\)/);
  assert.match(core, /isfinite\(material->E_ref\) && material->E_ref>0\.0/);
  assert.match(core, /G0 == RIVA_REFERENCE_G0 \? p->E_ref/);
  assert.match(gpu, /header\.riva_kernel_revision = RIVA_G0_CHECKPOINT_REVISION/);
  assert.match(gpu, /header\.riva_kernel_revision != RIVA_G0_CHECKPOINT_REVISION/);
  const input = readFileSync(root + "nonlinear_riva_sand_input.h", "utf8");
  const research = readFileSync(root + "nonlinear_riva_sand_intermediate_bias_research.h", "utf8");
  assert.match(cpu, /riva_read_material_rows\(fp, "material_properties_list"/);
  assert.match(cpu, /riva_read_material_rows\(fp,"material_point_constants"/);
  assert.match(input, /count!=1\+with_vs && count!=2\+with_vs/);
  assert.match(input, /automatic>=0 && automatic!=mode/);
  assert.match(cpu, /riva_target_shear\[nl_eindex\] = automatic \? mu : 0\.0/);
  assert.match(gpu, /riva_fit_initial_shear\(riva_parameters,&point_material/);
  assert.match(gpu, /riva_ib_fit_initial_shear\(riva_parameters,&point_material/);
  assert.equal((gpu.match(/point_material\.E_ref=state\.riva_E_refs\[point\]/g) ?? []).length, 2);
  assert.match(core, /RIVA_G0_CHECKPOINT_REVISION 4u/);
  assert.match(research, /10\+3\*bias_volume\+reversal_type-1/);
  for (const verb of ["read", "write"])
    assert.match(gpu, new RegExp(`nonlinear_checkpoint_${verb}_device\\(\\s*fp,g_nl\\.riva_E_refs`));
});


test("documents complete CPML products and the plain-CPML default", async () => {
  const html = await (await render("/docs/physics")).text();
  assert.match(html, /Complete second-order operator/);
  assert.match(html, /cpml_mpml_cross_fraction is now 0/);
  assert.match(html, /cpml_acoustic_sponge_fraction is an optional experimental acoustic bulk absorber/);
  assert.match(html, /defaults to 0\.3; set 0 to disable it/);
  assert.match(html, /checkpoint reads and writes/);
  const sourceRoot = fileURLToPath(new URL("../../quake/forward/", import.meta.url));
  if (!existsSync(sourceRoot)) return;
  const psolve = readFileSync(`${sourceRoot}psolve.cc`, "utf8");
  assert.match(psolve, /theCpmlMpmlCrossFraction = 0\.0/);
  assert.match(psolve, /theCpmlAcousticSpongeFraction = 0\.3/);
  assert.match(psolve, /parsetext\(fp, "cpml_acoustic_sponge_fraction"/);
  assert.match(psolve, /!isfinite\(cpml_acoustic_sponge_fraction\)/);
  assert.match(psolve, /if \(use_checkpoint \|\| checkpointing_rate > 0\)/);
});

test("keeps requested CPML defaults and alpha Hz units synchronized with the parser", async () => {
  const html = await (await render("/docs/parameters-damping-cpml")).text();
  for (const [key, value] of [
    ["cpml_reflection_target", "0.1"],
    ["cpml_kappa_max", "10"],
    ["cpml_acoustic_sponge_fraction", "0.3"],
  ]) {
    const row = html.split("</tr>").find((part) => part.includes(`>${key}</td>`));
    assert.ok(row?.includes(`Default: ${value};`), `${key} default`);
  }
  const alpha = html.split("</tr>").find((part) => part.includes(">cpml_alpha_max_hz</td>"));
  assert.match(alpha, /Default: simulation_wave_max_freq_hz \/ 2/);
  assert.match(alpha, /Omitted or nonpositive values/);
  assert.match(alpha, /positive numeric value overrides/);
  const source = new URL("../../quake/forward/psolve.cc", import.meta.url);
  if (!existsSync(source)) return;
  const solver = readFileSync(source, "utf8");
  assert.match(solver, /\.theCpmlReflectionTarget = 0\.1,/);
  assert.match(solver, /\.theCpmlKappaMax = 10\.0,/);
  assert.match(solver, /\.theCpmlAlphaMaxHz = -1\.0,/);
  assert.match(solver, /\.theCpmlAcousticSpongeFraction = 0\.3,/);
  assert.match(solver, /g_cpmlBuild\.alpha_max = \(Param\.theCpmlAlphaMaxHz > 0\.0\)\s*\? \(2\.0 \* M_PI \* Param\.theCpmlAlphaMaxHz\)\s*: \(M_PI \* Param\.theFreq\)/);
  assert.match(solver, /auto: simulation_wave_max_freq_hz \/ 2\)\\n",\s*0\.5 \* Param\.theFreq/);
  assert.match(solver, /\.theAcousticInterfaceMethod = ACOUSTIC_INTERFACE_FIELDS,/);
  assert.match(solver, /acoustic_interface_method_explicit =\s*\(parsetext\(fp, "acoustic_interface_method"/);
  assert.match(solver, /if \(!acoustic_interface_method_explicit &&\s*\(Param\.includeAcoustic != YES \|\| Param\.acousticInterfaceCoupling != YES \|\|\s*Param\.acousticPressureFormulation != YES\)\)\s*\{\s*Param\.theAcousticInterfaceMethod = ACOUSTIC_INTERFACE_ORIGINAL;/);
});

test("explains CPML material extension in model coordinates and mesh regeneration", async () => {
  const html = await (await render("/docs/physics")).text();
  assert.match(html, /Material extension clamps in absolute model coordinates/);
  assert.match(html, /depth below the surface at the clamped X\/Y location/);
  assert.match(html, /Building surface offsets/);
  assert.match(html, /input depth of 7500 m/);
  assert.match(html, /depth 7000 m below that surface/);
  assert.match(html, /Rebuild or recorrect an imported mesh/);
  assert.match(html, /corrected topographic bottom\/corner extension/);
});


test("documents CPML automatic thickness excluding air", async () => {
  const physics = await (await render("/docs/physics")).text();
  assert.match(physics, /largest non-air element edge/);
  assert.match(physics, /Air elements \(negative Vp\) are excluded/);
  assert.match(physics, /global fallback/);
  assert.match(physics, /not an exact cell count everywhere/);
  const parameters = await (await render("/docs/parameters-damping-cpml")).text();
  const thickness = parameters.split("</tr>").find((row) => row.includes("cpml_thickness_m"));
  const layers = parameters.split("</tr>").find((row) => row.includes(">cpml_layers</td>"));
  assert.ok(thickness);
  assert.ok(layers);
  assert.match(thickness, /largest non-air element edge/);
  assert.match(thickness, /Air is also excluded from the global fallback/);
  assert.match(layers, /finer regions contain more than this many cells/);
});


test("documents independent domain/material rotations and output frame conversion", async () => {
  const html = await (await render("/docs/input-files")).text();
  for (const phrase of ["Rotated domains and material grids", "3D_velocity_model_azimuth_deg",
    "coordinate_frame.txt", "CPML extension clamps in solver coordinates", "Rebuild persisted meshes",
    "true north", "vector_to_utm", "solver-aligned material grid"])
    assert.ok(html.includes(phrase), phrase);
  const params = await (await render("/docs/parameters-material-mesh")).text();
  const materialAngle = params.split("</tr>").find((row) => row.includes(">3D_velocity_model_azimuth_deg</td>"));
  assert.ok(materialAngle);
  assert.match(materialAngle, /Default: region_azimuth_leftface_deg; finite/);
  assert.match(materialAngle, /Explicit values override this default: 0/);
  assert.match(html, /If 3D_velocity_model_azimuth_deg is omitted, it inherits region_azimuth_leftface_deg/);
  assert.match(html, /Omit the next line for a material grid aligned with the domain/);

  const source = new URL("../../quake/forward/psolve.cc", import.meta.url);
  if (existsSync(source)) {
    const solver = readFileSync(source, "utf8");
    assert.match(solver, /Param\.the3DVelocityModelAzimuth = region_azimuth_leftface_deg;\s*const int material_azimuth_status =\s*hu_config_get_double_opt\(fp, "3D_velocity_model_azimuth_deg"/);
    assert.match(solver, /material_azimuth_status == 1 && Param\.the3DVelocityModelAzimuth != 0\.0/);
  }
});


test("documents large 3-D material counts and per-rank memory", async () => {
  const response = await render("/docs/input-files");
  assert.equal(response.status, 200);
  const html = await response.text();
  assert.match(html, /sample counts and profile start offsets are signed 64-bit/);
  assert.match(html, /replicated in host memory on every solver rank/);
  assert.match(html, /65\.3 GiB per solver rank/);
});


test("documents deferred CPML partitioning and the limits of GPU memory snapshots", async () => {
  const html = await (await render("/docs/physics")).text();
  for (const phrase of ["Completed-mesh partitioning and GPU memory",
    "ordinary partitions are still used throughout progressive construction",
    "local acoustic-interface matching", "equal-rank unified run",
    "preserves Q and geostatic data", "minimum free memory",
    "only a subset of the required storage", "later lazy allocations",
    "No diagnostics are added inside the timestep loop"])
    assert.ok(html.includes(phrase), phrase);
  const parameters = await (await render("/docs/parameters-damping-cpml")).text();
  const weight = parameters.split("</tr>").find((row) => row.includes(">cpml_partition_weight</td>"));
  assert.ok(weight);
  assert.match(weight, /Default: 1; &gt;= 1/);
  assert.match(weight, /Construction and mesh-only runs use ordinary partitions/);
  assert.match(weight, /not a GPU-memory cap/);
});
