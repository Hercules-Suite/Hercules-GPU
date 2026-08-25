import assert from "node:assert/strict";
import { existsSync, mkdtempSync, readdirSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { spawnSync } from "node:child_process";
import { fileURLToPath } from "node:url";
import test from "node:test";

async function render(path = "/") {
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

test("server-renders the Hercules documentation home", async () => {
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
  assert.match(html, /src="\/hercules-octree-mesh-mark\.png\?v=graded-centered-256"/);
  assert.doesNotMatch(html, /_vinext\/image\?url=.*hercules-octree-mesh-mark/);
  assert.ok(existsSync(join(fileURLToPath(new URL("../public/", import.meta.url)), "hercules-octree-mesh-mark.png")));
  assert.doesNotMatch(html, /codex-preview|react-loading-skeleton/);
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

  const assetsDirectory = fileURLToPath(new URL("../dist/client/assets/", import.meta.url));
  const clientJavaScript = readdirSync(assetsDirectory)
    .filter((name) => name.endsWith(".js"))
    .map((name) => readFileSync(join(assetsDirectory, name), "utf8"))
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
  assert.match(html, /all 207 unique keys currently consumed from parameters\.in/);
  assert.match(html, /Source-specific controls.*source\.in/);
});

test("renders the categorized parameter reference", async () => {
  const response = await render("/docs/parameters-coupled-physics");
  assert.equal(response.status, 200);

  const html = await response.text();
  assert.match(html, /Acoustic and poroelastic parameters/);
  assert.match(html, /acoustic_interface_face_filter/);
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

test("documents every parameters.in key parsed by the current source tree", () => {
  const referencePath = fileURLToPath(new URL("../app/docs/parameter-reference.ts", import.meta.url));
  const reference = readFileSync(referencePath, "utf8");
  const documented = [...reference.matchAll(/\bp\("([^"]+)"/g)].map((match) => match[1]);
  assert.equal(new Set(documented).size, 228);
  assert.equal(documented.length, 228, "parameter definitions must not be duplicated");

  const sourceRoot = fileURLToPath(new URL("../../quake/forward/", import.meta.url));
  if (!existsSync(sourceRoot)) return;

  const parserFiles = [
    "psolve.cc",
    "topography.cc",
    "nonlinear.cc",
    "buildings.cc",
    "drm.cc",
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
  const parserPattern = /(?:\b(?:hu_)?parsetext(?:_v)?|\bhu_config_get_[A-Za-z0-9_]+|\bparsedarray2?)\s*\(\s*[^,]+,\s*"([^"]+)"/gs;
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
  assert.match(physics, /Qs=C\+alpha×Vs_km_s\^beta/);
  assert.match(physics, /Baseline CPML example/);
  assert.match(physics, /cpml_reflection_target = 1\.0e-3/);
  assert.match(physics, /cpml_pml_material_damping = rayleigh/);
  assert.match(physics, /benchmark_inputs\/taper_topography\.py/);

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
  assert.match(parameters, /RIVASand: N × 2/);
  assert.match(parameters, /RIVASandCustom: N × 12/);
  assert.match(parameters, /Vs_limit Dr M kd h m zeta e_max e_min Q R n_G/);
  assert.match(parameters, /This is RIVA-Sand, not the experimental pre-name V8.1 h_p0 extension/);
  assert.match(parameters, /9585f0c155c9444885c5115e7753a1a5d97c783e2c685938a49a65435c9e8f83/);
  assert.match(parameters, /surface_count is an integer from 2 through 40/);
  assert.match(parameters, /Vs_limit G0 e_init Mc c lambda_c/);
  assert.match(parameters, /maximum absolute yield-surface residual/);
  assert.match(parameters, /documented MohrCoulomb, DruckerPrager, ManzariDafalias, PDMY03, RIVASand, and RIVASandCustom/);

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
  assert.match(gpu, /twenty-fixture suite/);
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
  assert.match(cpu, /property_columns\s*=\s*[\s\S]*MANZARI_DAFALIAS[\s\S]*\? 17 :[\s\S]*PRESSURE_DEPEND_MULTIYIELD_03[\s\S]*\? 22 :[\s\S]*RIVA_SAND_CUSTOM[\s\S]*\? 12 :[\s\S]*RIVA_SAND[\s\S]*\? 2 : 16/);
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
  assert.match(cpu, /PRESSURE_DEPEND_MULTIYIELD_03[\s\S]*\? 22 :[\s\S]*RIVA_SAND_CUSTOM[\s\S]*\? 12 :[\s\S]*RIVA_SAND[\s\S]*\? 2 : 16/);
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
  const gpu = readFileSync(`${sourceRoot}nonlinear/nonlinear_gpu_implementation.cuh`, "utf8");
  const gpuUnit = readFileSync(`${sourceRoot}nonlinear/nonlinear_gpu_riva_sand.cu`, "utf8");
  const checkpoint = readFileSync(`${sourceRoot}io_checkpoint.cc`, "utf8");
  const makefile = readFileSync(`${sourceRoot}Makefile`, "utf8");
  const fixture = readFileSync(`${benchmarkRoot}riva_sand_brickup_synthetic.in`, "utf8");
  const customFixture = readFileSync(`${benchmarkRoot}riva_sand_custom_equivalence.in`, "utf8");
  const lowBiasFixture = readFileSync(`${benchmarkRoot}riva_sand_low_bias_synthetic.in`, "utf8");
  const geostaticFixture = readFileSync(
    `${benchmarkRoot}riva_sand_ramped_geostatic_admission_synthetic.in`,
    "utf8",
  );
  const suite = readFileSync(`${benchmarkRoot}run_suite.sh`, "utf8");
  const nativeDriver = readFileSync(`${benchmarkRoot}riva_sand_native_driver.cc`, "utf8");
  const nativeStateTest = readFileSync(`${benchmarkRoot}riva_sand_native_state_test.cc`, "utf8");
  const frozenSha = "9585f0c155c9444885c5115e7753a1a5d97c783e2c685938a49a65435c9e8f83";
  const commonConstants = header.match(
    /typedef struct nlconstants_t \{([\s\S]*?)\} nlconstants_t;/,
  );

  assert.match(header, /RIVA_SAND_CUSTOM/);
  assert.match(header, /riva_material_parameters_t/);
  assert.ok(commonConstants);
  assert.doesNotMatch(
    commonConstants[1],
    /manzari_dafalias_constants_t|pdmy03_constants_t|riva_material|riva_void_ratio/,
  );
  assert.match(header, /typedef struct nlsolver_t \{[\s\S]*md_materials[\s\S]*pdmy03_materials[\s\S]*riva_materials[\s\S]*riva_void_ratios/);
  assert.match(cpu, /RIVASand[\s\S]*RIVA-Sand/);
  assert.match(cpu, /RIVASandCustom[\s\S]*RIVA-Sand-Custom/);
  assert.doesNotMatch(cpu, /PJLiq|PJ[-_ ]?Liq/);
  assert.match(cpu, /RIVA_SAND_CUSTOM[\s\S]*\? 12 :[\s\S]*RIVA_SAND[\s\S]*\? 2 : 16/);
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
  assert.match(gpu, /MATERIAL_MODEL_IS_RIVA_SAND\(g_config\.material_model\)[\s\S]*nonlinear_gpu_launch_riva_sand/);
  assert.match(gpu, /riva_initialize_geostatic_material/);
  assert.match(gpu, /riva_update_material\(riva_parameters, riva_material/);
  assert.match(gpu, /host_state->riva_materials[\s\S]*host_state->riva_void_ratios/);
  assert.match(gpu, /local_material=%.1f MiB/);
  assert.match(gpu, /full_effective_stress/);
  assert.match(gpuUnit, /NL_GPU_RIVA_SAND_ONLY/);
  assert.match(makefile, /nonlinear\/nonlinear_gpu_riva_sand\.o/);
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
  assert.match(customFixture, /Dr M kd h m zeta e_max e_min Q R n_G/);
  assert.match(suite, /riva_sand_custom_equivalence\.in/);
  assert.match(suite, /RIVASandCustom reference row differs from frozen RIVA-Sand/);
  assert.match(nativeDriver, /RIVA_STANDALONE/);
  assert.doesNotMatch(nativeDriver, /expected.*stress/i);
  assert.match(nativeStateTest, /frozen\/custom equivalence, zeta exposure/);
  assert.match(nativeStateTest, /relative-density conversion failed/);
  assert.match(nativeStateTest, /fluid mass balance before pressure is combined with[\s\S]*the[\s\S]*effective skeleton stress exactly once/);

  assert.match(checkpoint, /nonlinear_gpu_checkpoint_write/);
  assert.match(checkpoint, /nonlinear_gpu_checkpoint_read/);
  assert.match(gpu, /g_nonlinear_checkpoint_magic/);
  assert.match(gpu, /pore_pressure_tm1/);
  assert.match(gpu, /pore_pressure_tm2/);
  assert.match(gpu, /g_nl\.riva_states/);
  assert.match(gpu, /RIVA_PARAMETER_SHA256/);
  assert.match(gpu, /riva_kernel_revision = RIVA_KERNEL_REVISION/);
  assert.match(gpu, /header\.riva_kernel_revision != RIVA_KERNEL_REVISION/);
  assert.match(gpu, /Checkpoint lacks a compatible BrickUP\/RIVA-Sand ["\s]*extension/);
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
  assert.match(html, /extended_number_of_kinks/);
  assert.match(html, /number_of_point_sources/);
  assert.match(html, /x north, y east, and z positive downward/);
});

test("provides configuration and reader examples for output files", async () => {
  const response = await render("/docs/output-files");
  assert.equal(response.status, 200);

  const html = await response.text();
  assert.match(html, /np\.loadtxt/);
  assert.match(html, /np\.fromfile/);
  assert.match(html, /raw\.reshape/);
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
  assert.match(psolve, /station_rewrite_mediums_file_from_current_state\(\)/);
  assert.match(psolve, /return \(Param\.printStationPressures == YES\) \|\|[\s\S]{0,100}\(Param\.printStationPorePressureRatios == YES\)/);
  assert.match(nonlinear, /-myNonlinSolver->stresses\[nl_eindex\]\.qp\[4\]\.zz/);
  assert.match(header, /initial_vertical_effective_stress/);
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
  assert.match(html, /CSV is only a preparation format|CSV directly/);
  assert.match(html, /profile_id = j \* nx \+ i/);
  assert.match(html, /y as the outer loop and x as the inner loop/i);
  assert.match(html, /index\.in contains starts only/);
  assert.match(html, /index\.in = \[0, 4, 8, 11\]/);
  assert.match(html, /\[8, 11\)/);
  assert.match(html, /unequal layer counts: 4, 4, 3, and 4/);
  assert.match(html, /piecewise constant with depth/);
  assert.match(html, /bilinearly interpolates Vs, Vp, density/);
  assert.match(html, /benchmark_inputs\/convert_3d_velocity_model\.py/);
  assert.match(html, /csv\.DictReader/);
  assert.match(html, /--output-dir/);
  assert.match(html, /both qs and qp/);
  assert.match(html, /complete rectangular x-y grid/);
  assert.match(html, /index\.in does not define one ordered profile/);
  assert.match(html, /topography\.in — 3 × 3 grid/);
  assert.match(html, /flat_index = x_index \* n_east \+ y_index/);
  assert.match(html, /np\.meshgrid/);
  assert.match(html, /values\.reshape\(nx, ny\)/);
  assert.match(html, /perimeter range/);
  assert.match(html, /benchmark_inputs\/taper_topography\.py/);
  assert.match(html, /--flat-margin 1500/);
  assert.match(html, /--pml-thickness.*does not configure Hercules/s);
  assert.match(html, /max \|z_boundary - z_ref\| = 0/);
  assert.match(html, /cpml_3d_material_extension/);
});

test("converts a row-wise CSV into Hercules 3-D model files", () => {
  const repositoryRoot = fileURLToPath(new URL("../../", import.meta.url));
  const script = join(repositoryRoot, "benchmark_inputs", "convert_3d_velocity_model.py");
  if (!existsSync(script)) return;

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
