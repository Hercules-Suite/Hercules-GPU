import assert from "node:assert/strict";
import { existsSync, readFileSync } from "node:fs";
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
  assert.match(html, /one.*parameters\.in.*file/i);
  assert.match(html, /Search documentation/);
  assert.doesNotMatch(html, /codex-preview|react-loading-skeleton/);
});

test("server-renders a structured documentation page", async () => {
  const response = await render("/docs/getting-started");
  assert.equal(response.status, 200);

  const html = await response.text();
  assert.match(html, /Run your first Hercules simulation/);
  assert.match(html, /Before you begin/);
  assert.match(html, /Current Hercules GPU code path/);
});

test("documents the current single parameter-file interface", async () => {
  const response = await render("/docs/parameter-file");
  assert.equal(response.status, 200);

  const html = await response.text();
  assert.match(html, /The parameters\.in file/);
  assert.match(html, /psolve parameters\.in/);
  assert.match(html, /not the current psolve interface/);
  assert.match(html, /all 206 unique keys currently consumed from parameters\.in/);
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

test("documents every parameters.in key parsed by the current source tree", () => {
  const referencePath = fileURLToPath(new URL("../app/docs/parameter-reference.ts", import.meta.url));
  const reference = readFileSync(referencePath, "utf8");
  const documented = [...reference.matchAll(/\bp\("([^"]+)"/g)].map((match) => match[1]);
  assert.equal(new Set(documented).size, 206);
  assert.equal(documented.length, 206, "parameter definitions must not be duplicated");

  const sourceRoot = fileURLToPath(new URL("../../quake/forward/", import.meta.url));
  if (!existsSync(sourceRoot)) return;

  const parserFiles = [
    "psolve.cu",
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
  assert.match(html, /output_displacement_file/);
});
