import assert from "node:assert/strict";
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
