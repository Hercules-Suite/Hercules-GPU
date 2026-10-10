# Hercules Documentation

Technical documentation for the Hercules adaptive finite-element seismic wave
simulation code. The site covers installation, parameter and input files,
GPU/MPI execution, linear and nonlinear/plastic physical models, outputs,
validation, and troubleshooting.

The site source is maintained in `docs-site/` alongside the Hercules solver so
that implementation and documentation changes can be reviewed together. The
repository-level documentation workflow builds and tests the site without
compiling the scientific code.

## Quick start

Requires Node.js 22.13 or newer.

```bash
pnpm install
pnpm run dev
pnpm run build
```

## GitHub Pages

The same source also exports a static website for GitHub Pages, including the
current layout, browser-side search, theme switch, equations, and parameter
tables. The existing `build`, `start`, and `test` commands continue to use the
Sites/Vinext Worker build; `.openai/hosting.json` still identifies the university
Sites deployment.

From `docs-site/`, build and verify the static export with:

```bash
pnpm install --frozen-lockfile
pnpm run test:pages
```

`build:pages` runs the installed Next.js compiler in static-export mode and
replaces the generated `out/` directory. `test:pages` first builds that output,
then runs the existing rendered-content and parser audits against the exported
HTML, plus checks for missing pages, broken internal links, images, stylesheets,
JavaScript assets, and incorrect repository URL prefixes. No server or Sites
credentials are needed to serve the exported files. Keep `out/` and `.next/`
untracked. Both build paths require Node.js 22.13 or newer and the existing pnpm
lockfile; neither compiles the solver.

The default public address is `https://zwyll.github.io/Hercules_GPU/`. These
build-time controls allow a different hosting location:

| Control | Default for `build:pages` | Meaning |
| --- | --- | --- |
| `DOCS_STATIC_EXPORT` | `1`, set by the command | Select static metadata and Next.js export. Leave unset for the Sites build. |
| `NEXT_PUBLIC_DOCS_BASE_PATH` | `/Hercules_GPU` | URL prefix, beginning with `/` and without a trailing slash; use an empty string for a domain-root site. Next.js prefixes navigation links; the site prefixes public image URLs. |
| `DOCS_SITE_URL` | `https://zwyll.github.io/Hercules_GPU/` | Absolute public site URL used for social-preview metadata. Set this together with the URL prefix when changing host or repository. |

For example, to export for a custom domain at its root:

```bash
NEXT_PUBLIC_DOCS_BASE_PATH='' DOCS_SITE_URL=https://docs.example.org/ pnpm run test:pages
```

The documentation workflow builds and audits both targets on pull requests and
updates to `main`, including while the repository is private. Publication is
disabled until all of these activation steps are complete:

1. Make the repository public when ready. Changing its visibility is a separate
   maintainer action.
2. In **Settings → Pages → Build and deployment**, choose **GitHub Actions** as
   the source.
3. In **Settings → Secrets and variables → Actions → Variables**, create the
   repository variable `DOCS_PAGES_ENABLED` with the value `true`.
4. Run **Actions → Documentation site → Run workflow** on `main`.

After activation, successful documentation workflow runs on `main` publish the
static export through the `github-pages` environment. Pull requests and private
repositories never publish through this workflow. The workflow derives the
default URL and prefix from the repository owner/name; update its export-step
environment values before using a custom domain. Check navigation, search,
images, and a direct link to a nested documentation page on the first live
deployment, then update the repository's primary documentation link. Setting
`DOCS_PAGES_ENABLED` to `false` stops future deployments; it does not remove an
already published site.

## Project structure

- `app/` contains the public site, documentation routes, navigation, and styles.
- `content/upstream-rst/` preserves the authorized upstream reStructuredText
  sources while they are rewritten into the new information architecture.
- `public/og.png` is the branded social-preview image.
- `.openai/hosting.json` stores the Sites project identifier used for
  deployment.

## Content workflow

Changes to parsers, input formats, source or physical models, output formats,
or benchmark workflows should update the affected site content in the same
pull request. The documentation workflow builds the site and audits the
documented `parameters.in` keys against the current parser sources.

Documentation navigation lives in `app/docs/navigation.ts`; page content lives
in `app/docs/content.ts`. Add a navigation entry and a matching content record
to publish another guide through the shared documentation layout.

The repository's [Technical README](../doc/Technical%20README.md) is the concise
operational guide and variable index. Put lengthy engineering explanations,
examples, numerical derivations, and implementation/validation notes in
[Extended README](../doc/Extended%20README.md), keeping its current input
descriptions synchronized with the site and parser. Omit retired or unused
controls from both references. Label historical benchmark results with their
implementation scope.
