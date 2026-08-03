# Hercules Documentation

Technical documentation for the Hercules adaptive finite-element seismic wave
simulation code. The site covers installation, parameter and input files,
GPU/MPI execution, physical models, outputs, validation, and troubleshooting.

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

## Project structure

- `app/` contains the public site, documentation routes, navigation, and styles.
- `content/upstream-rst/` preserves the authorized upstream reStructuredText
  sources while they are rewritten into the new information architecture.
- `public/og.png` is the branded social-preview image.
- `.openai/hosting.json` stores the Sites project identifier used for
  deployment.

## Content workflow

Documentation navigation lives in `app/docs/navigation.ts`; page content lives
in `app/docs/content.ts`. Add a navigation entry and a matching content record
to publish another guide through the shared documentation layout.
