# Hercules Documentation

Modern documentation for the Hercules adaptive finite-element earthquake
simulation framework. The site covers GPU-resident execution, coupled physics,
input and output reference, validation, and troubleshooting.

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
- `.openai/hosting.json` stores metadata for the review deployment.

## Content workflow

Documentation navigation lives in `app/docs/navigation.ts`; page content lives
in `app/docs/content.ts`. Add a navigation entry and a matching content record
to publish another guide through the shared documentation layout.
