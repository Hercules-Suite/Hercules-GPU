import { docPages, type DocTable } from "./content";

type ParameterTarget = { href: string };

function tableKeys(table: DocTable, row: string[]): string[] {
  if (table.headers[0] === "Parameter") return [row[0]];
  const keys = ["Key", "Variable", "Setting", "File or key"].includes(table.headers[0])
    ? row[0]
    : table.headers[1] === "Additional keys or file" ? row[1] : "";
  // Expand compact vector-key labels such as base_force_direction_x/y/z.
  return keys
    .replace(/\b([A-Za-z0-9]+(?:_[A-Za-z0-9]+)*_)([xyz])\/([xyz])(?:\/([xyz]))?\b/g,
      (_, prefix: string, first: string, second: string, third?: string) =>
        [first, second, third].filter(Boolean).map((axis) => `${prefix}${axis}`).join(" "))
    .split(/[\s,/]+/)
    .filter((key) => /^[A-Za-z0-9]+(?:_[A-Za-z0-9]+)+$/.test(key));
}

const targets = new Map<string, ParameterTarget>();
const rowAnchors = new Map<string, string>();
const rowKey = (slug: string, sectionId: string, rowIndex: number) => `${slug}/${sectionId}/${rowIndex}`;

// The parameters.in catalog is authoritative even when a guide repeats a key.
// Source-file and launcher/runtime controls keep their own canonical tables.
const pages = Object.entries(docPages).sort(([, left], [, right]) =>
  Number(right.eyebrow === "Parameter reference") - Number(left.eyebrow === "Parameter reference"));
for (const [slug, page] of pages) {
  for (const section of page.sections) {
    const table = section.table;
    if (!table) continue;
    table.rows.forEach((row, rowIndex) => {
      const keys = [...new Set(tableKeys(table, row))].filter((key) => !targets.has(key));
      if (!keys.length) return;
      // Grouped keys share one row; repeated definitions do not duplicate IDs.
      const anchor = `parameter-${keys[0]}`;
      rowAnchors.set(rowKey(slug, section.id, rowIndex), anchor);
      for (const key of keys) targets.set(key, { href: `/docs/${slug}#${anchor}` });
    });
  }
}

export const parameterTargets: ReadonlyMap<string, ParameterTarget> = targets;

export function parameterRowAnchor(slug: string, sectionId: string, rowIndex: number) {
  return rowAnchors.get(rowKey(slug, sectionId, rowIndex));
}
