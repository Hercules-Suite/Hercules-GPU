import type { ReactNode } from "react";
import { docPages } from "../docs/content";
import { documentedParameterNames } from "../docs/parameter-reference";

// Reuse the catalogs rather than guessing that every underscored word is an
// input parameter. Source-file keys and launcher controls have their own tables.
const parameterNames = new Set([
  ...documentedParameterNames,
  ...Object.values(docPages).flatMap((page) => page.sections.flatMap(({ table }) =>
    table && ["Key", "Variable", "Setting", "File or key"].includes(table.headers[0])
      ? table.rows.flatMap((row) => row[0].split(/[\s,/]+/).filter((key) => /^[A-Za-z0-9]+(?:_[A-Za-z0-9]+)+$/.test(key)))
      : table?.headers[1] === "Additional keys or file"
        ? table.rows.flatMap((row) => row[1].split(/\s+/).filter((key) => /^[A-Za-z0-9]+(?:_[A-Za-z0-9]+)+$/.test(key)))
        : [],
  )),
]);

export function ParameterText({ text }: { text: string }) {
  const parts: ReactNode[] = [];
  let offset = 0;
  for (const match of text.matchAll(/[A-Za-z0-9_]+(?:\.[A-Za-z0-9_]+)*/g)) {
    if (!parameterNames.has(match[0])) continue;
    if (match.index > offset) parts.push(text.slice(offset, match.index));
    parts.push(<code className="parameterName" key={match.index}>{match[0]}</code>);
    offset = match.index + match[0].length;
  }
  if (offset < text.length) parts.push(text.slice(offset));
  return <>{parts}</>;
}
