import type { ReactNode } from "react";
import { parameterTargets } from "../docs/parameter-links";
import { documentedParameters, srfhDepthReferenceParameter } from "../docs/parameter-reference";
import { assetPath } from "../site-path";

// Reuse the catalogs rather than guessing that every underscored word is an
// input parameter. Source-file keys and launcher controls have their own tables.
const parameterNames = new Set(parameterTargets.keys());

// Underscored options are unambiguous in prose. Common words such as "surface"
// need explicit backticks so physical surfaces remain ordinary text.
const distinctiveValues = new Set(
  [...documentedParameters, srfhDepthReferenceParameter]
    .filter(({ type }) => /String|keyword/.test(type))
    .flatMap(({ values }) => values.match(/[A-Za-z0-9]+(?:_[A-Za-z0-9]+)+/g) ?? []),
);

export function ParameterText({ text, linkParameters = true }: { text: string; linkParameters?: boolean }) {
  const renderLiteral = (literal: string, key: string | number) => {
    const target = parameterTargets.get(literal);
    const code = <code className={target ? "parameterName" : "parameterValue"} key={key}>{literal}</code>;
    return target && linkParameters
      ? <a className="parameterLink" href={assetPath(target.href)} title="View parameter definition" key={key}>{code}</a>
      : code;
  };
  const parts: ReactNode[] = [];
  let offset = 0;
  for (const match of text.matchAll(/`([^`\n]+)`|[A-Za-z0-9_]+(?:\.[A-Za-z0-9_]+)*/g)) {
    if (match.index < offset) continue;
    const literal = match[1] ?? match[0];
    const isParameter = parameterNames.has(literal);
    if (!match[1] && !isParameter && !distinctiveValues.has(literal)) continue;
    if (match.index > offset) parts.push(text.slice(offset, match.index));
    parts.push(renderLiteral(literal, match.index));
    offset = match.index + match[0].length;

    // Only a documented key introduces an assignment; equations such as Z=0
    // remain prose/math. Include numeric, scientific-notation and path values.
    if (!isParameter) continue;
    const assignment = text.slice(offset).match(/^(\s*=\s*)(`[^`\n]+`|"[^"\n]*"|'[^'\n]*'|[^\s,;:()[\]{}<>="'`]+)/);
    if (!assignment) continue;
    const value = assignment[2].replace(/[.!?]+$/, "");
    if (!value) continue;
    const displayValue = value.startsWith("`") ? value.slice(1, -1) : value;
    parts.push(assignment[1]);
    parts.push(renderLiteral(displayValue, `${offset}-value`));
    offset += assignment[1].length + value.length;
  }
  if (offset < text.length) parts.push(text.slice(offset));
  return <>{parts}</>;
}
