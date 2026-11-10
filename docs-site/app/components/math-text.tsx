import katex from "katex";
import type { ReactNode } from "react";
import { ParameterText } from "./parameter-text";

type MathFormulaProps = {
  expression: string;
  display?: boolean;
};

export function MathFormula({ expression, display = false }: MathFormulaProps) {
  const markup = katex.renderToString(expression, {
    displayMode: display,
    output: "htmlAndMathml",
    strict: "error",
    throwOnError: true,
    trust: false,
  });

  return (
    <span
      className={`mathExpression ${display ? "mathDisplay" : "mathInline"}`}
      dangerouslySetInnerHTML={{ __html: markup }}
    />
  );
}

const mathSegment = /(\\\[[\s\S]*?\\\]|\\\([\s\S]*?\\\))/g;

export function MathText({ text }: { text: string }) {
  const parts = text.split(mathSegment);
  const rendered: ReactNode[] = parts.map((part, index) => {
    const inline = part.startsWith("\\(") && part.endsWith("\\)");
    const display = part.startsWith("\\[") && part.endsWith("\\]");

    if (!inline && !display) return <ParameterText text={part} key={index} />;

    // A standalone parameter/assignment is code, not a mathematical symbol.
    // Keep actual formulas in KaTeX, including their accessible MathML.
    const literal = inline && part.slice(2, -2).trim().match(/^\\mathtt\{([^{}]+)\}$/);
    if (literal) return <ParameterText text={literal[1].replaceAll("\\_", "_")} key={index} />;

    return (
      <MathFormula
        display={display}
        expression={part.slice(2, -2).trim()}
        key={`${index}-${part}`}
      />
    );
  });

  return <>{rendered}</>;
}
