import katex from "katex";
import type { ReactNode } from "react";

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

    if (!inline && !display) return part;

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
