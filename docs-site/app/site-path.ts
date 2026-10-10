// Next.js prefixes Link destinations itself; public image URLs need the prefix.
export function assetPath(path: string): string {
  return `${process.env.NEXT_PUBLIC_DOCS_BASE_PATH ?? ""}${path}`;
}
