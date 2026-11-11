// Next.js prefixes Link destinations itself; plain links and public assets need it here.
export function assetPath(path: string): string {
  return `${process.env.NEXT_PUBLIC_DOCS_BASE_PATH ?? ""}${path}`;
}
