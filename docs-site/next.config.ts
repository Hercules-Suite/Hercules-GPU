import type { NextConfig } from "next";

const staticExport = process.env.DOCS_STATIC_EXPORT === "1";
const basePath = staticExport
  ? (process.env.NEXT_PUBLIC_DOCS_BASE_PATH ?? "/Hercules-GPU")
  : "";

const nextConfig: NextConfig = {
  ...(staticExport ? { output: "export", trailingSlash: true } : {}),
  basePath,
  env: { NEXT_PUBLIC_DOCS_BASE_PATH: basePath },
  images: { unoptimized: true },
};

export default nextConfig;
