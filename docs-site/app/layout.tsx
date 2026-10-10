import type { Metadata } from "next";
import { headers } from "next/headers";
import { assetPath } from "./site-path";
import "katex/dist/katex.min.css";
import "./globals.css";

export async function generateMetadata(): Promise<Metadata> {
  let base: URL;
  if (process.env.DOCS_STATIC_EXPORT === "1") {
    base = new URL(process.env.DOCS_SITE_URL ?? "https://zwyll.github.io/Hercules_GPU/");
  } else {
    const requestHeaders = await headers();
    const host = requestHeaders.get("x-forwarded-host") ?? requestHeaders.get("host") ?? "localhost:3000";
    const protocol = requestHeaders.get("x-forwarded-proto") ?? (host.startsWith("localhost") ? "http" : "https");
    base = new URL(`${protocol}://${host}`);
  }
  const socialImage = new URL(assetPath("/og.png"), base).toString();

  return {
    title: {
      default: "Hercules Documentation",
      template: "%s · Hercules Documentation",
    },
    description:
      "Technical documentation for Hercules linear, nonlinear/plastic, and coupled-physics seismic simulations.",
    icons: {
      icon: assetPath("/hercules-octree-mesh-mark.png?v=graded-centered-256"),
      apple: assetPath("/hercules-octree-mesh-mark.png?v=graded-centered-256"),
    },
    metadataBase: base,
    openGraph: {
      title: "Hercules Documentation",
      description: "Build, configure, run, and validate Hercules linear, nonlinear/plastic, and coupled-physics simulations.",
      type: "website",
      images: [{ url: socialImage, width: 1732, height: 909, alt: "Hercules Documentation" }],
    },
    twitter: {
      card: "summary_large_image",
      title: "Hercules Documentation",
      description: "Build, configure, run, and validate Hercules linear, nonlinear/plastic, and coupled-physics simulations.",
      images: [socialImage],
    },
  };
}

export default function RootLayout({ children }: Readonly<{ children: React.ReactNode }>) {
  return (
    <html lang="en" suppressHydrationWarning>
      <body>{children}</body>
    </html>
  );
}
