import type { Metadata } from "next";
import { headers } from "next/headers";
import "./globals.css";

export async function generateMetadata(): Promise<Metadata> {
  const requestHeaders = await headers();
  const host = requestHeaders.get("x-forwarded-host") ?? requestHeaders.get("host") ?? "localhost:3000";
  const protocol = requestHeaders.get("x-forwarded-proto") ?? (host.startsWith("localhost") ? "http" : "https");
  const base = new URL(`${protocol}://${host}`);
  const socialImage = new URL("/og.png", base).toString();

  return {
    title: {
      default: "Hercules Documentation",
      template: "%s · Hercules Documentation",
    },
    description:
      "Technical documentation for Hercules linear, nonlinear/plastic, and coupled-physics seismic simulations.",
    icons: {
      icon: "/hercules-mark.png",
      apple: "/hercules-mark.png",
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
