"use client";

import { usePathname } from "next/navigation";
import { DocsSidebar } from "./sidebar";

export function DocsShell({ children }: { children: React.ReactNode }) {
  const pathname = usePathname();
  const slug = pathname.split("/").filter(Boolean).at(-1);

  return (
    <div className="docsShell">
      <DocsSidebar currentSlug={slug} />
      {children}
    </div>
  );
}
