import { SiteHeader } from "../components/site-header";
import { DocsShell } from "./docs-shell";

export default function DocumentationLayout({ children }: { children: React.ReactNode }) {
  return (
    <>
      <SiteHeader />
      <DocsShell>{children}</DocsShell>
    </>
  );
}
