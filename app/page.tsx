import { SiteHeader } from "./components/site-header";
import { DocsShell } from "./docs/docs-shell";
import DocsOverview from "./docs/page";

export default function Home() {
  return (
    <>
      <SiteHeader />
      <DocsShell>
        <DocsOverview />
      </DocsShell>
    </>
  );
}
