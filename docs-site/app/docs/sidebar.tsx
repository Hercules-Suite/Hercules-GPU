import Link from "next/link";
import { navGroups, navItemForSlug } from "./navigation";

function NavigationGroups({ currentSlug }: { currentSlug?: string }) {
  return <>
    <Link className="sidebarHome" href="/docs">Documentation home</Link>
    {navGroups.map((group) => <section className="sidebarGroup" key={group.label}>
      <h2>{group.label}</h2>
      {group.items.map((item) => item.children?.length ? (
        <details className="sidebarTopic" key={`${item.slug}-${currentSlug}`} open={item.slug === currentSlug || item.children.some((child) => child.slug === currentSlug)}>
          <summary>{item.title}</summary>
          <div className="sidebarChildren">
            {[{ ...item, title: `${item.title} overview` }, ...item.children].map((child) => <Link
              className={child.slug === currentSlug ? "active" : undefined}
              aria-current={child.slug === currentSlug ? "page" : undefined}
              href={`/docs/${child.slug}`} key={child.slug}
            >{child.title}</Link>)}
          </div>
        </details>
      ) : <Link
        className={item.slug === currentSlug ? "active" : undefined}
        aria-current={item.slug === currentSlug ? "page" : undefined}
        href={`/docs/${item.slug}`} key={item.slug}
      >{item.title}</Link>)}
    </section>)}
  </>;
}

export function DocsSidebar({ currentSlug }: { currentSlug?: string }) {
  return <aside className="docsSidebar" aria-label="Documentation navigation">
    <div className="sidebarScroll"><NavigationGroups currentSlug={currentSlug} /></div>
    <details className="mobileDocsContents" key={currentSlug}>
      <summary>Browse documentation <span>{navItemForSlug(currentSlug ?? "")?.title ?? "Contents"}</span></summary>
      <div className="mobileDocsTree"><NavigationGroups currentSlug={currentSlug} /></div>
    </details>
  </aside>;
}
