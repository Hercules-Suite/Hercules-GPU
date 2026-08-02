import Link from "next/link";
import { navGroups } from "./navigation";

export function DocsSidebar({ currentSlug }: { currentSlug?: string }) {
  return (
    <aside className="docsSidebar" aria-label="Documentation navigation">
      <div className="sidebarScroll">
        <Link className="sidebarHome" href="/docs">
          Documentation home
        </Link>
        {navGroups.map((group) => (
          <section className="sidebarGroup" key={group.label}>
            <h2>{group.label}</h2>
            {group.items.map((item) => (
              <Link
                className={item.slug === currentSlug ? "active" : undefined}
                href={`/docs/${item.slug}`}
                key={item.slug}
              >
                {item.title}
              </Link>
            ))}
          </section>
        ))}
      </div>
    </aside>
  );
}
