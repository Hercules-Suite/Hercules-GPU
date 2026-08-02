import Link from "next/link";
import { navGroups } from "./navigation";

export function DocsSidebar({ currentSlug }: { currentSlug?: string }) {
  return (
    <aside className="docsSidebar" aria-label="Documentation navigation">
      <div className="sidebarScroll">
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
        <div className="sidebarStatus">
          <span className="statusDot" />
          <div>
            <strong>GPU documentation</strong>
            <small>Continuously evolving</small>
          </div>
        </div>
      </div>
    </aside>
  );
}
