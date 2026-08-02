import Link from "next/link";
import { navGroups } from "./navigation";

export default function DocsOverview() {
  return (
    <main className="docsMain docsOverview">
      <div className="docsArticle wide">
        <span className="docEyebrow">Documentation</span>
        <h1>Build, run, and understand Hercules</h1>
        <p className="docLead">
          Follow a task-oriented guide or jump directly to the configuration and file
          references used by production simulations.
        </p>
        <div className="docsOverviewGrid">
          {navGroups.map((group) => (
            <section key={group.label}>
              <h2>{group.label}</h2>
              {group.items.map((item) => (
                <Link href={`/docs/${item.slug}`} key={item.slug}>
                  <strong>{item.title}</strong>
                  <span>{item.description}</span>
                </Link>
              ))}
            </section>
          ))}
        </div>
      </div>
    </main>
  );
}
