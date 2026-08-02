import type { Metadata } from "next";
import Link from "next/link";
import { notFound } from "next/navigation";
import { docPages } from "../content";
import { allNavItems, navItemForSlug } from "../navigation";

export function generateStaticParams() {
  return allNavItems.map((item) => ({ slug: item.slug }));
}

export async function generateMetadata({
  params,
}: {
  params: Promise<{ slug: string }>;
}): Promise<Metadata> {
  const { slug } = await params;
  const page = docPages[slug];
  return page ? { title: page.title, description: page.summary } : {};
}

export default async function DocPage({ params }: { params: Promise<{ slug: string }> }) {
  const { slug } = await params;
  const page = docPages[slug];
  if (!page) notFound();

  const index = allNavItems.findIndex((item) => item.slug === slug);
  const previous = index > 0 ? allNavItems[index - 1] : undefined;
  const next = index < allNavItems.length - 1 ? allNavItems[index + 1] : undefined;

  return (
    <main className="docsMain">
      <article className="docsArticle">
        <div className="docBreadcrumbs">
          <Link href="/docs">Docs</Link>
          <span>/</span>
          <span>{navItemForSlug(slug)?.title}</span>
        </div>
        <header className="docHeader">
          <span className="docEyebrow">{page.eyebrow}</span>
          <h1>{page.title}</h1>
          <p className="docLead">{page.summary}</p>
          <div className="docUpdated">
            Current Hercules GPU code path
          </div>
        </header>

        {page.sections.map((section) => (
          <section className="docSection" id={section.id} key={section.id}>
            <h2>{section.title}</h2>
            {section.paragraphs?.map((paragraph) => <p key={paragraph}>{paragraph}</p>)}
            {section.bullets && (
              <ul>
                {section.bullets.map((bullet) => <li key={bullet}>{bullet}</li>)}
              </ul>
            )}
            {section.table && (
              <div className="tableWrap">
                <table>
                  <thead>
                    <tr>
                      {section.table.headers.map((header) => <th key={header}>{header}</th>)}
                    </tr>
                  </thead>
                  <tbody>
                    {section.table.rows.map((row) => (
                      <tr key={row.join("-")}>
                        {row.map((cell, cellIndex) => (
                          <td key={`${cell}-${cellIndex}`}>{cell}</td>
                        ))}
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            )}
            {section.code && (
              <div className="codeBlock">
                <div className="codeHeader">
                  <small>{section.codeLabel ?? "configuration"}</small>
                </div>
                <pre><code>{section.code}</code></pre>
              </div>
            )}
            {section.note && (
              <aside className="docNote">
                <strong>Note</strong>
                <p>{section.note}</p>
              </aside>
            )}
          </section>
        ))}

        <nav className="docPagination" aria-label="Adjacent documentation pages">
          {previous ? (
            <Link href={`/docs/${previous.slug}`}>
              <small>Previous</small>
              <strong>← {previous.title}</strong>
            </Link>
          ) : <span />}
          {next && (
            <Link href={`/docs/${next.slug}`}>
              <small>Next</small>
              <strong>{next.title} →</strong>
            </Link>
          )}
        </nav>
      </article>

      <aside className="onThisPage" aria-label="On this page">
        <strong>On this page</strong>
        {page.sections.map((section) => (
          <a href={`#${section.id}`} key={section.id}>{section.title}</a>
        ))}
        <div className="tocRule" />
        <a href="https://github.com/zwyll/GPU_new" target="_blank" rel="noreferrer">
          Edit on GitHub ↗
        </a>
      </aside>
    </main>
  );
}
