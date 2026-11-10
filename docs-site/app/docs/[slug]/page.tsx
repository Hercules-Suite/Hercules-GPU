import type { Metadata } from "next";
import Link from "next/link";
import { notFound } from "next/navigation";
import { MathFormula, MathText } from "../../components/math-text";
import { docPages } from "../content";
import { TableOfContents } from "../table-of-contents";
import { LegacySectionRedirect } from "../legacy-section-redirect";
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

  const tocSections = page.sections.map((section) => ({ id: section.id, title: section.tocTitle ?? section.title, group: section.tocGroup }));
  if (page.legacySections?.length) tocSections.push({ id: "full-topic-index", title: "Full topic index", group: "All topics" });
  const parent = page.parent ? navItemForSlug(page.parent) : undefined;
  const index = allNavItems.findIndex((item) => item.slug === slug);
  const previous = index > 0 ? allNavItems[index - 1] : undefined;
  const next = index < allNavItems.length - 1 ? allNavItems[index + 1] : undefined;

  return (
    <main className={`docsMain ${page.eyebrow === "Parameter reference" ? "parameterReferenceMain" : ""}`}>
      <article className="docsArticle">
        {page.legacySections?.length ? <LegacySectionRedirect destinations={Object.fromEntries(page.legacySections.map((section) => [section.id, section.href]))} /> : null}
        <div className="docBreadcrumbs">
          <Link href="/docs">Docs</Link>
          <span>/</span>
          {parent && <><Link href={`/docs/${parent.slug}`}>{parent.title}</Link><span>/</span></>}
          <span>{navItemForSlug(slug)?.title}</span>
        </div>
        <header className="docHeader">
          <span className="docEyebrow">{page.eyebrow}</span>
          <h1><MathText text={page.title} /></h1>
          <p className="docLead"><MathText text={page.summary} /></p>
          <div className="docUpdated">
            Current Hercules GPU code path
          </div>
        </header>

        <TableOfContents mobile sections={tocSections} />
        {page.sections.map((section) => (
          <section className="docSection" id={section.id} key={section.id}>
            <h2><MathText text={section.title} /></h2>
            {section.links && <div className="guideLinks">
              {section.links.map((link) => <Link href={link.href} key={link.href}>
                <strong><MathText text={link.title} /><span aria-hidden="true"> →</span></strong>
                <span><MathText text={link.description} /></span>
              </Link>)}
            </div>}
            {section.paragraphs?.map((paragraph) => <p key={paragraph}><MathText text={paragraph} /></p>)}
            {section.equations?.map((equation) => (
              <MathFormula display expression={equation} key={equation} />
            ))}
            {section.afterEquationsParagraphs?.map((paragraph) => (
              <p key={paragraph}><MathText text={paragraph} /></p>
            ))}
            {section.bullets && (
              <ul>
                {section.bullets.map((bullet) => <li key={bullet}><MathText text={bullet} /></li>)}
              </ul>
            )}
            {section.table && (
              <div className={`tableWrap ${section.table.headers[0] === "Parameter" ? "parameterTableWrap" : ""}`}>
                {section.table.headers[0] === "Parameter" && (
                  <div className="parameterTableHint">Scroll horizontally to see every column.</div>
                )}
                <table className={section.table.headers[0] === "Parameter" ? "parameterTable" : undefined}>
                  <thead>
                    <tr>
                      {section.table.headers.map((header) => <th key={header}><MathText text={header} /></th>)}
                    </tr>
                  </thead>
                  <tbody>
                    {section.table.rows.map((row) => (
                      <tr key={row.join("-")}>
                        {row.map((cell, cellIndex) => (
                          <td
                            data-label={section.table?.headers[0] === "Parameter" ? section.table.headers[cellIndex] : undefined}
                            data-requirement={cellIndex === 1 ? cell : undefined}
                            key={`${cell}-${cellIndex}`}
                          >
                            <MathText text={cell} />
                          </td>
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
                <p><MathText text={section.note} /></p>
              </aside>
            )}
          </section>
        ))}

        {page.legacySections?.length ? <section className="docSection" id="full-topic-index">
          <h2>Full topic index</h2>
          <details className="legacyTopicIndex">
            <summary>Browse all {page.legacySections.length} detailed topics</summary>
            <ul>{page.legacySections.map((section) => <li id={section.id} key={section.id}>
              <Link data-legacy-section={section.id} href={section.href}><MathText text={section.title} /></Link>
            </li>)}</ul>
          </details>
        </section> : null}
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

      <TableOfContents sections={tocSections} />
    </main>
  );
}
