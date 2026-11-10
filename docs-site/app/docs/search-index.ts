import { docPages, type DocSection } from "./content";
import { allNavItems } from "./navigation";
import { normalizeSearchText } from "./search-utils.mjs";

export type DocumentationSearchItem = {
  href: string;
  title: string;
  description: string;
  normalizedTitle: string;
  normalizedDescription: string;
  searchText: string;
  kind: "page" | "section";
};

function sectionValues(section: DocSection): string[] {
  return [
    section.title,
    section.tocTitle ?? "",
    ...(section.links?.flatMap((link) => [link.title, link.description]) ?? []),
    ...(section.figures?.flatMap((figure) => [figure.title, figure.description, figure.alt, figure.caption]) ?? []),
    ...(section.paragraphs ?? []),
    ...(section.equations ?? []),
    ...(section.afterEquationsParagraphs ?? []),
    ...(section.bullets ?? []),
    section.codeLabel ?? "",
    section.code ?? "",
    section.note ?? "",
    ...(section.table?.headers ?? []),
    ...(section.table?.rows.flat() ?? []),
  ];
}

function sectionPreview(section: DocSection): string {
  const firstDetail =
    section.paragraphs?.[0] ??
    section.bullets?.[0] ??
    section.note ??
    section.table?.rows[0]?.join(" — ") ??
    section.codeLabel ??
    "Open this section";

  return firstDetail.length > 118 ? `${firstDetail.slice(0, 115)}…` : firstDetail;
}

export const documentationSearchItems: DocumentationSearchItem[] = allNavItems.flatMap(
  (navItem) => {
    const page = docPages[navItem.slug];
    if (!page) return [];

    const pageItem: DocumentationSearchItem = {
      href: `/docs/${navItem.slug}`,
      title: page.title,
      description: page.summary,
      normalizedTitle: normalizeSearchText(page.title),
      normalizedDescription: normalizeSearchText(page.summary),
      searchText: normalizeSearchText(
        [page.eyebrow, page.title, page.summary, navItem.description].join(" "),
      ),
      kind: "page",
    };

    const sectionItems = page.sections.map(
      (section): DocumentationSearchItem => {
        const description = `${page.title} · ${sectionPreview(section)}`;
        return {
          href: `/docs/${navItem.slug}#${section.id}`,
          title: section.title,
          description,
          normalizedTitle: normalizeSearchText(section.title),
          normalizedDescription: normalizeSearchText(description),
          searchText: normalizeSearchText(sectionValues(section).join(" ")),
          kind: "section",
        };
      },
    );

    return [pageItem, ...sectionItems];
  },
);
