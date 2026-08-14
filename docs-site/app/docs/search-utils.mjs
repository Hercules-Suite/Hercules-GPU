/**
 * Normalize documentation text so parameter names, paths, and ordinary words
 * can be matched with the same space-separated query.
 *
 * @param {string} value
 */
export function normalizeSearchText(value) {
  return value
    .toLowerCase()
    .replace(/[_./-]+/g, " ")
    .replace(/[^a-z0-9+]+/g, " ")
    .replace(/\s+/g, " ")
    .trim();
}

/**
 * Match every query term, independent of its order in the document.
 *
 * @param {string} normalizedText
 * @param {string} query
 */
export function matchesSearchQuery(normalizedText, query) {
  const terms = normalizeSearchText(query).split(" ").filter(Boolean);
  return terms.length > 0 && terms.every((term) => normalizedText.includes(term));
}

/**
 * Prefer title matches, then summary/context matches, then body-only matches.
 *
 * @param {{ normalizedTitle: string; normalizedDescription: string; searchText: string; kind: "page" | "section" }} item
 * @param {string} query
 */
export function searchResultScore(item, query) {
  const terms = normalizeSearchText(query).split(" ").filter(Boolean);
  let score = item.kind === "page" ? 0 : 1;

  for (const term of terms) {
    if (item.normalizedTitle.includes(term)) score += 0;
    else if (item.normalizedDescription.includes(term)) score += 2;
    else if (item.searchText.includes(term)) score += 4;
  }

  return score;
}
