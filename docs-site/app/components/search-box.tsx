"use client";

import Link from "next/link";
import { useEffect, useRef, useState } from "react";
import { allNavItems } from "../docs/navigation";
import { documentationSearchItems, type DocumentationSearchItem } from "../docs/search-index";
import { matchesSearchQuery, searchResultScore } from "../docs/search-utils.mjs";

const popularResults: DocumentationSearchItem[] = allNavItems.slice(0, 5).map((item) => ({
  href: `/docs/${item.slug}`,
  title: item.title,
  description: item.description,
  normalizedTitle: "",
  normalizedDescription: "",
  searchText: "",
  kind: "page",
}));

export function SearchBox() {
  const [query, setQuery] = useState("");
  const [open, setOpen] = useState(false);
  const wrapperRef = useRef<HTMLDivElement>(null);
  const inputRef = useRef<HTMLInputElement>(null);

  const results = query.trim()
    ? documentationSearchItems
        .filter((item) => matchesSearchQuery(item.searchText, query))
        .sort((left, right) => searchResultScore(left, query) - searchResultScore(right, query))
        .slice(0, 8)
    : popularResults;

  useEffect(() => {
    function handleKeyboard(event: KeyboardEvent) {
      if ((event.metaKey || event.ctrlKey) && event.key.toLowerCase() === "k") {
        event.preventDefault();
        setOpen(true);
        inputRef.current?.focus();
      }
      if (event.key === "Escape") {
        setOpen(false);
        inputRef.current?.blur();
      }
    }

    function handleClick(event: MouseEvent) {
      if (!wrapperRef.current?.contains(event.target as Node)) {
        setOpen(false);
      }
    }

    window.addEventListener("keydown", handleKeyboard);
    document.addEventListener("mousedown", handleClick);
    return () => {
      window.removeEventListener("keydown", handleKeyboard);
      document.removeEventListener("mousedown", handleClick);
    };
  }, []);

  return (
    <div className="search" ref={wrapperRef}>
      <span className="searchIcon" aria-hidden="true" />
      <input
        ref={inputRef}
        value={query}
        onChange={(event) => {
          setQuery(event.target.value);
          setOpen(true);
        }}
        onFocus={() => setOpen(true)}
        placeholder="Search documentation"
        aria-label="Search documentation"
        role="combobox"
        aria-autocomplete="list"
        aria-controls={open ? "documentation-search-results" : undefined}
        aria-expanded={open}
      />
      <kbd>⌘K</kbd>
      {open && (
        <div
          className="searchResults"
          id="documentation-search-results"
          role="listbox"
          aria-label="Documentation pages"
        >
          <p className="searchResultsLabel">
            {query ? "Matching documentation" : "Popular pages"}
          </p>
          {results.length ? (
            results.map((item) => (
              <Link
                href={item.href}
                key={item.href}
                onClick={() => {
                  setOpen(false);
                  setQuery("");
                }}
              >
                <span>{item.title}</span>
                <small>{item.description}</small>
              </Link>
            ))
          ) : (
            <p className="searchEmpty">No matching documentation.</p>
          )}
        </div>
      )}
    </div>
  );
}
