"use client";

import Link from "next/link";
import { useEffect, useRef, useState } from "react";
import { allNavItems } from "../docs/navigation";

export function SearchBox() {
  const [query, setQuery] = useState("");
  const [open, setOpen] = useState(false);
  const wrapperRef = useRef<HTMLDivElement>(null);
  const inputRef = useRef<HTMLInputElement>(null);

  const results = query.trim()
    ? allNavItems
        .filter((item) =>
          `${item.title} ${item.description}`
            .toLowerCase()
            .includes(query.trim().toLowerCase()),
        )
        .slice(0, 6)
    : allNavItems.slice(0, 5);

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
        aria-expanded={open}
      />
      <kbd>⌘K</kbd>
      {open && (
        <div className="searchResults" role="listbox" aria-label="Documentation pages">
          <p className="searchResultsLabel">
            {query ? "Matching pages" : "Popular pages"}
          </p>
          {results.length ? (
            results.map((item) => (
              <Link
                href={`/docs/${item.slug}`}
                key={item.slug}
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
            <p className="searchEmpty">No matching page yet.</p>
          )}
        </div>
      )}
    </div>
  );
}
