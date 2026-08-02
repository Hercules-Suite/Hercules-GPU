import Link from "next/link";
import { SearchBox } from "./search-box";
import { ThemeToggle } from "./theme-toggle";

const topLinks = [
  { href: "/docs/getting-started", label: "Docs" },
  { href: "/docs/running", label: "Guides" },
  { href: "/docs/validation", label: "Benchmarks" },
];

export function SiteHeader() {
  return (
    <header className="siteHeader">
      <div className="headerInner">
        <Link className="brand" href="/" aria-label="Hercules documentation home">
          <span className="brandMark" aria-hidden="true">
            H
          </span>
          <span className="brandCopy">
            <strong>Hercules</strong>
            <small>Documentation</small>
          </span>
          <span className="versionTag">GPU</span>
        </Link>

        <nav className="topNav" aria-label="Primary navigation">
          {topLinks.map((link) => (
            <Link href={link.href} key={link.href}>
              {link.label}
            </Link>
          ))}
          <a href="https://github.com/zwyll/GPU_new" target="_blank" rel="noreferrer">
            GitHub <span aria-hidden="true">↗</span>
          </a>
        </nav>

        <div className="headerTools">
          <SearchBox />
          <ThemeToggle />
          <details className="mobileMenu">
            <summary aria-label="Open navigation">
              <span />
              <span />
              <span />
            </summary>
            <div className="mobileMenuPanel">
              {topLinks.map((link) => (
                <Link href={link.href} key={link.href}>
                  {link.label}
                </Link>
              ))}
              <a href="https://github.com/zwyll/GPU_new">GitHub ↗</a>
            </div>
          </details>
        </div>
      </div>
    </header>
  );
}
