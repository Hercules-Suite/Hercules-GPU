import Image from "next/image";
import Link from "next/link";
import { SearchBox } from "./search-box";
import { ThemeToggle } from "./theme-toggle";
import { assetPath } from "../site-path";

const topLinks = [
  { href: "/docs", label: "Documentation" },
  { href: "/docs/parameter-file", label: "Parameter reference" },
];

export function SiteHeader() {
  return (
    <header className="siteHeader">
      <div className="headerInner">
        <Link className="brand" href="/" aria-label="Hercules documentation home">
          <span className="brandMark" aria-hidden="true">
            <Image src={assetPath("/hercules-octree-mesh-mark.png?v=graded-centered-256")} alt="" width={34} height={34} priority unoptimized />
          </span>
          <span className="brandCopy">
            <strong>Hercules Documentation</strong>
          </span>
          <span className="versionTag">current</span>
        </Link>

        <nav className="topNav" aria-label="Primary navigation">
          {topLinks.map((link) => (
            <Link href={link.href} key={link.href}>
              {link.label}
            </Link>
          ))}
          <a href="https://github.com/Hercules-Suite/Hercules-GPU" target="_blank" rel="noreferrer">
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
              <a href="https://github.com/Hercules-Suite/Hercules-GPU">GitHub ↗</a>
            </div>
          </details>
        </div>
      </div>
    </header>
  );
}
