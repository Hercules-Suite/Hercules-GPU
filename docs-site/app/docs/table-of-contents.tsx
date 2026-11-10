"use client";

import { useEffect, useState } from "react";

export type TocSection = { id: string; title: string; group?: string };

export function TableOfContents({ sections, mobile = false }: { sections: TocSection[]; mobile?: boolean }) {
  const [active, setActive] = useState(sections[0]?.id ?? "");
  useEffect(() => {
    let frame = 0;
    const update = () => {
      frame = 0;
      let current = sections[0]?.id ?? "";
      for (const section of sections) {
        const element = document.getElementById(section.id);
        if (element && element.getBoundingClientRect().top <= 150) current = section.id;
      }
      if (window.scrollY > 0 && window.scrollY + window.innerHeight >= document.documentElement.scrollHeight - 2) {
        current = sections.at(-1)?.id ?? current;
      }
      setActive(current);
    };
    const schedule = () => { if (!frame) frame = requestAnimationFrame(update); };
    schedule();
    window.addEventListener("scroll", schedule, { passive: true });
    window.addEventListener("resize", schedule);
    window.addEventListener("hashchange", schedule);
    return () => {
      cancelAnimationFrame(frame);
      window.removeEventListener("scroll", schedule);
      window.removeEventListener("resize", schedule);
      window.removeEventListener("hashchange", schedule);
    };
  }, [sections]);

  const groups: { title?: string; sections: TocSection[] }[] = [];
  for (const section of sections) {
    const previous = groups.at(-1);
    if (previous && previous.title === section.group) previous.sections.push(section);
    else groups.push({ title: section.group, sections: [section] });
  }
  const links = <nav aria-label={mobile ? "Page contents" : "On this page"}>
    {groups.map((group, index) => <div className="tocGroup" key={`${group.title}-${index}`}>
      {group.title && <p className="tocGroupTitle">{group.title}</p>}
      <ul>{group.sections.map((section) => <li key={section.id}>
        <a href={`#${section.id}`} aria-current={active === section.id ? "location" : undefined}>{section.title}</a>
      </li>)}</ul>
    </div>)}
  </nav>;
  if (mobile) return <details className="mobilePageContents">
    <summary>On this page <span>{sections.length} sections</span></summary>
    {links}
  </details>;
  return <aside className="onThisPage">
    <div className="tocHeading"><strong>On this page</strong><span>{sections.length} sections</span></div>
    <div className="tocScroll">{links}</div>
    <a className="tocFooter" href="https://github.com/Hercules-Suite/Hercules-GPU" target="_blank" rel="noreferrer">View on GitHub ↗</a>
  </aside>;
}
