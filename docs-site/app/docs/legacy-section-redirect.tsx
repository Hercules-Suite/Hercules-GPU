"use client";

import { useEffect } from "react";
import { useRouter } from "next/navigation";

// URL fragments never reach the server. Resolve old section bookmarks in the
// browser; the page also renders ordinary links as a no-JavaScript fallback.
export function LegacySectionRedirect({ destinations }: { destinations: Record<string, string> }) {
  const router = useRouter();
  useEffect(() => {
    const redirect = () => {
      let id: string;
      try { id = decodeURIComponent(window.location.hash.slice(1)); }
      catch { return; }
      const destination = destinations[id];
      if (destination) router.replace(destination, { scroll: true });
    };
    redirect();
    window.addEventListener("hashchange", redirect);
    return () => window.removeEventListener("hashchange", redirect);
  }, [destinations, router]);
  return null;
}
