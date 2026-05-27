// Spec §4.1 — pure URL parser. No side effects.

// Compare URLs are owned by ComparePage / ComparePageEdit (which use
// useParams + useNavigate). useStateFromUrl only needs to know "this is
// a compare path so set activePage='compare' and don't 404 it." Numeric
// ids are not tracked in Zustand.
export type ParsedUrl =
  | { kind: "root" }
  | { kind: "compare"; view: "list" }
  | { kind: "stale"; raw: string };

export function parseLocation(pathname: string, search: string): ParsedUrl {
  // Normalize trailing slash and leading slash.
  const segs = pathname.replace(/\/+$/, "").split("/").filter(Boolean);
  if (segs.length === 0) return { kind: "root" };

  const [page, a] = segs;

  // I4.4 (#181): Index is retired. `/index*` is redirected by the router
  // (sampleless → `/samples`; a slug-bearing `/index/:exp/:sample` →
  // `/sample/:id` via `IndexSlugRedirect`) before this parser runs, so there
  // is no `index` parse arm — any stray /index path falls through to `stale`.
  //
  // I1.7 (#163): Inspect is retired. `/inspect*` is redirected to `/samples`
  // by the router before this parser runs, so there is no `inspect` parse arm
  // — any stray /inspect path falls through to `stale`.

  // Compare URLs in the actual app live under two roots:
  //   /experiments/:eid/compare(/new|/:id|/:id/edit)
  //   /compare/all(/new|/:id|/:id/edit)
  // The hooks don't track Compare numeric ids in Zustand — ComparePage
  // handles its own URL via useNavigate/useParams. Recognizing the shape
  // here just lets useStateFromUrl set activePage="compare" without
  // falsely flagging the path as "stale".
  const isExperimentCompare = page === "experiments" && segs[2] === "compare";
  const isCompareAll = page === "compare" && a === "all";
  const isLegacyCompareList = page === "compare" && segs.length === 1;

  if (isExperimentCompare || isCompareAll || isLegacyCompareList) {
    return { kind: "compare", view: "list" };
  }

  return { kind: "stale", raw: pathname + (search ?? "") };
}
