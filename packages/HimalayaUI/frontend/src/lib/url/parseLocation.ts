// Spec §4.1 — pure URL parser. No side effects.

// Compare URLs are owned by ComparePage / ComparePageEdit (which use
// useParams + useNavigate). useStateFromUrl only needs to know "this is
// a compare path so set activePage='compare' and don't 404 it." Numeric
// ids are not tracked in Zustand.
export type ParsedUrl =
  | { kind: "root" }
  | { kind: "index";   experiment: string | undefined; sample: string | undefined }
  | { kind: "inspect"; experiment: string | undefined; sample: string | undefined; exposure: string | undefined }
  | { kind: "compare"; view: "list" }
  | { kind: "stale"; raw: string };

const safeDecode = (s: string): string => {
  try { return decodeURIComponent(s); } catch { return s; }
};

export function parseLocation(pathname: string, search: string): ParsedUrl {
  // Normalize trailing slash and leading slash.
  const segs = pathname.replace(/\/+$/, "").split("/").filter(Boolean);
  if (segs.length === 0) return { kind: "root" };

  const [page, a, b, _c] = segs;
  const decode = (s: string | undefined): string | undefined =>
    s === undefined ? undefined : safeDecode(s);

  if (page === "index" && segs.length <= 3) {
    return { kind: "index", experiment: decode(a), sample: decode(b) };
  }

  if (page === "inspect" && segs.length <= 3) {
    const experiment = decode(a);
    const sample = decode(b);
    let exposure: string | undefined = undefined;
    if (experiment !== undefined && sample !== undefined && search) {
      const params = new URLSearchParams(search.startsWith("?") ? search.slice(1) : search);
      const e = params.get("exposure");
      if (e !== null && e !== "") exposure = e;
    }
    return { kind: "inspect", experiment, sample, exposure };
  }

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
