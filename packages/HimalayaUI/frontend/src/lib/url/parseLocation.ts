// Spec §4.1 — pure URL parser. No side effects.

// Every legacy AppShell surface is retired — Inspect (#163), Index (#181),
// Compare (#177) — and their URLs are redirected at the router before this
// parser runs. So the only kinds this hook needs are `root` (bare `/`) and
// `stale` (anything else that reaches the AppShell catch-all). There is no
// longer a slug-bearing `index`/`inspect`/`compare` arm.
export type ParsedUrl =
  | { kind: "root" }
  | { kind: "stale"; raw: string };

export function parseLocation(pathname: string, search: string): ParsedUrl {
  // Normalize trailing slash and leading slash.
  const segs = pathname.replace(/\/+$/, "").split("/").filter(Boolean);
  if (segs.length === 0) return { kind: "root" };

  // I1.7 (#163) / I4.4 (#181) / I3.6 (#177): Inspect, Index, and Compare are
  // all retired. `/inspect*`, `/index*`, and `/compare*` are redirected by the
  // router (Inspect → /samples; Index → /samples or /sample/:id; Compare →
  // /series) before this parser runs, so there are no parse arms for them —
  // any stray legacy path falls through to `stale`.
  return { kind: "stale", raw: pathname + (search ?? "") };
}
