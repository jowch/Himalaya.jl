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

  return { kind: "stale", raw: pathname + (search ?? "") };
}
