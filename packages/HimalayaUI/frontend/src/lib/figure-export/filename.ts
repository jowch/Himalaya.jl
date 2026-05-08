// filename.ts — slug-safe filename helpers for figure export.

/**
 * Lowercase, collapse non-alphanumeric runs into single dashes, trim leading
 * and trailing dashes. Empty / fully-stripped input returns the sentinel
 * "figure" so the resulting filename always has a visible stem segment.
 *
 * Apply this PER SEGMENT before joining with static dashes — otherwise the
 * concatenation boundary is ambiguous. See spec §Filenames.
 */
export function slugifyForFilename(s: string): string {
  if (!s) return "figure";
  const slug = s
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, "-")
    .replace(/^-+|-+$/g, "");
  return slug.length === 0 ? "figure" : slug;
}

/**
 * Append `-{YYYY-MM-DD}.{ext}` to a pre-resolved stem. Date is local time
 * (matches user expectation when exporting near midnight); locale is pinned
 * to en-CA to guarantee `YYYY-MM-DD` regardless of system locale (en-US
 * defaults produce `05/08/2026` — `/` is invalid on every OS).
 */
export function buildFilename(stem: string, ext: "png" | "svg"): string {
  const date = new Intl.DateTimeFormat("en-CA", {
    year: "numeric",
    month: "2-digit",
    day: "2-digit",
  }).format(new Date());
  return `${stem}-${date}.${ext}`;
}
