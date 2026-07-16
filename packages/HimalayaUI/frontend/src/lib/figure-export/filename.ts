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
 * Slugify the stem, then append `-{YYYY-MM-DD}.{ext}`. SLUGIFYING IS LOAD-BEARING:
 * a raw title carries spaces, colons, and dots ("1:2.5 ratio", "JC042.dat") into
 * the filename. Spaces/dots before the extension make some browsers and OSes
 * mis-detect (or drop) the `.svg`/`.png` extension, and the unsanitized name can
 * defeat the anchor `download` attribute. Slugging guarantees a single, trailing
 * `.{ext}`. Idempotent on already-clean stems ("himalaya-trace-jc23").
 *
 * Date is local time (matches user expectation when exporting near midnight);
 * locale is pinned to en-CA to guarantee `YYYY-MM-DD` regardless of system
 * locale (en-US defaults produce `05/08/2026` — `/` is invalid on every OS).
 */
export function buildFilename(stem: string, ext: "png" | "svg"): string {
  const date = new Intl.DateTimeFormat("en-CA", {
    year: "numeric",
    month: "2-digit",
    day: "2-digit",
  }).format(new Date());
  return `${slugifyForFilename(stem)}-${date}.${ext}`;
}
