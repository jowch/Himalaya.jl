import type { LoadSample } from "../api";

/** Fuzzy/glob match against a sample's name and exposure filenames.
 *  `*` -> any run, `?` -> one char (anchored at start); otherwise substring. */
export function matchSample(s: LoadSample, query: string): boolean {
  const q = query.trim().toLowerCase();
  if (!q) return true;
  const fields = [s.name.toLowerCase(), ...s.exposures.map((e) => e.filename.toLowerCase())];
  if (/[*?]/.test(q)) {
    const re = new RegExp(
      "^" + q.replace(/[.+^${}()|[\]\\]/g, "\\$&").replace(/\*/g, ".*").replace(/\?/g, "."),
    );
    return fields.some((f) => re.test(f));
  }
  return fields.some((f) => f.includes(q));
}
