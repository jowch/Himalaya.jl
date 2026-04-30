/**
 * units.ts — pretty-print unit strings for axis labels.
 *
 * Storage convention is ASCII: users type `A-1`, `nm-1`, or `1/A` in
 * `experiment.toml` to avoid Unicode-input pain. The display layer renders
 * those as `Å⁻¹`, `nm⁻¹`, and `Å⁻¹` respectively.
 *
 * Idempotent — already-prettified strings are returned unchanged.
 */

const SUPERSCRIPT: Record<string, string> = {
  "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴",
  "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹",
};

/**
 * Convert ASCII unit strings to their Unicode-pretty form.
 *
 * Rules (applied in order):
 *  1. `1/X` is canonicalised to `X-1` (so `1/A` and `A-1` end up the same).
 *  2. ASCII negative-power suffix `-N` (digits) becomes `⁻N` superscript.
 *  3. A standalone capital `A` immediately before a `⁻` superscript becomes `Å`
 *     (Angstrom symbol). Other letters (`nm`, `cm`, `mm`, …) are untouched.
 *
 * Examples:
 *  - `"A-1"` → `"Å⁻¹"`
 *  - `"nm-1"` → `"nm⁻¹"`
 *  - `"1/A"` → `"Å⁻¹"`
 *  - `"Å⁻¹"` → `"Å⁻¹"` (idempotent)
 */
export function prettifyUnits(s: string): string {
  if (!s) return s;
  // 1. "1/X" → "X-1"
  let out = s.replace(/^1\/(.+)$/, "$1-1");
  // 2. "-N" → "⁻N" (Unicode superscript)
  out = out.replace(/-(\d+)/g, (_match, digits: string) =>
    "⁻" + Array.from(digits).map((d) => SUPERSCRIPT[d] ?? d).join(""));
  // 3. Capital A before a superscript power → Å (only when it's a standalone
  //    Angstrom symbol, not part of a longer letter run like "Asomething").
  out = out.replace(/(^|[^A-Za-z])A(?=⁻)/g, "$1Å");
  return out;
}
