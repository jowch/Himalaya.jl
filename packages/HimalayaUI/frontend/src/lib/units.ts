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

/**
 * Real-space lattice unit derived from the q-units stored on the experiment.
 * `1/A`, `A-1`, `Å⁻¹` → `Å`. `nm-1`, `nm⁻¹` → `nm`. Defaults to `Å` when
 * unset (matches the backend's `A-1` default in `experiment.toml`).
 */
export function latticeUnitFromQUnits(qUnits: string | null | undefined): string {
  if (!qUnits) return "Å";
  const lower = qUnits.toLowerCase();
  if (lower.includes("nm")) return "nm";
  return "Å";
}

/** Inverse-square form of the lattice unit (e.g. "Å⁻²", "nm⁻²"). */
export function inverseSquareUnits(qUnits: string | null | undefined): string {
  return `${latticeUnitFromQUnits(qUnits)}⁻²`;
}

/**
 * Format a κ (Gaussian curvature) value for display.
 *
 * κ ranges over many orders of magnitude depending on the lattice unit
 * (≈10⁻² nm⁻² ≈ 10⁻⁴ Å⁻²); fall back to scientific notation when fixed-point
 * would round all the leading digits to zero.
 */
export function formatKappa(k: number): string {
  const abs = Math.abs(k);
  if (abs === 0) return "0";
  if (abs < 0.01) return k.toExponential(2);
  return k.toFixed(3);
}

/**
 * Figure-grade κ: like {@link formatKappa} but renders the scientific-notation
 * branch with a real `×10ⁿ` (superscript) mantissa instead of the `e-4` ASCII
 * form, for a clean journal/presentation label (e.g. `1.70×10⁻⁴`).
 */
export function formatKappaPretty(k: number): string {
  const s = formatKappa(k);
  const m = /^(-?\d+(?:\.\d+)?)e([+-]?\d+)$/.exec(s);
  if (!m) return s;
  const mantissa = m[1]!;
  const expDigits = m[2]!.replace("+", "");
  const sup = Array.from(expDigits)
    .map((ch) => (ch === "-" ? "⁻" : (SUPERSCRIPT[ch] ?? ch)))
    .join("");
  return `${mantissa}×10${sup}`;
}
