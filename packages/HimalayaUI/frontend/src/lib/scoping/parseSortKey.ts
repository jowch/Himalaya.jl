/**
 * Derive a numeric sort key from a parsed ordering value so the worksheet can
 * render members low→high (series-scoping.html: SERIES.sort by `key`). Pure;
 * never mutates row identity or the batch payload. Unparseable → null (those
 * rows sort last, stable).
 *
 * Rules, in order:
 *  - `a : b` ratio → the second term `b` (LL37 : lipid titration sorts by the
 *    lipid side, matching the mockup's `1 : 0 … 1 : 4`).
 *  - otherwise the first standalone number anywhere in the string.
 */
export function parseSortKey(value: string): number | null {
  const ratio = value.match(/^\s*-?\d*\.?\d+\s*:\s*(-?\d*\.?\d+)\s*$/);
  if (ratio) return Number(ratio[1]);
  const num = value.match(/-?\d*\.?\d+/);
  return num ? Number(num[0]) : null;
}
