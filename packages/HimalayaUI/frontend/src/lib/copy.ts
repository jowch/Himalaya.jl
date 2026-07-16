/**
 * Copy sanitizers for user-facing strings assembled from UPSTREAM data.
 *
 * The house copy law bans em dashes (and "--"), enforced at the source level by
 * test/print-copy-no-emdash.test.ts. That guard cannot see data that arrives at
 * runtime — an experiment named "SSRL April 2026 — 1p7m" in the DB would render
 * an em dash no static scan can catch. Run such strings through `sanitizeDashes`
 * before they reach copy.
 */

/**
 * Fold em/en dashes into the middot facet separator " · " (the subtitle's own
 * chip-separator voice), collapsing the surrounding whitespace. Idempotent.
 */
export function sanitizeDashes(s: string): string {
  return s
    .replace(/\s*[—–]\s*/g, " · ")
    .replace(/\s{2,}/g, " ")
    .trim();
}
