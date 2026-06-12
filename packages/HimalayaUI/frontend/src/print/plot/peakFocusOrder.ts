// peakFocusOrder.ts — PURE (no React) helper for keyboard focus re-anchoring
// after a destructive peak edit (WCAG 2.4.3 Focus Order).
//
// When a keyboard user activates a peak mark to REMOVE it, that mark unmounts on
// the next render and browser focus falls to <body>. To preserve the user's
// place we move focus to the nearest SURVIVING peak: prefer the next peak by
// q-order; if the removed one was the last, fall back to the previous; if no
// peak survives, the caller focuses the "+ Peak" toolbar button instead.
//
// This is a pure ordering decision so it can be unit-tested in isolation — the
// real-browser focus move (querySelector + .focus()) is the caller's job and
// needs a live render-verify (jsdom focus-after-unmount is unreliable).

export interface FocusOrderPeak {
  id: number;
  q: number;
}

/**
 * Given the peak list as it stands BEFORE the removal and the id being removed,
 * return the id of the peak that should receive focus afterwards, or null if no
 * peak will survive (the caller should then fall back to the "+ Peak" button).
 *
 * "Nearest surviving sibling" is resolved in q-order: the next-higher-q peak
 * wins; if the removed peak was the highest q, the next-lower-q peak wins.
 */
export function nextFocusPeakId(
  peaks: ReadonlyArray<FocusOrderPeak>,
  removedId: number,
): number | null {
  // Order by q so "next / previous" means the visual left→right neighbour, not
  // array insertion order (manual peaks land out of q-order in the array).
  const ordered = [...peaks].sort((a, b) => a.q - b.q);
  const idx = ordered.findIndex((p) => p.id === removedId);
  if (idx === -1) return null;
  const next = ordered[idx + 1];
  if (next) return next.id;
  const prev = ordered[idx - 1];
  if (prev) return prev.id;
  return null;
}
