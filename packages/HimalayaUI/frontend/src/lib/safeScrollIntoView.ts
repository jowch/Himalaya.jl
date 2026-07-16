/**
 * scrollIntoView guarded against jsdom's "Not implemented" throw — jsdom has no
 * layout engine, and the scroll is cosmetic-only, so a failure is swallowed.
 */
export function safeScrollIntoView(
  el: Element | null | undefined,
  opts: ScrollIntoViewOptions = { block: "nearest" },
): void {
  try {
    el?.scrollIntoView(opts);
  } catch {
    /* no layout engine (jsdom) — cosmetic only */
  }
}
