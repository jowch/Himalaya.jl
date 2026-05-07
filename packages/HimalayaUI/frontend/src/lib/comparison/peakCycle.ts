/**
 * Peak click cycle for the Compare page edit-mode plot
 * (Plan §Phase 8, Task 8.1; spec §Peak click semantics).
 *
 * The triangle markers on each member's trace are clickable in edit mode.
 * A regular click cycles the peak's display state for that member:
 *
 *   shown (default) → labeled → hidden → shown
 *
 * `alt+click` (modifier) jumps directly to `hidden` regardless of starting
 * state — convenience for users with many peaks who want to bulk-suppress.
 *
 * State is encoded into the member's `peak_display` JSON column:
 *
 *   { hidden: number[]; labeled: number[] }
 *
 * "shown" is the implicit default — a peak id absent from BOTH arrays is
 * shown. The pure function returns a brand-new object (with brand-new
 * arrays) so React/Zustand reference-equality dirty checks fire correctly.
 *
 * Identity is by peak id, not q-value. If a peak is reanalyzed and its id
 * disappears, the rendering falls back to "shown" gracefully (the lookup
 * just misses both sets) — no migration needed at the comparison layer.
 *
 * Pure / no-op safe: passing `undefined` or `null` for the prior state
 * works as the implicit "shown" baseline.
 */

export interface PeakDisplay {
  hidden: number[];
  labeled: number[];
}

/**
 * Compute the next peak_display state for a single peak click.
 *
 * @param prior     The current display state for the member, or undefined
 *                  when never set (everything implicitly shown).
 * @param peakId    The peak id that was clicked.
 * @param altKey    True for alt+click (jump-to-hidden); false for the
 *                  three-step cycle.
 */
export function cyclePeakDisplay(
  prior: PeakDisplay | null | undefined,
  peakId: number,
  altKey: boolean,
): PeakDisplay {
  const hidden = new Set<number>(prior?.hidden ?? []);
  const labeled = new Set<number>(prior?.labeled ?? []);

  if (altKey) {
    // Jump-to-hidden: ensure peak is in `hidden` and not in `labeled`,
    // regardless of starting state. Idempotent on already-hidden peaks.
    labeled.delete(peakId);
    hidden.add(peakId);
  } else {
    // Three-step cycle. Determine current bucket:
    //   - in labeled  → labeled → hidden
    //   - in hidden   → hidden  → shown
    //   - else (shown) → shown  → labeled
    if (labeled.has(peakId)) {
      labeled.delete(peakId);
      hidden.add(peakId);
    } else if (hidden.has(peakId)) {
      hidden.delete(peakId);
      // (already removed from labeled above implicitly)
    } else {
      labeled.add(peakId);
    }
  }

  return {
    hidden: Array.from(hidden),
    labeled: Array.from(labeled),
  };
}
