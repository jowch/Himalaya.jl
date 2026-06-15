import { useShortcuts } from "./useShortcuts";

/**
 * The unified Alt+↑/↓ reorder power-gesture, bound through the shared shortcut
 * registry (`reorderUp` / `reorderDown` → `Alt+ArrowUp` / `Alt+ArrowDown`). It
 * moves whichever reorderable row currently holds keyboard focus, so a keyboard
 * user can reorder from ANY control inside a row — not only its drag grip /
 * ▲▼ buttons. One gesture means the same thing on every list (Scoping worksheet,
 * Series builder).
 *
 * The focused row is resolved from the key event's `target` (the focused element
 * for a window keydown) via `rowSelector`; its visual position is read from
 * `data-reorder-index`. `move(index, delta, row)` performs the surface's
 * EXISTING reorder (Scoping `moveRow`, Builder's ▲▼ buttons), so focus-retention
 * and screen-reader announcements stay owned by the surface that already gets
 * them right. When focus is not inside a reorderable row the binding DECLINES
 * (returns `false`), leaving Alt+Arrow its native meaning (word-nav in inputs,
 * etc.).
 */
export function useReorderShortcuts(opts: {
  /** Selector for the reorderable row carrying `data-reorder-index`. */
  rowSelector: string;
  /** Perform the surface's reorder of the row at `index` by `delta` (-1 up / +1 down). */
  move: (index: number, delta: -1 | 1, row: HTMLElement) => void;
  enabled?: boolean;
}): void {
  const { rowSelector, move, enabled = true } = opts;
  const step = (e: KeyboardEvent, delta: -1 | 1): boolean => {
    const t = e.target;
    if (!(t instanceof HTMLElement)) return false;
    const row = t.closest<HTMLElement>(rowSelector);
    if (!row) return false;
    const index = Number(row.getAttribute("data-reorder-index"));
    if (!Number.isInteger(index)) return false;
    move(index, delta, row);
    return true;
  };
  useShortcuts(
    {
      reorderUp: (e) => (step(e, -1) ? undefined : false),
      reorderDown: (e) => (step(e, 1) ? undefined : false),
    },
    enabled,
  );
}
