import {
  createContext,
  useContext,
  useCallback,
  useLayoutEffect,
  useMemo,
  useRef,
  useState,
} from "react";
import type { KeyboardEvent, RefObject } from "react";
import { nextGridCoord, type GridCoord } from "./rovingGrid";

/**
 * Roving-tabindex grid (WAI-ARIA "Data Grid").
 *
 * `SheetTable` owns ONE active cell; every other cell carries `tabIndex={-1}`,
 * so the whole table is a single tab stop. Arrow / Home / End / Page keys move
 * the active cell through `nextGridCoord` (the pure reducer in rovingGrid.ts);
 * landing focuses the active cell's registered element. Enter / F2 on a
 * multi-widget cell enters "interaction mode" (the cell's own focusables become
 * naturally tabbable; Escape exits back to navigation mode).
 *
 * The focus() call is gated by a ref flag set ONLY by a user-driven coord
 * change (key handler / requestActivate / enterInteraction / exitInteraction)
 * and consumed in a useLayoutEffect — exactly the RepresentativeBox discipline,
 * so SSE / foreign re-renders never yank focus. Re-applying tabindexes every
 * render is idempotent and fine; only focus() is flag-gated.
 *
 * The context default is INERT: `tabIndexFor` returns `undefined`,
 * register/activate are no-ops, `isActive` is false. So a `SampleTableRow`
 * rendered outside a grid (its story, or any non-grid consumer) renders exactly
 * as before — no roving. The grid is "on" only when `SheetTable` supplies the
 * provider via `useRovingGrid`.
 */

export interface RovingGridValue {
  /** The active cell coordinate (row 0 = header). */
  activeCoord: GridCoord;
  /** Is (row,col) the active cell? */
  isActive: (row: number, col: number) => boolean;
  /** Roving tabindex for a cell/widget: `0` for the active cell, `-1` for the
   *  rest. Returns `undefined` from the INERT default (no grid) so a non-grid
   *  consumer renders unchanged. */
  tabIndexFor: (row: number, col: number) => 0 | -1 | undefined;
  /** Make (row,col) the active cell. By default also requests focus (keyboard /
   *  programmatic activation). POINTER parity passes `{ focus: false }`: the
   *  native click already focused the clicked element, so the active coord moves
   *  WITHOUT a programmatic re-focus that would yank focus to the cell div. */
  requestActivate: (row: number, col: number, opts?: { focus?: boolean }) => void;
  /** Register a cell's focus target element (the inner single widget, or the
   *  gridcell div for zero/multi-widget cells) under its (row,col). Pass `null`
   *  on unmount to drop it. Returns nothing — wire as a ref callback. */
  registerCellEl: (row: number, col: number, el: HTMLElement | null) => void;
  /** Container keydown handler — attach to the `role="grid"` element. */
  onContainerKeyDown: (e: KeyboardEvent) => void;
  /** Is the grid in interaction mode (a multi-widget cell entered)? */
  interaction: boolean;
  /** The cell currently in interaction mode (null when in navigation mode). */
  interactionCoord: GridCoord | null;
  /** Enter interaction mode on a multi-widget cell (Enter / F2). */
  enterInteraction: (row: number, col: number) => void;
  /** Exit interaction mode (Escape), returning focus to the gridcell. */
  exitInteraction: () => void;
  /** The grid container ref — attach to the `role="grid"` element. */
  gridRef: RefObject<HTMLDivElement>;
}

const key = (row: number, col: number): string => `${row},${col}`;

const INERT: RovingGridValue = {
  activeCoord: { row: 0, col: 0 },
  isActive: () => false,
  tabIndexFor: () => undefined,
  requestActivate: () => {},
  registerCellEl: () => {},
  onContainerKeyDown: () => {},
  interaction: false,
  interactionCoord: null,
  enterInteraction: () => {},
  exitInteraction: () => {},
  // A stable inert ref so consumers can read `.current` (always null here).
  gridRef: { current: null },
};

const RovingGridContext = createContext<RovingGridValue>(INERT);

/** Read the roving-grid context. Returns the INERT value outside a provider. */
export function useRovingGridContext(): RovingGridValue {
  return useContext(RovingGridContext);
}

export const RovingGridProvider = RovingGridContext.Provider;

export interface UseRovingGridOptions {
  /** Total rows INCLUDING the header row at index 0. */
  rows: number;
  /** Total columns. */
  cols: number;
  /** Columns (by index) that are MULTI-WIDGET cells — Enter/F2 enters
   *  interaction mode rather than just activating a single control. For the
   *  Samples grid: Exposures (2) and Tags (4). */
  interactionCols?: ReadonlySet<number>;
  /** Initial active cell. Default = the first data Sample cell `{row:1,col:1}`
   *  (the first naturally-focusable cell; the checkbox at col 0 only exists when
   *  a checkbox column is present, and Sample is the row identity). */
  initialCoord?: GridCoord;
}

/**
 * The grid owner hook. Instantiated once by `SheetTable` when `roving` is on;
 * the returned value is fed to `RovingGridProvider` so every descendant cell
 * reads its tabindex / registers its element / reports clicks by (row,col).
 */
export function useRovingGrid(opts: UseRovingGridOptions): RovingGridValue {
  const { rows, cols, interactionCols, initialCoord } = opts;
  // Plain default — do NOT clamp at init: SheetTable mounts during the loading
  // window with rows:1 (header only), so a once-only clamp would pin {1,1}→{0,1}
  // (header) and never recover when the data rows arrive (a regression that put
  // the default tab stop on the "Sample" sort button instead of the first data
  // row). The clamp instead happens ON READ each render (effRow/effCol below).
  const [activeCoord, setActiveCoord] = useState<GridCoord>(
    initialCoord ?? { row: 1, col: 1 },
  );
  const [interactionCoord, setInteractionCoord] = useState<GridCoord | null>(null);

  const maxRow = Math.max(0, rows - 1);
  const maxCol = Math.max(0, cols - 1);
  // Effective active cell = stored coord clamped into the CURRENT dims. While the
  // table is still loading (rows:1, header only) a stored {1,1} resolves to the
  // header {0,1} so the grid always has exactly one tab stop; once the data rows
  // arrive it resolves back to {1,1} (first data row) with NO stored-state change.
  // This is both the empty-table fix AND avoids the loading-mount race. The WRITE
  // paths keep operating on the raw `activeCoord` — nextGridCoord clamps its
  // results to dims, so stored state stays in-range after any user navigation;
  // the only stored/effective divergence is the transient loading window.
  const effRow = Math.min(activeCoord.row, maxRow);
  const effCol = Math.min(activeCoord.col, maxCol);

  const gridRef = useRef<HTMLDivElement>(null);
  const cellEls = useRef<Map<string, HTMLElement>>(new Map());
  // Focus-request flag — set ONLY by a user-driven coord change, consumed in the
  // layout effect. SSE / foreign re-renders never set it, so they never steal
  // focus (RepresentativeBox discipline).
  const wantFocus = useRef(false);
  // Latest activeCoord, mirrored into a ref each render so the STABLE
  // (`useCallback([])`) requestActivate can read it without being re-created on
  // every coord change. This is what lets requestActivate BAIL when the target
  // cell is already active — see the focus-yank fix below.
  const coordRef = useRef(activeCoord);
  coordRef.current = activeCoord;

  const dims = useMemo(() => ({ rows, cols }), [rows, cols]);

  const registerCellEl = useCallback(
    (row: number, col: number, el: HTMLElement | null): void => {
      const k = key(row, col);
      if (el === null) cellEls.current.delete(k);
      else cellEls.current.set(k, el);
    },
    [],
  );

  // Read against the EFFECTIVE (clamped) coord, not the raw stored one, so a
  // stored coord that is transiently out of bounds (loading window) still yields
  // exactly one in-grid tab stop. Deps are effRow/effCol so these recompute when
  // dims change (data loads) — not just when activeCoord changes.
  const isActive = useCallback(
    (row: number, col: number): boolean => row === effRow && col === effCol,
    [effRow, effCol],
  );

  const tabIndexFor = useCallback(
    (row: number, col: number): 0 | -1 | undefined =>
      row === effRow && col === effCol ? 0 : -1,
    [effRow, effCol],
  );

  const requestActivate = useCallback(
    (row: number, col: number, opts?: { focus?: boolean }): void => {
      // BAIL when the cell is already active (no-op): drops redundant re-renders.
      if (coordRef.current.row === row && coordRef.current.col === col) return;
      // POINTER parity passes { focus: false }: a click already focuses the right
      // element natively, so we update the active coord WITHOUT setting wantFocus
      // — otherwise the layout effect would re-focus the cell div and yank focus
      // off the thumbnail/widget the user just clicked. Keyboard / programmatic
      // activations (default) DO want focus moved to the new active cell.
      if (opts?.focus !== false) wantFocus.current = true;
      // Activating a DIFFERENT cell leaves any interaction mode (the user moved
      // away). Same-cell already bailed above.
      setInteractionCoord((prev) =>
        prev && (prev.row !== row || prev.col !== col) ? null : prev,
      );
      setActiveCoord({ row, col });
    },
    [],
  );

  const enterInteraction = useCallback((row: number, col: number): void => {
    wantFocus.current = true;
    setActiveCoord({ row, col });
    setInteractionCoord({ row, col });
  }, []);

  const exitInteraction = useCallback((): void => {
    wantFocus.current = true;
    setInteractionCoord(null);
  }, []);

  const onContainerKeyDown = useCallback(
    (e: KeyboardEvent): void => {
      // In interaction mode the grid's cross-cell arrow handler is suppressed —
      // let Tab and the cell's own keys work. Escape exits.
      if (interactionCoord !== null) {
        if (e.key === "Escape") {
          e.preventDefault();
          e.stopPropagation();
          exitInteraction();
        }
        return;
      }

      // Enter / F2 on a multi-widget BODY cell enters interaction mode.
      // Interaction mode is a body-cell affordance ONLY: the header is always
      // row 0, and its cells (even at an interaction column like Exposures/Tags)
      // hold a sort button or a static label — never multi-widget content. Guard
      // on row >= 1 so Enter on a sortable header is NOT swallowed here: it falls
      // through to nextGridCoord (→ null for "Enter", no preventDefault) and the
      // header button's native Enter→click fires (the sort).
      if (
        (e.key === "Enter" || e.key === "F2") &&
        activeCoord.row >= 1 &&
        interactionCols?.has(activeCoord.col)
      ) {
        e.preventDefault();
        enterInteraction(activeCoord.row, activeCoord.col);
        return;
      }

      const ctrl = e.ctrlKey || e.metaKey;
      const next = nextGridCoord(activeCoord, e.key, dims, { ctrl });
      if (next === null) return; // not a nav key — let it bubble (⌘K, X/K/Esc cull)
      // preventDefault ONLY on keys the grid consumes.
      e.preventDefault();
      if (next.row !== activeCoord.row || next.col !== activeCoord.col) {
        wantFocus.current = true;
        setActiveCoord(next);
      }
    },
    [activeCoord, dims, interactionCoord, interactionCols, enterInteraction, exitInteraction],
  );

  // Focus the active cell's element — but ONLY after a user-driven change
  // (wantFocus). Re-applying nothing here is fine; tabindexes are re-applied by
  // each cell's render. Before paint so focus never flashes.
  //
  // When the active cell IS the interaction cell, focus its FIRST naturally
  // tabbable descendant instead of the gridcell div (the design note: interaction
  // mode focuses "the first (or selected) inner widget"). The child cells'
  // useMultiWidgetInertness layout effects run BEFORE this parent effect, so by
  // now the first inner widget already carries tabIndex=0 and the rest are
  // reachable via Tab.
  useLayoutEffect(() => {
    if (!wantFocus.current) return;
    wantFocus.current = false;
    // Focus the EFFECTIVE active cell (clamped into current dims).
    const el = cellEls.current.get(key(effRow, effCol));
    if (el == null) return;
    const interacting =
      interactionCoord?.row === effRow && interactionCoord?.col === effCol;
    if (interacting) {
      const inner = el.querySelector<HTMLElement>(
        '[tabindex="0"], button:not([tabindex="-1"]), a[href], input:not([tabindex="-1"])',
      );
      (inner ?? el).focus();
    } else {
      el.focus();
    }
  }, [effRow, effCol, interactionCoord]);

  return useMemo(
    () => ({
      // Expose the EFFECTIVE (clamped) coord so any consumer reading it sees an
      // in-grid cell, consistent with tabIndexFor/isActive.
      activeCoord: { row: effRow, col: effCol },
      isActive,
      tabIndexFor,
      requestActivate,
      registerCellEl,
      onContainerKeyDown,
      interaction: interactionCoord !== null,
      interactionCoord,
      enterInteraction,
      exitInteraction,
      gridRef,
    }),
    [
      effRow,
      effCol,
      isActive,
      tabIndexFor,
      requestActivate,
      registerCellEl,
      onContainerKeyDown,
      interactionCoord,
      enterInteraction,
      exitInteraction,
    ],
  );
}
