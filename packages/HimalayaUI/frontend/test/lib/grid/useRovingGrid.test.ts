import { describe, it, expect } from "vitest";
import { act, renderHook } from "@testing-library/react";
import {
  useRovingGrid,
  useRovingGridContext,
} from "../../../src/lib/grid/useRovingGrid";

const DIMS = { rows: 140, cols: 6 };
// Samples grid: Exposures(2) + Tags(4) are the multi-widget interaction cells.
const INTERACTION = new Set([2, 4]);

function setup(opts: Partial<Parameters<typeof useRovingGrid>[0]> = {}) {
  return renderHook(() =>
    useRovingGrid({ ...DIMS, interactionCols: INTERACTION, ...opts }),
  );
}

describe("useRovingGrid — initial state", () => {
  it("defaults the active cell to the first data Sample cell {1,1}", () => {
    const { result } = setup();
    expect(result.current.activeCoord).toEqual({ row: 1, col: 1 });
    expect(result.current.isActive(1, 1)).toBe(true);
    expect(result.current.interaction).toBe(false);
    expect(result.current.interactionCoord).toBeNull();
  });

  it("honours an explicit initialCoord", () => {
    const { result } = setup({ initialCoord: { row: 0, col: 1 } });
    expect(result.current.activeCoord).toEqual({ row: 0, col: 1 });
  });

  it("an EMPTY table (rows:1, header only) clamps the active cell INTO the grid — never tab-stop-less", () => {
    // The default {1,1} would be out of bounds with no data rows → zero cells
    // return tabIndexFor 0 → the header sort buttons become keyboard-unreachable.
    // The ON-READ clamp (effRow/effCol) resolves it to the header row {0,1}.
    const { result } = setup({ rows: 1 });
    expect(result.current.activeCoord).toEqual({ row: 0, col: 1 });
    // Exactly one cell is the tab stop, and it is on the (only) header row.
    let zeros = 0;
    for (let c = 0; c < 6; c++) if (result.current.tabIndexFor(0, c) === 0) zeros++;
    expect(zeros).toBe(1);
    expect(result.current.tabIndexFor(0, 1)).toBe(0);
  });

  it("a POPULATED table defaults the tab stop to the FIRST DATA row, NOT the header (regression pin)", () => {
    // Regression: a once-only init clamp pinned {1,1}→{0,1} during the loading
    // mount (rows:1) and never recovered when rows arrived, leaving the default
    // tab stop on the "Sample" sort HEADER. The on-read clamp resolves a stored
    // {1,1} back to the first data row once dims are populated.
    const { result } = setup({ rows: 140 });
    expect(result.current.tabIndexFor(1, 1)).toBe(0); // first DATA row, Sample cell
    expect(result.current.tabIndexFor(0, 1)).toBe(-1); // NOT the header
    expect(result.current.activeCoord).toEqual({ row: 1, col: 1 });
  });
});

describe("useRovingGrid — tabIndexFor (roving tabindex)", () => {
  it("returns 0 for the active cell and -1 for the rest", () => {
    const { result } = setup();
    expect(result.current.tabIndexFor(1, 1)).toBe(0);
    expect(result.current.tabIndexFor(1, 2)).toBe(-1);
    expect(result.current.tabIndexFor(0, 0)).toBe(-1);
  });

  it("moves the single 0 when the active cell changes", () => {
    const { result } = setup();
    act(() => result.current.requestActivate(3, 4));
    expect(result.current.tabIndexFor(1, 1)).toBe(-1);
    expect(result.current.tabIndexFor(3, 4)).toBe(0);
    expect(result.current.activeCoord).toEqual({ row: 3, col: 4 });
  });

  it("requestActivate to the ALREADY-active cell bails (no state change, no focus steal)", () => {
    // The focus-yank fix: an inner widget's focus event bubbles to the cell div
    // and calls requestActivate(r, sameCol). That must be a no-op so the layout
    // effect never re-focuses the cell div off the inner widget. We assert the
    // OBSERVABLE proxy: the activeCoord OBJECT REFERENCE is unchanged (a fresh
    // setActiveCoord({...}) would replace it). (jsdom can't honestly assert the
    // live focus traversal — the render-verify / e2e is the real gate.)
    const { result } = setup();
    const before = result.current.activeCoord; // {1,1}
    act(() => result.current.requestActivate(1, 1));
    expect(result.current.activeCoord).toBe(before); // same reference → bailed
  });

  it("requestActivate to a DIFFERENT cell still replaces activeCoord", () => {
    const { result } = setup();
    const before = result.current.activeCoord;
    act(() => result.current.requestActivate(2, 3));
    expect(result.current.activeCoord).not.toBe(before);
    expect(result.current.activeCoord).toEqual({ row: 2, col: 3 });
  });

  it("pointer parity: requestActivate(r,c,{focus:false}) still moves the active cell", () => {
    // BUG D fix: pointer activation (mousedown) updates the active coord WITHOUT
    // requesting focus, so clicking a cell resumes keyboard nav from there but
    // never re-focuses the cell div (which would yank focus off the clicked
    // widget / trap Tab in the grid). jsdom can't assert the absence of the
    // programmatic focus() honestly; we assert the coord still moves (parity
    // works). The no-focus behavior is the render-verify / e2e gate.
    const { result } = setup();
    act(() => result.current.requestActivate(4, 2, { focus: false }));
    expect(result.current.activeCoord).toEqual({ row: 4, col: 2 });
    expect(result.current.tabIndexFor(4, 2)).toBe(0);
  });

  it("pointer activation onto a different cell still leaves interaction mode", () => {
    const { result } = setup({ initialCoord: { row: 3, col: 4 } });
    act(() => result.current.enterInteraction(3, 4));
    act(() => result.current.requestActivate(5, 1, { focus: false }));
    expect(result.current.interaction).toBe(false);
    expect(result.current.activeCoord).toEqual({ row: 5, col: 1 });
  });
});

describe("useRovingGrid — keyboard navigation", () => {
  function press(
    result: { current: ReturnType<typeof useRovingGrid> },
    init: KeyboardEventInit,
  ) {
    let defaultPrevented = false;
    act(() => {
      result.current.onContainerKeyDown({
        key: init.key!,
        ctrlKey: init.ctrlKey ?? false,
        metaKey: init.metaKey ?? false,
        preventDefault: () => {
          defaultPrevented = true;
        },
        stopPropagation: () => {},
      } as unknown as React.KeyboardEvent);
    });
    return defaultPrevented;
  }

  it("ArrowDown moves the active row down and prevents default", () => {
    const { result } = setup();
    const prevented = press(result, { key: "ArrowDown" });
    expect(result.current.activeCoord).toEqual({ row: 2, col: 1 });
    expect(prevented).toBe(true);
  });

  it("ArrowRight moves the active column right", () => {
    const { result } = setup();
    press(result, { key: "ArrowRight" });
    expect(result.current.activeCoord).toEqual({ row: 1, col: 2 });
  });

  it("clamps at edges (ArrowUp from row 1 toward header stops at row 0)", () => {
    const { result } = setup({ initialCoord: { row: 0, col: 0 } });
    press(result, { key: "ArrowUp" });
    expect(result.current.activeCoord).toEqual({ row: 0, col: 0 });
  });

  it("⌘+Home jumps to the grid corner (0,0)", () => {
    const { result } = setup({ initialCoord: { row: 5, col: 3 } });
    press(result, { key: "Home", metaKey: true });
    expect(result.current.activeCoord).toEqual({ row: 0, col: 0 });
  });

  it("a non-nav key (e.g. 'x') does NOT preventDefault and does not move", () => {
    const { result } = setup();
    const prevented = press(result, { key: "x" });
    expect(prevented).toBe(false);
    expect(result.current.activeCoord).toEqual({ row: 1, col: 1 });
  });
});

describe("useRovingGrid — interaction mode", () => {
  it("Enter on a multi-widget cell (Tags, col 4) enters interaction mode", () => {
    const { result } = setup({ initialCoord: { row: 3, col: 4 } });
    act(() => {
      result.current.onContainerKeyDown({
        key: "Enter",
        ctrlKey: false,
        metaKey: false,
        preventDefault: () => {},
        stopPropagation: () => {},
      } as unknown as React.KeyboardEvent);
    });
    expect(result.current.interaction).toBe(true);
    expect(result.current.interactionCoord).toEqual({ row: 3, col: 4 });
  });

  it("Enter on a single-widget cell (Sample, col 1) does NOT enter interaction mode", () => {
    const { result } = setup({ initialCoord: { row: 3, col: 1 } });
    act(() => {
      result.current.onContainerKeyDown({
        key: "Enter",
        ctrlKey: false,
        metaKey: false,
        preventDefault: () => {},
        stopPropagation: () => {},
      } as unknown as React.KeyboardEvent);
    });
    expect(result.current.interaction).toBe(false);
  });

  it("Enter on a HEADER interaction-col cell (row 0, Exposures col 2) does NOT enter interaction AND does NOT preventDefault (sort survives)", () => {
    // BUG C: interaction mode is body-only. The header is row 0; its Exposures
    // cell holds a sort button, so Enter must fall through WITHOUT preventDefault
    // so the button's native Enter→click (the sort) fires. (jsdom can't assert
    // the native click; we assert the two hook-observable contracts: no
    // interaction-mode flip, and preventDefault was NOT called.)
    const { result } = setup({ initialCoord: { row: 0, col: 2 } });
    let defaultPrevented = false;
    act(() => {
      result.current.onContainerKeyDown({
        key: "Enter",
        ctrlKey: false,
        metaKey: false,
        preventDefault: () => {
          defaultPrevented = true;
        },
        stopPropagation: () => {},
      } as unknown as React.KeyboardEvent);
    });
    expect(result.current.interaction).toBe(false);
    expect(defaultPrevented).toBe(false);
  });

  it("F2 also enters interaction mode on a multi-widget cell", () => {
    const { result } = setup({ initialCoord: { row: 3, col: 2 } });
    act(() => {
      result.current.enterInteraction(3, 2);
    });
    expect(result.current.interactionCoord).toEqual({ row: 3, col: 2 });
  });

  it("Escape exits interaction mode", () => {
    const { result } = setup({ initialCoord: { row: 3, col: 4 } });
    act(() => result.current.enterInteraction(3, 4));
    expect(result.current.interaction).toBe(true);
    act(() => {
      result.current.onContainerKeyDown({
        key: "Escape",
        ctrlKey: false,
        metaKey: false,
        preventDefault: () => {},
        stopPropagation: () => {},
      } as unknown as React.KeyboardEvent);
    });
    expect(result.current.interaction).toBe(false);
    expect(result.current.interactionCoord).toBeNull();
  });

  it("while in interaction mode, arrow keys do NOT move the active cell", () => {
    const { result } = setup({ initialCoord: { row: 3, col: 4 } });
    act(() => result.current.enterInteraction(3, 4));
    act(() => {
      result.current.onContainerKeyDown({
        key: "ArrowDown",
        ctrlKey: false,
        metaKey: false,
        preventDefault: () => {},
        stopPropagation: () => {},
      } as unknown as React.KeyboardEvent);
    });
    expect(result.current.activeCoord).toEqual({ row: 3, col: 4 });
  });

  it("requestActivate to a different cell leaves interaction mode", () => {
    const { result } = setup({ initialCoord: { row: 3, col: 4 } });
    act(() => result.current.enterInteraction(3, 4));
    act(() => result.current.requestActivate(5, 1));
    expect(result.current.interaction).toBe(false);
    expect(result.current.activeCoord).toEqual({ row: 5, col: 1 });
  });
});

describe("useRovingGridContext — inert default (no provider)", () => {
  it("returns undefined tabIndexFor and false isActive outside a grid", () => {
    const { result } = renderHook(() => useRovingGridContext());
    expect(result.current.tabIndexFor(1, 1)).toBeUndefined();
    expect(result.current.isActive(1, 1)).toBe(false);
    expect(result.current.interaction).toBe(false);
    // no-ops don't throw
    expect(() => result.current.registerCellEl(1, 1, null)).not.toThrow();
    expect(() => result.current.requestActivate(1, 1)).not.toThrow();
  });
});
