/**
 * Shared cursor-contract harness.
 *
 * Usage: import `runCursorContract` and call it at module scope from any
 * `*.contract.test.tsx` file.  The caller supplies a `makeWrapper` factory that
 * returns a `ui(capture)` render function.  The `Wrapper` component the caller
 * defines must:
 *   - call `useListCursor` (or equivalent) to obtain the cursor
 *   - spread `cursor.rowProps(id)` onto every row element
 *   - render a `<DockStepper {...cursor.stepperProps()} />`
 *   - call `capture({ cursor, ids, onActivate })` **in the render body** (not
 *     in a useEffect) so the harness always holds the latest cursor object
 *
 * Stale-cursor hazard: because `useListCursor` returns a new object on every
 * render, reading the cursor from a closure captured before an `act()` will
 * see stale values.  The solution: store the latest cursor in a `let m` that
 * is reassigned by `capture` on every render.  All assertions read `m.cursor`
 * AFTER the `act()` flush.
 */
import { describe, it, expect } from "vitest";
import type { Mock } from "vitest";
import { render, fireEvent, act } from "@testing-library/react";
import type { ListCursor } from "../../src/print/interaction/types";

export interface Mounted {
  cursor: ListCursor;
  ids: number[];
  onActivate: Mock;
}

export function runCursorContract(
  name: string,
  makeWrapper: () => {
    ui: (capture: (m: Mounted) => void) => JSX.Element;
  },
): void {
  describe(`cursor contract: ${name}`, () => {
    // ── 1. Roving tabindex ──────────────────────────────────────────────────
    it("roving: exactly one [role='row'][tabindex='0'] at all times", () => {
      let m!: Mounted;
      const { ui } = makeWrapper();
      const { container } = render(ui((captured) => { m = captured; }));

      // Initial render: exactly one tabbable row
      expect(container.querySelectorAll('[role="row"][tabindex="0"]').length).toBe(1);

      // After moving the cursor: still exactly one
      act(() => { m.cursor.moveBy(1); });
      expect(container.querySelectorAll('[role="row"][tabindex="0"]').length).toBe(1);
    });

    // ── 2. Input parity ─────────────────────────────────────────────────────
    it("input parity: arrow == stepper-next == row-click all land on ids[1]", () => {
      let m!: Mounted;
      const { ui } = makeWrapper();
      const { container, getByTestId } = render(ui((captured) => { m = captured; }));

      // Path 1: moveBy(1) from ids[0] → ids[1]
      act(() => { m.cursor.moveBy(1); });
      const viaArrow = m.cursor.cursorId;

      // Reset to ids[0]
      act(() => { m.cursor.setCursor(m.ids[0]!); });

      // Path 2: stepper "next" button → ids[1]
      // Read testIdBase from the live cursor (it reflects the latest render)
      const nextTestId = `dock-next-${m.cursor.stepperProps().testIdBase}`;
      act(() => { fireEvent.click(getByTestId(nextTestId)); });
      const viaStepper = m.cursor.cursorId;

      // Reset to ids[0]
      act(() => { m.cursor.setCursor(m.ids[0]!); });

      // Path 3: click row at DOM index 1 (which is ids[1])
      const rows = container.querySelectorAll<HTMLElement>('[role="row"]');
      act(() => { fireEvent.click(rows[1]!); });
      const viaClick = m.cursor.cursorId;

      expect(viaArrow).toBe(m.ids[1]);
      expect(viaStepper).toBe(m.ids[1]);
      expect(viaClick).toBe(m.ids[1]);
    });

    // ── 3. Activate parity ──────────────────────────────────────────────────
    it("activate parity: cursor.activate() == double-click both fire onActivate with cursorId", () => {
      let m!: Mounted;
      const { ui } = makeWrapper();
      const { container } = render(ui((captured) => { m = captured; }));

      // Park cursor on ids[1]
      act(() => { m.cursor.setCursor(m.ids[1]!); });

      // Path 1: keyboard Enter via cursor.activate()
      act(() => { m.cursor.activate(); });
      const viaEnter = m.onActivate.mock.calls.at(-1)?.[0] as number | undefined;

      // Path 2: double-click on the row at DOM index 1 (ids[1])
      const rows = container.querySelectorAll<HTMLElement>('[role="row"]');
      act(() => { fireEvent.doubleClick(rows[1]!); });
      const viaDouble = m.onActivate.mock.calls.at(-1)?.[0] as number | undefined;

      expect(viaEnter).toBe(m.ids[1]);
      expect(viaDouble).toBe(m.ids[1]);
    });

    // ── 4. Orthogonality ────────────────────────────────────────────────────
    it("orthogonality: toggleSelect never moves cursor; moveBy never clears selection", () => {
      let m!: Mounted;
      const { ui } = makeWrapper();
      render(ui((captured) => { m = captured; }));

      // Selecting ids[2] must not move the cursor from ids[0]
      act(() => { m.cursor.toggleSelect(m.ids[2]!); });
      expect(m.cursor.selected.has(m.ids[2]!)).toBe(true);
      expect(m.cursor.cursorId).toBe(m.ids[0]!);

      // Moving cursor to ids[1] must not clear the selection
      act(() => { m.cursor.moveBy(1); });
      expect(m.cursor.selected.has(m.ids[2]!)).toBe(true);
      expect(m.cursor.cursorId).toBe(m.ids[1]!);

      // Toggling another id must not move the cursor
      act(() => { m.cursor.toggleSelect(m.ids[0]!); });
      expect(m.cursor.cursorId).toBe(m.ids[1]!);
    });
  });
}
