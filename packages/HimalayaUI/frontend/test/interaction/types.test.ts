import { describe, it, expect } from "vitest";
import type { Action, ListCursor, CursorStepperProps } from "../../src/print/interaction/types";

describe("interaction types", () => {
  it("an Action is constructable with the documented fields", () => {
    const a: Action = {
      id: "openFocus",
      label: "Focus",
      keys: ["Enter"],
      group: "Navigate",
      enabled: () => true,
      run: () => {},
      dock: "primary",
    };
    expect(a.id).toBe("openFocus");
  });

  it("a ListCursor exposes the parity contract", () => {
    const stepper: CursorStepperProps = {
      label: "Sample", axis: "vertical", testIdBase: "sample",
      count: "1 / 3", onPrev: () => {}, onNext: () => {},
      prevDisabled: true, nextDisabled: false,
    };
    const c: ListCursor = {
      cursorId: null, selected: new Set<number>(),
      setCursor: () => {}, moveBy: () => {}, activate: () => {},
      toggleSelect: () => {}, rowProps: () => ({} as never),
      stepperProps: () => stepper,
    };
    expect(c.cursorId).toBeNull();
  });
});
