import { describe, it, expect } from "vitest";
import * as ui from "../../src/components/ui";

describe("ui barrel public surface", () => {
  it("exports all 12 adopted primitives as functions", () => {
    expect(typeof ui.Button).toBe("function");
    expect(typeof ui.Card).toBe("function");
    expect(typeof ui.Dot).toBe("function");
    expect(typeof ui.HintText).toBe("function");
    expect(typeof ui.IconButton).toBe("function");
    expect(typeof ui.Kicker).toBe("function");
    expect(typeof ui.ModalShell).toBe("function");
    expect(typeof ui.PhaseChip).toBe("function");
    expect(typeof ui.PhaseStrip).toBe("function");
    expect(typeof ui.ScoreBar).toBe("function");
    expect(typeof ui.SegmentedControl).toBe("function");
    expect(typeof ui.ToastContainer).toBe("function");
  });

  it("does NOT export the deleted primitives", () => {
    expect("Input" in ui).toBe(false);
    expect("SectionLabel" in ui).toBe(false);
  });
});
