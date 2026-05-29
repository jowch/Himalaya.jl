import { describe, it, expect } from "vitest";
import * as ui from "../../src/components/ui";

describe("ui barrel public surface", () => {
  it("exports the adopted primitives", () => {
    expect(typeof ui.Button).toBe("function");
    expect(typeof ui.Dot).toBe("function");
    expect(typeof ui.HintText).toBe("function");
    expect(typeof ui.ScoreBar).toBe("function");
    expect(typeof ui.ToastContainer).toBe("function");
  });

  it("does NOT export the deleted primitives", () => {
    expect("Input" in ui).toBe(false);
    expect("SectionLabel" in ui).toBe(false);
  });
});
