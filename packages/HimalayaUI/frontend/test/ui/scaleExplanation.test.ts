import { describe, expect, it } from "vitest";
import { explainScale, SCALE_REASONS } from "../../src/ui/scaleExplanation";
import type { CompareScale } from "../../src/workspace/compare";

const sources = Object.keys(SCALE_REASONS) as CompareScale["source"][];

describe("explainScale", () => {
  it("says a pair divided by I₂ is not rescaled", () => {
    expect(explainScale({ scale: 1, source: "normalised" })).toBe(
      "Not rescaled. Both traces are already divided by the transmitted beam (I₂), so they share a scale.",
    );
  });

  it("gives the flux ratio as the factor for raw counts", () => {
    expect(explainScale({ scale: 0.98234, source: "monitor" })).toBe(
      "Multiplied by 0.982, the ratio of the transmitted beam (I₂) the two images received.",
    );
  });

  it("says what a high-q match does before it says why", () => {
    expect(explainScale({ scale: 1.1, source: "capillary" })).toBe(
      "Multiplied by 1.10 so the two traces agree at high q, leaving out peaks. Each capillary has its own path length, so the transmitted beam can't set the scale here.",
    );
  });

  it.each(sources.filter((s) => s !== "normalised" && s !== "monitor"))("opens every %s match with the factor the chip shows", (source) => {
    expect(explainScale({ scale: 1.5, source })).toMatch(/^Multiplied by 1\.50 so the two traces agree at high q, leaving out peaks\. [A-Z].*\.$/);
  });

  it("gives every source its own reason", () => {
    expect(new Set(sources.map((s) => SCALE_REASONS[s])).size).toBe(sources.length);
  });
});
