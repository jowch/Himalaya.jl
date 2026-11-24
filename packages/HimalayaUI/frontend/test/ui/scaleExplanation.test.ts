import { describe, expect, it } from "vitest";
import { explainScale, SCALE_REASONS } from "../../src/ui/scaleExplanation";
import type { CompareScale } from "../../src/workspace/compare";

const matched = Object.keys(SCALE_REASONS) as Exclude<CompareScale["source"], "normalised" | "monitor">[];

describe("explainScale", () => {
  it("explains why a pair divided by I₂ in the same cell needs no rescaling", () => {
    expect(explainScale({ scale: 1, source: "normalised" })).toBe(
      "Both traces are divided by I₂, the beam transmitted through the sample. Every sample fills the same autosampler cell, so the beam crosses the same thickness of sample and the traces are already on one scale.",
    );
  });

  it("gives the transmitted-beam ratio as the factor for raw counts", () => {
    expect(explainScale({ scale: 0.98234, source: "monitor" })).toBe(
      "Scattered intensity grows with the beam that gets through the sample. This trace is scaled ×0.982, the ratio of that transmitted beam (I₂ × exposure time) between the two images.",
    );
  });

  it("explains the physics of a high-q match before the factor", () => {
    expect(explainScale({ scale: 1.1, source: "capillary" })).toBe(
      "Scattered intensity grows with how much material the beam crosses: the path length and how evenly the sample fills it. Both change from capillary to capillary, and I₂ corrects only for transmission. This trace is scaled ×1.10 to match at high q, where the solvent dominates. A sample that still scatters there, such as vesicles, throws this factor off.",
    );
  });

  it.each(matched)("ends every %s match with the factor the chip shows", (source) => {
    expect(explainScale({ scale: 1.5, source })).toMatch(/^[A-ZI].*\. This trace is scaled ×1\.50 to match at high q, where the solvent dominates\. A sample that still scatters there, such as vesicles, throws this factor off\.$/);
  });

  it("names the background when it is the one scaled", () => {
    expect(explainScale({ scale: 1.5, source: "beam_off" }, "The background")).toContain("The background is scaled ×1.50");
  });

  it("gives every high-q match its own reason", () => {
    expect(new Set(matched.map((s) => SCALE_REASONS[s])).size).toBe(matched.length);
  });
});
