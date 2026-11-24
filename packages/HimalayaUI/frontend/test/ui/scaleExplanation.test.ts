import { describe, expect, it } from "vitest";
import { explainScale, SCALE_REASONS } from "../../src/ui/scaleExplanation";
import type { CompareScale } from "../../src/workspace/compare";

const matched = Object.keys(SCALE_REASONS) as Exclude<CompareScale["source"], "normalised" | "monitor" | "beam_off">[];
const HIGH_Q = "to match at high q, where both traces should be mostly solvent and capillary scattering. A sample that still scatters there, from a bilayer's form factor or just more lipid, throws this factor off.";

describe("explainScale", () => {
  it("explains why a pair divided by I₂ in the same cell needs no rescaling", () => {
    expect(explainScale({ scale: 1, source: "normalised" })).toBe(
      "Both traces are divided by I₂, the beam transmitted through the sample. Every sample fills the same autosampler cell, so the beam crosses the same thickness of sample and the traces are already on one scale.",
    );
  });

  it("explains why the transmitted beam alone sets a raw pair's scale, without a formula the code doesn't use", () => {
    const text = explainScale({ scale: 0.98234, source: "monitor" });
    expect(text).toBe(
      "Scattered intensity grows with the beam that gets through the sample. Every sample fills the same autosampler cell, so the thickness matches and the beam alone sets the scale. This trace is scaled ×0.982, the ratio of the beam that got through each sample while it was measured.",
    );
  });

  it("explains the physics of a high-q match before the factor", () => {
    expect(explainScale({ scale: 1.1, source: "capillary" })).toBe(
      `Scattered intensity grows with how much material the beam crosses: the path length and how evenly the sample fills it. Both change from capillary to capillary. Dividing by I₂ corrects for the beam and for absorption, not for how much material the beam crossed. This trace is scaled ×1.10 ${HIGH_Q}`,
    );
  });

  it.each(matched)("ends every %s match with the factor the chip shows and the high-q assumption", (source) => {
    const text = explainScale({ scale: 1.5, source });
    expect(text.startsWith(SCALE_REASONS[source])).toBe(true);
    expect(text.endsWith(`This trace is scaled ×1.50 ${HIGH_Q}`)).toBe(true);
  });

  it("says a beam-off factor matches only dark counts, and never calls it a high-q match", () => {
    const text = explainScale({ scale: 0.9, source: "beam_off" });
    expect(text).toBe(
      "One of these images was taken with the beam off, so it holds only the detector's dark counts and no scattering. The ×0.900 factor matches noise and means nothing.",
    );
    expect(text).not.toContain("high q");
  });

  it("names the background when it is the one scaled", () => {
    expect(explainScale({ scale: 1.5, source: "capillary" }, "The background")).toContain("The background is scaled ×1.50");
  });

  it("gives every high-q match its own reason", () => {
    expect(new Set(matched.map((s) => SCALE_REASONS[s])).size).toBe(matched.length);
  });
});
