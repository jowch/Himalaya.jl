import type { CompareScale } from "../workspace/compare";

type MatchedSource = Exclude<CompareScale["source"], "normalised" | "monitor">;

export const SCALE_REASONS: Record<MatchedSource, string> = {
  capillary: "Scattered intensity grows with how much material the beam crosses: the path length and how evenly the sample fills it. Both change from capillary to capillary, and I₂ corrects only for transmission.",
  other_setup: "Detector distance and the beamstop diode differ between setups, so counts per pixel and I₂ don't compare.",
  i2_collapsed: "I₂ fell far below its usual value on one of these images, so a trace divided by it is inflated.",
  beam_off: "One of these images was taken with the beam off, so it has no transmitted beam to scale by.",
  unknown_mode: "The amount of material in the beam only matches in an autosampler, where every sample fills the same cell, and this setup isn't marked as one.",
  unknown_normalisation: "The integration settings weren't found, so it isn't known whether these traces were divided by I₂.",
  matched: "Without an I₂ reading there is no transmitted beam to scale by.",
};

export const formatScale = (scale: number): string => scale.toPrecision(3);

export function explainScale({ scale, source }: CompareScale, subject = "This trace"): string {
  if (source === "normalised") {
    return "Both traces are divided by I₂, the beam transmitted through the sample. Every sample fills the same autosampler cell, so the beam crosses the same thickness of sample and the traces are already on one scale.";
  }
  if (source === "monitor") {
    return `Scattered intensity grows with the beam that gets through the sample. ${subject} is scaled ×${formatScale(scale)}, the ratio of that transmitted beam (I₂ × exposure time) between the two images.`;
  }
  return `${SCALE_REASONS[source]} ${subject} is scaled ×${formatScale(scale)} to match at high q, where the solvent dominates. A sample that still scatters there, such as vesicles, throws this factor off.`;
}
