import type { CompareScale } from "../workspace/compare";

type MatchedSource = Exclude<CompareScale["source"], "normalised" | "monitor" | "beam_off">;

export const SCALE_REASONS: Record<MatchedSource, string> = {
  capillary: "Scattered intensity grows with how much material the beam crosses: the path length and how evenly the sample fills it. Both change from capillary to capillary. Dividing by I₂ corrects for the beam and for absorption, not for how much material the beam crossed.",
  other_setup: "Detector distance, the beamstop diode and the I₂ base the traces are divided by can differ between setups, so counts per pixel and I₂ may not compare.",
  i2_collapsed: "The beam reaching the beamstop diode fell far below its usual value on one of these images, as when the beam misses the diode, so I₂ no longer measures the beam that got through the sample.",
  unknown_mode: "The amount of material in the beam only matches in an autosampler, where every sample fills the same cell, and this setup isn't marked as one.",
  unknown_normalisation: "The integration settings weren't found, so it isn't known whether these traces were divided by I₂.",
  matched: "Without an I₂ reading there is no transmitted beam to scale by.",
};

export const formatScale = (scale: number): string => scale.toPrecision(3);

export function explainScale({ scale, source }: CompareScale, subject = "This trace"): string {
  const factor = `×${formatScale(scale)}`;
  if (source === "normalised") {
    return "Both traces are divided by I₂, the beam transmitted through the sample. Every sample fills the same autosampler cell, so the beam crosses the same thickness of sample and the traces are already on one scale.";
  }
  if (source === "monitor") {
    return `Scattered intensity grows with the beam that gets through the sample. Every sample fills the same autosampler cell, so the thickness matches and the beam alone sets the scale. ${subject} is scaled ${factor}, the ratio of the beam that got through each sample while it was measured.`;
  }
  if (source === "beam_off") {
    return `One of these images was taken with the beam off, so it holds only the detector's dark counts and no scattering. The ${factor} factor matches noise and means nothing.`;
  }
  return `${SCALE_REASONS[source]} ${subject} is scaled ${factor} to match at high q, where both traces should be mostly solvent and capillary scattering. A sample that still scatters there, from a bilayer's form factor or just more lipid, throws this factor off.`;
}
