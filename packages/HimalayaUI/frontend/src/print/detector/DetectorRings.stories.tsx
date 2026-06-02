import type { Meta, StoryObj } from "@storybook/react-vite";
import { useState } from "react";
import { DetectorRings } from "./DetectorRings";
import { DetectorImage } from "./DetectorImage";
import { buildRingPlacements } from "./detectorGeometry";
import { phaseColor } from "../../phases";
import thumb37 from "../fixtures/thumbs/37.png?url";

const meta: Meta<typeof DetectorRings> = { title: "detector/DetectorRings", component: DetectorRings };
export default meta;
type Story = StoryObj<typeof DetectorRings>;

// Fixture 37's real detector geometry: 981×1043 (portrait). The frame matches
// that aspect so the contained image fills it with no letterbox — which is what
// lets the rings (a unit-square SVG stretched to the frame) register on the image.
const IMG37 = { w: 981, h: 1043 };
const ASPECT37 = IMG37.w / IMG37.h; // width/height → ry correction so rings stay round
const frame = "relative aspect-[981/1043] w-[420px] overflow-hidden rounded border border-frame-edge bg-frame-edge";
// Real beam center for sample 37 (px → normalized). y is measured from the
// BOTTOM (detector origin), so flip it to the top-left screen origin the overlay
// uses — same convention as buildRingPlacements. Off-center: high and left.
const BEAM37 = { x: 421.409 / IMG37.w, y: 1 - 836.946 / IMG37.h }; // ≈ { 0.430, 0.198 }
// Ring colour comes from phaseColor(assignedPhase) — what the DetectorPanel
// composite will pass per ring. Here: a Pn3m series (amber) with one predicted-
// absent ghost, plus an unassigned leftover (no colour → neutral ink-faint).
const PN3M = phaseColor("Pn3m");
const QS = [
  { q: 0.082, color: PN3M }, { q: 0.116, color: PN3M },
  { q: 0.142, color: PN3M }, { q: 0.171, color: PN3M, ghost: true },
  { q: 0.205 },
];
// NOTE: ring RADII here are presentational (the null-calibration fallback). The
// real beam CENTER is known (above); the optical params that set true radii
// (sample-detector distance, energy, pixel size) come from the backend later —
// at which point buildRingPlacements gets a real DetectorCalibration instead.
const presentationalRings = buildRingPlacements(QS, null).rings;

// Centered presentational fallback (null calibration) — what the overlay shows
// before any beam-center data exists.
export const Vocabulary: Story = {
  render: () => (
    <div className={frame}><DetectorRings beamCenter={{ x: 0.5, y: 0.5 }} rings={presentationalRings} imageAspect={ASPECT37} /></div>
  ),
};

// Real off-center beam (sample 37) → concentric rings about (0.43, 0.80) spill
// past the frame as partial arcs, clipped by overflow:hidden.
export const OffCenterBeam: Story = {
  render: () => (
    <div className={frame}><DetectorRings beamCenter={BEAM37} rings={presentationalRings} imageAspect={ASPECT37} /></div>
  ),
};

// The real composition: the fixture image (filling the portrait frame) + rings
// about the real beam center + interactive hover.
export const OverImage: Story = {
  render: () => {
    const [hovered, setHovered] = useState<number | undefined>(undefined);
    return (
      <div className={frame}>
        <DetectorImage src={thumb37} size="full" className="h-full w-full" />
        <DetectorRings beamCenter={BEAM37} rings={presentationalRings} imageAspect={ASPECT37}
          {...(hovered !== undefined ? { hoveredQ: hovered } : {})}
          onHoverQ={setHovered} />
      </div>
    );
  },
};
