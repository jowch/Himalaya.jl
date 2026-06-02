import type { Meta, StoryObj } from "@storybook/react-vite";
import { useState } from "react";
import { DetectorRings } from "./DetectorRings";
import { DetectorImage } from "./DetectorImage";
import { buildRingPlacements } from "./detectorGeometry";
import thumb37 from "../fixtures/thumbs/37.png?url";

const meta: Meta<typeof DetectorRings> = { title: "detector/DetectorRings", component: DetectorRings };
export default meta;
type Story = StoryObj<typeof DetectorRings>;

const frame = "relative aspect-square w-[420px] overflow-hidden rounded border border-frame-edge bg-frame-edge";
const QS = [
  { q: 0.082, color: "var(--color-success)" }, { q: 0.116, color: "var(--color-success)" },
  { q: 0.142, color: "var(--color-success)" }, { q: 0.171, color: "var(--color-success)", ghost: true },
  { q: 0.205 },
];

// Centered presentational fallback (null calibration) over the dark window.
export const Vocabulary: Story = {
  render: () => {
    const { beamCenter, rings } = buildRingPlacements(QS, null);
    return <div className={frame}><DetectorRings beamCenter={beamCenter} rings={rings} /></div>;
  },
};

// Off-center beam -> partial arcs clipped by the frame. Demonstrates the geometry
// path with a hand-set calibration (beam near the bottom-right).
export const OffCenterBeam: Story = {
  render: () => {
    const { beamCenter, rings } = buildRingPlacements(QS, {
      beamCenterPx: { x: 880, y: 980 }, imageSizePx: { w: 1000, h: 1000 },
      sampleDistanceMm: 800, pixelSizeMm: 0.172, energyKeV: 12.4,
    });
    return <div className={frame}><DetectorRings beamCenter={beamCenter} rings={rings} /></div>;
  },
};

// The real composition: image + rings + interactive hover.
export const OverImage: Story = {
  render: () => {
    const [hovered, setHovered] = useState<number | undefined>(undefined);
    const { beamCenter, rings } = buildRingPlacements(QS, null);
    return (
      <div className={frame}>
        <DetectorImage src={thumb37} size="full" className="h-full w-full" />
        <DetectorRings beamCenter={beamCenter} rings={rings}
          {...(hovered !== undefined ? { hoveredQ: hovered } : {})}
          onHoverQ={setHovered} />
      </div>
    );
  },
};
