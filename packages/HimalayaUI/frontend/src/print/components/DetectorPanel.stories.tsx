// src/print/components/DetectorPanel.stories.tsx
import { useState } from "react";
import type { Meta, StoryObj } from "@storybook/react-vite";
import { DetectorPanel } from "./DetectorPanel";
import { ThumbnailGallery, type GalleryExposure } from "./ThumbnailGallery";
import type { RingInput } from "../detector/detectorGeometry";
import { phaseColor } from "../../phases";
import thumb37 from "../fixtures/thumbs/37.png?url";
import thumb65 from "../fixtures/thumbs/65.png?url";
import thumb66 from "../fixtures/thumbs/66.png?url";
import thumb67 from "../fixtures/thumbs/67.png?url";

const SRCS = [thumb65, thumb66, thumb67, thumb37];
const EXPOSURES: GalleryExposure[] = SRCS.map((src, i) => ({
  id: i,
  src,
  frameNo: 65 + i,
  representative: i === 0,
}));

// Assigned reflections, phase-coloured: a dominant Pn3m + a subtle Im3m (one
// predicted-absent ghost). q in Å⁻¹; null calibration → centered fallback.
const RINGS: RingInput[] = [
  { q: 0.045, color: phaseColor("Pn3m") },
  { q: 0.064, color: phaseColor("Pn3m") },
  { q: 0.078, color: phaseColor("Pn3m") },
  { q: 0.052, color: phaseColor("Im3m") },
  { q: 0.09, color: phaseColor("Im3m"), ghost: true },
];

const meta: Meta<typeof DetectorPanel> = {
  title: "components/DetectorPanel",
  component: DetectorPanel,
  parameters: { layout: "padded" },
};
export default meta;
type Story = StoryObj<typeof meta>;

function DetectorDemo() {
  const [cur, setCur] = useState(0);
  // hoveredQ ↔ onHoverQ: hovering a ring lights it; in the full Focus page this
  // q would be shared with the trace + comb so all three light together.
  const [hoveredQ, setHoveredQ] = useState<number | undefined>(undefined);
  return (
    <div style={{ width: 372 }}>
      <DetectorPanel
        src={SRCS[cur]!}
        rings={RINGS}
        {...(hoveredQ !== undefined ? { hoveredQ } : {})}
        onHoverQ={setHoveredQ}
        hint="The real source. Rings in phase colour; hover a peak, ring, or comb tooth — the triple lights up."
        tools={
          <ThumbnailGallery
            exposures={EXPOSURES}
            selectedId={cur}
            onSelect={setCur}
            size="xs"
          />
        }
      />
    </div>
  );
}

export const Default: Story = { render: () => <DetectorDemo /> };

/** Many exposures — the header strip scrolls horizontally in place. */
function OverflowDemo() {
  const [cur, setCur] = useState(0);
  const many: GalleryExposure[] = Array.from({ length: 14 }, (_, i) => ({
    id: i,
    src: SRCS[i % SRCS.length]!,
    frameNo: 65 + i,
    representative: i === 0,
  }));
  return (
    <div style={{ width: 372 }}>
      <DetectorPanel
        src={SRCS[cur % SRCS.length]!}
        rings={RINGS}
        tools={
          <ThumbnailGallery exposures={many} selectedId={cur} onSelect={setCur} size="xs" />
        }
      />
    </div>
  );
}

export const ManyExposures: Story = { render: () => <OverflowDemo /> };

export const Empty: Story = {
  args: { src: null, hint: "No exposure selected." },
};
