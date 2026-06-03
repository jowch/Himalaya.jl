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

// Real sample-37 detector geometry — the SAME values the DetectorRings
// storyboard uses, so the panel registers rings on the actual (off-center) beam
// rather than the centered fallback. 981×1043 portrait; beam high + left.
const IMG37 = { w: 981, h: 1043 };
const ASPECT37 = IMG37.w / IMG37.h; // ≈ 0.94 (portrait)
const BEAM37 = { x: 421.409 / IMG37.w, y: 1 - 836.946 / IMG37.h }; // ≈ {0.430, 0.198}

// Assigned reflections, phase-coloured: a dominant Pn3m + a subtle Im3m (one
// predicted-absent ghost). q in Å⁻¹; radii are the presentational fallback
// (real radii arrive with full optical calibration).
const RINGS: RingInput[] = [
  { q: 0.082, color: phaseColor("Pn3m") },
  { q: 0.116, color: phaseColor("Pn3m") },
  { q: 0.142, color: phaseColor("Pn3m") },
  { q: 0.103, color: phaseColor("Im3m") },
  { q: 0.171, color: phaseColor("Im3m"), ghost: true },
];

const meta: Meta<typeof DetectorPanel> = {
  title: "components/DetectorPanel",
  component: DetectorPanel,
  parameters: { layout: "padded" },
};
export default meta;
type Story = StoryObj<typeof meta>;

function DetectorDemo() {
  // Default to sample 37 (index 3) so the displayed image pairs exactly with the
  // real beam center below; selecting another exposure keeps the same demo beam.
  const [cur, setCur] = useState(3);
  // hoveredQ ↔ onHoverQ: hovering a ring lights it. The cross-panel trace↔ring↔
  // comb link is wired at PAGE assembly (the page owns the shared hovered-q).
  const [hoveredQ, setHoveredQ] = useState<number | undefined>(undefined);
  return (
    <div style={{ width: 372 }}>
      <DetectorPanel
        src={SRCS[cur]!}
        rings={RINGS}
        beamCenter={BEAM37}
        imageAspect={ASPECT37}
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
        beamCenter={BEAM37}
        imageAspect={ASPECT37}
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
