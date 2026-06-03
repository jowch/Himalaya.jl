// src/print/components/DetectorPanel.stories.tsx
import { useState } from "react";
import type { Meta, StoryObj } from "@storybook/react-vite";
import { DetectorPanel } from "./DetectorPanel";
import { ThumbnailGallery, type GalleryExposure } from "./ThumbnailGallery";
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

const meta: Meta<typeof DetectorPanel> = {
  title: "components/DetectorPanel",
  component: DetectorPanel,
  parameters: { layout: "padded" },
};
export default meta;
type Story = StoryObj<typeof meta>;

function DetectorDemo() {
  const [cur, setCur] = useState(0);
  return (
    <div style={{ width: 372 }}>
      <DetectorPanel
        src={SRCS[cur]!}
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

export const Empty: Story = {
  args: { src: null, hint: "No exposure selected." },
};
