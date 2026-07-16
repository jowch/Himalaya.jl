import { useState } from "react";
import type { Meta, StoryObj } from "@storybook/react-vite";
import { SeriesPlate, type SeriesScale } from "./SeriesPlate";
import { TRANSITION } from "../waterfall/waterfall.fixtures";

const meta: Meta<typeof SeriesPlate> = {
  title: "components/SeriesPlate",
  component: SeriesPlate,
  parameters: { layout: "padded" },
};
export default meta;
type Story = StoryObj<typeof meta>;

function OpenDemo() {
  const [scale, setScale] = useState<SeriesScale>("log");
  return (
    <div style={{ maxWidth: 1080 }}>
      <SeriesPlate
        kicker="Series"
        title="LL37 titration into Lipid 1-2"
        subtitle="7 exposures · variable: LL37 : lipid · SSRL Apr 2026"
        rows={TRANSITION}
        scale={scale}
        onScaleChange={setScale}
        legendPhases={["Pn3m", "Lamellar"]}
        footHint="peaks are light anchors — hover a trace to read its indexing"
        footNote="offset ×1.0 · log I"
      />
    </div>
  );
}

export const Hero: Story = { render: () => <OpenDemo /> };
