import { useState } from "react";
import type { Meta, StoryObj } from "@storybook/react-vite";
import { CombsPanel, type CombView } from "./CombsPanel";
import { PN3M, IM3M, LEFTOVER } from "../comb/comb.fixtures";

const meta: Meta<typeof CombsPanel> = {
  title: "components/CombsPanel",
  component: CombsPanel,
  parameters: { layout: "padded" },
};
export default meta;
type Story = StoryObj<typeof meta>;

function CombsDemo() {
  const [view, setView] = useState<CombView>("comb");
  const [hoveredQ, setHoveredQ] = useState<number | undefined>(undefined);
  const hover = hoveredQ !== undefined ? { hoveredQ } : {};
  return (
    <div style={{ width: 720, height: 416 }}>
      <CombsPanel
        assigned={[PN3M, IM3M]}
        leftover={LEFTOVER}
        view={view}
        onViewChange={setView}
        onHoverQ={setHoveredQ}
        className="h-full"
        {...hover}
      />
    </div>
  );
}

export const Default: Story = { render: () => <CombsDemo /> };
