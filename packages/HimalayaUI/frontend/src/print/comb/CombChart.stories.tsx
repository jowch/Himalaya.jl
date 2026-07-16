import type { Meta, StoryObj } from "@storybook/react-vite";
import { useState } from "react";
import { CombChart } from "./CombChart";
import { ResidualChart } from "./ResidualChart";
import { PN3M, IM3M, LEFTOVER } from "./comb.fixtures";

const meta: Meta<typeof CombChart> = { title: "comb/CombChart", component: CombChart };
export default meta;
type Story = StoryObj<typeof CombChart>;

const frame = "rounded border border-hair bg-plate p-3";

export const Wide: Story = {
  render: () => <div className={frame} style={{ width: 760 }}><CombChart assigned={[PN3M, IM3M]} leftover={LEFTOVER} maxWidth={760} /></div>,
};

export const NarrowScrolling: Story = {
  render: () => <div className={frame} style={{ width: 360 }}><CombChart assigned={[PN3M, IM3M]} leftover={LEFTOVER} maxWidth={360} /></div>,
};

export const WithPreviewRow: Story = {
  render: () => <div className={frame} style={{ width: 760 }}><CombChart assigned={[PN3M]} hovered={IM3M} leftover={LEFTOVER} maxWidth={760} /></div>,
};

export const LinkedWithResidual: Story = {
  render: () => {
    const [hoveredQ, setHoveredQ] = useState<number | undefined>(undefined);
    const linkProps = { ...(hoveredQ !== undefined ? { hoveredQ } : {}), onHoverQ: setHoveredQ };
    return (
      <div className="flex flex-col gap-3" style={{ width: 760 }}>
        <div className={frame}><CombChart assigned={[PN3M, IM3M]} leftover={LEFTOVER} maxWidth={760} {...linkProps} /></div>
        <div className={frame}><ResidualChart assigned={[PN3M, IM3M]} maxWidth={760} {...linkProps} /></div>
      </div>
    );
  },
};
