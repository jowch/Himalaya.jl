import type { Meta, StoryObj } from "@storybook/react-vite";
import { useState } from "react";
import { WaterfallChart } from "./WaterfallChart";
import { FULL, TRANSITION, MIXED_STATES } from "./waterfall.fixtures";

const meta: Meta<typeof WaterfallChart> = {
  title: "waterfall/WaterfallChart",
  component: WaterfallChart,
  parameters: { layout: "padded" },
};
export default meta;
type Story = StoryObj<typeof WaterfallChart>;

const frame = "rounded border border-hair bg-plate p-4";

// Wide: the full real Sample-9 set at the plate max width.
export const Wide: Story = {
  render: () => (
    <div className={frame} style={{ width: 980 }}>
      <WaterfallChart rows={FULL} maxWidth={920} />
    </div>
  ),
};

// Narrow: fit-to-width down to the floor (no horizontal scroll).
export const Narrow: Story = {
  render: () => (
    <div className={frame} style={{ width: 600 }}>
      <WaterfallChart rows={TRANSITION} maxWidth={560} minWidth={560} />
    </div>
  ),
};

// Linear q.
export const LinearQ: Story = {
  render: () => (
    <div className={frame} style={{ width: 980 }}>
      <WaterfallChart rows={TRANSITION} maxWidth={920} xType="linear" />
    </div>
  ),
};

// Mixed states: indexed + form-factor (no anchors, neutral line) + unindexed.
export const MixedStates: Story = {
  render: () => (
    <div className={frame} style={{ width: 980 }}>
      <WaterfallChart rows={MIXED_STATES} maxWidth={920} />
    </div>
  ),
};

// Controlled hover: a parent drives the hot row (e.g. from a member list).
export const ControlledHover: Story = {
  render: () => {
    // eslint-disable-next-line react-hooks/rules-of-hooks
    const [hot, setHot] = useState<string | undefined>(TRANSITION[0]?.key);
    const hoverProps = { ...(hot !== undefined ? { hoveredKey: hot } : {}), onHoverRow: setHot };
    return (
      <div className={frame} style={{ width: 980 }}>
        <WaterfallChart rows={TRANSITION} maxWidth={920} {...hoverProps} />
      </div>
    );
  },
};
