import type { Meta, StoryObj } from "@storybook/react-vite";
import { CardFigure } from "./CardFigure";
import { FULL, TRANSITION, MIXED_STATES } from "./waterfall.fixtures";

const meta: Meta<typeof CardFigure> = {
  title: "waterfall/CardFigure",
  component: CardFigure,
  parameters: { layout: "padded" },
};
export default meta;
type Story = StoryObj<typeof CardFigure>;

const frame = "rounded border border-hair bg-plate p-4";

// The canonical three-member transition hero, frozen and mini.
export const Transition: Story = {
  render: () => (
    <div className={frame} style={{ width: 360 }}>
      <CardFigure rows={TRANSITION} />
    </div>
  ),
};

// The full real Sample-9 set — taller stack, same frozen treatment.
export const Full: Story = {
  render: () => (
    <div className={frame} style={{ width: 360 }}>
      <CardFigure rows={FULL} />
    </div>
  ),
};

// Mixed states: indexed + form-factor (no anchors, neutral line) + unindexed.
export const Mixed: Story = {
  render: () => (
    <div className={frame} style={{ width: 360 }}>
      <CardFigure rows={MIXED_STATES} />
    </div>
  ),
};
