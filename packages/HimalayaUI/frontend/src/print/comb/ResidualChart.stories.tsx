import type { Meta, StoryObj } from "@storybook/react-vite";
import { ResidualChart } from "./ResidualChart";
import { PN3M, IM3M } from "./comb.fixtures";
import { phaseColor } from "../../phases";
import type { CombSeries } from "./combModel";

const meta: Meta<typeof ResidualChart> = { title: "comb/ResidualChart", component: ResidualChart };
export default meta;
type Story = StoryObj<typeof ResidualChart>;

const frame = "rounded border border-hair bg-plate p-3";

// A fit gone wrong: an in-range reflection, one residual well off the TOP (over-
// predicted, +24%) and one well off the BOTTOM (under-predicted, −17%) — each draws
// a sign-flipped off-scale chevron pinned at the track edge with its magnitude label,
// distinct from the marginal +4% case in OverflowSignal.
const WAY_OFF: CombSeries = {
  phase: "Fm3m", color: phaseColor("Fm3m"), rSquared: 0.612,
  teeth: [
    { q: 0.071, label: "√3", observed: true, residual: 0.006 },  // in range
    { q: 0.092, label: "√4", observed: true, residual: 0.24 },   // off the top
    { q: 0.118, label: "√8", observed: true, residual: -0.17 },  // off the bottom
  ],
};

export const Wide: Story = {
  render: () => <div className={frame} style={{ width: 760 }}><ResidualChart assigned={[PN3M, IM3M]} maxWidth={760} /></div>,
};

export const NarrowScrolling: Story = {
  render: () => <div className={frame} style={{ width: 360 }}><ResidualChart assigned={[PN3M, IM3M]} maxWidth={360} /></div>,
};

export const OverflowSignal: Story = {
  render: () => <div className={frame} style={{ width: 760 }}><ResidualChart assigned={[IM3M]} maxWidth={760} /></div>,
};

// Residuals that run well beyond the ±3% track: up chevron + "+24%" off the top,
// down chevron + "−17%" off the bottom.
export const WayOffScale: Story = {
  render: () => <div className={frame} style={{ width: 760 }}><ResidualChart assigned={[WAY_OFF]} maxWidth={760} /></div>,
};
