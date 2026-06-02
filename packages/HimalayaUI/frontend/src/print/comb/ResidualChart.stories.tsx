import type { Meta, StoryObj } from "@storybook/react-vite";
import { ResidualChart } from "./ResidualChart";
import { PN3M, IM3M } from "./comb.fixtures";
import { phaseColor } from "../../phases";
import type { CombSeries } from "./combModel";

const meta: Meta<typeof ResidualChart> = { title: "comb/ResidualChart", component: ResidualChart };
export default meta;
type Story = StoryObj<typeof ResidualChart>;

const frame = "rounded border border-hair bg-plate p-3";

// One reflection in range, one just past the top of the track, one just past the
// bottom — each off-scale point draws a clamped hollow circle (sign by which edge)
// with a small directional chevron. No magnitude: off-scale is a binary "this is off".
const OFF_SCALE: CombSeries = {
  phase: "Fm3m", color: phaseColor("Fm3m"), rSquared: 0.71,
  teeth: [
    { q: 0.071, label: "√3", observed: true, residual: 0.006 },   // in range (filled dot)
    { q: 0.092, label: "√4", observed: true, residual: 0.038 },   // off the top → hollow + up chevron
    { q: 0.118, label: "√8", observed: true, residual: -0.034 },  // off the bottom → hollow + down chevron
  ],
};

export const Wide: Story = {
  render: () => <div className={frame} style={{ width: 760 }}><ResidualChart assigned={[PN3M, IM3M]} maxWidth={760} /></div>,
};

export const NarrowScrolling: Story = {
  render: () => <div className={frame} style={{ width: 360 }}><ResidualChart assigned={[PN3M, IM3M]} maxWidth={360} /></div>,
};

// Off-scale residuals (just past ±3%): a hollow clamped circle + small chevron at the
// top edge (over-predicted) and the bottom edge (under-predicted).
export const OffScale: Story = {
  render: () => <div className={frame} style={{ width: 760 }}><ResidualChart assigned={[OFF_SCALE]} maxWidth={760} /></div>,
};
