import type { Meta, StoryObj } from "@storybook/react-vite";
import { ResidualChart } from "./ResidualChart";
import { ResidualLegend } from "../components/ResidualLegend";
import { PN3M, IM3M } from "./comb.fixtures";
import { phaseColor } from "../../phases";
import type { CombSeries } from "./combModel";

const meta: Meta<typeof ResidualChart> = { title: "comb/ResidualChart", component: ResidualChart };
export default meta;
type Story = StoryObj<typeof ResidualChart>;

const frame = "rounded border border-hair bg-plate p-3";

// All three tiers in one row: in tolerance (filled dot), out of tolerance but on the
// track (hollow dot at its real value), and off-scale beyond the track (a chevron
// clamped at the edge, sign by direction). No magnitude on the chevron — off-scale is
// a binary "ran off the track".
const TIERS: CombSeries = {
  phase: "Fm3m", color: phaseColor("Fm3m"), rSquared: 0.71,
  teeth: [
    { q: 0.071, label: "√3", observed: true, residual: 0.006 },   // in tolerance → filled dot
    { q: 0.083, label: "√4", observed: true, residual: 0.026 },   // out of tolerance, on-scale → hollow dot (high)
    { q: 0.100, label: "√6", observed: true, residual: -0.028 },  // out of tolerance, on-scale → hollow dot (low)
    { q: 0.118, label: "√8", observed: true, residual: 0.045 },   // off-scale → up chevron at top edge
    { q: 0.140, label: "√10", observed: true, residual: -0.05 },  // off-scale → down chevron at bottom edge
  ],
};

export const Wide: Story = {
  render: () => <div className={frame} style={{ width: 760 }}><ResidualChart assigned={[PN3M, IM3M]} maxWidth={760} /></div>,
};

export const NarrowScrolling: Story = {
  render: () => <div className={frame} style={{ width: 360 }}><ResidualChart assigned={[PN3M, IM3M]} maxWidth={360} /></div>,
};

// The three tiers: in-tolerance filled dot, out-of-tolerance hollow dot (at value),
// and off-scale chevron clamped at the edge (both directions).
export const ToleranceTiers: Story = {
  render: () => (
    // Legend shown as in CombsPanel: the residual vocabulary + the band/track note.
    <div className={frame} style={{ width: 760 }}>
      <ResidualChart assigned={[TIERS]} maxWidth={760} />
      <ResidualLegend className="mt-2.5" />
    </div>
  ),
};
