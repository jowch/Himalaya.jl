import type { Meta, StoryObj } from "@storybook/react-vite";
import { Sparkline } from "./Sparkline";
import { realTraces } from "../fixtures/realTraces";

// Use a slice of trace 65 (Ia3d) — keep story data small
const trace65 = realTraces[65]!;
const sliceTrace = (
  full: { q: number[]; I: number[] },
  start: number,
  end: number,
) => ({
  q: full.q.slice(start, end),
  I: full.I.slice(start, end),
});

// ~100-point slice from the middle of the trace for a clean sparkline
const smallTrace = sliceTrace(trace65, 100, 200);

const meta = {
  title: "plot/Sparkline",
  component: Sparkline,
  parameters: { layout: "padded" },
} satisfies Meta<typeof Sparkline>;

export default meta;
type Story = StoryObj<typeof meta>;

/** Indexed trace with Ia3d violet hue. */
export const Indexed: Story = {
  args: {
    trace: smallTrace,
    phase: "Ia3d",
  },
};

/** Unindexed trace — ink-faint neutral. */
export const Unindexed: Story = {
  args: {
    trace: smallTrace,
    phase: null,
  },
};

/** Lamellar — indigo hue. */
export const Lamellar: Story = {
  args: {
    trace: smallTrace,
    phase: "Lamellar",
  },
};

/** Row of several sparklines — filmstrip feel. */
export const Row: Story = {
  args: { trace: smallTrace },
  render: () => (
    <div style={{ display: "flex", gap: 8, alignItems: "center" }}>
      <Sparkline trace={sliceTrace(trace65, 50, 130)} phase="Pn3m" />
      <Sparkline trace={sliceTrace(trace65, 100, 200)} phase="Ia3d" />
      <Sparkline trace={sliceTrace(trace65, 200, 300)} phase="Lamellar" />
      <Sparkline trace={sliceTrace(trace65, 300, 400)} phase="Hexagonal" />
      <Sparkline trace={sliceTrace(trace65, 150, 250)} phase={null} />
    </div>
  ),
};
