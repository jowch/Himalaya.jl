import type { Meta, StoryObj } from "@storybook/react-vite";
import { PhaseBlock } from "./PhaseBlock";

// Storybook action: log the remove click without pulling in @storybook/test.
const noop = (): void => {};

const meta = {
  title: "components/PhaseBlock",
  component: PhaseBlock,
  args: {
    phase: "Pn3m",
    score: 0.92,
    meta: "a = 14.2 nm · 5 reflections",
    onRemove: noop,
  },
} satisfies Meta<typeof PhaseBlock>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {};

// Three assigned phases at high / mid / low scores to show the phase colors
// and bar fills in the assignment cart.
export const Variants: Story = {
  render: () => (
    <div className="w-80">
      <PhaseBlock
        phase="Pn3m"
        score={0.94}
        meta="a = 14.2 nm · 5 reflections"
        onRemove={noop}
      />
      <PhaseBlock
        phase="Im3m"
        score={0.71}
        meta="a = 9.1 nm · 4 reflections"
        onRemove={noop}
      />
      <PhaseBlock
        phase="Lamellar"
        score={0.42}
        meta="d = 6.0 nm · 2 reflections"
        onRemove={noop}
      />
    </div>
  ),
};
