import type { Meta, StoryObj } from "@storybook/react-vite";
import { PhaseStrip } from "./PhaseStrip";

const meta = {
  title: "ui/PhaseStrip",
  component: PhaseStrip,
  args: {
    segments: [
      { phase: "Pn3m" },
      { phase: "Pn3m" },
      { phase: "Im3m" },
      { phase: null },
    ],
  },
} satisfies Meta<typeof PhaseStrip>;

export default meta;
type Story = StoryObj<typeof meta>;

// A → B transition (two distinct indexed phases).
export const Default: Story = {};

// Single phase throughout.
export const Throughout: Story = {
  args: {
    segments: [{ phase: "Lamellar" }, { phase: "Lamellar" }, { phase: "Lamellar" }],
  },
};

// Coexistence gradient in one cell.
export const Coexistence: Story = {
  args: {
    segments: [
      { phase: "Pn3m" },
      { phase: "Pn3m", coexistWith: "Im3m" },
      { phase: "Im3m" },
    ],
  },
};

// No indexed segment → empty caption.
export const Unindexed: Story = {
  args: { segments: [{ phase: null }, { phase: null }] },
};

export const Small: Story = { args: { size: "sm" } };

export const Vertical: Story = {
  args: { orientation: "vertical" },
  render: (args) => (
    <div style={{ height: 120 }}>
      <PhaseStrip {...args} />
    </div>
  ),
};
