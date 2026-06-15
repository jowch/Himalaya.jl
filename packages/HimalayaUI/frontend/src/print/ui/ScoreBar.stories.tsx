import type { Meta, StoryObj } from "@storybook/react-vite";
import { ScoreBar } from "./ScoreBar";

const meta = {
  title: "ui/ScoreBar",
  component: ScoreBar,
  args: { value: 0.86, phase: "Pn3m" },
  render: (args) => (
    <div style={{ width: 200 }}>
      <ScoreBar {...args} />
    </div>
  ),
} satisfies Meta<typeof ScoreBar>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {};

export const High: Story = { args: { value: 0.95, phase: "Im3m" } };
export const Mid: Story = { args: { value: 0.5, phase: "Lamellar" } };
export const Low: Story = { args: { value: 0.18, phase: "Ia3d" } };

export const Compact: Story = {
  args: { size: "compact", value: 0.72, phase: "Pn3m" },
  render: (args) => <ScoreBar {...args} />,
};
