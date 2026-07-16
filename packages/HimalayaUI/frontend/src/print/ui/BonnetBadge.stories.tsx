import type { Meta, StoryObj } from "@storybook/react-vite";
import { BonnetBadge } from "./BonnetBadge";

const meta = {
  title: "ui/BonnetBadge",
  component: BonnetBadge,
} satisfies Meta<typeof BonnetBadge>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {};

export const InCandidateRow: Story = {
  render: () => (
    <span className="font-mono text-sm">
      q = 0.142 <BonnetBadge />
    </span>
  ),
};
