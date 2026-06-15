import type { Meta, StoryObj } from "@storybook/react-vite";
import { GripHandle } from "./GripHandle";

const meta = {
  title: "ui/GripHandle",
  component: GripHandle,
} satisfies Meta<typeof GripHandle>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {};

export const InRow: Story = {
  render: () => (
    <div>
      <div className="group flex items-center gap-2 p-2 hover:bg-paper-sunk">
        <GripHandle />
        <span>Sample 037</span>
      </div>
      <p className="mt-2 text-ink-faint">
        Hover the row above to brighten the grip (the reveal is driven by the
        consumer row's `group`).
      </p>
    </div>
  ),
};
