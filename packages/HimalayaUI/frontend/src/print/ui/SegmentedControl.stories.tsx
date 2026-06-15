import type { Meta, StoryObj } from "@storybook/react-vite";
import { SegmentedControl } from "./SegmentedControl";

const meta = {
  title: "ui/SegmentedControl",
  component: SegmentedControl,
  args: {
    "aria-label": "view",
    options: [
      { value: "comb", label: "Comb" },
      { value: "resid", label: "Resid" },
    ],
    value: "comb",
    onChange: () => {},
  },
} satisfies Meta<typeof SegmentedControl>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Small: Story = { args: { size: "sm" } };
export const Medium: Story = { args: { size: "md" } };
export const MiniXs: Story = { args: { size: "xs" } };
