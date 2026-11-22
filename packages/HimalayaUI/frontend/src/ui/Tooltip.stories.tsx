import type { Meta, StoryObj } from "@storybook/react-vite";
import { Tooltip } from "./Tooltip";
import { Button } from "./Button";

const meta = {
  title: "ui/Tooltip",
  component: Tooltip,
  args: { label: "Reanalyze stale indices" },
} satisfies Meta<typeof Tooltip>;

export default meta;
type Story = StoryObj<typeof meta>;

export const OnButton: Story = {
  args: {
    children: <Button variant="ghost">Reanalyze</Button>,
  },
};

export const Bottom: Story = {
  args: {
    side: "bottom",
    children: <Button variant="ghost">Reanalyze</Button>,
  },
};
