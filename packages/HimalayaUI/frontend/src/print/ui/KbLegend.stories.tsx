import type { Meta, StoryObj } from "@storybook/react-vite";
import { KbLegend } from "./KbLegend";

const meta = {
  title: "ui/KbLegend",
  component: KbLegend,
} satisfies Meta<typeof KbLegend>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {
  args: {
    shortcuts: [
      { keyLabel: "←/→", description: "navigate" },
      { keyLabel: "X", description: "reject" },
      { keyLabel: "Esc", description: "close" },
    ],
  },
};
