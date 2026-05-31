import type { Meta, StoryObj } from "@storybook/react-vite";
import { SignalBars } from "./SignalBars";

const meta = {
  title: "ui/SignalBars",
  component: SignalBars,
  args: { value: 4, max: 5 },
} satisfies Meta<typeof SignalBars>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Strong: Story = { args: { value: 4, max: 5 } };
export const Weak: Story = { args: { value: 1, max: 5 } };
export const Full: Story = { args: { value: 5, max: 5 } };
