import type { Meta, StoryObj } from "@storybook/react-vite";
import { KbKey } from "./KbKey";

const meta = {
  title: "ui/KbKey",
  component: KbKey,
  args: { children: "⌘K" },
} satisfies Meta<typeof KbKey>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = { args: { children: "⌘K" } };
export const Letter: Story = { args: { children: "X" } };
