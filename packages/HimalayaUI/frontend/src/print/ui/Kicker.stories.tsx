import type { Meta, StoryObj } from "@storybook/react-vite";
import { Kicker } from "./Kicker";

const meta = {
  title: "ui/Kicker",
  component: Kicker,
  args: { children: "Phase call" },
} satisfies Meta<typeof Kicker>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {};

export const Faint: Story = { args: { tone: "faint" } };
export const Accent: Story = { args: { tone: "accent" } };
