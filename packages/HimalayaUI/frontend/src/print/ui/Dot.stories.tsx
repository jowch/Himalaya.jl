import type { Meta, StoryObj } from "@storybook/react-vite";
import { Dot } from "./Dot";

const meta = {
  title: "ui/Dot",
  component: Dot,
  args: { size: "sm" },
} satisfies Meta<typeof Dot>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {};

export const Accent: Story = { args: { tone: "accent" } };
export const Success: Story = { args: { tone: "success" } };
export const Muted: Story = { args: { tone: "muted" } };
export const Neutral: Story = { args: { tone: "neutral" } };

// label promotes the dot to role=img with an accessible name.
export const Labeled: Story = { args: { tone: "success", label: "Indexed" } };
