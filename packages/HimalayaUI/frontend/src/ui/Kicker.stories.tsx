import type { Meta, StoryObj } from "@storybook/react-vite";
import { Kicker } from "./Kicker";

const meta = {
  title: "ui/Kicker",
  component: Kicker,
  args: { children: "BEAMTIME 2024-03" },
} satisfies Meta<typeof Kicker>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Faint: Story = { args: { tone: "faint" } };
export const Accent: Story = { args: { tone: "accent" } };
