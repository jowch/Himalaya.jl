import type { Meta, StoryObj } from "@storybook/react-vite";
import { PhaseChip } from "./PhaseChip";

const meta = {
  title: "ui/PhaseChip",
  component: PhaseChip,
  args: { phase: "Pn3m", variant: "tint", size: "sm" },
} satisfies Meta<typeof PhaseChip>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = { args: { phase: "Pn3m" } };
export const Solid: Story = { args: { phase: "Pn3m", variant: "solid" } };
export const Coexistence: Story = { args: { phase: "Pn3m", coexistWith: ["Lamellar"] } };
export const ThreePhase: Story = { args: { phase: "Pn3m", coexistWith: ["Lamellar", "Hexagonal"] } };
