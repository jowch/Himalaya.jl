import type { Meta, StoryObj } from "@storybook/react-vite";
import { Term } from "./Term";
import { explainScale } from "./scaleExplanation";

const meta = {
  title: "ui/Term",
  component: Term,
  args: { explain: explainScale({ scale: 1.1, source: "capillary" }), children: "×1.10", className: "text-caption text-ink-faint" },
} satisfies Meta<typeof Term>;

export default meta;
type Story = StoryObj<typeof meta>;

export const CompareScale: Story = {};

export const InProse: Story = {
  args: { explain: "Average Gaussian curvature of the minimal surface. Phases related by the Bonnet transform share it.", children: "κ", className: "text-body text-ink" },
};
