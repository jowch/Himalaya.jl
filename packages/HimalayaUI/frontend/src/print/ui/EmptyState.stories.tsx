import type { Meta, StoryObj } from "@storybook/react-vite";
import { EmptyState } from "./EmptyState";

const meta = {
  title: "ui/EmptyState",
  component: EmptyState,
  args: { title: "No series match" },
} satisfies Meta<typeof EmptyState>;

export default meta;
type Story = StoryObj<typeof meta>;

export const FilterNoMatch: Story = {
  args: {
    title: "No series match",
    body: "Try clearing a filter, or widen the beamtime range to see more.",
  },
};

export const TitleOnly: Story = { args: { title: "No samples yet" } };
