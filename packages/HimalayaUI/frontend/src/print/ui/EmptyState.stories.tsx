import type { Meta, StoryObj } from "@storybook/react-vite";
import { EmptyState } from "./EmptyState";
import { Button } from "./Button";

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

// Error states stop dead-ending: the body states what happened, the action
// slot offers the way forward (the consumer composes a ui Button; the
// primitive owns the gap below the body).
export const WithAction: Story = {
  args: {
    title: "Could not load the corpus",
    body: "The sample list request failed.",
    action: <Button variant="outline">Reload the corpus</Button>,
  },
};
