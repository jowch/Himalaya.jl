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

/** The representative-marker form: 9px accent circle with a 1.5px plate ring,
 *  shown here on an accent backdrop so the separator ring reads. */
export const BorderedMarker: Story = {
  args: { tone: "accent", size: "md", bordered: true },
  render: (args) => (
    <span className="inline-flex bg-print-accent p-2 rounded-sm">
      <Dot {...args} aria-hidden="true" />
    </span>
  ),
};

export const WithCaption: Story = {
  // `label` is the screen-reader name; the visible word is the sibling text.
  render: () => (
    <span className="inline-flex items-center gap-1.5">
      <Dot tone="success" label="Kept" />
      <span className="text-sm text-ink-soft">Kept</span>
    </span>
  ),
};
