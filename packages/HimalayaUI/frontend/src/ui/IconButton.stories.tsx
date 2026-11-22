import type { Meta, StoryObj } from "@storybook/react-vite";
import { IconButton } from "./IconButton";

const meta = {
  title: "ui/IconButton",
  component: IconButton,
  args: { label: "Dismiss", dismiss: true },
} satisfies Meta<typeof IconButton>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {};

export const Ghost: Story = { args: { tone: "ghost" } };
export const Accent: Story = { args: { tone: "accent" } };
export const Danger: Story = { args: { tone: "danger" } };

export const WithGlyph: Story = {
  args: { label: "Next", dismiss: false, children: "›" },
};

export const Disabled: Story = { args: { disabled: true } };

export const HoverIntent: Story = {
  render: () => (
    <div className="flex flex-col gap-2">
      <div className="flex items-center gap-3">
        <IconButton label="Dismiss" dismiss tone="ghost" />
        <IconButton label="Action" tone="accent">›</IconButton>
        <IconButton label="Delete" tone="danger" dismiss />
      </div>
      <p className="text-xs text-ink-faint">
        danger is error-red at rest (the destructive baseline); ghost→ink and accent→terracotta appear on hover.
      </p>
    </div>
  ),
};
