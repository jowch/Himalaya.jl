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
