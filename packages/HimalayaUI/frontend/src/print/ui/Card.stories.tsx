import type { Meta, StoryObj } from "@storybook/react-vite";
import { Card } from "./Card";

const meta = {
  title: "ui/Card",
  component: Card,
  args: { children: "Pn3m · 4 peaks" },
} satisfies Meta<typeof Card>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {};
export const Elevated: Story = { args: { elevated: true } };
export const Draft: Story = { args: { draft: true, children: "New recipe" } };
export const Padded: Story = { args: { elevated: true, padding: "md", children: "Padded plate body" } };
export const Selected: Story = { args: { selected: true, padding: "sm", children: "Pn3m · selected" } };
