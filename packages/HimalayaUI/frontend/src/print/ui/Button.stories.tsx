import type { Meta, StoryObj } from "@storybook/react-vite";
import { Button } from "./Button";

const meta = {
  title: "ui/Button",
  component: Button,
  args: { children: "Curate" },
} satisfies Meta<typeof Button>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Ghost: Story = { args: { variant: "ghost" } };
export const Solid: Story = { args: { variant: "solid" } };
export const Accent: Story = { args: { variant: "accent" } };
export const Danger: Story = { args: { variant: "danger", children: "Delete sample" } };
export const Armed: Story = { args: { variant: "ghost", armed: true, children: "+ Peak" } };
