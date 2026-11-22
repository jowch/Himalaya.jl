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
export const Success: Story = { args: { variant: "success", children: "Keep" } };
export const Danger: Story = { args: { variant: "danger", children: "Delete sample" } };
export const Outline: Story = { args: { variant: "outline", children: "Drop" } };
export const Armed: Story = { args: { variant: "ghost", armed: true, children: "+ Peak" } };
export const Disabled: Story = {
  args: { variant: "solid", children: "Confirm & build →" },
  render: (args) => (
    <div className="inline-flex gap-3">
      <Button {...args} disabled />
      <Button {...args} />
    </div>
  ),
};
export const GhostInverse: Story = {
  args: { variant: "ghostInverse", children: "Restore" },
  render: (args) => (
    <div className="bg-ink p-4 rounded-md inline-flex">
      <Button {...args} />
    </div>
  ),
};
