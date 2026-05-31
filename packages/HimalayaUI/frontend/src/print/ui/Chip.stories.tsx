import type { Meta, StoryObj } from "@storybook/react-vite";
import { Chip } from "./Chip";

const meta = {
  title: "ui/Chip",
  component: Chip,
  args: { children: "Pn3m" },
} satisfies Meta<typeof Chip>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Static: Story = { args: { variant: "static" } };
export const Removable: Story = {
  args: { variant: "removable", onRemove: () => {} },
};
export const Add: Story = { args: { variant: "add", children: "+ add" } };
export const ToggleOff: Story = {
  args: { variant: "toggle", active: false, children: "Kept" },
};
export const ToggleOn: Story = {
  args: { variant: "toggle", active: true, children: "Kept" },
};
export const Trigger: Story = {
  args: { variant: "trigger", children: "Beamtime" },
};
