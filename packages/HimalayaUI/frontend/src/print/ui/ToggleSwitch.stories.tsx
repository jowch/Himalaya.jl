import type { Meta, StoryObj } from "@storybook/react-vite";
import { ToggleSwitch } from "./ToggleSwitch";

const meta = {
  title: "ui/ToggleSwitch",
  component: ToggleSwitch,
  args: { checked: false, label: "Heatmap", onChange: () => {} },
} satisfies Meta<typeof ToggleSwitch>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Off: Story = { args: { checked: false, label: "Heatmap", onChange: () => {} } };
export const On: Story = { args: { checked: true, label: "Heatmap", onChange: () => {} } };
export const NoVisibleLabel: Story = {
  args: { checked: true, label: "Heatmap", hideLabel: true, onChange: () => {} },
};
