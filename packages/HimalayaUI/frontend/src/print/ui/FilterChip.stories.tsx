import type { Meta, StoryObj } from "@storybook/react-vite";
import { FilterChip } from "./FilterChip";

const meta = {
  title: "ui/FilterChip",
  component: FilterChip,
  args: { label: "Coexistence", active: false, onClick: () => {} },
} satisfies Meta<typeof FilterChip>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Off: Story = {
  args: { label: "Coexistence", active: false, onClick: () => {} },
};
export const On: Story = {
  args: { label: "Coexistence", active: true, onClick: () => {} },
};
