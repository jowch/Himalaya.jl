import type { Meta, StoryObj } from "@storybook/react-vite";
import { TagPill } from "./TagPill";

const meta = {
  title: "ui/TagPill",
  component: TagPill,
  args: { children: "lipid-A" },
} satisfies Meta<typeof TagPill>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = { args: { children: "lipid-A" } };
export const Removable: Story = {
  args: { children: "lipid-A", onRemove: () => {} },
};
