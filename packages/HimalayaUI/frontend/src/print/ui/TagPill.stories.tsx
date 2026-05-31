import type { Meta, StoryObj } from "@storybook/react-vite";
import { TagPill } from "./TagPill";

const meta = {
  title: "ui/TagPill",
  component: TagPill,
  args: { tag: { key: "LL37" } },
} satisfies Meta<typeof TagPill>;

export default meta;
type Story = StoryObj<typeof meta>;

export const KeyOnly: Story = { args: { tag: { key: "LL37" } } };
export const KeyValue: Story = {
  args: { tag: { key: "temperature", value: "37C" } },
};
export const Removable: Story = {
  args: { tag: { key: "buffer", value: "PBS" }, onRemove: () => {} },
};
