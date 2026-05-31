import type { Meta, StoryObj } from "@storybook/react-vite";
import { TagList } from "./TagList";

const meta = {
  title: "ui/TagList",
  component: TagList,
  args: { tags: ["lipid-A", "37C", "run-3"] },
} satisfies Meta<typeof TagList>;

export default meta;
type Story = StoryObj<typeof meta>;

export const WithTags: Story = {
  args: { tags: ["lipid-A", "37C", "run-3"], onAdd: () => {} },
};
export const Empty: Story = {
  args: { tags: [], onAdd: () => {} },
};
export const Editable: Story = {
  args: {
    tags: ["lipid-A", "37C"],
    editable: true,
    onAdd: () => {},
    onRemove: () => {},
  },
};
