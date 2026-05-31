import type { Meta, StoryObj } from "@storybook/react-vite";
import { TagList } from "./TagList";

const meta = {
  title: "ui/TagList",
  component: TagList,
  args: {
    tags: [
      { key: "LL37" },
      { key: "temperature", value: "37C" },
      { key: "buffer", value: "PBS" },
    ],
  },
} satisfies Meta<typeof TagList>;

export default meta;
type Story = StoryObj<typeof meta>;

export const WithTags: Story = {
  args: {
    tags: [
      { key: "LL37" },
      { key: "temperature", value: "37C" },
      { key: "buffer", value: "PBS" },
    ],
    onAdd: () => {},
  },
};
export const Empty: Story = {
  args: { tags: [], onAdd: () => {} },
};
export const Editable: Story = {
  args: {
    tags: [{ key: "lipid", value: "DOPC" }],
    editable: true,
    onAdd: () => {},
    onRemove: () => {},
  },
};
