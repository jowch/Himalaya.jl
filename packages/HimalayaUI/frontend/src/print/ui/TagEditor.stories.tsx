import type { Meta, StoryObj } from "@storybook/react-vite";
import { TagEditor } from "./TagEditor";

const meta = {
  title: "ui/TagEditor",
  component: TagEditor,
  args: { onCommit: (tag) => console.log("commit", tag) },
} satisfies Meta<typeof TagEditor>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {};

export const WithKnownKeys: Story = {
  args: { knownKeys: ["temperature", "buffer", "lipid"] },
};
