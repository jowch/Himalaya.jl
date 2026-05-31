import type { Meta, StoryObj } from "@storybook/react-vite";
import { MetaList } from "./MetaList";

const meta = {
  title: "ui/MetaList",
  component: MetaList,
} satisfies Meta<typeof MetaList>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {
  args: {
    entries: [
      { key: "frame", value: "037" },
      { key: "integration", value: "1.2 s" },
      { key: "collected", value: "2024-03-14" },
      // Placeholder string; in real use this row's value would host a
      // SignalBars element rather than a bare label.
      { key: "signal", value: "strong" },
    ],
  },
};
