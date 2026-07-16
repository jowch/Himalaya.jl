import type { Meta, StoryObj } from "@storybook/react-vite";
import { StatusCell } from "./StatusCell";

const meta = {
  title: "components/StatusCell",
  component: StatusCell,
} satisfies Meta<typeof StatusCell>;

export default meta;
type Story = StoryObj<typeof meta>;

export const IndexedPn3m: Story = {
  args: { phase: "Pn3m" },
};

export const IndexedIm3m: Story = {
  args: { phase: "Im3m" },
};

export const Unset: Story = {
  args: { phase: null },
};
