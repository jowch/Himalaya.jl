import type { Meta, StoryObj } from "@storybook/react-vite";
import { KeptCell } from "./KeptCell";

const meta = {
  title: "components/KeptCell",
  component: KeptCell,
  args: {
    kept: 4,
    total: 5,
  },
} satisfies Meta<typeof KeptCell>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Partial: Story = {
  args: { kept: 4, total: 5 },
};

export const AllKept: Story = {
  args: { kept: 5, total: 5 },
};
