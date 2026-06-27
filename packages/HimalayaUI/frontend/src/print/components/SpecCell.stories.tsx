import type { Meta, StoryObj } from "@storybook/react-vite";
import { SpecCell } from "./SpecCell";

const meta = {
  title: "components/SpecCell",
  component: SpecCell,
  args: {
    name: "Sample A",
    sampleId: "s-001",
  },
} satisfies Meta<typeof SpecCell>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {};
