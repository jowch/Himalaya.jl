import type { Meta, StoryObj } from "@storybook/react-vite";
import { CombLegend } from "./CombLegend";

const meta = {
  title: "components/CombLegend",
  component: CombLegend,
  decorators: [
    (Story) => (
      <div className="bg-plate p-4">
        <Story />
      </div>
    ),
  ],
} satisfies Meta<typeof CombLegend>;

export default meta;
type Story = StoryObj<typeof meta>;

/** All four vocabulary atoms: observed, manual, predicted-absent, excluded. */
export const Default: Story = {};

/** Subset — only the two most common entries. */
export const Subset: Story = {
  args: {
    items: ["auto", "manual"],
  },
};
