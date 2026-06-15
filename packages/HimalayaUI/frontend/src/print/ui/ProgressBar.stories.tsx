import type { Meta, StoryObj } from "@storybook/react-vite";
import { ProgressBar } from "./ProgressBar";

const meta = {
  title: "ui/ProgressBar",
  component: ProgressBar,
  decorators: [
    (Story) => (
      <div className="w-48">
        <Story />
      </div>
    ),
  ],
} satisfies Meta<typeof ProgressBar>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Partial: Story = {
  args: { value: 8, total: 12, label: "Screened 8 of 12" },
};
export const Full: Story = { args: { value: 12, total: 12 } };
export const Empty: Story = { args: { value: 0, total: 12 } };
