import type { Meta, StoryObj } from "@storybook/react-vite";
import { Slider } from "./Slider";

const meta = {
  title: "ui/Slider",
  component: Slider,
  args: { value: 0.4, min: 0, max: 1.4, step: 0.1, onChange: () => {}, label: "Offset" },
  decorators: [
    (Story) => (
      <div className="w-64">
        <Story />
      </div>
    ),
  ],
} satisfies Meta<typeof Slider>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {};
export const WithValue: Story = { args: { valueDisplay: "0.40" } };
