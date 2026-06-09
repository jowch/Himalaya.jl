import type { Meta, StoryObj } from "@storybook/react-vite";
import { Checkbox } from "./Checkbox";

const meta = {
  title: "ui/Checkbox",
  component: Checkbox,
  decorators: [
    (Story) => (
      <div className="inline-flex items-center gap-3 bg-paper p-4 rounded-sm">
        <Story />
      </div>
    ),
  ],
} satisfies Meta<typeof Checkbox>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Unchecked: Story = {
  args: { checked: false, "aria-label": "select sample" },
};

export const Checked: Story = {
  args: { checked: true, "aria-label": "select sample" },
};

export const Indeterminate: Story = {
  args: { indeterminate: true, "aria-label": "select sample" },
};

export const Disabled: Story = {
  args: { disabled: true, "aria-label": "select sample" },
};

export const DisabledChecked: Story = {
  args: { checked: true, disabled: true, "aria-label": "select sample" },
};
