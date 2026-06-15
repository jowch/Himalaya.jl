import type { Meta, StoryObj } from "@storybook/react-vite";
import { RailBack } from "./RailBack";

const meta = {
  title: "components/RailBack",
  component: RailBack,
  args: {
    label: "Compose",
  },
} satisfies Meta<typeof RailBack>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {
  args: {
    onClick: () => {},
  },
  render: (args) => (
    <div className="bg-paper min-h-screen relative">
      <p className="text-sm text-ink-soft p-4">
        Rail is collapsed — the tab floats on the right edge to restore it.
      </p>
      <RailBack {...args} />
    </div>
  ),
};

export const CustomLabel: Story = {
  args: {
    label: "Edit",
    onClick: () => {},
  },
  render: (args) => (
    <div className="bg-paper min-h-screen relative">
      <RailBack {...args} />
    </div>
  ),
};
