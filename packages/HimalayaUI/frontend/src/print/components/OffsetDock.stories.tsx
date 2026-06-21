import type { Meta, StoryObj } from "@storybook/react-vite";
import { useState } from "react";
import { OffsetDock } from "./OffsetDock";

const meta = {
  title: "components/OffsetDock",
  component: OffsetDock,
  args: {
    offset: 1.2,
  },
} satisfies Meta<typeof OffsetDock>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {
  args: {
    onOffsetChange: () => {},
  },
  render: (args) => (
    <div className="bg-paper min-h-screen relative">
      <p className="text-sm text-ink-soft p-4">
        Full-bleed reading mode — the offset dock floats bottom-right.
      </p>
      <OffsetDock {...args} />
    </div>
  ),
};

function InteractiveOffsetDock(): JSX.Element {
  const [offset, setOffset] = useState(1.2);
  return (
    <div className="bg-paper min-h-screen relative">
      <p className="text-sm text-ink-soft p-4">Drag the slider to change the trace offset.</p>
      <OffsetDock offset={offset} onOffsetChange={setOffset} />
    </div>
  );
}

export const Interactive: Story = {
  args: {
    onOffsetChange: () => {},
  },
  render: () => <InteractiveOffsetDock />,
};
