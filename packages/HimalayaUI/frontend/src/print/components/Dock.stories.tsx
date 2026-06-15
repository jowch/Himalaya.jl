import type { Meta, StoryObj } from "@storybook/react-vite";
import { useState } from "react";
import { Dock } from "./Dock";

const meta = {
  title: "components/Dock",
  component: Dock,
  args: {
    offset: 1.2,
  },
} satisfies Meta<typeof Dock>;

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
      <Dock {...args} />
    </div>
  ),
};

function InteractiveDock(): JSX.Element {
  const [offset, setOffset] = useState(1.2);
  return (
    <div className="bg-paper min-h-screen relative">
      <p className="text-sm text-ink-soft p-4">Drag the slider to change the trace offset.</p>
      <Dock offset={offset} onOffsetChange={setOffset} />
    </div>
  );
}

export const Interactive: Story = {
  args: {
    onOffsetChange: () => {},
  },
  render: () => <InteractiveDock />,
};
