import type { Meta, StoryObj } from "@storybook/react-vite";
import { useState } from "react";
import { CullBar } from "./CullBar";

const meta = {
  title: "components/CullBar",
  component: CullBar,
  args: {
    count: 0,
    show: false,
  },
} satisfies Meta<typeof CullBar>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Selected: Story = {
  args: {
    count: 3,
    show: true,
    onReject: () => {},
    onRestore: () => {},
    onClear: () => {},
  },
  render: (args) => (
    <div className="bg-paper h-[220px] relative">
      <CullBar {...args} />
    </div>
  ),
};

export const Hidden: Story = {
  args: {
    count: 0,
    show: false,
  },
  render: (args) => (
    <div className="bg-paper h-[220px] relative">
      <p className="text-sm text-ink-soft p-4">Bar is hidden (opacity-0). Set show=true to reveal.</p>
      <CullBar {...args} />
    </div>
  ),
};

function InteractiveCullBar(): JSX.Element {
  const [count, setCount] = useState(4);
  return (
    <div className="bg-paper h-[220px] relative">
      <p className="text-sm text-ink-soft p-4">
        {count} frames selected. Drop or Clear to hide bar.{" "}
        <button
          className="underline"
          onClick={() => setCount(4)}
        >
          Reset
        </button>
      </p>
      <CullBar
        count={count}
        show={count > 0}
        onReject={() => setCount(0)}
        onRestore={() => {}}
        onClear={() => setCount(0)}
      />
    </div>
  );
}

export const Interactive: Story = {
  render: () => <InteractiveCullBar />,
};
