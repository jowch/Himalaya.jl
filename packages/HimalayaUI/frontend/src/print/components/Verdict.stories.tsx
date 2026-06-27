import type { Meta, StoryObj } from "@storybook/react-vite";
import { useState } from "react";
import { Verdict } from "./Verdict";

const meta = {
  title: "components/Verdict",
  component: Verdict,
  args: {
    dropped: false,
  },
} satisfies Meta<typeof Verdict>;

export default meta;
type Story = StoryObj<typeof meta>;

/** Kept — the default verdict (not dropped). */
export const Kept: Story = {
  args: { dropped: false, onToggle: () => {} },
  render: (args) => (
    <div className="bg-paper p-4 w-96">
      <Verdict {...args} />
    </div>
  ),
};

export const Dropped: Story = {
  args: { dropped: true, onToggle: () => {} },
  render: (args) => (
    <div className="bg-paper p-4 w-96">
      <Verdict {...args} />
    </div>
  ),
};

// Stateful wrapper so the single Drop toggle flips the verdict in Storybook.
function VerdictToggle(): JSX.Element {
  const [dropped, setDropped] = useState(false);
  return (
    <div className="bg-paper p-4 w-96">
      <Verdict dropped={dropped} onToggle={() => setDropped((d) => !d)} />
    </div>
  );
}

export const Interactive: Story = {
  render: () => <VerdictToggle />,
};
