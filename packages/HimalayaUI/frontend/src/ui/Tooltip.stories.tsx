import type { Meta, StoryObj } from "@storybook/react-vite";
import { Tooltip } from "./Tooltip";
import { Button } from "./Button";

const meta = {
  title: "ui/Tooltip",
  component: Tooltip,
  args: { content: "Reanalyze stale indices", children: <Button variant="ghost">Reanalyze</Button> },
} satisfies Meta<typeof Tooltip>;

export default meta;
type Story = StoryObj<typeof meta>;

export const OnButton: Story = {};

export const Bottom: Story = { args: { side: "bottom" } };

export const WithShortcut: Story = {
  args: { content: "Keyboard shortcuts", shortcut: "?", children: <Button variant="ghost">?</Button> },
};

export const TwoLines: Story = {
  args: {
    content: (
      <>
        <span className="block font-semibold">MO_H2O_37C · SAXS 1.5 m</span>
        <span className="block text-ink-soft">Click to hide. ⌘-click to subtract from this sample.</span>
      </>
    ),
    children: <Button variant="ghost">MO_H2O_37C</Button>,
  },
};

export const OverflowOnly: Story = {
  args: {
    overflowOnly: true,
    content: "A_very_long_sample_name_that_does_not_fit_S2401",
    children: <span className="block w-32 truncate">A_very_long_sample_name_that_does_not_fit_S2401</span>,
  },
};

export const Sweep: Story = {
  render: () => (
    <div className="flex">
      <Tooltip content="First"><Button variant="ghost">One</Button></Tooltip>
      <Tooltip content="Second"><Button variant="ghost">Two</Button></Tooltip>
      <Tooltip content="Third"><Button variant="ghost">Three</Button></Tooltip>
    </div>
  ),
};
