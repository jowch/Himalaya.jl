import { useState } from "react";
import type { Meta, StoryObj } from "@storybook/react-vite";
import { Button, SegmentedControl } from "../ui";
import { ToolBar } from "./ToolBar";

const meta = {
  title: "components/ToolBar",
  component: ToolBar,
  // `children` is required in ToolBarProps; stories use render() overrides,
  // so supply a null placeholder here to satisfy the CSF3 args constraint.
  args: { children: null },
} satisfies Meta<typeof ToolBar>;

export default meta;
type Story = StoryObj<typeof meta>;

/** Canonical Focus composition: scale toggle + Auto-fit action + Peak armed toggle. */
function FocusToolBarDemo() {
  const [scale, setScale] = useState<"log" | "lin">("log");
  const [peakArmed, setPeakArmed] = useState(false);
  return (
    <ToolBar>
      <SegmentedControl
        options={[
          { value: "log", label: "log q" },
          { value: "lin", label: "linear q" },
        ]}
        value={scale}
        onChange={setScale}
        aria-label="q scale"
      />
      <Button variant="ghost" onClick={() => {}}>Auto-fit</Button>
      <Button armed={peakArmed} onClick={() => setPeakArmed((p) => !p)}>+ Peak</Button>
    </ToolBar>
  );
}

export const Default: Story = {
  render: () => <FocusToolBarDemo />,
};

/** Two ghost action buttons only — no segmented control. */
function JustButtonsDemo() {
  return (
    <ToolBar>
      <Button variant="ghost">Export</Button>
      <Button variant="ghost">Share</Button>
    </ToolBar>
  );
}

export const JustButtons: Story = {
  render: () => <JustButtonsDemo />,
};

/** Scale toggle only — e.g. a stripped-down mini toolbar. */
function ScaleOnlyDemo() {
  const [scale, setScale] = useState<"log" | "lin">("log");
  return (
    <ToolBar>
      <SegmentedControl
        options={[
          { value: "log", label: "log q" },
          { value: "lin", label: "linear q" },
        ]}
        value={scale}
        onChange={setScale}
        aria-label="q scale"
      />
    </ToolBar>
  );
}

export const ScaleOnly: Story = {
  render: () => <ScaleOnlyDemo />,
};
