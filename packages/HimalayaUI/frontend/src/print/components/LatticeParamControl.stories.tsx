import type { Meta, StoryObj } from "@storybook/react-vite";
import { useState } from "react";
import { LatticeParamControl } from "./LatticeParamControl";

const meta: Meta<typeof LatticeParamControl> = {
  title: "Print/Components/LatticeParamControl",
  component: LatticeParamControl,
  parameters: { layout: "padded" },
};

export default meta;
type Story = StoryObj<typeof LatticeParamControl>;

function InteractiveHarness(): JSX.Element {
  const [value, setValue] = useState("252");
  return (
    <div className="flex flex-col gap-4 w-80">
      <LatticeParamControl
        value={value}
        min={120}
        max={360}
        step={1}
        onValueChange={setValue}
        unit="Å"
        label="lattice a"
      />
      <div className="text-caption text-ink-faint font-mono">value: {value}</div>
    </div>
  );
}

export const Default: Story = {
  render: () => <InteractiveHarness />,
};
