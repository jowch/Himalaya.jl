import type { Meta, StoryObj } from "@storybook/react-vite";
import { Swatch } from "./Swatch";

const meta = {
  title: "ui/Swatch",
  component: Swatch,
  args: { phase: "Pn3m" },
} satisfies Meta<typeof Swatch>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {};

const PHASES = ["Pn3m", "Im3m", "Ia3d", "Lamellar"];

export const AllShapes: Story = {
  render: () => (
    <div className="inline-flex flex-col gap-2">
      <div className="inline-flex items-center gap-2">
        {PHASES.map((p) => (
          <Swatch key={p} phase={p} shape="square" />
        ))}
      </div>
      <div className="inline-flex items-center gap-2">
        {PHASES.map((p) => (
          <Swatch key={p} phase={p} shape="circle" />
        ))}
      </div>
    </div>
  ),
};
