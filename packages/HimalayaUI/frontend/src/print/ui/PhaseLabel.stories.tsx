import type { Meta, StoryObj } from "@storybook/react-vite";
import { PhaseLabel } from "./PhaseLabel";

const meta = {
  title: "ui/PhaseLabel",
  component: PhaseLabel,
  args: { phase: "Pn3m", children: "Pn3m" },
} satisfies Meta<typeof PhaseLabel>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {};

const PHASES = ["Pn3m", "Im3m", "Ia3d", "Lamellar"];

export const AllPhases: Story = {
  render: () => (
    <div className="inline-flex flex-col gap-2">
      {PHASES.map((p) => (
        <PhaseLabel key={p} phase={p} className="text-data-strong font-bold">
          {p}
        </PhaseLabel>
      ))}
    </div>
  ),
};

export const Roles: Story = {
  render: () => (
    <div className="inline-flex flex-col gap-2">
      <PhaseLabel phase="Pn3m" className="text-display">
        Pn3m
      </PhaseLabel>
      <PhaseLabel phase="Pn3m" className="text-data-strong font-bold">
        Pn3m
      </PhaseLabel>
    </div>
  ),
};
