import type { Meta, StoryObj } from "@storybook/react";
import { ReadingRow } from "./ReadingRow";

const meta: Meta<typeof ReadingRow> = {
  title: "components/ReadingRow",
  component: ReadingRow,
};
export default meta;

type Story = StoryObj<typeof ReadingRow>;

export const Default: Story = {
  args: {
    phase: "Pn3m",
    span: "1:0 → 1:0.75",
    lattice: "a 205 → 195 Å",
  },
};

export const Panel: Story = {
  render: () => (
    <div className="flex flex-col gap-2">
      <ReadingRow phase="Pn3m" span="1:0 → 1:0.75" lattice="a 205 → 195 Å" />
      <ReadingRow phase="Im3m" span="1:0.25 → 1:0.9" lattice="a 318 → 296 Å" />
      <ReadingRow phase="Lamellar" span="1:0.5 → 1:1" lattice="d 64 → 58 Å" />
    </div>
  ),
};
