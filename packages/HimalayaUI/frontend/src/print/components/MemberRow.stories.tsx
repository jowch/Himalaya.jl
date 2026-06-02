import type { Meta, StoryObj } from "@storybook/react";
import { MemberRow } from "./MemberRow";

const meta: Meta<typeof MemberRow> = {
  title: "components/MemberRow",
  component: MemberRow,
};
export default meta;

type Story = StoryObj<typeof MemberRow>;

export const Default: Story = {
  args: {
    name: "S-014 · 25 °C",
    sub: "sample-014 · 1:0.5",
    phase: "Pn3m",
  },
};

export const Coexistence: Story = {
  args: {
    name: "S-022 · 40 °C",
    sub: "sample-022 · 1:0.8",
    phase: "Pn3m",
    coexistWith: ["Im3m"],
  },
};

export const List: Story = {
  render: () => (
    <div className="flex flex-col gap-0.5">
      <MemberRow name="S-014 · 25 °C" sub="sample-014 · 1:0.25" phase="Pn3m" />
      <MemberRow name="S-022 · 32 °C" sub="sample-022 · 1:0.5" phase="Im3m" />
      <MemberRow name="S-031 · 40 °C" sub="sample-031 · 1:0.75" phase="Ia3d" />
      <MemberRow name="S-047 · 55 °C" sub="sample-047 · 1:1" phase="Lamellar" />
    </div>
  ),
};
