import type { Meta, StoryObj } from "@storybook/react-vite";
import { ScopeSampleRow } from "./ScopeSampleRow";
import { realTraces } from "../fixtures/realTraces";

const noop = (): void => {};

const meta: Meta<typeof ScopeSampleRow> = {
  title: "components/ScopeSampleRow",
  component: ScopeSampleRow,
};

export default meta;
type Story = StoryObj<typeof meta>;

export const Ok: Story = {
  args: {
    name: "Lipid 1-2 + LL37 1:0.25",
    sampleId: "smp_07",
    trace: realTraces[37]!,
    phase: "Pn3m",
    value: "1 : 0.25",
    onToggleFlag: noop,
  },
};

export const Flagged: Story = {
  args: {
    name: "Lipid 1-2, no LL37",
    sampleId: "smp_04",
    trace: realTraces[65]!,
    phase: "Pn3m",
    value: "1 : 0",
    flagged: true,
    onToggleFlag: noop,
  },
};

export const Coexistence: Story = {
  args: {
    name: "Lipid 1-2 + LL37 1:0.5",
    sampleId: "smp_09",
    trace: realTraces[66]!,
    phase: "Im3m",
    value: "1 : 0.5",
    onToggleFlag: noop,
  },
};

export const List: Story = {
  render: () => (
    <div className="max-w-[420px]">
      <ScopeSampleRow
        name="Lipid 1-2, no LL37"
        sampleId="smp_04"
        trace={realTraces[65]!}
        phase="Pn3m"
        value="1 : 0"
        onToggleFlag={noop}
      />
      <ScopeSampleRow
        name="Lipid 1-2 + LL37 1:0.25"
        sampleId="smp_07"
        trace={realTraces[37]!}
        phase="Pn3m"
        value="1 : 0.25"
        flagged
        onToggleFlag={noop}
      />
      <ScopeSampleRow
        name="Lipid 1-2 + LL37 1:0.5"
        sampleId="smp_09"
        trace={realTraces[66]!}
        phase="Im3m"
        value="1 : 0.5"
        onToggleFlag={noop}
      />
      <ScopeSampleRow
        name="Lipid 1-2 + LL37 1:1"
        sampleId="smp_12"
        trace={realTraces[93]!}
        phase={null}
        value="1 : 1"
        onToggleFlag={noop}
      />
    </div>
  ),
};
