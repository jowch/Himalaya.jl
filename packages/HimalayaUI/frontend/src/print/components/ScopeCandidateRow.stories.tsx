import type { Meta, StoryObj } from "@storybook/react-vite";
import { ScopeCandidateRow } from "./ScopeCandidateRow";
import { realTraces } from "../fixtures/realTraces";

const meta = {
  title: "components/ScopeCandidateRow",
  component: ScopeCandidateRow,
  args: {
    name: "Lipid 1-1 + LL37 1:1",
    sampleId: "smp_37",
    phase: "Pn3m",
    trace: realTraces[37]!,
    why: (
      <>
        has LL37, but the{" "}
        <strong className="text-accent font-semibold">1-1 lipid line</strong>{" "}
        — its own series?
      </>
    ),
  },
} satisfies Meta<typeof ScopeCandidateRow>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {};
