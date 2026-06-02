import type { Meta, StoryObj } from "@storybook/react-vite";
import { AutoGroup } from "./AutoGroup";

const meta = {
  title: "components/AutoGroup",
  component: AutoGroup,
} satisfies Meta<typeof AutoGroup>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Summary: Story = {
  args: {
    variant: "summary",
    children: (
      <>
        You selected <strong className="text-ink font-semibold">6 samples</strong> on the
        contact sheet. Himalaya grouped them from their names and read the order from the{" "}
        <strong className="text-ink font-semibold">LL37 : lipid ratio</strong> — five ratios
        parsed cleanly, one needs a look.
      </>
    ),
  },
};

export const Compose: Story = {
  args: {
    variant: "compose",
    title: "Auto-grouped",
    children: (
      <>
        Himalaya read <strong className="text-ink font-semibold">6 samples</strong> as one
        series from their names, ordered by{" "}
        <strong className="text-ink font-semibold">LL37 : lipid ratio</strong>.
      </>
    ),
    actions: [{ label: "Confirm series" }, { label: "Adjust", muted: true }],
  },
};
