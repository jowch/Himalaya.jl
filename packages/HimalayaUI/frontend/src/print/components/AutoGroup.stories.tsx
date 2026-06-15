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
        contact sheet. Himalaya grouped them by their stored{" "}
        <strong className="text-ink font-semibold">LL37 : lipid ratio</strong> value. Click any
        value to skip that sample from the write.
      </>
    ),
  },
};

// The builder's compose card as actually used (BU-AUTOGROUP-STALE): no
// "Auto-grouped" title — this is a user-owned saved series, not a machine
// suggestion — a neutral grouping body, and real verbs (Save changes / Edit).
export const Compose: Story = {
  args: {
    variant: "compose",
    children: (
      <>
        <strong className="text-ink font-semibold">6 samples</strong>, ordered by{" "}
        <strong className="text-ink font-semibold">LL37 : lipid ratio</strong>.
      </>
    ),
    actions: [{ label: "Save changes" }, { label: "Edit", muted: true }],
  },
};
