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
