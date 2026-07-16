import { useState } from "react";
import type { Meta, StoryObj } from "@storybook/react-vite";
import { FlagButton } from "./FlagButton";

const meta = {
  title: "ui/FlagButton",
  component: FlagButton,
  args: { value: "1 : 0.25" },
} satisfies Meta<typeof FlagButton>;

export default meta;
type Story = StoryObj<typeof meta>;

/** Kept read (will be committed): hover reveals the dotted affordance; clicking skips it. */
export const Default: Story = {};

/** Skipped read: accent value, dashed underline, "skipped" caption; clicking restores it. */
export const Flagged: Story = { args: { flagged: true } };

/** Click the value to flip between the skipped and kept looks. */
export const Toggling: Story = {
  render: () => {
    const [flagged, setFlagged] = useState(true);
    return (
      <FlagButton
        value="1 : 0.25"
        flagged={flagged}
        onClick={() => setFlagged((f) => !f)}
      />
    );
  },
};
