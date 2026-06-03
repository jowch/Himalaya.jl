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

/** Resolved, re-openable value: hover reveals the dotted editable affordance. */
export const Default: Story = {};

/** Uncertain parse: accent value, dashed underline, "check the read" caption. */
export const Flagged: Story = { args: { flagged: true } };

/** Click the value to flip between the flagged and resolved looks. */
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
