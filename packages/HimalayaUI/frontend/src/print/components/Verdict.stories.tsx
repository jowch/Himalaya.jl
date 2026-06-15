import type { Meta, StoryObj } from "@storybook/react-vite";
import { useState } from "react";
import { Verdict } from "./Verdict";

const meta = {
  title: "components/Verdict",
  component: Verdict,
  args: {
    dropped: false,
  },
} satisfies Meta<typeof Verdict>;

export default meta;
type Story = StoryObj<typeof meta>;

/** Neither verdict yet — the honest default (SA-SCREENED). */
export const Unscreened: Story = {
  args: { dropped: false, onToggle: () => {}, onToggleKeep: () => {} },
  render: (args) => (
    <div className="bg-paper p-4 w-96">
      <Verdict {...args} />
    </div>
  ),
};

/** Explicitly accepted via the Keep verb. */
export const Kept: Story = {
  args: { dropped: false, kept: true, onToggle: () => {}, onToggleKeep: () => {} },
  render: (args) => (
    <div className="bg-paper p-4 w-96">
      <Verdict {...args} />
    </div>
  ),
};

export const Dropped: Story = {
  args: { dropped: true, onToggle: () => {}, onToggleKeep: () => {} },
  render: (args) => (
    <div className="bg-paper p-4 w-96">
      <Verdict {...args} />
    </div>
  ),
};

// Stateful wrapper so the tri-state toggles work interactively in Storybook.
// Mirrors the LoupePage keyboard semantics: last verb wins, restore nulls.
function VerdictToggle(): JSX.Element {
  const [status, setStatus] = useState<"accepted" | "rejected" | null>(null);
  return (
    <div className="bg-paper p-4 w-96">
      <Verdict
        dropped={status === "rejected"}
        kept={status === "accepted"}
        onToggle={() => setStatus((s) => (s === "rejected" ? null : "rejected"))}
        onToggleKeep={() => setStatus((s) => (s === "accepted" ? null : "accepted"))}
      />
    </div>
  );
}

export const Interactive: Story = {
  render: () => <VerdictToggle />,
};
