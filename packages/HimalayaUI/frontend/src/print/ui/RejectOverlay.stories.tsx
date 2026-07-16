import type { Meta, StoryObj } from "@storybook/react-vite";
import { RejectOverlay } from "./RejectOverlay";

const meta = {
  title: "ui/RejectOverlay",
  component: RejectOverlay,
} satisfies Meta<typeof RejectOverlay>;

export default meta;
type Story = StoryObj<typeof meta>;

// NOTE: rejection reads via the ✕ SHAPE + the consumer's frame dimming, never hue
// alone (checklist A). The dimming below (bg-ink-soft) is supplied by the consumer,
// not the overlay; the rejected STATE itself is announced by sibling text elsewhere
// (checklist A) so the svg is aria-hidden. Accent is the rationed reject grease-pencil
// mark (checklist B). Purely decorative, so C/D (interaction states / touch target) N/A.
export const OverThumbnail: Story = {
  render: () => (
    <div className="relative w-[62px] h-[62px] bg-ink-soft">
      <RejectOverlay />
    </div>
  ),
};

// Larger relative parent proves the ✕ scales to fill its container.
export const OverFrame: Story = {
  render: () => (
    <div className="relative w-[200px] h-[200px] bg-ink-soft">
      <RejectOverlay />
    </div>
  ),
};
