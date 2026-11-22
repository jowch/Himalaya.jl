import type { Meta, StoryObj } from "@storybook/react-vite";
import { ModalShell } from "./ModalShell";

// Rendered OPEN with a placeholder title + body so the frame is visible in
// Storybook. ModalShell renders into the normal React tree (a fixed-position
// scrim + frame, no portal), so the story preview shows it directly.
const meta = {
  title: "ui/ModalShell",
  component: ModalShell,
  args: {
    open: true,
    onClose: () => {},
    "aria-label": "Example dialog",
    children: (
      <div className="p-5">
        <h2 className="text-h3 font-semibold text-ink">Example dialog</h2>
        <p className="mt-2 text-body text-ink-soft">
          Placeholder body content. Press Escape or click the backdrop to close
          (no-op in this story).
        </p>
      </div>
    ),
  },
} satisfies Meta<typeof ModalShell>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {};

export const Small: Story = { args: { size: "sm" } };
export const Large: Story = { args: { size: "lg" } };
export const TopAligned: Story = { args: { align: "top" } };
export const Drawer: Story = { args: { variant: "drawer" } };
