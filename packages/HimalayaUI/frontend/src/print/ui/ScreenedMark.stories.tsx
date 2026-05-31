import type { Meta, StoryObj } from "@storybook/react-vite";
import { ScreenedMark } from "./ScreenedMark";

const meta = {
  title: "ui/ScreenedMark",
  component: ScreenedMark,
  // Render on a small light surface so the empty (transparent) disc and its
  // hairline ring are visible against the page.
  decorators: [
    (Story) => (
      <div className="inline-flex items-center gap-2 bg-paper p-3 rounded-sm">
        <Story />
      </div>
    ),
  ],
} satisfies Meta<typeof ScreenedMark>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Screened: Story = { args: { done: true } };
export const Unscreened: Story = { args: { done: false } };
