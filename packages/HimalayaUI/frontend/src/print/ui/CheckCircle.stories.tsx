import type { Meta, StoryObj } from "@storybook/react-vite";
import { CheckCircle } from "./CheckCircle";

const meta = {
  title: "ui/CheckCircle",
  component: CheckCircle,
  // Render on a small light surface so the empty (transparent) disc and its
  // hairline ring are visible against the page.
  decorators: [
    (Story) => (
      <div className="inline-flex items-center gap-2 bg-paper p-3 rounded-sm">
        <Story />
      </div>
    ),
  ],
} satisfies Meta<typeof CheckCircle>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Selected: Story = { args: { checked: true } };
export const Unselected: Story = { args: { checked: false } };
export const ScreenedStatus: Story = { args: { checked: true, label: "Screened" } };
