import type { Meta, StoryObj } from "@storybook/react-vite";
import { Popover } from "./Popover";

const meta = {
  title: "ui/Popover",
  component: Popover,
  parameters: {
    docs: {
      description: {
        component:
          "A click-triggered LIGHT-plate popover for ARBITRARY (possibly interactive) content — the ARIA counterpart to the hover Tooltip, which may hold only non-interactive text. Opens on click, closes on Escape / outside-pointerdown / re-click.",
      },
    },
  },
  decorators: [
    (Story) => (
      <div className="flex justify-center p-16">
        <Story />
      </div>
    ),
  ],
} satisfies Meta<typeof Popover>;

export default meta;
type Story = StoryObj<typeof meta>;

/** Click the trigger to reveal interactive content below it. */
export const Default: Story = {
  args: {
    label: "Details",
    trigger: (
      <button className="text-ink-faint text-xs hover:text-ink cursor-pointer">
        show more
      </button>
    ),
    children: (
      <div className="flex flex-col gap-1">
        <button className="text-ink-soft hover:text-ink text-sm text-left">
          Action one
        </button>
        <button className="text-ink-soft hover:text-ink text-sm text-left">
          Action two
        </button>
      </div>
    ),
  },
};
