import type { Meta, StoryObj } from "@storybook/react-vite";
import { Popover } from "./Popover";
import { TagPill } from "./TagPill";

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

/** The popover can hold real chips (this is what TagList's overflow uses). */
export const WithTagPills: Story = {
  args: {
    label: "More tags",
    side: "bottom",
    trigger: (
      <button className="text-ink-faint text-xs hover:text-ink cursor-pointer">
        +3 more
      </button>
    ),
    children: (
      <div className="flex flex-wrap gap-1.5">
        <TagPill tag={{ key: "buffer", value: "PBS" }} />
        <TagPill tag={{ key: "lipid", value: "DOPC" }} />
        <TagPill tag={{ key: "pH", value: "7.4" }} />
      </div>
    ),
  },
};
