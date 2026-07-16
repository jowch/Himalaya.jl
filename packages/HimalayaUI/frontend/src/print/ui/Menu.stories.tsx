import type { Meta, StoryObj } from "@storybook/react-vite";
import { Menu } from "./Menu";

const meta = {
  title: "ui/Menu",
  component: Menu,
  args: {
    "aria-label": "Sort by",
    open: true,
    onSelect: () => {},
    onClose: () => {},
    options: [
      { value: "score", label: "Score" },
      { value: "name", label: "Name" },
      { value: "date", label: "Date added" },
    ],
  },
  decorators: [
    (Story) => (
      <div className="relative h-40">
        <Story />
      </div>
    ),
  ],
} satisfies Meta<typeof Menu>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Open: Story = { args: { activeValue: "score" } };

export const WithDisabled: Story = {
  args: {
    options: [
      { value: "score", label: "Score" },
      { value: "name", label: "Name" },
      { value: "date", label: "Date added", disabled: true },
    ],
    activeValue: "name",
  },
};
