import type { Meta, StoryObj } from "@storybook/react-vite";
import { SuggestionList } from "./SuggestionList";

const meta = {
  title: "ui/SuggestionList",
  component: SuggestionList,
  args: {
    id: "folders",
    label: "Folders",
    active: 1,
    onPick: () => {},
    items: [
      { key: "a", name: "beamtime-2026-10", hint: "beamtime, 1 setup · added as October", strong: true },
      { key: "b", name: "beamtime-2026-11", hint: "beamtime, 3 setups", strong: true },
      { key: "c", name: "notes" },
    ],
  },
} satisfies Meta<typeof SuggestionList>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Folders: Story = {};
export const WithNote: Story = { args: { note: "More folders match. Keep typing to narrow the list." } };
