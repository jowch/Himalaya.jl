import type { Meta, StoryObj } from "@storybook/react-vite";
import { FolioHeader } from "./FolioHeader";

const meta = {
  title: "components/FolioHeader",
  component: FolioHeader,
  args: {
    kicker: "Folio",
    title: "Saved series",
    subtitle: "A wall of saved series across every experiment — experiment is a filter facet, not a container.",
    count: 5,
    countLabel: "series in the folio",
  },
} satisfies Meta<typeof FolioHeader>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {};

export const NoSubtitle: Story = {
  args: {
    subtitle: undefined,
  },
};
