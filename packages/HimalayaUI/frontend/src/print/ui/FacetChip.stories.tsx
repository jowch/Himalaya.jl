import type { Meta, StoryObj } from "@storybook/react-vite";
import { FacetChip } from "./FacetChip";

const meta = {
  title: "ui/FacetChip",
  component: FacetChip,
  args: { label: "Beamtime", onClick: () => {} },
} satisfies Meta<typeof FacetChip>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {};
