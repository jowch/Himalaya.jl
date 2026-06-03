import type { Meta, StoryObj } from "@storybook/react-vite";
import { FitMetadata } from "./FitMetadata";

const meta = {
  title: "components/FitMetadata",
  component: FitMetadata,
  args: {
    landed: 7,
    total: 9,
    paramName: "d",
    paramValue: "81.2",
    unit: "Å",
  },
} satisfies Meta<typeof FitMetadata>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Snapped: Story = {
  args: {
    snapped: true,
  },
};

export const NotSnapped: Story = {
  args: {
    snapped: false,
  },
};
