import type { Meta, StoryObj } from "@storybook/react-vite";
import { DockStepper } from "./DockStepper";

const meta = {
  title: "ui/DockStepper",
  component: DockStepper,
} satisfies Meta<typeof DockStepper>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Sample: Story = {
  args: {
    label: "Sample",
    axis: "vertical",
    testIdBase: "sample",
    onPrev: () => {},
    onNext: () => {},
    count: "3 / 12",
  },
};

export const Frame: Story = {
  args: {
    label: "Frame",
    axis: "horizontal",
    testIdBase: "frame",
    onPrev: () => {},
    onNext: () => {},
    count: "2 / 5",
    countWidthClass: "min-w-[2.75rem]",
    prevDisabled: true,
  },
};

export const NoReadout: Story = {
  args: {
    label: "Sample",
    axis: "vertical",
    testIdBase: "sample",
    onPrev: () => {},
    onNext: () => {},
  },
};
