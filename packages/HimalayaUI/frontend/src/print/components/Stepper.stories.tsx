import type { Meta, StoryObj } from "@storybook/react-vite";
import { Stepper } from "./Stepper";

const meta = {
  title: "components/Stepper",
  component: Stepper,
  args: {
    label: "Lipid 1‑2 + LL37 1:0.5",
    position: "sample 4 of 9",
  },
} satisfies Meta<typeof Stepper>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {
  args: {
    onPrev: () => {},
    onNext: () => {},
  },
};

export const FirstSample: Story = {
  args: {
    label: "DPPC neat",
    position: "sample 1 of 9",
    prevDisabled: true,
    onNext: () => {},
  },
};

export const LastSample: Story = {
  args: {
    label: "DPPC + cholesterol 3:1",
    position: "sample 9 of 9",
    onPrev: () => {},
    nextDisabled: true,
  },
};

export const NoSubtitle: Story = {
  args: {
    label: "Lipid 1‑2 + LL37 1:0.5",
    position: undefined,
    onPrev: () => {},
    onNext: () => {},
  },
};
