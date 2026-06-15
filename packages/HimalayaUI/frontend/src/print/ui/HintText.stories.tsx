import type { Meta, StoryObj } from "@storybook/react-vite";
import { HintText } from "./HintText";

const meta = {
  title: "ui/HintText",
  component: HintText,
  args: { children: "Peaks are fit automatically; drag to nudge a marker." },
} satisfies Meta<typeof HintText>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {};
