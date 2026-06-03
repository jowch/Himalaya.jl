import type { Meta, StoryObj } from "@storybook/react-vite";
import { Field } from "./Field";

const meta = {
  title: "ui/Field",
  component: Field,
  args: { value: "LL37 : lipid ratio" },
} satisfies Meta<typeof Field>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {};

export const Placeholder: Story = {
  args: { value: "", placeholder: "Choose a variable" },
};
