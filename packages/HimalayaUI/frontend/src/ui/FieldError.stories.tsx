import type { Meta, StoryObj } from "@storybook/react-vite";
import { FieldError } from "./FieldError";

const meta = {
  title: "ui/FieldError",
  component: FieldError,
  args: { children: "Use {name} exactly once, for example {name}.tif." },
} satisfies Meta<typeof FieldError>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {};
