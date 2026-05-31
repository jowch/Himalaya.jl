import type { Meta, StoryObj } from "@storybook/react-vite";
import { Badge } from "./Badge";
import { Button } from "./Button";

const meta = {
  title: "ui/Badge",
  component: Badge,
  args: { children: "1" },
} satisfies Meta<typeof Badge>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = { args: { children: "1" } };

export const WithinButton: Story = {
  render: () => (
    <Button>
      Notes <Badge>3</Badge>
    </Button>
  ),
};
