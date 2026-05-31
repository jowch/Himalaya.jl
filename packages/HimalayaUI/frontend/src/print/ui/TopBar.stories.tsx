import type { Meta, StoryObj } from "@storybook/react-vite";
import { TopBar } from "./TopBar";
import { Button } from "./Button";

const meta = {
  title: "ui/TopBar",
  component: TopBar,
  args: {
    wordmark: <span className="font-serif">Himalaya</span>,
    children: <span>tabs</span>,
    rightSlot: <Button>Account</Button>,
  },
} satisfies Meta<typeof TopBar>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {};
