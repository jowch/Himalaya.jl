import type { Meta, StoryObj } from "@storybook/react-vite";
import { TopBar } from "./TopBar";
import { Button } from "./Button";
import { Wordmark } from "./Wordmark";

const meta = {
  title: "ui/TopBar",
  component: TopBar,
  args: {
    wordmark: <Wordmark tail="SAXS">Himalaya</Wordmark>,
    children: <span>tabs</span>,
    rightSlot: <Button>Account</Button>,
  },
} satisfies Meta<typeof TopBar>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {};
