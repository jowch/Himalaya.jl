import type { Meta, StoryObj } from "@storybook/react-vite";
import { Wordmark } from "./Wordmark";

const meta = {
  title: "ui/Wordmark",
  component: Wordmark,
  args: { children: "Himalaya" },
} satisfies Meta<typeof Wordmark>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = { args: { children: "Himalaya" } };
export const WithTail: Story = { args: { children: "Himalaya", tail: "SAXS" } };
