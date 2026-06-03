import type { Meta, StoryObj } from "@storybook/react-vite";
import { NoticePill } from "./NoticePill";

const meta: Meta<typeof NoticePill> = {
  title: "ui/NoticePill",
  component: NoticePill,
};

export default meta;

type Story = StoryObj<typeof NoticePill>;

export const New: Story = {
  render: () => <NoticePill tone="new">+2 new match</NoticePill>,
};

export const Draft: Story = {
  render: () => <NoticePill tone="draft">Draft</NoticePill>,
};
