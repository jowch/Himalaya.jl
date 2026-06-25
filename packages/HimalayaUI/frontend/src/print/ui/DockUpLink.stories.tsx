import type { Meta, StoryObj } from "@storybook/react-vite";
import { DockUpLink } from "./DockUpLink";

const meta = {
  title: "ui/DockUpLink",
  component: DockUpLink,
} satisfies Meta<typeof DockUpLink>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Button: Story = {
  args: { label: "Corpus", onClick: () => {} },
};

export const Link: Story = {
  args: { label: "Experiments", href: "/experiments", onClick: () => {} },
};
