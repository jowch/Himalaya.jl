import type { Meta, StoryObj } from "@storybook/react-vite";
import { StageTabs } from "./StageTabs";

const meta = {
  title: "ui/StageTabs",
  component: StageTabs,
  args: { active: "samples", onChange: () => {} },
} satisfies Meta<typeof StageTabs>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Samples: Story = { args: { active: "samples", onChange: () => {} } };
export const Index: Story = { args: { active: "index", onChange: () => {} } };
export const Series: Story = { args: { active: "series", onChange: () => {} } };
