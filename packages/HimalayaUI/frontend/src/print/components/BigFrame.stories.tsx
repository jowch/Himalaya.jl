import type { Meta, StoryObj } from "@storybook/react-vite";
import { BigFrame } from "./BigFrame";
import thumb37 from "../fixtures/thumbs/37.png?url";

const meta: Meta<typeof BigFrame> = {
  title: "components/BigFrame",
  component: BigFrame,
};

export default meta;
type Story = StoryObj<typeof meta>;

export const Kept: Story = {
  render: () => (
    <div className="bg-paper p-6 w-[520px]">
      <BigFrame src={thumb37} caption="frame 37 · 0.40 s" />
    </div>
  ),
};

export const Dropped: Story = {
  render: () => (
    <div className="bg-paper p-6 w-[520px]">
      <BigFrame src={thumb37} caption="frame 37 · 0.40 s" rejected />
    </div>
  ),
};

export const Empty: Story = {
  render: () => (
    <div className="bg-paper p-6 w-[520px]">
      <BigFrame src={null} caption="no exposure" />
    </div>
  ),
};
