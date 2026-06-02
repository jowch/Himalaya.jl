import type { Meta, StoryObj } from "@storybook/react-vite";
import { ResidualChart } from "./ResidualChart";
import { PN3M, IM3M } from "./comb.fixtures";

const meta: Meta<typeof ResidualChart> = { title: "comb/ResidualChart", component: ResidualChart };
export default meta;
type Story = StoryObj<typeof ResidualChart>;

const frame = "rounded border border-hair bg-plate p-3";

export const Wide: Story = {
  render: () => <div className={frame} style={{ width: 760 }}><ResidualChart assigned={[PN3M, IM3M]} maxWidth={760} /></div>,
};

export const NarrowScrolling: Story = {
  render: () => <div className={frame} style={{ width: 360 }}><ResidualChart assigned={[PN3M, IM3M]} maxWidth={360} /></div>,
};

export const OverflowSignal: Story = {
  render: () => <div className={frame} style={{ width: 760 }}><ResidualChart assigned={[IM3M]} maxWidth={760} /></div>,
};
