import type { Meta, StoryObj } from "@storybook/react-vite";
import { AssignmentCart } from "./AssignmentCart";
import { PhaseBlock } from "./PhaseBlock";

const noop = () => {};

const meta: Meta<typeof AssignmentCart> = {
  title: "components/AssignmentCart",
  component: AssignmentCart,
  parameters: { layout: "padded" },
};
export default meta;

type Story = StoryObj<typeof AssignmentCart>;

export const Empty: Story = {
  render: () => (
    <div className="bg-paper-sunk p-4" style={{ width: 300 }}>
      <AssignmentCart />
    </div>
  ),
};

export const SinglePhase: Story = {
  render: () => (
    <div className="bg-paper-sunk p-4" style={{ width: 300 }}>
      <AssignmentCart>
        <PhaseBlock
          phase="Im3m"
          score={0.84}
          meta="a = 14.2 nm · 5 reflections"
          series="QII alkane"
          onRemove={noop}
        />
      </AssignmentCart>
    </div>
  ),
};

export const Coexistence: Story = {
  render: () => (
    <div className="bg-paper-sunk p-4" style={{ width: 300 }}>
      <AssignmentCart>
        <PhaseBlock
          phase="Ia3d"
          score={0.71}
          meta="a = 11.4 nm · 4 reflections"
          onRemove={noop}
        />
        <PhaseBlock
          phase="Im3m"
          score={0.66}
          meta="a = 9.1 nm · 3 reflections"
          onRemove={noop}
        />
      </AssignmentCart>
    </div>
  ),
};
