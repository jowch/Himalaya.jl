import type { Meta, StoryObj } from "@storybook/react-vite";
import { AssignmentRail } from "./AssignmentRail";
import { AssignmentCart } from "./AssignmentCart";
import { PhaseBlock } from "./PhaseBlock";
import { CandidateRow, CandidateList } from "./CandidateRow";

const noop = (): void => {};

const meta: Meta<typeof AssignmentRail> = {
  title: "components/AssignmentRail",
  component: AssignmentRail,
};

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {
  render: () => (
    <div style={{ width: 340, height: 620 }}>
      <AssignmentRail
        assignmentCount={2}
        assignment={
          <AssignmentCart onCustomIndex={noop}>
            <PhaseBlock
              phase="Pn3m"
              score={0.91}
              meta="a = 197 Å · 3 reflections"
              onRemove={noop}
            />
            <PhaseBlock
              phase="Im3m"
              score={0.74}
              meta="a = 252 Å · 2 reflections"
              onRemove={noop}
            />
          </AssignmentCart>
        }
        candidates={
          <CandidateList>
            <CandidateRow
              phase="Pn3m"
              score={0.91}
              why="explains 3 peaks · in the call"
              selected
            />
            <CandidateRow
              phase="Im3m"
              score={0.74}
              why="explains 2 peaks"
            />
            <CandidateRow
              phase="Ia3d"
              score={0.55}
              why="explains 2 peaks"
              bonnet
            />
          </CandidateList>
        }
        candidatesNote={
          <>
            Ranked by fit; <b>Bonnet</b> bumps a coexisting cubic whose lattice
            matches the Gauss–Bonnet ratio. Add the ones that make sense.
          </>
        }
      />
    </div>
  ),
};
