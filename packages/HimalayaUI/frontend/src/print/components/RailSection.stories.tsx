import type { Meta, StoryObj } from "@storybook/react-vite";
import { RailSection } from "./RailSection";

const meta: Meta<typeof RailSection> = {
  title: "components/RailSection",
  component: RailSection,
  parameters: { layout: "padded" },
};
export default meta;

type Story = StoryObj<typeof RailSection>;

export const Default: Story = {
  render: () => (
    <div className="bg-paper-sunk p-4" style={{ width: 280 }}>
      <RailSection label="Assignment" count={2} note="Ranked by coverage × consistency score.">
        <div className="bg-plate border border-hair rounded p-3 text-body text-ink">
          Phase block placeholder
        </div>
      </RailSection>
    </div>
  ),
};

export const NoCountOrNote: Story = {
  render: () => (
    <div className="bg-paper-sunk p-4" style={{ width: 280 }}>
      <RailSection label="Candidates">
        <div className="bg-plate border border-hair rounded p-3 text-body text-ink">
          Candidate block placeholder
        </div>
        <div className="bg-plate border border-hair rounded p-3 text-body text-ink">
          Candidate block placeholder
        </div>
      </RailSection>
    </div>
  ),
};
