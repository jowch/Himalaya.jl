import type { Meta, StoryObj } from "@storybook/react-vite";
import { ReadingPanel } from "./ReadingPanel";

const meta: Meta<typeof ReadingPanel> = {
  title: "components/ReadingPanel",
  component: ReadingPanel,
};
export default meta;

type Story = StoryObj<typeof ReadingPanel>;

const READINGS = [
  { phase: "Pn3m", span: "1:0 → 1:0.5", lattice: "a 205 → 195 Å" },
  { phase: "Lamellar", span: "1:0.5 → 1:1", lattice: "d 60 → 55 Å" },
];

export const Default: Story = {
  args: {
    readings: READINGS,
    coexistenceNote: "coexistence at 1:0.5",
    formFactorNote: "form factor only at 1:1.5",
  },
  render: (args) => (
    <div style={{ width: 320 }}>
      <ReadingPanel {...args} />
    </div>
  ),
};

export const NoNotes: Story = {
  args: {
    readings: READINGS,
  },
  render: (args) => (
    <div style={{ width: 320 }}>
      <ReadingPanel {...args} />
    </div>
  ),
};
