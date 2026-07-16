import type { Meta, StoryObj } from "@storybook/react-vite";
import { SeriesMemberRow } from "./SeriesMemberRow";

const meta: Meta<typeof SeriesMemberRow> = {
  title: "components/SeriesMemberRow",
  component: SeriesMemberRow,
};
export default meta;

type Story = StoryObj<typeof SeriesMemberRow>;

const noop = () => {};

export const Gallery: Story = {
  render: () => (
    <div className="flex flex-col gap-0.5">
      <SeriesMemberRow
        phases={["Pn3m"]}
        variableValue="1:0"
        dataLine="a = 205 Å · q₁ 0.061 Å⁻¹"
        onHover={noop}
        onLeave={noop}
      />
      <SeriesMemberRow
        phases={["Pn3m", "Lamellar"]}
        variableValue="1:0.5"
        dataLine="a 195 · d 60 Å · q₁ 0.057 Å⁻¹"
        onHover={noop}
        onLeave={noop}
      />
      <SeriesMemberRow
        phases={[]}
        variableValue="1:1.5"
        dataLine="no Bragg peaks · q₁ —"
        onHover={noop}
        onLeave={noop}
      />
    </div>
  ),
};

export const Hot: Story = {
  args: {
    phases: ["Pn3m"],
    variableValue: "1:0",
    dataLine: "a = 205 Å · q₁ 0.061 Å⁻¹",
    hot: true,
  },
};
