import { useState } from "react";
import type { Meta, StoryObj } from "@storybook/react-vite";
import { MemberList, type MemberDatum } from "./MemberList";

const meta: Meta<typeof MemberList> = {
  title: "components/MemberList",
  component: MemberList,
};

export default meta;
type Story = StoryObj<typeof MemberList>;

const MEMBERS: MemberDatum[] = [
  { key: "1:1.5",  phases: [],                   variableValue: "1:1.5",  dataLine: "no Bragg peaks · q₁ —" },
  { key: "1:1",    phases: ["Lamellar"],          variableValue: "1:1",    dataLine: "d = 55 Å · q₁ 0.114 Å⁻¹" },
  { key: "1:0.75", phases: ["Lamellar"],          variableValue: "1:0.75", dataLine: "d = 57 Å · q₁ 0.110 Å⁻¹" },
  { key: "1:0.5",  phases: ["Pn3m", "Lamellar"],  variableValue: "1:0.5",  dataLine: "a 195 · d 60 Å · q₁ 0.057 Å⁻¹" },
  { key: "1:0.25", phases: ["Pn3m"],              variableValue: "1:0.25", dataLine: "a = 198 Å · q₁ 0.063 Å⁻¹" },
  { key: "1:0.1",  phases: ["Pn3m"],              variableValue: "1:0.1",  dataLine: "a = 202 Å · q₁ 0.062 Å⁻¹" },
  { key: "1:0",    phases: ["Pn3m"],              variableValue: "1:0",    dataLine: "a = 205 Å · q₁ 0.061 Å⁻¹" },
];

function OpenDemo() {
  const [hoveredKey, setHoveredKey] = useState<string | undefined>(undefined);
  return (
    <div style={{ width: 320 }}>
      <MemberList
        members={MEMBERS}
        {...(hoveredKey !== undefined ? { hoveredKey } : {})}
        onHoverMember={setHoveredKey}
      />
    </div>
  );
}

export const Default: Story = {
  render: () => <OpenDemo />,
  args: { members: MEMBERS },
};
