import type { Meta, StoryObj } from "@storybook/react-vite";
import { useState } from "react";
import { CandidateRow, CandidateList } from "./CandidateRow";

const meta = {
  title: "components/CandidateRow",
  component: CandidateRow,
  args: {
    phase: "Pn3m",
    score: 0.91,
    why: "dominant — three sharp reflections",
  },
} satisfies Meta<typeof CandidateRow>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {};

export const Selected: Story = {
  args: { selected: true },
};

export const WithBonnet: Story = {
  args: {
    phase: "Im3m",
    score: 0.84,
    bonnet: true,
    why: "explains the two subtle peaks; a matches Bonnet",
  },
};

export const Custom: Story = {
  args: {
    score: null,
    why: "hand-drawn index — no auto-fit score",
  },
};

// Stateful wrapper so toggling works interactively in Storybook.
function ListDemo(): JSX.Element {
  const CANDIDATES: Array<{
    phase: string;
    score: number | null;
    why: string;
    bonnet?: boolean;
  }> = [
    { phase: "Pn3m",       score: 0.91, why: "dominant — three sharp reflections" },
    { phase: "Im3m",       score: 0.84, why: "explains the two subtle peaks; a matches Bonnet", bonnet: true },
    { phase: "Ia3d",       score: 0.61, why: "partial fit — two of three expected peaks present" },
    { phase: "Hexagonal",  score: 0.48, why: "single broad reflection; weak evidence" },
  ];
  const [selected, setSelected] = useState<Set<string>>(new Set(["Pn3m"]));

  return (
    <CandidateList className="w-72">
      {CANDIDATES.map((c) => (
        <CandidateRow
          key={c.phase}
          phase={c.phase}
          score={c.score}
          why={c.why}
          {...(c.bonnet ? { bonnet: true } : {})}
          selected={selected.has(c.phase)}
          onToggle={() =>
            setSelected((prev) => {
              const next = new Set(prev);
              if (next.has(c.phase)) { next.delete(c.phase); } else { next.add(c.phase); }
              return next;
            })
          }
        />
      ))}
    </CandidateList>
  );
}

export const List: Story = {
  render: () => <ListDemo />,
};
