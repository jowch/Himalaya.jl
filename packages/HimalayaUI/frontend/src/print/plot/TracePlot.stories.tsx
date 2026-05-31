import type { Meta, StoryObj } from "@storybook/react-vite";
import { TracePlot, type TraceModel } from "./TracePlot";
import type { PlotPeak } from "./marks/PlotPeaks";
import { realTraces } from "../fixtures/realTraces";
import { realMembers, transitionSeries } from "../fixtures/realSeriesMembers";
import type { SeriesMember } from "../../api";

function modelFor(member: SeriesMember): TraceModel {
  const trace = realTraces[member.exposure_id as number]!;
  const peaks: PlotPeak[] = (member.snapshot?.effective_peaks ?? []).map(
    (p) => ({
      id: p.id,
      q: p.q,
      intensity: p.intensity,
      source: p.source,
    }),
  );
  return {
    trace,
    peaks,
    phase: member.snapshot?.confirmed_index?.phase ?? null,
  };
}

const heroModel = modelFor(realMembers[0]!); // exp 65, Ia3d

const meta = {
  title: "plot/TracePlot",
  component: TracePlot,
  parameters: { layout: "padded" },
} satisfies Meta<typeof TracePlot>;

export default meta;
type Story = StoryObj<typeof meta>;

// Focus hero: axes on, full interaction (zoom + click-to-add/select peaks).
export const Hero: Story = {
  args: {
    traces: [heroModel],
    height: 360,
    interaction: {
      onXDomain: (d) => console.log("xDomain", d),
      onAddPeak: (q) => console.log("addPeak", q),
      onClickPeak: (id, alt) => console.log("clickPeak", id, alt),
    },
    paperColor: "var(--color-paper)",
  },
};

// Mini: axes + interaction off, small — cheap enough for a masonry of cards.
export const Mini: Story = {
  args: {
    traces: [heroModel],
    height: 64,
    width: 200,
    axes: false,
    interaction: false,
  },
};

// Overlay: the three-member Sample-9 transition in shared scales (a preview of
// the Plan #2 band layout — here simply overlaid).
export const TransitionOverlay: Story = {
  args: {
    traces: transitionSeries.map(modelFor),
    height: 360,
  },
};
