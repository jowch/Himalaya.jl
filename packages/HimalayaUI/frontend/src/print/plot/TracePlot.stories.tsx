import { useState } from "react";
import type { Meta, StoryObj } from "@storybook/react-vite";
import { TracePlot, type TraceModel } from "./TracePlot";
import type { PlotPeak } from "./marks/PlotPeaks";
import { realTraces } from "../fixtures/realTraces";
import { realMembers } from "../fixtures/realSeriesMembers";
import type { SeriesMember } from "../../api";
import { Card, Kicker, SegmentedControl, Button } from "../ui";
import { phaseColor } from "../../phases";

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
    trace: heroModel,
    height: 360,
    interaction: {
      onXDomain: (d) => console.log("xDomain", d),
      onAddPeak: (q) => console.log("addPeak", q),
      onClickPeak: (id) => console.log("clickPeak", id),
    },
    paperColor: "var(--color-paper)",
  },
};

// Mini: axes + interaction off, small — cheap enough for a masonry of cards.
export const Mini: Story = {
  args: {
    trace: heroModel,
    height: 64,
    width: 200,
    axes: false,
    interaction: false,
  },
};

// annotatedModel: heroModel peaks annotated with Ia3d colour + Miller labels.
const annotatedModel: TraceModel = {
  ...heroModel,
  peaks: heroModel.peaks.map((p, i) => ({
    ...p,
    color: phaseColor("Ia3d"),
    label: (["110", "111", "200", "211", "220", "221", "?"] as const)[i] ?? "?",
  })),
};

// HeroPlate: the full Print figure plate — closed-look primitives composing the
// publication framing around the redesigned TracePlot, matching the mockup at
// docs/redesign-mockups/2026-05-29-focus-plot.html.
export const HeroPlate: Story = {
  // args satisfies the Story constraint; render overrides the entire output.
  args: { trace: heroModel, height: 360 },
  render: () => {
    // eslint-disable-next-line react-hooks/rules-of-hooks
    const [scale, setScale] = useState<"log" | "linear">("log");
    return (
      <div style={{ maxWidth: 640, margin: "0 auto" }}>
        <Card elevated padding="md">
          {/* Header rail */}
          <div
            style={{
              display: "flex",
              alignItems: "flex-start",
              justifyContent: "space-between",
              marginBottom: 12,
              gap: 12,
            }}
          >
            {/* Left: kicker + title + subtitle stack */}
            <div style={{ display: "flex", flexDirection: "column", gap: 2 }}>
              <Kicker tone="accent" as="div">Integration</Kicker>
              <div className="text-display" style={{ lineHeight: 1.15 }}>
                exp 65 · Ia3d
              </div>
              <div className="text-data" style={{ opacity: 0.7 }}>
                smp_65 · SSRL Apr 2026
              </div>
            </div>
            {/* Right: scale toggle + add-peak button */}
            <div style={{ display: "flex", alignItems: "center", gap: 8, flexShrink: 0 }}>
              <SegmentedControl
                aria-label="Scale"
                options={[
                  { value: "log", label: "Log" },
                  { value: "linear", label: "Linear" },
                ]}
                value={scale}
                onChange={(v) => setScale(v)}
                size="xs"
              />
              <Button variant="ghost" onClick={() => console.log("addPeak")}>
                + Peak
              </Button>
            </div>
          </div>
          {/* Plot body */}
          <TracePlot
            trace={heroModel}
            height={360}
            paperColor="var(--color-paper)"
            yType={scale}
            interaction={{
              onXDomain: () => {},
              onAddPeak: () => {},
              onClickPeak: () => {},
            }}
          />
        </Card>
      </div>
    );
  },
};

// Annotated: live layer toggles — peaks, labels, confidence band.
// All three layers can be toggled independently via Storybook controls.
export const Annotated: StoryObj<{ showPeaks: boolean; showLabels: boolean; showBand: boolean }> = {
  args: { showPeaks: true, showLabels: true, showBand: true },
  argTypes: {
    showPeaks: { control: "boolean" },
    showLabels: { control: "boolean" },
    showBand: { control: "boolean" },
  },
  render: ({ showPeaks, showLabels, showBand }) => (
    <TracePlot
      trace={annotatedModel}
      height={360}
      show={{ peaks: showPeaks, labels: showLabels, band: showBand }}
      interaction={{ onXDomain: () => {}, onAddPeak: () => {}, onClickPeak: () => {} }}
      paperColor="var(--color-paper)"
    />
  ),
};

// Hover: manual-fidelity story — pointer or keyboard-focus a peak to see the
// terracotta ring, q-line, and axis chip; other peaks do NOT dim.
export const Hover: Story = {
  args: {
    trace: annotatedModel,
    height: 360,
    show: { labels: true },
    interaction: {
      onXDomain: () => {},
      onAddPeak: () => {},
      onClickPeak: () => {},
    },
    paperColor: "var(--color-paper)",
  },
  parameters: {
    docs: {
      description: {
        story:
          "Pointer or keyboard-focus a peak: it grows + gains the terracotta ring, the q-line shows, and a q-value chip appears at the axis foot. Other peaks do NOT dim.",
      },
    },
  },
};

// PhaseHighlight: a subset of peaks stays at full colour; the rest fade to neutral.
export const PhaseHighlight: Story = {
  args: {
    trace: annotatedModel,
    height: 360,
    show: { labels: true },
    highlightPeakIds: new Set(annotatedModel.peaks.slice(0, 3).map((p) => p.id)),
  },
};
