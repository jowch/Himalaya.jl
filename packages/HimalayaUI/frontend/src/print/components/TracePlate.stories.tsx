import { useState } from "react";
import type { Meta, StoryObj } from "@storybook/react-vite";
import { TracePlate, type TraceScale } from "./TracePlate";
import type { TraceModel } from "../plot/TracePlot";
import type { PlotPeak } from "../plot/marks/PlotPeaks";
import { realTraces } from "../fixtures/realTraces";
import { realMembers } from "../fixtures/realSeriesMembers";
import type { SeriesMember } from "../../api";

function modelFor(member: SeriesMember): TraceModel {
  const trace = realTraces[member.exposure_id as number]!;
  const peaks: PlotPeak[] = (member.snapshot?.effective_peaks ?? []).map((p) => ({
    id: p.id,
    q: p.q,
    intensity: p.intensity,
    source: p.source,
  }));
  return { trace, peaks, phase: member.snapshot?.confirmed_index?.phase ?? null };
}

const heroModel = modelFor(realMembers[0]!); // exp 65

const meta: Meta<typeof TracePlate> = {
  title: "components/TracePlate",
  component: TracePlate,
  parameters: { layout: "padded" },
};
export default meta;
type Story = StoryObj<typeof meta>;

function HeroDemo() {
  const [scale, setScale] = useState<TraceScale>("log");
  const [armed, setArmed] = useState(false);
  // Scroll-to-zoom round-trip: the wheel emits a window via onXDomain; we store
  // it and feed it straight back as the controlled `xDomain` so the zoom renders.
  // Double-click (TracePlot reset) emits null → back to full extent. Auto-fit
  // here clears the zoom too.
  const [zoom, setZoom] = useState<[number, number] | null>(null);
  return (
    <div style={{ maxWidth: 1180 }}>
      <TracePlate
        kicker="Integration"
        title="Lipid 1-2 + LL37 1:0.5"
        subtitle="smp_09 · SSRL Apr 2026 · representative exposure smp_09_e03"
        trace={heroModel}
        scale={scale}
        onScaleChange={setScale}
        onAutoFit={() => setZoom(null)}
        addPeakArmed={armed}
        onToggleAddPeak={() => setArmed((p) => !p)}
        xDomain={zoom}
        interaction={{ onXDomain: setZoom }}
      />
    </div>
  );
}

export const Hero: Story = { render: () => <HeroDemo /> };

/** A pre-zoomed window. The trace + any peaks/labels whose q falls outside
 *  [0.04, 0.08] clip cleanly at the axes instead of overdrawing the spines and
 *  tick labels (peak labels keep their headroom above the curve). */
export const Zoomed: Story = {
  render: () => (
    <div style={{ maxWidth: 1180 }}>
      <TracePlate
        kicker="Integration"
        title="Lipid 1-2 + LL37 1:0.5"
        subtitle="zoomed to q ∈ [0.04, 0.08] — annotations clip at the axes"
        trace={heroModel}
        scale="log"
        onScaleChange={() => {}}
        xDomain={[0.04, 0.08]}
        interaction={{ onXDomain: () => {} }}
      />
    </div>
  ),
};
