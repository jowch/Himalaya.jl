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
  return (
    <div style={{ maxWidth: 1180 }}>
      <TracePlate
        kicker="Integration"
        title="Lipid 1-2 + LL37 1:0.5"
        subtitle="smp_09 · SSRL Apr 2026 · representative exposure smp_09_e03"
        trace={heroModel}
        scale={scale}
        onScaleChange={setScale}
        onAutoFit={() => {}}
        addPeakArmed={armed}
        onToggleAddPeak={() => setArmed((p) => !p)}
        interaction={{ onXDomain: () => {} }}
      />
    </div>
  );
}

export const Hero: Story = { render: () => <HeroDemo /> };
